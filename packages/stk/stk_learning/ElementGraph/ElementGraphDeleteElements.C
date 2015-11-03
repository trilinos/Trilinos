#include <gtest/gtest.h>
#include <vector>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include "ElementGraphTester.hpp"
namespace
{

typedef std::map<stk::mesh::EntityId, stk::mesh::impl::IdViaSidePair> MockGraphEdges;

class ElemGraphDeleteElementsTester : public stk::unit_test_util::MeshFixture
{
protected:
    ElemGraphDeleteElementsTester() : elementGraph(nullptr) { }
    ~ElemGraphDeleteElementsTester()
    {
        delete elementGraph;
    }

    void create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                                    stk::mesh::EntityIdVector elemIdsToDelete,
                                                                    size_t goldNumGlobalConnectionsAfterDeletes,
                                                                    MockGraphEdges &goldConnections)
    {
        setup_mesh("generated:1x1x4", auraOption);
        elementGraph = new ElemElemGraphTester(get_bulk());
        test(elemIdsToDelete, goldNumGlobalConnectionsAfterDeletes, goldConnections);
    }

    void test(const stk::mesh::EntityIdVector &elemIdsToDelete, size_t goldNumGlobalConnectionsAfterDeletes, MockGraphEdges &goldConnections)
    {
        remove_elements(elemIdsToDelete);

        check_graph_connections_updated_after_deletions(goldConnections, goldNumGlobalConnectionsAfterDeletes);
    }

    void remove_elements(const stk::mesh::EntityIdVector &elemIdsToDelete)
    {
        remove_inactive_elements(elemIdsToDelete);
        stk::mesh::impl::DeletedElementInfoVector elementsToDelete;
        fill_elements_to_delete(elemIdsToDelete, elementsToDelete);
        destroy_elements(elementsToDelete);
        elementGraph->delete_elements(elementsToDelete);
    }

    void expect_local_elem_counts_are_equal_to_active_elements(stk::mesh::Selector locallyOwned, size_t numActiveElements)
    {
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(get_bulk(), counts, &locallyOwned);
        EXPECT_EQ(numActiveElements, counts[stk::topology::ELEM_RANK]);
    }

    void expect_correct_number_of_connections_in_graph_after_deletes(stk::mesh::Selector locallyOwned,
                                                                     size_t goldNumGlobalConnectionsAfterDeletes,
                                                                     stk::mesh::EntityVector &elems)
    {
        size_t localNumConnections = 0;
        for(stk::mesh::Entity elem : elems)
        {
            localNumConnections += elementGraph->get_num_connected_elems(elem);
        }
        size_t globalNumConnections = 0;
        stk::all_reduce_sum(get_comm(), &localNumConnections, &globalNumConnections, 1);
        EXPECT_EQ(goldNumGlobalConnectionsAfterDeletes, globalNumConnections);
    }

    void expect_connection_is_in_gold_connections(const MockGraphEdges &goldConnections, stk::mesh::EntityId elemId, stk::mesh::impl::IdViaSidePair connection)
    {
        MockGraphEdges::const_iterator it = goldConnections.find(elemId);
        if(it != goldConnections.end())
        {
            stk::mesh::impl::IdViaSidePair goldConnection = it->second;
            EXPECT_EQ(goldConnection.id, connection.id);
            EXPECT_EQ(goldConnection.side, connection.side);
        }
    }

    void remove_inactive_elements(const stk::mesh::EntityIdVector &elemIdsToDelete)
    {
        for(stk::mesh::EntityId elemId : elemIdsToDelete)
        {
            activeElements.erase(elemId);
        }
    }

    void fill_elements_to_delete(const stk::mesh::EntityIdVector &elemIdsToDelete, stk::mesh::impl::DeletedElementInfoVector &elements_to_delete)
    {
        for(stk::mesh::EntityId elemId : elemIdsToDelete)
        {
            stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
            if(get_bulk().is_valid(elem) && get_bulk().bucket(elem).owned())
            {
                elements_to_delete.push_back( {elem, elemId, get_bulk().bucket(elem).topology().is_shell()});
            }
        }
    }

    void destroy_elements(stk::mesh::impl::DeletedElementInfoVector &elementsToDelete)
    {
        get_bulk().modification_begin();
        for(stk::mesh::impl::DeletedElementInfo elem : elementsToDelete)
        {
            get_bulk().destroy_entity(elem.entity);
            stk::mesh::Entity elemCheck = get_bulk().get_entity(stk::topology::ELEM_RANK, elem.identifier);
            EXPECT_FALSE(get_bulk().is_valid(elemCheck));
        }
        get_bulk().modification_end();
    }

    void check_graph_connections_updated_after_deletions(const MockGraphEdges &goldConnections, size_t goldNumGlobalConnectionsAfterDeletes)
    {
        expect_local_elem_counts_are_equal_to_active_elements(get_meta().locally_owned_part(), activeElements.size());

        stk::mesh::EntityVector elems;
        stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elems);
        expect_correct_number_of_connections_in_graph_after_deletes(get_meta().locally_owned_part(), goldNumGlobalConnectionsAfterDeletes, elems);

        for(stk::mesh::Entity elem : elems)
            for(size_t i = 0; i < elementGraph->get_num_connected_elems(elem); i++)
                expect_connection_is_in_gold_connections(goldConnections, get_bulk().identifier(elem), get_connection(elem, i));
    }

    stk::mesh::impl::IdViaSidePair get_connection(stk::mesh::Entity elem, size_t sideOffset)
    {
        if(elementGraph->is_connected_elem_locally_owned(elem, sideOffset))
            return get_local_connection_information(elem, sideOffset);
        else
            return elementGraph->get_connected_remote_id_and_via_side(elem, sideOffset);
    }

    stk::mesh::impl::IdViaSidePair get_local_connection_information(stk::mesh::Entity elem, size_t sideOffset)
    {
        stk::mesh::impl::ElementViaSidePair conn = elementGraph->get_connected_element_and_via_side(elem, sideOffset);
        return {get_bulk().identifier(elem), conn.side};
    }

protected:
    ElemElemGraphTester *elementGraph;
    std::set<stk::mesh::EntityId> activeElements = {1, 2, 3, 4};
};

TEST_F(ElemGraphDeleteElementsTester, deleteElements1and4WithAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {{4, {3, 2}}, {5, {2, 3}}};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::AUTO_AURA, {1, 4}, 2, goldConnections);
    }
}
TEST_F(ElemGraphDeleteElementsTester, deleteElements1and4WithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {{4, {3, 2}}, {5, {2, 3}}};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::NO_AUTO_AURA, {1, 4}, 2, goldConnections);
    }
}
TEST_F(ElemGraphDeleteElementsTester, deleteElements1and3WithAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::AUTO_AURA, {1, 3}, 0, goldConnections);
    }
}
TEST_F(ElemGraphDeleteElementsTester, deleteElements1and3WithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::NO_AUTO_AURA, {1, 3}, 0, goldConnections);
    }
}
TEST_F(ElemGraphDeleteElementsTester, deleteEveryElementWithAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::AUTO_AURA, {1, 2, 3, 4}, 0, goldConnections);
    }
}
TEST_F(ElemGraphDeleteElementsTester, deleteEveryElementWithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) < 4)
    {
        MockGraphEdges goldConnections = {};
        create_and_expect_graph_updated_after_elements_are_deleted(stk::mesh::BulkData::NO_AUTO_AURA, {1, 2, 3, 4}, 0, goldConnections);
    }
}
//        create_and_expect_graph_updated_after_elements_are_deleted({2, 3}, stk::mesh::BulkData::AUTO_AURA);
//        create_and_expect_graph_updated_after_elements_are_deleted({1, 3}, stk::mesh::BulkData::AUTO_AURA);
//        create_and_expect_graph_updated_after_elements_are_deleted({2, 4}, stk::mesh::BulkData::AUTO_AURA);
//        create_and_expect_graph_updated_after_elements_are_deleted({}, stk::mesh::BulkData::AUTO_AURA);
//        create_and_expect_graph_updated_after_elements_are_deleted({1, 2, 3, 4}, stk::mesh::BulkData::AUTO_AURA);

}
