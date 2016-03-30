#include "gtest/gtest.h"
#include <mpi.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>

class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
    ElemElemGraphTester(stk::mesh::BulkData& bulkData) :
        stk::mesh::ElemElemGraph(bulkData) { }

    stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
};


class ParallelGraphUpdate : public stk::unit_test_util::MeshFixture
{
public:
    ParallelGraphUpdate() : stk::unit_test_util::MeshFixture()
    {
        airPart = &get_meta().declare_part_with_topology("air", stk::topology::HEXAHEDRON_8);
        air = *airPart;

        skinPart = &get_meta().declare_part_with_topology("things to skin", stk::topology::HEXAHEDRON_8);
        skinSelector = *skinPart;

        make_part_an_element_block_by_setting_id();
    }

    void put_some_elements_into_some_parts()
    {
        stk::mesh::EntityVector elements = get_half_of_elements();
        add_part_to_elements(elements,airPart);
        add_part_to_elements(elements,skinPart);
    }

    void add_part_to_elements(const stk::mesh::EntityVector &elements, stk::mesh::Part* part)
    {
        get_bulk().modification_begin();
        for(stk::mesh::Entity element : elements )
            get_bulk().change_entity_parts(element, {part}, {});
        get_bulk().modification_end();
    }

    void remove_part_from_elements(const stk::mesh::EntityVector &elements, stk::mesh::Part* part)
    {
        get_bulk().modification_begin();
        for(stk::mesh::Entity element : elements )
            get_bulk().change_entity_parts(element, {}, {part});
        get_bulk().modification_end();
    }

    stk::mesh::Selector get_air_selector() const
    {
        return air;
    }

    stk::mesh::Selector get_skin_selector() const
    {
        return skinSelector;
    }

    void move_elements_out_of_parts()
    {
        remove_part_from_elements(get_half_of_elements(), airPart);
        remove_part_from_elements(get_half_of_elements(), skinPart);
    }

    void move_elements_into_parts_and_verify()
    {
        expect_no_elements_selected(get_air_selector());
        expect_no_elements_selected(get_skin_selector());
        put_some_elements_into_some_parts();
        expect_some_elements_selected(get_air_selector());
        expect_some_elements_selected(get_skin_selector());
    }

    void verify_selectors_not_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info)
    {
        for(const auto& item : parallel_info)
        {
            stk::mesh::impl::LocalId elemOnOtherProc = item.first.elem2;
            EXPECT_FALSE(remoteAirSelector[elemOnOtherProc]);
            EXPECT_FALSE(remoteSkinSelector[elemOnOtherProc]);
        }
    }

    void verify_part_ordinals_not_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
    {
        for(const auto& item : parallelPartInfo)
        {
            const stk::mesh::impl::PartOrdinals &partOrdinals = item.second;
            EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
            EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
        }
    }

    void verify_selectors_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info)
    {
        for(const auto& item : parallel_info)
        {
            stk::mesh::impl::LocalId elemOnOtherProc = item.first.elem2;
            if(elemOnOtherProc %2 == 0)
            {
                EXPECT_TRUE(remoteAirSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
                EXPECT_TRUE(remoteSkinSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
            }
            else
            {
                EXPECT_FALSE(remoteAirSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
                EXPECT_FALSE(remoteSkinSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
            }
        }
    }

    void verify_part_ordinals_are_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
    {
        for(const auto& item : parallelPartInfo)
        {
            const stk::mesh::impl::LocalId elemOnOtherProc = item.first;
            const stk::mesh::impl::PartOrdinals &partOrdinals = item.second;
            if(elemOnOtherProc%2 == 0)
            {
                EXPECT_TRUE(std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
                EXPECT_TRUE(std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
            }
            else
            {
                EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
                EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
            }
        }
    }


    void make_part_an_element_block_by_setting_id()
    {
        get_meta().set_part_id(*airPart, 100);
        get_meta().set_part_id(*skinPart, 101);
    }

    stk::mesh::EntityVector get_half_of_elements()
    {
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elements);
        return get_elements_that_have_even_identifiers(elements);
    }

    stk::mesh::EntityVector get_elements_that_have_even_identifiers(const stk::mesh::EntityVector& elements)
    {
        stk::mesh::EntityVector someElements;
        for(stk::mesh::Entity element : elements)
            if(get_bulk().identifier(element)%2 == 0)
                someElements.push_back(element);
        return someElements;
    }

    void expect_no_elements_selected(stk::mesh::Selector selector)
    {
        size_t num_elements = stk::mesh::count_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK));
        EXPECT_EQ(0u, num_elements);
    }

    void expect_some_elements_selected(stk::mesh::Selector selector)
    {
        size_t num_elements= stk::mesh::count_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK));
        EXPECT_TRUE(num_elements != 0u);
    }

    stk::mesh::Part* airPart = nullptr;
    stk::mesh::Part* skinPart = nullptr;
    stk::mesh::Selector air;
    stk::mesh::Selector skinSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteSkinSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteAirSelector;
};

// Need updater for part ordinals

TEST_F(ParallelGraphUpdate, updateAirSelector)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD)==4)
    {
        setup_mesh("generated:4x4x4", stk::mesh::BulkData::NO_AUTO_AURA);
        ElemElemGraphTester graph(get_bulk());

        stk::mesh::impl::ParallelPartInfo parallelPartInfo;
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
        stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);

        verify_selectors_not_in_parallel_graph(graph.get_parallel_info());
        verify_part_ordinals_not_in_parallel_graph(parallelPartInfo);

        move_elements_into_parts_and_verify();
        stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_air_selector(), remoteAirSelector);
        stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
        verify_selectors_in_parallel_graph(graph.get_parallel_info());
        verify_part_ordinals_are_in_parallel_graph(parallelPartInfo);

        move_elements_out_of_parts();
        stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_air_selector(), remoteAirSelector);
        stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
        verify_selectors_not_in_parallel_graph(graph.get_parallel_info());
        verify_part_ordinals_not_in_parallel_graph(parallelPartInfo);
    }
}

