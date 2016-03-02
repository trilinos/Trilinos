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
    ElemElemGraphTester(stk::mesh::BulkData& bulkData, const stk::mesh::Selector &selector, const stk::mesh::Selector *air = nullptr) :
        stk::mesh::ElemElemGraph(bulkData, selector, air) { }

    stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
};

void pack_data_for_air(stk::CommSparse &comm, ElemElemGraphTester& graph, const stk::mesh::BulkData& bulkData, stk::mesh::Selector sel)
{
    stk::mesh::impl::ParallelGraphInfo& parallel_info = graph.get_parallel_info();
    for(const auto& item : parallel_info)
    {
        const stk::mesh::GraphEdge &edge = item.first;
        const stk::mesh::impl::ParallelInfo &pinfo = item.second;
        stk::mesh::Entity local_element = graph.get_entity(edge.elem1);

        bool didPartChangesOccur = true;

        if(didPartChangesOccur)
        {
            int64_t id = bulkData.identifier(local_element);
            if(sel(bulkData.bucket(local_element)))
                comm.send_buffer(pinfo.get_proc_rank_of_neighbor()).pack<int64_t>(id);
            else
                comm.send_buffer(pinfo.get_proc_rank_of_neighbor()).pack<int64_t>(-id);
        }
    }
}

void unpack_and_update_air(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, ElemElemGraphTester& graph)
{
    std::vector<int64_t> idsTrue;
    std::vector<int64_t> idsFalse;
    for(int i=0;i<bulkData.parallel_size();++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            int64_t id;
            comm.recv_buffer(i).unpack<int64_t>(id);
            if(id>0)
                idsTrue.push_back(id);
            else
                idsFalse.push_back(-id);
        }
    }

    std::sort(idsTrue.begin(), idsTrue.end());
    std::sort(idsFalse.begin(), idsFalse.end());

    stk::mesh::impl::ParallelGraphInfo &parallel_info = graph.get_parallel_info();
    for(auto& item : parallel_info)
    {
        const stk::mesh::GraphEdge &edge = item.first;
        stk::mesh::impl::ParallelInfo &pinfo = item.second;
        int64_t id = -edge.elem2;
        if(std::binary_search(idsTrue.begin(), idsTrue.end(), id))
        {
            pinfo.set_is_in_air(true);
        }
        else if(std::binary_search(idsFalse.begin(), idsFalse.end(), id))
            pinfo.set_is_in_air(false);
    }
}

void update_parallel_graph_for_air_selector(ElemElemGraphTester& graph, const stk::mesh::BulkData& bulkData, stk::mesh::Selector air)
{
    stk::CommSparse comm(bulkData.parallel());
    pack_data_for_air(comm, graph, bulkData, air);
    comm.allocate_buffers();
    pack_data_for_air(comm, graph, bulkData, air);
    comm.communicate();
    unpack_and_update_air(comm, bulkData, graph);
}

class ParallelGraphUpdate : public stk::unit_test_util::MeshFixture
{
public:
    ParallelGraphUpdate() : stk::unit_test_util::MeshFixture()
    {
        part = &get_meta().declare_part_with_topology("air", stk::topology::HEXAHEDRON_8);
        air = *part;
        make_part_an_element_block_by_setting_id();
    }

    void make_some_elements_air()
    {
        stk::mesh::EntityVector elements = get_half_of_elements();
        add_air_to_elements(elements);
    }

    void add_air_to_elements(const stk::mesh::EntityVector &elements)
    {
        get_bulk().modification_begin();
        for(stk::mesh::Entity element : elements )
            get_bulk().change_entity_parts(element, {part}, {});
        get_bulk().modification_end();
    }

    void remove_air_from_elements(const stk::mesh::EntityVector &elements)
    {
        get_bulk().modification_begin();
        for(stk::mesh::Entity element : elements )
            get_bulk().change_entity_parts(element, {}, {part});
        get_bulk().modification_end();
    }

    stk::mesh::Selector get_air() const
    {
        return air;
    }

    void move_elements_out_of_air()
    {
        remove_air_from_elements(get_half_of_elements());
    }

    void move_elements_into_air_and_verify()
    {
        size_t num_elements_air = stk::mesh::count_selected_entities(get_air(), get_bulk().buckets(stk::topology::ELEM_RANK));
        EXPECT_TRUE(num_elements_air==0u);
        make_some_elements_air();
        num_elements_air = stk::mesh::count_selected_entities(get_air(), get_bulk().buckets(stk::topology::ELEM_RANK));
        EXPECT_TRUE(num_elements_air!=0u);
    }

    void verify_no_air_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info) const
    {
        for(const auto& item : parallel_info)
        {
            const stk::mesh::impl::ParallelInfo &pinfo = item.second;
            EXPECT_TRUE(pinfo.is_considered_air() == false);
        }
    }

    void verify_no_air_part_ordinal_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
    {
        for(const auto& item : parallelPartInfo)
        {
            const stk::mesh::impl::PartOrdinals &partOrdinals = item.second;
            EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), part->mesh_meta_data_ordinal()));
        }
    }

    void verify_air_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info) const
    {
        for(const auto& item : parallel_info)
        {
            const stk::mesh::GraphEdge &edge = item.first;
            const stk::mesh::impl::ParallelInfo& pinfo = item.second;
            if(edge.elem2 %2 == 0)
                EXPECT_TRUE(pinfo.is_considered_air() == true) << " for element " << -edge.elem2;
            else
                EXPECT_TRUE(pinfo.is_considered_air() == false) << " for element " << -edge.elem2;
        }
    }

    void verify_air_part_ordinal_is_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
    {
        for(const auto& item : parallelPartInfo)
        {
            const stk::mesh::GraphEdge &edge = item.first;
            const stk::mesh::impl::PartOrdinals &partOrdinals = item.second;
            if(edge.elem2%2 == 0)
                EXPECT_TRUE(std::binary_search(partOrdinals.begin(), partOrdinals.end(), part->mesh_meta_data_ordinal()));
            else
                EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), part->mesh_meta_data_ordinal()));
        }
    }

private:

    void make_part_an_element_block_by_setting_id()
    {
        get_meta().set_part_id(*part, 100);
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

    stk::mesh::Part* part = nullptr;
    stk::mesh::Selector air;
};

// Need updater for part ordinals

TEST_F(ParallelGraphUpdate, updateAirSelector)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD)==4)
    {
        setup_mesh("generated:4x4x4", stk::mesh::BulkData::NO_AUTO_AURA);
        ElemElemGraphTester graph(get_bulk(), get_meta().locally_owned_part());

        stk::mesh::impl::ParallelPartInfo parallelPartInfo;
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);

        verify_no_air_in_parallel_graph(graph.get_parallel_info());
        verify_no_air_part_ordinal_in_parallel_graph(parallelPartInfo);

        move_elements_into_air_and_verify();
        update_parallel_graph_for_air_selector(graph, get_bulk(), get_air()); // Note: would work for skinned selector also!
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
        verify_air_in_parallel_graph(graph.get_parallel_info());
        verify_air_part_ordinal_is_in_parallel_graph(parallelPartInfo);

        move_elements_out_of_air();
        update_parallel_graph_for_air_selector(graph, get_bulk(), get_air());
        stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
        verify_no_air_in_parallel_graph(graph.get_parallel_info());
        verify_no_air_part_ordinal_in_parallel_graph(parallelPartInfo);
    }
}

