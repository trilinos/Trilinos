#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>

#include <unit_tests/BulkDataTester.hpp>

namespace
{

class FaceCreator : public stk::unit_test_util::MeshFixture
{
protected:

    void test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary()
    {
        each_proc_make_face_on_proc_boundary();
        test_that_num_sides_is_expected_value(1);
    }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        bulkData = new stk::mesh::unit_test::BulkDataFaceSharingTester(metaData, communicator, auraOption);
    }

    void test_elem_elem_graph_for_face_connection_info()
    {
        stk::mesh::ElemElemGraph egraph(get_bulk(), get_meta().universal_part());
        test_sides_using_graph(egraph);
    }

private:

    void test_sides_using_graph(const stk::mesh::ElemElemGraph & egraph)
    {
        stk::mesh::EntityVector sides = get_shared_sides();
        for(size_t i = 0;i < sides.size();++i)
            test_graph_info_on_side(egraph, sides[i]);
    }

    void each_proc_make_face_on_proc_boundary()
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        stk::mesh::EntityVector nodes_of_face = get_nodes_of_face_for_this_proc();
        create_faces(elem, nodes_of_face);
    }

    stk::mesh::EntityVector get_nodes_of_face_for_this_proc()
    {
        //        std::vector<unsigned> face_node_ids = { 5, 6, 8, 7 };
        std::vector<unsigned> face_node_ids = { 8, 7, 5, 6 };// 6, 5, 7, 8
        return get_nodes_for_proc(face_node_ids);
    }

    void create_faces(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        get_bulk().modification_begin();
        create_face_per_proc(element, nodes_of_face);
        get_bulk().modification_end();
    }

    void create_face_per_proc(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity side = stk::mesh::declare_element_to_sub_topology_with_nodes(get_bulk(), element, nodes_of_face, id, stk::topology::FACE_RANK,
                get_meta().get_topology_root_part(stk::topology::QUAD_4_2D));
        EXPECT_TRUE(get_bulk().is_valid(side));
        test_that_num_sides_is_expected_value(2);
    }

    unsigned get_permuted_index(unsigned i)
    {
        std::vector<std::vector<unsigned> > index_for_proc = {
                {0, 1, 2, 3},
                {3, 2, 1, 0}
        };
        return index_for_proc[get_bulk().parallel_rank()][i];
    }

    stk::mesh::EntityVector get_nodes_for_proc(const std::vector<unsigned>& face_node_ids)
    {
        stk::mesh::EntityVector nodes(face_node_ids.size());
        for(size_t n = 0; n < nodes.size(); ++n)
            nodes[n] = get_bulk().get_entity(stk::topology::NODE_RANK, face_node_ids[get_permuted_index(n)]);
        return nodes;
    }


    void test_that_num_sides_is_expected_value(size_t num_sides_gold)
    {
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(get_bulk(), counts);
        EXPECT_EQ(num_sides_gold, counts[get_meta().side_rank()]);
    }

    stk::mesh::EntityVector get_shared_sides()
    {
        stk::mesh::EntityVector sides;
        stk::mesh::get_selected_entities(get_meta().globally_shared_part(), get_bulk().buckets(get_meta().side_rank()), sides);
        EXPECT_EQ(1u, sides.size());
        return sides;
    }

    void test_graph_info_on_side(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity side)
    {
        unsigned num_elements = get_bulk().num_elements(side);
        const stk::mesh::Entity* elements = get_bulk().begin_elements(side);
        for(unsigned j=0;j<num_elements;++j)
            check_if_remote_sides_have_consistent_parallel_info(egraph, elements[j], side);
    }

    void check_if_remote_sides_have_consistent_parallel_info(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity element, stk::mesh::Entity side)
    {
        if(get_bulk().bucket(element).owned())
        {
            int num_connected_elems = egraph.get_num_connected_elems(element);
            for(int k=0;k<num_connected_elems;++k)
                check_if_remote_side_has_consistent_parallel_info(egraph, element, k, side);
        }
    }

    void check_if_remote_side_has_consistent_parallel_info(const stk::mesh::ElemElemGraph& egraph, stk::mesh::Entity element, int element_index, stk::mesh::Entity side)
    {
        if(!egraph.is_connected_elem_locally_owned(element, element_index))
        {
            stk::mesh::impl::IdViaSidePair sidePair = egraph.get_connected_remote_id_and_via_side(element, element_index);
            const stk::mesh::impl::parallel_info &p_info = egraph.get_const_parallel_edge_info(element, sidePair.id);
            test_parallel_info_for_side(p_info, element, side, egraph.get_side_from_element1_to_remote_element2(element,sidePair.id));
        }
    }

    void test_parallel_info_for_side(const stk::mesh::impl::parallel_info &p_info, stk::mesh::Entity element, stk::mesh::Entity side, int side_ordinal)
    {
        test_parallel_info(p_info);
        stk::mesh::EntityVector side_nodes(get_bulk().begin_nodes(side),get_bulk().end_nodes(side));
        stk::mesh::EntityVector permuted_element_side_nodes = get_permuted_element_side_nodes(element, side_ordinal, p_info.m_other_proc, p_info.m_permutation);
        EXPECT_TRUE(side_nodes == permuted_element_side_nodes);
    }

    stk::mesh::EntityVector get_permuted_element_side_nodes(stk::mesh::Entity element, int side_ordinal, int other_proc_id, unsigned permutation_other_proc)
    {
        stk::mesh::EntityVector element_nodes(get_bulk().begin_nodes(element),get_bulk().end_nodes(element));
        stk::topology side_topology = get_side_topology(element, side_ordinal);
        stk::mesh::EntityVector element_side_nodes = get_element_side_nodes(element, element_nodes, side_ordinal);
        return get_permuted_nodes(side_topology, element_side_nodes, other_proc_id, permutation_other_proc);
    }

    stk::mesh::EntityVector get_element_side_nodes(stk::mesh::Entity element, const stk::mesh::EntityVector& element_nodes, int side_ordinal)
    {
        stk::topology element_topology = get_bulk().bucket(element).topology();
        stk::mesh::EntityVector element_side_nodes(get_side_topology(element, side_ordinal).num_nodes());
        element_topology.side_nodes(element_nodes, side_ordinal, element_side_nodes.begin());
        return element_side_nodes;
    }

    stk::topology get_side_topology(stk::mesh::Entity element, int side_ordinal)
    {
        stk::topology element_topology = get_bulk().bucket(element).topology();
        return element_topology.side_topology(side_ordinal);
    }

    stk::mesh::EntityVector get_permuted_nodes(stk::topology side_topology, const stk::mesh::EntityVector &element_side_nodes, int other_proc_id, unsigned permutation_other_proc)
    {
        stk::mesh::EntityVector permuted_element_side_nodes(side_topology.num_nodes());
        unsigned permutation = get_bulk().parallel_rank()>other_proc_id ? permutation_other_proc : 0;
        side_topology.permutation_nodes(element_side_nodes, permutation, permuted_element_side_nodes.begin());
        return permuted_element_side_nodes;
    }

    void test_parallel_info(const stk::mesh::impl::parallel_info &p_info)
    {
        EXPECT_EQ(1-get_bulk().parallel_rank(), p_info.m_other_proc);
        EXPECT_EQ(4, p_info.m_permutation);
        EXPECT_EQ(stk::topology::HEXAHEDRON_8, p_info.m_remote_element_toplogy);
        EXPECT_EQ(2u, p_info.m_chosen_side_id);
    }
};

TEST_F(FaceCreator, twoHexesTwoProcsCreateTwoFacesWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreator, DISABLED_testFaceDataUsingElemElemGraphWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
        test_elem_elem_graph_for_face_connection_info();
    }
}

TEST_F(FaceCreator, DISABLED_testFaceDataUsingElemElemGraphWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
        test_elem_elem_graph_for_face_connection_info();
    }
}

}
