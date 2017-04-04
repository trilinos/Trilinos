#ifndef FACECREATORFIXTURE_HPP_
#define FACECREATORFIXTURE_HPP_

#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/FaceTestingUtils.hpp>

class FaceCreatorFixture : public stk::unit_test_util::MeshFixture
{
protected:

    FaceCreatorFixture(unsigned spatial_dim = 3) : MeshFixture(spatial_dim) {}

    virtual ~FaceCreatorFixture() { }

    void test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary()
    {
        each_proc_make_face_on_proc_boundary();
        test_that_num_sides_is_expected_value(1);
    }

    void each_proc_make_face_on_proc_boundary()
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        stk::mesh::EntityVector nodes_of_face = get_nodes_of_face_for_this_proc();
        create_faces(elem, nodes_of_face);
    }

    virtual stk::mesh::EntityVector get_nodes_of_face_for_this_proc()
    {
        std::vector<unsigned> face_node_ids = { 8, 7, 5, 6 };
        return get_nodes_for_proc(face_node_ids);
    }

    void create_faces(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        get_bulk().modification_begin();
        create_face_per_proc(element, nodes_of_face);
        test_that_num_sides_is_expected_value(2);
        get_bulk().modification_end();
    }

    virtual void test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face()
    {
        only_proc_0_makes_a_face();
        test_that_num_sides_is_expected_value(1);
        test_that_each_proc_has_num_sides_with_expected_value(1);
    }

    void only_proc_0_makes_a_face()
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        stk::mesh::EntityVector nodes_of_face = get_nodes_of_face_for_this_proc();
        create_faces_only_one_proc(elem, nodes_of_face);
    }

    virtual void create_faces_only_one_proc(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        get_bulk().modification_begin();
        if(get_bulk().parallel_rank()==0)
        {
            create_face_per_proc(element, nodes_of_face);
        }
        test_that_num_sides_is_expected_value(1);
        get_bulk().modification_end();
    }

    void create_face_per_proc(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::topology side_topology = get_bulk().bucket(element).topology().side_topology();
        stk::mesh::Entity side = stk::unit_test_util::declare_element_side_with_nodes(get_bulk(), element, nodes_of_face, id, get_meta().get_topology_root_part(side_topology));
        EXPECT_TRUE(get_bulk().is_valid(side));
    }

    virtual unsigned get_permuted_index(unsigned i)
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

    void test_that_each_proc_has_num_sides_with_expected_value(unsigned expected_num_sides)
    {
        unsigned num_local_sides = stk::mesh::count_selected_entities(get_bulk().mesh_meta_data().globally_shared_part(), get_bulk().buckets(get_bulk().mesh_meta_data().side_rank()));
        EXPECT_EQ(expected_num_sides, num_local_sides);
    }
};

#endif /* FACECREATORFIXTURE_HPP_ */
