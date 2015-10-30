#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include "FaceCreatorFixture.hpp"

namespace
{

class FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester : public FaceCreatorFixture
{
protected:
    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        bulkData = new stk::mesh::unit_test::BulkDataElemGraphFaceSharingTester(metaData, communicator, auraOption);
    }

    void test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face()
    {
        only_proc_0_makes_a_face();
        test_that_num_sides_is_expected_value(1);
        test_that_each_proc_has_num_sides_with_expected_value(1);
    }

    virtual void create_faces(stk::mesh::Entity element, stk::mesh::EntityVector& nodes_of_face)
    {
        get_bulk().modification_begin();
        if(get_bulk().parallel_rank()==0)
        {
            create_face_per_proc(element, nodes_of_face);
        }
        test_that_num_sides_is_expected_value(1);
        get_bulk().modification_end();
    }

    void only_proc_0_makes_a_face()
    {
        unsigned id = get_bulk().parallel_rank()+1;
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        stk::mesh::EntityVector nodes_of_face = get_nodes_of_face_for_this_proc();
        create_faces(elem, nodes_of_face);
    }

    void test_that_each_proc_has_num_sides_with_expected_value(unsigned expected_num_sides)
    {
        unsigned num_local_sides = stk::mesh::count_selected_entities(get_bulk().mesh_meta_data().globally_shared_part(), get_bulk().buckets(get_bulk().mesh_meta_data().side_rank()));
        EXPECT_EQ(expected_num_sides, num_local_sides);
    }
};


TEST_F(FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}

TEST_F(FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}


}
