#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include "FaceCreatorFixture.hpp"

namespace
{

class FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester : public FaceCreatorFixture
{
protected:

    FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester() : FaceCreatorFixture(2) {}

    void setup_2x1_2d_mesh(stk::mesh::BulkData::AutomaticAuraOption aura_option)
    {
        bulkData = new stk::mesh::unit_test::BulkDataElemGraphFaceSharingTester(metaData, get_comm(), aura_option);
        unsigned numX = 2, numY = 1;
        stk::unit_test_util::convert_quad_fixture_to_my_bulk_data_flavor(numX, numY, bulkData);
    }

    virtual stk::mesh::EntityVector get_nodes_of_face_for_this_proc()
    {
        std::vector<unsigned> face_node_ids = { 2, 5 };
        return get_nodes_for_proc(face_node_ids);
    }

    virtual unsigned get_permuted_index(unsigned i)
    {
        std::vector<std::vector<unsigned> > index_for_proc = {
                {0, 1},
                {1, 0}
        };
        return index_for_proc[get_bulk().parallel_rank()][i];
    }
};


TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_2x1_2d_mesh(stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_2x1_2d_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_2x1_2d_mesh(stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithoutAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_2x1_2d_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
    }
}

}
