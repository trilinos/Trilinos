#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "FaceCreatorFixture.hpp"
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace
{

class FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester : public FaceCreatorFixture
{
protected:
    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        set_bulk(new stk::unit_test_util::BulkDataElemGraphFaceSharingTester(get_meta(), get_comm(), auraOption));
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
