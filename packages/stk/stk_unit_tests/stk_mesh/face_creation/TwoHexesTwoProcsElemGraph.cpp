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
};


TEST_F(FaceCreatorElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
        test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
    }
}


}
