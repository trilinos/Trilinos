#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "FaceCreatorFixture.hpp"

namespace
{

class TwoHexGeneratedMeshWithFaceAdjacentElementGraph : public FaceCreatorFixture
{
protected:
  void setup_mesh_and_create_graph(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x2", auraOption);
  }
};


TEST_F(TwoHexGeneratedMeshWithFaceAdjacentElementGraph, onlyOneFaceExistsAfterBothProcsCreateFaceOnProcBoundary_aura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh_and_create_graph(stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }
}

TEST_F(TwoHexGeneratedMeshWithFaceAdjacentElementGraph, onlyOneFaceExistsAfterBothProcsCreateFaceOnProcBoundary_noAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh_and_create_graph(stk::mesh::BulkData::NO_AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }
}

TEST_F(TwoHexGeneratedMeshWithFaceAdjacentElementGraph, faceExistsOnBothProcsAfterOnlyOneProcCreatesFace_aura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh_and_create_graph(stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
  }
}

TEST_F(TwoHexGeneratedMeshWithFaceAdjacentElementGraph, faceExistsOnBothProcsAfterOnlyOneProcCreatesFace_noAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_mesh_and_create_graph(stk::mesh::BulkData::NO_AUTO_AURA);
    test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
  }
}


}
