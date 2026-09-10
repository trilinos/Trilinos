#include <gtest/gtest.h>

#include <Akri_CDMesh_Utils.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SmallClusterCleanup.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <Akri_Unit_InterfaceGeometry.hpp>

namespace krino {

template <typename MESHSPEC>
class SmallClusterCleanupFixture : public DecompositionFixture<MESHSPEC, LSPerPhasePolicy, 3>
{
public:
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::write_mesh;
  using DecompositionFixture<MESHSPEC, LSPerPhasePolicy, 3>::cdmesh;

protected:
  void set_levelsets_and_phases()
  {
    this->setup_ls_fields();
    mySurfaceIdentifiers.clear();
    for (auto & ls : this->levelset_fields())
      mySurfaceIdentifiers.push_back(ls.identifier);
  }
  void decompose_with_element_phases(const std::vector<int> & elemPhases)
  {
    PhasePerElementInterfaceGeometry interfaceGeometry(mySurfaceIdentifiers);
    interfaceGeometry.set_element_phases(mMesh, mBuilder.get_assigned_element_global_ids(), elemPhases);
    cdmesh->decompose_mesh(interfaceGeometry, -1);
  }
  void cleanup_small_clusters_for_levelsets(const Surface_Identifier fromNegLS, const Surface_Identifier toNegLS, const size_t smallClusterSize)
  {
    NamedPhase fromPhase(std::to_string(fromNegLS.get()));
    fromPhase.tag().add(fromNegLS, -1);
    NamedPhase toPhase(std::to_string(toNegLS.get()));
    toPhase.tag().add(toNegLS, -1);
    const Phase_Support & phaseSupport = Phase_Support::get(mMesh.mesh_meta_data());
    cleanup_small_clusters_of_phase(mMesh, phaseSupport, fromPhase, toPhase, smallClusterSize);
  }

  void test_elements_in_phases(const std::vector<int> & elemPhases)
  {
    const Phase_Support & phaseSupport = Phase_Support::get(mMesh.mesh_meta_data());
    const auto & elemIds = mBuilder.get_assigned_element_global_ids();
    STK_ThrowRequire(elemIds.size() == elemPhases.size());
    for (size_t n=0; n<mBuilder.get_assigned_element_global_ids().size(); ++n)
    {
      stk::mesh::Entity elem = mMesh.get_entity(stk::topology::ELEMENT_RANK, elemIds[n]);
      if (mMesh.is_valid(elem))
      {
        const PhaseTag elemPhase =  determine_phase_for_entity(mMesh, elem, phaseSupport);
        EXPECT_TRUE(elemPhase.contain(mySurfaceIdentifiers[elemPhases[n]], -1));
      }
    }
  }

  void test_phases_after_decomposition_and_hinge_removal(const std::vector<int> & initialElemPhases, const std::vector<int> & goldElemPhases)
  {
    this->reset_cdmesh(); // new interface conforming mesh each time
    decompose_with_element_phases(initialElemPhases);

    const size_t smallClusterSize = 3;
    cleanup_small_clusters_for_levelsets(mySurfaceIdentifiers[1], mySurfaceIdentifiers[0], smallClusterSize);
    cleanup_small_clusters_for_levelsets(mySurfaceIdentifiers[2], mySurfaceIdentifiers[0], smallClusterSize);

    test_elements_in_phases(goldElemPhases);
  }

  std::vector<Surface_Identifier> mySurfaceIdentifiers;
};

class SmallClusterCleanupPatchOfRegularTrisAroundNode : public SmallClusterCleanupFixture<PatchOfRegularTrisAroundNode>
{
public:
  SmallClusterCleanupPatchOfRegularTrisAroundNode()
  {
    set_valid_proc_sizes_for_test({1,2,3,4});
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1}, {0,0,0,0,0,0});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1}, {0,0,1,1,0,0});
    else if(stk::parallel_machine_size(mComm) == 3)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1}, {0,1,2,0,1,2});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1}, {0,1,2,3,0,1});
    set_levelsets_and_phases();
  }
};

TEST_F(SmallClusterCleanupPatchOfRegularTrisAroundNode, allElementsInPhase0_smallClusterCleanup_doesNothing)
{
  if(is_valid_proc_size_for_test())
  {
    test_phases_after_decomposition_and_hinge_removal({0,0,0,0,0,0}, {0,0,0,0,0,0});
  }
}

TEST_F(SmallClusterCleanupPatchOfRegularTrisAroundNode, hingeElementsInPhase1_smallClusterCleanup_hingesAreRemoved)
{
  if(is_valid_proc_size_for_test())
  {
    test_phases_after_decomposition_and_hinge_removal({1,0,0,0,0,0}, {0,0,0,0,0,0});
    test_phases_after_decomposition_and_hinge_removal({1,1,0,0,0,0}, {0,0,0,0,0,0});
    test_phases_after_decomposition_and_hinge_removal({1,0,1,0,1,0}, {0,0,0,0,0,0});

    test_phases_after_decomposition_and_hinge_removal({1,1,1,0,0,0}, {1,1,1,0,0,0});
    test_phases_after_decomposition_and_hinge_removal({1,1,1,0,1,0}, {1,1,1,0,0,0});
  }
}

TEST_F(SmallClusterCleanupPatchOfRegularTrisAroundNode, hingeElementsInPhase1and2_smallClusterCleanup_hingesAreRemoved)
{
  if(is_valid_proc_size_for_test())
  {
    test_phases_after_decomposition_and_hinge_removal({1,2,0,0,0,0}, {0,0,0,0,0,0});
    test_phases_after_decomposition_and_hinge_removal({1,2,2,0,0,0}, {0,0,0,0,0,0});
    test_phases_after_decomposition_and_hinge_removal({1,2,1,2,1,2}, {0,0,0,0,0,0});
  }
}

}


