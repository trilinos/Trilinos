#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace krino {

class RegularQuadRefinement : public RefinementFixture<RegularQuad>
{
public:
  RegularQuadRefinement()
  {
    StkMeshQuadFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
protected:
};

TEST_F(RegularQuadRefinement, meshAfter3LevelsOfUMRViaUniformMarker_have85Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement_with_uniform_marker(3);

    EXPECT_EQ(81u, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
    EXPECT_EQ(85u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    // Check that boundary sides were correctly propagated
    stk::mesh::create_all_sides(mMesh, mMesh.mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(mBuilder.check_boundary_sides());
  }
}

class UMRRegularQuadRefinement : public RefinementFixture<UMRRegularQuad>
{
public:
  UMRRegularQuadRefinement()
  {;
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,1,0});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});
  }
protected:
};

TEST_F(UMRRegularQuadRefinement, meshAfterEachLevelOfRefineHasCorrectNumElemsAndNodesInParallel)
{
  const std::vector<int> testActiveProcs{1,2,4};
  const int parallelSize = stk::parallel_machine_size(mComm);
  if (std::find(testActiveProcs.begin(), testActiveProcs.end(), parallelSize) == testActiveProcs.end())
    return;

  const std::vector<size_t> goldNumElementsByRefinementLevel = {20, 84, 340, 1364};
  const std::vector<size_t> goldNumNodesByRefinementLevel = {25, 81, 289, 1089};
  for (size_t i=0; i<goldNumElementsByRefinementLevel.size(); ++i)
  {
    perform_iterations_of_uniform_refinement_with_uniform_marker(1);

    const size_t numElems = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
    const size_t numNodes = get_global_num_entities(mMesh, stk::topology::NODE_RANK);
    EXPECT_EQ(goldNumElementsByRefinementLevel[i], numElems);
    EXPECT_EQ(goldNumNodesByRefinementLevel[i], numNodes);

    if (0 == stk::parallel_machine_rank(mComm))
    {
      std::cout << "Refinement Level " << i+1 << ", num elements = " << numElems << ", num nodes = " << numNodes << std::endl;
    }
  }
}

}
