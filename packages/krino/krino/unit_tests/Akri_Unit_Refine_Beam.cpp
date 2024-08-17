#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>

namespace krino {

class RegularBeamRefinement : public RefinementFixture<RegularBeam>
{
public:
  RegularBeamRefinement()
  {
    StkMeshBeamFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
protected:
};

TEST_F(RegularBeamRefinement, meshAfter3LevelsOfUMRViaUniformMarker_have15Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement_with_uniform_marker(3);

    EXPECT_EQ(9u, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
    EXPECT_EQ(15u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

class UMRRegularBeamRefinement : public RefinementFixture<UMRRegularBeam>
{
public:
  UMRRegularBeamRefinement()
  {;
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});
  }
protected:
};

TEST_F(UMRRegularBeamRefinement, meshAfterEachLevelOfRefineHasCorrectNumElemsAndNodesInParallel)
{
  const std::vector<int> testActiveProcs{1,2};
  const int parallelSize = stk::parallel_machine_size(mComm);
  if (std::find(testActiveProcs.begin(), testActiveProcs.end(), parallelSize) == testActiveProcs.end())
    return;

  const std::vector<size_t> goldNumElementsByRefinementLevel = {6, 14, 30, 62};
  const std::vector<size_t> goldNumNodesByRefinementLevel = {5, 9, 17, 33};
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
