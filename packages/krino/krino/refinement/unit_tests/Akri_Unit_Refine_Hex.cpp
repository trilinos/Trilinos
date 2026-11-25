#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <Akri_Unit_RefinementFixture_Edges.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace krino {

class RegularHexRefinement : public RefinementFixture<RegularHex>
{
public:
  RegularHexRefinement()
  {
    StkMeshHexFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
protected:
};

TEST_F(RegularHexRefinement, meshAfter3LevelsOfUMRViaUniformMarker_have585Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement_with_uniform_marker(3);

    EXPECT_EQ(729u, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
    EXPECT_EQ(585u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    // Check that boundary sides were correctly propagated
    stk::mesh::create_all_sides(mMesh, mMesh.mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(mBuilder.check_boundary_sides());
  }
}

class UMRRegularHexRefinement : public RefinementFixture<UMRRegularHex>
{
public:
  UMRRegularHexRefinement()
  {;
    set_valid_proc_sizes_for_test({1,2,4,8});
    this->build_mesh_from_spec(meshSpec, {{2, {0,1,1,0,0,1,1,0}}});
  }
protected:
};

TEST_F(UMRRegularHexRefinement, meshAfterEachLevelOfRefineHasCorrectNumElemsAndNodesInParallel)
{
  const std::vector<int> testActiveProcs{1,2,4,8};
  const int parallelSize = stk::parallel_machine_size(mComm);
  if (std::find(testActiveProcs.begin(), testActiveProcs.end(), parallelSize) == testActiveProcs.end())
    return;

  const std::vector<size_t> goldNumElementsByRefinementLevel = {72, 584, 4680};
  const std::vector<size_t> goldNumNodesByRefinementLevel = {125, 729, 4913};
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

class UMRRegularHexRefinementWithEdges : public WithExplicitEdges<UMRRegularHexRefinement>
{
public:
  UMRRegularHexRefinementWithEdges()
  {
    if (is_valid_proc_size_for_test()) create_edges();
  }
};

TEST_F(UMRRegularHexRefinementWithEdges,
    uniformRefinementThenUnrefinement_edgesAttachedToAllActiveElements)
{
  if (!is_valid_proc_size_for_test()) return;

  const size_t numOriginalEdges = num_edges_in(mMesh.mesh_meta_data().universal_part());
  EXPECT_EQ(num_original_edges(), numOriginalEdges);

  // Each original edge is split in two at each level, and the edges are checked after each level.
  perform_iterations_of_uniform_refinement_with_uniform_marker(1);
  EXPECT_EQ(numOriginalEdges * (1 + 2), num_edges_in(mMesh.mesh_meta_data().universal_part()));
  perform_iterations_of_uniform_refinement_with_uniform_marker(1);
  EXPECT_EQ(numOriginalEdges * (1 + 2 + 4), num_edges_in(mMesh.mesh_meta_data().universal_part()));

  // Marker unrefinement does not support hex face and centroid nodes, so the mesh is fully
  // unrefined.
  unrefine_mesh();
  EXPECT_EQ(numOriginalEdges, num_edges_in(mMesh.mesh_meta_data().universal_part()));
}
}
