#include <Akri_Unit_RefinementFixture_Tet.hpp>
#include <Akri_Edge.hpp>
#include <Akri_OutputUtils.hpp>
#include <stk_util/diag/PrintTimer.hpp>

namespace krino {

TEST_F(RegularTetRefinement, singleTetWithAllEdgesMarked_afterEdgeRefinement_have8ChildElements)
{
  if(is_valid_proc_size_for_test())
  {
    std::vector<Edge> edgesToRefine;
    fill_entity_edges(mMesh, get_element(), edgesToRefine);
    mark_edges_for_refinement(edgesToRefine);
    do_edge_refinement();

    EXPECT_EQ(9u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    EXPECT_TRUE(myRefinement.is_parent(get_element()));
    EXPECT_EQ(8u, myRefinement.get_num_children(get_element()));

    const std::vector<stk::mesh::Entity> children = myRefinement.get_children(get_element());
    for (auto && child : children)
    {
      EXPECT_TRUE(myRefinement.is_child(child));
      EXPECT_EQ(get_element(), myRefinement.get_parent(child));
    }
  }
}

TEST_F(RegularTetRefinement, singleTetWith1EdgeMarked_afterEdgeRefinement_have2ChildElements)
{
  if(is_valid_proc_size_for_test())
  {
    std::vector<Edge> elemEdges;
    fill_entity_edges(mMesh, get_element(), elemEdges);
    mark_edges_for_refinement({elemEdges[0]});
    do_edge_refinement();

    EXPECT_EQ(3u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    EXPECT_TRUE(myRefinement.is_parent(get_element()));
    EXPECT_EQ(2u, myRefinement.get_num_children(get_element()));
  }
}

template<typename FIXTURE>
void test_edge_refinement(FIXTURE & fixture, const std::vector<size_t> & goldNumElementsByRefinementLevel)
{
  if(!fixture.is_valid_proc_size_for_test())
    return;

  const bool doWriteMesh = false;

  if (doWriteMesh)
    fixture.write_mesh("test.e");

  int count = 0;


  for (size_t i=0; i<goldNumElementsByRefinementLevel.size(); ++i)
  {
    fixture.mark_edges_crossing_x_value(0.01);

    if (doWriteMesh)
      fixture.refine_marked_edges(create_file_name("test", ++count));
    else
      fixture.refine_marked_edges();

    const size_t numElements = get_global_num_entities(fixture.mMesh, stk::topology::ELEMENT_RANK);
    EXPECT_EQ(goldNumElementsByRefinementLevel[i], numElements);

    if (0 == stk::parallel_machine_rank(fixture.mComm))
      std::cout << "Refinement Level " << i+1 << ", num elements = " << numElements << std::endl;
  }
}

TEST_F(UMRRegularTetRefinement, performanceEdgeRefinementTest)
{
  // As of 9/5/2024:
  // Krino times:
  // NP=1 Time:    682 ms

  const std::vector<size_t> goldNumElementsByRefinementLevel = {40, 140, 364, 760, 1376, 3392, 8318, 21976, 57398};
  test_edge_refinement(*this, goldNumElementsByRefinementLevel);

  stk::diag::printTimersTable(std::cout, sierra::Diag::sierraTimer(), stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, mComm);
}

TEST_F(UMRRegularTet10Refinement, performanceEdgeRefinementTest)
{
  // goldNumElementsByRefinementLevel differ from Tet4 case because quality metric is different
  const std::vector<size_t> goldNumElementsByRefinementLevel = {40, 140, 364, 760, 1376, 3384, 8294, 21924, 57282};
  test_edge_refinement(*this, goldNumElementsByRefinementLevel);
}

}
