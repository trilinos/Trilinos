/*
 * Akri_Unit_Refine_Tet.cpp
 *
 *  Created on: Dec 15, 2022
 *      Author: drnoble
 */
#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <Akri_OutputUtils.hpp>

namespace krino {

class RegularTetRefinement : public RefinementFixture<RegularTet>
{
public:
  RegularTetRefinement()
  {;
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element()
  {
    const std::vector<stk::mesh::Entity> ownedElements = get_owned_elements();
    return ownedElements[0];
  }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }
protected:
};


class UMRRegularTetRefinement : public RefinementFixture<UMRRegularTet>
{
public:
  UMRRegularTetRefinement()
  {;
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,1,0,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 8)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,4,5,6,7});
  }
protected:
};

class FourRightTetsRefinement : public RefinementFixture<FourRightTets>
{
public:
  FourRightTetsRefinement()
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,0,1,1});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,2});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});
    }
  }
protected:
};

class RightTetSurroundedByEdgeTetsRefinement : public RefinementFixture<RightTetSurroundedByEdgeTets>
{
public:
  RightTetSurroundedByEdgeTetsRefinement()
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,0,0,0,0,0,0});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,0,1,0,1,0});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,2,0,1,2,0});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,2,3,0,1,2});
    }
  }
protected:
};

class RightTetSurroundedByFaceTetsRefinement : public RefinementFixture<RightTetSurroundedByFaceTets>
{
public:
  RightTetSurroundedByFaceTetsRefinement()
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,0,0,0,0});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,0,1,0});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,2,0,1});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,2,3,0});
    }
  }
protected:
};

TEST_F(RegularTetRefinement, givenMeshWithSingleTetMarked_afterRefinement_have8ChildElements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    mark_elements_for_refinement({get_element()});
    do_refinement();

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

TEST_F(RegularTetRefinement, meshAfter3LevelsOfUMR_have585Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement(3);

    EXPECT_EQ(585u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RegularTetRefinement, perceptRefinement_meshAfter3LevelsOfUMR_have585Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    setup_percept_refinement();
    perform_iterations_of_percept_uniform_refinement(3);

    EXPECT_EQ(585u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(FourRightTetsRefinement, meshAfter2LevelsOfUMR_haveCorrectNumberOfElementsAndNodes)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    perform_iterations_of_uniform_refinement(2);

    const unsigned goldNumElems = 4 + 32 + 256;
    EXPECT_EQ(goldNumElems, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
    const unsigned goldNumNodes = 6 + 13 + 66;
    EXPECT_EQ(goldNumNodes, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
  }
}

TEST_F(FourRightTetsRefinement, meshAfter3LevelsOfUMR_have2340Elements)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    perform_iterations_of_uniform_refinement(3);

    EXPECT_EQ(2340u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
    const unsigned goldNumNodes = 6 + 13 + 66 + 404;
    EXPECT_EQ(goldNumNodes, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
  }
}

TEST_F(FourRightTetsRefinement, perceptRefinment_meshAfter3LevelsOfUMR_have2340Elements)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    setup_percept_refinement();
    perform_iterations_of_percept_uniform_refinement(3);

    EXPECT_EQ(2340u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
    const unsigned goldNumNodes = 6 + 13 + 66 + 404;
    EXPECT_EQ(goldNumNodes, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
  }
}

TEST_F(RegularTetRefinement, meshAfter2LevelsOfUMR_haveExpectedQuality)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement(2);

    const double quality = compute_mesh_quality();
    const double goldQuality = 0.7;
    EXPECT_GT(quality, goldQuality);
  }
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, checkAllPossibleRefinementsInParallel_getSameQualityAsProducedByPercept)
{
  // These gold values were generated by running Percept in 10/2022
  const std::array<unsigned,64> goldNumElementsByCase{7,27,17,33,17,33,26,39,17,33,26,39,26,39,35,45,17,33,26,39,27,39,36,45,26,39,35,45,36,45,45,51,17,33,26,39,26,39,35,45,27,39,36,45,36,45,45,51,26,39,35,45,36,45,45,51,36,45,45,51,45,51,54,57};
  const std::array<unsigned,64> goldNumNodesByCase{16,22,22,27,22,27,28,32,22,27,28,32,28,32,34,37,22,27,28,32,28,32,34,37,28,32,34,37,34,37,40,42,22,27,28,32,28,32,34,37,28,32,34,37,34,37,40,42,28,32,34,37,34,37,40,42,34,37,40,42,40,42,46,47};
  const std::array<double,64> goldQualityByCase{
    0.71,0.32,0.32,0.32,0.41,0.32,0.18,0.32,0.32,0.32,0.20,0.32,0.18,0.32,0.18,0.32,
    0.32,0.32,0.20,0.32,0.32,0.32,0.18,0.32,0.20,0.32,0.20,0.41,0.18,0.32,0.18,0.41,
    0.41,0.32,0.18,0.32,0.41,0.32,0.41,0.32,0.32,0.32,0.18,0.32,0.18,0.32,0.18,0.32,
    0.18,0.32,0.18,0.32,0.24,0.32,0.24,0.32,0.26,0.32,0.18,0.41,0.18,0.32,0.18,0.41};

    const unsigned numEdges = 6;
    for (int iCaseId=0; iCaseId<(1<<numEdges); ++iCaseId)
  {
    std::vector<unsigned> edgeElementsToRefine;
    for (unsigned iEdge=0; iEdge<numEdges; ++iEdge)
      if (iCaseId & (1<<iEdge))
        edgeElementsToRefine.push_back(iEdge);

    test_refinement_of_given_elements(false, edgeElementsToRefine, goldNumElementsByCase[iCaseId], goldNumNodesByCase[iCaseId], goldQualityByCase[iCaseId]);

    unrefine_mesh(false);
  }
}

TEST_F(RightTetSurroundedByFaceTetsRefinement, checkAllPossibleRefinementsInParallel_expectBoundarySidesAreCorrect)
{
  const unsigned numFaces = 4;
  const unsigned numCases = 1<<numFaces;
  for (unsigned i=0; i<numCases; ++i)
  {
    std::vector<unsigned> faceElementsToRefine;
    for (unsigned iFace=0; iFace<4; ++iFace)
      if (i & (1<<iFace))
        faceElementsToRefine.push_back(1+iFace);

    refine_elements_with_given_indices(false, faceElementsToRefine);

    stk::mesh::create_all_sides(mMesh, mMesh.mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(mBuilder.check_boundary_sides());
    EXPECT_TRUE(mBuilder.check_block_boundary_sides());

    unrefine_mesh(false);
  }
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, refineCenterElemAndThenChildOfCenterElem_noHangingNodes)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[0]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(false, elementsToRefine);

    elementsToRefine.clear();
    if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(false, centerElem);
      ASSERT_EQ(8u, childElems.size());
      elementsToRefine.assign({childElems[7]});
    }

    refine_elements(false, elementsToRefine);

    EXPECT_FALSE(myRefinement.have_any_hanging_refined_nodes());
  }
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, afterEachOfSixRoundsOfRefinementOfEdgeElements_centerElementHasCorrectNumberOfChildren)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[0]);

    const std::array<unsigned,6> elementsToRefineByRefinementLevel = {1,2,3,4,5,6};
    const std::array<unsigned,6> goldNumCenterElemChildrenByRefinementLevel = {2,3,4,6,7,8};

    for (unsigned i=0; i<6; ++i)
    {
      refine_elements_with_given_indices(false, {elementsToRefineByRefinementLevel[i]});

      if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
      {
        EXPECT_EQ(goldNumCenterElemChildrenByRefinementLevel[i], myRefinement.get_num_children(centerElem));
      }
    }
  }
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, markedAnyTransitionElementForEveryEdgeConfiguration_parentElementGetsFullyRefined)
{
  const bool usePercept = false;
  const int indexOfCentralElement = 0;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(usePercept, indexOfCentralElement);
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, DISABLED_BECAUSE_SLOW_percept_markedAnyTransitionElementForEveryEdgeConfiguration_parentElementGetsFullyRefined)
{
  const bool usePercept = true;
  const int indexOfCentralElement = 0;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(usePercept, indexOfCentralElement);
}

TEST_F(UMRRegularTetRefinement, performanceRefinementThenUnrefinementTest)
{
  // As of 12/13/2022:
  // Krino times:
  // NP=1 Time:    21718 ms
  // NP=2 Time:    11991 ms
  // NP=4 Time:     7063 ms
  // NP=8 Time:     6447 ms

  // Percept times:
  // NP=1 Percept: 91221 ms
  // NP=2 Percept: 51785 ms
  // NP=4 Percept: 32157 ms
  // NP=8 Percept: 28709 ms


  const std::vector<int> testActiveProcs{1,2,4,8};
  const int parallelSize = stk::parallel_machine_size(mComm);
  if (std::find(testActiveProcs.begin(), testActiveProcs.end(), parallelSize) == testActiveProcs.end())
    return;

  const bool doWriteMesh = false;
  const bool usePercept = false;

  if (doWriteMesh)
    write_mesh("test.e");

  int count = 0;

  const std::vector<size_t> goldNumElementsByRefinementLevel = {72, 512, 2892, 13908, 61146, 255890, 1045094};
  for (size_t i=0; i<goldNumElementsByRefinementLevel.size(); ++i)
  {
    mark_elements_spanning_x_equal_0();

    if (doWriteMesh)
      refine_marked_elements(usePercept, create_file_name("test", ++count));
    else
      refine_marked_elements(usePercept);

    const size_t numElements = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
    EXPECT_EQ(goldNumElementsByRefinementLevel[i], numElements);

    if (0 == stk::parallel_machine_rank(mComm))
      std::cout << "Refinement Level " << i+1 << ", num elements = " << numElements << std::endl;
  }

  const std::vector<size_t> goldNumElementsByUnrefinementLevel = {287808, 76440, 18164, 3648, 576, 72, 8};
  for (size_t i=0; i<goldNumElementsByUnrefinementLevel.size(); ++i)
  {
    mark_all_elements_for_unrefinement();

    if (doWriteMesh)
      refine_marked_elements(usePercept, create_file_name("test", ++count));
    else
      refine_marked_elements(usePercept);

    const size_t numElements = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
    EXPECT_EQ(goldNumElementsByUnrefinementLevel[i], numElements);

    if (0 == stk::parallel_machine_rank(mComm))
      std::cout << "Unrefinement Level " << i+1 << ", num elements = " << numElements << std::endl;
  }
}

TEST_F(RegularTetRefinement, meshRefinedTwiceWithParentAndChildrenMovedToNewProc_unrefinementCausesParentToReturnToOriginatingProc)
{
  if(stk::parallel_machine_size(mComm) > 1)
  {
    perform_iterations_of_uniform_refinement(2);

    move_owned_elements_with_given_ids_and_owned_attached_entities_to_processor({1004, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057}, 1);

    EXPECT_TRUE(mBuilder.check_boundary_sides());
    EXPECT_TRUE(check_face_and_edge_ownership(mMesh));

    mark_all_elements_for_unrefinement();
    refine_marked_elements();

    EXPECT_TRUE(mBuilder.check_boundary_sides());
    EXPECT_TRUE(check_face_and_edge_ownership(mMesh));

    stk::mesh::Entity movedParentElement = mMesh.get_entity(stk::topology::ELEMENT_RANK, 1004);
    if (mMesh.is_valid(movedParentElement))
    {
      EXPECT_EQ(0, mMesh.parallel_owner_rank(movedParentElement));
    }

    mark_all_elements_for_unrefinement();
    ASSERT_NO_THROW(refine_marked_elements());

    EXPECT_EQ(1u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    EXPECT_TRUE(mBuilder.check_boundary_sides());
    EXPECT_TRUE(check_face_and_edge_ownership(mMesh));
  }
}


}


