#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <Akri_OutputUtils.hpp>


namespace krino {

class RegularTriRefinement : public RefinementFixture<RegularTri>
{
public:
  RegularTriRefinement()
  {
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element()
  {
    const std::vector<stk::mesh::Entity> ownedElements = get_owned_elements();
    return ownedElements[0];
  }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }
protected:
};

class UMRRegularTriRefinement : public RefinementFixture<UMRRegularTri>
{
public:
  UMRRegularTriRefinement()
  {;
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,0,1});
    else if(stk::parallel_machine_size(mComm) == 3)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,0});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});
  }
protected:
};

class RightTriSurroundedByEdgeTrisRefinement : public RefinementFixture<RightTriSurroundedByEdgeTris>
{
public:
  RightTriSurroundedByEdgeTrisRefinement()
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {2,2,2,1}, {0,0,0,0});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {2,2,2,1}, {0,0,1,1});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {2,2,2,1}, {0,1,2,2});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {2,2,2,1}, {0,1,2,3});
    }
  }
protected:
};

TEST_F(RegularTriRefinement, givenMeshWithoutAnyRefinement_whenQueryingParentsAndChildren_noParentOrChildren)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const stk::mesh::Entity elem = get_element();
    EXPECT_FALSE(myRefinement.is_parent(elem));
    EXPECT_FALSE(myRefinement.is_child(elem));

    EXPECT_EQ(0u, myRefinement.get_num_children(elem));

    stk::mesh::Entity parent = myRefinement.get_parent(elem);
    EXPECT_FALSE(mMesh.is_valid(parent));
  }
}

TEST_F(RegularTriRefinement, givenMeshWithNoElementMarked_whenFindingEdgesToRefine_noEdgesToRefine)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    clear_refinement_marker();

    const TransitionElementEdgeMarker edgeMarker(mMesh, myRefinement, myElementMarkerField.name());
    myRefinement.find_edges_to_refine(edgeMarker);

    EXPECT_EQ(0u, myRefinement.get_num_edges_to_refine());
  }
}

TEST_F(RegularTriRefinement, givenMeshWithSingleTriMarked_whenFindingEdgesToRefine_3EdgesToRefine)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    clear_refinement_marker();
    mark_elements_for_refinement({get_element()});

    const TransitionElementEdgeMarker edgeMarker(mMesh, myRefinement, myElementMarkerField.name());
    myRefinement.find_edges_to_refine(edgeMarker);

    EXPECT_EQ(3u, myRefinement.get_num_edges_to_refine());
  }
}

TEST_F(RegularTriRefinement, givenMeshWithSingleTriMarked_afterRefinement_all3EdgesHaveRefineNodes)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    mark_elements_for_refinement({get_element()});
    do_refinement();
    EXPECT_EQ(6u, get_global_num_entities(mMesh, stk::topology::NODE_RANK));

    for (int iEdge=0; iEdge<3; ++iEdge)
    {
      const Edge edge = get_edge(iEdge);
      const stk::mesh::Entity childEdgeNode = myRefinement.get_edge_child_node(edge);
      EXPECT_TRUE(mMesh.is_valid(childEdgeNode));

      const auto & edgeNodes = get_edge_nodes(edge);
      const stk::math::Vector3d edgeNode0Coords = get_node_coordinates(edgeNodes[0]);
      const stk::math::Vector3d edgeNode1Coords = get_node_coordinates(edgeNodes[1]);
      const stk::math::Vector3d childNodeCoords = get_node_coordinates(childEdgeNode);
      const double error = (0.5*(edgeNode0Coords + edgeNode1Coords) - childNodeCoords).length();
      EXPECT_NEAR(0., error, 1.e-6);
    }
  }
}

TEST_F(RegularTriRefinement, givenMeshWithSingleTriMarked_afterRefinement_have4ChildElements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    mark_elements_for_refinement({get_element()});
    do_refinement();

    EXPECT_EQ(5u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    EXPECT_TRUE(myRefinement.is_parent(get_element()));
    EXPECT_EQ(4u, myRefinement.get_num_children(get_element()));

    const std::vector<stk::mesh::Entity> children = myRefinement.get_children(get_element());
    for (auto && child : children)
    {
      EXPECT_TRUE(myRefinement.is_child(child));
      EXPECT_EQ(get_element(), myRefinement.get_parent(child));
    }
  }
}

TEST_F(RegularTriRefinement, twoRoundsOfMarkingSameElement_secondRoundDoesNothing)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    refine_elements_with_given_indices(false, {0});

    EXPECT_EQ(5u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    refine_elements_with_given_indices(false, {0}); // marking parent element, should do nothing

    EXPECT_EQ(5u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RegularTriRefinement, meshAfter3LevelsOfUMR_have85Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    perform_iterations_of_uniform_refinement(3);

    EXPECT_EQ(85u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RegularTriRefinement, perceptRefinement_meshAfter3LevelsOfUMR_have85Elements)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    setup_percept_refinement();
    perform_iterations_of_percept_uniform_refinement(3);

    EXPECT_EQ(85u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}


TEST_F(RightTriSurroundedByEdgeTrisRefinement, refinementOfOneTriInParallel_expectEdgeIsRefinedAndCoordinatesAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    const unsigned indexOfElemToRefine = 0;
    const stk::mesh::Entity edgeNode0 = mMesh.get_entity(stk::topology::NODE_RANK, mBuilder.get_assigned_node_global_ids()[3]);
    const stk::mesh::Entity edgeNode1 = mMesh.get_entity(stk::topology::NODE_RANK, mBuilder.get_assigned_node_global_ids()[5]);

    refine_elements_with_given_indices(false, {indexOfElemToRefine});

    EXPECT_EQ(10u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
    EXPECT_EQ(9u, get_global_num_entities(mMesh, stk::topology::NODE_RANK));

    const bool isEdgeNodeOwnedOrShared0 = mMesh.is_valid(edgeNode0) && (mMesh.bucket(edgeNode0).owned() || mMesh.bucket(edgeNode0).shared());
    const bool isEdgeNodeOwnedOrShared1 = mMesh.is_valid(edgeNode1) && (mMesh.bucket(edgeNode1).owned() || mMesh.bucket(edgeNode1).shared());

    if (isEdgeNodeOwnedOrShared0 && isEdgeNodeOwnedOrShared1)
    {
      const Edge edge = edge_from_edge_nodes(mMesh, edgeNode0, edgeNode1);

      const stk::mesh::Entity childEdgeNode = myRefinement.get_edge_child_node(edge);
      ASSERT_TRUE(mMesh.is_valid(childEdgeNode));

      const stk::math::Vector3d edgeNode0Coords = get_node_coordinates(edgeNode0);
      const stk::math::Vector3d edgeNode1Coords = get_node_coordinates(edgeNode1);
      const stk::math::Vector3d childNodeCoords = get_node_coordinates(childEdgeNode);
      const double error = (0.5*(edgeNode0Coords + edgeNode1Coords) - childNodeCoords).length();
      EXPECT_NEAR(0., error, 1.e-6);
    }
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, checkAllPossibleRefinementsInParallel_expectBoundarySidesAreCorrect)
{
  for (int i=0; i<8; ++i)
  {
    std::vector<unsigned> edgeElementsToRefine;
    for (unsigned iEdge=0; iEdge<3; ++iEdge)
      if (i & (1<<iEdge))
        edgeElementsToRefine.push_back(iEdge);

    refine_elements_with_given_indices(false, edgeElementsToRefine);

    stk::mesh::create_all_sides(mMesh, mMesh.mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(mBuilder.check_boundary_sides());
    EXPECT_TRUE(mBuilder.check_block_boundary_sides());

    unrefine_mesh(false);
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, checkAllPossibleRefinementsInParallel_getSameQualityAsProducedByPercept)
{
  // These gold values were generated by running Percept in 10/2022
  const std::array<unsigned,8> goldNumElementsByCase{4,10,10,15,10,15,15,20};
  const std::array<unsigned,8> goldNumNodesByCase{6,9,9,12,9,12,12,15};
  const std::array<double,8> goldQualityByCase{0.82,0.36,0.36,0.36,0.82,0.82,0.82,0.82};

  for (int i=0; i<8; ++i)
  {
    std::vector<unsigned> edgeElementsToRefine;
    for (unsigned iEdge=0; iEdge<3; ++iEdge)
      if (i & (1<<iEdge))
        edgeElementsToRefine.push_back(iEdge);

    test_refinement_of_given_elements(false, edgeElementsToRefine, goldNumElementsByCase[i], goldNumNodesByCase[i], goldQualityByCase[i]);

    unrefine_mesh(false);
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, afterEachOfThreeRoundsOfRefinementOfEdgeElements_centerElementHasCorrectNumberOfChildren)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[3]);

    const std::array<unsigned,3> elementsToRefineByRefinementLevel = {0,1,2};
    const std::array<unsigned,3> goldNumCenterElemChildrenByRefinementLevel = {2,3,4};

    for (unsigned i=0; i<3; ++i)
    {
      refine_elements_with_given_indices(false, {elementsToRefineByRefinementLevel[i]});

      if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
      {
        EXPECT_EQ(goldNumCenterElemChildrenByRefinementLevel[i], myRefinement.get_num_children(centerElem));
      }
    }
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, refineCenterElemAndThenChildOfCenterElem_noHangingNodes)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[3]);

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
      ASSERT_EQ(4u, childElems.size());
      elementsToRefine.assign({childElems[0]});
    }

    refine_elements(false, elementsToRefine);

    EXPECT_FALSE(myRefinement.have_any_hanging_refined_nodes());

    const unsigned goldNumElems = 4 + 4+3*2 + 4+2+2*(6-2);
    EXPECT_EQ(goldNumElems, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, refineCenterElemAndThenMarkTransitionElement_parentOfTransitionElementGetsRefined)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[3]);
    stk::mesh::Entity edgeElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[0]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(false, elementsToRefine);

    elementsToRefine.clear();
    if (mMesh.is_valid(edgeElem) && mMesh.bucket(edgeElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(false, edgeElem);
      ASSERT_EQ(2u, childElems.size());
      elementsToRefine.assign({childElems[0]});
    }

    refine_elements(false, elementsToRefine);

    EXPECT_FALSE(myRefinement.have_any_hanging_refined_nodes());

    if (mMesh.is_valid(edgeElem) && mMesh.bucket(edgeElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(false, edgeElem);
      EXPECT_EQ(4u, childElems.size());
      for (auto && childElem : childElems)
      {
        EXPECT_EQ(0u, get_num_children(false, childElem));
      }
    }
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, refineCenterElemAndThenMarkAllCenterElementChildrenForUnrefinement_perceptAfterUnrefineOriginalMeshObtained)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[3]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(true, elementsToRefine);

    const unsigned goldNumElemsAfterRefinement = 4 + 4 + 6;
    EXPECT_EQ(goldNumElemsAfterRefinement, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));

    clear_refinement_marker();
    if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(true, centerElem);
      EXPECT_EQ(4u, childElems.size());
      mark_elements_for_unrefinement(childElems);
    }

    refine_marked_elements(true);

    const unsigned goldNumElemsAfterUnrefinement = 4;
    EXPECT_EQ(goldNumElemsAfterUnrefinement, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, markedAnyTransitionElementForEveryEdgeConfiguration_parentElementGetsFullyRefined)
{
  const bool usePercept = false;
  const int indexOfCentralElement = 3;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(usePercept, indexOfCentralElement);
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, percept_markedAnyTransitionElementForEveryEdgeConfiguration_parentElementGetsFullyRefined)
{
  const bool usePercept = true;
  const int indexOfCentralElement = 3;
  test_refinement_of_transition_element_leads_to_refinement_of_parent(usePercept, indexOfCentralElement);
}

TEST_F(UMRRegularTriRefinement, refinementThenUnrefinementTest)
{
  if(stk::parallel_machine_size(mComm) > 4)
    return;

  const bool doWriteMesh = false;
  const bool usePercept = false;

  if (doWriteMesh)
    write_mesh("test.e");

  int count = 0;

  const std::vector<size_t> goldNumElementsByRefinementLevel = {20, 68, 180, 420, 916, 1924, 3956};
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

  const std::vector<size_t> goldNumElementsByUnrefinementLevel = {1924, 916, 420, 180, 68, 20, 4};
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

TEST_F(RegularTriRefinement, meshRefinedTwiceWithParentAndChildrenMovedToNewProc_unrefinementCausesParentToReturnToOriginatingProc)
{
  if(stk::parallel_machine_size(mComm) > 1)
  {
    perform_iterations_of_uniform_refinement(2);

    move_owned_elements_with_given_ids_and_owned_attached_entities_to_processor({1004, 1010, 1011, 1012, 1013}, 1);

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

