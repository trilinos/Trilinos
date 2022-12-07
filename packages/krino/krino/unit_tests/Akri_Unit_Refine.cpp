#include <stk_mesh/base/MetaData.hpp>
#include <Akri_Refinement.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Edge.hpp>
#include <Akri_AuxMetaData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_AdaptivityInterface.hpp>
#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Quality.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MOAB_TetRefiner.hpp>
#include <Akri_TriRefiner.hpp>
#include <Akri_TransitionElementEdgeMarker.hpp>
#include <Akri_RefinementInterface.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

namespace krino {

void set_refinement_marker_field(FieldRef elementMarkerField, const int value)
{
  stk::mesh::field_fill(value, elementMarkerField, stk::mesh::selectField(elementMarkerField));
}

void clear_refinement_marker_field(FieldRef elementMarkerField)
{
  set_refinement_marker_field(elementMarkerField, Refinement::NOTHING);
}

template<int DIM>
class RefinementFixture : public StkMeshFixture<DIM>
{
public:
  RefinementFixture()
  : myRefinement(mesh().mesh_meta_data(), &this->get_aux_meta().active_part())
  {
    stk::mesh::MetaData & meta = mesh().mesh_meta_data();
    stk::mesh::FieldBase & field = meta.declare_field<int>(stk::topology::ELEMENT_RANK, myElementMarkerFieldName, 1);
    myElementMarkerField = FieldRef(field);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), 1, 1, nullptr);
  }
  stk::mesh::BulkData & mesh() { return this->mMesh; }
  StkMeshBuilder<DIM> & builder() { return this->mBuilder; }
  const stk::ParallelMachine & comm() { return this->mComm; }
  void clear_refinement_marker()
  {
    clear_refinement_marker_field(myElementMarkerField);
  }
  void mark_elements_for_refinement(const std::vector<stk::mesh::Entity> & elements)
  {
    for (auto && elem : elements)
    {
      int * elemMarker = field_data<int>(myElementMarkerField, elem);
      *elemMarker = Refinement::REFINE;
    }
  }
  void mark_elements_for_unrefinement(const std::vector<stk::mesh::Entity> & elements)
  {
    for (auto && elem : elements)
    {
      int * elemMarker = field_data<int>(myElementMarkerField, elem);
      *elemMarker = Refinement::COARSEN;
    }
  }
  void mark_nonparent_elements()
  {
    for ( auto && bucket : mesh().get_buckets( stk::topology::ELEMENT_RANK, mesh().mesh_meta_data().locally_owned_part() ) )
    {
      const int markerValue = myRefinement.is_parent(*bucket) ? Refinement::NOTHING : Refinement::REFINE;
      auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for (size_t i=0; i<bucket->size(); ++i)
        elemMarker[i] = markerValue;
    }
  }
  void mark_percept_nonparent_elements()
  {
    stk::mesh::Part & parentPart = myPerceptRefinement->parent_part();
    for ( auto && bucket : mesh().get_buckets( stk::topology::ELEMENT_RANK, mesh().mesh_meta_data().locally_owned_part() ) )
    {
      const int markerValue = bucket->member(parentPart) ? Refinement::NOTHING : Refinement::REFINE;
      auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for (size_t i=0; i<bucket->size(); ++i)
        elemMarker[i] = markerValue;
    }
  }
  void do_refinement()
  {
    const TransitionElementEdgeMarker edgeMarker(mesh(), myRefinement, myElementMarkerField.name());
    myRefinement.do_refinement(edgeMarker);
  }
  void perform_iterations_of_uniform_refinement(const int numIterationsOfUMR)
  {
    for (int iRefine=0; iRefine<numIterationsOfUMR; ++iRefine)
    {
      myTimer.start();
      mark_nonparent_elements();

      do_refinement();
      myTimer.stop();
      std::cout << "After " << iRefine+1 << " levels of refinement, there are " << get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK) << " elements, time = " << myTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
    }
  }
  void setup_percept_refinement()
  {
    if (myPerceptRefinement)
      return;

    stk::mesh::MetaData & meta = mesh().mesh_meta_data();

    meta.enable_late_fields();
    myPerceptRefinement = &create_percept_refinement(meta, myTimer);
    meta.disable_late_fields();
  }
  void perform_iterations_of_percept_uniform_refinement(const int numIterationsOfUMR)
  {
    for (int iRefine=0; iRefine<numIterationsOfUMR; ++iRefine)
    {
      myTimer.start();
      mark_percept_nonparent_elements();
      HAdapt::do_adaptive_refinement(mesh().mesh_meta_data(), myElementMarkerField.name());
      myTimer.stop();
      std::cout << "After " << iRefine+1 << " levels of percept refinement, there are " << get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK) << " elements, time = " << myTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
    }
  }
  void write_mesh(const std::string &fileName)
  {
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(mesh());

    Ioss::PropertyManager properties;
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));

    size_t outputFileIndex = stkIo.create_output_mesh(fileName, stk::io::WRITE_RESULTS, properties);
    stkIo.set_active_selector(this->get_aux_meta().active_part());
    stkIo.set_subset_selector(outputFileIndex, this->get_aux_meta().active_part());
    stkIo.write_output_mesh(outputFileIndex);
    stkIo.begin_output_step(outputFileIndex, 0.);
    stkIo.write_defined_output_fields(outputFileIndex);
    stkIo.end_output_step(outputFileIndex);
  }

  void refine_marked_elements(const bool usePercept, const std::string fileName = "")
  {
    if (usePercept)
      setup_percept_refinement();

    if (usePercept)
      HAdapt::do_adaptive_refinement(mesh().mesh_meta_data(), myElementMarkerField.name());
    else
      do_refinement();

    if (!fileName.empty())
      write_mesh(fileName);
  }

  bool element_spans_x_equal_0(const stk::mesh::Entity elem)
  {
    static constexpr double tol{1.e-9};
    bool hasNeg = false;
    bool hasPos = false;
    for (auto && node : StkMeshEntities{mesh().begin_nodes(elem), mesh().end_nodes(elem)})
    {
      const stk::math::Vector3d nodeCoords = this->get_node_coordinates(node);
      if (std::abs(nodeCoords[0]) < tol)
        return true;
      else if (nodeCoords[0] < 0)
        hasNeg = true;
      else
        hasPos = true;
    }
    return hasNeg && hasPos;
  }

  void mark_elements_spanning_x_equal_0()
  {
    clear_refinement_marker();

    for ( auto && bucket : mesh().get_buckets( stk::topology::ELEMENT_RANK, mesh().mesh_meta_data().locally_owned_part() ) )
    {
      int * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for ( size_t iElem=0; iElem<bucket->size(); ++iElem )
      {
        if (element_spans_x_equal_0((*bucket)[iElem]))
          elemMarker[iElem] = Refinement::REFINE;
      }
    }
  }

  std::vector<stk::mesh::Entity> get_elements_with_given_ids(const std::vector<stk::mesh::EntityId> & idsOfElems)
  {
    std::vector<stk::mesh::Entity> elementsToRefine;
    for (auto && idOfElem : idsOfElems)
    {
      stk::mesh::Entity elemToRefine = mesh().get_entity(stk::topology::ELEMENT_RANK, idOfElem);
      if (mesh().is_valid(elemToRefine) && mesh().bucket(elemToRefine).owned())
        elementsToRefine.push_back(elemToRefine);
    }
    return elementsToRefine;
  }

  void refine_elements(const bool usePercept, const std::vector<stk::mesh::Entity> & elemsToRefine, const std::string fileName = "")
  {
    clear_refinement_marker();
    mark_elements_for_refinement(elemsToRefine);
    refine_marked_elements(usePercept, fileName);
  }

  void refine_elements_with_given_ids(const bool usePercept, const std::vector<stk::mesh::EntityId> & idsOfElemsToRefine, const std::string fileName = "")
  {
    refine_elements(usePercept, get_elements_with_given_ids(idsOfElemsToRefine), fileName);
  }

  void refine_elements_with_given_indices(const bool usePercept, const std::vector<unsigned> & indicesOfElemsToRefine, const std::string fileName = "")
  {
    refine_elements_with_given_ids(usePercept, builder().get_ids_of_elements_with_given_indices(indicesOfElemsToRefine), fileName);
  }
  double compute_mesh_quality()
  {
    const ScaledJacobianQualityMetric qualityMetric;
    return krino::compute_mesh_quality(mesh(), this->get_aux_meta().active_part(), qualityMetric);
  }
  void test_refinement_of_given_elements(const bool usePercept, const std::vector<unsigned> & indicesOfElemsToRefine, const unsigned goldNumElements, const unsigned goldNumNodes, const double goldQuality, const std::string fileName = "")
  {
    if(stk::parallel_machine_size(comm()) <= 4)
    {
      refine_elements_with_given_indices(usePercept, indicesOfElemsToRefine, fileName);

      EXPECT_EQ(goldNumElements, get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK));
      EXPECT_EQ(goldNumNodes, get_global_num_entities(mesh(), stk::topology::NODE_RANK));

      const double quality = compute_mesh_quality();
      EXPECT_NEAR(quality, goldQuality, 0.02);
    }
  }
  void mark_all_elements_for_unrefinement()
  {
    set_refinement_marker_field(myElementMarkerField, Refinement::COARSEN);
  }
  void unrefine_mesh(const bool usePercept)
  {
    if (usePercept)
    {
      bool converged = false;
      while (!converged)
      {
        const unsigned numElementsBefore = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
        stk::mesh::Part & parentPart = myPerceptRefinement->parent_part();
        for ( auto && bucket : mesh().get_buckets( stk::topology::ELEMENT_RANK, mesh().mesh_meta_data().locally_owned_part() ) )
        {
          const int markerValue = bucket->member(parentPart) ? Refinement::NOTHING : Refinement::COARSEN;
          auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
          for (size_t i=0; i<bucket->size(); ++i)
            elemMarker[i] = markerValue;
        }
        HAdapt::do_adaptive_refinement(mesh().mesh_meta_data(), myElementMarkerField.name());
        const unsigned numElementsAfter = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
        converged = numElementsAfter == numElementsBefore;
      }
    }
    else
    {
      myRefinement.fully_unrefine_mesh();
    }
  }
  std::vector<stk::mesh::Entity> get_children(const bool usePercept, const stk::mesh::Entity elem)
  {
    if (usePercept)
    {
      std::vector<stk::mesh::Entity> children;
      myPerceptRefinement->fill_children(elem, children);
      return children;
    }
    return myRefinement.get_children(elem);
  }
  unsigned get_num_children(const bool usePercept, const stk::mesh::Entity elem)
  {
    if (usePercept)
    {
      myPerceptRefinement->get_num_children(elem);
    }
    return myRefinement.get_num_children(elem);
  }
  void test_refinement_of_transition_element_leads_to_refinement_of_parent(const bool usePercept, const int indexOfCenterElement)
  {
    const stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[indexOfCenterElement]);

    const unsigned numEdges = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK) - 1;
    for (int iCaseId=0; iCaseId<(1<<numEdges); ++iCaseId)
    {
      std::vector<unsigned> edgeElementsToRefine;
      for (unsigned iEdge=0; iEdge<numEdges; ++iEdge)
        if (iCaseId & (1<<iEdge))
          edgeElementsToRefine.push_back(iEdge);

      refine_elements_with_given_indices(usePercept, edgeElementsToRefine);

      std::vector<stk::mesh::Entity> transitionElements = get_children(usePercept, centerElem);
      const unsigned numTransitionElements = transitionElements.size();

      unrefine_mesh(usePercept);

      for (unsigned iTransitionElement=0; iTransitionElement<numTransitionElements; ++iTransitionElement)
      {
        refine_elements_with_given_indices(usePercept, edgeElementsToRefine);
        transitionElements = get_children(usePercept, centerElem);

        ASSERT_EQ(numTransitionElements, transitionElements.size()) << "Number of transition elements changed from " << numTransitionElements << " to " << transitionElements.size() << std::endl;

        refine_elements_with_given_ids(usePercept, {mesh().identifier(transitionElements[iTransitionElement])});
        if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
        {
          const unsigned numChildrenAfterRefinementOfTransition = (get_children(usePercept, centerElem)).size();
          EXPECT_EQ(myRefinement.get_num_children_when_fully_refined(centerElem), numChildrenAfterRefinementOfTransition);
        }

        unrefine_mesh(usePercept);
      }
    }
  }
protected:
  stk::diag::TimerSet myTimerSet{sierra::Diag::TIMER_ALL};
  stk::diag::Timer myTimer{stk::diag::createRootTimer("Refinement", myTimerSet)};
  Refinement myRefinement;
  std::string myElementMarkerFieldName{"ELEMENT_MARKER"};
  FieldRef myElementMarkerField;
  PerceptRefinement * myPerceptRefinement{nullptr};
};

class RegularTriRefinement : public RefinementFixture<2>
{
public:
  RegularTriRefinement()
  {
    RegularTri meshSpec;
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

class RegularTetRefinement : public RefinementFixture<3>
{
public:
  RegularTetRefinement()
  {
    RegularTet meshSpec;
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

class FourRightTetsRefinement : public RefinementFixture<3>
{
public:
  FourRightTetsRefinement()
  {
    FourRightTets meshSpec;
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

    const TransitionElementEdgeMarker edgeMarker(mesh(), myRefinement, myElementMarkerField.name());
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

    const TransitionElementEdgeMarker edgeMarker(mesh(), myRefinement, myElementMarkerField.name());
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

class RightTriSurroundedByEdgeTrisRefinement : public RefinementFixture<2>
{
public:
  RightTriSurroundedByEdgeTrisRefinement()
  {
    RightTriSurroundedByEdgeTris meshSpec;
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

TEST_F(RightTriSurroundedByEdgeTrisRefinement, refinementOfOneTriInParallel_expectEdgeIsRefinedAndCoordinatesAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    const unsigned indexOfElemToRefine = 0;
    const stk::mesh::Entity edgeNode0 = mesh().get_entity(stk::topology::NODE_RANK, builder().get_assigned_node_global_ids()[3]);
    const stk::mesh::Entity edgeNode1 = mesh().get_entity(stk::topology::NODE_RANK, builder().get_assigned_node_global_ids()[5]);

    refine_elements_with_given_indices(false, {indexOfElemToRefine});

    EXPECT_EQ(10u, get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK));
    EXPECT_EQ(9u, get_global_num_entities(mesh(), stk::topology::NODE_RANK));

    const bool isEdgeNodeOwnedOrShared0 = mesh().is_valid(edgeNode0) && (mesh().bucket(edgeNode0).owned() || mesh().bucket(edgeNode0).shared());
    const bool isEdgeNodeOwnedOrShared1 = mesh().is_valid(edgeNode1) && (mesh().bucket(edgeNode1).owned() || mesh().bucket(edgeNode1).shared());

    if (isEdgeNodeOwnedOrShared0 && isEdgeNodeOwnedOrShared1)
    {
      const Edge edge = edge_from_edge_nodes(mesh(), edgeNode0, edgeNode1);

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

    stk::mesh::create_all_sides(mesh(), mesh().mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(builder().check_boundary_sides());
    EXPECT_TRUE(builder().check_block_boundary_sides());

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

class RightTetSurroundedByEdgeTetsRefinement : public RefinementFixture<3>
{
public:
  RightTetSurroundedByEdgeTetsRefinement()
  {
    RightTetSurroundedByEdgeTets meshSpec;
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

class RightTetSurroundedByFaceTetsRefinement : public RefinementFixture<3>
{
public:
  RightTetSurroundedByFaceTetsRefinement()
  {
    RightTetSurroundedByFaceTets meshSpec;
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

TEST_F(RightTetSurroundedByFaceTetsRefinement, checkAllPossibleRefinementsInParallel_expectBoundarySidesAreCorrect)
{
  for (int i=0; i<64; ++i)
  {
    std::vector<unsigned> edgeElementsToRefine;
    for (unsigned iEdge=0; iEdge<6; ++iEdge)
      if (i & (1<<iEdge))
        edgeElementsToRefine.push_back(iEdge);

    refine_elements_with_given_indices(false, edgeElementsToRefine);

    stk::mesh::create_all_sides(mesh(), mesh().mesh_meta_data().universal_part(), stk::mesh::PartVector{}, false);
    EXPECT_TRUE(builder().check_boundary_sides());
    EXPECT_TRUE(builder().check_block_boundary_sides());

    unrefine_mesh(false);
  }
}

static stk::topology get_refinement_topology(const stk::topology baseTopology)
{
  switch(baseTopology())
  {
    case stk::topology::TRIANGLE_3_2D:
      return stk::topology::TRIANGLE_6_2D;
    case stk::topology::TETRAHEDRON_4:
      return stk::topology::TETRAHEDRON_10;
    default:
        ThrowRuntimeError("Element topology not found in get_refinement_topology: " << baseTopology.name());
  }
}

void fill_permutations_and_permuted_case_ids(const stk::topology baseTopo, std::vector<unsigned> & permutations, std::vector<unsigned> & permutedCaseIds)
{
  const stk::topology topo = get_refinement_topology(baseTopo);
  const unsigned numEdges = topo.num_edges();
  const unsigned numNodes = topo.num_nodes();

  std::cout << "Building permutation tables for " << topo.name() << " with " << numEdges << " edges." << std::endl;

  const unsigned numCases = 1<<numEdges;
  permutations.resize(numCases);
  permutedCaseIds.resize(numCases);

  const unsigned numBaseNodes = numNodes-numEdges;
  std::vector<unsigned> permutedNodeIds(numNodes);

  for (unsigned i=0; i<numCases; ++i)
  {
    unsigned lowestCaseId = numCases;
    unsigned lowestPermutation = topo.num_positive_permutations();
    for (unsigned iPerm=0; iPerm<topo.num_positive_permutations(); ++iPerm)
    {
      topo.permutation_node_ordinals(iPerm, permutedNodeIds.data());

      unsigned permutedCaseId = 0;
      for (unsigned iEdge=0; iEdge<numEdges; ++iEdge)
      {
        const int permutedEdge = permutedNodeIds[iEdge+numBaseNodes]-numBaseNodes;
        if (i & (1<<permutedEdge))
          permutedCaseId += 1<<iEdge;
      }

      if (permutedCaseId < lowestCaseId)
      {
        lowestCaseId = permutedCaseId;
        lowestPermutation = iPerm;
      }
    }

    permutations[i] = lowestPermutation;
    permutedCaseIds[i] = lowestCaseId;
  }

  std::cout << "  Permutations = ";
  for (auto perm : permutations) std::cout << perm << ",";
  std::cout << std::endl;

  std::cout << "  Permuted case ids = ";
  for (auto permutedCaseId : permutedCaseIds) std::cout << permutedCaseId << ",";
  std::cout << std::endl;
}

TEST(Refinement,permutationTables)
{
  stk::ParallelMachine comm{MPI_COMM_WORLD};
  if(stk::parallel_machine_size(comm) == 1)
  {
    std::vector<unsigned> permutations;
    std::vector<unsigned> permutedCaseIds;
    fill_permutations_and_permuted_case_ids(stk::topology::TRIANGLE_3_2D, permutations, permutedCaseIds);
    for (unsigned iCase=0; iCase<permutations.size(); ++iCase)
    {
      EXPECT_EQ(permutations[iCase], TriRefiner::determine_permutation_tri3(iCase));
      EXPECT_EQ(permutedCaseIds[iCase], TriRefiner::determine_permuted_case_id_tri3(iCase));
    }
    fill_permutations_and_permuted_case_ids(stk::topology::TETRAHEDRON_4, permutations, permutedCaseIds);
    for (unsigned iCase=0; iCase<permutations.size(); ++iCase)
    {
      EXPECT_EQ(permutations[iCase], moab::SimplexTemplateRefiner::determine_permutation_tet4(iCase));
      EXPECT_EQ(permutedCaseIds[iCase], moab::SimplexTemplateRefiner::determine_permuted_case_id_tet4(iCase));
    }
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, afterEachOfThreeRoundsOfRefinementOfEdgeElements_centerElementHasCorrectNumberOfChildren)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[3]);

    const std::array<unsigned,3> elementsToRefineByRefinementLevel = {0,1,2};
    const std::array<unsigned,3> goldNumCenterElemChildrenByRefinementLevel = {2,3,4};

    for (unsigned i=0; i<3; ++i)
    {
      refine_elements_with_given_indices(false, {elementsToRefineByRefinementLevel[i]});

      if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
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
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[3]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(false, elementsToRefine);

    elementsToRefine.clear();
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(false, centerElem);
      ASSERT_EQ(4u, childElems.size());
      elementsToRefine.assign({childElems[0]});
    }

    refine_elements(false, elementsToRefine);

    EXPECT_FALSE(myRefinement.have_any_hanging_refined_nodes());

    const unsigned goldNumElems = 4 + 4+3*2 + 4+2+2*(6-2);
    EXPECT_EQ(goldNumElems, get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RightTriSurroundedByEdgeTrisRefinement, refineCenterElemAndThenMarkTransitionElement_parentOfTransitionElementGetsRefined)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[3]);
    stk::mesh::Entity edgeElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[0]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(false, elementsToRefine);

    elementsToRefine.clear();
    if (mesh().is_valid(edgeElem) && mesh().bucket(edgeElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(false, edgeElem);
      ASSERT_EQ(2u, childElems.size());
      elementsToRefine.assign({childElems[0]});
    }

    refine_elements(false, elementsToRefine);

    EXPECT_FALSE(myRefinement.have_any_hanging_refined_nodes());

    if (mesh().is_valid(edgeElem) && mesh().bucket(edgeElem).owned())
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
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[3]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(true, elementsToRefine);

    const unsigned goldNumElemsAfterRefinement = 4 + 4 + 6;
    EXPECT_EQ(goldNumElemsAfterRefinement, get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK));

    clear_refinement_marker();
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      std::vector<stk::mesh::Entity> childElems = get_children(true, centerElem);
      EXPECT_EQ(4u, childElems.size());
      mark_elements_for_unrefinement(childElems);
    }

    refine_marked_elements(true);

    const unsigned goldNumElemsAfterUnrefinement = 4;
    EXPECT_EQ(goldNumElemsAfterUnrefinement, get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK));
  }
}

TEST_F(RightTetSurroundedByEdgeTetsRefinement, refineCenterElemAndThenChildOfCenterElem_noHangingNodes)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[0]);

    std::vector<stk::mesh::Entity> elementsToRefine;
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
    {
      elementsToRefine.assign({centerElem});
    }

    refine_elements(false, elementsToRefine);

    elementsToRefine.clear();
    if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
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
    stk::mesh::Entity centerElem = mesh().get_entity(stk::topology::ELEMENT_RANK, builder().get_assigned_element_global_ids()[0]);

    const std::array<unsigned,6> elementsToRefineByRefinementLevel = {1,2,3,4,5,6};
    const std::array<unsigned,6> goldNumCenterElemChildrenByRefinementLevel = {2,3,4,6,7,8};

    for (unsigned i=0; i<6; ++i)
    {
      refine_elements_with_given_indices(false, {elementsToRefineByRefinementLevel[i]});

      if (mesh().is_valid(centerElem) && mesh().bucket(centerElem).owned())
      {
        EXPECT_EQ(goldNumCenterElemChildrenByRefinementLevel[i], myRefinement.get_num_children(centerElem));
      }
    }
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

std::string get_filename(const std::string & baseName, int fileIndex)
{
  char counterChar[32] = {'\0'};
  std::sprintf(counterChar, "%.4d", fileIndex+1);
  std::string filename = baseName + ".e-s" + std::string(counterChar);
  return filename;
}

TEST_F(RegularTriRefinement, refinementThenUnrefinementTest)
{
  std::vector<size_t> goldNumElementsByRefinementLevel = {5, 21, 69, 181, 421, 917, 1925};



  if(stk::parallel_machine_size(mComm) <= 4)
  {
    for (size_t i=0; i<goldNumElementsByRefinementLevel.size(); ++i)
    {
      mark_elements_spanning_x_equal_0();

      refine_marked_elements(false);

      const size_t numElements = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
      EXPECT_EQ(goldNumElementsByRefinementLevel[i], numElements);

      std::cout << "Refinement Level " << i+1 << ", num elements = " << numElements << std::endl;
    }

    //write_mesh("test.e");

    const bool doUnrefinement = true;
    if (doUnrefinement)
    {
      std::vector<size_t> goldNumElementsByUnrefinementLevel = {917, 421, 181, 69, 21, 5, 1};

      for (size_t i=0; i<goldNumElementsByUnrefinementLevel.size(); ++i)
      {
        mark_all_elements_for_unrefinement();

        refine_marked_elements(false);

        const size_t numElements = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
        EXPECT_EQ(goldNumElementsByUnrefinementLevel[i], numElements);

        std::cout << "Unrefinement Level " << i+1 << ", num elements = " << numElements << std::endl;
      }
    }
  }
}

TEST_F(RegularTetRefinement, DISABLED_performanceRefinementThenUnrefinementTest)
{
  // As of 11/15/2022 (first time with refinement followed by unrefinement)
  // Time: 19544 ms
  // Percept: 86114 ms
  std::vector<size_t> goldNumElementsByRefinementLevel = {9, 73, 513, 2893, 13909, 61147, 255891, 1045095};

  const bool doWriteMesh = false;
  const bool usePercept = false;

  if (doWriteMesh)
    write_mesh("test.e");

  int count = 0;
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    for (size_t i=0; i<goldNumElementsByRefinementLevel.size(); ++i)
    {
      mark_elements_spanning_x_equal_0();

      if (doWriteMesh)
        refine_marked_elements(usePercept, get_filename("test", ++count));
      else
        refine_marked_elements(usePercept);

      const size_t numElements = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
      EXPECT_EQ(goldNumElementsByRefinementLevel[i], numElements);

      std::cout << "Refinement Level " << i+1 << ", num elements = " << numElements << std::endl;
    }

    const bool doUnrefinement = true;
    if (doUnrefinement)
    {
      std::vector<size_t> goldNumElementsByUnrefinementLevel = {287809, 76441, 18165, 3649, 577, 73, 9, 1};

      for (size_t i=0; i<goldNumElementsByUnrefinementLevel.size(); ++i)
      {
        mark_all_elements_for_unrefinement();

        if (doWriteMesh)
          refine_marked_elements(usePercept, get_filename("test", ++count));
        else
          refine_marked_elements(usePercept);

        const size_t numElements = get_global_num_entities(mesh(), stk::topology::ELEMENT_RANK);
        EXPECT_EQ(goldNumElementsByUnrefinementLevel[i], numElements);

        std::cout << "Unrefinement Level " << i+1 << ", num elements = " << numElements << std::endl;
      }
    }
  }
}

}
