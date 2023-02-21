#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_

#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_AdaptivityInterface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Refinement.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_TransitionElementEdgeMarker.hpp>
#include <Akri_AuxMetaData.hpp>

namespace krino {

inline void set_refinement_marker_field(FieldRef elementMarkerField, const int value)
{
  stk::mesh::field_fill(value, elementMarkerField, stk::mesh::selectField(elementMarkerField));
}

inline void clear_refinement_marker_field(FieldRef elementMarkerField)
{
  set_refinement_marker_field(elementMarkerField, Refinement::NOTHING);
}


template <typename MESHSPEC>
class RefinementFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  RefinementFixture()
  : myRefinement(mMesh.mesh_meta_data(), &this->get_aux_meta().active_part())
  {
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    stk::mesh::FieldBase & field = meta.declare_field<int>(stk::topology::ELEMENT_RANK, myElementMarkerFieldName, 1);
    myElementMarkerField = FieldRef(field);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), 1, 1, nullptr);
  }

  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;

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
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() ) )
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
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() ) )
    {
      const int markerValue = bucket->member(parentPart) ? Refinement::NOTHING : Refinement::REFINE;
      auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for (size_t i=0; i<bucket->size(); ++i)
        elemMarker[i] = markerValue;
    }
  }
  void do_refinement()
  {
    const TransitionElementEdgeMarker edgeMarker(mMesh, myRefinement, myElementMarkerField.name());
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
      std::cout << "After " << iRefine+1 << " levels of refinement, there are " << get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK) << " elements, time = " << myTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
    }
  }
  void setup_percept_refinement()
  {
    if (myPerceptRefinement)
      return;

    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();

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
      HAdapt::do_adaptive_refinement(mMesh.mesh_meta_data(), myElementMarkerField.name());
      myTimer.stop();
      std::cout << "After " << iRefine+1 << " levels of percept refinement, there are " << get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK) << " elements, time = " << myTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
    }
  }
  void write_mesh(const std::string &fileName)
  {
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(mMesh);

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

  void refine_marked_elements(const bool usePercept=false, const std::string fileName = "")
  {
    if (usePercept)
      setup_percept_refinement();

    if (usePercept)
      HAdapt::do_adaptive_refinement(mMesh.mesh_meta_data(), myElementMarkerField.name());
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
    for (auto && node : StkMeshEntities{mMesh.begin_nodes(elem), mMesh.end_nodes(elem)})
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

    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() ) )
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
      stk::mesh::Entity elemToRefine = mMesh.get_entity(stk::topology::ELEMENT_RANK, idOfElem);
      if (mMesh.is_valid(elemToRefine) && mMesh.bucket(elemToRefine).owned())
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
    refine_elements_with_given_ids(usePercept, mBuilder.get_ids_of_elements_with_given_indices(indicesOfElemsToRefine), fileName);
  }

  double compute_mesh_quality()
  {
    const ScaledJacobianQualityMetric qualityMetric;
    return krino::compute_mesh_quality(mMesh, this->get_aux_meta().active_part(), qualityMetric);
  }

  void unrefine_mesh(const bool usePercept)
  {
    if (usePercept)
    {
      bool converged = false;
      while (!converged)
      {
        const unsigned numElementsBefore = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
        stk::mesh::Part & parentPart = myPerceptRefinement->parent_part();
        for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() ) )
        {
          const int markerValue = bucket->member(parentPart) ? Refinement::NOTHING : Refinement::COARSEN;
          auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
          for (size_t i=0; i<bucket->size(); ++i)
            elemMarker[i] = markerValue;
        }
        HAdapt::do_adaptive_refinement(mMesh.mesh_meta_data(), myElementMarkerField.name());
        const unsigned numElementsAfter = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
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

  void clear_refinement_marker()
  {
    clear_refinement_marker_field(myElementMarkerField);
  }

  void mark_all_elements_for_unrefinement()
  {
    set_refinement_marker_field(myElementMarkerField, Refinement::COARSEN);
  }

  void test_refinement_of_transition_element_leads_to_refinement_of_parent(const bool usePercept, const int indexOfCenterElement)
  {
    const stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[indexOfCenterElement]);

    const unsigned numEdges = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK) - 1;
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

        refine_elements_with_given_ids(usePercept, {mMesh.identifier(transitionElements[iTransitionElement])});
        if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
        {
          const unsigned numChildrenAfterRefinementOfTransition = (get_children(usePercept, centerElem)).size();
          EXPECT_EQ(myRefinement.get_num_children_when_fully_refined(centerElem), numChildrenAfterRefinementOfTransition);
        }

        unrefine_mesh(usePercept);
      }
    }
  }

  void test_refinement_of_given_elements(const bool usePercept, const std::vector<unsigned> & indicesOfElemsToRefine, const unsigned goldNumElements, const unsigned goldNumNodes, const double goldQuality, const std::string fileName = "")
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      refine_elements_with_given_indices(usePercept, indicesOfElemsToRefine, fileName);

      EXPECT_EQ(goldNumElements, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
      EXPECT_EQ(goldNumNodes, get_global_num_entities(mMesh, stk::topology::NODE_RANK));

      const double quality = compute_mesh_quality();
      EXPECT_NEAR(quality, goldQuality, 0.02);
    }
  }

  void move_owned_elements_with_given_ids_and_owned_attached_entities_to_processor(const std::vector<stk::mesh::EntityId> & elemIds, const int proc)
  {
    std::vector<stk::mesh::EntityProc> entitiesToMove;
    for (auto elemId : elemIds)
    {
      stk::mesh::Entity elem = mMesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
      if (mMesh.is_valid(elem) && mMesh.parallel_owner_rank(elem) == mMesh.parallel_rank())
      {
        entitiesToMove.emplace_back(elem, proc);
        for (auto && relRank : {stk::topology::EDGE_RANK, stk::topology::FACE_RANK})
          for (auto edgeOrFace : StkMeshEntities{mMesh.begin(elem,relRank), mMesh.end(elem,relRank)})
            if (mMesh.parallel_owner_rank(edgeOrFace) == mMesh.parallel_rank())
              entitiesToMove.emplace_back(edgeOrFace, proc);
      }
    }
    mMesh.change_entity_owner(entitiesToMove);
  }


protected:
  MESHSPEC meshSpec;
  stk::diag::Timer myTimer{"Refinement", sierra::Diag::sierraTimer()};
  Refinement myRefinement;
  std::string myElementMarkerFieldName{"ELEMENT_MARKER"};
  FieldRef myElementMarkerField;
  PerceptRefinement * myPerceptRefinement{nullptr};
};

}




#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_ */
