#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_

#include <Akri_StkMeshFixture.hpp>
#include <random>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_AllReduce.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Refinement.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_TransitionElementEdgeMarker.hpp>
#include <Akri_AuxMetaData.hpp>

namespace krino {

inline void set_refinement_marker_field(FieldRef elementMarkerField, const Refinement::RefinementMarker value)
{
  stk::mesh::field_fill(static_cast<int>(value), elementMarkerField, stk::mesh::selectField(elementMarkerField));
}

inline void clear_refinement_marker_field(FieldRef elementMarkerField)
{
  set_refinement_marker_field(elementMarkerField, Refinement::RefinementMarker::NOTHING);
}


template <typename MESHSPEC>
class RefinementFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  RefinementFixture()
  : myRefinement(mMesh.mesh_meta_data(), &this->get_aux_meta().active_part(), sierra::Diag::sierraTimer())
  {
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    stk::mesh::FieldBase & elemMarkerField = meta.declare_field<int>(stk::topology::ELEMENT_RANK, myElementMarkerFieldName, 1);
    stk::mesh::FieldBase & elemField = meta.declare_field<double>(stk::topology::ELEMENT_RANK, myElemFieldName, 1);
    myElementMarkerField = FieldRef(elemMarkerField);
    myElemField = FieldRef(elemField);
    stk::mesh::put_field_on_mesh(elemMarkerField, meta.universal_part(), 1, 1, nullptr);
    stk::mesh::put_field_on_mesh(elemField, meta.universal_part(), 1, 1, nullptr);
    mMesh.set_automatic_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::write_mesh;

  void mark_elements_for_refinement(const std::vector<stk::mesh::Entity> & elements)
  {
    for (auto && elem : elements)
    {
      int * elemMarker = field_data<int>(myElementMarkerField, elem);
      *elemMarker = static_cast<int>(Refinement::RefinementMarker::REFINE);
    }
  }
  void mark_elements_for_unrefinement(const std::vector<stk::mesh::Entity> & elements)
  {
    for (auto && elem : elements)
    {
      int * elemMarker = field_data<int>(myElementMarkerField, elem);
      *elemMarker = static_cast<int>(Refinement::RefinementMarker::COARSEN);
    }
  }
  void mark_nonparent_elements()
  {
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() ) )
    {
      const auto markerValue = myRefinement.is_parent(*bucket) ? Refinement::RefinementMarker::NOTHING : Refinement::RefinementMarker::REFINE;
      auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for (size_t i=0; i<bucket->size(); ++i)
        elemMarker[i] = static_cast<int>(markerValue);
    }
  }

  bool do_refinement()
  {
    const TransitionElementEdgeMarker edgeMarker(mMesh, myRefinement, myElementMarkerField.name());
    return myRefinement.do_refinement(edgeMarker);
  }

  void perform_iterations_of_uniform_refinement_with_general_element_marker(const int numIterationsOfUMR)
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

  void perform_iterations_of_uniform_refinement_with_uniform_marker(const int numIterationsOfUMR)
  {
    myTimer.start();
    myRefinement.do_uniform_refinement(numIterationsOfUMR);
    myTimer.stop();
    const size_t numElems = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK);
    if (0 == stk::parallel_machine_rank(mComm))
      std::cout << "After " << numIterationsOfUMR << " levels of uniform refinement, there are " << numElems << " elements, time = " << myTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
  }

  void refine_marked_elements(const std::string fileName = "")
  {
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

  std::pair<bool, double> tag_element_spans_z_equal_0_and_compute_element_field(const stk::mesh::Entity elem, bool flip)
  {
    static constexpr double tol{1.e-9};
    bool hasNeg = false;
    bool hasPos = false;
    bool hasNodeOnZEqual0 = false;
    double centroidValue = 0.0;
    for (auto && node : StkMeshEntities{mMesh.begin_nodes(elem), mMesh.end_nodes(elem)})
    {
      const stk::math::Vector3d nodeCoords = this->get_node_coordinates(node);
      centroidValue += nodeCoords[2];

      if (std::abs(nodeCoords[2]) < tol)
        hasNodeOnZEqual0 = true;
      else if (nodeCoords[2] < 0)
        hasNeg = true;
      else
        hasPos = true;
    }
    centroidValue /= (double)mMesh.num_nodes(elem);
    if (hasNodeOnZEqual0) return std::pair<bool, double>{true,centroidValue};
    return std::pair<bool, double>{hasNeg && hasPos, centroidValue};
  }

  void mark_elements_spanning_z_equal_0_and_populate_elem_field(bool flip)
  {
    clear_refinement_marker();
    AuxMetaData & auxMeta = AuxMetaData::get(mMesh.mesh_meta_data());
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() & auxMeta.active_part()) )
    {
      int * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      double * elemField = field_data<double>(myElemField, *bucket);
      for ( size_t iElem=0; iElem<bucket->size(); ++iElem )
      {
        auto doMarkAndFieldVal = tag_element_spans_z_equal_0_and_compute_element_field((*bucket)[iElem], flip);
        auto doMark = std::get<0>(doMarkAndFieldVal);
        auto centroidValue = std::get<1>(doMarkAndFieldVal);
        if (doMark) elemMarker[iElem] = static_cast<int>(Refinement::RefinementMarker::REFINE);

        if (centroidValue <= 0 ) elemField[iElem] = -1;
        else elemField[iElem] = 1;

        if(flip) elemField[iElem] *= -1.0;
      }
    }
  }

  void check_elem_field_values_spanning_z_equal_0(bool flip)
  {
    clear_refinement_marker();
    AuxMetaData & auxMeta = AuxMetaData::get(mMesh.mesh_meta_data());
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() & auxMeta.active_part()) )
    {
      double * elemField = field_data<double>(myElemField, *bucket);
      for ( size_t iElem=0; iElem<bucket->size(); ++iElem )
      {
        auto doMarkAndFieldVal = tag_element_spans_z_equal_0_and_compute_element_field((*bucket)[iElem], flip);
        auto centroidValue = std::get<1>(doMarkAndFieldVal);
        double goldValue;
        if (centroidValue < 0 ) goldValue = -1.0;
        else goldValue = 1.0;

        if(flip) goldValue *= -1.0;
        EXPECT_EQ(goldValue, elemField[iElem]);
      }
    }
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
          elemMarker[iElem] = static_cast<int>(Refinement::RefinementMarker::REFINE);
      }
    }
  }

  void randomly_select_children(std::mt19937 & rand_gen,
      const double select_fraction,
      std::vector<stk::mesh::Entity> & elems_out)
  {
    elems_out.clear();

    std::uniform_real_distribution<> rand_dist(0., 1.);
    const auto & meta = mMesh.mesh_meta_data();
    const stk::mesh::Selector selector = meta.locally_owned_part() & myRefinement.child_part();
    for (auto && bucket : mMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
    {
      for (size_t iElem = 0; iElem < bucket->size(); ++iElem)
      {
        const double rand_val = rand_dist(rand_gen);
        if (rand_val < select_fraction)
        {
          elems_out.push_back((*bucket)[iElem]);
        }
      }
    }
  }

  void randomly_mark_elements(std::mt19937 &rand_gen)
  {
    clear_refinement_marker();

    std::uniform_real_distribution<> rand_dist(0., 1.);
    const double refine_prob = 0.1;
    const double unrefine_prob = 0.95;
    const auto & meta = mMesh.mesh_meta_data();
    const stk::mesh::Selector owned_leaf = meta.locally_owned_part() & !myRefinement.parent_part();
    const stk::mesh::Selector owned_parent = meta.locally_owned_part() & myRefinement.parent_part();
    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, owned_leaf ) )
    {
      int * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for ( size_t iElem=0; iElem<bucket->size(); ++iElem )
      {
        const double rand_val = rand_dist(rand_gen);
        if(rand_val < refine_prob)
        {
          elemMarker[iElem] = static_cast<int>(Refinement::RefinementMarker::REFINE);
        }
      }
    }

    for ( auto && bucket : mMesh.get_buckets( stk::topology::ELEMENT_RANK, owned_parent ) )
    {
      int * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for ( size_t iElem=0; iElem<bucket->size(); ++iElem )
      {
        const double rand_val = rand_dist(rand_gen);
        if(rand_val < unrefine_prob)
        {
          elemMarker[iElem] = static_cast<int>(Refinement::RefinementMarker::COARSEN);
          for(auto && child : myRefinement.get_children((*bucket)[iElem]))
          {
            *field_data<int>(myElementMarkerField, child) = static_cast<int>(Refinement::RefinementMarker::COARSEN);
          }
        }
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

  void refine_elements(const std::vector<stk::mesh::Entity> & elemsToRefine, const std::string fileName = "")
  {
    clear_refinement_marker();
    mark_elements_for_refinement(elemsToRefine);
    refine_marked_elements(fileName);
  }

  void refine_elements_with_given_ids(const std::vector<stk::mesh::EntityId> & idsOfElemsToRefine, const std::string fileName = "")
  {
    refine_elements(get_elements_with_given_ids(idsOfElemsToRefine), fileName);
  }

  void refine_elements_with_given_indices(const std::vector<unsigned> & indicesOfElemsToRefine, const std::string fileName = "")
  {
    refine_elements_with_given_ids(mBuilder.get_ids_of_elements_with_given_indices(indicesOfElemsToRefine), fileName);
  }

  double compute_mesh_quality()
  {
    const ScaledJacobianQualityMetric qualityMetric;
    return krino::compute_mesh_quality(mMesh, this->get_aux_meta().active_part(), qualityMetric);
  }

  void unrefine_mesh()
  {
    myRefinement.fully_unrefine_mesh();
  }
  std::vector<stk::mesh::Entity> get_children(const stk::mesh::Entity elem)
  {
    return myRefinement.get_children(elem);
  }
  unsigned get_num_children(const stk::mesh::Entity elem)
  {
    return myRefinement.get_num_children(elem);
  }

  void clear_refinement_marker()
  {
    clear_refinement_marker_field(myElementMarkerField);
  }

  void mark_all_elements_for_unrefinement()
  {
    set_refinement_marker_field(myElementMarkerField, Refinement::RefinementMarker::COARSEN);
  }

  void mark_all_elements_for_refinement()
  {
    set_refinement_marker_field(myElementMarkerField, Refinement::RefinementMarker::REFINE);
  }

  void test_refinement_of_transition_element_leads_to_refinement_of_parent(const int indexOfCenterElement)
  {
    const stk::mesh::Entity centerElem = mMesh.get_entity(stk::topology::ELEMENT_RANK, mBuilder.get_assigned_element_global_ids()[indexOfCenterElement]);
    std::vector<stk::mesh::Entity> transitionElements;

    const unsigned numEdges = get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK) - 1;
    for (int iCaseId=0; iCaseId<(1<<numEdges); ++iCaseId)
    {
      std::vector<unsigned> edgeElementsToRefine;
      for (unsigned iEdge=0; iEdge<numEdges; ++iEdge)
        if (iCaseId & (1<<iEdge))
          edgeElementsToRefine.push_back(iEdge);

      refine_elements_with_given_indices(edgeElementsToRefine);

      unsigned numTransitionElements = 0;
      if (mMesh.is_valid(centerElem))
        numTransitionElements = get_num_children(centerElem);
      all_reduce_max(mMesh.parallel(), numTransitionElements);

      unrefine_mesh();

      for (unsigned iTransitionElement=0; iTransitionElement<numTransitionElements; ++iTransitionElement)
      {
        refine_elements_with_given_indices(edgeElementsToRefine);
        std::vector<stk::mesh::Entity> elementsToRefine;
        if (mMesh.is_valid(centerElem))
        {
          transitionElements = get_children(centerElem);

          ASSERT_EQ(numTransitionElements, transitionElements.size()) << "Number of transition elements changed from " << numTransitionElements << " to " << transitionElements.size() << std::endl;

          elementsToRefine.push_back(transitionElements[iTransitionElement]);
        }

        refine_elements(elementsToRefine);
        if (mMesh.is_valid(centerElem) && mMesh.bucket(centerElem).owned())
        {
          const unsigned numChildrenAfterRefinementOfTransition = (get_children(centerElem)).size();
          EXPECT_EQ(myRefinement.get_num_children_when_fully_refined(centerElem), numChildrenAfterRefinementOfTransition);
        }

        unrefine_mesh();
      }
    }
  }

  void test_refinement_of_given_elements(const std::vector<unsigned> & indicesOfElemsToRefine, const unsigned goldNumElements, const unsigned goldNumNodes, const double goldQuality, const std::string fileName = "")
  {
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      refine_elements_with_given_indices(indicesOfElemsToRefine, fileName);

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
    fix_node_owners_to_assure_active_owned_element_for_node(mMesh, this->get_aux_meta().active_part());
  }


protected:
  MESHSPEC meshSpec;
  stk::diag::Timer myTimer{"Refinement", sierra::Diag::sierraTimer()};
  Refinement myRefinement;
  std::string myElementMarkerFieldName{"ELEMENT_MARKER"};
  std::string myElemFieldName{"ELEMENT_FIELD"};
  FieldRef myElementMarkerField;
  FieldRef myElemField;
};

}




#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_HPP_ */
