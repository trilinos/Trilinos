/*
 * Akri_RefinementInterface.cpp
 *
 *  Created on: Oct 31, 2022
 *      Author: drnoble
 */
#include "Akri_RefinementInterface.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <Akri_AuxMetaData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include "Akri_DiagWriter.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_ReportHandler.hpp"
#include "Akri_TransitionElementEdgeMarker.hpp"
#include "stk_util/environment/Env.hpp"

namespace krino {

void clear_refinement_marker(const RefinementInterface & refinement)
{
  FieldRef markerField = refinement.get_marker_field_and_sync_to_host();
  stk::mesh::field_fill(static_cast<int>(Refinement::RefinementMarker::NOTHING), markerField, stk::mesh::selectField(markerField));
}

void mark_selected_elements_for_refinement(const RefinementInterface & refinement, const stk::mesh::Selector & selector)
{
  FieldRef markerField = refinement.get_marker_field_and_sync_to_host();
  clear_refinement_marker(refinement);
  stk::mesh::field_fill(static_cast<int>(Refinement::RefinementMarker::REFINE), markerField, selector);
}

void mark_selected_elements_for_refinement(const RefinementInterface & refinement, const int current_refinement_level, const int max_refinement_levels, const stk::mesh::Selector & selector)
{
  const bool doAdapt = current_refinement_level < max_refinement_levels;
  if (doAdapt)
    mark_selected_elements_for_refinement(refinement, selector);
  else
    clear_refinement_marker(refinement);
}

void mark_elements_for_refinement(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & elemsToRefine)
{
  clear_refinement_marker(refinement);
  FieldRef markerField = refinement.get_marker_field_and_sync_to_host();
  for (auto && elem : elemsToRefine)
  {
    int * elemMarker = field_data<int>(markerField, elem);
    *elemMarker = static_cast<int>(Refinement::RefinementMarker::REFINE);
  }
}

void mark_elements_with_given_ids_for_refinement(const stk::mesh::BulkData & mesh, const RefinementInterface & refinement, const std::vector<stk::mesh::EntityId> & idsOfElemsToRefine)
{
  std::vector<stk::mesh::Entity> elementsToRefine;
  for (auto && idOfElemToRefine : idsOfElemsToRefine)
  {
    stk::mesh::Entity elemToRefine = mesh.get_entity(stk::topology::ELEMENT_RANK, idOfElemToRefine);
    if (mesh.is_valid(elemToRefine) && mesh.bucket(elemToRefine).owned())
      elementsToRefine.push_back(elemToRefine);
  }

  mark_elements_for_refinement(refinement, elementsToRefine);
}

void mark_based_on_indicator_field(const stk::mesh::BulkData & mesh,
    const RefinementInterface & refinement,
     const std::string & indicatorFieldName,
    const int maxRefinementLevel,
    const int currentRefinementLevel,
    const uint64_t targetElemCount)
{
  const auto & meta = mesh.mesh_meta_data();

  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  if (!auxMeta.has_field(stk::topology::ELEMENT_RANK, indicatorFieldName))
  {
    krinolog << "Did not find registered indicator field " << indicatorFieldName << ".  Marking no elements for refinement" << stk::diag::dendl;
    clear_refinement_marker(refinement);
    return;
  }

  const FieldRef indicatorField = auxMeta.get_field(stk::topology::ELEMENT_RANK, indicatorFieldName);
  const FieldRef markerField = refinement.get_marker_field_and_sync_to_host();
  const auto & parentPart = refinement.parent_part();

  const auto & activeBuckets =
      mesh.get_buckets(stk::topology::ELEMENT_RANK,
          stk::mesh::selectField(indicatorField) &
          !parentPart &
          meta.locally_owned_part());

  const int bucketsToSample = std::min(1000ul, activeBuckets.size());
  std::vector<double> sampleValues;
  const int default_bucket_capacity = 512;
  sampleValues.reserve(default_bucket_capacity*bucketsToSample);
  for(int i=0; i < bucketsToSample; ++i)
  {
    const double * indicatorData = field_data<double>(indicatorField, *activeBuckets[i]);
    sampleValues.insert(sampleValues.end(), indicatorData, indicatorData + activeBuckets[i]->size());
  }

  // This is not scalable
  std::vector<double> globalVals;
  stk::parallel_vector_concat(mesh.parallel(), sampleValues, globalVals);
  std::sort(globalVals.begin(), globalVals.end());

  const auto currentActiveElements = globalVals.size();

  const double remaining_levels = maxRefinementLevel - currentRefinementLevel;
  const double fraction_to_refine =
      currentActiveElements >= targetElemCount ? 0:
          1./12. * (
              std::pow(
                  static_cast<double>(targetElemCount) / static_cast<double>(currentActiveElements),
                  1./remaining_levels)
              - 1.);

  const int global_fraction_index = (1.-fraction_to_refine) * currentActiveElements;
  const double threshold_val = globalVals[global_fraction_index];

  for(auto && bucketPtr : activeBuckets)
  {
    const auto & bucket = *bucketPtr;
    const double * indicatorData = field_data<double>(indicatorField, bucket);
    int * markerData = field_data<int>(markerField, bucket);
    const int size = bucket.size();
    for(int i=0; i < size; ++i)
    {
      markerData[i] = indicatorData[i] > threshold_val 
        ? static_cast<int>(Refinement::RefinementMarker::REFINE)
        : static_cast<int>(Refinement::RefinementMarker::NOTHING);
    }
  }
}

void fill_leaf_children(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & leaf_children)
{
  std::vector<stk::mesh::Entity> grand_children;
  for (auto&& child : children)
  {
    refinement.fill_children(child, grand_children);
    if (grand_children.empty())
    {
      leaf_children.push_back(child);
    }
    else
    {
      fill_leaf_children(refinement, grand_children, leaf_children);
    }
  }
}

void fill_leaf_children(const RefinementInterface & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & leaf_children)
{
  leaf_children.clear();
  std::vector<stk::mesh::Entity> children;
  refinement.fill_children(entity, children);
  fill_leaf_children(refinement, children, leaf_children);
}

void fill_all_children(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & all_children)
{
  std::vector<stk::mesh::Entity> grand_children;
  for (auto&& child : children)
  {
    refinement.fill_children(child, grand_children);
    all_children.insert(all_children.end(), grand_children.begin(), grand_children.end());
    fill_all_children(refinement, grand_children, all_children);
  }
}

void fill_all_children(const RefinementInterface & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & all_children)
{
  std::vector<stk::mesh::Entity> children;
  refinement.fill_children(entity, children);
  all_children = children;
  fill_all_children(refinement, children, all_children);
}

KrinoRefinement &
KrinoRefinement::get(const stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  STK_ThrowRequireMsg(nullptr != refinement, "Refinement not found on MetaData.");
  return *refinement;
}

KrinoRefinement &
KrinoRefinement::create(stk::mesh::MetaData & meta, stk::diag::Timer & timer)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  STK_ThrowRequireMsg(nullptr == refinement, "KrinoRefinement::create should be called only once per MetaData.");
  if (nullptr == refinement)
  {
    AuxMetaData & auxMeta = AuxMetaData::get_or_create(meta);
    refinement = new KrinoRefinement(meta,
        &auxMeta.active_part(),
        false /*auxMeta.get_force_64bit_flag()*/,
        auxMeta.get_assert_32bit_flag(),
        timer);
    meta.declare_attribute_with_delete<KrinoRefinement>(refinement);
  }
  return *refinement;
}

KrinoRefinement &
KrinoRefinement::create(stk::mesh::MetaData & meta)
{
  return KrinoRefinement::create(meta, sierra::Diag::sierraTimer());
}

bool
KrinoRefinement::is_created(const stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  return refinement != nullptr;
}

KrinoRefinement &
KrinoRefinement::get_or_create(stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  if (refinement)
    return *refinement;
  return create(meta);
}

KrinoRefinement &
KrinoRefinement::get_or_create(stk::mesh::MetaData & meta, stk::diag::Timer & timer)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  if (refinement)
    return *refinement;

  return create(meta, timer);
}


void
KrinoRefinement::register_parts_and_fields_via_aux_meta_for_fmwk(stk::mesh::MetaData & meta)
{
  // FIXME: Ugly workaround for fmwk
  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  if (auxMeta.using_fmwk())
  {
    stk::mesh::Part & RefinedEdgeNodePart = auxMeta.declare_io_part_with_topology("Refinement_Edge_Node", stk::topology::NODE, true);
    auxMeta.register_field("REFINEMENT_REFINED_EDGE_NODE_PARENTS_IDS", krino::FieldType::UNSIGNED_INTEGER_64, stk::topology::NODE_RANK, 1, 2, RefinedEdgeNodePart);
    auxMeta.register_field("REFINEMENT_LEVEL", krino::FieldType::INTEGER, stk::topology::ELEMENT_RANK, 1, 1, meta.universal_part()); // needed everywhere for restart
    auxMeta.register_field("REFINEMENT_PARENT_ELEMENT_ID", krino::FieldType::UNSIGNED_INTEGER_64, stk::topology::ELEMENT_RANK, 1, 1, meta.universal_part()); // needed everywhere for restart
    auxMeta.register_field("ORIGINATING_PROC_FOR_PARENT_ELEMENT", krino::FieldType::INTEGER, stk::topology::ELEMENT_RANK, 1, 1, meta.universal_part()); // needed everywhere for restart
  }
}

KrinoRefinement::KrinoRefinement(stk::mesh::MetaData & meta,
    stk::mesh::Part * activePart,
    const bool force64Bit,
    const bool assert32Bit,
    stk::diag::Timer & parent_timer)
    : myMeta(meta),
     myRefinementTimer("Refinement", parent_timer),
     myRefinement(meta, activePart, force64Bit, assert32Bit, myRefinementTimer)
{
  myElementMarkerField = meta.declare_field<int>(stk::topology::ELEMENT_RANK, "REFINEMENT_ELEMENT_MARKER", 1);
  stk::mesh::put_field_on_mesh(myElementMarkerField.field(), meta.universal_part(), 1, 1, nullptr);
}

KrinoRefinement::~KrinoRefinement()
{
}

bool KrinoRefinement::is_child(const stk::mesh::Entity entity) const
{
  return myRefinement.is_child(entity);
}

bool KrinoRefinement::is_parent(const stk::mesh::Entity entity) const
{
  return myRefinement.is_parent(entity);
}

bool KrinoRefinement::is_parent_side(const stk::mesh::Entity side) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  bool haveAttachedParentElement = false;
  for (auto element : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
  {
    if (is_parent(element))
      haveAttachedParentElement = true;
    else
      return false;
  }
  return haveAttachedParentElement;
}

stk::mesh::Part & KrinoRefinement::parent_part() const
{
  return myRefinement.parent_part();
}

stk::mesh::Part & KrinoRefinement::child_part() const
{
  return myRefinement.child_part();
}

stk::mesh::Entity KrinoRefinement::get_parent(const stk::mesh::Entity child) const
{
  return myRefinement.get_parent(child);
}

std::pair<stk::mesh::EntityId,int> KrinoRefinement::get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const
{
  return myRefinement.get_parent_id_and_parallel_owner_rank(child);
}

unsigned KrinoRefinement::get_num_children(const stk::mesh::Entity elem) const
{
  return myRefinement.get_num_children(elem);
}

void KrinoRefinement::fill_children(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children) const
{
  myRefinement.fill_children(parent, children);
}

void KrinoRefinement::fill_child_element_ids(const stk::mesh::Entity parent, std::vector<stk::mesh::EntityId> & childElemIds) const
{
  myRefinement.fill_child_element_ids(parent, childElemIds);
}

std::string KrinoRefinement::locally_check_leaf_children_have_parents_on_same_proc() const
{
  return myRefinement.locally_check_leaf_children_have_parents_on_same_proc();
}

void KrinoRefinement::fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const
{
  myRefinement.fill_child_elements_that_must_stay_on_same_proc_as_parent(parent, dependents);
}

bool KrinoRefinement::has_rebalance_constraint(const stk::mesh::Entity entity) const
{
  return myRefinement.has_parallel_owner_rebalance_constraint(entity);
}

void KrinoRefinement::update_element_rebalance_weights_incorporating_parallel_owner_constraints(stk::mesh::Field<double> & elemWtField) const
{
  return myRefinement.update_element_rebalance_weights_incorporating_parallel_owner_constraints(elemWtField);
}

TransitionElementEdgeMarker & KrinoRefinement::setup_marker() const
{
  // This is delayed to make sure the bulkdata has been created.  It might be better to just pass in the bulkdata when the marker is used.
  if (!myMarker)
  {
    myMarker = std::make_unique<TransitionElementEdgeMarker>(myMeta.mesh_bulk_data(), const_cast<Refinement&>(myRefinement), myElementMarkerField.name());
  }
  return *myMarker;
}

void KrinoRefinement::set_marker_field(const std::string & markerFieldName)
{
  myElementMarkerField = myMeta.get_field(stk::topology::ELEMENT_RANK, markerFieldName);
}

FieldRef KrinoRefinement::get_marker_field_and_sync_to_host() const
{
  setup_marker();
  myElementMarkerField.sync_to_host();

  return myElementMarkerField;
}

int KrinoRefinement::fully_refined_level(const stk::mesh::Entity elem) const
{
  const int refineLevel = myRefinement.refinement_level(elem);
  if (is_transition(elem))
    return refineLevel-1;
  return refineLevel;
}

int KrinoRefinement::partially_refined_level(const stk::mesh::Entity elem) const
{
  const int refineLevel = myRefinement.refinement_level(elem);
  return refineLevel;
}

bool KrinoRefinement::is_transition(const stk::mesh::Entity elem) const
{
  return get_marker().is_transition(elem);
}

TransitionElementEdgeMarker & KrinoRefinement::get_marker() const
{
  setup_marker();
  STK_ThrowRequireMsg(myMarker, "Logic error.  Marker should have already been setup.");
  return *myMarker;
}

std::array<unsigned, 2> KrinoRefinement::get_marked_element_counts() const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const FieldRef markerField = get_marker_field_and_sync_to_host();
  const stk::mesh::Selector selector = stk::mesh::selectField(markerField) & AuxMetaData::get(myMeta).active_part() & myMeta.locally_owned_part() & !myRefinement.parent_part();

  unsigned numRefine = 0;
  unsigned numUnrefine = 0;

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const int * marker = field_data<int>(markerField, *bucket);
    const int length = bucket->size();

    for (int i = 0; i < length; ++i)
    {
      if (marker[i] == static_cast<int>(Refinement::RefinementMarker::COARSEN))
        ++numUnrefine;
      else if (marker[i] == static_cast<int>(Refinement::RefinementMarker::REFINE))
        ++numRefine;
    }
  }

  const std::array<unsigned, 2> localNum{numRefine, numUnrefine};
  std::array<unsigned, 2> globalNum{0,0};
  stk::all_reduce_sum(mesh.parallel(), localNum.data(), globalNum.data(), 2);

  return globalNum;
}

bool KrinoRefinement::is_supported_uniform_refinement_element() const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const FieldRef markerField = get_marker_field_and_sync_to_host();
  const stk::mesh::Selector selector = stk::mesh::selectField(markerField) &
      AuxMetaData::get(myMeta).active_part() & myMeta.locally_owned_part() &
      !myRefinement.parent_part();

  bool is_supported = true;

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const auto topology = bucket->topology();

    is_supported &= (topology == stk::topology::TRI_3 || topology == stk::topology::TRI_3_2D ||
        topology == stk::topology::TETRAHEDRON_4 || topology == stk::topology::HEX_8 ||
        topology == stk::topology::QUAD_4 || topology == stk::topology::QUAD_4_2D ||
        topology == stk::topology::BEAM_2 || topology == stk::topology::BEAM_3);
  }

  const int is_supported_int = static_cast<int>(is_supported);
  int reduced_result = 0;
  
  MPI_Allreduce(&is_supported_int, &reduced_result, 1, MPI_INT, MPI_LAND, mesh.parallel());
  
  return static_cast<bool>(reduced_result);
}

bool KrinoRefinement::do_refinement(const int debugLevel)
{ 
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  const auto markerCounts = get_marked_element_counts();
  const unsigned numRefine = markerCounts[0];
  const unsigned numUnrefine = markerCounts[1];

  krinolog << "Number of elements marked for refinement = " << numRefine << "\n";
  krinolog << "Number of elements marked for unrefinement = " << numUnrefine << stk::diag::dendl;

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  //if it can be uniformly refined, and all elements are marked for refinement just do uniform refinement
  if (is_supported_uniform_refinement_element() && counts[3] == numRefine && (counts[3]+numRefine) > 0)
  {
    return do_uniform_refinement(1);
  }

  krinolog << "Adaptive refinement: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  bool didMakeAnyChanges = false;
  {
    stk::diag::TimeBlock timer_(myRefinementTimer);
    didMakeAnyChanges = myRefinement.do_refinement(get_marker());
  }

  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Adapt: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));

  return didMakeAnyChanges;
}

bool KrinoRefinement::do_uniform_refinement(const int numUniformRefinementLevels)
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Uniform refinement: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  const bool didMakeAnyChanges = myRefinement.do_uniform_refinement(numUniformRefinementLevels);

  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Uniform refinement: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));

  return didMakeAnyChanges;
}

void KrinoRefinement::restore_after_restart()
{
  myRefinement.restore_after_restart();

  ParallelThrowAssert(myMeta.mesh_bulk_data().parallel(), check_face_and_edge_ownership(myMeta.mesh_bulk_data()));
  ParallelThrowAssert(myMeta.mesh_bulk_data().parallel(), check_face_and_edge_relations(myMeta.mesh_bulk_data()));
}

void check_leaf_children_have_parents_on_same_proc(const stk::ParallelMachine comm, const RefinementInterface & refinement)
{
  RequireEmptyErrorMsg(comm, refinement.locally_check_leaf_children_have_parents_on_same_proc(), "Leaf child without parent owned on same proc.");
}

}
