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
#include "Akri_DiagWriter.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_ReportHandler.hpp"
#include "Akri_TransitionElementEdgeMarker.hpp"

namespace krino {

void clear_refinement_marker(const RefinementInterface & refinement)
{
  FieldRef markerField = refinement.get_marker_field();
  stk::mesh::field_fill(static_cast<int>(Refinement::NOTHING), markerField, stk::mesh::selectField(markerField));
}

void mark_selected_elements_for_refinement(const RefinementInterface & refinement, const stk::mesh::Selector & selector)
{
  FieldRef markerField = refinement.get_marker_field();
  clear_refinement_marker(refinement);
  stk::mesh::field_fill(static_cast<int>(Refinement::REFINE), markerField, selector);
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
  FieldRef markerField = refinement.get_marker_field();
  for (auto && elem : elemsToRefine)
  {
    int * elemMarker = field_data<int>(markerField, elem);
    *elemMarker = Refinement::REFINE;
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
  const FieldRef markerField = refinement.get_marker_field();
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
      markerData[i] = indicatorData[i] > threshold_val ? Refinement::REFINE : Refinement::NOTHING;
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

template<class RefinementClass>
void check_leaf_children_have_parents_on_same_proc(const stk::mesh::BulkData & mesh, const RefinementClass * refinement)
{
  bool error = false;
  std::string msg = "Error: leaf child without parent owned by same processor";
  auto buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part());
  for(auto bucket : buckets)
  {
    for(auto && elem : *bucket)
    {
      if(refinement->is_child(elem) && !refinement->is_parent(elem))
      {
        auto parent = refinement->get_parent(elem);
        if(!mesh.is_valid(parent))
        {
          error = true;
          msg +=  " " + debug_entity(mesh, elem);
          break;
        }
        else if(!mesh.bucket(parent).owned())
        {
          error = true;
          msg +=  " " + debug_entity(mesh, elem);
          break;
        }
      }
    }
    if(error == true) break;
  }
  ParallelThrowRequireMsg(mesh.parallel(), !error, msg);
}

PerceptRefinement &
PerceptRefinement::get(const stk::mesh::MetaData & meta)
{
  PerceptRefinement * refinement = const_cast<PerceptRefinement*>(meta.get_attribute<PerceptRefinement>());
  ThrowRequireMsg(nullptr != refinement, "Refinement not found on MetaData.");
  return *refinement;
}

bool
PerceptRefinement::is_created(const stk::mesh::MetaData & meta)
{
  PerceptRefinement * refinement = const_cast<PerceptRefinement*>(meta.get_attribute<PerceptRefinement>());
  return refinement != nullptr;
}

PerceptRefinement &
PerceptRefinement::create(stk::mesh::MetaData & meta)
{
  PerceptRefinement * refinement = const_cast<PerceptRefinement*>(meta.get_attribute<PerceptRefinement>());
  ThrowRequireMsg(nullptr == refinement, "PerceptRefinement::create should be called only once per MetaData.");
  if (nullptr == refinement)
  {
    refinement = new PerceptRefinement(meta);
    meta.declare_attribute_with_delete<PerceptRefinement>(refinement);
  }
  return *refinement;
}

PerceptRefinement::PerceptRefinement(stk::mesh::MetaData & meta) :
    myMeta(meta)
{
  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  if (auxMeta.has_field(stk::topology::ELEMENT_RANK, "refine_level"))
    myRefinementLevelField = auxMeta.get_field(stk::topology::ELEMENT_RANK, "refine_level");
  const std::string transitionElementFieldName = (myMeta.spatial_dimension() == 2) ?
    "transition_element" : "transition_element_3";
  if (auxMeta.has_field(stk::topology::ELEMENT_RANK, transitionElementFieldName))
    myTransitionElementField = auxMeta.get_field(stk::topology::ELEMENT_RANK, transitionElementFieldName);

  myElementMarkerField = meta.declare_field<int>(stk::topology::ELEMENT_RANK, "REFINEMENT_ELEMENT_MARKER", 1);
  stk::mesh::put_field_on_mesh(myElementMarkerField.field(), meta.universal_part(), 1, 1, nullptr);
}

bool PerceptRefinement::is_child(const stk::mesh::Entity entity) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);
  const unsigned num_family_trees = mesh.num_connectivity(entity, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(entity, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity ft = family_trees[ifamily];
    const stk::mesh::Entity* family_tree_entities = mesh.begin(ft, entity_rank);
    if(family_tree_entities[0] != entity) return true;
  }
  return false;
}

bool PerceptRefinement::is_parent(const stk::mesh::Entity parent) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(parent);
  const unsigned num_family_trees = mesh.num_connectivity(parent, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(parent, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity* family_tree_entities = mesh.begin(family_trees[ifamily], entity_rank);
    const stk::mesh::Entity tree_parent = family_tree_entities[0]; // 0th entry in the family_tree is the parent
    if (parent == tree_parent) // I am the parent
    {
      return true;
    }
  }
  return false;
}

bool PerceptRefinement::is_parent_side(const stk::mesh::Entity side) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  return mesh.num_connectivity(side, stk::topology::CONSTRAINT_RANK) > 0; // Not a great answer but was working before and this is temporary
}

stk::mesh::Entity PerceptRefinement::get_parent(const stk::mesh::Entity entity) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);
  const unsigned num_family_trees = mesh.num_connectivity(entity, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(entity, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity ft = family_trees[ifamily];
    const stk::mesh::Entity* family_tree_entities = mesh.begin(ft, entity_rank);
    if(family_tree_entities[0] != entity) return family_tree_entities[0];
  }
  return stk::mesh::Entity();
}

std::pair<stk::mesh::EntityId,int> PerceptRefinement::get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const
{
  stk::mesh::Entity parent = get_parent(child);
  ThrowAssert(myMeta.mesh_bulk_data().is_valid(parent));
  return {myMeta.mesh_bulk_data().identifier(parent), myMeta.mesh_bulk_data().parallel_owner_rank(parent)};
}

static stk::mesh::Part & get_percept_refinement_inactive_part(const stk::mesh::MetaData & meta, stk::mesh::EntityRank rank)
{
  const std::string inactive_part_name = "refine_inactive_elements_part_"+std::to_string((int)rank);
  stk::mesh::Part* inactive_part = meta.get_part(inactive_part_name);
  ThrowRequireMsg(nullptr != inactive_part, "Inactive (parent) part not found: " << inactive_part_name);
  return *inactive_part;
}

static stk::mesh::Part & get_percept_refinement_active_part(const stk::mesh::MetaData & meta, stk::mesh::EntityRank rank)
{
  const std::string active_part_name = "refine_active_elements_part_"+std::to_string((int)rank);
  stk::mesh::Part* active_part = meta.get_part(active_part_name);
  ThrowRequireMsg(nullptr != active_part, "Active (child) part not found: " << active_part_name);
  return *active_part;
}

stk::mesh::Part & PerceptRefinement::parent_part() const
{
  return get_percept_refinement_inactive_part(myMeta, stk::topology::ELEMENT_RANK);
}

stk::mesh::Part & PerceptRefinement::child_part() const
{
  return get_percept_refinement_active_part(myMeta, stk::topology::ELEMENT_RANK);
}

unsigned PerceptRefinement::get_num_children(const stk::mesh::Entity parent) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(parent);
  const unsigned num_family_trees = mesh.num_connectivity(parent, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(parent, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity* family_tree_entities = mesh.begin(family_trees[ifamily], entity_rank);
    const stk::mesh::Entity tree_parent = family_tree_entities[0]; // 0th entry in the family_tree is the parent
    if (parent == tree_parent) // I am the parent
    {
      const unsigned num_family_tree_entities = mesh.num_connectivity(family_trees[ifamily], entity_rank);
      return num_family_tree_entities-1;
    }
  }
  return 0;
}

void PerceptRefinement::fill_children(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  children.clear();
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(parent);
  const unsigned num_family_trees = mesh.num_connectivity(parent, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(parent, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity* family_tree_entities = mesh.begin(family_trees[ifamily], entity_rank);
    const stk::mesh::Entity tree_parent = family_tree_entities[0]; // 0th entry in the family_tree is the parent
    if (parent == tree_parent) // I am the parent
    {
      const unsigned num_family_tree_entities = mesh.num_connectivity(family_trees[ifamily], entity_rank);
      for (unsigned ichild=1; ichild < num_family_tree_entities; ++ichild)
      {
        children.push_back(family_tree_entities[ichild]);
      }
    }
  }
}

void PerceptRefinement::fill_child_element_ids(const stk::mesh::Entity parent, std::vector<stk::mesh::EntityId> & childElementIds) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  childElementIds.clear();
  const unsigned num_family_trees = mesh.num_connectivity(parent, stk::topology::CONSTRAINT_RANK);
  const stk::mesh::Entity* family_trees = mesh.begin(parent, stk::topology::CONSTRAINT_RANK);
  for (unsigned ifamily=0; ifamily < num_family_trees; ++ifamily)
  {
    const stk::mesh::Entity* family_tree_entities = mesh.begin_elements(family_trees[ifamily]);
    const stk::mesh::Entity tree_parent = family_tree_entities[0]; // 0th entry in the family_tree is the parent
    if (parent == tree_parent) // I am the parent
    {
      const unsigned num_family_tree_entities = mesh.num_elements(family_trees[ifamily]);
      for (unsigned ichild=1; ichild < num_family_tree_entities; ++ichild)
      {
        childElementIds.push_back(mesh.identifier(family_tree_entities[ichild]));
      }
    }
  }
}

void PerceptRefinement::fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const
{
  dependents.clear();
  if(is_child(parent)) return;
  fill_all_children(*this, parent, dependents);
}

bool PerceptRefinement::has_rebalance_constraint(const stk::mesh::Entity entity) const
{
  return is_child(entity);
}

int PerceptRefinement::fully_refined_level(const stk::mesh::Entity elem) const
{
  const int refineLevel = *field_data<int>(myRefinementLevelField, elem);
  const int transitionElement = *field_data<int>(myTransitionElementField, elem);
  return refineLevel-transitionElement;
}

FieldRef PerceptRefinement::get_marker_field() const
{
  ThrowRequireMsg(myElementMarkerField.valid(), "PerceptRefinement created without setting the hadapt function, which is needed in get_marker_field().  Is this a unit test?");
  return myElementMarkerField;
}

bool PerceptRefinement::is_transition(const stk::mesh::Entity elem) const
{
  const int isTransitionElement = *field_data<int>(myTransitionElementField, elem);
  return isTransitionElement;
}

void PerceptRefinement::set_adaptive_refinement_function(const std::function<void(const std::string &, int)> & adaptiveRefinement)
{
  myAdaptiveRefinement = adaptiveRefinement;
}

void PerceptRefinement::set_uniform_refinement_function(const std::function<void(int)> & uniformRefinement)
{
  myUniformRefinement = uniformRefinement;
}

void PerceptRefinement::do_refinement(const int debugLevel)
{
  ThrowRequireMsg(myAdaptiveRefinement && myElementMarkerField.valid(), "PerceptRefinement created without calling set_adaptive_refinement_function, which is needed in do_adaptive_refinement().");
  myAdaptiveRefinement(myElementMarkerField.name(), debugLevel);
}

void PerceptRefinement::do_uniform_refinement(const int numUniformRefinementLevels)
{
  ThrowRequireMsg(myAdaptiveRefinement && myElementMarkerField.valid(), "PerceptRefinement created without calling set_uniform_refinement_function, which is needed in do_uniform_refinement().");
  myUniformRefinement(numUniformRefinementLevels);
}

KrinoRefinement &
KrinoRefinement::get(const stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  ThrowRequireMsg(nullptr != refinement, "Refinement not found on MetaData.");
  return *refinement;
}

KrinoRefinement &
KrinoRefinement::create(stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  ThrowRequireMsg(nullptr == refinement, "KrinoRefinement::create should be called only once per MetaData.");
  if (nullptr == refinement)
  {
    AuxMetaData & auxMeta = AuxMetaData::get(meta);
    refinement = new KrinoRefinement(meta, &auxMeta.active_part(), false/*auxMeta.get_force_64bit_flag()*/, auxMeta.get_assert_32bit_flag());
    meta.declare_attribute_with_delete<KrinoRefinement>(refinement);
  }
  return *refinement;
}

bool
KrinoRefinement::is_created(const stk::mesh::MetaData & meta)
{
  KrinoRefinement * refinement = const_cast<KrinoRefinement*>(meta.get_attribute<KrinoRefinement>());
  return refinement != nullptr;
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

KrinoRefinement::KrinoRefinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart, const bool force64Bit, const bool assert32Bit) :
    myMeta(meta),
    myRefinement(meta, activePart, force64Bit, assert32Bit)
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

void KrinoRefinement::fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const
{
  dependents.clear();
  std::vector<stk::mesh::Entity> children;
  fill_children(parent, children);
  for(auto && child : children)
  {
    if(child == stk::mesh::Entity()) continue;
    if(is_parent(child)) continue;
    dependents.push_back(child);
  }
}

bool KrinoRefinement::has_rebalance_constraint(const stk::mesh::Entity entity) const
{
  if(!is_parent(entity))
  {
    //if not a parent but is a child, must be leaf element, constrained
    if(is_child(entity)) return true;
    //if not a parent or child, must be completed unadapted element. No constraint
    else return false;
  }
  //if a parent, check if any children are parents or invalid elements (already moved)
  //if so, constrained. If not, not constrained
  std::vector<stk::mesh::Entity> children;
  fill_children(entity, children);
  for(auto && child : children)
  {
    if(child == stk::mesh::Entity())
      return true;
    if(is_parent(child))
      return true;
  }
  return false;
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

FieldRef KrinoRefinement::get_marker_field() const
{
  setup_marker();
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
  ThrowRequireMsg(myMarker, "Logic error.  Marker should have already been setup.");
  return *myMarker;
}

std::pair<unsigned,unsigned> KrinoRefinement::get_marked_element_counts() const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const FieldRef markerField = get_marker_field();
  const stk::mesh::Selector selector = stk::mesh::selectField(markerField) & AuxMetaData::get(myMeta).active_part() & myMeta.locally_owned_part() & !myRefinement.parent_part();

  unsigned numRefine = 0;
  unsigned numUnrefine = 0;

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const int * marker = field_data<int>(markerField, *bucket);
    const int length = bucket->size();

    for (int i = 0; i < length; ++i)
    {
      if (marker[i] == Refinement::COARSEN)
        ++numUnrefine;
      else if (marker[i] == Refinement::REFINE)
        ++numRefine;
    }
  }

  const std::array<unsigned, 2> localNum{numRefine, numUnrefine};
  std::array<unsigned, 2> globalNum{0,0};
  stk::all_reduce_sum(mesh.parallel(), localNum.data(), globalNum.data(), 2);

  return std::make_pair(globalNum[0], globalNum[1]);
}

void KrinoRefinement::do_refinement(const int debugLevel)
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  const auto markerCounts = get_marked_element_counts();
  const unsigned numRefine = markerCounts.first;
  const unsigned numUnrefine = markerCounts.second;
  krinolog << "Number of elements marked for refinement = " << numRefine << "\n";
  krinolog << "Number of elements marked for unrefinement = " << numUnrefine << stk::diag::dendl;

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Adapt: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  myRefinement.do_refinement(get_marker());

  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Adapt: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));

}

void KrinoRefinement::do_uniform_refinement(const int numUniformRefinementLevels)
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Uniform refinement: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  myRefinement.do_uniform_refinement(numUniformRefinementLevels);

  stk::mesh::comm_mesh_counts(mesh, counts);

  krinolog << "Uniform refinement: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << stk::diag::dendl;

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));
}

void KrinoRefinement::restore_after_restart()
{
  myRefinement.restore_after_restart();

  ParallelThrowAssert(myMeta.mesh_bulk_data().parallel(), check_face_and_edge_ownership(myMeta.mesh_bulk_data()));
  ParallelThrowAssert(myMeta.mesh_bulk_data().parallel(), check_face_and_edge_relations(myMeta.mesh_bulk_data()));
}

template void check_leaf_children_have_parents_on_same_proc(const stk::mesh::BulkData & mesh, const RefinementInterface * refinement);
template void check_leaf_children_have_parents_on_same_proc(const stk::mesh::BulkData & mesh, const Refinement * refinement);
}
