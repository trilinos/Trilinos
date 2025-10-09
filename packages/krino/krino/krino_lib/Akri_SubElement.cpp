// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourUtils.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_Utility.hpp>
#include <stk_math/StkVector.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

#include <cmath>
#include <memory>

namespace krino{

bool
SubElementNode::on_common_edge(const SubElementNode * other) const
{
  NodeSet ancestors;
  get_ancestors(ancestors);
  if (ancestors.size() <= 2)
  {
    other->get_ancestors(ancestors);
  }
  if (ancestors.size() <= 2)
  {
    return true;
  }
  return false;
}

SubElementMidSideNode::SubElementMidSideNode( const Mesh_Element * owner,
    const SubElementNode *parent1,
    const SubElementNode *parent2)
  : SubElementNode(owner),
    my_is_mesh_node(false),
    my_parent1(parent1),
    my_parent2(parent2)
{
  // fill base class data
  my_cached_owner_coords = compute_owner_coords( owner );
  my_global_coords = 0.5*(my_parent1->coordinates()) + 0.5*(my_parent2->coordinates());
}

SubElementMidSideNode::SubElementMidSideNode( const Mesh_Element * owner,
    const SubElementNode *parent1,
    const SubElementNode *parent2,
    stk::mesh::Entity meshNode,
    stk::mesh::EntityId meshNodeId)
  : SubElementMidSideNode(owner, parent1, parent2)
{
  my_is_mesh_node = true;
  set_entity(meshNode, meshNodeId);
}

SubElementChildNode::SubElementChildNode( const Mesh_Element * in_owner,
    const NodeVec & parents,
    const std::vector<double> & weights )
  : SubElementNode(in_owner),
    my_parents(parents),
    my_weights(weights)
{
  // fill base class data
  my_cached_owner_coords = compute_owner_coords( in_owner );
  my_global_coords = in_owner->coordinates( my_cached_owner_coords );
}

SubElementMeshNode::SubElementMeshNode( const Mesh_Element * in_owner,
    stk::mesh::Entity nodeEntity,
    stk::mesh::EntityId nodeEntityId,
    const stk::math::Vector3d & in_owner_coords,
    const stk::math::Vector3d & in_global_coords )
  : SubElementNode(in_owner)
{
  // fill base class data
  set_entity(nodeEntity, nodeEntityId);
  my_cached_owner_coords = in_owner_coords;
  my_global_coords = in_global_coords;
}

void
SubElementNode::fill_parent_entity_pointers(std::vector<stk::mesh::Entity*> & parentEntities) const
{
  NodeVec parents = get_parents();

  const unsigned numParents = parents.size();
  parentEntities.resize(numParents);

  for (unsigned i=0; i<numParents; ++i)
    parentEntities[i] = &parents[i]->entity();
}

static bool is_on_multiple_blocks(const stk::mesh::BulkData& mesh, stk::mesh::Entity node)
{
  bool foundVolumePart = false;
  for (auto && part : mesh.bucket(node).supersets())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
       !stk::mesh::is_auto_declared_part(*part) &&
       part->subsets().empty())
    {
      if (foundVolumePart) return true;
      foundVolumePart = true;
    }
  }
  return false;
}

bool SubElementChildNode::needs_to_be_ale_prolonged(const CDMesh & mesh) const
{
  if (mesh.get_prolongation_model() == INTERPOLATION)
    return false;

  const bool is_initial_mesh = !mesh.was_mesh_previously_decomposed();
  if (is_initial_mesh)
    return false;

  if (is_on_multiple_blocks(mesh.stk_bulk(), entity()))
    return true;

  // relatively unusual case of an edge node that is not on a block-block boundary.
  // this is currently handled by using interpolation.  This possibly needs further
  // testing/development to treat these like mesh nodes where we see if they have
  // changed phase.
  return false;
}

stk::math::Vector3d SubElementChildNode::compute_owner_coords( const Mesh_Element * in_owner ) const
{
  stk::math::Vector3d calcOwnerCoords{stk::math::Vector3d::ZERO};
  STK_ThrowAssert(my_parents.size() == my_weights.size());
  for (size_t i=0; i<my_parents.size(); ++i)
    calcOwnerCoords += my_weights[i]*my_parents[i]->owner_coords(in_owner);
  return calcOwnerCoords;
}

void
SubElementChildNode::prolongate_fields(const CDMesh & mesh) const
{
  for (auto && parent : my_parents)
    if (!parent->is_prolonged())
      parent->prolongate_fields(mesh);

  if (my_is_prolonged_flag) return;
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "SubElementEdgeNode::prolongate_fields for node#" << entityId() << "\n";
  my_is_prolonged_flag = true;

  ProlongationQuery prolongQuery;
  if (needs_to_be_ale_prolonged(mesh))
    prolongQuery = mesh.find_prolongation_node(*this);
  const ProlongationPointData * prolong_node = prolongQuery.get_prolongation_point_data();

  prolong_cdfem_displacements(mesh, prolong_node);

  prolong_zeroed_fields(mesh, nullptr);

  prolong_ale_fields(mesh, prolong_node);

  prolong_edge_interpolation_fields(mesh.get_edge_interpolation_fields());

  prolong_interpolation_fields(mesh);
}

bool SubElementMidSideNode::is_mesh_node_that_needs_to_be_prolonged(const CDMesh & mesh) const
{

  STK_ThrowRequire(my_is_mesh_node);

  const SubElementMeshNode * parent1 = dynamic_cast<const SubElementMeshNode *>(my_parent1);
  const SubElementMeshNode * parent2 = dynamic_cast<const SubElementMeshNode *>(my_parent2);
  const int num_ale_prolonged_parents = (parent1->needs_to_be_ale_prolonged(mesh) ? 1 : 0) + (parent2->needs_to_be_ale_prolonged(mesh) ? 1 : 0);

  if (num_ale_prolonged_parents == 0) return false;
  if (num_ale_prolonged_parents == 2) return true;

  // 1 parent node needed to be ALE prolonged and this node is active (so the edge is not cut).
  // This means the interface was cutting this edge, but now is not -> prolong OR
  // the interface is passing through one of the parents of this uncut edge -> do not prolong.

  const bool have_or_did_have_interface = my_cached_owner->have_interface() || mesh.fetch_prolong_element(my_cached_owner->entityId())->have_subelements();

  return have_or_did_have_interface && nullptr == mesh.fetch_prolong_node(entityId());
}

void
SubElementMidSideNode::prolongate_fields(const CDMesh & mesh) const
{
  if (!my_parent1->is_prolonged())
  {
    my_parent1->prolongate_fields(mesh);
  }
  if (!my_parent2->is_prolonged())
  {
    my_parent2->prolongate_fields(mesh);
  }
  if (my_is_prolonged_flag) return;
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "SubElementMidSideNode::prolongate_fields for node#" << entityId() << "\n";
  my_is_prolonged_flag = true;

  const bool needs_to_be_prolonged = !my_is_mesh_node || is_mesh_node_that_needs_to_be_prolonged(mesh);
  if (needs_to_be_prolonged)
  {
    // Note: CDFEM displacement is not present on midside nodes
    prolong_zeroed_fields(mesh, nullptr);

    prolong_edge_interpolation_fields(mesh.get_edge_interpolation_fields());

    prolong_interpolation_fields(mesh);

    prolong_ale_fields(mesh);
  }
}

void
SubElementMidSideNode::prolong_interpolation_fields(const CDMesh & mesh) const
{
  const ProlongationElementData * interpolationElem = nullptr;
  stk::math::Vector3d interp_elem_p_coords;
  const ProlongationElementData * prolongElem =  mesh.fetch_prolong_element(my_cached_owner->entityId());
  STK_ThrowRequire(prolongElem);
  prolongElem->find_subelement_and_parametric_coordinates_at_point(coordinates(), interpolationElem, interp_elem_p_coords);

  for(auto && field : mesh.get_interpolation_fields())
  {
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    interpolationElem->evaluate_prolongation_field(mesh.get_cdfem_support(), field, field_length, interp_elem_p_coords, val);
  }
}

static void prolongate_edge_interpolation_fields_for_node(const FieldSet & edgeInterpFields,
    const NodeVec & parents,
    const std::vector<double> & weights,
    const stk::mesh::Entity node)
{
  for(auto && field : edgeInterpFields)
  {
    double * val = field_data<double>(field, node);
    if (nullptr == val) continue;

    const unsigned fieldLength = field.length();
    for (unsigned i=0; i<fieldLength; ++i)
      val[i] = 0.0;

    const size_t numParents = parents.size();
    for(size_t iParent=0; iParent<numParents; ++iParent)
    {
      double * parentVal = field_data<double>(field, parents[iParent]->entity());
      STK_ThrowRequireMsg(parentVal, "All parents must have edge interpolation field if child has field.");

      for (unsigned i=0; i<fieldLength; ++i)
        val[i] += weights[iParent] * parentVal[i];
    }
  }
}

void
SubElementChildNode::prolong_edge_interpolation_fields(const FieldSet & edgeInterpFields) const
{
  prolongate_edge_interpolation_fields_for_node(edgeInterpFields, get_parents(), get_parent_weights(), my_entity);
}

void
SubElementMidSideNode::prolong_edge_interpolation_fields(const FieldSet & edgeInterpFields) const
{
  prolongate_edge_interpolation_fields_for_node(edgeInterpFields, get_parents(), get_parent_weights(), my_entity);
}

void
SubElementMidSideNode::prolong_ale_fields(const CDMesh & mesh) const
{
  // simply average parent nodes
  for(auto && field : mesh.get_ale_prolongation_fields())
  {
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (nullptr != val)
    {
      double * val1 = field_data<double>(field, my_parent1->entity());
      double * val2 = field_data<double>(field, my_parent2->entity());
      STK_ThrowRequire(val1 && val2);
      for (unsigned i=0; i<field_length; ++i)
      {
        val[i] = 0.5*val1[i] + 0.5*val2[i];
      }
    }
  }
}

void
SubElementSteinerNode::prolongate_fields(const CDMesh & mesh) const
{
  for (auto && parent : get_parents())
  {
    if (!parent->is_prolonged())
    {
      parent->prolongate_fields(mesh);
    }
  }
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "SubElementInternalNode::prolongate_fields for node#" << entityId() << "\n";
  my_is_prolonged_flag = true;

  const stk::mesh::BulkData& stk_mesh = mesh.stk_bulk();
  const stk::mesh::FieldVector & all_fields = stk_mesh.mesh_meta_data().get_fields();
  for ( stk::mesh::FieldVector::const_iterator it = all_fields.begin(); it != all_fields.end() ; ++it )
  {
    const FieldRef field = (const FieldRef)(**it);

    // Do not try to prolong non-real variables
    if( field.entity_rank()!=stk::topology::NODE_RANK || !field.type_is<double>() ) continue;

    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    for (unsigned i=0; i<field_length; ++i)
    {
      val[i] = 0.0;
    }

    bool parent_error = false;
    double tot_wt = 0.0;
    for (unsigned p=0; p<get_parents().size(); ++p)
    {
      const double * parent_val = field_data<double>(field, get_parents()[p]->entity());
      if (NULL == parent_val)
      {
        parent_error = true;
      }
      else
      {
        tot_wt += get_parent_weights()[p];
        for (unsigned i=0; i<field_length; ++i)
        {
          val[i] += get_parent_weights()[p] * parent_val[i];
        }
      }
    }

    if (parent_error)
    {
      for (unsigned i=0; i<field_length; ++i)
      {
        val[i] /= tot_wt;
      }
      krinolog << "Error prolongating internal node field for node " << debug_entity(stk_mesh, entity()) << stk::diag::dendl;
      for (unsigned p=0; p<get_parents().size(); ++p)
      {
        krinolog << "  Parent " << p << ": " << debug_entity(stk_mesh, get_parents()[p]->entity()) << stk::diag::dendl;
      }
    }
  }
}

bool on_interface_or_io_parts_have_changed(const CDMesh & mesh, const Phase_Support & phaseSupport, stk::mesh::Entity node, const ProlongationNodeData & oldProlongNode)
{
  const auto newParts = PartAndFieldCollections::determine_io_parts(mesh.stk_bulk().bucket(node));
  const auto & oldParts = mesh.get_prolong_part_and_field_collections().get_parts(oldProlongNode.get_part_collection_id());
  if (newParts != oldParts)
    return true;
  for (auto && partOrdinal : newParts)
    if (phaseSupport.is_interface(mesh.stk_meta().get_part(partOrdinal)))
      return true;
  return false;
}

bool SubElementMeshNode::needs_to_be_ale_prolonged(const CDMesh & mesh) const
{
  const ProlongationNodeData * old_prolong_node = NULL;
  old_prolong_node = mesh.fetch_prolong_node(entityId());
  const bool is_initial_mesh = !mesh.was_mesh_previously_decomposed();
  return !is_initial_mesh && nullptr != old_prolong_node && on_interface_or_io_parts_have_changed(mesh, mesh.get_phase_support(), entity(), *old_prolong_node);
}

void
SubElementMeshNode::prolongate_fields(const CDMesh & mesh) const
{

  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "SubElementMeshNode::prolongate_fields for node#" << entityId() << "\n";
  my_is_prolonged_flag = true;

  const ProlongationNodeData * old_prolong_node = nullptr;
  old_prolong_node = mesh.fetch_prolong_node(entityId());

  const bool needsToBeALEProlonged = needs_to_be_ale_prolonged(mesh);
  ProlongationQuery prolongQuery;
  if (mesh.get_prolongation_model() != INTERPOLATION && needsToBeALEProlonged)
    prolongQuery = mesh.find_prolongation_node(*this);
  const ProlongationPointData * prolong_data = prolongQuery.get_prolongation_point_data();

  if( !old_prolong_node && !prolong_data )
  {
    return;
  }

  prolong_cdfem_displacements(mesh, prolong_data, false);

  const ProlongationNodeData * nodeToExamineForPreExistingField = needsToBeALEProlonged ? nullptr : old_prolong_node;
  prolong_zeroed_fields(mesh, nodeToExamineForPreExistingField);

  prolong_ale_fields(mesh, prolong_data, old_prolong_node);

  prolong_interpolation_fields(mesh, old_prolong_node);
}

void
SubElementNode::prolong_zeroed_fields(const CDMesh & mesh, const ProlongationNodeData * nodeToExamineForPreExistingField) const
{
  const FieldSet & zeroed_fields = mesh.get_zeroed_fields();
  for(auto&& field : zeroed_fields)
  {
    if (!nodeToExamineForPreExistingField || !nodeToExamineForPreExistingField->get_field_data(field)) // If this node existed before and had this field, leave it alone
    {
      double * val = field_data<double>(field, my_entity);
      if (nullptr != val) std::fill(val, val+field.length(), 0.);
    }
  }
}

void SubElementNode::prolong_cdfem_displacements(const CDMesh & mesh,
    const ProlongationPointData * prolong_data,
    const bool zero_if_no_prolong_data) const
{
  const FieldRef field = mesh.get_cdfem_displacements_field();
  if( !field.valid()) return;

  const unsigned field_length = field.length();
  for ( unsigned is = 0; is < field.number_of_states(); ++is )
  {
    const stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(is);
    const FieldRef state_field = field.field_state(state);

    double * val = field_data<double>(state_field, my_entity);
    if(val == NULL) continue;

    if(!prolong_data)
    {
      if (zero_if_no_prolong_data)
      {
        std::fill(val, val+field_length, 0);
      }
    }
    else
    {
      const double * prolong_field = prolong_data->get_field_data(state_field);
      STK_ThrowRequire(NULL != prolong_field);
      std::copy(prolong_field, prolong_field+field_length, val);

      if (state == stk::mesh::StateNew)
      {
        const stk::math::Vector3d & coords = coordinates();
        const stk::math::Vector3d & old_coords = prolong_data->get_previous_coordinates();
        for (unsigned i=0; i<field_length; ++i)
        {
          val[i] += coords[i] - old_coords[i];
        }
      }
    }
  }
}

const SubElementNode *
SubElementNode::common_child( const NodeVec & parents )
{
  const size_t numParents = parents.size();
  STK_ThrowAssert(numParents > 0);

  for (auto && child : parents[0]->my_children)
  {
    if (child->get_num_parents() == numParents)
    {
      bool childOfAllParents = true;
      for (size_t iParent=1; iParent<numParents; ++iParent)
      {
        if (!parents[iParent]->have_child(child))
        {
          childOfAllParents = false;
          break;
        }
      }
      if (childOfAllParents)
        return child;
    }
  }
  return nullptr;
}

bool SubElementNode::have_child() const
{
  return !my_children.empty();
}

bool SubElementNode::have_child(const SubElementNode* child) const
{
  return (std::find(my_children.begin(), my_children.end(), child) != my_children.end());
}

void
SubElementChildNode::prolong_ale_fields(const CDMesh & mesh, const ProlongationPointData * prolong_data) const
{
  const ProlongationElementData * interpolationElem = nullptr;
  stk::math::Vector3d interp_elem_p_coords;
  const FieldSet & ale_prolongation_fields = mesh.get_ale_prolongation_fields();
  for(FieldSet::const_iterator it = ale_prolongation_fields.begin(); it != ale_prolongation_fields.end(); ++it)
  {
    const FieldRef field = *it;
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    if(prolong_data)
    {
      // this node has changed phase
      // prolong based on prolong_node
      const double * prolong_field = prolong_data->get_field_data(field);
      // We cannot yet handle the case where a prolongation field is not defined on the prolongation node that
      // was found. Throw here to avoid the possibility of not prolonging a prolongation field that then
      // has its time derivative screwed up because it has a mesh velocity associated with it.
      // We think this should only occur in problems with multiple different level sets.
      STK_ThrowRequire(prolong_field);
      std::copy(prolong_field, prolong_field+field_length, val);
    }
    else
    {
      if (nullptr == interpolationElem)
      {
        const ProlongationElementData * prolongElem = mesh.fetch_prolong_element(my_cached_owner->entityId());
        STK_ThrowRequire(prolongElem);
        prolongElem->find_subelement_and_parametric_coordinates_at_point(coordinates(), interpolationElem, interp_elem_p_coords);
      }

      interpolationElem->evaluate_prolongation_field(mesh.get_cdfem_support(), field, field_length, interp_elem_p_coords, val);
    }
  }
}

void
SubElementMeshNode::prolong_ale_fields(const CDMesh & mesh,
    const ProlongationPointData * prolong_data,
    const ProlongationNodeData * old_node) const
{
  const FieldSet & ale_prolongation_fields = mesh.get_ale_prolongation_fields();
  for(FieldSet::const_iterator it = ale_prolongation_fields.begin(); it != ale_prolongation_fields.end(); ++it)
  {
    const FieldRef field = *it;
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    if(prolong_data)
    {
      // this node has changed phase
      // prolong based on prolong_node
      const double * prolong_field = prolong_data->get_field_data(field);
      // We cannot yet handle the case where a prolongation field is not defined on the prolongation node that
      // was found. Throw here to avoid the possibility of not prolonging a prolongation field that then
      // has its time derivative screwed up because it has a mesh velocity associated with it.
      // We think this should only occur in problems with multiple different level sets.
      STK_ThrowRequire(prolong_field);
      std::copy(prolong_field, prolong_field+field_length, val);
    }
    else if(old_node)
    {
      const double * old_field = old_node->get_field_data(field);
      if(!old_field)
      {
        const FieldRef initial_field = mesh.get_cdfem_support().get_initial_prolongation_field( field );
        if(initial_field.valid())
        {
          const double * initial_data = old_node->get_field_data(initial_field);
          std::copy(initial_data, initial_data+field_length, val);
        }
        else
        {
          std::fill(val, val+field_length, 0.);
        }
      }
    }
  }
}

void
SubElementChildNode::prolong_interpolation_fields(const CDMesh & mesh) const
{
  const ProlongationElementData * interpolationElem = nullptr;
  stk::math::Vector3d interpElemParamCoords;
  const ProlongationElementData * prolongElem =  mesh.fetch_prolong_element(my_cached_owner->entityId());
  STK_ThrowRequire(prolongElem);
  prolongElem->find_subelement_and_parametric_coordinates_at_point(coordinates(), interpolationElem, interpElemParamCoords);

  const FieldSet & interpolation_fields = mesh.get_interpolation_fields();
  for(FieldSet::const_iterator it = interpolation_fields.begin(); it != interpolation_fields.end(); ++it)
  {
    const FieldRef field = *it;
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    interpolationElem->evaluate_prolongation_field(mesh.get_cdfem_support(), field, field_length, interpElemParamCoords, val);
  }
}

void
SubElementMeshNode::prolong_interpolation_fields(const CDMesh & mesh, const ProlongationNodeData * old_node) const
{
  const FieldSet & interpolation_fields = mesh.get_interpolation_fields();
  for(FieldSet::const_iterator it = interpolation_fields.begin(); it != interpolation_fields.end(); ++it)
  {
    const FieldRef field = *it;
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, my_entity);
    if (NULL == val) continue;

    if(old_node)
    {
      const double * old_field = old_node->get_field_data(field);
      if(!old_field)
      {
        const FieldRef initial_field = mesh.get_cdfem_support().get_initial_prolongation_field( field );
        if(initial_field.valid())
        {
          const double * initial_data = old_node->get_field_data(initial_field);
          std::copy(initial_data, initial_data+field_length, val);
        }
        else
        {
          std::fill(val, val+field_length, 0.);
        }
      }
    }
  }
}

void removeParts(stk::mesh::PartVector & parts, const stk::mesh::PartVector & parts_to_remove)
{
  stk::mesh::PartVector::iterator parts_begin = parts.begin();
  stk::mesh::PartVector::iterator begin_parts_to_erase = parts.end();
  for(stk::mesh::PartVector::const_iterator it = parts_to_remove.begin(); it != parts_to_remove.end(); ++it)
  {
    begin_parts_to_erase = std::remove(parts_begin, begin_parts_to_erase, *it);
  }
  parts.erase(begin_parts_to_erase, parts.end());
}

std::vector<unsigned> SubElementNode::prolongation_node_fields(const CDMesh & mesh) const
{
  const FieldSet & ale_prolongation_fields = mesh.get_ale_prolongation_fields();

  std::vector<unsigned> ale_prolongation_fields_on_node;
  ale_prolongation_fields_on_node.reserve(ale_prolongation_fields.size());
  for(auto && field : ale_prolongation_fields)
    if (nullptr != field_data<double>(field, my_entity))
      ale_prolongation_fields_on_node.push_back(field.field().mesh_meta_data_ordinal());

  std::sort(ale_prolongation_fields_on_node.begin(), ale_prolongation_fields_on_node.end());
  return ale_prolongation_fields_on_node;
}

static bool float_less(double a, double b)
{
  return static_cast<float>(a) < static_cast<float>(b);
}

bool SubElementNode::higher_priority_by_score_then_ancestry(const SubElementNode & a, const SubElementNode & b, const bool globalIDsAreParallelConsistent)
{
  // higher score wins
  if (float_less(b.get_node_score(), a.get_node_score())) return true;
  if (float_less(a.get_node_score(), b.get_node_score())) return false;

  if (globalIDsAreParallelConsistent)
    return SubElementNodeAncestry::compare(a.get_ancestry(), b.get_ancestry(), SubElementNode::less_by_entity_id);
  return SubElementNodeAncestry::compare(a.get_ancestry(), b.get_ancestry(), SubElementNode::less_by_coordinates_then_by_entity_id);
}

bool SubElementNode::less_by_entity_id(const SubElementNode & a, const SubElementNode & b)
{
  // lower id wins (this can be an issue because it is sensitive to the global ids provided by percept, which are dependent on the number of procs)
  return a.entityId() < b.entityId();
}

bool SubElementNode::less_by_coordinates_then_by_entity_id(const SubElementNode & a, const SubElementNode & b)
{
  const stk::math::Vector3d & aCoord = a.coordinates();
  const stk::math::Vector3d & bCoord = b.coordinates();
  if (float_less(aCoord[0], bCoord[0])) return true;
  if (float_less(bCoord[0], aCoord[0])) return false;
  if (float_less(aCoord[1], bCoord[1])) return true;
  if (float_less(bCoord[1], aCoord[1])) return false;
  if (float_less(aCoord[2], bCoord[2])) return true;
  if (float_less(bCoord[2], aCoord[2])) return false;

  // as a last resort, lower id wins (this can be an issue because it is sensitive to the global ids provided by percept, which are dependent on the number of procs)
  return a.entityId() < b.entityId();
}

bool SubElementNode::captures_intersection_point_domains(const std::vector<int> & intersectionPointDomains) const
{
  return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(my_sorted_node_domains, intersectionPointDomains);
}

bool SubElementNode::captures_interface(const InterfaceID & interface) const
{
  if (interface.first_ls() == interface.second_ls())
    return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(my_sorted_node_domains, {interface.first_ls()});
  return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(my_sorted_node_domains, {interface.first_ls(),interface.second_ls()});
}

void SubElementNode::insert_node_domains(const std::vector<int> & domainsToAdd) const
{
  my_sorted_node_domains.insert(my_sorted_node_domains.end(), domainsToAdd.begin(), domainsToAdd.end());
  stk::util::sort_and_unique(my_sorted_node_domains);
}

void
SubElementNode::get_ancestors(NodeSet & ancestors) const
{
  if (is_mesh_node())
  {
    ancestors.insert(this);
    return;
  }
  const NodeVec parents = get_parents();
  for(auto&& parent : parents)
  {
    parent->get_ancestors(ancestors);
  }
}

SubElementNodeAncestry
SubElementNode::get_ancestry() const
{
  return SubElementNodeAncestry(this);
}

void
SubElementNode::build_stencil(std::map<const SubElementNode *, double> & stencil, const double self_weight) const
{
  if (is_mesh_node())
  {
    stencil[this] += self_weight;
    return;
  }
  const NodeVec parents = get_parents();
  const std::vector<double> parent_weights = get_parent_weights();
  for(unsigned i=0; i<parents.size(); ++i)
  {
    parents[i]->build_stencil(stencil, self_weight * parent_weights[i]);
  }
}

void
SubElementNode::build_constraint_stencil(const FieldRef field, std::vector<stk::mesh::Entity> & entities, std::vector<double> & weights) const
{
  STK_ThrowRequire(!is_mesh_node());
  static const double wt_min = 1.e-9;
  typedef std::tuple<stk::mesh::EntityId, stk::mesh::Entity, double> EntityAndWeight;
  std::vector<EntityAndWeight> entitiesAndWeights;

  const MasterElement* master_elem = my_cached_owner->get_evaluation_master_element(field);

  const unsigned npe = master_elem->get_topology().num_nodes();
  std::vector<double> shapefcn (npe,0.);
  master_elem->shape_fcn(1, my_cached_owner_coords.data(), shapefcn.data());

  const auto & nodes = my_cached_owner->get_nodes();
  for (unsigned n=0; n<npe; ++n)
  {
    if (shapefcn[n] > wt_min)
    {
      entitiesAndWeights.push_back(std::make_tuple(nodes[n]->entityId(), nodes[n]->entity(), shapefcn[n]));
    }
  }

  std::sort(entitiesAndWeights.begin(), entitiesAndWeights.end(), [](const EntityAndWeight & a, const EntityAndWeight & b) { return std::get<0>(a) > std::get<0>(b); });

  entities.clear();
  weights.clear();

  entities.push_back(entity());
  weights.push_back(-1.0);

  for (auto && entityAndWeight : entitiesAndWeights)
  {
    entities.push_back(std::get<1>(entityAndWeight));
    weights.push_back(std::get<2>(entityAndWeight));
  }
}

// Notes on coordinate systems.
// (1) real coordinates
//     This can be obtained using Element::coordinates( Vector3d of owner coordinates ).
// (2) owner coordinates
//     This is the parametric coordinates of the finite element.
//     This can be obtained using SubElement::owner_coordinates( Vector3d of subelement coordinates ).
//     This is the type of coordinates stored in SubElement::my_coords, no matter what level of subelement.
//     In places where both owner coordinates and subelement coordinates are used, this is kept in the var owner_coords.
// (3) subelement coordinates
//     This is the parametric coordinates within the subelement.
//     In places where both owner coordinates and subelement coordinates are used, this is kept in the var local_coords.

SubElement::SubElement( const stk::topology topo,
               const NodeVec & nodes,
               const std::vector<int> & side_ids,
               const Mesh_Element * owner)
    : ElementObj( topo, nodes),
      my_parent_side_ids( side_ids ),
      my_owner( owner )
{
  STK_ThrowAssert( nodes.size() == topology().num_nodes() );
  STK_ThrowAssert( my_parent_side_ids.size() == topology().num_sides() );

  set_permutation();
}

void
SubElement::set_permutation()
{
  // For true subelements (not just coincident with the owning mesh element), permute the element
  // nodes and sides in a consistent way.  This will produce consistent node and side ordering
  // from decomposition to decomposition. In this way, repeated decompositions will produce identical results.

  bool coincident_with_owner = true;
  for (size_t n=0; n<my_nodes.size(); ++n)
  {
    if (my_nodes[n] != my_owner->get_nodes()[n])
    {
       coincident_with_owner = false;
       break;
    }
  }
  if (coincident_with_owner) return;

  std::vector<SubElementNodeAncestry> node_ancestries;
  node_ancestries.reserve(my_nodes.size());
  for (auto&& node : my_nodes)
  {
    node_ancestries.emplace_back(node);
  }
  const unsigned permutation_index = topology().lexicographical_smallest_permutation(node_ancestries.data(), true); // only consider positive permutations (true means this)

  if (permutation_index == stk::mesh::DEFAULT_PERMUTATION) return;

  // permute nodes
  NodeVec permuted_nodes(my_nodes.size());
  topology().permutation_nodes(my_nodes.data(), permutation_index, permuted_nodes.data());
  my_nodes = permuted_nodes;

  // permute sides
  std::vector<unsigned> side_permutation = get_side_permutation(topology(), static_cast<stk::mesh::Permutation>(permutation_index));

  std::vector<int> permuted_parent_side_ids(my_parent_side_ids.size());
  for (unsigned iside=0; iside<my_parent_side_ids.size(); ++iside)
  {
    permuted_parent_side_ids[iside] = my_parent_side_ids[side_permutation[iside]];
  }
  my_parent_side_ids = permuted_parent_side_ids;
}

bool
SubElement::check_entity_nodes(const stk::mesh::BulkData & stkMesh) const
{
  // Before using a consistent permutation, this could be a problem
  const stk::mesh::Entity * elem_nodes = stkMesh.begin_nodes(entity());

  for (unsigned inode=0; inode<my_nodes.size(); ++inode)
  {
    if (elem_nodes[inode] != my_nodes[inode]->entity()) return false;
  }
  return true;
}

void
SubElement::get_owner_coord_transform(double * dOwnerdSub) const
{
//  const NodeVec & owner_nodes = my_owner->get_nodes();
//  for ( int n = 0; n < my_num_nodes; n++ )
//  {
//    krinolog << "node " << n << ", owner node phys_coords = " << owner_nodes[n]->coordinates()[0] << "," << owner_nodes[n]->coordinates()[1] << stk::diag::dendl;
//  }

  std::vector<stk::math::Vector3d> owner_coords;
  fill_node_owner_coords(my_owner, owner_coords);

  // Hard coded for linear simplex elements with constant transformations
  const unsigned nnodes = num_nodes();
  STK_ThrowAssert(my_master_elem.num_intg_pts() == nnodes);
  const double * d_shape = my_master_elem.shape_fcn_deriv();

  const int dim = spatial_dim();
  for ( int i = 0; i < dim; i++ )
  {
    for ( int j = 0; j < dim; j++ )
    {
      double & dOwnerdSub_ij = dOwnerdSub[i*dim + j];

      dOwnerdSub_ij = 0.0;
      for ( unsigned n = 0; n < nnodes; n++ )
      {
        dOwnerdSub_ij += owner_coords[n][i] * d_shape[n*dim + j];
      }

      //krinolog << "dOwnerdSub[" << i << "][" << j << "] = " << dOwnerdSub_ij << stk::diag::dendl;
    }
  }
}

void
SubElement::determine_decomposed_elem_phase(const std::vector<Surface_Identifier> & surfaceIDs)
{
  if(have_subelements())
  {
    for(auto && subelem : my_subelements)
    {
      subelem->determine_decomposed_elem_phase(surfaceIDs);
    }
    // Phase for SubElement with subelements is left empty
    return;
  }

  if(!my_owner->have_interface())
  {
    set_phase(my_owner->get_phase());
  }
  else
  {
    const PhaseTag & startPhase = my_phase.empty() ? my_owner->get_phase() : my_phase;
    my_phase = update_phase(surfaceIDs, startPhase, my_owner->get_sorted_cutting_interfaces(), myInterfaceSigns);
  }

  if(krinolog.shouldPrint(LOG_DEBUG))
  {
    krinolog << "SubElement with nodes ";
    for (auto && node : my_nodes)
      krinolog << node->get_ancestry() << " ";
    krinolog << "has phase " << my_phase << "\n";
  }
}

std::vector<int> SubElement::subelement_interface_signs(const InterfaceID interface, const int sign) const
{
  if (have_interface(interface) && sign != 0)
  {
    std::vector<int> interfaceSigns = myInterfaceSigns;
    interfaceSigns[my_owner->get_interface_index(interface)] = sign;
    return interfaceSigns;
  }
  return myInterfaceSigns;
}

void SubElement::initialize_interface_signs()
{
  myInterfaceSigns.assign(my_owner->get_num_interfaces(), 0);
}

void SubElement::set_interface_signs(const std::vector<int> & interfaceSigns)
{
  myInterfaceSigns = interfaceSigns;
}

void SubElement::update_interface_signs(const InterfaceID interface, const int sign)
{
  set_interface_signs(subelement_interface_signs(interface, sign));
}

double
SubElement::relative_volume() const
{
  // This is a relative volume compared to the owner volume.
  // Actually this is a relative volume if the "parametric" volume of the element is unity.
  // Otherwise, it is off by a factor.
  const int nelem = 1;
  const int dim   = spatial_dim();
  const int nint  = my_master_elem.num_intg_pts();
  std::vector<double> coords(my_nodes.size() * dim, 0.);
  std::vector<double> det_J(nint, 0.);
  double error = 0.;

  // integration weights
  const double * intg_weights = my_master_elem.intg_weights();

  // load coords
  int count = 0;
  for ( auto && node : my_nodes )
  {
    const stk::math::Vector3d & owner_coords = node->owner_coords(my_owner);
    for ( int j = 0; j < dim; j++ )
      {
        coords[count++] = owner_coords[j];
      }
  }

  // determinant at integration points
  my_master_elem.determinant( dim, nelem, coords.data(), det_J.data(), &error );

  double elem_volume = 0.;
  for ( int ip = 0; ip < nint; ip++ )
    {
      elem_volume += det_J[ip] * intg_weights[ip];
    }

  return elem_volume;
}

double
SubElement::maximum_relative_angle() const
{
  // Find the maximum angle formed at the vertices in parametric coordinates.
  // These are obviously differentt than that maximum angle in physical coordinates due to
  // the shape of the owning element.
  double max_angle = 0;

  const stk::topology topol = topology();
  const unsigned num_edges = topol.num_edges();
  for ( unsigned edge0 = 0; edge0 < num_edges; edge0++ )
  {
    const unsigned * lnn0 = get_edge_node_ordinals(topol, edge0);
    STK_ThrowAssert(
        2 == topol.edge_topology(edge0).num_nodes() || 3 == topol.edge_topology(edge0).num_nodes());

    for ( unsigned edge1 = edge0+1; edge1 < num_edges; edge1++ )
    {
      const unsigned * lnn1 = get_edge_node_ordinals(topol, edge1);
      STK_ThrowAssert(2 == topol.edge_topology(edge1).num_nodes() ||
          3 == topol.edge_topology(edge1).num_nodes());

      int node0 = -1;
      int node1 = -1;
      if (lnn0[0] == lnn1[0])
      {
        node0 = 0; node1 = 0;
      }
      else if (lnn0[0] == lnn1[1])
      {
        node0 = 0; node1 = 1;
      }
      else if (lnn0[1] == lnn1[0])
      {
        node0 = 1; node1 = 0;
      }
      else if (lnn0[1] == lnn1[1])
      {
        node0 = 1; node1 = 1;
      }
      else
      {
        continue;
      }
      const stk::math::Vector3d vec0 = my_nodes[lnn0[1-node0]]->owner_coords(my_owner) - my_nodes[lnn0[node0]]->owner_coords(my_owner);
      const stk::math::Vector3d vec1 = my_nodes[lnn1[1-node1]]->owner_coords(my_owner) - my_nodes[lnn1[node1]]->owner_coords(my_owner);
      const double angle = std::acos( Dot(vec0.unit_vector(),vec1.unit_vector()) );

      //if (angle > 2.4)
      //{
      //  krinolog << "DEBUG: bad angle = " << angle << " between edges=" << edge0 << "," << edge1 << " at node=" << lnn0[node0] << stk::diag::dendl;
      //}

      if (angle > max_angle)
      {
        max_angle = angle;
      }
    }
  }

  return max_angle;
}

void
SubElement::decompose_edges(CDMesh & /*mesh*/, const InterfaceID /*interface_key*/)
{
  const std::string & owner_type = my_owner->topology().name();
  const std::string & sub_type = topology().name();
  ThrowRuntimeError("Subelement decomposition for subelement of type '" << sub_type
      << "' which was generated from owning element of type '" << owner_type
      << "' is missing the capability to generate conformal facets.");
}

void
SubElement::find_refined_edges(std::vector<unsigned> & refined_edges) const
{

  const stk::topology topol = topology();
  const unsigned num_edges = topol.num_edges();
  refined_edges.reserve(num_edges);
  for ( unsigned edge = 0; edge < num_edges; edge++ )
  {
    const unsigned * edge_node_ordinals = get_edge_node_ordinals(topol, edge);

    const int num_edge_nodes = topol.edge_topology(edge).num_nodes();
    STK_ThrowRequire(2 == num_edge_nodes || 3 == num_edge_nodes);

    if ((2 == num_edge_nodes &&
         NULL != SubElementNode::common_child({my_nodes[edge_node_ordinals[0]], my_nodes[edge_node_ordinals[1]]})) ||
        (3 == num_edge_nodes &&
         (NULL != SubElementNode::common_child({my_nodes[edge_node_ordinals[0]], my_nodes[edge_node_ordinals[2]]}) ||
          NULL != SubElementNode::common_child({my_nodes[edge_node_ordinals[1]], my_nodes[edge_node_ordinals[2]]}))))
    {
      refined_edges.push_back(edge);
    }
  }
}

int
SubElement::find_longest_bad_edge(std::vector<unsigned> & bad_edges) const
{

  const stk::topology topol = topology();

  const unsigned num_bad_edges = bad_edges.size();
  std::vector<stk::math::Vector3d> edge_midpt(num_bad_edges,stk::math::Vector3d::ZERO);

  if (0 == num_bad_edges) return -1;

  double max_length = 0;
  unsigned longest_bad_edge_index = 0;
  for ( unsigned index = 0; index < num_bad_edges; index++ )
  {
    const unsigned edge = bad_edges[index];

    const unsigned * const lnn = get_edge_node_ordinals(topol, edge);
    const int num_edge_nodes = topol.edge_topology(edge).num_nodes();

    const SubElementNode * const node0 = my_nodes[lnn[0]];
    const SubElementNode * const node1 = my_nodes[lnn[1]];

    const stk::math::Vector3d & coord0 = node0->coordinates();
    const stk::math::Vector3d & coord1 = node1->coordinates();

    const double edge_straight_length = (coord0 - coord1).length();
    STK_ThrowRequire(edge_straight_length > 0.0);

    if (2 == num_edge_nodes)
    {
      edge_midpt[index] = 0.5*(coord0 + coord1);
    }
    else
    {
      STK_ThrowAssert(3 == num_edge_nodes);
      edge_midpt[index] = my_nodes[lnn[0]]->coordinates();
    }

    // we need an absolute mechanism for selecting the edge to bisect so that all elements that share
    // common edges will make the same decisions
    if (utility::is_more(edge_straight_length,max_length))
    {
      longest_bad_edge_index = index;
      max_length = edge_straight_length;
    }
    else if (!utility::is_less(edge_straight_length,max_length)) // tie breaker
    {
      const stk::math::Vector3d & edge_midside_coords = edge_midpt[index];
      // note that it is safe to assume that longest_bad_edge is already assigned if edge_length == max_length
      const stk::math::Vector3d longest_edge_midside_coords = edge_midpt[longest_bad_edge_index];

      STK_ThrowAssert((utility::is_not_equal(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
                    utility::is_not_equal(edge_midside_coords[1],longest_edge_midside_coords[1])));

      if (utility::is_more(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
          (!utility::is_less(edge_midside_coords[0],longest_edge_midside_coords[0]) &&
            (utility::is_more(edge_midside_coords[1],longest_edge_midside_coords[1]))))
      {
        longest_bad_edge_index = index;
        max_length = edge_straight_length;
      }
    }
  }
  return bad_edges[longest_bad_edge_index];
}

int
SubElement::parent_side_id(const int iside) const
{
  return my_parent_side_ids[iside];
}

void
SubElement::debug_subelements(const NodeVec & lnodes, const InterfaceID & interface, const int case_id) const
{
  krinolog << "owner_id=" << my_owner->entityId() << ", after cutting with interface " << interface << ", case_id=" << case_id << stk::diag::dendl;

  for (unsigned n=0; n<lnodes.size(); ++n)
  {
    const SubElementNode * node = lnodes[n];
    if (NULL != node)
    {
      krinolog << "  Node[" << n << "]: coords=" << node->coordinates() << ": ["<< node->get_ancestry() << "]\n";
    }
  }
  krinolog << " interface signs = ";
  const std::vector<InterfaceID> interfaces = my_owner->get_sorted_cutting_interfaces();
  STK_ThrowRequire(interfaces.size() == myInterfaceSigns.size());
  for (unsigned i=0; i<interfaces.size(); ++i) krinolog << interfaces[i] << "@" << myInterfaceSigns[i] << " ";
  krinolog << "\n";
  for (unsigned subid=0; subid<my_subelements.size(); ++subid)
  {
    krinolog << "Subelement[" << subid << "]:\n";
    my_subelements[subid]->debug();
  }
  krinolog << stk::diag::dendl;
}

void
SubElement::debug() const
{
  const double sub_vol = relative_volume();
  krinolog << "  owner_id=" << my_owner->entityId() << ", relative_volume=" << sub_vol << ", interface signs = ";
  const std::vector<InterfaceID> interfaces = my_owner->get_sorted_cutting_interfaces();
  STK_ThrowRequire(interfaces.size() == myInterfaceSigns.size());
  for (unsigned i=0; i<interfaces.size(); ++i) krinolog << interfaces[i] << "@" << myInterfaceSigns[i] << " ";
  krinolog << "\n";
  if (true || sub_vol < 1.e-10)
  {
    for (unsigned n=0; n<my_nodes.size(); ++n)
    {
      const SubElementNode * node = my_nodes[n];
      if (NULL != node)
      {
        krinolog << "  SubNode[" << n << "]: coords=" << node->coordinates() << ": ["<< node->get_ancestry() << "] with domains { ";
        for (int domain : node->get_sorted_node_domains()) krinolog << domain << " ";
        krinolog << "}\n";
      }
    }
    for (unsigned n=0; n<num_sides(); ++n)
    {
      krinolog << "  SubSide[" << n << "]: on parents side " << parent_side_id(n) << "\n";
    }
  }
}

bool
SubElement::have_interface(const InterfaceID & interface) const
{
  return my_owner->have_interface(interface) && (myInterfaceSigns[my_owner->get_interface_index(interface)] == 0);
}

SubElement_Tri_6::SubElement_Tri_6(
  const NodeVec & nodes,
  const std::vector<int> & parent_side_ids,
  const Mesh_Element * owner)
    : SubElement( stk::topology::TRIANGLE_6_2D,
                   nodes,
                   parent_side_ids,
                   owner)
{
}

SubElement_Tet_10::SubElement_Tet_10(
  const NodeVec & nodes,
  const std::vector<int> & parent_side_ids,
  const Mesh_Element * owner)
    : SubElement( stk::topology::TETRAHEDRON_10,
                   nodes,
                   parent_side_ids,
                   owner)
{
}

SubElement_Tri_3::SubElement_Tri_3(
  const NodeVec & nodes,
  const std::vector<int> & parent_side_ids,
  const Mesh_Element * owner)
    : SubElement( stk::topology::TRIANGLE_3_2D,
                   nodes,
                   parent_side_ids,
                   owner)
{
}

void
SubElement_Tri_3::build_quadratic_subelements(CDMesh & mesh)
{

  if ( my_subelements.size() > 0 )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->build_quadratic_subelements(mesh);
    }
    return;
  }

  // create 1, 6-noded tri
  NodeVec sub_nodes = my_nodes;
  sub_nodes.resize(6,(SubElementNode *)NULL);

  sub_nodes[3] = mesh.create_midside_node(my_owner, my_nodes[0], my_nodes[1]);
  sub_nodes[4] = mesh.create_midside_node(my_owner, my_nodes[1], my_nodes[2]);
  sub_nodes[5] = mesh.create_midside_node(my_owner, my_nodes[2], my_nodes[0]);

  std::unique_ptr<SubElement> sub = std::make_unique<SubElement_Tri_6>( sub_nodes, my_parent_side_ids, my_owner );
  sub->set_interface_signs(get_interface_signs());
  add_subelement( std::move(sub) );
}

void SubElement_Tri_3::cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains)
{
  const std::vector<double> weights{1.-pCoords[0]-pCoords[1], pCoords[0], pCoords[1]};

  bool badCut = false;
  for (auto && weight : weights)
  {
    if(weight < mesh.get_snapper().get_edge_tolerance())
    {
      if (krinolog.shouldPrint(LOG_DEBUG))
        krinolog << "Skipping cut of interior intersection point because of quality." << stk::diag::dendl;
      badCut = true;
      break;
    }
  }

  if (!badCut)
  {
    const SubElementNode * cutNode = mesh.create_child_internal_or_face_node( my_owner, my_nodes, weights );
    cutNode->set_node_domains(sortedDomains);

    NodeVec lnodes = my_nodes;
    lnodes.push_back(cutNode);

    handle_tri(lnodes, get_interface_signs(), 0,1,3, 0,-1,-1, false,false,false);
    handle_tri(lnodes, get_interface_signs(), 1,2,3, 1,-1,-1, false,false,false);
    handle_tri(lnodes, get_interface_signs(), 2,0,3, 2,-1,-1, false,false,false);
  }
}

void
SubElement_Tri_3::determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key)
{
  bool setSignOnEdge = false;
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 0,1 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 1,2 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 2,0 );

  if (!setSignOnEdge && myInterfaceSigns[my_owner->get_interface_index(interface_key)]==0)
    set_node_signs_on_uncrossed_subelement(interface_key);

}

void
SubElement_Tri_3::determine_node_scores(const CDMesh & /*mesh*/, const InterfaceID /*interface_key*/)
{
  // No-op, node scores are not used on tris
}

void
SubElement_Tri_3::decompose_edges(CDMesh & mesh, const InterfaceID interface_key)
{
  process_edge( mesh, interface_key, 0,1 );
  process_edge( mesh, interface_key, 1,2 );
  process_edge( mesh, interface_key, 2,0 );
}

void
SubElement_Tri_3::fix_hanging_children(CDMesh & mesh, const InterfaceID & interface, const std::vector<int> & edges_with_children)
{
  int edge_case_id = 0;
  for (auto edge_with_children : edges_with_children)
  {
    edge_case_id += 1<<edge_with_children;
  }
  // It is invalid to have all three edges with hanging children.
  STK_ThrowErrorMsgIf(edge_case_id == 7, "Found Tri 3, with invalid configuration of edges with hanging children.");

  std::array<int,3> node_signs = {{0, 0, 0}};

  if (edge_case_id == 0) // uncut subelement
  {
    if (my_nodes[0]->node_sign_is_set()) node_signs[0] = my_nodes[0]->get_node_sign();
    if (my_nodes[1]->node_sign_is_set()) node_signs[1] = my_nodes[1]->get_node_sign();
    if (my_nodes[2]->node_sign_is_set()) node_signs[2] = my_nodes[2]->get_node_sign();

    if (!(node_signs[0] >= 0 && node_signs[1] >= 0 && node_signs[2] >= 0) &&
        !(node_signs[0] <= 0 && node_signs[1] <= 0 && node_signs[2] <= 0))
    {
      return;
    }
  }
  else
  {
    for (auto && edge_with_child : edges_with_children)
    {
      const unsigned * edge_node_ordinals = get_edge_node_ordinals(topology(), edge_with_child);
      const int i0 = edge_node_ordinals[0];
      const int i1 = edge_node_ordinals[1];
      node_signs[i0] = my_nodes[i0]->get_node_sign();
      node_signs[i1] = my_nodes[i1]->get_node_sign();
    }
  }

  perform_decomposition(mesh, interface, node_signs);
}

void
SubElement_Tri_3::perform_decomposition(CDMesh & mesh, const InterfaceID interface_key, const std::array<int,3> & node_signs)
{
  const int caseId = ContourTri::compute_case_id(node_signs);

  if (caseId == 0  || // ls[0]<0 && ls[1]<0 && ls[2]<0
      caseId == 13 || // ls[0]=0 && ls[1]=0 && ls[2]=0
      caseId == 26)   // ls[0]>0 && ls[1]>0 && ls[2]>0
  {
    const int sign = (caseId==0) ? -1 : 1;
    update_interface_signs(interface_key, sign);
    return;
  }

  NodeVec lnodes = my_nodes;
  lnodes.resize(6,(SubElementNode *)NULL);

  const std::array<unsigned,6> & permute = ContourTri::get_permuted_node_ordinals(caseId);
  const int permutedCaseId = ContourTri::get_permuted_case_id(caseId);

  const unsigned i0 = permute[0];
  const unsigned i1 = permute[1];
  const unsigned i2 = permute[2];
  const unsigned i3 = permute[3];
  const unsigned i5 = permute[5];

  // nodes and sides permute the same way
  const unsigned s0 = permute[0];
  const unsigned s1 = permute[1];
  const unsigned s2 = permute[2];

  const Simplex_Generation_Method simplexMethod = mesh.get_cdfem_support().get_simplex_generation_method();

  switch (permutedCaseId)
  {
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0
    {
      update_interface_signs(interface_key, -1);
    }
    break;

    case 25: // ls[0]=0 && ls[1]>0 && ls[2]>0
    {
      update_interface_signs(interface_key, +1);
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0
    case 24: // ls[0]<0 && ls[1]>0 && ls[2]>0
    {
      lnodes[i3] = SubElementNode::common_child({lnodes[i0], lnodes[i1]});
      lnodes[i5] = SubElementNode::common_child({lnodes[i2], lnodes[i0]});
      STK_ThrowAssert(nullptr != lnodes[i3] && nullptr != lnodes[i5]);

      const bool diag = determine_diagonal_for_cut_triangle(simplexMethod, lnodes, i0, i1, i2, i3, i5);

      const int sign = (permutedCaseId==2) ? 1 : -1;
      handle_tri(lnodes, subelement_interface_signs(interface_key, sign), i0,i3,i5, s0,-1,s2, false,true,false);
      handle_quad(mesh, lnodes, subelement_interface_signs(interface_key, -sign), i3,i1,i2,i5, s0,s1,s2,-1, false,false,false,true, diag);
    }
    break;

    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0
    case 22: // ls[0]=0 && ls[1]=0 && ls[2]>0
    {
      const int sign = (permutedCaseId==4) ? -1 : 1;
      update_interface_signs(interface_key, sign);
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0
    case 21: // ls[0]<0 && ls[1]=0 && ls[2]>0
    {
      lnodes[i5] = SubElementNode::common_child({lnodes[i2], lnodes[i0]});
      STK_ThrowAssert(nullptr != lnodes[i5]);

      const int sign = (permutedCaseId==5) ? 1 : -1;
      handle_tri(lnodes, subelement_interface_signs(interface_key,  sign), i1,i5,i0, -1,s2,s0, true,false,false);
      handle_tri(lnodes, subelement_interface_signs(interface_key, -sign), i5,i1,i2, -1,s1,s2, true,false,false);
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error. caseId, permutedCaseId=" << caseId << "," << permutedCaseId);
  }

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    debug_subelements(lnodes, interface_key, caseId);
  }
}

bool
SubElement_Tri_3::determine_diagonal_for_cut_triangle(const Simplex_Generation_Method & simplexMethod, const NodeVec & lnodes, const int /*i0*/, const int i1, const int i2, const int i3, const int i5)
{
  /*
   *     2 o
   *      / \
   *   5 o   \
   *    / \   \
   *   o---o---o
   *   0   3   1
   */

  // true:  connect nodes 3 and 2
  // false: connect nodes 5 and 1

  if (simplexMethod == CUT_QUADS_BY_GLOBAL_IDENTIFIER)
  {
    SubElementNodeAncestry ancestry1 = lnodes[i1]->get_ancestry();
    SubElementNodeAncestry ancestry2 = lnodes[i2]->get_ancestry();

    return ancestry2 < ancestry1;
  }
  else
  {
    STK_ThrowAssert(simplexMethod == CUT_QUADS_BY_DEFAULT_METHOD || simplexMethod == CUT_QUADS_BY_LARGEST_ANGLE);

    // Angle-based scheme
    // Select diagonal that cuts largest angle in quad.  Since there isn't an issue with
    // conforming, this is always possible (unlike tet).
    return ElementObj::will_cutting_quad_from_0to2_cut_largest_angle(lnodes[i3],lnodes[i1],lnodes[i2],lnodes[i5]);
  }
}

void
SubElement_Tri_3::handle_tri( NodeVec & lnodes,
                               const std::vector<int> & subInterfaceSigns,
                               const int i0, const int i1, const int i2,
                               const int s0, const int s1, const int s2,
                               const bool /*is_interface0*/, const bool /*is_interface1*/, const bool /*is_interface2*/)
{
  STK_ThrowAssert(!is_degenerate(lnodes, i0,i1,i2));

  NodeVec sub_nodes(3,(SubElementNode *)NULL);
  std::vector<int> sub_parent_ids(3);

  sub_nodes[0] = lnodes[i0];
  sub_nodes[1] = lnodes[i1];
  sub_nodes[2] = lnodes[i2];
  sub_parent_ids[0] = (s0>= 0) ? my_parent_side_ids[s0] : -1;
  sub_parent_ids[1] = (s1>= 0) ? my_parent_side_ids[s1] : -1;
  sub_parent_ids[2] = (s2>= 0) ? my_parent_side_ids[s2] : -1;

  std::unique_ptr<SubElement> sub = std::make_unique<SubElement_Tri_3>( sub_nodes, sub_parent_ids, my_owner);
  sub->set_interface_signs(subInterfaceSigns);
  add_subelement( std::move(sub) );
}

void
SubElement_Tri_3::handle_quad( CDMesh & /*mesh*/,
    NodeVec & lnodes,
    const std::vector<int> & subInterfaceSigns,
    const int i0, const int i1, const int i2, const int i3,
    const int s0, const int s1, const int s2, const int s3,
    const bool is_interface0, const bool is_interface1, const bool is_interface2, const bool is_interface3,
    const bool face )
{
#if 1
  if (face)
  {
    handle_tri( lnodes, subInterfaceSigns, i0,i1,i2, s0,s1,-1, is_interface0,is_interface1,false );
    handle_tri( lnodes, subInterfaceSigns, i2,i3,i0, s2,s3,-1, is_interface2,is_interface3,false );
  }
  else
  {
    handle_tri( lnodes, subInterfaceSigns, i0,i1,i3, s0,-1,s3, is_interface0,false,is_interface3 );
    handle_tri( lnodes, subInterfaceSigns, i2,i3,i1, s2,-1,s1, is_interface2,false,is_interface1 );
  }
#else
  NodeVec quad_nodes;
  std::vector<double> weights;
  const double x0 = 0;
  const double y0 = 0;

  quad_nodes.push_back(lnodes[i0]);
  weights.push_back(0.25*(1-x0)*(1-y0));
  quad_nodes.push_back(lnodes[i1]);
  weights.push_back(0.25*(1+x0)*(1-y0));
  quad_nodes.push_back(lnodes[i2]);
  weights.push_back(0.25*(1+x0)*(1+y0));
  quad_nodes.push_back(lnodes[i3]);
  weights.push_back(0.25*(1-x0)*(1+y0));

  lnodes[6] = mesh.create_internal_node( my_owner, quad_nodes, weights );

  handle_tri(lnodes, ls_index, sign, i0,i1,6, s0,-1,-1);
  handle_tri(lnodes, ls_index, sign, i1,i2,6, s1,-1,-1);
  handle_tri(lnodes, ls_index, sign, i2,i3,6, s2,-1,-1);
  handle_tri(lnodes, ls_index, sign, i3,i0,6, s3,-1,-1);
#endif
}

bool
SubElement_Tri_3::is_degenerate( NodeVec & lnodes,
                                  const int i0, const int i1, const int i2 )
{

  if ( lnodes[i0] == lnodes[i1] ||
       lnodes[i0] == lnodes[i2] ||
       lnodes[i1] == lnodes[i2] )
  {
    // this tri is degenerate with two coincident nodes
    return true;
  }

  return false;
}

SubElement_Tet_4::SubElement_Tet_4(
  const NodeVec & nodes,
  const std::vector<int> & parent_side_ids,
  const Mesh_Element * owner)
    : SubElement( stk::topology::TETRAHEDRON_4,
                   nodes,
                   parent_side_ids,
                   owner)
{
}

void SubElement_Tet_4::cut_face_intersection_point_with_permutation(CDMesh & mesh, const std::array<int,4> & permuteNodes, const std::array<int,4> & permuteSides, const std::vector<double> & faceNodeWeights, const std::vector<int> & sortedDomains)
{
  const SubElementNode * cutNode = mesh.create_child_internal_or_face_node(my_owner,
      {my_nodes[permuteNodes[0]], my_nodes[permuteNodes[1]], my_nodes[permuteNodes[2]]},
      {faceNodeWeights[0], faceNodeWeights[1], faceNodeWeights[2]});
  const auto & previousDomains = cutNode->get_sorted_node_domains();
  STK_ThrowRequire(previousDomains.empty() || sortedDomains == previousDomains);
  cutNode->set_node_domains(sortedDomains);

  NodeVec lnodes;
  lnodes.reserve(5);
  for (int i=0; i<4; ++i)
    lnodes.push_back(my_nodes[permuteNodes[i]]);
  lnodes.push_back(cutNode);

  handle_tet(lnodes, get_interface_signs(), 0,1,4,3, permuteSides[0],-1,-1,permuteSides[3]);
  handle_tet(lnodes, get_interface_signs(), 1,2,4,3, permuteSides[1],-1,-1,permuteSides[3]);
  handle_tet(lnodes, get_interface_signs(), 2,0,4,3, permuteSides[2],-1,-1,permuteSides[3]);
}

void SubElement_Tet_4::cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains)
{
  const std::vector<double> weights{1.-pCoords[0]-pCoords[1]-pCoords[2], pCoords[0], pCoords[1], pCoords[2]};

  bool badCut = false;
  for (auto && weight : weights)
  {
    if(weight < mesh.get_snapper().get_edge_tolerance())
    {
      if (krinolog.shouldPrint(LOG_DEBUG))
        krinolog << "Skipping cut of interior intersection point because of quality." << stk::diag::dendl;
      badCut = true;
      break;
    }
  }

  if (!badCut)
  {
    const SubElementNode * cutNode = mesh.create_child_internal_or_face_node( my_owner, my_nodes, weights );
    cutNode->set_node_domains(sortedDomains);

    NodeVec lnodes = my_nodes;
    lnodes.push_back(cutNode);

    handle_tet(lnodes, get_interface_signs(), 0,3,1,4, -1,-1,-1,0);
    handle_tet(lnodes, get_interface_signs(), 1,3,2,4, -1,-1,-1,1);
    handle_tet(lnodes, get_interface_signs(), 0,2,3,4, -1,-1,-1,2);
    handle_tet(lnodes, get_interface_signs(), 0,1,2,4, -1,-1,-1,3);
  }
}

void
SubElement_Tet_4::build_quadratic_subelements(CDMesh & mesh)
{

  if ( my_subelements.size() > 0 )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->build_quadratic_subelements(mesh);
    }
    return;
  }

  // create 1, 10-noded tet
  NodeVec sub_nodes = my_nodes;
  sub_nodes.resize(10,(SubElementNode *)NULL);

  const stk::topology tet10_topology = stk::topology::TETRAHEDRON_10;

  for (unsigned edge_i=0; edge_i<tet10_topology.num_edges(); ++edge_i)
  {
    const unsigned * tet10_edge_lnn = get_edge_node_ordinals(tet10_topology, edge_i);
    sub_nodes[tet10_edge_lnn[2]] = mesh.create_midside_node(my_owner, my_nodes[tet10_edge_lnn[0]], my_nodes[tet10_edge_lnn[1]]);
  }

  std::unique_ptr<SubElement> sub = std::make_unique<SubElement_Tet_10>( sub_nodes, my_parent_side_ids, my_owner );
  sub->set_interface_signs(get_interface_signs());
  add_subelement( std::move(sub) );
}

struct FaceIntersection
{
  FaceIntersection(const int iFace, const stk::math::Vector3d & coords, const std::vector<int> & domains)
  : face(iFace),
    parametricCoords(coords),
    sortedDomains(domains) {}

  int face;
  stk::math::Vector3d parametricCoords;
  std::vector<int> sortedDomains;
};

void
SubElement_Tet_4::cut_face_interior_intersection_points(CDMesh & mesh, int level)
{

  if ( my_subelements.size() > 0 )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->cut_face_interior_intersection_points(mesh, level);
    }
    return;
  }

  static constexpr std::array<std::array<int,4>,4> permuteNodes{{ {{0,3,1,2}}, {{1,3,2,0}}, {{0,2,3,1}}, {{0,1,2,3}} }};
  static constexpr std::array<std::array<int,4>,4> permuteSides{{ {{2,1,3,0}}, {{0,2,3,1}}, {{3,1,0,2}}, {{0,1,2,3}} }};

  std::vector<const SubElementNode*> faceNodes(3);
  // Note: this recursively cuts the first identified face interior intersection point
  for (unsigned iFace=0; iFace<topology().num_faces(); ++iFace)
  {
    for (int i=0; i<3; ++i)
      faceNodes[i] = my_nodes[permuteNodes[iFace][i]];

    std::vector<ElementIntersection> faceIntersections;
    my_owner->fill_face_interior_intersections(faceNodes, faceIntersections);

    for (auto && faceIntersection : faceIntersections)
    {
      const stk::math::Vector3d & faceCoords = faceIntersection.parametricCoords;
      const std::vector<double> faceNodeWeights{1.-faceCoords[0]-faceCoords[1], faceCoords[0], faceCoords[1]};
      bool badCut = false;
      for (auto && weight : faceNodeWeights)
      {
        if (weight < mesh.get_snapper().get_edge_tolerance())
        {
          if (krinolog.shouldPrint(LOG_DEBUG))
            krinolog << "Skipping cut of interior face intersection point because of quality." << stk::diag::dendl;
          badCut = true;
          break;
        }
      }
      if (!badCut)
      {
        STK_ThrowRequireMsg(level < 8, "Face cut recursion level exceeded.");
        cut_face_intersection_point_with_permutation(mesh, permuteNodes[iFace], permuteSides[iFace], faceNodeWeights, faceIntersection.sortedDomains);
        cut_face_interior_intersection_points(mesh, ++level);
        return;
      }
    }
  }
}

double SubElement_Tet_4::tet_volume(const std::array<stk::math::Vector3d,4> & nodes)
{
  return Dot(nodes[3]-nodes[0],Cross(nodes[1]-nodes[0], nodes[2]-nodes[0]))/6.0;
}

void
SubElement_Tet_4::fix_hanging_children(CDMesh & mesh, const InterfaceID & interface_key, const std::vector<int> & edges_with_children)
{
  int edge_case_id = 0;
  for (unsigned i=0; i<edges_with_children.size(); ++i)
  {
    edge_case_id += 1<<edges_with_children[i];
  }

  // This is figured out in a table that covers all possible combinations of hanging children.
  // Invalid cases occur when the combination of hanging children doesn't make sense.  For example,
  // it is invalid to have all three edges of a face with hanging children.
  static const int case_id_from_edge_case_id [] =
    { 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,  // 0-9
      2, 3, 1, 1, 4, 0, 1, 1, 1, 1,  // 10-19
      2, 4, 3, 0, 1, 0, 4, 0, 3, 0,  // 20-29
      1, 0, 1, 2, 1, 4, 1, 3, 1, 0,  // 30-39
      1, 4, 3, 1, 0, 0, 0, 0, 1, 3,  // 40-49
      0, 0, 4, 1, 0, 0, 1, 0, 0, 0,  // 50-59
      0, 0, 0, 0                     // 60-63
    };
  const int case_id = case_id_from_edge_case_id[edge_case_id];
  STK_ThrowErrorMsgIf(case_id == 0, "Found Tet 4, with invalid configuration of edges with hanging children.");

  if (case_id == 1)
  {
    std::array<int,4> node_signs = {{0, 0, 0, 0}};

    if (edge_case_id == 0) // uncut subelement
    {
      if (my_nodes[0]->node_sign_is_set()) node_signs[0] = my_nodes[0]->get_node_sign();
      if (my_nodes[1]->node_sign_is_set()) node_signs[1] = my_nodes[1]->get_node_sign();
      if (my_nodes[2]->node_sign_is_set()) node_signs[2] = my_nodes[2]->get_node_sign();
      if (my_nodes[3]->node_sign_is_set()) node_signs[3] = my_nodes[3]->get_node_sign();

      if (!(node_signs[0] >= 0 && node_signs[1] >= 0 && node_signs[2] >= 0 && node_signs[3] >= 0) &&
          !(node_signs[0] <= 0 && node_signs[1] <= 0 && node_signs[2] <= 0 && node_signs[3] <= 0))
      {
        return;
      }
    }
    else
    {
      for (auto && edge_with_child : edges_with_children)
      {
        const unsigned * edge_node_ordinals = get_edge_node_ordinals(topology(), edge_with_child);
        const int i0 = edge_node_ordinals[0];
        const int i1 = edge_node_ordinals[1];
        node_signs[i0] = my_nodes[i0]->get_node_sign();
        node_signs[i1] = my_nodes[i1]->get_node_sign();
      }
    }

    perform_decomposition(mesh, interface_key, node_signs);
  }
  else
  {
    STK_ThrowAssert(case_id == 2 || case_id == 3 || case_id == 4);

    static const int edge_case_permutations [] =
      { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // 0-9
         2,  6, -1, -1,  5, -1, -1, -1, -1, -1,  // 10-19
         0,  8,  5, -1, -1, -1,  0, -1,  0, -1,  // 20-29
        -1, -1, -1,  1, -1,  4, -1,  3, -1, -1,  // 30-39
        -1,  2,  2, -1, -1, -1, -1, -1, -1,  1,  // 40-49
        -1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  // 50-59
        -1, -1, -1, -1                           // 60-63
      };

    const int permutation = edge_case_permutations[edge_case_id];
    STK_ThrowRequire(permutation >= 0);

    const std::array<unsigned, 10> & permutedNodeOrdinals = ContourTet::get_permuted_node_ordinals_for_permutation(permutation);
    const unsigned i0 = permutedNodeOrdinals[0];
    const unsigned i1 = permutedNodeOrdinals[1];
    const unsigned i2 = permutedNodeOrdinals[2];
    const unsigned i3 = permutedNodeOrdinals[3];
    const unsigned i5 = permutedNodeOrdinals[5];
    const unsigned i6 = permutedNodeOrdinals[6];
    const unsigned i7 = permutedNodeOrdinals[7];
    const unsigned i8 = permutedNodeOrdinals[8];

    const std::array<unsigned, 4> & permutedSideOrdinals = ContourTet::get_permuted_side_ordinals_for_permutation(permutation);
    const unsigned s0 = permutedSideOrdinals[0];
    const unsigned s1 = permutedSideOrdinals[1];
    const unsigned s2 = permutedSideOrdinals[2];
    const unsigned s3 = permutedSideOrdinals[3];

    NodeVec lnodes = my_nodes;
    lnodes.resize(10,(SubElementNode *)NULL);

    const int arbitrary_sign = 0;
    const auto subInterfaceSigns = subelement_interface_signs(interface_key, arbitrary_sign);
    const Simplex_Generation_Method simplexMethod = mesh.get_cdfem_support().get_simplex_generation_method();
    const bool globalIDsAreParallelConsistent = mesh.get_cdfem_support().get_global_ids_are_parallel_consistent();

    if (case_id == 2)
    {
      lnodes[i6] = SubElementNode::common_child({lnodes[i0], lnodes[i2]});
      lnodes[i8] = SubElementNode::common_child({lnodes[i1], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i6] && nullptr != lnodes[i8]);


      handle_tet( lnodes, subInterfaceSigns, i6,i2,i3,i8, -1,s1,-1,s2 );
      handle_tet( lnodes, subInterfaceSigns, i3,i0,i6,i8, s0,-1,-1,s2 );
      handle_tet( lnodes, subInterfaceSigns, i6,i0,i1,i8, -1,s0,-1,s3 );
      handle_tet( lnodes, subInterfaceSigns, i1,i2,i6,i8, s1,-1,-1,s3 );
    }
    else if (case_id == 3)
    {
      lnodes[i6] = SubElementNode::common_child({lnodes[i0], lnodes[i2]});
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      lnodes[i8] = SubElementNode::common_child({lnodes[i1], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i6] && nullptr != lnodes[i7] && nullptr != lnodes[i8]);

      const bool face0 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i3, i0, i1, i7, i8);
      const bool face2 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i3, i2, i7, i6);

      // Connect 6-8
      handle_tet( lnodes, subInterfaceSigns, i1,i2,i6,i8, s1,-1,-1,s3 );
      handle_pyramid( lnodes, subInterfaceSigns, i1,i0,i7,i8,i6, s3,s2,-1,-1,s0, face0 );
      handle_pyramid( lnodes, subInterfaceSigns, i2,i3,i7,i6,i8, s1,s0,-1,-1,s2, face2 );
    }
    else if (case_id == 4)
    {
      lnodes[i5] = SubElementNode::common_child({lnodes[i1], lnodes[i2]});
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      lnodes[i8] = SubElementNode::common_child({lnodes[i1], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i5] && nullptr != lnodes[i7] && nullptr != lnodes[i8]);

      const bool face0 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i3, i0, i1, i7, i8);
      const bool face1 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i1, i2, i3, i5, i8);

      // Connect 5-7
      handle_tet( lnodes, subInterfaceSigns, i0,i5,i2,i7, -1,-1,s2,s3 );
      handle_pyramid( lnodes, subInterfaceSigns, i1,i0,i7,i8,i5, s3,-1,-1,s1,s0, face0 );
      handle_pyramid( lnodes, subInterfaceSigns, i5,i8,i3,i2,i7, -1,s0,s2,-1,s1, face1 );
    }

    if (krinolog.shouldPrint(LOG_DEBUG))
    {
      debug_subelements(lnodes, interface_key, case_id);
    }
  }
}

void
SubElement_Tet_4::determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key)
{
  STK_ThrowAssert(!have_subelements());

  bool setSignOnEdge = false;
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 0,1 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 1,2 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 0,2 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 0,3 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 1,3 );
  setSignOnEdge |= determine_node_signs_on_edge( mesh, interface_key, 2,3 );

  if (!setSignOnEdge && myInterfaceSigns[my_owner->get_interface_index(interface_key)]==0)
    set_node_signs_on_uncrossed_subelement(interface_key);
}

void
SubElement_Tet_4::determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key)
{
  STK_ThrowAssert(!have_subelements());

  determine_node_scores_on_face( mesh, interface_key, 0,1,3 );
  determine_node_scores_on_face( mesh, interface_key, 1,2,3 );
  determine_node_scores_on_face( mesh, interface_key, 2,0,3 );
  determine_node_scores_on_face( mesh, interface_key, 0,2,1 );
}


static std::pair<double,double> get_quad_angle_measures(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & x3)
{
  std::array<stk::math::Vector3d,4> sides{x1-x0, x2-x1, x3-x2, x0-x3};
  for (auto && side : sides) side.unitize();

  // return measure02,measure13 where measureAB=std::max(-cos(A),-cos(B))
  return {std::max(Dot(sides[3],sides[0]), Dot(sides[1],sides[2])), std::max(Dot(sides[0],sides[1]), Dot(sides[2],sides[3]))};
}

static void determine_node_scores_on_triangle_face( const Simplex_Generation_Method simplexMethod, const std::array<const SubElementNode *,5> & faceNodes )
{
  /*
   *     2 o
   *      / \
   *     /   o 4
   *    /   / \
   *   o---o---o
   *   0   3   1
   */

  if (simplexMethod == CUT_QUADS_BY_DEFAULT_METHOD || simplexMethod == CUT_QUADS_BY_NEAREST_EDGE_CUT)
  {
    //  nodal edge cut length based criterion
    //  Use globally consistent comparison at nodes based on the shortest relative edge length for the cut edges that use the node.

    //  The general idea is to prefer edges that emanate away from nodes that have nearby cuts.  This, for perfectly shaped elements,
    //  will cut the largest angles.

    for (int iChild=3; iChild<5; ++iChild)
    {
      const SubElementNode * child = faceNodes[iChild];
      const NodeVec & parents = child->get_parents();
      const std::vector<double> parent_weights = child->get_parent_weights();
      for (size_t i=0; i<parents.size(); ++i)
        if (parents[i] == faceNodes[0] || parents[i] == faceNodes[2])
          parents[i]->set_node_score(1.-parent_weights[i]);
    }
  }
  else if (simplexMethod == CUT_QUADS_BY_LARGEST_ANGLE)
  {
    // nodal face angle criterion
    // Use globally consistent comparison at nodes based on the largest angle for the cut faces that use the node.
    // This does not rely on perfectly shaped elements to cut the largest angles.

    const std::pair<double,double> measure23Andmeasure04 = get_quad_angle_measures(faceNodes[2]->coordinates(), faceNodes[0]->coordinates(), faceNodes[3]->coordinates(), faceNodes[4]->coordinates());
    faceNodes[2]->set_node_score(measure23Andmeasure04.first);
    faceNodes[0]->set_node_score(measure23Andmeasure04.second);
  }
}

void
SubElement_Tet_4::determine_node_scores_on_face( const CDMesh & mesh, const InterfaceID /*interface*/, const int i0, const int i1, const int i2 )
{
  const Simplex_Generation_Method simplexMethod = mesh.get_cdfem_support().get_simplex_generation_method();

  const SubElementNode * parent0 = my_nodes[i0];
  const SubElementNode * parent1 = my_nodes[i1];
  const SubElementNode * parent2 = my_nodes[i2];

  if (parent0->get_node_on_interface() || parent1->get_node_on_interface() || parent2->get_node_on_interface()) return;

  const SubElementNode * child0 = SubElementNode::common_child({parent0, parent1});
  const SubElementNode * child1 = SubElementNode::common_child({parent1, parent2});
  const SubElementNode * child2 = SubElementNode::common_child({parent2, parent0});

  const int caseId =
      ((child0 == nullptr) ? 0 : 1) +
      ((child1 == nullptr) ? 0 : 2) +
      ((child2 == nullptr) ? 0 : 4);

  if (caseId == 3)
    determine_node_scores_on_triangle_face(simplexMethod, {parent0, parent1, parent2, child0, child1});
  else if (caseId == 5)
    determine_node_scores_on_triangle_face(simplexMethod, {parent2, parent0, parent1, child2, child0});
  else if (caseId == 6)
    determine_node_scores_on_triangle_face(simplexMethod, {parent1, parent2, parent0, child1, child2});
}

void
SubElement_Tet_4::decompose_edges(CDMesh & mesh, const InterfaceID interface_key)
{

  process_edge( mesh, interface_key, 0,1 );
  process_edge( mesh, interface_key, 1,2 );
  process_edge( mesh, interface_key, 0,2 );
  process_edge( mesh, interface_key, 0,3 );
  process_edge( mesh, interface_key, 1,3 );
  process_edge( mesh, interface_key, 2,3 );
}

void
SubElement_Tet_4::perform_decomposition(CDMesh & mesh, const InterfaceID interface_key, const std::array<int,4> & node_signs)
{
  // create between 4 to 6 conforming tetrahedral subelements

  const int caseId = ContourTet::compute_case_id(node_signs);

  if (caseId == 0  || // ls[0]<0 && ls[1]<0 && ls[2]<0 && ls[3]<0
      caseId == 40 || // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]=0
      caseId == 80)   // ls[0]>0 && ls[1]>0 && ls[2]>0 && ls[3]>0
  {
    const int sign = (caseId==0) ? -1 : 1;
    update_interface_signs(interface_key, sign);
    return;
  }

  NodeVec lnodes = my_nodes;
  lnodes.resize(10,(SubElementNode *)NULL);

  const int permutedCaseId = ContourTet::get_permuted_case_id(caseId);
  const std::array<unsigned,10> & permutedNodeOrdinals = ContourTet::get_permuted_node_ordinals(caseId);

  const unsigned i0 = permutedNodeOrdinals[0];
  const unsigned i1 = permutedNodeOrdinals[1];
  const unsigned i2 = permutedNodeOrdinals[2];
  const unsigned i3 = permutedNodeOrdinals[3];
  const unsigned i4 = permutedNodeOrdinals[4];
  const unsigned i5 = permutedNodeOrdinals[5];
  const unsigned i6 = permutedNodeOrdinals[6];
  const unsigned i7 = permutedNodeOrdinals[7];
  const unsigned i8 = permutedNodeOrdinals[8];

  const std::array<unsigned, 4> & permutedSideOrdinals = ContourTet::get_permuted_side_ordinals(caseId);
  const unsigned s0 = permutedSideOrdinals[0];
  const unsigned s1 = permutedSideOrdinals[1];
  const unsigned s2 = permutedSideOrdinals[2];
  const unsigned s3 = permutedSideOrdinals[3];

  const Simplex_Generation_Method simplexMethod = mesh.get_cdfem_support().get_simplex_generation_method();
  const bool globalIDsAreParallelConsistent = mesh.get_cdfem_support().get_global_ids_are_parallel_consistent();

  switch (permutedCaseId)
  {
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    {
      update_interface_signs(interface_key, -1);
    }
    break;

    case 76: // ls[0]=0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    case 79: // ls[0]=0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      update_interface_signs(interface_key, +1);
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 78: // ls[0]<0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      const int sign = (permutedCaseId==2) ? -1 : 1;

      lnodes[i4] = SubElementNode::common_child({lnodes[i0], lnodes[i1]});
      lnodes[i6] = SubElementNode::common_child({lnodes[i0], lnodes[i2]});
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i4] && nullptr != lnodes[i6] && nullptr != lnodes[i7]);

      // face0: true: connect 4 and 3, false: connect 7 and 1
      // face2: true: connect 7 and 2, false: connect 6 and 3
      // face3: true: connect 6 and 1, false: connect 4 and 2
      const bool face0 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i1, i3, i4, i7);
      const bool face2 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i3, i2, i7, i6);
      const bool face3 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i2, i1, i6, i4);

      handle_tet( lnodes, subelement_interface_signs(interface_key, -sign), i6,i4,i7,i0, s3,s0,s2,-1 );
      handle_wedge( mesh, lnodes, subelement_interface_signs(interface_key, sign), i4,i7,i6,i1,i3,i2, s0,s2,s3,-1,s1, face0,face2,!face3 );
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 75: // ls[0]<0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    {
      const int sign = (permutedCaseId==5) ? -1 : 1;

      lnodes[i6] = SubElementNode::common_child({lnodes[i0], lnodes[i2]});
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i6] && nullptr != lnodes[i7]);

      // face2: true: connect 7 and 2, false: connect 6 and 3
      const bool face2 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i3, i2, i7, i6);

      handle_tet( lnodes, subelement_interface_signs(interface_key, -sign), i6,i1,i7,i0, s3,s0,s2,-1 );
      handle_pyramid( lnodes, subelement_interface_signs(interface_key, sign), i7,i6,i2,i3,i1, -1,s3,s1,s0,s2, face2 );
    }
    break;

    case 8:  // ls[0]>0 && ls[1]>0 && ls[2]<0 && ls[3]<0
    {
      lnodes[i5] = SubElementNode::common_child({lnodes[i1], lnodes[i2]});
      lnodes[i6] = SubElementNode::common_child({lnodes[i0], lnodes[i2]});
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      lnodes[i8] = SubElementNode::common_child({lnodes[i1], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i5] && nullptr != lnodes[i6] && nullptr != lnodes[i7] && nullptr != lnodes[i8]);

      // face0: true: connect 7 and 1, false: connect 8 and 0
      // face1: true: connect 5 and 3, false: connect 8 and 2
      // face2: true: connect 7 and 2, false: connect 6 and 3
      // face3: true: connect 5 and 0, false: connect 6 and 1
      const bool face0 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i3, i0, i1, i7, i8);
      const bool face1 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i1, i2, i3, i5, i8);
      const bool face2 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i0, i3, i2, i7, i6);
      const bool face3 = determine_diagonal_for_cut_triangular_face(simplexMethod, globalIDsAreParallelConsistent, lnodes, i2, i1, i0, i5, i6);

      // face 4: true: connect 6 and 8, false: connect 7 and 5
      const bool face4 = ElementObj::determine_diagonal_for_internal_quad_of_cut_tet_from_edge_nodes(simplexMethod, lnodes[i8], lnodes[i5], lnodes[i6], lnodes[i7],
          face0, face1, face2, face3);

      handle_wedge( mesh, lnodes, subelement_interface_signs(interface_key, -1), i8,i3,i7,i5,i2,i6, s1,s2,-1,s0,s3, !face1,!face2,face4 );
      handle_wedge( mesh, lnodes, subelement_interface_signs(interface_key,  1), i8,i1,i5,i7,i0,i6, s0,s3,-1,s1,s2, !face0,!face3,face4 );
    }
    break;

    case 13: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    case 67: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]>0
    {
      const int sign = (permutedCaseId==13) ? -1 : 1;
      update_interface_signs(interface_key, sign);
    }
    break;

    case 14: // ls[0]>0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    {
      lnodes[i7] = SubElementNode::common_child({lnodes[i0], lnodes[i3]});
      STK_ThrowRequire(nullptr != lnodes[i7]);

      handle_tet( lnodes, subelement_interface_signs(interface_key,  1), i1,i7,i2,i0, s0,s2,s3,-1 );
      handle_tet( lnodes, subelement_interface_signs(interface_key, -1), i2,i7,i1,i3, s2,s0,s1,-1 );
    }
    break;

    default: ThrowRuntimeError("Subelement decomposition error.");
  }

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    debug_subelements(lnodes, interface_key, caseId);
  }
}

bool SubElement_Tet_4::determine_diagonal_for_cut_triangular_face(const Simplex_Generation_Method & simplexMethod, const bool globalIDsAreParallelConsistent, const NodeVec & lnodes, const int /*i0*/, const int i1, const int i2, const int /*i3*/, const int /*i5*/)
{
  /*
   *     2 o
   *      / \
   *   5 o   \
   *    / \   \
   *   o---o---o
   *   0   3   1
   */

  // true:  connect nodes 3 and 2
  // false: connect nodes 5 and 1

  if (simplexMethod == CUT_QUADS_BY_GLOBAL_IDENTIFIER)
  {
    return lnodes[i2]->get_ancestry() < lnodes[i1]->get_ancestry();
  }
  else
  {
    STK_ThrowAssert(simplexMethod == CUT_QUADS_BY_DEFAULT_METHOD  || simplexMethod == CUT_QUADS_BY_LARGEST_ANGLE || simplexMethod == CUT_QUADS_BY_NEAREST_EDGE_CUT);
    return SubElementNode::higher_priority_by_score_then_ancestry(*lnodes[i2],*lnodes[i1], globalIDsAreParallelConsistent);
  }
}

void
SubElement_Tet_4::handle_pyramid( NodeVec & lnodes,
                               const std::vector<int> & subInterfaceSigns,
                               const int i0, const int i1, const int i2, const int i3, const int i4,
                               const int s0, const int s1, const int s2, const int s3, const int s4,
                               const bool face4 )
{
  if (face4)
  {
    handle_tet( lnodes, subInterfaceSigns, i0,i1,i2,i4, s0,s1,-1,s4 );
    handle_tet( lnodes, subInterfaceSigns, i2,i3,i0,i4, s2,s3,-1,s4 );
  }
  else
  {
    handle_tet( lnodes, subInterfaceSigns, i0,i1,i3,i4, s0,-1,s3,s4 );
    handle_tet( lnodes, subInterfaceSigns, i2,i3,i1,i4, s2,-1,s1,s4 );
  }
}

void
SubElement_Tet_4::handle_tet( NodeVec & lnodes,
                               const std::vector<int> & subInterfaceSigns,
                               const int i0, const int i1, const int i2, const int i3,
                               const int s0, const int s1, const int s2, const int s3)
{
  STK_ThrowAssert(!is_degenerate(lnodes, i0,i1,i2,i3));

  NodeVec sub_nodes(4,(SubElementNode *)NULL);
  std::vector<int> sub_parent_ids(4);

  sub_nodes[0] = lnodes[i0];
  sub_nodes[1] = lnodes[i1];
  sub_nodes[2] = lnodes[i2];
  sub_nodes[3] = lnodes[i3];
  sub_parent_ids[0] = (s0>= 0) ? my_parent_side_ids[s0] : -1;
  sub_parent_ids[1] = (s1>= 0) ? my_parent_side_ids[s1] : -1;
  sub_parent_ids[2] = (s2>= 0) ? my_parent_side_ids[s2] : -1;
  sub_parent_ids[3] = (s3>= 0) ? my_parent_side_ids[s3] : -1;

  std::unique_ptr<SubElement> sub = std::make_unique<SubElement_Tet_4>( sub_nodes, sub_parent_ids, my_owner);
  sub->set_interface_signs(subInterfaceSigns);
  add_subelement( std::move(sub) );
}

void
SubElement_Tet_4::handle_wedge( CDMesh & mesh,
    NodeVec & lnodes,
    const std::vector<int> & subInterfaceSigns,
    const int i0, const int i1, const int i2, const int i3, const int i4, const int i5,
    const int s0, const int s1, const int s2, const int s3, const int s4,
    const bool face0, const bool face1, const bool face2 )
{
/*
 *                                   PARENT Linear 6-Node Wedge Nodes
 *                       5           (SPACE_DIM=3!)
 *                     . o
 *                    . / \
 *                   . /   \         Face_Quad_4_3D()  0-1-4-3
 *                  . /     \        Face_Quad_4_3D()  1-2-5-4
 *                 . /       \       Face_Quad_4_3D()  0-3-5-2
 *              2 . o---------o 4    Face_Tri_3_3D()   0-2-1
 *               o . 3      .        Face_Tri_3_3D()   3-4-5
 *              /.\        .
 *             /.  \      .
 *            /.    \    .
 *           /.      \  .
 *          o---------o
 *         0           1
 *
 */
 // face0: true: connect 0 and 4, false: connect 1 and 3
 // face1: true: connect 1 and 5, false: connect 2 and 4
 // face2: true: connect 0 and 5, false: connect 2 and 3

  std::vector<int> wedge_nodes(6);
  wedge_nodes[0] = i0;
  wedge_nodes[1] = i1;
  wedge_nodes[2] = i2;
  wedge_nodes[3] = i3;
  wedge_nodes[4] = i4;
  wedge_nodes[5] = i5;

  std::vector<int> wedge_sides(5);
  wedge_sides[0] = s0;
  wedge_sides[1] = s1;
  wedge_sides[2] = s2;
  wedge_sides[3] = s3;
  wedge_sides[4] = s4;

  // Each of the 3 quad faces can be subdivided in 2 ways, giving a total of 8 possible decompositions.
  // 2 of these are illegal, however, since they don't result in tetrahedra.

  static const unsigned tet_nodes0[] = { 0,5,4,3, 0,2,1,5, 0,1,4,5 };
  static const unsigned tet_nodes1[] = { 1,3,5,4, 1,0,2,5, 1,0,5,3 };
  static const unsigned tet_nodes2[] = { 0,5,4,3, 0,2,1,4, 0,2,4,5 };
  static const unsigned tet_nodes3[] = { 0,0,0,0, 0,0,0,0, 0,0,0,0 }; // illegal
  static const unsigned tet_nodes4[] = { 0,0,0,0, 0,0,0,0, 0,0,0,0 }; // illegal
  static const unsigned tet_nodes5[] = { 1,3,5,4, 1,0,2,3, 1,2,5,3 };
  static const unsigned tet_nodes6[] = { 2,4,3,5, 2,1,0,4, 2,0,3,4 };
  static const unsigned tet_nodes7[] = { 2,4,3,5, 2,1,0,3, 2,1,3,4 };
  static const unsigned * tet_node_map[] = { tet_nodes0 , tet_nodes1 , tet_nodes2, tet_nodes3 , tet_nodes4 , tet_nodes5, tet_nodes6, tet_nodes7 };

  static const int tet_sides0[] = { 2, 4, 0,-1,  2, 1,-1, 3, -1, 1,-1, 0 };
  static const int tet_sides1[] = { 0, 4, 1,-1, -1, 2, 1, 3,  0, 2,-1,-1 };
  static const int tet_sides2[] = { 2, 4, 0,-1, -1, 1, 0, 3,  2, 1,-1,-1 };
  static const int tet_sides3[] = {-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 }; // illegal
  static const int tet_sides4[] = {-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 }; // illegal
  static const int tet_sides5[] = { 0, 4, 1,-1,  0, 2,-1, 3, -1, 2,-1, 1 };
  static const int tet_sides6[] = { 1, 4, 2,-1,  1, 0,-1, 3, -1, 0,-1, 2 };
  static const int tet_sides7[] = { 1, 4, 2,-1, -1, 0, 2, 3,  1, 0,-1,-1 };
  static const int * tet_side_map[] = { tet_sides0 , tet_sides1 , tet_sides2 , tet_sides3 , tet_sides4 , tet_sides5, tet_sides6, tet_sides7 };

  const unsigned case_id =
      (face0 ? 0 : 1) +
      (face1 ? 0 : 2) +
      (face2 ? 0 : 4);

  //krinolog << "Wedge case_id = " << case_id << stk::diag::dendl;

  if (case_id < 3 || case_id > 4)
  {
    const unsigned * tet_nodes = tet_node_map[case_id];
    const int * tet_sides = tet_side_map[case_id];

    unsigned lnn[12];
    int lsn[12];

    for (int n=0; n<12; ++n)
    {
      lnn[n] = wedge_nodes[tet_nodes[n]];
      lsn[n] = (tet_sides[n]<0) ? -1 : wedge_sides[tet_sides[n]];
    }

    handle_tet(lnodes, subInterfaceSigns, lnn[0], lnn[1], lnn[2], lnn[3],   lsn[0], lsn[1], lsn[2], lsn[3]);
    handle_tet(lnodes, subInterfaceSigns, lnn[4], lnn[5], lnn[6], lnn[7],   lsn[4], lsn[5], lsn[6], lsn[7]);
    handle_tet(lnodes, subInterfaceSigns, lnn[8], lnn[9], lnn[10], lnn[11], lsn[8], lsn[9], lsn[10], lsn[11]);
  }
  else
  {
    STK_ThrowRequire(case_id == 3 || case_id == 4);
    if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "Schonhardt's polyhedron formed, Adding Steiner point." << "\n";

    // Schonhardt's polyhedra should never be forced now that diagonals are cut using a globally consistent nodal criterion
    STK_ThrowRequireMsg(false, "Schonhardt polyhedron found.  This should not happen.");

    NodeVec wedge_nodevec(6);
    wedge_nodevec[0] = lnodes[i0];
    wedge_nodevec[1] = lnodes[i1];
    wedge_nodevec[2] = lnodes[i2];
    wedge_nodevec[3] = lnodes[i3];
    wedge_nodevec[4] = lnodes[i4];
    wedge_nodevec[5] = lnodes[i5];

    const stk::math::Vector3d & v0 = lnodes[i0]->owner_coords(my_owner);
    const stk::math::Vector3d & v1 = lnodes[i1]->owner_coords(my_owner);
    const stk::math::Vector3d & v2 = lnodes[i2]->owner_coords(my_owner);
    const stk::math::Vector3d & v3 = lnodes[i3]->owner_coords(my_owner);
    const stk::math::Vector3d & v4 = lnodes[i4]->owner_coords(my_owner);
    const stk::math::Vector3d & v5 = lnodes[i5]->owner_coords(my_owner);

    std::vector<double> weights(6);

    // Use the centroid of 1 of the elements that will be within the volume of the wedge depending
    // on which way that face2 is cut (which is on the internal quad, which may be non-planar).
    // This is not guaranteed, however, to always produce tets that stay within the deformed wedge.
    if (face2)
    {
      weights[0] = 1./4.;
      weights[1] = 1./4.;
      weights[4] = 1./4.;
      weights[5] = 1./4.;
    }
    else
    {
      weights[1] = 1./4.;
      weights[2] = 1./4.;
      weights[3] = 1./4.;
      weights[4] = 1./4.;
    }

    const SubElementNode * centroid = mesh.create_steiner_node( my_owner, wedge_nodevec, weights );

    const int i6 = lnodes.size();
    lnodes.push_back(centroid);
    const stk::math::Vector3d & v6 = lnodes[i6]->owner_coords(my_owner);

    // create 8 tets
    handle_tet(lnodes, subInterfaceSigns, i0, i2, i1, i6, -1, -1, -1, s3);
    handle_tet(lnodes, subInterfaceSigns, i3, i4, i5, i6, -1, -1, -1, s4);

    if (tet_volume({{v0, v2, v1, v6}}) < 0. || tet_volume({{v3, v4, v5, v6}}) < 0.)
    {
      krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
    }

    if (face0)
    {
      handle_tet(lnodes, subInterfaceSigns, i0, i1, i4, i6, -1, -1, -1, s0);
      handle_tet(lnodes, subInterfaceSigns, i0, i4, i3, i6, -1, -1, -1, s0);

      if (tet_volume({{v0, v1, v4, v6}}) < 0. || tet_volume({{v0, v4, v3, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }
    else
    {
      handle_tet(lnodes, subInterfaceSigns, i0, i1, i3, i6, -1, -1, -1, s0);
      handle_tet(lnodes, subInterfaceSigns, i1, i4, i3, i6, -1, -1, -1, s0);

      if (tet_volume({{v0, v1, v3, v6}}) < 0. || tet_volume({{v1, v4, v3, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }

    if (face1)
    {
      handle_tet(lnodes, subInterfaceSigns, i1, i2, i5, i6, -1, -1, -1, s1);
      handle_tet(lnodes, subInterfaceSigns, i1, i5, i4, i6, -1, -1, -1, s1);

      if (tet_volume({{v1, v2, v5, v6}}) < 0. || tet_volume({{v1, v5, v4, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }
    else
    {
      handle_tet(lnodes, subInterfaceSigns, i1, i2, i4, i6, -1, -1, -1, s1);
      handle_tet(lnodes, subInterfaceSigns, i2, i5, i4, i6, -1, -1, -1, s1);

      if (tet_volume({{v1, v2, v4, v6}}) < 0. || tet_volume({{v2, v5, v4, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }

    if (face2)
    {
      handle_tet(lnodes, subInterfaceSigns, i0, i3, i5, i6, -1, -1, -1, s2);
      handle_tet(lnodes, subInterfaceSigns, i0, i5, i2, i6, -1, -1, -1, s2);

      if (tet_volume({{v0, v3, v5, v6}}) < 0. || tet_volume({{v0, v5, v2, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }
    else
    {
      handle_tet(lnodes, subInterfaceSigns, i0, i3, i2, i6, -1, -1, -1, s2);
      handle_tet(lnodes, subInterfaceSigns, i3, i5, i2, i6, -1, -1, -1, s2);

      if (tet_volume({{v0, v3, v2, v6}}) < 0. || tet_volume({{v3, v5, v2, v6}}) < 0.)
      {
        krinolog << "Error: Schonhard polyhedron decomposition includes inverted tets." << stk::diag::dendl;
      }
    }
  }
}

double
SubElement::parametric_distance( const SubElementNode * node0, const SubElementNode * node1 )
{
  if (node0->is_mesh_node() && node1->is_mesh_node())
    return 1.0;

  std::map<const krino::SubElementNode *, double> node0_stencil;
  std::map<const krino::SubElementNode *, double> node1_stencil;
  node0->build_stencil(node0_stencil);
  node1->build_stencil(node1_stencil);

  // parametric distance = sqrt(0.5*sum(dx_i^2))

  double sum_sqr_dist = 0.;
  for (auto && entry : node0_stencil)
  {
    const SubElementNode * parent = entry.first;
    const double wt0 = entry.second;

    auto it = node1_stencil.find(parent);
    const double wt1 = (it == node1_stencil.end()) ? 0.0 : (it->second);

    sum_sqr_dist += (wt0 - wt1)*(wt0 - wt1);
  }
  for (auto && entry : node1_stencil)
  {
    const SubElementNode * parent = entry.first;
    const double wt1 = entry.second;

    auto it = node0_stencil.find(parent);
    if (it == node0_stencil.end())
    {
      sum_sqr_dist += wt1*wt1;
    }
  }

  return std::sqrt(0.5*sum_sqr_dist);
}

void set_node_signs_for_edge(const InterfaceID /*interface*/, const SubElementNode * node1, const SubElementNode * node2, const int node1Sign, const int node2Sign, const double crossingLocation, const CDFEM_Snapper & snapper)
{
  if (node1Sign != node2Sign)
  {
    // determine tolerance for this edge that may have already been cut
    const bool always_snap = snapper.always_snap();
    const double edge_tol = always_snap ? 0.0 : snapper.get_edge_tolerance()/SubElement::parametric_distance(node1, node2);

    if (always_snap || edge_tol > 0.5)
    {
      if (crossingLocation < 0.5)
      {
        node1->set_node_sign(0);
      }
      else
      {
        node2->set_node_sign(0);
      }
    }
    else
    {
      if (crossingLocation < edge_tol)
      {
        node1->set_node_sign(0);
      }
      else if (crossingLocation > 1.-edge_tol)
      {
        node2->set_node_sign(0);
      }
    }
  }

  node1->set_node_sign(node1Sign);
  node2->set_node_sign(node2Sign);
}

bool
SubElement::determine_node_signs_on_edge( const CDMesh & mesh, const InterfaceID interface, const int i0, const int i1 )
{
  bool setSignOnEdge = false;
  if (my_owner->have_interface(interface))
  {
    const SubElementNode * parent1 = my_nodes[i0];
    const SubElementNode * parent2 = my_nodes[i1];

    if (parent1->captures_interface(interface))
      parent1->set_node_sign(0);
    if (parent2->captures_interface(interface))
      parent2->set_node_sign(0);

    const int sign = myInterfaceSigns[my_owner->get_interface_index(interface)];
    if (sign == 0)
    {
      const auto [crossingSign, position] = my_owner->interface_edge_crossing_sign_and_position(interface, parent1, parent2);
      if (crossingSign != 0)
      {
        set_node_signs_for_edge(interface, parent1, parent2, -crossingSign, crossingSign, position, mesh.get_snapper());
        setSignOnEdge = true;
      }
    }
  }
  return setSignOnEdge;
}

void SubElement::set_node_signs_on_uncrossed_subelement( const InterfaceID interface )
{
  std::vector<stk::math::Vector3d> nodeOwnerCoords;
  fill_node_owner_coords(my_owner, nodeOwnerCoords);
  const int interfaceSign = my_owner->get_interface_sign_for_uncrossed_subelement(interface, nodeOwnerCoords);
  for (auto && node : my_nodes)
    node->set_node_sign(interfaceSign);
}

void
SubElement::process_edge( CDMesh & mesh, const InterfaceID interface, const int i0, const int i1 )
{
  const SubElementNode * parent1 = my_nodes[i0];
  const SubElementNode * parent2 = my_nodes[i1];

  if (parent1->get_node_on_interface() || parent2->get_node_on_interface()) return;

  const SubElementNode * subnode = SubElementNode::common_child({parent1, parent2});

  if (subnode == nullptr &&
      parent1->node_sign_is_set() &&
      parent2->node_sign_is_set() &&
      parent1->get_node_sign() == -parent2->get_node_sign() &&
      have_interface(interface))
  {
    const auto [crossingSign, position] = my_owner->interface_edge_crossing_sign_and_position(interface, parent1, parent2);
    STK_ThrowRequireMsg(crossingSign!=0 && position > 0. && position < 1., "Error process_edge " << crossingSign << " " << position << " " << parent1->get_node_sign() << " " << parent2->get_node_sign() << " " << parent1->get_ancestry() << " " << parent2->get_ancestry());

    std::unique_ptr<SubElementNode> newNode = std::make_unique<SubElementEdgeNode>(my_owner, position, parent1, parent2);
    subnode = mesh.add_managed_node(std::move(newNode));

    // add subnode as child to parents
    parent1->add_child(subnode);
    parent2->add_child(subnode);
  }
}

void
SubElement::determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->SubElement::determine_node_signs(mesh, interface_key);
    }
    return;
  }

  determine_node_signs(mesh, interface_key);
}

void
SubElement::determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->SubElement::determine_node_scores(mesh, interface_key);
    }
    return;
  }

  determine_node_scores(mesh, interface_key);
}

void
SubElement::decompose(CDMesh & mesh, const InterfaceID interface_key)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->decompose(mesh, interface_key);
    }
    return;
  }

  decompose_edges(mesh, interface_key);
}

std::vector<int>
SubElement::get_edges_with_children(const InterfaceID & /*interface*/) const
{
  // Iterate edges looking for any common children of the edge nodes
  const stk::topology Top = topology();
  const int num_edges = Top.num_edges();
  std::vector<int> edgesWithChildren;

  for ( int edge = 0; edge < num_edges; ++edge )
  {
    const unsigned * edge_node_ordinals = get_edge_node_ordinals(Top, edge);

    const SubElementNode * node0 = my_nodes[edge_node_ordinals[0]];
    const SubElementNode * node1 = my_nodes[edge_node_ordinals[1]];
    const SubElementNode * child = SubElementNode::common_child({node0, node1});
    if( child )
    {
      edgesWithChildren.reserve(num_edges);
      if(krinolog.shouldPrint(LOG_DEBUG))
      {
        krinolog << "Found hanging node on edge " << edge << " of element id=" << entityId() << "\n";
        krinolog << " Ancestry: " << child->get_ancestry() << "\n";
      }
      edgesWithChildren.push_back(edge);
    }
  }
  return edgesWithChildren;
}

void
SubElement::handle_hanging_children(CDMesh & mesh, const InterfaceID & interface)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->handle_hanging_children(mesh, interface);
    }
    return;
  }

  const std::vector<int> edgesWithChildren = get_edges_with_children(interface);

  fix_hanging_children(mesh, interface, edgesWithChildren);

  for ( auto && subelem : my_subelements )
  {
    STK_ThrowRequire(!subelem->have_edges_with_children());
  }
}

void
SubElement::build_quadratic_subelements(CDMesh & mesh)
{
  STK_ThrowRequireMsg(have_subelements(), "This subelement type does not support quadratic subelements.");
  for ( auto && subelem : my_subelements )
  {
    subelem->build_quadratic_subelements(mesh);
  }
}

void
SubElement::cut_face_interior_intersection_points(CDMesh & /*mesh*/, int /*level*/)
{

  if (topology().num_faces() == 0)
    return;

  STK_ThrowRequireMsg(false, "This subelement type does not support cut_face_interior_intersection_points: " << topology());
}

bool
SubElement_Tet_4::is_degenerate( NodeVec & lnodes,
                                  const int i0, const int i1, const int i2, const int i3 )
{

  if ( lnodes[i0] == lnodes[i1] ||
       lnodes[i0] == lnodes[i2] ||
       lnodes[i0] == lnodes[i3] ||
       lnodes[i1] == lnodes[i2] ||
       lnodes[i1] == lnodes[i3] ||
       lnodes[i2] == lnodes[i3] )
  {
    // this tet is degenerate with two coincident nodes
    return true;
  }

  return false;
}

} // namespace krino
