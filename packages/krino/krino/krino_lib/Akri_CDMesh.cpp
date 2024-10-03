// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDMesh.hpp>

#include <stk_mesh/base/CommunicateMeshTypes.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/environment/LogWithTimeAndMemory.hpp>
#include <stk_util/parallel/ParallelReduce.hpp> // Needed for all_reduce_max
#include <stk_util/parallel/CommSparse.hpp>

#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_ChildNodeCreator.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_Element.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDMesh_Debug.hpp>
#include <Akri_CDMesh_Refinement.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_ChildNodeStencil.hpp>
#include <Akri_ConformingPhaseParts.hpp>
#include <Akri_ReportHandler.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>
#include <Akri_DecompositionHasChanged.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_Snap.hpp>
#include <Akri_SnapToNode.hpp>
#include <Akri_SubElementChildNodeAncestry.hpp>
#include <Akri_Surface_Identifier.hpp>
#include <stk_math/StkVector.hpp>

#include <math.h>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <limits>

#include <Akri_RefinementSupport.hpp>
#include <Akri_ParallelErrorMessage.hpp>
namespace krino{

std::unique_ptr<CDMesh> CDMesh::the_new_mesh;

//--------------------------------------------------------------------------------

CDMesh::CDMesh( stk::mesh::BulkData & mesh )
  : my_meta(mesh.mesh_meta_data()),
    my_aux_meta(AuxMetaData::get(my_meta)),
    my_entity_id_pool(my_meta),
    my_spatial_dim( my_meta.spatial_dimension() ),
    my_cdfem_support(CDFEM_Support::get(my_meta)),
    my_phase_support(Phase_Support::get(my_meta)),
    myRefinementSupport(RefinementSupport::get(my_meta)),
    my_stash_step_count(-1),
    my_missing_remote_prolong_facets(false),
    my_timer_decompose("Decompose", my_cdfem_support.get_timer_cdfem()),
    my_timer_decomposition_has_changed("Need CDFEM", my_cdfem_support.get_timer_cdfem()),
    my_timer_snap("Snapping", my_timer_decompose),
    my_timer_stash_field_data("Stash Field Data", my_timer_decompose),
    my_timer_modify_mesh("Modify Mesh", my_cdfem_support.get_timer_cdfem()),
    my_timer_prolongation("Prolongation", my_cdfem_support.get_timer_cdfem()),
    my_timer_compute_CFL("Compute CFL", my_cdfem_support.get_timer_cdfem())
{
  stk::mesh::insert(my_attribute_parts, my_aux_meta.active_part());
  stk::mesh::insert(my_attribute_parts, my_aux_meta.exposed_boundary_part());
  stk::mesh::insert(my_attribute_parts, my_aux_meta.block_boundary_part());
}

CDMesh::~CDMesh()
{
  clear();
}

SubElementNode * CDMesh::add_managed_node(std::unique_ptr<SubElementNode> node)
{
  nodes.emplace_back(std::move(node));
  return nodes.back().get();
}

//--------------------------------------------------------------------------------


stk::mesh::Part & CDMesh::get_locally_owned_part() const { return my_meta.locally_owned_part(); }
stk::mesh::Part & CDMesh::get_globally_shared_part() const { return my_meta.globally_shared_part(); }
stk::mesh::Part & CDMesh::get_active_part() const { return my_aux_meta.active_part(); }
stk::mesh::Part & CDMesh::get_block_boundary_part() const { return my_aux_meta.block_boundary_part(); }

//--------------------------------------------------------------------------------

void CDMesh::add_periodic_node_pair(stk::mesh::Entity node1, stk::mesh::Entity node2)
{
  my_periodic_node_id_map[stk_bulk().identifier(node1)].push_back(stk_bulk().identifier(node2));
  my_periodic_node_id_map[stk_bulk().identifier(node2)].push_back(stk_bulk().identifier(node1));
}

const std::vector<InterfaceID> & CDMesh::all_interface_ids(const std::vector<Surface_Identifier> & surfaceIdentifiers) const
{
  if(crossing_keys.empty())
  {
    const int numSurfaces = surfaceIdentifiers.size();

    if(numSurfaces < 2 || !my_phase_support.has_one_levelset_per_phase())
    {
      crossing_keys.resize(numSurfaces);
      for(int i=0; i < numSurfaces; ++i)
      {
        crossing_keys[i] = InterfaceID(i,i);
      }
    }
    else
    {
      for(int i=0; i < numSurfaces; ++i)
      {
        for(int j=i+1; j < numSurfaces; ++j)
        {
          crossing_keys.push_back(InterfaceID(i,j));
        }
      }
    }
  }
  return crossing_keys;
}

std::vector<InterfaceID> CDMesh::active_interface_ids(const std::vector<Surface_Identifier> & surfaceIdentifiers) const
{
  const std::vector<InterfaceID> all_interfaces = all_interface_ids(surfaceIdentifiers);
  if (all_interfaces.size() == 1) return all_interfaces;

  std::vector<int> id_is_active_locally(all_interfaces.size(), false);
  for (const auto & elem : elements)
  {
    for (auto && elemInterface : elem->get_sorted_cutting_interfaces())
    {
      const auto lower = std::lower_bound(all_interfaces.begin(), all_interfaces.end(), elemInterface);
      STK_ThrowAssert(*lower == elemInterface);
      id_is_active_locally[std::distance(all_interfaces.begin(), lower)] = true;
    }
  }

  std::vector<int> id_is_active_globally(all_interfaces.size());
  stk::all_reduce_sum(stk_bulk().parallel(), id_is_active_locally.data(), id_is_active_globally.data(), id_is_active_locally.size());

  std::vector<InterfaceID> active_ids;
  for (size_t id=0; id<id_is_active_globally.size(); ++id)
    if (id_is_active_globally[id])
      active_ids.push_back(all_interfaces[id]);

  return active_ids;
}

void
CDMesh::handle_possible_failed_time_step( stk::mesh::BulkData & mesh, const int step_count )
{
  const bool haveEverPerformedDecomposition = nullptr != the_new_mesh.get();
  if (haveEverPerformedDecomposition)
  {
    the_new_mesh->rebuild_after_rebalance_or_failed_step();
  }
}

static void interpolate_nodal_field(const FieldRef field,
    stk::mesh::Entity node,
    const std::vector<stk::mesh::Entity> & interpNodes,
    const std::vector<double> & interpWeights)
{
  const unsigned fieldLength = field.length();

  double * val = field_data<double>(field, node);
  if (nullptr == val) return;

  for (unsigned i=0; i<fieldLength; ++i)
    val[i] = 0.;

  for (size_t iNode=0; iNode<interpNodes.size(); ++iNode)
  {
    const double * nodeVal = field_data<double>(field, interpNodes[iNode]);
    STK_ThrowRequire(nullptr != nodeVal);

    for (unsigned i=0; i<fieldLength; ++i)
      val[i] += interpWeights[iNode] * nodeVal[i];
  }
}

static bool any_node_was_snapped(const std::vector<stk::mesh::Entity> & nodes,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  for (auto && node : nodes)
  {
    const auto iter = nodesToCapturedDomains.find(node);
    if (iter != nodesToCapturedDomains.end() && !iter->second.empty())
      return true;
  }
  return false;
}

const SubElementNode *
CDMesh::find_new_node_with_common_ancestry_as_existing_node_with_given_id(const stk::mesh::EntityId nodeId) const
{
  const SubElementNode * cdmeshNode = get_mesh_node(nodeId);
  if (cdmeshNode != nullptr)
    return cdmeshNode;

  const stk::mesh::Entity node = stk_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
  STK_ThrowAssert(stk_bulk().is_valid(node));
  if (stk_bulk().bucket(node).member(get_child_edge_node_part()))
    return find_new_node_with_common_ancestry_as_existing_child_node(node);

  return nullptr;
}

const SubElementNode *
CDMesh::find_new_node_with_common_ancestry_as_existing_child_node(const stk::mesh::Entity node) const
{
  STK_ThrowAssert(stk_bulk().bucket(node).member(get_child_edge_node_part()));

  const auto parentIds = get_child_node_parent_ids(stk_bulk(), get_parent_node_ids_field(), node);
  std::vector<const SubElementNode *> parents;
  parents.reserve(parentIds.size());
  for (auto && parentId : parentIds)
  {
    const SubElementNode * parent = find_new_node_with_common_ancestry_as_existing_node_with_given_id(parentId);
    if (nullptr == parent)
      return nullptr;
    parents.push_back(parent);
  }
  return SubElementNode::common_child(parents);
}

static void apply_snapping_to_children_of_snapped_nodes(const std::vector<ChildNodeStencil> & childNodeStencils,
    const FieldSet & snapFields,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  for(const auto & childNodeStencil : childNodeStencils)
    if (any_node_was_snapped(childNodeStencil.parentNodes, nodesToCapturedDomains))
      for (auto && field : snapFields)
        interpolate_nodal_field(field, childNodeStencil.childNode, childNodeStencil.parentNodes, childNodeStencil.parentWeights);
}

static void undo_previous_snaps_without_interpolation(const stk::mesh::BulkData & mesh, const FieldRef coordsField, FieldRef cdfemSnapField)
{
  FieldRef oldSnapDisplacements = cdfemSnapField.field_state(stk::mesh::StateOld);
  stk::mesh::field_axpby(-1.0, oldSnapDisplacements, +1.0, coordsField);
  FieldRef modelCoords(mesh.mesh_meta_data().coordinate_field());
  if (modelCoords != coordsField)
    stk::mesh::field_axpby(-1.0, oldSnapDisplacements, +1.0, modelCoords);
}

void CDMesh::prepare_for_resnapping(const stk::mesh::BulkData & mesh, const InterfaceGeometry & interfaceGeometry)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  FieldRef cdfemSnapField = cdfemSupport.get_cdfem_snap_displacements_field();

  if (cdfemSnapField.valid())
  {
    const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());

    if (cdfemSupport.get_resnap_method() == krino::RESNAP_AFTER_USING_INTERPOLATION_TO_UNSNAP)
    {
      undo_previous_snaps_using_interpolation(mesh, auxMeta.active_part(), cdfemSupport.get_coords_field(), cdfemSnapField, cdfemSupport.get_snap_fields());
    }
    else if (cdfemSupport.get_resnap_method() == krino::RESNAP_USING_INTERFACE_ON_PREVIOUS_SNAPPED_MESH ||
        cdfemSupport.get_resnap_method() == RESNAP_USING_INTERPOLATION)
    {
      if (cdfemSupport.get_resnap_method() == krino::RESNAP_USING_INTERFACE_ON_PREVIOUS_SNAPPED_MESH)
        cdfemSnapField.field().rotate_multistate_data();

      interfaceGeometry.set_do_update_geometry_when_mesh_changes(true);
      interfaceGeometry.prepare_to_intersect_elements(mesh);
      interfaceGeometry.set_do_update_geometry_when_mesh_changes(false);

      undo_previous_snaps_without_interpolation(mesh, cdfemSupport.get_coords_field(), cdfemSnapField);
    }
  }
}

void CDMesh::snap_and_update_fields_and_captured_domains(const InterfaceGeometry & interfaceGeometry,
    NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  stk::diag::TimeBlock timer__(my_timer_snap);

  const FieldSet snapFields = my_cdfem_support.get_snap_fields();

  FieldRef cdfemSnapField = my_cdfem_support.get_cdfem_snap_displacements_field();

  if (cdfemSnapField.valid())
    stk::mesh::field_copy(my_cdfem_support.get_coords_field(), cdfemSnapField);

  std::vector<ChildNodeStencil> childNodeStencils;
  if (cdfemSnapField.valid())
    fill_child_node_stencils(stk_bulk(), get_child_edge_node_part(), get_parent_node_ids_field(), get_parent_node_weights_field(), childNodeStencils);

  const double minIntPtWeightForEstimatingCutQuality = get_snapper().get_edge_tolerance();

  const stk::mesh::Selector parentElementSelector = get_cdfem_parent_element_selector(get_active_part(), my_cdfem_support, my_phase_support);
  nodesToCapturedDomains = snap_as_much_as_possible_while_maintaining_quality(stk_bulk(),
      parentElementSelector,
      snapFields,
      interfaceGeometry,
      my_cdfem_support.get_global_ids_are_parallel_consistent(),
      my_cdfem_support.get_snapping_sharp_feature_angle_in_degrees(),
      minIntPtWeightForEstimatingCutQuality,
      my_cdfem_support.get_max_edge_snap());

  if (cdfemSnapField.valid())
  {
    apply_snapping_to_children_of_snapped_nodes(childNodeStencils, snapFields, nodesToCapturedDomains);
    stk::mesh::field_axpby(+1.0, my_cdfem_support.get_coords_field(), -1.0, cdfemSnapField);
  }

  if (cdfemSnapField.valid() && my_cdfem_support.get_resnap_method() == krino::RESNAP_USING_INTERPOLATION)
  {
    const FieldSet & postSnapInterpFields = my_cdfem_support.get_interpolation_fields();
    snap_fields_using_interpolation(stk_bulk(), my_aux_meta.active_part(), my_cdfem_support.get_coords_field(), cdfemSnapField, postSnapInterpFields);
    apply_snapping_to_children_of_snapped_nodes(childNodeStencils, postSnapInterpFields, nodesToCapturedDomains);
  }
}

int
CDMesh::decompose_mesh(stk::mesh::BulkData & mesh,
      const InterfaceGeometry & interfaceGeometry,
      const int stepCount,
      const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>> & periodic_node_pairs)
{
  stk::diag::TimeBlock root_timer__(CDFEM_Support::get(mesh.mesh_meta_data()).get_timer_cdfem());

  stk::log_with_time_and_memory(mesh.parallel(), "Begin Mesh Decomposition.");

  CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());

  const bool wasPreviouslyDecomposed = nullptr != the_new_mesh;
  if (!wasPreviouslyDecomposed)
  {
    // FIXME: This can cause problems for shells.
    attach_sides_to_elements(mesh);
  }

  krinolog << "Decomposing mesh for region into phase conformal elements." << stk::diag::dendl;
  NodeToCapturedDomainsMap nodesToCapturedDomains;

  {
    the_new_mesh = std::make_unique<CDMesh>(mesh);

    fix_node_owners_to_assure_active_owned_element_for_node(mesh, the_new_mesh->get_active_part());

    for(auto && pair : periodic_node_pairs)
    {
      the_new_mesh->add_periodic_node_pair(pair.first, pair.second);
    }

    // Not sure if this is krino's responsibility or the driving application.  If we have
    // elemental death fields, these need to be parallel consistent on aura elements.
    the_new_mesh->parallel_communicate_elemental_death_fields();

    stk::diag::TimeBlock timer__(the_new_mesh->my_timer_decompose);

    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      the_new_mesh->snap_and_update_fields_and_captured_domains(interfaceGeometry, nodesToCapturedDomains);

    interfaceGeometry.prepare_to_decompose_elements(the_new_mesh->stk_bulk(), nodesToCapturedDomains);
  }

  {
    stk::diag::TimeBlock timer__(the_new_mesh->my_timer_decompose);

    the_new_mesh->generate_nonconformal_elements();
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      the_new_mesh->snap_nearby_intersections_to_nodes(interfaceGeometry, nodesToCapturedDomains);
    the_new_mesh->set_phase_of_uncut_elements(interfaceGeometry);
    the_new_mesh->triangulate(interfaceGeometry);
    the_new_mesh->decompose(interfaceGeometry);
  }

  const int stashStepCount = wasPreviouslyDecomposed ? stepCount : (-1);
  the_new_mesh->stash_field_data(stashStepCount);

  const bool mesh_modified = the_new_mesh->modify_mesh();

  the_new_mesh->prolongation();

  // debugging
  if ( krinolog.shouldPrint(LOG_DEBUG) )
  {
    the_new_mesh->debug_output();
  }

  {
    const ScaledJacobianQualityMetric qualityMetric;
    krinolog << "After cutting quality is " << compute_mesh_quality(mesh, the_new_mesh->get_active_part(), qualityMetric) << stk::diag::dendl;
  }


  if (!the_new_mesh->aux_meta().using_fmwk())
  {
    the_new_mesh->print_conformal_volumes_and_surface_areas();
  }

  stk::log_with_time_and_memory(mesh.parallel(), "End Mesh Decomposition.");

  const int status = mesh_modified ? (COORDINATES_MAY_BE_MODIFIED | MESH_MODIFIED) : COORDINATES_MAY_BE_MODIFIED;
  return status;
}

static void fill_subelement_node_parent_entities(const SubElementNode & node, std::vector<stk::mesh::Entity> & nodeParents)
{
  nodeParents.clear();
  for (auto * parent : node.get_parents())
    nodeParents.push_back(parent->entity());
}

void CDMesh::store_child_node_parent_ids_and_weights_fields() const
{
  if (get_parent_node_ids_field().valid() && get_parent_node_weights_field().valid())
  {
    std::vector<stk::mesh::Entity> nodeParents;
    for (auto && node : nodes)
    {
      if (nullptr != dynamic_cast<const SubElementChildNode *>(node.get()))
      {
        fill_subelement_node_parent_entities(*node, nodeParents);
        store_child_node_parent_ids_and_weights(stk_bulk(), get_parent_node_ids_field(), get_parent_node_weights_field(), node->entity(), nodeParents, node->get_parent_weights());
      }
    }
    const std::vector<const stk::mesh::FieldBase *> parentFields{&get_parent_node_ids_field().field(), &get_parent_node_weights_field().field()};
    stk::mesh::communicate_field_data(stk_bulk(), parentFields);
  }
}

bool
CDMesh::modify_mesh()
{
  stk::diag::TimeBlock timer__(my_timer_modify_mesh);

  ParallelThrowAssert(stk_bulk().parallel(), check_face_and_edge_ownership(stk_bulk()));
  ParallelThrowAssert(stk_bulk().parallel(), check_face_and_edge_relations(stk_bulk()));

  set_entities_for_child_nodes_with_common_ancestry_as_existing_child_nodes();
  const bool all_elems_are_set_and_correct = set_entities_for_existing_child_elements();

  std::vector< stk::mesh::Entity> ownedUnusedOldChildElems = get_owned_unused_old_child_elements_and_clear_child_elements();

  const bool modificationIsNeeded = (myRefinementSupport.get_interface_maximum_refinement_level() > 0) || stk::is_true_on_any_proc(stk_bulk().parallel(), !all_elems_are_set_and_correct || !ownedUnusedOldChildElems.empty());

  if (modificationIsNeeded)
  {
    stk::mesh::toggle_sideset_updaters(stk_bulk(), false);
    stk_bulk().modification_begin();
    create_node_entities();
    std::vector<SideDescription> side_requests;
    create_element_and_side_entities(side_requests);
    destroy_custom_ghostings(stk_bulk());
    stk::mesh::destroy_elements_no_mod_cycle(stk_bulk(), ownedUnusedOldChildElems, stk_meta().universal_part());
    stk_bulk().modification_end();
    ParallelThrowAssert(stk_bulk().parallel(), check_shared_entity_nodes(stk_bulk()));

    add_possible_interface_sides(side_requests);
    batch_create_sides(stk_bulk(), side_requests);

    stk::mesh::toggle_sideset_updaters(stk_bulk(), true);
    activate_selected_entities_touching_active_elements(stk_bulk(), stk::topology::NODE_RANK, stk_meta().universal_part(), aux_meta().active_part()); // we should be able to skip this step if there are no higher order elements
    update_element_side_parts();

    ParallelThrowAssert(stk_bulk().parallel(), check_element_side_connectivity(stk_bulk(), aux_meta().exposed_boundary_part(), aux_meta().active_part()));
    ParallelThrowAssert(stk_bulk().parallel(), check_element_side_parts());

    aux_meta().induce_topology_nodesets(aux_meta().active_locally_owned_selector());

    ParallelThrowAssert(stk_bulk().parallel(), check_induced_parts(stk_bulk()));
  }

  store_child_node_parent_ids_and_weights_fields();

  return modificationIsNeeded;
}

void
CDMesh::set_entities_for_child_nodes_with_common_ancestry_as_existing_child_nodes()
{
  if (!was_mesh_previously_decomposed()) return;

  const stk::mesh::BulkData & mesh = stk_bulk();
  const stk::mesh::Selector ownedOrSharedChildNodeSelector = get_child_edge_node_part() & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedOrSharedChildNodeSelector))
  {
    for(const auto childNode : *bucketPtr)
    {
      const SubElementNode * newNode = find_new_node_with_common_ancestry_as_existing_child_node(childNode);
      if (nullptr != newNode)
        newNode->set_entity(mesh, childNode);
    }
  }
}

void
CDMesh::parallel_communicate_elemental_death_fields() const
{
  std::vector< const stk::mesh::FieldBase *> element_fields;
  for (auto && field : my_cdfem_support.get_levelset_fields())
  {
    if (field.valid() && field.entity_rank() == stk::topology::ELEMENT_RANK)
    {
      element_fields.push_back(&field.field());
    }
  }
  stk::mesh::communicate_field_data(stk_bulk(), element_fields);
}

bool
CDMesh::set_entities_for_existing_child_elements()
{
  std::vector<stk::mesh::Entity> subelem_node_entities;
  std::vector<stk::mesh::Entity> existing_elems;

  bool all_element_entities_are_set_and_correct = true;
  for (const auto & elem : elements)
  {
    if (elem->have_subelements())
    {
      std::vector<const SubElement *> conformal_subelems;
      elem->get_subelements(conformal_subelems);

      for (auto && subelem : conformal_subelems)
      {
        // If all nodes are set, look for existing element using the nodes.
        subelem_node_entities.clear();
        for (auto&& subelem_node : subelem->get_nodes())
        {
          stk::mesh::Entity node = subelem_node->entity();
          if (stk_bulk().is_valid(node)) subelem_node_entities.push_back(node);
        }
        existing_elems.clear();
        if (subelem_node_entities.size() == subelem->get_nodes().size())
        {
          stk::mesh::get_entities_through_relations(stk_bulk(), subelem_node_entities, stk::topology::ELEMENT_RANK, existing_elems);
          STK_ThrowAssert(existing_elems.size() <= 1);
        }

        if (existing_elems.empty())
        {
          all_element_entities_are_set_and_correct = false;
        }
        else
        {
          subelem->set_entity(stk_bulk(), existing_elems[0]);
          STK_ThrowAssert(subelem->check_entity_nodes(stk_bulk()));
          if (all_element_entities_are_set_and_correct && elem_io_part_changed(*subelem)) all_element_entities_are_set_and_correct = false;
        }
      }
    }
    else
    {
      if (all_element_entities_are_set_and_correct && elem_io_part_changed(*elem)) all_element_entities_are_set_and_correct = false;
    }
  }
  return all_element_entities_are_set_and_correct;
}

std::vector<stk::mesh::Entity>
CDMesh::get_owned_unused_old_child_elements_and_clear_child_elements()
{
  std::vector<stk::mesh::Entity> ownedUnusedOldChildElems;
  stk::mesh::Selector selector = get_child_part() & get_locally_owned_part();
  std::vector<stk::mesh::Entity> old_child_elems;
  stk::mesh::get_selected_entities( selector, stk_bulk().buckets( stk::topology::ELEMENT_RANK ), old_child_elems );

  ownedUnusedOldChildElems.reserve(old_child_elems.size());

  for (auto&& old_child_elem : old_child_elems)
  {
    const SubElement * subelem = find_child_element(old_child_elem);
    if (subelem == nullptr)
    {
      ownedUnusedOldChildElems.push_back(old_child_elem);
    }
  }
  child_elements.clear(); // reset child element vector

  return ownedUnusedOldChildElems;
}

bool
CDMesh::decomposition_needs_update(const InterfaceGeometry & interfaceGeometry,
      const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>> & periodic_node_pairs)
{
  return !the_new_mesh || the_new_mesh->decomposition_has_changed(interfaceGeometry);
}

void
CDMesh::mark_interface_elements_for_adaptivity(stk::mesh::BulkData & mesh, const FieldRef coordsField, const RefinementSupport & refinementSupport, const InterfaceGeometry & interfaceGeometry, const int num_refinements)
{
  krino::mark_interface_elements_for_adaptivity(mesh,
      refinementSupport.get_non_interface_conforming_refinement(),
      interfaceGeometry,
      refinementSupport,
      coordsField,
      num_refinements);
}

void
CDMesh::nonconformal_adaptivity(stk::mesh::BulkData & mesh, const FieldRef coordsField, const InterfaceGeometry & interfaceGeometry)
{
  const auto & refinementSupport = RefinementSupport::get(mesh.mesh_meta_data());
  stk::diag::TimeBlock timer__(refinementSupport.get_timer());

  stk::log_with_time_and_memory(mesh.parallel(), "Begin Nonconformal Adaptivity.");

  auto & refinement = refinementSupport.get_non_interface_conforming_refinement();

  std::function<void(int)> markerFunction;
  if (refinementSupport.has_refinement_interval())
  {
    markerFunction = [&mesh, &refinementSupport, &interfaceGeometry](int num_refinements)
    {
      constexpr bool isDefaultCoarsen = true;
      mark_elements_that_intersect_interval(mesh,
          refinementSupport.get_non_interface_conforming_refinement(),
          interfaceGeometry,
          refinementSupport.get_refinement_interval(),
          refinementSupport.get_interface_minimum_refinement_level(),
          isDefaultCoarsen);
    };
  }
  else
  {
    markerFunction = [&mesh, &coordsField, &refinementSupport, &interfaceGeometry](int num_refinements)
    {
      mark_interface_elements_for_adaptivity(mesh, coordsField, refinementSupport, interfaceGeometry, num_refinements);
    };
  }

  perform_multilevel_adaptivity(refinement, mesh, markerFunction, refinementSupport.get_do_not_refine_or_unrefine_selector());

  stk::log_with_time_and_memory(mesh.parallel(), "End Nonconformal Adaptivity.");
}

void
CDMesh::rebuild_after_rebalance_or_failed_step()
{
  clear();
  generate_nonconformal_elements();
  restore_subelements();
}

static bool side_is_adaptivity_or_cdfem_parent(stk::mesh::BulkData & mesh, const RefinementSupport & refinementSupport, const stk::mesh::Entity side, const stk::mesh::Part & cdfemParentPart)
{
  if (refinementSupport.has_non_interface_conforming_refinement() && refinementSupport.get_non_interface_conforming_refinement().is_parent_side(side))
    return true;
  for (auto element : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
    if (mesh.bucket(element).member(cdfemParentPart))
      return true;
  return false;
}

static void delete_extraneous_inactive_sides(stk::mesh::BulkData & mesh, const RefinementSupport & refinementSupport, const stk::mesh::Part & cdfemParentPart, const stk::mesh::Part & activePart)
{
  stk::mesh::Selector notActive = !activePart;

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_selected_entities(notActive, mesh.buckets(mesh.mesh_meta_data().side_rank()), sides, false);

  mesh.modification_begin();

  for (auto && side : sides)
    if (!side_is_adaptivity_or_cdfem_parent(mesh, refinementSupport, side, cdfemParentPart))
      STK_ThrowRequireMsg(disconnect_and_destroy_entity(mesh, side), "Could not destroy entity " << mesh.entity_key(side));

  mesh.modification_end();
}

void
CDMesh::rebuild_from_restart_mesh(stk::mesh::BulkData & mesh)
{
  ParallelThrowRequire(mesh.parallel(), !the_new_mesh);

  the_new_mesh = std::make_unique<CDMesh>(mesh);
  the_new_mesh->rebuild_child_part();
  the_new_mesh->rebuild_parent_and_active_parts_using_nonconformal_and_child_parts();
  the_new_mesh->generate_nonconformal_elements();
  the_new_mesh->restore_subelements();

  activate_selected_entities_touching_active_elements(the_new_mesh->stk_bulk(), stk::topology::NODE_RANK, the_new_mesh->stk_meta().universal_part(), the_new_mesh->aux_meta().active_part()); // we should be able to skip this step if there are no higher order elements
  the_new_mesh->update_element_side_parts(); // rebuild conformal side parts

  delete_extraneous_inactive_sides(mesh, the_new_mesh->myRefinementSupport, the_new_mesh->get_parent_part(), the_new_mesh->get_active_part());

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));
}

static bool is_child_elem(const stk::mesh::BulkData & mesh, const stk::mesh::Part & childEdgeNodePart, stk::mesh::Entity elem)
{
  for (auto elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
    if(mesh.bucket(elemNode).member(childEdgeNodePart))
      return true;
  return false;
}

static void batch_change_entity_parts(stk::mesh::BulkData & mesh, 
    const stk::mesh::EntityVector & entitiesWithWrongParts,
    const stk::mesh::PartVector & addParts,
    const stk::mesh::PartVector & removeParts)
{
  if (stk::is_true_on_any_proc(mesh.parallel(), !entitiesWithWrongParts.empty()))
  {
    mesh.batch_change_entity_parts(entitiesWithWrongParts, addParts, removeParts);
  }
}

void
CDMesh::rebuild_child_part()
{
  auto & child_part = get_child_part();
  auto & mesh = stk_bulk();

  // Need to iterate all locally owned elements to find child elements,
  // which are identified by detecting that they use child edge nodes

  const auto & childEdgeNodePart = get_child_edge_node_part();
  const stk::mesh::Selector activeLocallyOwnedNotChild = get_active_part() & get_locally_owned_part() & !child_part;
  
  stk::mesh::EntityVector entitiesWithWrongParts;
  
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwnedNotChild))
    for(const auto & elem : *bucketPtr)
      if(is_child_elem(mesh, childEdgeNodePart, elem))
        entitiesWithWrongParts.push_back(elem);

  batch_change_entity_parts(mesh, entitiesWithWrongParts, stk::mesh::PartVector{&child_part}, stk::mesh::PartVector{});
}

void
CDMesh::rebuild_parent_and_active_parts_using_nonconformal_and_child_parts()
{
  auto & parent_part = get_parent_part();
  auto & mesh = stk_bulk();

  stk::mesh::EntityVector entitiesWithWrongParts;

  // Find parents of child elements are active or do not have parent part
  stk::mesh::Selector locallyOwnedChild = get_locally_owned_part() & get_child_part();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, locallyOwnedChild))
  {
    for(const auto & elem : *bucketPtr)
    {
      auto parentElem = get_parent_element(elem);
      const stk::mesh::Bucket & parentElemBucket = mesh.bucket(parentElem);
      if (parentElemBucket.member(get_active_part()) || !parentElemBucket.member(parent_part))
        entitiesWithWrongParts.push_back(parentElem);
    }
  }
  stk::util::sort_and_unique(entitiesWithWrongParts);
  
  batch_change_entity_parts(mesh, entitiesWithWrongParts, stk::mesh::PartVector{&parent_part}, stk::mesh::PartVector{&get_active_part()});

  // Also remove active part from nonconformal sides
  entitiesWithWrongParts.clear();
  auto sideRank = mesh.mesh_meta_data().side_rank();
  const stk::mesh::Selector activeLocallyOwnedNonConformal = get_active_part() & get_locally_owned_part() &
      stk::mesh::selectUnion(my_phase_support.get_nonconformal_parts_of_rank(sideRank));
  for(const auto & bucketPtr : mesh.get_buckets(sideRank, activeLocallyOwnedNonConformal))
    for (auto && side : *bucketPtr)
      entitiesWithWrongParts.push_back(side);
      
  batch_change_entity_parts(mesh, entitiesWithWrongParts, stk::mesh::PartVector{}, stk::mesh::PartVector{&get_active_part()});
}

const SubElementNode *
CDMesh::find_or_build_subelement_node_with_id(const stk::mesh::EntityId nodeId, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode )
{
  const auto & iter = idToSubElementNode.find(nodeId);
  if (iter != idToSubElementNode.end())
    return iter->second;
  return build_subelement_child_node(stk_bulk().get_entity(stk::topology::NODE_RANK, nodeId), ownerMeshElem, idToSubElementNode);
}

const SubElementNode *
CDMesh::find_or_build_subelement_node(const stk::mesh::Entity node, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode )
{
  const auto & iter = idToSubElementNode.find(stk_bulk().identifier(node));
  if (iter != idToSubElementNode.end())
    return iter->second;
  return build_subelement_child_node(node, ownerMeshElem, idToSubElementNode);
}

void
CDMesh::find_or_build_midside_nodes(const stk::topology & elemTopo, const Mesh_Element & ownerMeshElem, const stk::mesh::Entity * elemNodes, const NodeVec & subelemNodes )
{
  if (elemTopo.num_nodes() > elemTopo.base().num_nodes())
  {
    for (unsigned iEdge=0; iEdge<elemTopo.num_edges(); ++iEdge)
    {
      const unsigned * edgeLNN = get_edge_node_ordinals(elemTopo, iEdge);
      create_midside_node(&ownerMeshElem, subelemNodes[edgeLNN[0]], subelemNodes[edgeLNN[1]], elemNodes[edgeLNN[2]]);
    }
  }
}

const SubElementNode *
CDMesh::build_subelement_child_node(const stk::mesh::Entity node, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode )
{
  const auto & mesh = stk_bulk();

  const auto parentIdsAndWts = get_child_node_parent_ids_and_weights(mesh, get_parent_node_ids_field(), get_parent_node_weights_field(), node);

  std::vector<const SubElementNode *> immediateParents;
  std::vector<double> immediateParentWts;
  immediateParents.reserve(parentIdsAndWts.size());
  immediateParentWts.reserve(parentIdsAndWts.size());
  for (const auto & [parentId, parentWt] : parentIdsAndWts)
  {
    immediateParents.push_back(find_or_build_subelement_node_with_id(parentId, ownerMeshElem, idToSubElementNode));
    immediateParentWts.push_back(parentWt);
  }

  const SubElementNode * childNode = nullptr;
  if (immediateParents.size() == 2)
    childNode = create_edge_node(&ownerMeshElem, immediateParents[0], immediateParents[1], immediateParentWts[1]);
  else
    childNode = create_child_internal_or_face_node(&ownerMeshElem, immediateParents, immediateParentWts);

  childNode->set_entity(stk_bulk(), node);

  idToSubElementNode[mesh.identifier(node)] = childNode;

  return childNode;
}

void
CDMesh::restore_subelements()
{
  stk::mesh::Selector selector = get_locally_owned_part() & get_child_part();
  const auto & mesh = stk_bulk();

  std::map<stk::mesh::EntityId, const SubElementNode*> idToSubElementNode;
  for (auto && node : nodes)
    idToSubElementNode[node->entityId()] = node.get();

  bool error = false;

  const auto & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, selector);
  NodeVec subelem_nodes;
  for(const auto & b_ptr : buckets)
  {
    const stk::topology & topo = b_ptr->topology();
    const unsigned num_base_nodes = topo.base().num_nodes();
    subelem_nodes.reserve(num_base_nodes);
    for(const auto & elem : *b_ptr)
    {
      const stk::mesh::Entity parent = get_parent_element(elem);
      STK_ThrowRequire(mesh.is_valid(parent) && parent != elem);

      auto parentMeshElem = find_mesh_element(mesh.identifier(parent));
      if (!parentMeshElem)
      {
        krinolog << "Could not find Mesh_Element for CDFEM parent element " << mesh.identifier(parent) << " for CDFEM child element " << mesh.identifier(elem) << stk::diag::dendl;
        krinolog << debug_entity_1line(mesh,parent) << stk::diag::dendl;
        krinolog << debug_entity_1line(mesh,elem) << stk::diag::dendl;
        error = true;
        continue;
      }

      subelem_nodes.clear();
      const auto * elem_nodes = mesh.begin_nodes(elem);
      for(unsigned i=0; i < num_base_nodes; ++i)
      {
        const SubElementNode * node = find_or_build_subelement_node(elem_nodes[i], *parentMeshElem, idToSubElementNode);
        subelem_nodes.push_back(node);
      }
      find_or_build_midside_nodes(topo, *parentMeshElem, elem_nodes, subelem_nodes);

      std::unique_ptr<SubElement> subelem;
      switch(topo)
      {
      case stk::topology::TRIANGLE_3_2D:
      case stk::topology::TRIANGLE_6_2D:
        subelem = std::make_unique<SubElement_Tri_3>(subelem_nodes, std::vector<int>{-1, -1, -1}, parentMeshElem);
        break;
      case stk::topology::TETRAHEDRON_4:
      case stk::topology::TETRAHEDRON_10:
        subelem = std::make_unique<SubElement_Tet_4>(subelem_nodes, std::vector<int>{-1, -1, -1, -1}, parentMeshElem);
        break;
      default:
        ThrowRuntimeError("At present only Tri3, Tri6, Tet4 and Tet10 topologies are supported for restart of CDFEM problems.");
      }

      if (topo == stk::topology::TRIANGLE_6_2D || topo == stk::topology::TETRAHEDRON_10)
      {
        subelem->build_quadratic_subelements(*this);
        std::vector<SubElement *> highOrderSubElems;
        subelem->get_subelements( highOrderSubElems );
        STK_ThrowRequire(highOrderSubElems.size() == 1);
        highOrderSubElems[0]->set_entity(stk_bulk(), elem);
      }
      else
      {
        subelem->set_entity(stk_bulk(), elem);
      }

      const_cast<Mesh_Element *>(parentMeshElem)->add_subelement(std::move(subelem));
    }
  }

  ParallelThrowRequire(stk_bulk().parallel(), !error);

  std::vector<SubElement *> subelems;
  for (auto && element : elements)
  {
    element->get_subelements( subelems );
    if (subelems.size() > 1)
    {
      element->set_have_interface();
    }
  }
}

void
CDMesh::fixup_adapted_element_parts(stk::mesh::BulkData & mesh)
{/* %TRACE[SPEC]% */ Tracespec trace__("Mesh::fixup_adapted_element_parts(CDFEM_Support & cdfem_support)"); /* %TRACE% */
  // Fixup volume parts that currently can be messed up by adaptivity.
  // There are two types of fixes:
  // 1. CDFEM parent elements that are activated by Encore adaptivity (this is expected, but needs to be fixed).
  // 2. Conformal elements that have somehow picked up the non-conformal part (this probably shouldn't happen).

  Phase_Support & phase_support = Phase_Support::get(mesh.mesh_meta_data());
  CDFEM_Support & cdfem_support = CDFEM_Support::get(mesh.mesh_meta_data());
  AuxMetaData & aux_meta = AuxMetaData::get(mesh.mesh_meta_data());
  stk::mesh::Selector cdfem_parent_selector = cdfem_support.get_parent_part();

  stk::mesh::Selector locally_owned_selector(mesh.mesh_meta_data().locally_owned_part());
  std::vector<stk::mesh::Entity> entities;

  std::vector<stk::mesh::PartVector> remove_parts;
  stk::mesh::PartVector bucket_remove_parts;
  const stk::mesh::BucketVector & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, locally_owned_selector);
  for (auto&& bucket_ptr : buckets)
  {
    unsigned num_volume_parts = 0;
    stk::mesh::Part * extraneous_nonconformal_part = nullptr;
    const stk::mesh::PartVector & bucket_parts = bucket_ptr->supersets();
    bucket_remove_parts.clear();
    for(auto&& part : bucket_parts)
    {
      if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && part->subsets().empty() && part->topology() != stk::topology::INVALID_TOPOLOGY)
      {
        ++num_volume_parts;
        if (phase_support.is_nonconformal(part))
        {
          extraneous_nonconformal_part = part;
        }
      }
    }
    if (num_volume_parts > 1)
    {
      bucket_remove_parts.push_back(extraneous_nonconformal_part);
    }
    if (cdfem_parent_selector(*bucket_ptr))
    {
      bucket_remove_parts.push_back(&aux_meta.active_part());
    }
    if (!bucket_remove_parts.empty())
    {
      entities.insert(entities.end(), bucket_ptr->begin(), bucket_ptr->end());
      remove_parts.insert(remove_parts.end(), bucket_ptr->size(), bucket_remove_parts);
    }
  }
  stk::mesh::PartVector empty;
  std::vector<stk::mesh::PartVector> add_parts(entities.size(), empty);

  mesh.batch_change_entity_parts(entities, add_parts, remove_parts);
}

//--------------------------------------------------------------------------------

void
CDMesh::stash_field_data(const int stepCount) const
{
  stk::diag::TimeBlock timer__(my_timer_stash_field_data);
  my_stash_step_count = stepCount;
  clear_prolongation_data();

  myProlongPartAndFieldCollections.build(stk_bulk());

  stash_nodal_field_data();
  stash_elemental_field_data();
}

//--------------------------------------------------------------------------------

static void fill_nodes_of_elements_with_subelements_or_changed_phase(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<std::unique_ptr<Mesh_Element>> & meshElements,
    std::set<stk::mesh::Entity> & nodesOfElements)
{
  nodesOfElements.clear();

  for (auto && element : meshElements)
  {
    bool haveSubelementsOrChangedPhase = false;
    if (element->have_subelements())
    {
      haveSubelementsOrChangedPhase = true;
    }
    else
    {
      PhaseTag currentPhase = determine_phase_for_entity(mesh, element->entity(), phaseSupport);
      if (element->get_phase() != currentPhase)
        haveSubelementsOrChangedPhase = true;
    }

    if (haveSubelementsOrChangedPhase)
      for (auto && node : element->get_nodes())
        nodesOfElements.insert(node->entity());
  }
}

static std::set<stk::mesh::Entity> get_nodes_of_elements_with_subelements_or_have_changed_phase(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<std::unique_ptr<Mesh_Element>> & meshElements)
{
  std::set<stk::mesh::Entity> nodesOfElements;
  fill_nodes_of_elements_with_subelements_or_changed_phase(mesh, phaseSupport, meshElements, nodesOfElements);

  communicate_shared_nodes_to_sharing_procs(mesh, nodesOfElements);

  return nodesOfElements;
}

//--------------------------------------------------------------------------------

void
CDMesh::stash_nodal_field_data() const
{

  ProlongationPointData::set_coords_fields(my_spatial_dim, get_coords_field(), get_cdfem_support().get_cdfem_snap_displacements_field());

  // stash child nodes
  {
    stk::mesh::Selector selector = get_locally_owned_part() & get_child_part();

    stk::mesh::BucketVector const& buckets = stk_bulk().get_buckets( stk::topology::ELEMENT_RANK, selector );
    for ( auto&& bucket_ptr : buckets )
    {
      const stk::mesh::Bucket & b = *bucket_ptr;
      const size_t length = b.size();
      for (size_t ielem = 0; ielem < length; ++ielem)
      {
        stk::mesh::Entity elem = b[ielem];
        const unsigned num_elem_nodes = stk_bulk().num_nodes(elem);
        const stk::mesh::Entity* elem_nodes = stk_bulk().begin_nodes(elem);

        for (unsigned inode = 0; inode < num_elem_nodes; ++inode)
        {
          stk::mesh::Entity node = elem_nodes[inode];
          STK_ThrowAssert((stk_bulk().bucket(node).member(get_active_part())));
          ProlongationNodeData *& node_data = my_prolong_node_map[stk_bulk().identifier(node)];
          if (nullptr == node_data)
          {
            const bool communicate_me_to_all_sharers = stk_bulk().bucket(node).member(get_globally_shared_part());
            node_data = new ProlongationNodeData(*this, myProlongPartAndFieldCollections, node, communicate_me_to_all_sharers);
          }
        }
      }
    }
  }

  // Stash all nodes of elements that have child elements or have changed phase.
  // Due to hanging nodes, etc, this is more than just the cut elements.
  for (auto&& node : get_nodes_of_elements_with_subelements_or_have_changed_phase(stk_bulk(), my_phase_support, elements))
  {
    if (stk_bulk().bucket(node).member(get_active_part())) // Don't stash inactive midside nodes
    {
      STK_ThrowAssert(stk_bulk().is_valid(node));
      ProlongationNodeData *& node_data = my_prolong_node_map[stk_bulk().identifier(node)];
      if (nullptr == node_data)
      {
        const bool communicate_me_to_all_sharers = stk_bulk().bucket(node).member(get_globally_shared_part());
        node_data = new ProlongationNodeData(*this, myProlongPartAndFieldCollections, node, communicate_me_to_all_sharers);
      }
    }
  }

  // stash all inter-block nodes
  if (need_nodes_for_prolongation())
  {
    stk::mesh::Selector active_not_ghost_selector(get_active_part() & (get_locally_owned_part() | get_globally_shared_part()));
    const stk::mesh::BucketVector & buckets = stk_bulk().get_buckets(stk::topology::NODE_RANK, active_not_ghost_selector);
    for (auto&& bucket_ptr : buckets)
    {
      unsigned num_conformal_parts = 0;
      const stk::mesh::PartVector & node_parts = bucket_ptr->supersets();
      for(auto&& node_part_ptr : node_parts)
      {
        // This is designed to catch side with block_2 + block_1_air, block_1_air + block_1_solid, etc.
        // These are included so that we can prolongate a node on the block_1_air + block_1_solid + block_2
        // from a node on that same part ownership.  (This is needed in cases where block_2 has other vars).
        if (node_part_ptr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
            !my_phase_support.is_nonconformal(node_part_ptr) &&
            stk::io::is_part_io_part(*node_part_ptr))
        {
          ++num_conformal_parts;
        }
      }

      if (num_conformal_parts > 1)
      {
        for (auto&& node : *bucket_ptr)
        {
          ProlongationNodeData *& node_data = my_prolong_node_map[stk_bulk().identifier(node)];
          if (nullptr == node_data)
          {
            node_data = new ProlongationNodeData(*this, myProlongPartAndFieldCollections, node, false);
          }
        }
      }
    }
  }

  // build facets if needed
  if (need_facets_for_prolongation())
  {
    stk::mesh::Selector active_locally_owned_selector(get_active_part() & get_locally_owned_part());

    const stk::mesh::BucketVector & buckets = stk_bulk().get_buckets(stk_bulk().mesh_meta_data().side_rank(), active_locally_owned_selector);
    for (auto&& bucket_ptr : buckets)
    {
      unsigned num_conformal_parts = 0;
      const stk::mesh::PartVector & side_parts = bucket_ptr->supersets();
      for(auto&& side_part_ptr : side_parts)
      {
        // This is designed to catch sides like block_1_air + block_1_solid, etc and not block_2 + block_1_air.
        // If we include the nondecomposed blocks like block_2, this could result in prolongation of a node
        // on the interface (block_1_air + block_1_solid) from a node on the boundary of the undecomposed block
        // (block_1_air + block_2).
        if (side_part_ptr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
            my_phase_support.is_conformal(side_part_ptr) &&
            stk::io::is_part_io_part(*side_part_ptr))
        {
          ++num_conformal_parts;
        }
      }

      if (num_conformal_parts > 1)
      {
        for (auto&& side : *bucket_ptr)
        {
          STK_ThrowAssert( stk_bulk().num_elements(side) > 0 );

          ProlongationFacet * prolong_facet = new ProlongationFacet(*this, side);
          my_prolong_facets.push_back(prolong_facet);

          if (spatial_dim() == 3)
            prolong_facet->build_and_append_edge_facets(my_prolong_facets);
        }
      }
    }
  }
}

std::vector<std::vector<stk::mesh::Entity>> CDMesh::get_subelements_for_CDFEM_parents(const std::vector<stk::mesh::Entity> & sortedCdfemParentElems) const
{
  std::vector<std::vector<stk::mesh::Entity>> childrenForParents;
  childrenForParents.resize(sortedCdfemParentElems.size());

  stk::mesh::Selector selector = get_locally_owned_part() & get_child_part();

  std::ostringstream errLog;
  const auto & buckets = stk_bulk().get_buckets(stk::topology::ELEMENT_RANK, selector);
  for(const auto & bucketPtr : buckets)
  {
    for(const auto & elem : *bucketPtr)
    {
      const stk::mesh::Entity parent = get_parent_element(elem);
      if (stk_bulk().is_valid(parent))
      {
        auto iter = std::lower_bound(sortedCdfemParentElems.begin(), sortedCdfemParentElems.end(), parent, stk::mesh::EntityLess(stk_bulk()));
        if (iter != sortedCdfemParentElems.end() && *iter == parent)
        {
          const size_t index = std::distance(sortedCdfemParentElems.begin(), iter);
          childrenForParents[index].push_back(elem);
        }
        else
        {
          errLog << "Failed to find parent element:\n " << debug_entity_1line(stk_bulk(), parent) << " for child element " << debug_entity_1line(stk_bulk(), elem);
        }
      }
      else
      {
        errLog << "Parent element:\n " << debug_entity_1line(stk_bulk(), parent) << "is not valid for child element " << debug_entity_1line(stk_bulk(), elem);
      }
    }
  }

  RequireEmptyErrorMsg(stk_bulk().parallel(), errLog.str(), "Error in get_subelements_for_CDFEM_parents:");
  return childrenForParents;
}

//--------------------------------------------------------------------------------

void
CDMesh::stash_elemental_field_data() const
{
  const bool haveElemFields = !get_element_fields().empty();

  const std::vector<stk::mesh::Entity> nonconformalElems = get_nonconformal_elements();
  const auto subelemsForParents = get_subelements_for_CDFEM_parents(nonconformalElems);

  for (size_t iParent = 0; iParent < nonconformalElems.size(); ++iParent)
  {
    const stk::mesh::Entity parent = nonconformalElems[iParent];
    const stk::mesh::EntityId parentId = stk_bulk().identifier(parent);
    const std::vector<stk::mesh::Entity> & subelems = subelemsForParents[iParent];

    if (subelems.empty())
    {
      ProlongationElementData * elem_data = new ProlongationLeafElementData(*this, myProlongPartAndFieldCollections, parent);
      STK_ThrowAssert(0 == my_prolong_element_map.count(parentId));
      my_prolong_element_map[parentId] = elem_data;
    }
    else
    {
      const unsigned numSubelems = subelems.size();
      std::vector<const ProlongationElementData *> subelemsData(numSubelems);

      for (unsigned iSub=0; iSub<numSubelems; ++iSub)
      {
        const stk::mesh::EntityId subelemId = stk_bulk().identifier(subelems[iSub]);
        ProlongationElementData * subElemData = new ProlongationLeafElementData(*this, myProlongPartAndFieldCollections, subelems[iSub]);
        STK_ThrowAssertMsg(0 == my_prolong_element_map.count(subelemId), "Duplicate subelement entityId " << subelemId);
        my_prolong_element_map[subelemId] = subElemData;
        subelemsData[iSub] = subElemData;
      }

      const bool single_coincident_subelement = (numSubelems == 1);
      if (!single_coincident_subelement)
      {
        ProlongationElementData * elemData = new ProlongationParentElementData(*this, parent, subelemsData, haveElemFields);
        STK_ThrowAssert(0 == my_prolong_element_map.count(parentId));
        my_prolong_element_map[parentId] = elemData;
      }
    }
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::clear_prolongation_trees() const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::clear_prolongation_trees() const"); /* %TRACE% */
  my_phase_prolong_tree_map.clear();
}

//--------------------------------------------------------------------------------

void
CDMesh::build_prolongation_trees() const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::build_prolongation_trees() const"); /* %TRACE% */

  clear_prolongation_trees();

  if (need_facets_for_prolongation())
  {
    std::map<std::vector<unsigned>,std::vector<const ProlongationFacet*>> phase_prolong_facet_map;

    for ( unsigned n=0; n<my_prolong_facets.size(); ++n )
    {
      const ProlongationFacet * prolong_facet = my_prolong_facets[n];
      phase_prolong_facet_map[prolong_facet->compute_common_fields()].push_back(prolong_facet);
    }

    for (auto && entry : phase_prolong_facet_map)
    {
      const std::vector<unsigned> & fields = entry.first;
      std::vector<const ProlongationFacet *> & facets = entry.second;

      my_phase_prolong_tree_map[fields] = std::make_unique<SearchTree<const ProlongationFacet*>>(facets, ProlongationFacet::get_centroid, ProlongationFacet::insert_into_bounding_box);
      STK_ThrowRequire(!my_phase_prolong_tree_map[fields]->empty());
    }
  }
}

//--------------------------------------------------------------------------------

std::string print_fields(const stk::mesh::MetaData & meta, const std::vector<unsigned> & fieldOrdinals)
{
  const stk::mesh::FieldVector & all_fields = meta.get_fields();
  std::ostringstream os;
  os << "Fields { ";
  for (unsigned fieldOrdinal : fieldOrdinals)
    os << all_fields[fieldOrdinal]->name() << " ";
  os << "}";
  return os.str();
}

//--------------------------------------------------------------------------------

static
void pack_facet_fields_for_all_other_procs(const PhaseProlongTreeMap & phaseProlongTreeMap,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      buffer.pack(phaseProlongTreeMap.size());

      for (auto && entry : phaseProlongTreeMap)
      {
        const std::vector<unsigned> & facetFields = entry.first;
        buffer.pack(facetFields.size());
        for (unsigned field : facetFields)
          buffer.pack(field);
      }
    }
  });
}

static
void receive_facet_fields_from_other_procs(PhaseProlongTreeMap & phaseProlongTreeMap,
    stk::CommSparse &commSparse)
{
  std::vector<unsigned> facetFields;

  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    size_t numEntries = 0;
    buffer.unpack(numEntries);

    for (size_t ientry = 0; ientry<numEntries; ++ientry)
    {
      size_t numFields = 0;
      buffer.unpack(numFields);

      facetFields.resize(numFields);
      buffer.unpack(facetFields.data(), numFields);

      if (phaseProlongTreeMap.find(facetFields) == phaseProlongTreeMap.end())
      {
        phaseProlongTreeMap[facetFields].reset();
      }
    }
  });
}

void
CDMesh::communicate_prolongation_facet_fields() const
{
  if (1 ==  stk_bulk().parallel_size() || !need_facets_for_prolongation())
    return;

  stk::CommSparse commSparse(stk_bulk().parallel());
  pack_facet_fields_for_all_other_procs(my_phase_prolong_tree_map, commSparse);
  receive_facet_fields_from_other_procs(my_phase_prolong_tree_map, commSparse);
}

//--------------------------------------------------------------------------------

static void find_nearest_matching_prolong_facet(const stk::mesh::BulkData & mesh,
    const PhaseProlongTreeMap & phaseProlongTreeMap,
    const std::vector<unsigned> & requiredFields,
    const SubElementNode & targetNode,
    const ProlongationFacet *& nearestProlongFacet,
    FacetDistanceQuery<Facet> & nearestFacetQuery,
    bool & haveMissingRemoteProlongFacets)
{
  const stk::math::Vector3d & targetCoordinates = targetNode.coordinates();
  bool haveMatchingButEmptyTree = false;
  for (auto && entry : phaseProlongTreeMap)
  {
    const std::vector<unsigned> & treeFields = entry.first;
    SearchTree<const ProlongationFacet*> * facetTree = entry.second.get();

    if (std::includes(treeFields.begin(), treeFields.end(), requiredFields.begin(), requiredFields.end()))
    {
      if (nullptr == facetTree)
      {
        haveMatchingButEmptyTree = true;
      }
      else
      {
        std::vector<const ProlongationFacet*> nearest_prolong_facets;
        facetTree->find_closest_entities( targetCoordinates, nearest_prolong_facets );
        STK_ThrowAssert(!nearest_prolong_facets.empty());

        for (auto && prolong_facet : nearest_prolong_facets)
        {
          FacetDistanceQuery<Facet> facet_query(*prolong_facet->get_facet(), targetCoordinates);
          if (nearestFacetQuery.empty() || facet_query.distance_squared() < nearestFacetQuery.distance_squared())
          {
            nearestProlongFacet = prolong_facet;
            nearestFacetQuery = facet_query;
          }
        }
      }
    }
  }

  if (nullptr == nearestProlongFacet && haveMatchingButEmptyTree)
  {
    haveMissingRemoteProlongFacets = true;
    if (krinolog.shouldPrint(LOG_DEBUG))
    {
      krinolog << "Found missing remote prolong facet for node for " << targetNode.entityId() << stk::diag::dendl;
      krinolog << "  with required part fields = " << print_fields(mesh.mesh_meta_data(), requiredFields) << stk::diag::dendl;
    }
  }
}

static void write_diagnostics_for_node_not_found(const stk::mesh::BulkData & mesh,
    const PhaseProlongTreeMap & phaseProlongTreeMap,
    const SubElementNode & targetNode,
    const std::vector<unsigned> & requiredFields)
{
  krinolog << "Failed to find prolongation node for node#" << targetNode.entityId() << stk::diag::dendl;
  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    krinolog << "  with  required part fields=" << print_fields(mesh.mesh_meta_data(), requiredFields) << stk::diag::dendl;
    krinolog << "  with parts=";
    const stk::mesh::PartVector & parts = mesh.bucket(targetNode.entity()).supersets();
    for(stk::mesh::PartVector::const_iterator part_iter = parts.begin(); part_iter != parts.end(); ++part_iter)
    {
      const stk::mesh::Part * const part = *part_iter;
      krinolog << "\"" << part->name() << "\"" << " ";
    }
    krinolog << stk::diag::dendl;
    const unsigned num_dst_node_elements = mesh.num_elements(targetNode.entity());
    const stk::mesh::Entity* dst_node_elements = mesh.begin_elements(targetNode.entity());
    for (unsigned dst_node_elem_index=0; dst_node_elem_index<num_dst_node_elements; ++dst_node_elem_index)
    {
      stk::mesh::Entity elem = dst_node_elements[dst_node_elem_index];

      krinolog << "  Elem: id=" << mesh.identifier(elem) << stk::diag::dendl;
      krinolog << "    Mesh parts=";
      const stk::mesh::PartVector & elem_parts = mesh.bucket(targetNode.entity()).supersets();
      for(stk::mesh::PartVector::const_iterator part_iter = elem_parts.begin(); part_iter != elem_parts.end(); ++part_iter)
      {
        const stk::mesh::Part * const part = *part_iter;
        krinolog << "\"" << part->name() << "\"" << " ";
      }
      krinolog << stk::diag::dendl;
    }

    krinolog << "Candidate prolongation facets:" << stk::diag::dendl;
    for (auto && entry : phaseProlongTreeMap)
    {
      const std::vector<unsigned> & tree_fields = entry.first;
      krinolog << "  matching fields=" << std::includes(tree_fields.begin(), tree_fields.end(), requiredFields.begin(), requiredFields.end())
               << ", tree fields=" << print_fields(mesh.mesh_meta_data(), tree_fields)
               << stk::diag::dendl;
    }
  }
}

ProlongationQuery
CDMesh::find_prolongation_node(const SubElementNode & targetNode) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::find_prolongation_node(const SubElementNode & dst_node) const"); /* %TRACE% */

  const std::vector<unsigned> requiredFields = targetNode.prolongation_node_fields(*this);

  const ProlongationFacet * nearestProlongFacet = nullptr;
  FacetDistanceQuery<Facet> nearestFacetQuery;

  find_nearest_matching_prolong_facet(stk_bulk(), my_phase_prolong_tree_map, requiredFields, targetNode, nearestProlongFacet, nearestFacetQuery, my_missing_remote_prolong_facets);

  if(nullptr != nearestProlongFacet)
  {
    ProlongationQuery prolongationQuery(*nearestProlongFacet, nearestFacetQuery);
    if (krinolog.shouldPrint(LOG_DEBUG))
    {
      const std::vector<const ProlongationNodeData *> & facetNodes = nearestProlongFacet->get_prolongation_nodes();
      krinolog << "Prolongation facet for " << targetNode.entityId() << " with ancestry " << targetNode.get_ancestry() << " has nodes ";
      for (auto&& node : facetNodes)
      {
        krinolog << node->entityId() << " " << node->get_previous_coordinates() << "  ";
      }
      krinolog << *(nearestProlongFacet->get_facet()) << stk::diag::dendl;
      krinolog << "  with required fields " << print_fields(stk_meta(), requiredFields) << stk::diag::dendl;
      krinolog << "Prolongation data for node#" << stk_bulk().identifier(targetNode.entity()) << " (" << targetNode.coordinates() << ")"
               << " will be point at location (" << prolongationQuery.get_prolongation_point_data()->get_previous_coordinates() << ")" << stk::diag::dendl;
    }
    return prolongationQuery;
  }

  // Search for facet failed.  Now try nodes.  This will handle triple points.  Something better that handles an actual edge search might be better in 3d.
  const stk::math::Vector3d & targetCoordinates = targetNode.coordinates();
  const ProlongationNodeData * closest_node = nullptr;
  double closest_dist2 = std::numeric_limits<double>::max();
  for (auto && entry : my_prolong_node_map)
  {
    const ProlongationNodeData * node = entry.second;
    const std::vector<unsigned> & tree_fields = myProlongPartAndFieldCollections.get_fields(node->get_field_collection_id());
    if (std::includes(tree_fields.begin(), tree_fields.end(), requiredFields.begin(), requiredFields.end()))
    {
      const double dist2 = (node->get_previous_coordinates()-targetCoordinates).length_squared();
      if (dist2 < closest_dist2)
      {
        closest_node = node;
        closest_dist2 = dist2;
      }
    }
  }

  ProlongationQuery prolongationQuery(closest_node);

  if (nullptr == closest_node)
    write_diagnostics_for_node_not_found(stk_bulk(),  my_phase_prolong_tree_map, targetNode, requiredFields);
  else if (krinolog.shouldPrint(LOG_DEBUG))
  {
    krinolog << "Prolongation node for " << targetNode.entityId() << " is " << closest_node->entityId() << stk::diag::dendl;
    krinolog << "Prolongation data for node#" << stk_bulk().identifier(targetNode.entity()) << " (" << targetNode.coordinates() << ")"
               << " will be point at location (" << prolongationQuery.get_prolongation_point_data()->get_previous_coordinates() << ")" << stk::diag::dendl;

    krinolog << "With required fields " << print_fields(stk_meta(), requiredFields) << stk::diag::dendl;
    krinolog << "Prolongation data for node#" << stk_bulk().identifier(targetNode.entity()) << " (" << targetNode.coordinates() << ")"
             << " will be point at location (" << prolongationQuery.get_prolongation_point_data()->get_previous_coordinates() << ") with distance " << std::sqrt(closest_dist2) << stk::diag::dendl;
  }

  return prolongationQuery;
}

//--------------------------------------------------------------------------------

static bool entity_has_any_node_in_selector(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const stk::mesh::Selector & selector)
{
  const unsigned numNodes = mesh.num_nodes(entity);
  stk::mesh::Entity const* entityNodes = mesh.begin_nodes(entity);
  for (unsigned n=0; n<numNodes; ++n)
    if (selector(mesh.bucket(entityNodes[n]))) return true;
  return false;
}

//--------------------------------------------------------------------------------

static std::vector<stk::mesh::Entity> get_elements_that_might_get_decomposed(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & ownedPossibleParentElementSelector,
    const stk::mesh::Selector & allDecomposedBlocksSelector)
{
  std::vector<stk::mesh::Entity> elems;
  for (auto&& bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedPossibleParentElementSelector))
  {
    if (Mesh_Element::is_supported_topology(bucket->topology()))
    {
      if (allDecomposedBlocksSelector(bucket))
      {
        elems.insert(elems.end(), bucket->begin(), bucket->end());
      }
      else
      {
        for (auto&& elem : *bucket)
          if (entity_has_any_node_in_selector(mesh, elem, allDecomposedBlocksSelector))
            elems.push_back(elem);
      }
    }
  }

  return elems;
}

std::vector<stk::mesh::Entity>
CDMesh::get_nonconformal_elements() const
{
  stk::mesh::Selector possibleParentElementSelector = get_locally_owned_part() & (get_parent_part() | (get_active_part() & !get_child_part()));

  std::vector<stk::mesh::Entity> elems = get_elements_that_might_get_decomposed(stk_bulk(), possibleParentElementSelector, my_phase_support.get_all_decomposed_blocks_selector());
  std::sort(elems.begin(), elems.end(), stk::mesh::EntityLess(stk_bulk()));

  return elems;
}

//--------------------------------------------------------------------------------

void
CDMesh::generate_nonconformal_elements()
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::generate_nonconformal_elements()"); /* %TRACE% */
  ParallelThrowRequire(stk_bulk().parallel(), nodes.empty());
  ParallelThrowRequire(stk_bulk().parallel(), elements.empty());

  const std::vector<stk::mesh::Entity> nonconformalElems = get_nonconformal_elements();
  elements.reserve(nonconformalElems.size());
  for(auto && elem : nonconformalElems)
  {
    auto mesh_elem = std::make_unique<Mesh_Element>(*this, elem);
    elements.emplace_back(std::move(mesh_elem));
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::clear()
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::clear()"); /* %TRACE% */

  nodes.clear();
  elements.clear();

  clear_prolongation_data();

  mesh_node_map.clear();
  child_elements.clear();
}

//--------------------------------------------------------------------------------

void
CDMesh::clear_prolongation_data() const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::clear()"); /* %TRACE% */
  for (auto && map_entry : my_prolong_node_map)
  {
    delete map_entry.second;
  }
  my_prolong_node_map.clear();
  for (auto && map_entry : my_prolong_element_map)
  {
    delete map_entry.second;
  }
  my_prolong_element_map.clear();

  for (auto && prolong_facet : my_prolong_facets)
  {
    delete prolong_facet;
  }
  my_prolong_facets.clear();

  clear_prolongation_trees();
}

//--------------------------------------------------------------------------------

PhaseTag
CDMesh::determine_entity_phase(stk::mesh::Entity entity) const
{
  return determine_phase_for_entity(stk_bulk(), entity, my_phase_support);
}

//--------------------------------------------------------------------------------

bool
CDMesh::elem_io_part_changed(const ElementObj & elem) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::verify_elem_part(const Mesh_Element * elem) const"); /* %TRACE% */
  const stk::mesh::Part & current_elem_io_part = find_element_part(stk_bulk(),elem.entity());
  const stk::mesh::Part * const conformal_elem_io_part = my_phase_support.find_conformal_io_part(current_elem_io_part, elem.get_phase());
  return (&current_elem_io_part != conformal_elem_io_part || !stk_bulk().bucket(elem.entity()).member(get_active_part()));
}

//--------------------------------------------------------------------------------

void
CDMesh::determine_nonconformal_parts(
    stk::mesh::Entity entity,
    stk::mesh::PartVector & add_parts,
    stk::mesh::PartVector & remove_parts) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::determine_nonconformal_parts(stk::mesh::Entity entity, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const"); /* %TRACE% */

  add_parts.clear();
  remove_parts.clear();

  const auto & all_decomposed_blocks_selector = my_phase_support.get_all_decomposed_blocks_selector();
  stk::mesh::EntityRank entity_rank = stk_bulk().entity_rank(entity);
  const stk::mesh::PartVector & current_parts = stk_bulk().bucket(entity).supersets();
  for(stk::mesh::PartVector::const_iterator part_iter = current_parts.begin(); part_iter != current_parts.end(); ++part_iter)
  {
    stk::mesh::Part & part = **part_iter;
    if( part.primary_entity_rank() == entity_rank && all_decomposed_blocks_selector(part) )
    {
      stk::mesh::Part * nonconformal_io_part = const_cast<stk::mesh::Part *>(my_phase_support.find_nonconformal_part(part));
      if (nullptr != nonconformal_io_part && nonconformal_io_part != &part)
      {
        add_parts.push_back(nonconformal_io_part);
        remove_parts.push_back(&part);

        for(stk::mesh::PartVector::const_iterator sup_it = part.supersets().begin(); sup_it != part.supersets().end(); ++sup_it)
        {
          stk::mesh::Part & superset = **sup_it;
          if (!stk::mesh::is_auto_declared_part(superset))
          {
            remove_parts.push_back(&superset);
          }
        }
      }
    }
  }

  // Set to inactive
  remove_parts.push_back(&aux_meta().active_part());

  if (entity_rank == stk::topology::ELEMENT_RANK)
  {
    add_parts.push_back(&get_parent_part());
    remove_parts.push_back(&get_child_part());
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::determine_conformal_parts(
    const stk::mesh::PartVector & current_parts,
    const stk::mesh::EntityRank entity_rank,
    const PhaseTag & phase,
    stk::mesh::PartVector & add_parts,
    stk::mesh::PartVector & remove_parts) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::determine_conformal_parts(stk::mesh::Entity entity, const PhaseTag & phase, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const"); /* %TRACE% */

  const auto & all_decomposed_blocks_selector = my_phase_support.get_all_decomposed_blocks_selector();
  for(stk::mesh::PartVector::const_iterator part_iter = current_parts.begin(); part_iter != current_parts.end(); ++part_iter)
  {
    stk::mesh::Part & part = **part_iter;
    if( part.primary_entity_rank() == entity_rank &&
        (stk::io::is_part_io_part(part) || all_decomposed_blocks_selector(&part)) )
    {
      stk::mesh::Part * conformal_elem_io_part = const_cast<stk::mesh::Part *>(my_phase_support.find_conformal_io_part(part, phase));
      if (nullptr != conformal_elem_io_part && conformal_elem_io_part != &part)
      {
        add_parts.push_back(conformal_elem_io_part);
        remove_parts.push_back(&part);

        for(stk::mesh::PartVector::const_iterator sup_it = part.supersets().begin(); sup_it != part.supersets().end(); ++sup_it)
        {
          stk::mesh::Part & superset = **sup_it;
          if (!stk::mesh::is_auto_declared_part(superset))
          {
            remove_parts.push_back(&superset);
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::determine_conformal_parts(
    stk::mesh::Entity entity,
    const PhaseTag & phase,
    stk::mesh::PartVector & add_parts,
    stk::mesh::PartVector & remove_parts) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::determine_conformal_parts(stk::mesh::Entity entity, const PhaseTag & phase, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const"); /* %TRACE% */

  add_parts.clear();
  remove_parts.clear();

  STK_ThrowAssert(stk_bulk().is_valid(entity));

  stk::mesh::EntityRank entity_rank = stk_bulk().entity_rank(entity);
  const stk::mesh::PartVector & current_parts = stk_bulk().bucket(entity).supersets();
  determine_conformal_parts(current_parts, entity_rank, phase, add_parts, remove_parts);
}

//--------------------------------------------------------------------------------

void
CDMesh::determine_child_conformal_parts(
    stk::topology topology,
    const stk::mesh::PartVector & parent_parts,
    const PhaseTag & phase,
    stk::mesh::PartVector & child_parts) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::determine_child_conformal_parts(stk::mesh::Entity entity, const PhaseTag & phase, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const"); /* %TRACE% */

  child_parts.clear();

  const auto & all_decomposed_blocks_selector = my_phase_support.get_all_decomposed_blocks_selector();
  stk::mesh::EntityRank entity_rank = topology.rank();
  for(stk::mesh::PartVector::const_iterator part_iter = parent_parts.begin(); part_iter != parent_parts.end(); ++part_iter)
  {
    stk::mesh::Part & part = **part_iter;
    if( part.primary_entity_rank() == entity_rank &&
        (stk::io::is_part_io_part(part) || all_decomposed_blocks_selector(&part)) )
    {
      stk::mesh::Part * conformal_elem_io_part = const_cast<stk::mesh::Part *>(my_phase_support.find_conformal_io_part(part, phase));
      if (nullptr != conformal_elem_io_part && !my_phase_support.is_interface(&part))
      {
        child_parts.push_back(conformal_elem_io_part);
      }
    }
    else if (stk::mesh::contain(my_attribute_parts, part))
    {
        child_parts.push_back(&part);
    }
  }

  child_parts.push_back(&stk_meta().get_topology_root_part(topology));

  if (entity_rank == stk::topology::ELEMENT_RANK)
  {
    child_parts.push_back(&get_child_part());
  }

  // Set to active
  child_parts.push_back(&aux_meta().active_part());
}

//--------------------------------------------------------------------------------

bool
CDMesh::triangulate(const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::triangulate(InterfaceGeometry & interfaceGeometry)"); /* %TRACE% */
  bool made_changes = false;
  for (auto && elem : elements)
  {
    made_changes |= elem->triangulate(*this, interfaceGeometry);
  }
  return made_changes;
}

//--------------------------------------------------------------------------------

void
CDMesh::cut_sharp_features(const InterfaceGeometry & interfaceGeometry)
{
  if (!interfaceGeometry.might_have_interior_or_face_intersections())
    return;

  for (auto && elem : elements)
  {
    elem->cut_interior_intersection_points(*this);
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::set_phase_of_uncut_elements(const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::snap_nearby_intersections_to_nodes(void)"); /* %TRACE% */

  const std::vector<Surface_Identifier> & surfaceIDs = interfaceGeometry.get_surface_identifiers();
  const bool oneLSPerPhase = my_phase_support.has_one_levelset_per_phase();
  for (auto && entry : interfaceGeometry.get_phase_for_uncut_elements())
  {
    Mesh_Element * elem = find_mesh_element(stk_bulk().identifier(entry.first));
    if (elem)
    {
      PhaseTag elemPhase;
      if (oneLSPerPhase)
      {
        elemPhase.add(surfaceIDs[entry.second], -1);
        elem->set_phase(elemPhase);
      }
      else
      {
        STK_ThrowRequire(1 == surfaceIDs.size());
        elemPhase.add(surfaceIDs[0], entry.second);
        elem->set_phase(elemPhase);
      }
      if (false)
        krinolog << "Set phase for elem " << stk_bulk().identifier(entry.first) << ".\n" << elem->visualize(*this) << stk::diag::dendl;
    }
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::snap_nearby_intersections_to_nodes(const InterfaceGeometry & interfaceGeometry, NodeToCapturedDomainsMap & domainsAtNodes)
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::snap_nearby_intersections_to_nodes(void)"); /* %TRACE% */
  stk::diag::TimeBlock timer__(my_timer_snap);

  const stk::mesh::Selector parentElementSelector = get_cdfem_parent_element_selector(get_active_part(), my_cdfem_support, my_phase_support);
  snap_to_node(stk_bulk(), parentElementSelector, interfaceGeometry, get_snapper(), domainsAtNodes);
  for (auto && entry : domainsAtNodes)
  {
    const SubElementNode * node = get_mesh_node(stk_bulk().identifier(entry.first));
    if (node)
      node->set_node_domains(entry.second);
  }

  domainsAtNodes.clear(); // done using this
}

//--------------------------------------------------------------------------------

void
CDMesh::decompose(const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::decompose(void)"); /* %TRACE% */

  if (my_cdfem_support.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
    cut_sharp_features(interfaceGeometry);

  const std::vector<Surface_Identifier> & surfaceIDs = interfaceGeometry.get_surface_identifiers();

  // TODO: N^2 in number of phases
  for (auto && interface : active_interface_ids(surfaceIDs))
  {
    determine_node_signs(interface);
    decompose_edges(interface);
    determine_node_scores(interface);
    handle_hanging_children(interface);
  }
  for (auto && elem : elements)
  {
    elem->build_quadratic_subelements(*this);
  }
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << stk::diag::dendl;

  for (auto && elem : elements )
  {
    if(krinolog.shouldPrint(LOG_DEBUG))
    {
      krinolog << "Determining subelement phases for Mesh_Element local_id=" << " identifier=" << elem->entityId();
      krinolog << "\n";
    }
    elem->determine_decomposed_elem_phase(surfaceIDs);
    if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "\n";
  }
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << stk::diag::dendl;
}

void
CDMesh::build_parallel_hanging_edge_nodes()
{ /* %TRACE[ON]% */ Trace trace__("krino::CDMesh::build_parallel_hanging_edge_nodes(void)"); /* %TRACE% */
  stk::mesh::BulkData & mesh = stk_bulk();

  if (mesh.parallel_size() < 2) return;

  // Get all cut edges in the mesh that are parallel shared. Processing by edge nodes should be cheaper than processing
  // by elements since we don't have to deal with duplicates.
  std::vector<SubElementChildNodeAncestry> shared_edge_nodes;
  for (auto&& node : nodes)
  {
    const SubElementEdgeNode * edge_node = dynamic_cast<const SubElementEdgeNode *>( node.get() );
    if (nullptr != edge_node)
    {
      if (SubElementChildNodeAncestry::is_shared(mesh, edge_node))
      {
        shared_edge_nodes.emplace_back(edge_node);
      }
    }
  }

  std::vector<int> sharing_procs;
  std::vector<stk::mesh::EntityKey> edge_node_keys;

  stk::CommSparse comm_spec(mesh.parallel());

  for (int phase=0;phase<2;++phase)
  {
    for (auto&& shared_edge_node : shared_edge_nodes)
    {
      shared_edge_node.get_parent_node_keys(edge_node_keys);
      stk_bulk().shared_procs_intersection(edge_node_keys, sharing_procs);

      for (auto&& other_proc : sharing_procs)
      {
        if (other_proc != mesh.parallel_rank())
        {
          shared_edge_node.pack_into_buffer(comm_spec.send_buffer(other_proc));
        }
      }
    }

    if ( phase == 0 )
    {
      comm_spec.allocate_buffers();
    }
    else
    {
      comm_spec.communicate();
    }
  }

  for(int i = 0; i < mesh.parallel_size(); ++i)
  {
    if(i != mesh.parallel_rank())
    {
      while(comm_spec.recv_buffer(i).remaining())
      {
        SubElementChildNodeAncestry shared_child_node(comm_spec.recv_buffer(i));
        shared_child_node.build_missing_child_nodes(*this);
      }
    }
  }
}

void
CDMesh::determine_node_signs(const InterfaceID & interface)
{
  for (auto && node : nodes)
  {
    node->clear_node_sign();
  }
  for (auto && elem : elements)
  {
    elem->determine_node_signs(*this, interface);
  }
  sync_node_signs_on_constrained_nodes();
  parallel_sync_node_signs_on_shared_nodes();
}

void
CDMesh::decompose_edges(const InterfaceID & interface)
{
  for (auto && elem : elements)
  {
    if(krinolog.shouldPrint(LOG_DEBUG))
    {
      krinolog << "Decomposing Mesh_Element local_id=" << " identifier=" << elem->entityId();
      krinolog << "\n";
    }
    elem->decompose(*this, interface);
    if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "\n";
  }
}

void
CDMesh::determine_node_scores(const InterfaceID & interface)
{
  for (auto && node : nodes)
  {
    node->clear_node_score();
  }
  for (auto && elem : elements)
  {
    elem->determine_node_scores(*this, interface);
  }
  sync_node_scores_on_constrained_nodes();
  parallel_sync_node_scores_on_shared_nodes();
}

template <typename T>
std::vector<std::vector<int>> determine_owning_procs_of_nodes_in_ancestries(const stk::mesh::BulkData & mesh, const std::vector<std::pair<SubElementChildNodeAncestry,T>> & constrainedNodesAndData)
{
  std::vector<stk::mesh::EntityKey> edgeNodeKeys;

  std::vector<std::vector<int>> owningProcsOfNodesInAncestries;
  owningProcsOfNodesInAncestries.reserve(constrainedNodesAndData.size());
  for (auto&& constrainedNodeAndData : constrainedNodesAndData)
  {
    owningProcsOfNodesInAncestries.emplace_back();
    std::vector<int> & owningProcs = owningProcsOfNodesInAncestries.back();

    const auto & nodeAncestry = constrainedNodeAndData.first;
    nodeAncestry.get_parent_node_keys(edgeNodeKeys);
    for (auto && edgeNodeKey : edgeNodeKeys)
      owningProcs.push_back(mesh.parallel_owner_rank(mesh.get_entity(edgeNodeKey))); //Expensive?

    stk::util::sort_and_unique(owningProcs);
  }

  return owningProcsOfNodesInAncestries;
}

template <typename T>
std::vector<std::vector<int>> determine_sharing_procs_of_nodes_in_ancestries(const stk::mesh::BulkData & mesh, const std::vector<std::pair<SubElementChildNodeAncestry,T>> & sharedNodesAndData)
{
  std::vector<stk::mesh::EntityKey> edgeNodeKeys;

  std::vector<std::vector<int>> sharingProcsOfNodesInAncestries;
  sharingProcsOfNodesInAncestries.reserve(sharedNodesAndData.size());
  for (auto&& sharedNodeAndData : sharedNodesAndData)
  {
    sharingProcsOfNodesInAncestries.emplace_back();
    std::vector<int> & sharingProcs = sharingProcsOfNodesInAncestries.back();

    const auto & nodeAncestry = sharedNodeAndData.first;
    nodeAncestry.get_parent_node_keys(edgeNodeKeys);
    mesh.shared_procs_intersection(edgeNodeKeys, sharingProcs);
  }

  return sharingProcsOfNodesInAncestries;
}

template <typename T>
void pack_node_data_for_node_ancestries(const stk::mesh::BulkData & mesh, const std::vector<std::pair<SubElementChildNodeAncestry,T>> & nodeAncestriesAndData, const std::vector<std::vector<int>> & destinationProcs, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    STK_ThrowAssert(nodeAncestriesAndData.size() == destinationProcs.size());
    std::vector<stk::mesh::EntityKey> edgeNodeKeys;

    for (size_t i=0; i<nodeAncestriesAndData.size(); ++i)
    {
      const auto & nodeAncestry = nodeAncestriesAndData[i].first;
      const T nodeData = nodeAncestriesAndData[i].second;
      nodeAncestry.get_parent_node_keys(edgeNodeKeys);

      for (auto&& otherProc : destinationProcs[i])
      {
        if (otherProc != commSparse.parallel_rank()) // please avoid talking to yourself
        {
          auto & buffer = commSparse.send_buffer(otherProc);
          nodeAncestry.pack_into_buffer(buffer);
          buffer.pack(nodeData);
        }
      }
    }
  });
}

template <typename T>
void set_node_sign_or_score(const SubElementNode * node, const T & signOrScore) { STK_ThrowRequireMsg(false, "Unsupported type in set_node_sign_or_score."); }

template <>
void set_node_sign_or_score(const SubElementNode * node, const int & sign) { node->set_node_sign(sign); }

template <>
void set_node_sign_or_score(const SubElementNode * node, const double & score) { node->set_node_score(score); }

template <typename T>
T get_node_sign_or_score(const SubElementNode * node) { STK_ThrowRequireMsg(false, "Unsupported type in get_node_sign_or_score."); }

template <>
int get_node_sign_or_score(const SubElementNode * node) { return node->get_node_sign(); }

template <>
double get_node_sign_or_score(const SubElementNode * node) { return node->get_node_score(); }

template <typename T>
bool node_sign_or_score_is_set(const SubElementNode * node) { STK_ThrowRequireMsg(false, "Unsupported type in node_sign_or_score_is_set."); return false; }

template <>
bool node_sign_or_score_is_set<int>(const SubElementNode * node) { return node->node_sign_is_set(); }

template <>
bool node_sign_or_score_is_set<double>(const SubElementNode * node) { return node->node_score_is_set(); }

template <typename T>
std::vector<std::pair<SubElementChildNodeAncestry,T>> gather_constrained_node_ancestries_and_sign_or_score(const std::vector<std::unique_ptr<SubElementNode>> & nodes,
    const std::unordered_map<stk::mesh::EntityId, std::vector<stk::mesh::EntityId> > & periodicNodeIDMap)
{
  std::vector<std::pair<SubElementChildNodeAncestry,T>> ancestriesAndNodeOrScore;
  for (auto && node : nodes)
  {
    if (node_sign_or_score_is_set<T>(node.get()))
    {
      SubElementChildNodeAncestry nodeAncestry(node.get());
      const auto & constrainedNodeAncestries = nodeAncestry.get_constrained_node_ancestries(periodicNodeIDMap);
      for (auto && constrainedNodeAncestry : constrainedNodeAncestries)
        ancestriesAndNodeOrScore.emplace_back(constrainedNodeAncestry,get_node_sign_or_score<T>(node.get()));
    }
  }
  return ancestriesAndNodeOrScore;
}

template <typename T>
std::vector<std::pair<SubElementChildNodeAncestry,T>> gather_shared_node_ancestries_and_sign_or_score(const stk::mesh::BulkData & mesh,
    const std::vector<std::unique_ptr<SubElementNode>> & nodes)
{
  // Get all cut edges in the mesh that are parallel shared. Processing by edge nodes should be cheaper than processing
  // by elements since we don't have to deal with duplicates.

  std::vector<std::pair<SubElementChildNodeAncestry,T>> sharedNodeAncestriesAndSignOrScore;
  for (auto && node : nodes)
    if (node_sign_or_score_is_set<T>(node.get()) && SubElementChildNodeAncestry::is_shared(mesh, node.get()))
      sharedNodeAncestriesAndSignOrScore.emplace_back(node.get(),get_node_sign_or_score<T>(node.get()));

  return sharedNodeAncestriesAndSignOrScore;
}

template <typename T>
void receive_node_sign_or_score(CDMesh & cdmesh, stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&commSparse, &cdmesh](int procId)
  {
    auto & buffer = commSparse.recv_buffer(procId);
    SubElementChildNodeAncestry shared_edge_node(buffer);
    T signOrScore = 0;
    buffer.unpack(signOrScore);

    const SubElementNode * node = shared_edge_node.find_subelement_node(cdmesh);
    if (node)
      set_node_sign_or_score(node, signOrScore);
  });
}

template <typename T>
void sync_node_sign_or_score_on_local_constrained_nodes(CDMesh & cdmesh, const std::vector<std::pair<SubElementChildNodeAncestry,T>> & constrainedNodesAndSignOrScore, const std::vector<std::vector<int>> & owningProcsOfNodesInAncestries)
{
  STK_ThrowAssert(constrainedNodesAndSignOrScore.size() == owningProcsOfNodesInAncestries.size());

  for (size_t i=0; i<constrainedNodesAndSignOrScore.size(); ++i)
  {
    if (std::binary_search(owningProcsOfNodesInAncestries[i].begin(), owningProcsOfNodesInAncestries[i].end(), cdmesh.stk_bulk().parallel_rank()))
    {
      const auto & nodeAncestry = constrainedNodesAndSignOrScore[i].first;
      const T signOrScore = constrainedNodesAndSignOrScore[i].second;

      const SubElementNode * node = nodeAncestry.find_subelement_node(cdmesh);
      if (node)
        set_node_sign_or_score(node, signOrScore);
    }
  }
}

template <typename T>
void sync_node_sign_or_score_on_constrained_nodes(CDMesh & cdmesh,
    const stk::mesh::BulkData & mesh,
    const std::vector<std::unique_ptr<SubElementNode>> & nodes,
    const std::unordered_map<stk::mesh::EntityId, std::vector<stk::mesh::EntityId> > & periodicNodeIDMap)
{
  if (stk::is_true_on_all_procs(cdmesh.stk_bulk().parallel(), periodicNodeIDMap.empty()))
    return;

  const std::vector<std::pair<SubElementChildNodeAncestry,T>> constrainedNodeAncestriesAndSignOrScore = gather_constrained_node_ancestries_and_sign_or_score<T>(nodes, periodicNodeIDMap);
  const std::vector<std::vector<int>> owningProcsOfNodesInAncestries = determine_owning_procs_of_nodes_in_ancestries(mesh, constrainedNodeAncestriesAndSignOrScore);

  sync_node_sign_or_score_on_local_constrained_nodes(cdmesh, constrainedNodeAncestriesAndSignOrScore, owningProcsOfNodesInAncestries);

  if (mesh.parallel_size() < 2) return;

  stk::CommSparse commSparse(mesh.parallel());
  pack_node_data_for_node_ancestries(mesh, constrainedNodeAncestriesAndSignOrScore, owningProcsOfNodesInAncestries, commSparse);
  receive_node_sign_or_score<T>(cdmesh, commSparse);
}

template <typename T>
void sync_node_sign_or_score_on_shared_nodes(CDMesh & cdmesh,
    const stk::mesh::BulkData & mesh,
    const std::vector<std::unique_ptr<SubElementNode>> & nodes)
{
  if (mesh.parallel_size() < 2) return;

  const std::vector<std::pair<SubElementChildNodeAncestry,T>> sharedNodeAncestriesAndSignOrScore = gather_shared_node_ancestries_and_sign_or_score<T>(mesh, nodes);
  const std::vector<std::vector<int>> sharingProcsOfNodesInAncestries = determine_sharing_procs_of_nodes_in_ancestries(mesh, sharedNodeAncestriesAndSignOrScore);

  stk::CommSparse commSparse(mesh.parallel());
  pack_node_data_for_node_ancestries(mesh, sharedNodeAncestriesAndSignOrScore, sharingProcsOfNodesInAncestries, commSparse);
  receive_node_sign_or_score<T>(cdmesh, commSparse);
}

void
CDMesh::sync_node_signs_on_constrained_nodes()
{
  sync_node_sign_or_score_on_constrained_nodes<int>(*this, stk_bulk(), nodes, my_periodic_node_id_map);
}

void
CDMesh::sync_node_scores_on_constrained_nodes()
{
  sync_node_sign_or_score_on_constrained_nodes<double>(*this, stk_bulk(), nodes, my_periodic_node_id_map);
}

void
CDMesh::parallel_sync_node_signs_on_shared_nodes()
{
  sync_node_sign_or_score_on_shared_nodes<int>(*this, stk_bulk(), nodes);
}

void
CDMesh::parallel_sync_node_scores_on_shared_nodes()
{
  sync_node_sign_or_score_on_shared_nodes<double>(*this, stk_bulk(), nodes);
}

void
CDMesh::handle_hanging_children(const InterfaceID & interface)
{ /* %TRACE[ON]% */ Trace trace__("krino::CDMesh::handle_hanging_children(const InterfaceID & interface)"); /* %TRACE% */

  build_parallel_hanging_edge_nodes();

  for (auto && elem : elements)
  {
    if(krinolog.shouldPrint(LOG_DEBUG))
    {
      krinolog << "Handling hanging children Mesh_Element identifier=" << elem->entityId();
      krinolog << "\n";
    }
    elem->handle_hanging_children(*this, interface);
    if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "\n";
  }
}

//--------------------------------------------------------------------------------

const SubElementMeshNode *
CDMesh::get_mesh_node( stk::mesh::EntityId node_id ) const
{
  NodeMap::const_iterator it = mesh_node_map.find( node_id );
  return ( it != mesh_node_map.end() ) ? it->second : nullptr;
}

const SubElementMeshNode *
CDMesh::get_mesh_node(const SubElementNode * new_node) const
{
  const SubElementMeshNode * mesh_node = dynamic_cast<const SubElementMeshNode *>( new_node );
  if (nullptr != mesh_node)
  {
    return get_mesh_node(new_node->entityId());
  }
  return nullptr;
}

SubElementMeshNode *
CDMesh::add_managed_mesh_node( std::unique_ptr<SubElementMeshNode> node )
{
  SubElementMeshNode * ptr = node.get();
  mesh_node_map.insert( NodeMap::value_type( node->entityId(), ptr ) );
  add_managed_node(std::move(node));
  return ptr;
}

//--------------------------------------------------------------------------------

const SubElementNode *
CDMesh::create_mesh_node(
    const Mesh_Element * owner,
    const int lnn,
    stk::mesh::Entity nodeEntity )
{
  const stk::mesh::EntityId nodeId = stk_bulk().identifier(nodeEntity);
  const SubElementMeshNode * subnode  = get_mesh_node(nodeId);

  if ( nullptr == subnode )
  {
    const stk::math::Vector3d owner_coords = owner->get_node_parametric_coords( lnn );
    const double * global_coords_ptr = field_data<double>(get_coords_field(), nodeEntity);

    stk::math::Vector3d global_coords(global_coords_ptr, my_spatial_dim);

    std::unique_ptr<SubElementMeshNode> meshNode = std::make_unique<SubElementMeshNode>( owner, nodeEntity, nodeId, owner_coords, global_coords );
    subnode = add_managed_mesh_node(std::move(meshNode));
  }

  return subnode;
}

//--------------------------------------------------------------------------------

const SubElementNode *
CDMesh::create_edge_node( const Mesh_Element * owner,
                        const SubElementNode * parent1,
                        const SubElementNode * parent2,
                        const double position)
{
  const SubElementNode * subnode = SubElementNode::common_child({parent1, parent2});

  if (subnode == nullptr)
  {
    std::unique_ptr<SubElementNode> newNode = std::make_unique<SubElementEdgeNode>(owner, position, parent1, parent2);
    subnode = add_managed_node(std::move(newNode));

    parent1->add_child(subnode);
    parent2->add_child(subnode);
  }

  return subnode;
}

//--------------------------------------------------------------------------------

const SubElementNode *
CDMesh::create_midside_node( const Mesh_Element * owner,
                        const SubElementNode * parent1,
                        const SubElementNode * parent2,
                        stk::mesh::Entity entity)
{
  std::pair<const SubElementNode *,const SubElementNode *> parents =
      parent1<parent2 ? std::make_pair(parent1,parent2) : std::make_pair(parent2,parent1);

  const SubElementNode * subnode = nullptr;

  auto iter = my_midside_node_map.find(parents);
  if (iter != my_midside_node_map.end())
  {
    subnode = iter->second;
  }
  else
  {
    std::unique_ptr<SubElementNode> newNode = (stk_bulk().is_valid(entity)) ?
        std::make_unique<SubElementMidSideNode>(owner, parent1, parent2, entity, stk_bulk().identifier(entity)) :
        std::make_unique<SubElementMidSideNode>(owner, parent1, parent2);
    subnode = add_managed_node(std::move(newNode));
    my_midside_node_map[parents] = subnode;
  }

  return subnode;
}

//--------------------------------------------------------------------------------

const SubElementNode *
CDMesh::create_steiner_node( const Mesh_Element * owner,
    const NodeVec & parents,
    const std::vector<double> & weights )
{
  std::unique_ptr<SubElementNode> newNode = std::make_unique<SubElementSteinerNode>( owner, parents, weights );
  return add_managed_node(std::move(newNode));
}

//--------------------------------------------------------------------------------

const SubElementNode *
CDMesh::create_child_internal_or_face_node( const Mesh_Element * owner,
    const NodeVec & parents,
    const std::vector<double> & weights )
{
  const SubElementNode * subnode = SubElementNode::common_child(parents);

  if (subnode == nullptr)
  {
    std::unique_ptr<SubElementNode> newNode = std::make_unique<SubElementChildNode>( owner, parents, weights );
    subnode = add_managed_node(std::move(newNode));

    for (auto && parent : parents)
      parent->add_child(subnode);
  }

  return subnode;
}

//--------------------------------------------------------------------------------

void
CDMesh::create_subelement_mesh_entities(
    const Mesh_Element & elem,
    const std::vector<const SubElement *> conformal_subelems)
{
  stk::mesh::Entity parent_elem = elem.entity();
  for (auto && subelem : conformal_subelems)
  {
    const stk::mesh::PartVector & parent_parts = stk_bulk().bucket(parent_elem).supersets();
    const stk::topology parent_topology = stk_bulk().bucket(parent_elem).topology();
    stk::mesh::PartVector subelem_parts;
    determine_child_conformal_parts(parent_topology, parent_parts, subelem->get_phase(), subelem_parts);

    if (0 == subelem->entityId())
    {
      const stk::mesh::EntityId new_id = my_entity_id_pool.get_EntityId(stk::topology::ELEMENT_RANK);
      STK_ThrowAssert(!stk_bulk().is_valid(stk_bulk().get_entity(stk::topology::ELEMENT_RANK, new_id)));
      stk::mesh::Entity subelem_entity = stk_bulk().declare_element(new_id, subelem_parts);
      subelem->set_entity( stk_bulk(), subelem_entity );
      STK_ThrowAssert(stk_bulk().bucket(subelem_entity).topology() != stk::topology::INVALID_TOPOLOGY);

      const NodeVec & elem_nodes = subelem->get_nodes();
      for (unsigned n=0; n<elem_nodes.size(); ++n)
      {
        stk::mesh::Entity node = elem_nodes[n]->entity();
        stk_bulk().declare_relation( subelem_entity, node , n );
      }
    }
    else
    {
      stk::mesh::Entity subelem_entity = subelem->entity();
      stk_bulk().change_entity_parts(subelem_entity, subelem_parts, get_removable_parts(stk_bulk(), subelem_entity));
    }
  }
}

void
CDMesh::attach_existing_and_identify_missing_subelement_sides(
    const Mesh_Element & elem,
    const std::vector<const SubElement *> conformal_subelems,
    std::vector<SideDescription> & side_requests)
{
  stk::mesh::BulkData & stk_mesh = stk_bulk();
  const bool build_internal_sides = my_cdfem_support.use_internal_face_stabilization();

  for (auto && subelem : conformal_subelems)
  {
    const stk::topology topology = subelem->topology();
    const stk::mesh::Entity * elem_nodes = stk_bulk().begin_nodes(subelem->entity());

    for (unsigned s=0; s<topology.num_sides(); ++s)
    {
      const stk::topology side_topology = topology.side_topology(s);
      std::vector<stk::mesh::Entity> side_nodes(side_topology.num_nodes());
      topology.side_nodes(elem_nodes, s, side_nodes.data());

      std::vector<stk::mesh::Entity> sides;
      stk::mesh::get_entities_through_relations(stk_mesh, side_nodes, stk_meta().side_rank(), sides);

      if (sides.empty())
      {
        stk::mesh::Entity parent_side = find_entity_by_ordinal(stk_bulk(), elem.entity(), stk_meta().side_rank(), subelem->parent_side_id(s));
        const bool have_parent_side = stk_bulk().is_valid(parent_side);
        const bool is_internal_side = subelem->parent_side_id(s) == -1;

        if (have_parent_side || (is_internal_side && build_internal_sides))
        {
          static stk::mesh::PartVector empty_parts;
          const stk::mesh::PartVector & parent_parts = have_parent_side ? stk_bulk().bucket(parent_side).supersets() : empty_parts;

          // We have to make sure that pre-existing sideset parts are added to the side so that we
          // can figure out the correct conformal side parts during the second modification pass.
          stk::mesh::PartVector side_parts;
          determine_child_conformal_parts(side_topology, parent_parts, subelem->get_phase(), side_parts);
          if (is_internal_side)
          {
            side_parts.push_back(&get_internal_side_part());
          }

          side_requests.emplace_back(subelem->entity(), s, side_parts);
        }
      }
      else
      {
        STK_ThrowRequire(sides.size() == 1);
        attach_entity_to_element(stk_bulk(), stk_meta().side_rank(), sides[0], subelem->entity());
      }
    }
  }
}

bool
CDMesh::check_element_side_parts() const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::check_element_side_parts()"); /* %TRACE% */
  // This method requires aura to work correctly.
  if (!stk_bulk().is_automatic_aura_on())
  {
    // Skip check if we don't have aura
    return true;
  }

  bool success = true;
  stk::mesh::Selector active_locally_owned = aux_meta().active_locally_owned_selector();
  stk::mesh::BucketVector const& buckets = stk_bulk().get_buckets(stk::topology::ELEMENT_RANK, active_locally_owned);

  std::vector<stk::mesh::Entity> side_nodes;

  for (auto&& bucket : buckets)
  {
    const stk::topology topology = bucket->topology();
    const unsigned num_sides = topology.num_sides();
    for (auto&& elem : *bucket)
    {
      auto elem_nodes = stk_bulk().begin(elem, stk::topology::NODE_RANK);
      for (unsigned s=0; s<num_sides; ++s)
      {
        auto side_topology = topology.side_topology(s);
        side_nodes.resize(side_topology.num_nodes());
        topology.side_nodes(elem_nodes, s, side_nodes.data());

        if (!check_element_side_parts(side_nodes))
        {
          krinolog << "Side nodes: ";
          for(auto && node : side_nodes) krinolog << debug_entity(stk_bulk(), node) << stk::diag::dendl;

          krinolog << "Elements connected to side nodes: ";
          std::vector<stk::mesh::Entity> elems;
          stk::mesh::get_entities_through_relations(stk_bulk(), side_nodes, stk::topology::ELEMENT_RANK, elems);
          for(auto && touching_elem : elems) krinolog << debug_entity(stk_bulk(), touching_elem) << stk::diag::dendl;

          success = false;
        }
      }
    }
  }

  return success;
}

std::vector<unsigned> get_conformal_volume_part_ordinals(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, stk::mesh::Entity entity)
{
  std::vector<unsigned> conformalVolumeParts;

  for(auto && part : mesh.bucket(entity).supersets())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      if (phaseSupport.is_conformal(part))
      {
        conformalVolumeParts.push_back(part->mesh_meta_data_ordinal());
      }
    }
  }

  return conformalVolumeParts;
}

bool have_multiple_conformal_volume_parts_in_common(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const std::vector<stk::mesh::Entity> & sideNodes)
{
  const int numSideNodes = sideNodes.size();
  STK_ThrowRequire(numSideNodes > 0);

  std::vector<unsigned> commonConformalVolumeParts = get_conformal_volume_part_ordinals(mesh, phaseSupport, sideNodes[0]);

  for (int n=1; n<numSideNodes; ++n)
  {
    const std::vector<unsigned> nodeConformalVolumeParts = get_conformal_volume_part_ordinals(mesh, phaseSupport, sideNodes[n]);

    std::vector<unsigned> workingSet;
    workingSet.swap(commonConformalVolumeParts);
    std::set_intersection(workingSet.begin(),workingSet.end(),nodeConformalVolumeParts.begin(),nodeConformalVolumeParts.end(),std::back_inserter(commonConformalVolumeParts));

    if (commonConformalVolumeParts.empty())
      return false;
  }
  return commonConformalVolumeParts.size() > 1;
}

void
CDMesh::add_possible_interface_sides(std::vector<SideDescription> & sideRequests) const
{
  // This will add sides that *might be* interface sides.
  // Because this probes the nodes, it will add "keyhole" sides that aren't actually on an interface
  // This should be harmless, however, and avoids extra communication or requiring aura.

  stk::mesh::Selector active_locally_owned = aux_meta().active_locally_owned_selector();
  stk::mesh::BucketVector const& buckets = stk_bulk().get_buckets(stk::topology::ELEMENT_RANK, active_locally_owned);

  std::vector<stk::mesh::Entity> sideNodes;

  for (auto&& bucket : buckets)
  {
    const stk::topology topology = bucket->topology();
    const unsigned num_sides = topology.num_sides();
    for (auto&& elem : *bucket)
    {
      auto elem_nodes = stk_bulk().begin(elem, stk::topology::NODE_RANK);
      for (unsigned s=0; s<num_sides; ++s)
      {
        auto sideTopology = topology.side_topology(s);
        sideNodes.resize(sideTopology.num_nodes());
        topology.side_nodes(elem_nodes, s, sideNodes.data());

        const bool possibleInterfaceSide = have_multiple_conformal_volume_parts_in_common(stk_bulk(), my_phase_support, sideNodes);
        if (possibleInterfaceSide)
        {
          stk::mesh::PartVector sideParts{&stk_meta().get_topology_root_part(sideTopology)};
          sideRequests.emplace_back(elem, s, sideParts);
        }
      }
    }
  }
}


bool
CDMesh::check_element_side_parts(const std::vector<stk::mesh::Entity> & side_nodes) const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::check_element_side_parts(const std::vector<stk::mesh::Entity> & side_nodes)"); /* %TRACE% */

  // This method requires aura.
  STK_ThrowRequire(stk_bulk().is_automatic_aura_on());

  std::vector<stk::mesh::Entity> elems;
  stk::mesh::get_entities_through_relations(stk_bulk(), side_nodes, stk::topology::ELEMENT_RANK, elems);

  std::vector<const stk::mesh::Part *> conformal_volume_parts;
  for (auto&& elem : elems)
  {
    if (!stk_bulk().bucket(elem).member(get_active_part()))
    {
      continue;
    }
    auto& elem_parts = stk_bulk().bucket(elem).supersets();
    for(auto&& part : elem_parts)
    {
      if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && my_phase_support.is_conformal(part))
      {
        if (std::find(conformal_volume_parts.begin(), conformal_volume_parts.end(), part) == conformal_volume_parts.end())
        {
          conformal_volume_parts.push_back(part);
        }
      }
    }
  }

  if (conformal_volume_parts.empty())
  {
    return true;
  }

  if (conformal_volume_parts.size() > 2)
  {
    krinolog << "Expected to find 1 or 2 conformal side parts when examining side nodes: ";
    for (auto&& side_node : side_nodes) krinolog << stk_bulk().identifier(side_node) << " ";
    krinolog << " but instead found the parts: ";
    for (auto&& part : conformal_volume_parts) krinolog << part->name() << " ";
    krinolog << stk::diag::dendl;
    return false;
  }

  std::vector<PhaseTag> side_phases(conformal_volume_parts.size());
  for (unsigned iphase = 0; iphase<side_phases.size(); ++iphase)
  {
    side_phases[iphase] = my_phase_support.get_iopart_phase(*conformal_volume_parts[iphase]);
    STK_ThrowRequire(!side_phases[iphase].empty());
  }

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_entities_through_relations(stk_bulk(), side_nodes, stk_meta().side_rank(), sides);

  if (conformal_volume_parts.size() == 2 && side_phases[0] != side_phases[1])
  {
    stk::mesh::PartVector conformal_side_parts;
    const stk::mesh::Part * conformal_side_part = my_phase_support.find_interface_part(*conformal_volume_parts[0], *conformal_volume_parts[1]);
    if (nullptr != conformal_side_part) conformal_side_parts.push_back(const_cast<stk::mesh::Part *>(conformal_side_part));
    conformal_side_part = my_phase_support.find_interface_part(*conformal_volume_parts[1], *conformal_volume_parts[0]);
    if (nullptr != conformal_side_part) conformal_side_parts.push_back(const_cast<stk::mesh::Part *>(conformal_side_part));

    if (!conformal_side_parts.empty())
    {
      // Check that side exists and has conformal side parts
      if (sides.size() != 1)
      {
        krinolog << "Expected to find 1 conformal side, but instead found " << sides.size() << " when examining side nodes: ";
        for (auto&& side_node : side_nodes) krinolog << stk_bulk().identifier(side_node) << " ";
        krinolog << " with conformal volume parts: ";
        for (auto&& part : conformal_volume_parts) krinolog << part->name() << " ";
        krinolog << stk::diag::dendl;
        return false;
      }
      else
      {
        auto& side_bucket = stk_bulk().bucket(sides[0]);
        if (!side_bucket.member_all(conformal_side_parts))
        {
          krinolog << "Side " << stk_bulk().identifier(sides[0]) << " is missing at least one of the conformal side parts: ";
          for (auto&& part : conformal_side_parts) krinolog << part->name() << " ";
          krinolog << ", actual parts: ";
          auto& side_parts = side_bucket.supersets();
          for (auto&& part : side_parts) krinolog << part->name() << " ";
          krinolog << stk::diag::dendl;
          return false;
        }
      }
    }
  }
  else
  {
    // Check that if the side exists, then it does not have any interface sides
    if (sides.size() > 1)
    {
      krinolog << "Expected to find 0 or 1 side, but instead found " << sides.size() << " when examining side nodes: ";
      for (auto&& side_node : side_nodes) krinolog << stk_bulk().identifier(side_node) << " ";
      krinolog << " with conformal volume parts: ";
      for (auto&& part : conformal_volume_parts) krinolog << part->name() << " ";
      krinolog << stk::diag::dendl;
      return false;
    }
    else if (sides.size() == 1)
    {
      const stk::mesh::PartVector & existing_side_parts = stk_bulk().bucket(sides[0]).supersets();
      for(auto && sidePart : existing_side_parts)
      {
        if(sidePart->primary_entity_rank() == stk_meta().side_rank() && my_phase_support.is_interface(sidePart))
        {
          krinolog << "Side " << stk_bulk().identifier(sides[0]) << " has an erroneous interface part " << sidePart->name() << "." << stk::diag::dendl;
          return false;
        }
      }
    }
  }

  return true;
}

void
CDMesh::update_element_side_parts()
{
  // This method makes sure the correct conformal side parts are on the element sides
  stk::mesh::Selector locally_owned = get_locally_owned_part();

  std::vector< stk::mesh::Entity> sides;
  stk::mesh::get_selected_entities( locally_owned, stk_bulk().buckets( stk_bulk().mesh_meta_data().side_rank() ), sides );

  stk::mesh::PartVector addParts;
  stk::mesh::PartVector removeParts;
  std::vector<stk::mesh::PartVector> batchAddParts;
  std::vector<stk::mesh::PartVector> batchRemoveParts;
  batchAddParts.reserve(sides.size());
  batchRemoveParts.reserve(sides.size());

  for (auto && side : sides)
  {
    determine_element_side_parts(side, addParts, removeParts);
    batchAddParts.push_back(addParts);
    batchRemoveParts.push_back(removeParts);
  }

  stk_bulk().batch_change_entity_parts(sides, batchAddParts, batchRemoveParts);
}

void
CDMesh::determine_element_side_parts(const stk::mesh::Entity side, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const
{
  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    krinolog << "Analyzing side " << debug_entity_1line(stk_bulk(), side) << "\n";
    for (auto sideNode : StkMeshEntities{stk_bulk().begin_nodes(side), stk_bulk().end_nodes(side)})
      krinolog << "  " << debug_entity_1line(stk_bulk(), sideNode, true) << "\n";
    for (auto sideElem : StkMeshEntities{stk_bulk().begin_elements(side), stk_bulk().end_elements(side)})
      krinolog << "  " << debug_entity_1line(stk_bulk(), sideElem, true) << "\n";
    krinolog << stk::diag::dendl;
  }

  add_parts.clear();
  remove_parts.clear();

  std::vector<const stk::mesh::Part *> volume_parts;
  std::vector<const stk::mesh::Part *> conformal_volume_parts;
  std::vector<const stk::mesh::Part *> nonconformal_volume_parts;
  const stk::mesh::PartVector & existing_side_parts = stk_bulk().bucket(side).supersets();
  for(stk::mesh::PartVector::const_iterator part_iter = existing_side_parts.begin(); part_iter != existing_side_parts.end(); ++part_iter)
  {
    const stk::mesh::Part * const side_part = *part_iter;
    if (side_part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      if (my_phase_support.is_conformal(side_part))
      {
        conformal_volume_parts.push_back(side_part);
      }
      if (my_phase_support.is_nonconformal(side_part))
      {
        nonconformal_volume_parts.push_back(side_part);
      }
      else if (stk::io::is_part_io_part(*side_part) &&
          !stk::io::is_part_assembly_io_part(*side_part))
      {
        volume_parts.push_back(side_part);
      }
    }
  }

  if (volume_parts.size() > 2)
  {
    krinolog << "Found side with more than 2 volume parts:" << stk::diag::dendl;
    krinolog << debug_entity_1line(stk_bulk(), side) << stk::diag::dendl;
    for (auto elem : StkMeshEntities{stk_bulk().begin_elements(side), stk_bulk().end_elements(side)})
      krinolog << " " << debug_entity_1line(stk_bulk(), elem) << stk::diag::dendl;
  }

  STK_ThrowRequire(volume_parts.size() <= 2); // Can be zero for inactive elements supporting a face

  if (conformal_volume_parts.empty())
  {
    /* There are two possible cases where no conformal volume parts are found:
     *   1) This side is part of a surface that does not touch any blocks that are being decomposed.
     *      Only the active parts for these sides should be updated.
     *   2) This side is a parent side that should be deactivated and moved to the nonconformal part.
     *      These sides will have at least 1 nonconformal volume part from the parent volume element.
     */
    if(nonconformal_volume_parts.empty())
    {
      if(element_side_should_be_active(side))
      {
        add_parts.push_back(&aux_meta().active_part());
      }
      else
      {
        remove_parts.push_back(&aux_meta().active_part());
      }
    }
    else
    {
      determine_nonconformal_parts(side, add_parts, remove_parts);
    }
  }

  if (volume_parts.size() == 2)
  {
    add_parts.push_back(&get_block_boundary_part());
  }
  else
  {
    remove_parts.push_back(&get_block_boundary_part());
  }

  if (conformal_volume_parts.empty())
  {
    return;
  }

  STK_ThrowRequire(conformal_volume_parts.size() == 1 || conformal_volume_parts.size() == 2);

  std::vector<PhaseTag> side_phases(conformal_volume_parts.size());
  for (unsigned iphase = 0; iphase<side_phases.size(); ++iphase)
  {
    side_phases[iphase] = my_phase_support.get_iopart_phase(*conformal_volume_parts[iphase]);
    STK_ThrowRequire(!side_phases[iphase].empty());
  }

  if (conformal_volume_parts.size() == 2 && side_phases[0] != side_phases[1])
  {
    // interface side, add interface parts
    stk::mesh::Part * conformal_side0 = const_cast<stk::mesh::Part *>(my_phase_support.find_interface_part(*conformal_volume_parts[0], *conformal_volume_parts[1]));
    if (nullptr != conformal_side0) add_parts.push_back(conformal_side0);
    stk::mesh::Part * conformal_side1 = const_cast<stk::mesh::Part *>(my_phase_support.find_interface_part(*conformal_volume_parts[1], *conformal_volume_parts[0]));
    if (nullptr != conformal_side1) add_parts.push_back(conformal_side1);
  }

  for (auto && side_phase : side_phases)
  {
    determine_conformal_parts(existing_side_parts, stk_meta().side_rank(), side_phase, add_parts, remove_parts);
  }

  if(element_side_should_be_active(side))
  {
    add_parts.push_back(&aux_meta().active_part());
  }
  else
  {
    remove_parts.push_back(&aux_meta().active_part());
  }
}

bool
CDMesh::element_side_should_be_active(const stk::mesh::Entity side) const
{
  auto num_elems = stk_bulk().num_connectivity(side, stk::topology::ELEMENT_RANK);
  const auto * touching_elems = stk_bulk().begin(side, stk::topology::ELEMENT_RANK);
  auto & active_part = aux_meta().active_part();
  bool active = false;
  for(unsigned i=0; i < num_elems; ++i)
  {
    if(stk_bulk().bucket(touching_elems[i]).member(active_part))
    {
      active = true;
      break;
    }
  }
  return active;
}
void
CDMesh::handle_single_coincident_subelement(const Mesh_Element & elem, const SubElement * subelem, std::vector<SideDescription> & side_requests)
{
  stk::mesh::Entity elem_entity = elem.entity();
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "single coincident subelement for elem #" << stk_bulk().identifier(elem_entity) << " with phase " << subelem->get_phase() << stk::diag::dendl;
  subelem->set_entity( stk_bulk(), elem_entity );
  stk::mesh::PartVector add_parts;
  stk::mesh::PartVector remove_parts;
  determine_conformal_parts(elem_entity, subelem->get_phase(), add_parts, remove_parts);

  add_parts.push_back(&get_active_part());
  remove_parts.push_back(&get_parent_part());
  remove_parts.push_back(&get_child_part());

  stk_bulk().change_entity_parts(elem_entity, add_parts, remove_parts);

  std::vector<const SubElement *> subelem_vec(1, subelem);
  attach_existing_and_identify_missing_subelement_sides(elem, subelem_vec, side_requests);
}

//--------------------------------------------------------------------------------

void
CDMesh::generate_sorted_child_elements()
{
  child_elements.clear();

  for (const auto & elem : elements)
  {
    if (elem->have_subelements())
    {
      std::vector<const SubElement *> conformal_subelems;
      elem->get_subelements(conformal_subelems);

      for (auto && subelem : conformal_subelems)
      {
        child_elements.push_back(subelem);
      }
    }
  }

  std::sort(child_elements.begin(),child_elements.end(), ElementObj::is_less);
}

//--------------------------------------------------------------------------------

const SubElement *
CDMesh::find_child_element(stk::mesh::Entity elem_mesh_obj) const
{
  // Ugh
  CDMesh *const mesh = const_cast<CDMesh*>(this);

  if (child_elements.empty())
  {
    mesh->generate_sorted_child_elements();
  }

  const stk::mesh::EntityId elem_id = stk_bulk().identifier(elem_mesh_obj);
  auto lb_cmp = [](const ElementObj * elem, stk::mesh::EntityId target_id) { return elem->entityId() < target_id; };
  auto first = std::lower_bound(child_elements.begin(),child_elements.end(), elem_id, lb_cmp);

  if (first != child_elements.end() && (*first)->entityId() == elem_id)
  {
    return dynamic_cast<const SubElement *>(*first);
  }
  return nullptr;
}

//--------------------------------------------------------------------------------

stk::mesh::Entity
CDMesh::get_parent_element(stk::mesh::Entity elem_entity) const
{
  std::set<stk::mesh::Entity> parent_elem_node_set;

  const stk::mesh::Entity * elem_nodes = stk_bulk().begin_nodes(elem_entity);
  const unsigned num_base_elem_nodes = stk_bulk().bucket(elem_entity).topology().base().num_nodes();

  for (unsigned inode=0; inode<num_base_elem_nodes; ++inode)
    recursively_fill_parent_nodes(stk_bulk(), elem_nodes[inode], get_parent_node_ids_field(), parent_elem_node_set);

  const std::vector<stk::mesh::Entity> parent_elem_nodes(parent_elem_node_set.begin(), parent_elem_node_set.end());
  std::vector<stk::mesh::Entity> parent_elems;
  stk::mesh::get_entities_through_relations(stk_bulk(), parent_elem_nodes, stk::topology::ELEMENT_RANK, parent_elems);

  STK_ThrowAssert(parent_elems.size() <= 1);

  if (parent_elems.empty())
  {
    krinolog << "Did not find parent element for element \n" << debug_entity(stk_bulk(), elem_entity) << stk::diag::dendl;
    return stk::mesh::Entity();
  }
  else
  {
    return parent_elems.front();
  }
}

//--------------------------------------------------------------------------------

bool
CDMesh::get_parent_child_coord_transformation(stk::mesh::Entity elem_mesh_obj, double * dParentdChild) const
{
  STK_ThrowAssert(my_cdfem_support.use_nonconformal_element_size());

  const SubElement * subelem = find_child_element(elem_mesh_obj);

  if (nullptr == subelem)
  {
    krinolog << "did not find element " << stk_bulk().identifier(elem_mesh_obj) << stk::diag::dendl;
    return false;
  }

  subelem->get_owner_coord_transform(dParentdChild);
  return true;
}

void
CDMesh::get_parent_nodes_and_weights(stk::mesh::Entity child, stk::mesh::Entity & parent0, stk::mesh::Entity & parent1, double & position) const
{
  // Really slow!
  auto id = stk_bulk().identifier(child);
  auto find_existing = std::find_if(nodes.begin(), nodes.end(),
      [id](const std::unique_ptr<krino::SubElementNode> & compare)
      { return compare->entityId() == id; });
  STK_ThrowAssert(find_existing != nodes.end());
  STK_ThrowAssert(dynamic_cast<krino::SubElementEdgeNode *>(find_existing->get()) != nullptr);
  const krino::SubElementEdgeNode& edge_node = dynamic_cast<krino::SubElementEdgeNode &>(*find_existing->get());
  krino::NodeVec edge_node_parents = edge_node.get_parents();
  position = edge_node.get_position();
  parent0 = edge_node_parents[0]->entity();
  parent1 = edge_node_parents[1]->entity();
}

//--------------------------------------------------------------------------------

std::function<double(stk::mesh::Entity)> build_get_local_length_scale_for_side_function(const CDMesh & cdmesh)
{
  const stk::mesh::Selector elementSelector = selectUnion(cdmesh.get_phase_support().get_conformal_parts()) & cdmesh.get_active_part() & cdmesh.get_locally_owned_part();

  auto get_length_scale_for_side =
      [&cdmesh,elementSelector](stk::mesh::Entity side)
      {
        const stk::mesh::BulkData & mesh = cdmesh.stk_bulk();
        double minElemVolume = 0.;
        for (auto elem : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
        {
          if (elementSelector(mesh.bucket(elem)))
          {
            stk::mesh::Entity volumeElement = cdmesh.get_cdfem_support().use_nonconformal_element_size() ? cdmesh.get_parent_element(elem) : elem;
            STK_ThrowRequire(cdmesh.stk_bulk().is_valid(volumeElement));
            const double elemVol = ElementObj::volume( mesh, volumeElement, cdmesh.get_coords_field() );
            if (minElemVolume == 0. || elemVol < minElemVolume)
              minElemVolume = elemVol;
          }
        }
        double lengthScale = 0.;
        if (minElemVolume > 0.)
        {
          const double invDim = 1.0 / mesh.mesh_meta_data().spatial_dimension();
          lengthScale = std::pow(minElemVolume, invDim);
        }

        return lengthScale;
      };
  return get_length_scale_for_side;
}

std::function<double(stk::mesh::Entity)> build_get_constant_length_scale_for_side_function(const double lengthScale)
{
  auto get_length_scale_for_side =
      [lengthScale](stk::mesh::Entity side)
      {
        return lengthScale;
      };
  return get_length_scale_for_side;
}

std::vector<stk::mesh::Entity> get_unique_owned_volume_elements_using_sides(const CDMesh & cdmesh, const stk::mesh::Selector & interfaceSideSelector)
{
  // Not exactly cheap
  const stk::mesh::BulkData & mesh = cdmesh.stk_bulk();
  const stk::mesh::Selector elementSelector = selectUnion(cdmesh.get_phase_support().get_conformal_parts()) & cdmesh.get_active_part() & cdmesh.get_locally_owned_part();

  std::vector<stk::mesh::Entity> volumeElements;
  for( auto&& bucket : mesh.get_buckets(mesh.mesh_meta_data().side_rank(), interfaceSideSelector) )
  {
    for (auto && side : *bucket)
    {
      for (auto elem : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
      {
        if (elementSelector(mesh.bucket(elem)))
        {
          stk::mesh::Entity volumeElement = cdmesh.get_cdfem_support().use_nonconformal_element_size() ? cdmesh.get_parent_element(elem) : elem;
          volumeElements.push_back(volumeElement);
        }
      }
    }
  }
  stk::util::sort_and_unique(volumeElements);
  return volumeElements;
}

double compute_L1_norm_of_side_length_scales(const CDMesh & cdmesh, const stk::mesh::Selector & interfaceSideSelector)
{
  const std::vector<stk::mesh::Entity> elementsInNorm = get_unique_owned_volume_elements_using_sides(cdmesh, interfaceSideSelector);

  const double invDim = 1.0 / cdmesh.spatial_dim();

  double sumLengths = 0.;
  for (auto elem : elementsInNorm)
  {
    const double elemVolume = ElementObj::volume( cdmesh.stk_bulk(), elem, cdmesh.get_coords_field() );
    sumLengths += std::pow(elemVolume, invDim);
  }

  const double sumCount = elementsInNorm.size();

  const std::array<double,2> localSum{sumLengths, sumCount};
  std::array<double,2> globalSum;
  stk::all_reduce_sum(cdmesh.stk_bulk().parallel(), localSum.data(), globalSum.data(), localSum.size());
  return globalSum[0]/globalSum[1];
}

stk::math::Vector3d get_side_average_of_vector(const stk::mesh::BulkData& mesh,
    const FieldRef vectorField,
    const stk::mesh::Entity side)
{
  const int spatialDim = mesh.mesh_meta_data().spatial_dimension();

  stk::math::Vector3d avg{stk::math::Vector3d::ZERO};
  int numNodes = 0;
  for (auto node : StkMeshEntities{mesh.begin_nodes(side), mesh.end_nodes(side)})
  {
    double * vectorPtr = field_data<double>(vectorField, node);
    if(nullptr != vectorPtr)
    {
      const stk::math::Vector3d vec(vectorPtr, spatialDim);
      avg += vec;
      ++numNodes;
    }
  }
  if (numNodes > 0)
    avg /= numNodes;

  return avg;
}

stk::math::Vector3d get_side_average_of_vector_difference(const stk::mesh::BulkData& mesh,
    const FieldRef addVectorField,
    const FieldRef subtractVectorField,
    const stk::mesh::Entity side)
{
  const int spatialDim = mesh.mesh_meta_data().spatial_dimension();

  stk::math::Vector3d avg{stk::math::Vector3d::ZERO};
  int numNodes = 0;
  for (auto node : StkMeshEntities{mesh.begin_nodes(side), mesh.end_nodes(side)})
  {
    double * addVectorPtr = field_data<double>(addVectorField, node);
    double * subtractVectorPtr = field_data<double>(subtractVectorField, node);
    if(nullptr != addVectorPtr && nullptr != subtractVectorPtr)
    {
      const stk::math::Vector3d addVec(addVectorPtr, spatialDim);
      const stk::math::Vector3d subtractVec(subtractVectorPtr, spatialDim);
      avg += addVec - subtractVec;
      ++numNodes;
    }
  }
  if (numNodes > 0)
    avg /= numNodes;

  return avg;
}

std::function<stk::math::Vector3d(stk::mesh::Entity)> build_get_side_displacement_from_cdfem_displacements_function(const stk::mesh::BulkData& mesh, const FieldRef cdfemDisplacementsField)
{
  auto get_element_size =
      [&mesh, cdfemDisplacementsField](stk::mesh::Entity side)
      {
        return get_side_average_of_vector(mesh, cdfemDisplacementsField, side);
      };
  return get_element_size;
}

std::function<stk::math::Vector3d(stk::mesh::Entity)> build_get_side_displacement_from_change_in_cdfem_displacements_function(const stk::mesh::BulkData& mesh, const FieldRef newCdfemDisplacementsField)
{
  auto get_element_size =
      [&mesh, newCdfemDisplacementsField](stk::mesh::Entity side)
      {
        return get_side_average_of_vector_difference(mesh, newCdfemDisplacementsField, newCdfemDisplacementsField.field_state(stk::mesh::StateOld), side);
      };
  return get_element_size;
}

std::function<stk::math::Vector3d(stk::mesh::Entity)> build_get_side_displacement_from_velocity_function(const stk::mesh::BulkData& mesh, const FieldRef velocity, const double dt)
{
  auto get_element_size =
      [&mesh, velocity, dt](stk::mesh::Entity side)
      {
        return dt*get_side_average_of_vector(mesh, velocity, side);
      };
  return get_element_size;
}

double get_side_cdfem_cfl(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    const std::function<stk::math::Vector3d(stk::mesh::Entity)> & get_side_displacement,
    const std::function<double(stk::mesh::Entity)> & get_length_scale_for_side,
    stk::mesh::Entity side)
{
  const stk::math::Vector3d sideCDFEMDisplacement = get_side_displacement(side);
  const stk::math::Vector3d sideNormal = get_side_normal(mesh, coordsField, side);
  const double sideNormalDisplacement = std::abs(Dot(sideCDFEMDisplacement, sideNormal));

  const double sideLengthScale = get_length_scale_for_side(side);

  return (sideLengthScale == 0.) ? 0. : sideNormalDisplacement/sideLengthScale;
}

double CDMesh::compute_cdfem_cfl(const Interface_CFL_Length_Scale lengthScaleType, const std::function<stk::math::Vector3d(stk::mesh::Entity)> & get_side_displacement) const
{
  stk::diag::TimeBlock timer__(my_timer_compute_CFL);

  const stk::mesh::Selector interfaceSideSelector = my_phase_support.get_all_conformal_surfaces_selector();

  get_coords_field().field().sync_to_host();

  std::function<double(stk::mesh::Entity)> get_length_scale_for_side;
  if (lengthScaleType == CONSTANT_LENGTH_SCALE)
  {
    get_length_scale_for_side = build_get_constant_length_scale_for_side_function(my_cdfem_support.get_constant_length_scale_for_interface_CFL());
  }
  else if (lengthScaleType == LOCAL_LENGTH_SCALE)
  {
    get_length_scale_for_side = build_get_local_length_scale_for_side_function(*this);
  }
  else
  {
    STK_ThrowRequire(lengthScaleType == L1_NORM_LENGTH_SCALE);
    const double lengthScaleNorm = compute_L1_norm_of_side_length_scales(*this, interfaceSideSelector);
    krinolog << "Using L1 Norm length scale " << lengthScaleNorm << " to compute Interface CFL." << stk::diag::dendl;
    get_length_scale_for_side = build_get_constant_length_scale_for_side_function(lengthScaleNorm);
  }

  double cfl = 0.;
  for( auto&& bucket : stk_bulk().get_buckets(stk_bulk().mesh_meta_data().side_rank(), interfaceSideSelector) )
  {
    for (auto && side : *bucket)
    {
      const double sideCFL = get_side_cdfem_cfl(stk_bulk(), get_coords_field(), get_side_displacement, get_length_scale_for_side, side);
      if (sideCFL > 0.)
        cfl = std::max(cfl, sideCFL);
    }
  }

  const double localCFL = cfl;
  stk::all_reduce_max(stk_bulk().parallel(), &localCFL, &cfl, 1);

  return cfl;
}

double CDMesh::compute_cdfem_displacement_cfl() const
{
  get_cdfem_displacements_field().field().sync_to_host();

  auto get_side_displacement = build_get_side_displacement_from_cdfem_displacements_function(stk_bulk(), get_cdfem_displacements_field());

  return compute_cdfem_cfl(my_cdfem_support.get_length_scale_type_for_interface_CFL(), get_side_displacement);
}

double CDMesh::compute_non_rebased_cdfem_displacement_cfl() const
{
  auto get_side_displacement = build_get_side_displacement_from_change_in_cdfem_displacements_function(stk_bulk(), get_cdfem_displacements_field());

  return compute_cdfem_cfl(LOCAL_LENGTH_SCALE, get_side_displacement);
}

double CDMesh::compute_interface_velocity_cfl(const FieldRef velocityField, const double dt) const
{
  velocityField.field().sync_to_host();

  auto get_side_displacement = build_get_side_displacement_from_velocity_function(stk_bulk(), velocityField, dt);

  return compute_cdfem_cfl(my_cdfem_support.get_length_scale_type_for_interface_CFL(), get_side_displacement);
}

void
CDMesh::update_adaptivity_parent_entities()
{
  // The CDFEM child and parent entities have been updated with the new interface locations and phases.
  // The inactive adaptivity parent entities now need to be updated similarly.  Otherwise, when they
  // are unrefined, they would become active with the wrong block parts.  Similarly, when we update the
  // side parts, we want the elements (even the inactive ones) to have the correct parts.

  if (myRefinementSupport.get_interface_maximum_refinement_level() <= 0 || !myRefinementSupport.has_non_interface_conforming_refinement())
  {
    return;
  }

  const RefinementInterface & refinement = myRefinementSupport.get_non_interface_conforming_refinement();
  stk::mesh::BulkData& mesh = stk_bulk();

  stk::mesh::PartVector add_parts;
  stk::mesh::PartVector remove_parts;

  const std::vector<std::pair<stk::mesh::Entity,unsigned>> parentsAndElementPart = get_owned_adaptivity_parents_and_their_element_part(mesh, refinement, my_phase_support);

  for( auto&& parentAndElementPart : parentsAndElementPart )
  {
    const stk::mesh::Entity parent = parentAndElementPart.first;
    const stk::mesh::Part & targetElementPart = stk_meta().get_part(parentAndElementPart.second);
    const stk::mesh::Part & currentElementPart = find_element_part(mesh, parent);

    if (currentElementPart.mesh_meta_data_ordinal() != targetElementPart.mesh_meta_data_ordinal())
    {
      if (my_phase_support.is_nonconformal(&targetElementPart))
      {
        determine_nonconformal_parts(parent, add_parts, remove_parts);
        const auto& add_it = std::find(add_parts.begin(), add_parts.end(), &get_parent_part());
        if (add_it != add_parts.end())
        {
          add_parts.erase(add_it);
        }
      }
      else
      {
        const PhaseTag & parent_phase = my_phase_support.get_iopart_phase(targetElementPart);
        determine_conformal_parts(parent, parent_phase, add_parts, remove_parts);
        remove_parts.push_back(&get_parent_part());
        remove_parts.push_back(&get_child_part());
      }
      mesh.change_entity_parts(parent, add_parts, remove_parts);
    }
  }
}

void
CDMesh::update_uncut_element(const Mesh_Element & elem)
{
  const stk::mesh::Entity elem_entity = elem.entity();
  stk::mesh::BulkData & stk_mesh = stk_bulk();
  if (stk_mesh.bucket(elem_entity).member(get_parent_part()) || elem_io_part_changed(elem))
  {
    stk::mesh::PartVector add_parts;
    stk::mesh::PartVector remove_parts;
    determine_conformal_parts(elem_entity, elem.get_phase(), add_parts, remove_parts);
    add_parts.push_back(&get_active_part());
    remove_parts.push_back(&get_parent_part());
    remove_parts.push_back(&get_child_part());

    stk_mesh.change_entity_parts(elem_entity, add_parts, remove_parts);
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::create_node_entities()
{
  stk::mesh::BulkData& stk_mesh = stk_bulk();

  std::vector<stk::mesh::Entity*> nodeParents;
  std::vector<ChildNodeRequest> childNodeRequests;
  std::vector<ChildNodeRequest> midSideNodeRequests;
  for (auto && node : nodes)
  {
    if (!stk_bulk().is_valid(node->entity()))
    {
      node->fill_parent_entity_pointers(nodeParents);
      if (nullptr != dynamic_cast<const SubElementMidSideNode *>(node.get()))
        midSideNodeRequests.push_back(ChildNodeRequest(nodeParents, &node->entity()));
      else
        childNodeRequests.push_back(ChildNodeRequest(nodeParents, &node->entity()));
    }
  }

  auto generate_new_ids = [&](stk::topology::rank_t entityRank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds)
  {
     EntityIdPool::generate_new_ids(stk_mesh, entityRank, numIdsNeeded, requestedIds, my_aux_meta.get_assert_32bit_flag(), my_aux_meta.get_force_64bit_flag());
  };

  stk::mesh::PartVector node_parts = {&aux_meta().active_part(),
    &get_child_edge_node_part(),
    &stk_meta().get_topology_root_part(stk::topology::NODE)
  };
  batch_create_child_nodes(stk_mesh, childNodeRequests, node_parts, generate_new_ids);

  stk::mesh::PartVector higher_order_node_parts = {&aux_meta().active_part(),
    &stk_meta().get_topology_root_part(stk::topology::NODE)
  };
  batch_create_child_nodes(stk_mesh, midSideNodeRequests, higher_order_node_parts, generate_new_ids);

  // Since batch_create_child_nodes took pointers to the entities, the entityIds were not updated, Ugh.
  for (auto&& node : nodes)
  {
    node->set_entityId_from_entity(stk_mesh);
    if(krinolog.shouldPrint(LOG_DEBUG))
    {
      if (!node->is_mesh_node())
        krinolog << "NODE ID : " << node->entityId() << ": ancestry: [" << node->get_ancestry() << "]" << stk::diag::dendl;
    }
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::create_element_and_side_entities(std::vector<SideDescription> & side_requests)
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::create_element_and_side_entities(void)"); /* %TRACE% */

  // Count how many we need to set pool size
  unsigned num_local_subelems = 0;
  for (auto && elem : elements)
  {
    if (elem->have_subelements())
    {
      std::vector<const SubElement *> conformal_subelems;
      elem->get_subelements(conformal_subelems);
      for (auto && subelem : conformal_subelems)
        if (0 == subelem->entityId())
          ++num_local_subelems;
    }
  }

  my_entity_id_pool.reserve(stk::topology::ELEMENT_RANK, num_local_subelems, my_aux_meta.get_assert_32bit_flag(), my_aux_meta.get_force_64bit_flag());

  stk::mesh::PartVector add_parts;
  stk::mesh::PartVector remove_parts;

  for (const auto & elem : elements)
  {
    if (elem->have_subelements())
    {
      std::vector<const SubElement *> conformal_subelems;
      elem->get_subelements(conformal_subelems);

      // check for corner case of a single subelement that is coincident with parent
      if (elem->is_single_coincident())
      {
        handle_single_coincident_subelement(*elem, conformal_subelems[0], side_requests);
      }
      else
      {
        create_subelement_mesh_entities(*elem, conformal_subelems);
        attach_existing_and_identify_missing_subelement_sides(*elem, conformal_subelems, side_requests);

        determine_nonconformal_parts(elem->entity(), add_parts, remove_parts);
        stk_bulk().change_entity_parts(elem->entity(), add_parts, remove_parts);
      }
    }
    else
    {
      update_uncut_element(*elem);
    }
  }

  update_adaptivity_parent_entities();
}

static void delete_all_child_elements(stk::mesh::BulkData & mesh)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());

  stk::mesh::Selector childSelector = cdfemSupport.get_child_part() & !cdfemSupport.get_parent_part();
  std::vector<stk::mesh::Entity> childElems;
  stk::mesh::get_selected_entities( childSelector, mesh.buckets( stk::topology::ELEMENT_RANK ), childElems );

  stk::mesh::destroy_elements(mesh, childElems, mesh.mesh_meta_data().universal_part());
}

void append_part_changes_to_reset_entities_to_original_undecomposed_state(const stk::mesh::BulkData & mesh,
  const stk::mesh::EntityRank entityRank,
  const Phase_Support & phaseSupport,
  stk::mesh::Part & childPart,
  stk::mesh::Part & parentPart,
  stk::mesh::Part & activePart,
  std::vector<stk::mesh::Entity> & batchEntities,
  std::vector<stk::mesh::PartVector> & batchAddParts,
  std::vector<stk::mesh::PartVector> & batchRemoveParts)
{
  stk::mesh::PartVector bucketAddParts;
  stk::mesh::PartVector bucketRemoveParts;

  for (auto * bucketPtr : mesh.get_buckets(entityRank, mesh.mesh_meta_data().locally_owned_part()))
  {
    determine_original_undecomposed_part_changes_for_entities(mesh, *bucketPtr, phaseSupport, childPart, parentPart, activePart, bucketAddParts, bucketRemoveParts);
    if (!bucketAddParts.empty() || !bucketRemoveParts.empty())
    {
      batchEntities.insert(batchEntities.end(), bucketPtr->begin(), bucketPtr->end());
      batchAddParts.insert(batchAddParts.end(), bucketPtr->size(), bucketAddParts);
      batchRemoveParts.insert(batchRemoveParts.end(), bucketPtr->size(), bucketRemoveParts);
    }
  }
}

void
CDMesh::reset_mesh_to_original_undecomposed_state(stk::mesh::BulkData & mesh)
{
  delete_all_child_elements(mesh);

  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());

  std::vector<stk::mesh::Entity> batchEntities;
  std::vector<stk::mesh::PartVector> batchAddParts;
  std::vector<stk::mesh::PartVector> batchRemoveParts;

  append_part_changes_to_reset_entities_to_original_undecomposed_state(mesh, stk::topology::ELEMENT_RANK, phaseSupport, cdfemSupport.get_child_part(), cdfemSupport.get_parent_part(), auxMeta.active_part(), batchEntities, batchAddParts, batchRemoveParts);
  append_part_changes_to_reset_entities_to_original_undecomposed_state(mesh, mesh.mesh_meta_data().side_rank(), phaseSupport, cdfemSupport.get_child_part(), cdfemSupport.get_parent_part(), auxMeta.active_part(), batchEntities, batchAddParts, batchRemoveParts);

  mesh.batch_change_entity_parts(batchEntities, batchAddParts, batchRemoveParts);

  if (the_new_mesh)
    the_new_mesh.reset();
}

void CDMesh::determine_processor_prolongation_bounding_box(const bool guessAndCheckProcPadding, const double maxCFLGuess, BoundingBox & procBbox) const
{
  procBbox.clear();
  if (stk_bulk().parallel_size() == 1 || nodes.empty())
  {
    procBbox.accommodate(stk::math::Vector3d::ZERO);
    return;
  }

  const stk::mesh::Selector cdfemParentOrActiveSelector = my_cdfem_support.get_parent_part() | get_active_part();

  for (auto && node : nodes)
  {
    if (node->needs_to_be_ale_prolonged(*this))
    {
      BoundingBox nodeBBox;
      nodeBBox.accommodate(node->coordinates());
      if (guessAndCheckProcPadding)
      {
        const double maxElemSizeForNode = compute_maximum_size_of_selected_elements_using_node(stk_bulk(), cdfemParentOrActiveSelector, node->entity());
        const double nodeBoxPadding = maxCFLGuess*maxElemSizeForNode;
        nodeBBox.pad(nodeBoxPadding);
        procBbox.accommodate(nodeBBox);
      }
    }
  }

  if (!guessAndCheckProcPadding)
  {
    procBbox.pad_epsilon();
  }
}

//--------------------------------------------------------------------------------

void
CDMesh::prolongation()
{
  stk::diag::TimeBlock timer__(my_timer_prolongation);

  const bool guessAndCheckProcPadding =
      was_mesh_previously_decomposed() &&
      stk_bulk().parallel_size() > 1 &&
      get_cdfem_displacements_field().valid();
  BoundingBox proc_target_bbox;

  double maxCFLGuess = 3.0;

  bool done = false;
  while (!done)
  {
    done = true;

    determine_processor_prolongation_bounding_box(guessAndCheckProcPadding, maxCFLGuess, proc_target_bbox);
    std::vector<BoundingBox> proc_target_bboxes;
    BoundingBox::gather_bboxes( proc_target_bbox, proc_target_bboxes );

    const size_t facet_precomm_size = my_prolong_facets.size();
    ProlongationFacet::communicate(*this, my_prolong_facets, my_prolong_node_map, proc_target_bboxes);

    build_prolongation_trees();

    communicate_prolongation_facet_fields();

    const stk::mesh::Part & active_part = get_active_part();

    // update nodal fields
    my_missing_remote_prolong_facets = false;
    for (auto && node : nodes)
    {
      if(!(node->is_prolonged()) && stk_bulk().bucket(node->entity()).member(active_part))
      {
        node->prolongate_fields(*this);
      }
    }

    if (guessAndCheckProcPadding)
    {
      const double maxCFL = compute_non_rebased_cdfem_displacement_cfl();
      const bool missingFacetsOnSomeProc = stk::is_true_on_any_proc(stk_bulk().parallel(), my_missing_remote_prolong_facets);
      const bool mustRedoGhosting = maxCFL > maxCFLGuess || missingFacetsOnSomeProc;
      if (mustRedoGhosting)
      {
        const double maxAcceptableCFLGuess = 1.e4;
        STK_ThrowRequireMsg(maxCFLGuess < maxAcceptableCFLGuess, "Error in prolongation.  Communication did not succeed with CFL guess of " << maxCFLGuess);

        const double CFLFactorOfSafety = 1.5;
        const double CFLGrowthMultiplier = 2.0;
        maxCFLGuess = std::max(CFLFactorOfSafety*maxCFL, CFLGrowthMultiplier*maxCFLGuess);
        krinolog << "Must redo ghosting for prolongation. Missing facets = " << missingFacetsOnSomeProc
            << ", current CFL estimate = " << maxCFL
            << ", New max CFL guess = " << maxCFLGuess << stk::diag::dendl;

        done = false;

        for (size_t i = facet_precomm_size; i < my_prolong_facets.size(); i++ )
          delete my_prolong_facets[i];
        my_prolong_facets.resize(facet_precomm_size);

        for (auto && node : nodes)
        {
          node->set_prolonged_flag(false);
        }
      }
    }
  }

  rebase_cdfem_displacements();

  // prolongate element fields
  for (const auto & elem : elements)
  {
    if (elem->have_subelements())
    {
      std::vector<const SubElement *> conformal_subelems;
      elem->get_subelements(conformal_subelems);

      for (auto && subelem : conformal_subelems)
      {
        subelem->prolongate_fields(*this);
      }
    }
    else
    {
      elem->prolongate_fields(*this);
    }
  }

  // We might want to check what causes any parallel discrepencies, but sync everything here
  const stk::mesh::FieldVector & all_fields = stk_bulk().mesh_meta_data().get_fields();
  const std::vector<const stk::mesh::FieldBase *> const_fields(all_fields.begin(), all_fields.end());
  for (auto && f : all_fields)
  {
    f->sync_to_host();
    f->modify_on_host();
  }
  stk::mesh::communicate_field_data(stk_bulk(), const_fields);
}

void
CDMesh::rebase_cdfem_displacements()
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::rebase_cdfem_displacements(void)"); /* %TRACE% */
  // rebase cdfem mesh displacements such that STATE_OLD is zero
  const FieldRef cdfem_displacements_field = get_cdfem_displacements_field();
  if (cdfem_displacements_field.valid())
  {
    const unsigned field_length = cdfem_displacements_field.length();
    std::vector<FieldRef> stk_fields(cdfem_displacements_field.number_of_states());
    for ( unsigned is = 0 ; is < cdfem_displacements_field.number_of_states(); ++is )
    {
      const stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(is);
      stk_fields[is] = cdfem_displacements_field.field_state(state);
    }

    std::vector< stk::mesh::Entity> objs;
    stk::mesh::get_selected_entities( stk::mesh::selectField(cdfem_displacements_field), stk_bulk().buckets( stk::topology::NODE_RANK ), objs );

    std::vector<double> old_displacement(field_length);

    const unsigned len = objs.size();
    for ( unsigned iObj(0); iObj < len; ++iObj )
    {
      stk::mesh::Entity node = objs[iObj];

      const double * old_data = field_data<double>( stk_fields[stk::mesh::StateOld], node);
      STK_ThrowRequire(nullptr != old_data);
      for (unsigned d=0; d<field_length; ++d)
      {
        old_displacement[d] = old_data[d];
      }

      for ( unsigned is = 0 ; is < cdfem_displacements_field.number_of_states(); ++is )
      {
        double * displacement = field_data<double>( stk_fields[is], node);
        STK_ThrowRequire(nullptr != displacement);
        for (unsigned d=0; d<field_length; ++d)
        {
          displacement[d] -= old_displacement[d];
        }
      }
    }
  }
}

double
CDMesh::get_maximum_cdfem_displacement() const
{ /* %TRACE[ON]% */ Trace trace__("krino::Mesh::get_maximum_cdfem_displacement(void)"); /* %TRACE% */
  double max_sqr_displacement = 0.0;
  const FieldRef cdfem_displacements_field = get_cdfem_displacements_field();
  if (cdfem_displacements_field.valid())
  {
    const auto & buckets = stk_bulk().get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(cdfem_displacements_field));

    for( auto&& b : buckets )
    {
      double * cdfem_displacements = field_data<double>(cdfem_displacements_field, *b);
      STK_ThrowAssert(nullptr != cdfem_displacements);

      const unsigned field_length = cdfem_displacements_field.length(*b);

      const int num_nodes = b->size();
      for(int n=0; n<num_nodes; ++n)
      {
        double displacement_sqrmag = 0.0;
        for (unsigned d=0; d<field_length; ++d)
        {
          const double displacements_comp = cdfem_displacements[n*field_length+d];
          displacement_sqrmag += displacements_comp*displacements_comp;
        }
        max_sqr_displacement = std::max(max_sqr_displacement, displacement_sqrmag);
      }
    }
  }

  const double local_max_sqr_displacement = max_sqr_displacement;
  stk::all_reduce_max(stk_bulk().parallel(), &local_max_sqr_displacement, &max_sqr_displacement, 1);

  return std::sqrt(max_sqr_displacement);
}

bool
CDMesh::decomposition_has_changed(const InterfaceGeometry & interfaceGeometry)
{
  stk::diag::TimeBlock timer_(my_timer_decomposition_has_changed);
  return krino::decomposition_has_changed(stk_bulk(), interfaceGeometry, my_aux_meta.active_part(), my_cdfem_support, my_phase_support);
}

void
CDMesh::print_conformal_volumes_and_surface_areas() const
{
  const stk::mesh::PartVector all_conformal_parts = my_phase_support.get_conformal_parts();
  stk::mesh::PartVector volume_conformal_parts;
  stk::mesh::PartVector side_conformal_parts;
  stk::mesh::PartVector interfacial_conformal_parts;

  for (auto && conformal_part : all_conformal_parts)
  {
    if (conformal_part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      volume_conformal_parts.push_back(conformal_part);
    }
    else if (my_phase_support.is_interface(conformal_part))
    {
      interfacial_conformal_parts.push_back(conformal_part);
    }
    else if (conformal_part->primary_entity_rank() == stk_meta().side_rank())
    {
      side_conformal_parts.push_back(conformal_part);
    }
  }

  print_volume_or_surface_area(stk_bulk(), stk::topology::ELEMENT_RANK, get_active_part(), volume_conformal_parts);
  print_volume_or_surface_area(stk_bulk(), stk_meta().side_rank(), get_active_part(), interfacial_conformal_parts);
  if ( krinolog.shouldPrint(LOG_PARTS) )
    print_volume_or_surface_area(stk_bulk(), stk_meta().side_rank(), get_active_part(), side_conformal_parts);
}

void
CDMesh::debug_output() const
{
  for (const auto & elem : elements)
    debug_elem_parts_and_relations(stk_bulk(), *elem);
  krinolog << stk::diag::dendl;

  for (auto && node : nodes)
    debug_nodal_parts_and_fields(stk_bulk(), node.get());
  krinolog << stk::diag::dendl;

  if (false)
  {
    debug_sides(stk_bulk(), get_active_part());
    krinolog << stk::diag::dendl;
  }
}

const Mesh_Element * CDMesh::find_mesh_element(stk::mesh::EntityId elemId, const std::vector<std::unique_ptr<Mesh_Element>> & searchElements)
{
  auto lb_cmp = [](const std::unique_ptr<Mesh_Element> & elem, stk::mesh::EntityId target_id) { return elem->entityId() < target_id; };
  auto first = std::lower_bound(searchElements.begin(),searchElements.end(), elemId, lb_cmp);

  if (first != searchElements.end() && (*first)->entityId() == elemId)
  {
    return (*first).get();
  }

  return nullptr;
}

} // namespace krino


