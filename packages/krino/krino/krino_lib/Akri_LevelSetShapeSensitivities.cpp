#include "Akri_Snap.hpp"
#include <Akri_LevelSetShapeSensitivities.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ChildNodeStencil.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalGradient.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_CDMesh.hpp>

namespace krino {

static std::string copied_levelset_field_name_for_snapping(const std::string & origLSFieldName)
{
  return origLSFieldName + "_COPYFORSNAPPING";
}

std::vector<LS_Field> create_levelset_copies_and_set_to_use_as_snap_fields(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
{
  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  std::vector<LS_Field> copiedLSFields;
  for(auto & lsField : origLSFields)
  {
    const FieldRef origField = lsField.isovar;
    const stk::mesh::Selector fieldSelector = stk::mesh::selectField(origField);
    const unsigned numStates = 1;
    FieldType fieldType = FieldType::REAL;
    FieldRef copyField = auxMeta.declare_field(copied_levelset_field_name_for_snapping(origField.name()), FieldType::REAL, stk::topology::NODE_RANK, numStates);
    stk::mesh::put_field_on_mesh(copyField.field(), fieldSelector, nullptr);
    copiedLSFields.emplace_back(copyField.name(), lsField.identifier, copyField, lsField.isoval, lsField.ptr, lsField.deathPtr);
    cdfemSupport.add_interpolation_field(copyField);
    cdfemSupport.add_interpolation_field(origField);
    cdfemSupport.add_do_not_snap_field(origField);
  }
  return copiedLSFields;
}

std::vector<LS_Field> get_levelset_copy_fields(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
{
  std::vector<LS_Field> copiedLSFields;
  for(auto & lsField : origLSFields)
  {
    const std::string copyFieldName = copied_levelset_field_name_for_snapping(lsField.isovar.name());
    stk::mesh::FieldBase * copyFieldBase = meta.get_field(stk::topology::NODE_RANK, copyFieldName);
    STK_ThrowRequireMsg(nullptr != copyFieldBase, "Did not find Levelset copy field " << copyFieldName);
    copiedLSFields.emplace_back(copyFieldName, lsField.identifier, FieldRef(copyFieldBase), lsField.isoval, lsField.ptr, lsField.deathPtr);
  }
  return copiedLSFields;
}

std::vector<LS_Field> update_levelset_copies_to_prepare_for_snapping(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
{
  std::vector<LS_Field> copiedLSFields = get_levelset_copy_fields(meta, origLSFields);
  for(size_t i=0; i<origLSFields.size(); ++i)
    stk::mesh::field_copy(origLSFields[i].isovar.field(), copiedLSFields[i].isovar.field());
  return copiedLSFields;
}

std::string output_sensitivity(const LevelSetShapeSensitivity & sens)
{
  std::ostringstream os;
  os << "Sensitivities for node " << sens.interfaceNodeId << ":\n";
  for (size_t i=0; i<sens.parentNodeIds.size(); ++i)
    os << "  dCoordParentLevelSet for parent node " << sens.parentNodeIds[i] << " = " << sens.dCoordsdParentLevelSets[i] << "\n";
  return os.str();
}

static stk::math::Vector3d compute_edge_vector(const unsigned dim, const FieldRef coordsField, const std::array<stk::mesh::Entity,2> & parentEdgeNodes)
{
  const stk::math::Vector3d x0(field_data<double>(coordsField, parentEdgeNodes[0]), dim);
  const stk::math::Vector3d x1(field_data<double>(coordsField, parentEdgeNodes[1]), dim);
  return x1-x0;
}

static stk::math::Vector3d compute_dXdLs_along_edge(const stk::math::Vector3d & edgeVec0to1,
  const double ls0, const double ls1)
{
  const double diff = ls0-ls1;
  const double factor = (diff == 0.) ? 0. : (1/diff);
  return edgeVec0to1 * factor;
}

static stk::math::Vector3d compute_dXdLs_along_mesh_edge(const unsigned dim,
    const FieldRef coordsField,
    const FieldRef levelSetField,
    const stk::mesh::Entity edgeNode0,
    const stk::mesh::Entity edgeNode1)
{
  const stk::math::Vector3d dx = compute_edge_vector(dim, coordsField, {{edgeNode0, edgeNode1}});
  const double ls0 = *field_data<double>(levelSetField, edgeNode0);
  const double ls1 = *field_data<double>(levelSetField, edgeNode1);
  return compute_dXdLs_along_edge(dx, ls0, ls1);
}

static void append_edge_sensitivity(const stk::mesh::BulkData & mesh,
  const stk::mesh::Entity interfaceNode,
  const stk::mesh::Entity edgeNode0,
  const stk::mesh::Entity edgeNode1,
  const stk::math::Vector3d & dXdLs,
  const double edgeWt,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const std::vector<stk::mesh::EntityId> parentNodeIds = {mesh.identifier(edgeNode0), mesh.identifier(edgeNode1)};
  const std::vector<stk::math::Vector3d> dCoordsdParentLevelSets = {(1.-edgeWt)*dXdLs, edgeWt*dXdLs};

  shapeSensitivities.emplace_back(mesh.identifier(interfaceNode), parentNodeIds, dCoordsdParentLevelSets);
}

static void append_sensitivity_for_stencil(const stk::mesh::BulkData & mesh,
  const stk::mesh::Entity childNode,
  const std::vector<stk::mesh::Entity> & parentNodes,
  const std::vector<double> & parentWeights,
  const stk::math::Vector3d & dXdLs,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;
  fill_node_ids_for_nodes(mesh, parentNodes, parentNodeIds);

  dCoordsdParentLevelSets.clear();
  dCoordsdParentLevelSets.reserve(parentWeights.size());
  for (unsigned i=0; i<parentWeights.size(); ++i)
    dCoordsdParentLevelSets.push_back(parentWeights[i] * dXdLs);

  shapeSensitivities.emplace_back(mesh.identifier(childNode), parentNodeIds, dCoordsdParentLevelSets);
}

static void append_sensitivity_for_stencil(const stk::mesh::BulkData & mesh,
  const ChildNodeStencil & stencil,
  const stk::math::Vector3d & dXdLs,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  append_sensitivity_for_stencil(mesh, stencil.childNode, stencil.parentNodes, stencil.parentWeights, dXdLs, shapeSensitivities);
}

static void append_sensitivity_for_background_node_on_interface(const stk::mesh::BulkData & mesh,
  stk::mesh::Entity interfaceNode,
  const stk::math::Vector3d & dXdLs,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const std::vector<stk::mesh::EntityId> parentNodeIds = {mesh.identifier(interfaceNode)};
  const std::vector<stk::math::Vector3d> dCoordsdParentLevelSets = {dXdLs};

  shapeSensitivities.emplace_back(mesh.identifier(interfaceNode), parentNodeIds, dCoordsdParentLevelSets);
}

static stk::math::Vector3d compute_dXdLs_along_immediate_parent_edge(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const unsigned dim,
    const FieldRef coordsField,
    const FieldRef levelSetField,
    const stk::mesh::Entity childNode)
{
  // multiple ls leading to a non-parent edge that is cut by this level set
  // Find dXdLs along this edge
  const auto nodeImmediateNodeParentsAndWeights = get_child_node_parents_and_weights(mesh,
      cdfemSupport.get_parent_node_ids_field(),
      cdfemSupport.get_parent_node_weights_field(),
      childNode);
  STK_ThrowRequire(2 == nodeImmediateNodeParentsAndWeights.size());
  return compute_dXdLs_along_mesh_edge(dim, coordsField, levelSetField, nodeImmediateNodeParentsAndWeights[0].first, nodeImmediateNodeParentsAndWeights[1].first);
}

void append_sensitivities_for_child_nodes(const stk::mesh::BulkData & mesh, const LS_Field & lsField, std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  const FieldRef coordsField = cdfemSupport.get_coords_field();
  std::vector<ChildNodeStencil> childNodeStencils;
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  fill_child_node_stencils_for_interface(mesh, cdfemSupport.get_child_node_part(), interfaceSelector, cdfemSupport.get_parent_node_ids_field(), cdfemSupport.get_parent_node_weights_field(), childNodeStencils);

  shapeSensitivities.reserve(childNodeStencils.size());
  for (auto & stencil : childNodeStencils)
  {
    STK_ThrowAssert(stencil.parentNodes.size() >= 2 && stencil.parentNodes.size() <=4);
    const stk::math::Vector3d dXdLs = (2 == stencil.parentNodes.size()) ?
        compute_dXdLs_along_mesh_edge(dim, coordsField, lsField.isovar, stencil.parentNodes[0], stencil.parentNodes[1]) :
        compute_dXdLs_along_immediate_parent_edge(mesh, cdfemSupport, dim, coordsField, lsField.isovar, stencil.childNode);
    append_sensitivity_for_stencil(mesh, stencil, dXdLs, shapeSensitivities);
  }
}

static void fill_neighbor_nodes_of_selected_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Selector & elementSelector, std::vector<stk::mesh::Entity> & neighbors)
{
  //NOTE: This is for linear simplex elements
  STK_ThrowAssert(mesh.is_automatic_aura_on());
  neighbors.clear();
  for (auto && elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    if (is_entity_selected(mesh, elementSelector, elem))
      for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
        if (elemNode != node)
          neighbors.push_back(elemNode);
  stk::util::sort_and_unique(neighbors);
}

static std::vector<stk::mesh::Entity> get_neighbor_nodes_of_selected_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Selector & elementSelector)
{
  std::vector<stk::mesh::Entity> neighbors;
  fill_neighbor_nodes_of_selected_elements(mesh, node, elementSelector, neighbors);
  return neighbors;
}

static bool first_wins_tiebreaker(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity nodeA, stk::mesh::Entity nodeB)
{
  STK_ThrowAssert(nodeA.is_local_offset_valid());
  if (!nodeB.is_local_offset_valid())
    return true;
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  return is_less_than_in_x_then_y_then_z(get_vector_field(mesh, coordsField, nodeA, dim), get_vector_field(mesh, coordsField, nodeB, dim));
}

std::pair<double, double> snapped_node_edge_position_and_square_error(const unsigned dim,
    const stk::math::Vector3d & snapDisp,
    const stk::math::Vector3d & unsnappedEdgeLoc0,
    const stk::math::Vector3d & unsnappedEdgeLoc1)
{
  const unsigned bestDim = get_edge_longest_dimension(dim, unsnappedEdgeLoc0.data(), unsnappedEdgeLoc1.data());
  const double interpPos = snapDisp[bestDim]/(unsnappedEdgeLoc1[bestDim] - unsnappedEdgeLoc0[bestDim]);
  const double pos = std::min(1., std::max(0., interpPos));
  const stk::math::Vector3d edgeLocErr = snapDisp - pos*(unsnappedEdgeLoc1-unsnappedEdgeLoc0);
  return std::make_pair(pos, edgeLocErr.length_squared());
}

std::pair<stk::mesh::Entity, double> find_cut_edge_for_snapped_node(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & neighborNodes,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const stk::mesh::Entity snappedNode)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  stk::mesh::Entity otherEdgeNode;
  double cutEdgeLoc = 0.;

  const stk::math::Vector3d nodeSnapDisp = get_vector_field(mesh, snapDisplacements, snappedNode, dim);
  if (!nodeSnapDisp.zero_length())
  {
    const stk::math::Vector3d nodeUnsnappedCoords = get_vector_field(mesh, coordsField, snappedNode, dim) - nodeSnapDisp;
    double minError = std::numeric_limits<double>::max();

    for (auto & nbrNode : neighborNodes)
    {
      const stk::math::Vector3d nbrNodeUnsnappedCoords = get_vector_field(mesh, coordsField, nbrNode, dim) - get_vector_field(mesh, snapDisplacements, nbrNode, dim);
      const auto & [edgeLoc, edgeLocErr] = snapped_node_edge_position_and_square_error(dim, nodeSnapDisp, nodeUnsnappedCoords, nbrNodeUnsnappedCoords);
      if (edgeLocErr < minError)
      {
        otherEdgeNode = nbrNode;
        cutEdgeLoc = edgeLoc;
        minError = edgeLocErr;
      }
    }
  }

  return std::make_pair(otherEdgeNode, cutEdgeLoc);
}

std::pair<stk::mesh::Entity, double> find_cut_edge_for_snapped_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const stk::mesh::Entity snappedNode)
{
  const std::vector<stk::mesh::Entity> neighborNodes = get_neighbor_nodes_of_selected_elements(mesh, snappedNode, elemSelector);
  return find_cut_edge_for_snapped_node(mesh, neighborNodes, coordsField, snapDisplacements, snappedNode);
}

stk::mesh::Entity find_best_arbitrary_neighbor_edge_node(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & neighborNodes,
    const FieldRef coordsField,
    const FieldRef isovar,
    stk::mesh::Entity node)
{
  const double nodeLs = *field_data<double>(isovar, node);

  stk::mesh::Entity bestArbitraryOtherNode;
  double maxDeltaLS = 0.;

  for (auto & nbrNode : neighborNodes)
  {
    const double nbrNodeLs = *field_data<double>(isovar, nbrNode);
    if (nbrNodeLs * nodeLs <= 0)
    {
      const double deltaLS = std::abs(nbrNodeLs-nodeLs);
      if((deltaLS > maxDeltaLS) || (deltaLS == maxDeltaLS && first_wins_tiebreaker(mesh, coordsField, nbrNode, bestArbitraryOtherNode)))
      {
        bestArbitraryOtherNode = nbrNode;
        maxDeltaLS = deltaLS;
      }
    }
  }

  STK_ThrowRequireMsg(bestArbitraryOtherNode.is_local_offset_valid(), "find_best_arbitrary_other_parent_node: couldn't find node to pick");
  return bestArbitraryOtherNode;
}

stk::mesh::Entity find_best_arbitrary_neighbor_edge_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector,
    const FieldRef coordsField,
    const FieldRef isovar,
    stk::mesh::Entity node)
{
  const std::vector<stk::mesh::Entity> neighborNodes = get_neighbor_nodes_of_selected_elements(mesh, node, elemSelector);
  return find_best_arbitrary_neighbor_edge_node(mesh, neighborNodes, coordsField, isovar, node);
}

void append_sensitivities_for_snapped_nodes_from_parents(const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & elemSelector,
  const FieldRef coordsField,
  const FieldRef snapDisplacements,
  const FieldRef lsFieldToUse,
  stk::mesh::Entity snappedNode,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const std::vector<stk::mesh::Entity> neighborNodes = get_neighbor_nodes_of_selected_elements(mesh, snappedNode, elemSelector);
  auto [otherEdgeNode, cutEdgeLoc] = find_cut_edge_for_snapped_node(mesh, neighborNodes, coordsField, snapDisplacements, snappedNode);

  if(!otherEdgeNode.is_local_offset_valid())
  {
    otherEdgeNode = find_best_arbitrary_neighbor_edge_node(mesh, neighborNodes, coordsField, lsFieldToUse, snappedNode);
    cutEdgeLoc = 0.;
  }

  const double ls0 = *field_data<double>(lsFieldToUse, snappedNode);
  const double ls1 = *field_data<double>(lsFieldToUse, otherEdgeNode);

  const stk::math::Vector3d thisNodeUnsnappedLoc = get_vector_field(mesh, coordsField, snappedNode, dim) -
    get_vector_field(mesh, snapDisplacements, snappedNode, dim);
  const stk::math::Vector3d otherNodeUnsnappedLoc = get_vector_field(mesh, coordsField, otherEdgeNode, dim) -
    get_vector_field(mesh, snapDisplacements, otherEdgeNode, dim);
  const stk::math::Vector3d dXdLs = compute_dXdLs_along_edge(otherNodeUnsnappedLoc - thisNodeUnsnappedLoc, ls0, ls1);

  append_edge_sensitivity(mesh, snappedNode, snappedNode, otherEdgeNode, dXdLs, cutEdgeLoc, shapeSensitivities);
}

void append_sensitivities_for_background_node_on_interface(const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & elemSelector,
  const FieldRef coordsField,
  const FieldRef lsFieldToUse,
  stk::mesh::Entity node,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const stk::mesh::Entity arbitraryOtherNode = find_best_arbitrary_neighbor_edge_node(mesh, elemSelector, coordsField, lsFieldToUse, node);
  const stk::math::Vector3d dXdLs = compute_dXdLs_along_mesh_edge(mesh.mesh_meta_data().spatial_dimension(), coordsField, lsFieldToUse, node, arbitraryOtherNode);
  append_edge_sensitivity(mesh, node, node, arbitraryOtherNode, dXdLs, 0., shapeSensitivities);
}

void append_sensitivities_for_background_nodes_on_interface(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  const stk::mesh::Selector interfaceNotChildSelector = (!cdfemSupport.get_child_node_part()) & interfaceSelector;
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  const stk::mesh::Selector ownershipSelector = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
  const stk::mesh::Selector cutElementSelector = phaseSupport.get_levelset_decomposed_blocks_selector(lsField.identifier);
  const FieldRef snapDisplacements = cdfemSupport.get_cdfem_snap_displacements_field();

  const bool isSnappingActive = (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE) &&
    cdfemSupport.get_cdfem_snap_displacements_field().valid();

  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, ownershipSelector & interfaceNotChildSelector ) )
  {
    for (auto node : *bucket)
    {
      if(isSnappingActive)
        append_sensitivities_for_snapped_nodes_from_parents(mesh, cutElementSelector, coordsField, snapDisplacements, lsField.isovar, node, shapeSensitivities);
      else
        append_sensitivities_for_background_node_on_interface(mesh, cutElementSelector, coordsField, lsField.isovar, node, shapeSensitivities);
    }
  }
}

static stk::math::Vector3d compute_closest_point_dXdLs_using_nodal_gradient(const stk::mesh::BulkData & mesh,
    const unsigned dim,
    const FieldRef nodalGradLevelSetField,
    const stk::mesh::Entity interfaceNode)
{
  stk::math::Vector3d dXdLs = -get_vector_field(mesh, nodalGradLevelSetField, interfaceNode, dim);
  const double sqrMag = dXdLs.length_squared();
  if (sqrMag > 0.)
    dXdLs /= sqrMag;
  return dXdLs;
}

void append_closest_point_sensitivities_for_child_nodes(const stk::mesh::BulkData & mesh, const LS_Field & lsField, std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<ChildNodeStencil> childNodeStencils;
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  fill_child_node_stencils_for_interface(mesh, cdfemSupport.get_child_node_part(), interfaceSelector, cdfemSupport.get_parent_node_ids_field(), cdfemSupport.get_parent_node_weights_field(), childNodeStencils);

  const FieldRef nodalGradLevelSetField = get_nodal_gradient_for_scalar_field(mesh.mesh_meta_data(), lsField.isovar);

  shapeSensitivities.reserve(childNodeStencils.size());
  for (auto & stencil : childNodeStencils)
  {
    const stk::math::Vector3d dXdLs = compute_closest_point_dXdLs_using_nodal_gradient(mesh, dim, nodalGradLevelSetField, stencil.childNode);
    append_sensitivity_for_stencil(mesh, stencil, dXdLs, shapeSensitivities);
  }
}


void append_closest_point_sensitivities_for_background_nodes_on_interface(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  const stk::mesh::Selector interfaceNotChildSelector = (!cdfemSupport.get_child_node_part()) & interfaceSelector;
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  const stk::mesh::Selector ownershipSelector = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
  const stk::mesh::Selector cutElementSelector = phaseSupport.get_levelset_decomposed_blocks_selector(lsField.identifier);
  const FieldRef snapDisplacements = cdfemSupport.get_cdfem_snap_displacements_field();
  const FieldRef nodalGradLevelSetField = get_nodal_gradient_for_scalar_field(mesh.mesh_meta_data(), lsField.isovar);
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> parentNodes;
  std::vector<double> parentNodeWts;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;

  const bool isSnappingActive = (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE) &&
    cdfemSupport.get_cdfem_snap_displacements_field().valid();

  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, ownershipSelector & interfaceNotChildSelector ) )
  {
    for (auto node : *bucket)
    {
      const stk::math::Vector3d dXdLs = compute_closest_point_dXdLs_using_nodal_gradient(mesh, dim, nodalGradLevelSetField, node);
      if(isSnappingActive)
      {
        const auto & [otherEdgeNode, cutEdgeLoc] = find_cut_edge_for_snapped_node(mesh, cutElementSelector, coordsField, snapDisplacements, node);
        if(otherEdgeNode.is_local_offset_valid())
          append_edge_sensitivity(mesh, node, node, otherEdgeNode, dXdLs, cutEdgeLoc, shapeSensitivities);
        else
          append_sensitivity_for_background_node_on_interface(mesh, node, dXdLs, shapeSensitivities);
      }
      else
      {
        append_sensitivity_for_background_node_on_interface(mesh, node, dXdLs, shapeSensitivities);
      }
    }
  }
}

void fill_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities,
    const bool doComputeClosestPointSensitivities)
{
  shapeSensitivities.clear();
  if (doComputeClosestPointSensitivities)
  {
    append_closest_point_sensitivities_for_child_nodes(mesh, lsField, shapeSensitivities);
    append_closest_point_sensitivities_for_background_nodes_on_interface(mesh, lsField, shapeSensitivities);
  }
  else
  {
    append_sensitivities_for_child_nodes(mesh, lsField, shapeSensitivities);
    append_sensitivities_for_background_nodes_on_interface(mesh, lsField, shapeSensitivities);
  }
}

std::vector<LevelSetShapeSensitivity> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const LS_Field & lsField, const bool doComputeClosestPointSensitivities)
{
  std::vector<LevelSetShapeSensitivity> shapeSensitivities;

  fill_levelset_shape_sensitivities(mesh, lsField, shapeSensitivities, doComputeClosestPointSensitivities);

  return shapeSensitivities;
}

std::vector<std::vector<LevelSetShapeSensitivity>> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & lsFields, const bool doComputeClosestPointSensitivities)
{
  std::vector<std::vector<LevelSetShapeSensitivity>> shapeSensitivities(lsFields.size());
  for (unsigned i=0; i<lsFields.size(); ++i)
  {
    fill_levelset_shape_sensitivities(mesh, lsFields[i], shapeSensitivities[i], doComputeClosestPointSensitivities);
  }

  return shapeSensitivities;
}

std::vector<stk::mesh::Entity> get_owned_interface_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & interfaceSelector)
{
  std::vector<stk::mesh::Entity> interfaceSides;

  for ( auto * bucket : mesh.get_buckets( mesh.mesh_meta_data().side_rank(), mesh.mesh_meta_data().locally_owned_part() & interfaceSelector ) )
    interfaceSides.insert(interfaceSides.end(), bucket->begin(), bucket->end());

  return interfaceSides;
}

size_t get_sensitivity_index(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
{
  const stk::mesh::EntityId nodeId = mesh.identifier(node);
  const auto iter = std::lower_bound(sortedShapeSensitivities.begin(), sortedShapeSensitivities.end(), nodeId,
    [](const LevelSetShapeSensitivity & sens, const stk::mesh::EntityId searchNodeId) {return sens.interfaceNodeId < searchNodeId;});
  STK_ThrowRequire(iter != sortedShapeSensitivities.end() && iter->interfaceNodeId == nodeId);
  return std::distance(sortedShapeSensitivities.begin(), iter);
}

template <size_t NNODES>
std::array<size_t,NNODES> get_facet_sensitivity_indices(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & negativeSideElementSelector,
    const stk::mesh::Entity side,
    const std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);

template<>
std::array<size_t,2> get_facet_sensitivity_indices<2>(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & negativeSideElementSelector,
    const stk::mesh::Entity side,
    const std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
{
  const std::array<stk::mesh::Entity,2> sideNodes = get_oriented_line_side_nodes(mesh, negativeSideElementSelector, side);
  std::array<size_t,2> facetSensIndices =
    {
      get_sensitivity_index(mesh, sideNodes[0], sortedShapeSensitivities),
      get_sensitivity_index(mesh, sideNodes[1], sortedShapeSensitivities)
    };
  return facetSensIndices;
}

template<>
std::array<size_t,3> get_facet_sensitivity_indices<3>(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & negativeSideElementSelector,
    const stk::mesh::Entity side,
    const std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
{
  const std::array<stk::mesh::Entity,3> sideNodes = get_oriented_triangle_side_nodes(mesh, negativeSideElementSelector, side);
  std::array<size_t,3> facetSensIndices =
    {
      get_sensitivity_index(mesh, sideNodes[0], sortedShapeSensitivities),
      get_sensitivity_index(mesh, sideNodes[1], sortedShapeSensitivities),
      get_sensitivity_index(mesh, sideNodes[2], sortedShapeSensitivities)
    };
  return facetSensIndices;
}

template <size_t NNODES>
void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
     const LS_Field & lsField,
     const bool doComputeClosestPointSensitivities,
     std::vector<std::array<size_t,NNODES>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
{
  fill_levelset_shape_sensitivities(mesh, lsField, sortedShapeSensitivities, doComputeClosestPointSensitivities);

  std::sort(sortedShapeSensitivities.begin(), sortedShapeSensitivities.end(),
      [](const LevelSetShapeSensitivity & a, const LevelSetShapeSensitivity & b) {return a.interfaceNodeId < b.interfaceNodeId;});

  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  const stk::mesh::Selector negativeSideElementSelector = phaseSupport.get_negative_levelset_block_selector(lsField.identifier);

  const std::vector<stk::mesh::Entity> interfaceSides = get_owned_interface_sides(mesh, interfaceSelector);

  facetSensitivityIndices.clear();
  facetSensitivityIndices.reserve(interfaceSides.size());

  for (auto & interfaceSide : interfaceSides)
  {
    const std::array<size_t,NNODES> facetSensIndices = get_facet_sensitivity_indices<NNODES>(mesh, negativeSideElementSelector, interfaceSide, sortedShapeSensitivities);
    facetSensitivityIndices.push_back(facetSensIndices);
  }
}

stk::mesh::Entity get_facet_node(const stk::mesh::BulkData & mesh, const LevelSetShapeSensitivity &sens)
{
  stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId);
  STK_ThrowAssert(mesh.is_valid(node));
  return node;
}

template <size_t NNODES>
std::array<const LevelSetShapeSensitivity*,NNODES> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,NNODES> & facetSensIndices)
{
  std::array<const LevelSetShapeSensitivity*,NNODES> facetNodeSens;
  for (size_t i=0; i<NNODES; ++i)
    facetNodeSens[i] = &sensitivities[facetSensIndices[i]];
  return facetNodeSens;
}

template <size_t NNODES>
std::array<stk::mesh::Entity,NNODES> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,NNODES> & facetNodeSens)
{
  std::array<stk::mesh::Entity,NNODES> facetNodes;
  for (size_t i=0; i<NNODES; ++i)
    facetNodes[i] = get_facet_node(mesh, *facetNodeSens[i]);
  return facetNodes;
}

void setup_fields_for_conforming_decomposition(stk::mesh::MetaData & meta, const std::vector<krino::LS_Field> & lsFields, const bool doSetupSnapping, const bool doSetupNodalLevelsetGradient)
{
  krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meta);
  const krino::FieldRef coordsField = meta.coordinate_field();

  if (doSetupSnapping)
  {
    cdfemSupport.register_cdfem_snap_displacements_field();
    if (!doSetupNodalLevelsetGradient)
      krino::create_levelset_copies_and_set_to_use_as_snap_fields(meta, lsFields);
  }

  if (doSetupNodalLevelsetGradient)
  {
    for(auto & lsField : lsFields)
    {
      krino::FieldRef nodalGrad = krino::register_nodal_gradient_for_scalar_field(meta, lsField.isovar);
      cdfemSupport.add_interpolation_field(nodalGrad);
    }
  }

  cdfemSupport.set_coords_field(coordsField);
  cdfemSupport.register_parent_node_ids_field();
  cdfemSupport.finalize_fields();
}

void decompose_mesh_to_conform_to_levelsets(stk::mesh::BulkData & mesh, const std::vector<krino::LS_Field> & lsFields, const bool doSnapping, const bool doSetupNodalLevelsetGradient)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  krino::AuxMetaData & auxMeta = krino::AuxMetaData::get(meta);
  krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meta);
  krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);
  if(doSnapping) 
  {
    cdfemSupport.set_cdfem_edge_degeneracy_handling(krino::Edge_Degeneracy_Handling::SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  }

  std::vector<krino::LS_Field> lsFieldsToUse = lsFields;
  if (doSetupNodalLevelsetGradient)
  {
    for(auto & lsField : lsFields)
      krino::update_nodal_gradient(mesh, cdfemSupport.get_coords_field(), lsField.isovar);
  }
  else if(doSnapping)
  {
    lsFieldsToUse = krino::update_levelset_copies_to_prepare_for_snapping(meta, lsFields);
  }
  std::unique_ptr<krino::InterfaceGeometry> interfaceGeometry = krino::create_levelset_geometry(meta.spatial_dimension(), auxMeta.active_part(), cdfemSupport, phaseSupport, lsFieldsToUse);
  
  auxMeta.clear_force_64bit_flag();
  krino::CDMesh::decompose_mesh(mesh, *interfaceGeometry);
}

// Explicit template instantiation

template void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
     const LS_Field & lsField,
     const bool doComputeClosestPointSensitivities,
     std::vector<std::array<size_t,2>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);
template void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
     const LS_Field & lsField,
     const bool doComputeClosestPointSensitivities,
     std::vector<std::array<size_t,3>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);

template std::array<const LevelSetShapeSensitivity*,2> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,2> & facetSensIndices);
template std::array<const LevelSetShapeSensitivity*,3> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,3> & facetSensIndices);

template std::array<stk::mesh::Entity,2> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,2> & facetNodeSens);
template std::array<stk::mesh::Entity,3> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,3> & facetNodeSens);

}

