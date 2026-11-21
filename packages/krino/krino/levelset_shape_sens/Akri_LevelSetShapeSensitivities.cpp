#include <Akri_Snap.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_ChildNodeStencil.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceInverseGradientArea.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalMotionConstraint.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <Akri_SharpFeature.hpp>
#include <Akri_SurfaceNodalNormal.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_CDMesh.hpp>

#include <algorithm>
#include <cmath>
#include <memory>

namespace krino {

static std::string copied_levelset_field_name_for_snapping(const std::string & origLSFieldName)
{
  return origLSFieldName + "_COPYFORSNAPPING";
}

static std::vector<LS_Field> create_levelset_copies_and_set_to_use_as_snap_fields(
    const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
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

static std::vector<LS_Field> get_levelset_copy_fields(
    const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
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

static std::vector<LS_Field> update_levelset_copies_to_prepare_for_snapping(
    const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields)
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
  auto coordsData = coordsField.data<double>();
  const stk::math::Vector3d x0(coordsData.entity_values(parentEdgeNodes[0]), dim);
  const stk::math::Vector3d x1(coordsData.entity_values(parentEdgeNodes[1]), dim);
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
  auto levelSetData = levelSetField.data<double>();
  const double ls0 = levelSetData.entity_values(edgeNode0)();
  const double ls1 = levelSetData.entity_values(edgeNode1)();
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

static stk::math::Vector3d compute_dXdLs_along_immediate_parent_edge(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const unsigned dim,
    const FieldRef coordsField,
    const FieldRef levelSetField,
    const stk::mesh::Entity childNode)
{
  // With multiple level sets the node can lie on an edge between nodes created by earlier cuts.
  const auto nodeImmediateNodeParentsAndWeights = get_child_node_parents_and_weights(mesh,
      cdfemSupport.get_parent_node_ids_field(),
      cdfemSupport.get_parent_node_weights_field(),
      childNode);
  STK_ThrowRequire(2 == nodeImmediateNodeParentsAndWeights.size());
  return compute_dXdLs_along_mesh_edge(dim, coordsField, levelSetField, nodeImmediateNodeParentsAndWeights[0].first, nodeImmediateNodeParentsAndWeights[1].first);
}

static void append_sensitivities_for_child_nodes(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSelector,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  const FieldRef coordsField = cdfemSupport.get_coords_field();
  std::vector<ChildNodeStencil> childNodeStencils;
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

static std::pair<stk::mesh::Entity, double> find_cut_edge_for_snapped_node(const stk::mesh::BulkData & mesh,
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

static stk::mesh::Entity find_best_arbitrary_neighbor_edge_node(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & neighborNodes,
    const FieldRef coordsField,
    const FieldRef isovar,
    stk::mesh::Entity node)
{
  auto isovarData = isovar.data<double>();
  const double nodeLs = isovarData.entity_values(node)();

  stk::mesh::Entity bestArbitraryOtherNode;
  double maxDeltaLS = 0.;

  for (auto & nbrNode : neighborNodes)
  {
    const double nbrNodeLs = isovarData.entity_values(nbrNode)();
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

  STK_ThrowRequireMsg(bestArbitraryOtherNode.is_local_offset_valid(),
      "find_best_arbitrary_neighbor_edge_node: no neighbor with opposite or zero level set");
  return bestArbitraryOtherNode;
}

static std::string surface_nodal_normal_field_name(const std::string & lsFieldName)
{
  return "SURFNORMAL_" + lsFieldName;
}

static FieldRef update_interface_surface_nodal_normal(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSelector)
{
  // On the current (for snapping, snapped) coordinates: the interface the objective is evaluated
  // on. Each side is counted once, from its negative-phase element.
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const FieldRef normalField = mesh.mesh_meta_data().get_field(
      stk::topology::NODE_RANK, surface_nodal_normal_field_name(lsField.isovar.name()));
  STK_ThrowRequireMsg(normalField.valid(),
      "The interface normal of "
          << lsField.isovar.name()
          << " was not registered; construct a ConformingDecompositionWithSensitivities before "
             "committing the mesh.");
  update_surface_nodal_normal(mesh,
      cdfemSupport.get_coords_field(),
      normalField,
      interfaceSelector,
      phaseSupport.get_negative_levelset_block_selector(lsField.identifier) &
          AuxMetaData::get(mesh.mesh_meta_data()).active_part());
  return normalField;
}

// The interface fields for the velocity of background nodes on the interface with no cut edge.
struct InterfaceNormalAndInverseGradientArea
{
  FieldRef normal;
  FieldRef inverseGradientArea;
};

// The velocity of a background node on the interface with no cut edge: along the surface nodal
// normal N, at the speed that makes the enclosed volume's derivative exact.  False if N is zero.
static bool compute_dXdLs_along_surface_normal(const stk::mesh::BulkData & mesh,
    const InterfaceNormalAndInverseGradientArea & interfaceFields,
    const stk::mesh::Entity node,
    stk::math::Vector3d & dXdLs)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::math::Vector3d normal = get_vector_field(mesh, interfaceFields.normal, node, dim);
  const double sqrMag = normal.length_squared();
  if (sqrMag == 0.) return false;
  const double inverseGradientArea =
      interfaceFields.inverseGradientArea.data<double>().entity_values(node)();
  dXdLs = -(inverseGradientArea / sqrMag) * normal;
  return true;
}

// Whether snapping would let the node snap toward the neighbor, along their edge: a node on a
// domain boundary only along that boundary (side part compatibility), and a node on a sharp feature
// only along it.
static bool can_snap_toward_neighbor(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const stk::mesh::Entity node,
    const stk::mesh::Entity nbr)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  if (!which_intersection_point_nodes_are_compatible_for_snapping(
          mesh, auxMeta, phaseSupport, {node, nbr})[0])
    return false;
  return !sharpFeatureInfo ||
      is_intersection_point_node_compatible_for_snapping_based_on_sharp_features(
          *sharpFeatureInfo, node, {node, nbr});
}

// The edge a node on the interface (zero snap displacement) moves along when it is snapped with
// predictable snap selection and the level set increases: among the edges to neighbors with
// negative level set (the ones cut once the node's level set is positive) that snapping would allow
// the node to snap along, the one along which the intersection point moves slowest, with ties
// broken by direction (x, then y, then z).  Like snapping, it allows a node on a domain boundary
// only edges on that boundary, and a node on a sharp feature only edges along it.  Uses the
// unsnapped coordinates and the level set that the snapping used.
static stk::mesh::Entity
find_predictable_snap_edge_for_node_on_interface(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & neighborNodes,
    const SharpFeatureInfo * sharpFeatureInfo,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const FieldRef isovar,
    const stk::mesh::Entity node)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::mesh::Part & childNodePart =
      CDFEM_Support::get(mesh.mesh_meta_data()).get_child_node_part();
  auto isovarData = isovar.data<double>();
  const double nodeLs = isovarData.entity_values(node)();
  const stk::math::Vector3d nodeLoc = get_vector_field(mesh, coordsField, node, dim) -
      get_vector_field(mesh, snapDisplacements, node, dim);

  stk::mesh::Entity bestNode;
  float bestSpeed = 0.f;
  stk::math::Vector3d bestDirection;
  for (auto nbr : neighborNodes)
  {
    if (mesh.bucket(nbr).member(childNodePart))
      continue; // snapping only moves nodes along background edges
    const double nbrLs = isovarData.entity_values(nbr)();
    if (nbrLs >= 0. || nbrLs == nodeLs) continue;
    if (!can_snap_toward_neighbor(mesh, sharpFeatureInfo, node, nbr)) continue;
    const stk::math::Vector3d nbrLoc = get_vector_field(mesh, coordsField, nbr, dim) -
        get_vector_field(mesh, snapDisplacements, nbr, dim);
    const stk::math::Vector3d velocity = (nbrLoc - nodeLoc) / (nodeLs - nbrLs);
    const float speed = static_cast<float>(velocity.length());
    const stk::math::Vector3d direction = velocity.unit_vector();
    if (!bestNode.is_local_offset_valid() || speed < bestSpeed ||
        (speed == bestSpeed && is_less_than_in_x_then_y_then_z(direction, bestDirection)))
    {
      bestNode = nbr;
      bestSpeed = speed;
      bestDirection = direction;
    }
  }
  return bestNode;
}

static void append_sensitivities_for_snapped_nodes_from_parents(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector,
    const SharpFeatureInfo * sharpFeatureInfo,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const FieldRef lsFieldToUse,
    const InterfaceNormalAndInverseGradientArea & interfaceFields,
    stk::mesh::Entity snappedNode,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const std::vector<stk::mesh::Entity> neighborNodes = get_neighbor_nodes_of_selected_elements(mesh, snappedNode, elemSelector);
  auto [otherEdgeNode, cutEdgeLoc] = find_cut_edge_for_snapped_node(mesh, neighborNodes, coordsField, snapDisplacements, snappedNode);

  const bool isPredictable =
      CDFEM_Support::get(mesh.mesh_meta_data()).get_snap_selection() == Snap_Selection::PREDICTABLE;
  if (!otherEdgeNode.is_local_offset_valid() && isPredictable)
  {
    // A node on the interface: snapping moves it along a single edge as the level set increases.
    otherEdgeNode = find_predictable_snap_edge_for_node_on_interface(mesh,
        neighborNodes,
        sharpFeatureInfo,
        coordsField,
        snapDisplacements,
        lsFieldToUse,
        snappedNode);
    cutEdgeLoc = 0.;
  }

  if(!otherEdgeNode.is_local_offset_valid())
  {
    otherEdgeNode = find_best_arbitrary_neighbor_edge_node(mesh, neighborNodes, coordsField, lsFieldToUse, snappedNode);
    stk::math::Vector3d dXdLs;
    if (compute_dXdLs_along_surface_normal(mesh, interfaceFields, snappedNode, dXdLs))
    {
      append_edge_sensitivity(
          mesh, snappedNode, snappedNode, otherEdgeNode, dXdLs, 0., shapeSensitivities);
      return;
    }
    cutEdgeLoc = 0.;
  }

  auto lsFieldToUseData = lsFieldToUse.data<double>();
  const double ls0 = lsFieldToUseData.entity_values(snappedNode)();
  const double ls1 = lsFieldToUseData.entity_values(otherEdgeNode)();

  const stk::math::Vector3d thisNodeUnsnappedLoc = get_vector_field(mesh, coordsField, snappedNode, dim) -
    get_vector_field(mesh, snapDisplacements, snappedNode, dim);
  const stk::math::Vector3d otherNodeUnsnappedLoc = get_vector_field(mesh, coordsField, otherEdgeNode, dim) -
    get_vector_field(mesh, snapDisplacements, otherEdgeNode, dim);
  const stk::math::Vector3d dXdLs = compute_dXdLs_along_edge(otherNodeUnsnappedLoc - thisNodeUnsnappedLoc, ls0, ls1);

  append_edge_sensitivity(mesh, snappedNode, snappedNode, otherEdgeNode, dXdLs, cutEdgeLoc, shapeSensitivities);
}

// The parts of a side that snapping keeps nodes on (the exposed boundary, io and nonconformal side
// parts), with the phase parts of a decomposed part replaced by its nonconformal part, as snapping
// does.
static std::vector<const stk::mesh::Part *>
get_side_parts_to_check(const stk::mesh::BulkData & mesh, const stk::mesh::Entity side)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  std::vector<const stk::mesh::Part *> parts;
  for (auto * part : mesh.bucket(side).supersets())
    if (is_part_to_check_for_snapping_compatibility(
            phaseSupport, auxMeta, mesh.mesh_meta_data().side_rank(), *part))
      parts.push_back(
          phaseSupport.is_decomposed(*part) ? &phaseSupport.get_nonconforming_part(*part) : part);
  stk::util::sort_and_unique(parts);
  return parts;
}

static stk::mesh::Entity get_side_node_shared_with_other_side(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1)
{
  for (auto n0 : StkMeshEntities{mesh.begin_nodes(side0), mesh.end_nodes(side0)})
    if (n0 != node)
      for (auto n1 : StkMeshEntities{mesh.begin_nodes(side1), mesh.end_nodes(side1)})
        if (n0 == n1) return n0;
  return stk::mesh::Entity();
}

// Two boundary sides at the node are on the same wall unless their side parts differ or, with the
// sharp feature angle set, they meet at a sharp edge (2D: at a sharp angle).
static bool are_boundary_sides_on_same_wall(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const double cosSharpFeatureAngle,
    const stk::mesh::Entity node,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1)
{
  if (get_side_parts_to_check(mesh, side0) != get_side_parts_to_check(mesh, side1)) return false;
  if (cosSharpFeatureAngle >= 1.) return true;
  if (mesh.mesh_meta_data().spatial_dimension() == 2)
    return !is_sharp_between_line_sides(
        mesh, coordsField, node, side0, side1, cosSharpFeatureAngle);
  const stk::mesh::Entity sharedNode =
      get_side_node_shared_with_other_side(mesh, node, side0, side1);
  return !is_sharp_between_tri_sides(
      mesh, coordsField, {{node, sharedNode}}, side0, side1, cosSharpFeatureAngle);
}

static void append_sensitivities_for_background_node_on_interface(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector,
    const double cosSharpFeatureAngle,
    const stk::mesh::Selector & boundarySideSelector,
    const stk::mesh::Selector & parentElementSelector,
    const FieldRef coordsField,
    const FieldRef lsFieldToUse,
    const InterfaceNormalAndInverseGradientArea & interfaceFields,
    stk::mesh::Entity node,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const std::vector<stk::mesh::Entity> neighborNodes =
      get_neighbor_nodes_of_selected_elements(mesh, node, elemSelector);
  // The other edge node only names the second parent; its level set does not affect the node.
  const stk::mesh::Entity arbitraryOtherNode =
      find_best_arbitrary_neighbor_edge_node(mesh, neighborNodes, coordsField, lsFieldToUse, node);
  stk::math::Vector3d dXdLs;
  if (compute_dXdLs_along_surface_normal(mesh, interfaceFields, node, dXdLs))
  {
    // On the domain boundary, the node follows the interface's contact with the boundary, with the
    // same normal speed so the volume stays exact.
    const NodalMotionConstraint constraint = compute_zero_levelset_boundary_node_constraint(mesh,
        coordsField,
        lsFieldToUse,
        node,
        boundarySideSelector,
        parentElementSelector,
        [&](const stk::mesh::Entity side0, const stk::mesh::Entity side1)
        {
          return are_boundary_sides_on_same_wall(
              mesh, coordsField, cosSharpFeatureAngle, node, side0, side1);
        });
    dXdLs = constrain_keeping_normal_speed(constraint, dXdLs);
  }
  else
  {
    dXdLs = compute_dXdLs_along_mesh_edge(mesh.mesh_meta_data().spatial_dimension(),
        coordsField,
        lsFieldToUse,
        node,
        arbitraryOtherNode);
  }
  append_edge_sensitivity(mesh, node, node, arbitraryOtherNode, dXdLs, 0., shapeSensitivities);
}

static bool is_snapping_active(const CDFEM_Support & cdfemSupport)
{
  return cdfemSupport.get_cdfem_edge_degeneracy_handling() ==
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE &&
      cdfemSupport.get_cdfem_snap_displacements_field().valid();
}

static void append_sensitivities_for_background_nodes_on_interface(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSelector,
    const InterfaceNormalAndInverseGradientArea & interfaceFields,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interfaceNotChildSelector = (!cdfemSupport.get_child_node_part()) & interfaceSelector;
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  const stk::mesh::Selector ownershipSelector = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
  const stk::mesh::Selector cutElementSelector = phaseSupport.get_levelset_decomposed_blocks_selector(lsField.identifier);
  const FieldRef snapDisplacements = cdfemSupport.get_cdfem_snap_displacements_field();

  const bool isSnappingActive = is_snapping_active(cdfemSupport);

  // The boundary and sharp features snapping uses: on the undecomposed (parent) elements, at the
  // unsnapped coordinates.
  const stk::mesh::Selector parentElementSelector = get_potential_cdfem_parent_element_selector(
      AuxMetaData::get(mesh.mesh_meta_data()).active_part(), cdfemSupport);
  const stk::mesh::Selector boundarySideSelector = get_snapping_compatibility_side_selector(mesh);
  std::unique_ptr<SharpFeatureInfo> sharpFeatureInfo;
  const double sharpFeatureAngle = cdfemSupport.get_snapping_sharp_feature_angle_in_degrees();
  const double cosSharpFeatureAngle =
      (sharpFeatureAngle > 0.) ? std::cos(sharpFeatureAngle * M_PI / 180.) : 1.;
  if (isSnappingActive && sharpFeatureAngle > 0.)
  {
    sharpFeatureInfo = std::make_unique<SharpFeatureInfo>();
    const bool includeSharedNodes = true;
    sharpFeatureInfo->find_sharp_features(mesh,
        coordsField,
        parentElementSelector,
        boundarySideSelector,
        cosSharpFeatureAngle,
        includeSharedNodes,
        snapDisplacements);
  }

  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, ownershipSelector & interfaceNotChildSelector ) )
  {
    for (auto node : *bucket)
    {
      if(isSnappingActive)
        append_sensitivities_for_snapped_nodes_from_parents(mesh,
            cutElementSelector,
            sharpFeatureInfo.get(),
            coordsField,
            snapDisplacements,
            lsField.isovar,
            interfaceFields,
            node,
            shapeSensitivities);
      else
        append_sensitivities_for_background_node_on_interface(mesh,
            cutElementSelector,
            cosSharpFeatureAngle,
            boundarySideSelector,
            parentElementSelector,
            coordsField,
            lsField.isovar,
            interfaceFields,
            node,
            shapeSensitivities);
    }
  }
}

static void project_sensitivities_onto_nodal_normal(const stk::mesh::BulkData & mesh,
    const FieldRef normalField,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  for (auto & sens : shapeSensitivities)
  {
    const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId);
    const stk::math::Vector3d normal = get_vector_field(mesh, normalField, node, dim);
    const double sqrMag = normal.length_squared();
    if (sqrMag == 0.) continue;
    for (auto & dCoords : sens.dCoordsdParentLevelSets)
      dCoords = (Dot(dCoords, normal) / sqrMag) * normal;
  }
}

static void fill_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSelector,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities,
    const ShapeSensitivityMethod method)
{
  shapeSensitivities.clear();
  const InterfaceNormalAndInverseGradientArea interfaceFields{
      update_interface_surface_nodal_normal(mesh, lsField, interfaceSelector),
      update_interface_inverse_gradient_area(mesh, lsField, interfaceSelector)};
  append_sensitivities_for_child_nodes(mesh, lsField, interfaceSelector, shapeSensitivities);
  append_sensitivities_for_background_nodes_on_interface(
      mesh, lsField, interfaceSelector, interfaceFields, shapeSensitivities);
  if (method == ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT)
    project_sensitivities_onto_nodal_normal(mesh, interfaceFields.normal, shapeSensitivities);
}

static std::vector<stk::mesh::Entity> get_owned_interface_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & interfaceSelector)
{
  std::vector<stk::mesh::Entity> interfaceSides;

  // Inactive sides (of inactive parent elements) can keep their interface parts.
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  for (auto * bucket : mesh.get_buckets(meta.side_rank(),
           meta.locally_owned_part() & AuxMetaData::get(meta).active_part() & interfaceSelector))
    interfaceSides.insert(interfaceSides.end(), bucket->begin(), bucket->end());

  return interfaceSides;
}

static size_t get_sensitivity_index(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
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

static void sort_by_interface_node_id(std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  std::sort(shapeSensitivities.begin(),
      shapeSensitivities.end(),
      [](const LevelSetShapeSensitivity & a, const LevelSetShapeSensitivity & b)
      { return a.interfaceNodeId < b.interfaceNodeId; });
}

static stk::mesh::Entity
get_facet_node(const stk::mesh::BulkData & mesh, const LevelSetShapeSensitivity & sens)
{
  stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId);
  STK_ThrowAssert(mesh.is_valid(node));
  return node;
}

template <size_t NNODES>
std::array<const LevelSetShapeSensitivity *, NNODES>
InterfaceFacetsWithSensitivities<NNODES>::node_sensitivities(const size_t facet) const
{
  std::array<const LevelSetShapeSensitivity *, NNODES> facetNodeSens;
  for (size_t i=0; i<NNODES; ++i)
    facetNodeSens[i] = &sensitivities[facetSensitivityIndices[facet][i]];
  return facetNodeSens;
}

template <size_t NNODES>
std::array<stk::mesh::Entity, NNODES> InterfaceFacetsWithSensitivities<NNODES>::nodes(
    const stk::mesh::BulkData & mesh, const size_t facet) const
{
  std::array<stk::mesh::Entity, NNODES> facetNodes;
  for (size_t i=0; i<NNODES; ++i)
    facetNodes[i] = get_facet_node(mesh, sensitivities[facetSensitivityIndices[facet][i]]);
  return facetNodes;
}

ConformingDecompositionWithSensitivities::ConformingDecompositionWithSensitivities(
    stk::mesh::MetaData & meta,
    const std::vector<LS_Field> & lsFields,
    const bool doSnapping,
    const stk::mesh::Selector * interfaceRestriction)
    : myMeta(meta), myLSFields(lsFields), myDoSnapping(doSnapping)
{
  if (interfaceRestriction != nullptr) myInterfaceRestriction = *interfaceRestriction;

  CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  const FieldRef coordsField = meta.coordinate_field();

  if (doSnapping)
  {
    // The sensitivities read the original level set after a snapped decomposition, so the
    // decomposition snaps a copy; predictable snap selection makes the snapped sensitivities exact.
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.set_snap_selection(Snap_Selection::PREDICTABLE);
    create_levelset_copies_and_set_to_use_as_snap_fields(meta, lsFields);
  }
  else
  {
    // A cut within the edge tolerance of a node is moved onto the node, which then does not move
    // for small level set perturbations, so use the smallest tolerance the parser allows rather
    // than the default 1e-3.  Callers can still set another tolerance after this.
    constexpr double edgeTolWithoutSnapping = 1.e-12;
    cdfemSupport.set_cdfem_edge_tol(edgeTolWithoutSnapping);
  }

  const Phase_Support & phaseSupport = Phase_Support::get(meta);
  for (auto & lsField : lsFields)
  {
    const stk::mesh::Selector interfaceSelector =
        phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
    register_surface_nodal_normal_field(
        meta, surface_nodal_normal_field_name(lsField.isovar.name()), interfaceSelector);
    register_interface_inverse_gradient_area_field(meta, lsField);
  }

  cdfemSupport.set_coords_field(coordsField);
  cdfemSupport.register_parent_node_ids_field();
  cdfemSupport.finalize_fields();
}

void ConformingDecompositionWithSensitivities::decompose()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  AuxMetaData & auxMeta = AuxMetaData::get(myMeta);
  CDFEM_Support & cdfemSupport = CDFEM_Support::get(myMeta);
  Phase_Support & phaseSupport = Phase_Support::get(myMeta);
  if (myDoSnapping)
  {
    cdfemSupport.set_cdfem_edge_degeneracy_handling(
        Edge_Degeneracy_Handling::SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  }

  const std::vector<LS_Field> lsFieldsToUse = myDoSnapping ?
      update_levelset_copies_to_prepare_for_snapping(myMeta, myLSFields) :
      myLSFields;
  std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_levelset_geometry(
      myMeta.spatial_dimension(), auxMeta.active_part(), cdfemSupport, phaseSupport, lsFieldsToUse);

  auxMeta.clear_force_64bit_flag();
  CDMesh::decompose_mesh(mesh, *interfaceGeometry);
}

void ConformingDecompositionWithSensitivities::reset_to_undecomposed_mesh()
{
  CDMesh::reset_mesh_to_original_undecomposed_state(myMeta.mesh_bulk_data());
}

stk::mesh::Selector
ConformingDecompositionWithSensitivities::interface_selector(const LS_Field & lsField) const
{
  const stk::mesh::Selector interfaceSelector =
      Phase_Support::get(myMeta).get_negative_levelset_interface_selector(lsField.identifier);
  return myInterfaceRestriction ? (interfaceSelector & *myInterfaceRestriction) : interfaceSelector;
}

const LS_Field &
ConformingDecompositionWithSensitivities::levelset_field(const unsigned lsIndex) const
{
  STK_ThrowRequireMsg(lsIndex < myLSFields.size(),
      "Level set index " << lsIndex << " is out of range (" << myLSFields.size()
                         << " level sets).");
  return myLSFields[lsIndex];
}

std::vector<LevelSetShapeSensitivity> ConformingDecompositionWithSensitivities::sensitivities(
    const unsigned lsIndex, const ShapeSensitivityMethod method) const
{
  std::vector<LevelSetShapeSensitivity> shapeSensitivities;
  const LS_Field & lsField = levelset_field(lsIndex);
  fill_levelset_shape_sensitivities(
      myMeta.mesh_bulk_data(), lsField, interface_selector(lsField), shapeSensitivities, method);
  return shapeSensitivities;
}

std::vector<std::vector<LevelSetShapeSensitivity>>
ConformingDecompositionWithSensitivities::all_sensitivities(
    const ShapeSensitivityMethod method) const
{
  std::vector<std::vector<LevelSetShapeSensitivity>> shapeSensitivities(myLSFields.size());
  for (unsigned i = 0; i < myLSFields.size(); ++i)
    fill_levelset_shape_sensitivities(myMeta.mesh_bulk_data(),
        myLSFields[i],
        interface_selector(myLSFields[i]),
        shapeSensitivities[i],
        method);
  return shapeSensitivities;
}

template <size_t NNODES>
InterfaceFacetsWithSensitivities<NNODES>
ConformingDecompositionWithSensitivities::facets_and_sensitivities(
    const unsigned lsIndex, const ShapeSensitivityMethod method) const
{
  STK_ThrowRequireMsg(NNODES == myMeta.spatial_dimension(),
      "Facets of " << NNODES << " nodes requested for a mesh of dimension "
                   << myMeta.spatial_dimension() << ".");
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const LS_Field & lsField = levelset_field(lsIndex);

  const stk::mesh::Selector interfaceSelector = interface_selector(lsField);
  InterfaceFacetsWithSensitivities<NNODES> facets;
  fill_levelset_shape_sensitivities(mesh, lsField, interfaceSelector, facets.sensitivities, method);
  sort_by_interface_node_id(facets.sensitivities);

  const Phase_Support & phaseSupport = Phase_Support::get(myMeta);
  const stk::mesh::Selector negativeSideElementSelector =
      phaseSupport.get_negative_levelset_block_selector(lsField.identifier);

  const std::vector<stk::mesh::Entity> interfaceSides =
      get_owned_interface_sides(mesh, interfaceSelector);
  facets.facetSensitivityIndices.reserve(interfaceSides.size());
  for (auto & interfaceSide : interfaceSides)
    facets.facetSensitivityIndices.push_back(get_facet_sensitivity_indices<NNODES>(
        mesh, negativeSideElementSelector, interfaceSide, facets.sensitivities));

  return facets;
}

// Explicit template instantiation

template struct InterfaceFacetsWithSensitivities<2>;
template struct InterfaceFacetsWithSensitivities<3>;

template InterfaceFacetsWithSensitivities<2>
ConformingDecompositionWithSensitivities::facets_and_sensitivities<2>(
    const unsigned lsIndex, const ShapeSensitivityMethod method) const;
template InterfaceFacetsWithSensitivities<3>
ConformingDecompositionWithSensitivities::facets_and_sensitivities<3>(
    const unsigned lsIndex, const ShapeSensitivityMethod method) const;

} // namespace krino
