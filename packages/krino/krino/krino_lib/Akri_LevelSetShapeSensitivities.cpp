#include "Akri_Snap.hpp"
#include <Akri_LevelSetShapeSensitivities.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ChildNodeStencil.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <Akri_QualityMetric.hpp>

namespace krino {

std::string output_sensitivity(const LevelSetShapeSensitivity & sens)
{
  std::ostringstream os;
  os << "Sensitivities for node " << sens.interfaceNodeId << ":\n";
  for (size_t i=0; i<sens.parentNodeIds.size(); ++i)
    os << "  dCoordParentLevelSet for parent node " << sens.parentNodeIds[i] << " = " << sens.dCoordsdParentLevelSets[i] << "\n";
  return os.str();
}

static std::array<double,2> compute_d_loc_d_levelsets_for_edge(const FieldRef levelSetField, const std::array<stk::mesh::Entity,2> & parentEdgeNodes)
{
  const double ls0 = *field_data<double>(levelSetField, parentEdgeNodes[0]);
  const double ls1 = *field_data<double>(levelSetField, parentEdgeNodes[1]);
  const double sqrLo = (ls0-ls1)*(ls0-ls1);
  return {{(-ls1/sqrLo), (ls0/sqrLo)}};
}

static stk::math::Vector3d compute_edge_vector(const FieldRef coordsField, const std::array<stk::mesh::Entity,2> & parentEdgeNodes)
{
  const stk::math::Vector3d x0(field_data<double>(coordsField, parentEdgeNodes[0]));
  const stk::math::Vector3d x1(field_data<double>(coordsField, parentEdgeNodes[1]));
  return x1-x0;
}

static void fill_d_coords_d_levelsets(const FieldRef coordsField, const FieldRef levelSetField, const std::vector<stk::mesh::Entity> & parentNodes, const std::vector<double> & /*parentWeights*/, std::vector<stk::math::Vector3d> & dCoordsdParentLevelSets)
{
  STK_ThrowRequireMsg(2 == parentNodes.size(), "Currently only edge intersections are supported.");
  dCoordsdParentLevelSets.clear();
  const std::array<double,2> dLocdLs = compute_d_loc_d_levelsets_for_edge(levelSetField, {{parentNodes[0], parentNodes[1]}});
  const stk::math::Vector3d dx = compute_edge_vector(coordsField, {{parentNodes[0], parentNodes[1]}});
  dCoordsdParentLevelSets.push_back(dLocdLs[0] * dx);
  dCoordsdParentLevelSets.push_back(dLocdLs[1] * dx);
}

void append_sensitivities_for_child_nodes(const stk::mesh::BulkData & mesh, const FieldRef levelSetField, std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  FieldRef lsFieldToUse = levelSetField;
  auto stashedLSFld = cdfemSupport.get_stashed_levelsets().find(levelSetField);
  if(stashedLSFld != cdfemSupport.get_stashed_levelsets().end()) lsFieldToUse = stashedLSFld->second;

  const FieldRef coordsField = cdfemSupport.get_coords_field();
  std::vector<ChildNodeStencil> childNodeStencils;
  fill_child_node_stencils(mesh, cdfemSupport.get_child_node_part(), cdfemSupport.get_parent_node_ids_field(), cdfemSupport.get_parent_node_weights_field(), childNodeStencils);

  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;

  shapeSensitivities.reserve(childNodeStencils.size());
  for (auto & stencil : childNodeStencils)
  {
    fill_node_ids_for_nodes(mesh, stencil.parentNodes, parentNodeIds);
    fill_d_coords_d_levelsets(coordsField, lsFieldToUse, stencil.parentNodes, stencil.parentWeights, dCoordsdParentLevelSets);
    shapeSensitivities.emplace_back(mesh.identifier(stencil.childNode), parentNodeIds, dCoordsdParentLevelSets);
  }
}

static void fill_neighbor_nodes_not_on_interface(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Selector & interfaceSelector, std::vector<stk::mesh::Entity> & neighbors)
{
  //NOTE: This is for linear simplex elements
  STK_ThrowAssert(mesh.is_automatic_aura_on());
  neighbors.clear();
  for (auto && elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)}) //NOTE: This does not limit elements to active ones or anything like that
      if (elemNode != node && !interfaceSelector(mesh.bucket(elemNode)))
        neighbors.push_back(elemNode);
  stk::util::sort_and_unique(neighbors);
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

static bool first_wins_tiebreaker(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity nodeA, stk::mesh::Entity nodeB)
{
  STK_ThrowAssert(nodeA.is_local_offset_valid());
  if (!nodeB.is_local_offset_valid())
    return true;
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  return is_less_than_in_x_then_y_then_z(get_vector_field(mesh, coordsField, nodeA, dim), get_vector_field(mesh, coordsField, nodeB, dim));
}

static std::pair<bool, double> is_node_on_edge_and_location(const unsigned dim,
    const stk::math::Vector3d & loc,
    const stk::math::Vector3d & edgeLoc0,
    const stk::math::Vector3d & edgeLoc1,
    const double tol)
{
  const double pos = compute_child_position(dim, loc, edgeLoc0, edgeLoc1);
  const stk::math::Vector3d edgeLoc = (1.-pos)*edgeLoc0 + pos*edgeLoc1;
  const bool isOnEdge = (pos>=0. && pos<=1.) && ((edgeLoc-loc).length_squared() < tol*tol*(edgeLoc1-edgeLoc0).length_squared());
  return std::make_pair(isOnEdge, pos);
}

stk::mesh::Entity find_other_parent_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const FieldRef isovar,
    stk::mesh::Entity node)
{
  const double tol = 1e-8;
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const double nodeLs = *field_data<double>(isovar, node);
  const stk::math::Vector3d nodeCurrentCoords = get_vector_field(mesh, coordsField, node, dim);
  const stk::math::Vector3d nodeUnsnappedCoords = nodeCurrentCoords - get_vector_field(mesh, snapDisplacements, node, dim);

  std::vector<stk::mesh::Entity> neighborNodes;
  fill_neighbor_nodes_of_selected_elements(mesh, node, elemSelector, neighborNodes);

  stk::mesh::Entity otherNode;
  stk::mesh::Entity bestArbitraryOtherNode;
  double maxDeltaLS = 0.;

  for (auto & nbrNode : neighborNodes)
  {
    const stk::math::Vector3d nbrNodeUnsnappedCoords = get_vector_field(mesh, coordsField, nbrNode, dim) - get_vector_field(mesh, snapDisplacements, nbrNode, dim);

    const auto & [isOnEdge, edgeLoc] = is_node_on_edge_and_location(dim, nodeCurrentCoords, nodeUnsnappedCoords, nbrNodeUnsnappedCoords, tol);
    if (isOnEdge)
    {
      if (edgeLoc > tol)
      {
        STK_ThrowRequireMsg(!otherNode.is_local_offset_valid(), "find_other_parent_node: more than 2 nodes found with nonzero weight");
        otherNode = nbrNode;
      }
      else
      {
        const double nbrNodeLs = *field_data<double>(isovar, nbrNode);
        const double deltaLS = std::abs(nbrNodeLs-nodeLs);
        if((deltaLS > maxDeltaLS) || (deltaLS == maxDeltaLS && first_wins_tiebreaker(mesh, coordsField, nbrNode, bestArbitraryOtherNode)))
        {
          bestArbitraryOtherNode = nbrNode;
          maxDeltaLS = deltaLS;
        }
      }
    }
  }

  if(!otherNode.is_local_offset_valid())
  {
    STK_ThrowRequireMsg(bestArbitraryOtherNode.is_local_offset_valid(), "find_other_parent_node: couldn't find node to pick");
    otherNode = bestArbitraryOtherNode;
  }
  return otherNode;
}

void calculate_sensitivities_for_snapped_nodes_from_parents(const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & elemSelector,
  const FieldRef coordsField,
  const FieldRef snapDisplacements,
  const FieldRef isovar,
  stk::mesh::Entity node,
  std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const auto & LSStash = cdfemSupport.get_stashed_levelsets();
  STK_ThrowRequireMsg(LSStash.find(isovar) != LSStash.end(), 
    "calculate_sensitivities_for_snapped_nodes_from_parents: no stashed levelset field found. Please ensure "
    "CDFEM_Support::setup_levelset_field_stash is called before decomposition");

  const auto & stashedLS = LSStash.at(isovar);

  stk::mesh::Entity otherNode = find_other_parent_node(mesh, elemSelector, coordsField, snapDisplacements, stashedLS, node);

  std::vector<stk::mesh::EntityId> parentNodeIds{mesh.identifier(node), mesh.identifier(otherNode)};
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets(2);

  const std::array<double,2> dLocdLs = compute_d_loc_d_levelsets_for_edge(stashedLS, {{node, otherNode}});
  const stk::math::Vector3d thisNodeUnsnappedLoc = get_vector_field(mesh, coordsField, node, dim) -
    get_vector_field(mesh, snapDisplacements, node, dim);
  const stk::math::Vector3d otherNodeUnsnappedLoc = get_vector_field(mesh, coordsField, otherNode, dim) -
    get_vector_field(mesh, snapDisplacements, otherNode, dim);
  const stk::math::Vector3d dx = otherNodeUnsnappedLoc - thisNodeUnsnappedLoc;

  dCoordsdParentLevelSets[0] = (dLocdLs[0] * dx);
  dCoordsdParentLevelSets[1] = (dLocdLs[1] * dx);

  shapeSensitivities.emplace_back(mesh.identifier(node), parentNodeIds, dCoordsdParentLevelSets);
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

  std::vector<stk::mesh::EntityId> parentNodeIds(1);
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets(1);
  std::vector<stk::mesh::Entity> neighbors;

  const bool calcSensFromParents = (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE) &&
    cdfemSupport.get_cdfem_snap_displacements_field().valid();

  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, ownershipSelector & interfaceNotChildSelector ) )
  {
    for (auto node : *bucket)
    {
      const stk::math::Vector3d x0(field_data<double>(coordsField, node));
      if(calcSensFromParents)
      {
        calculate_sensitivities_for_snapped_nodes_from_parents(mesh, cutElementSelector, coordsField, snapDisplacements, lsField.isovar, node, shapeSensitivities);
      }
      else 
      {
        fill_neighbor_nodes_not_on_interface(mesh, node, interfaceSelector, neighbors);

        dCoordsdParentLevelSets[0] = stk::math::Vector3d::ZERO;
        for (auto nbr : neighbors)
        {
          const stk::math::Vector3d x1(field_data<double>(coordsField, nbr));
          const double ls1 = *field_data<double>(lsField.isovar, nbr);
          dCoordsdParentLevelSets[0] += -1./ls1 * (x1-x0);
        }

        dCoordsdParentLevelSets[0] /= neighbors.size();
        stk::mesh::EntityId nodeId = mesh.identifier(node);
        parentNodeIds[0] = nodeId;
        shapeSensitivities.emplace_back(nodeId, parentNodeIds, dCoordsdParentLevelSets);
      }
    }
  }
}

void fill_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  shapeSensitivities.clear();
  append_sensitivities_for_child_nodes(mesh, lsField.isovar, shapeSensitivities);
  append_sensitivities_for_background_nodes_on_interface(mesh, lsField, shapeSensitivities);
}

std::vector<LevelSetShapeSensitivity> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & lsFields)
{
  std::vector<LevelSetShapeSensitivity> shapeSensitivities;

  STK_ThrowRequireMsg(1 == lsFields.size(), "Currently only one level set is supported.");
  fill_levelset_shape_sensitivities(mesh, lsFields[0], shapeSensitivities);

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
     std::vector<std::array<size_t,NNODES>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities)
{
  fill_levelset_shape_sensitivities(mesh, lsField, sortedShapeSensitivities);

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

// Explicit template instantiation

template void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
     const LS_Field & lsField,
     std::vector<std::array<size_t,2>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);
template void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
     const LS_Field & lsField,
     std::vector<std::array<size_t,3>> & facetSensitivityIndices,
     std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);

template std::array<const LevelSetShapeSensitivity*,2> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,2> & facetSensIndices);
template std::array<const LevelSetShapeSensitivity*,3> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,3> & facetSensIndices);

template std::array<stk::mesh::Entity,2> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,2> & facetNodeSens);
template std::array<stk::mesh::Entity,3> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,3> & facetNodeSens);

}

