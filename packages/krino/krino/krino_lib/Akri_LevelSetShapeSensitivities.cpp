#include "Akri_Snap.hpp"
#include <Akri_LevelSetShapeSensitivities.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ChildNodeStencil.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>

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
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  std::vector<ChildNodeStencil> childNodeStencils;
  fill_child_node_stencils(mesh, cdfemSupport.get_child_node_part(), cdfemSupport.get_parent_node_ids_field(), cdfemSupport.get_parent_node_weights_field(), childNodeStencils);

  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;

  shapeSensitivities.reserve(childNodeStencils.size());
  for (auto & stencil : childNodeStencils)
  {
    fill_node_ids_for_nodes(mesh, stencil.parentNodes, parentNodeIds);
    fill_d_coords_d_levelsets(coordsField, levelSetField, stencil.parentNodes, stencil.parentWeights, dCoordsdParentLevelSets);
    shapeSensitivities.emplace_back(mesh.identifier(stencil.childNode), parentNodeIds, dCoordsdParentLevelSets);
  }
}

static void fill_neighbor_nodes_not_on_interface(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Selector & interfaceSelector, std::vector<stk::mesh::Entity> & neighbors)
{
  STK_ThrowAssert(mesh.is_automatic_aura_on());
  neighbors.clear();
  for (auto && elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)}) //NOTE: This does not limit elements to active ones or anything like that
      if (elemNode != node && !interfaceSelector(mesh.bucket(elemNode)))
        neighbors.push_back(elemNode);
  stk::util::sort_and_unique(neighbors);
}

stk::mesh::Entity find_other_parent_node(stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & interpNodes,
  const std::vector<double> & interpWeights, FieldRef isovar)
{
  const double tol = 1e-8;
  const double ls0 = *field_data<double>(isovar, node);

  stk::mesh::Entity otherNode;
  int arbitraryNodeToPick = -1;
  double maxDeltaLS = 0.;
  for(unsigned i=0; i<interpNodes.size(); i++)
  {
    if(interpNodes[i] == node)
    {
      STK_ThrowRequireMsg(interpWeights[i] > tol, "find_other_parent_node: weight of current node should be > 0");
      continue;
    }
    else if(interpWeights[i] > tol)
    {
      STK_ThrowRequireMsg(!otherNode.is_local_offset_valid(), "find_other_parent_node: more than 2 nodes found with nonzero weight");
      otherNode = interpNodes[i];
    }
    else
    {
      const double ls1 = *field_data<double>(isovar, interpNodes[i]);
      const double deltaLS = std::fabs(ls1-ls0);
      if(deltaLS > maxDeltaLS) 
      {
        arbitraryNodeToPick = i;
        maxDeltaLS = deltaLS;
      }
    }
  }
  if(!otherNode.is_local_offset_valid()) 
  {
    STK_ThrowRequireMsg(arbitraryNodeToPick != -1, "find_other_parent_node: couldn't find node to pick");
    otherNode = interpNodes[arbitraryNodeToPick];
  }
  return otherNode;
}

void calcuate_sensitivities_for_snapped_nodes_from_parents(const stk::mesh::BulkData & mesh, const CDFEM_Support & cdfemSupport, 
  const LS_Field & lsField, stk::mesh::Entity node, std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  const FieldRef snapDisplacements = cdfemSupport.get_cdfem_snap_displacements_field();
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const stk::math::Vector3d snapDisp = get_vector_field(mesh, snapDisplacements, node, dim);
  const stk::math::Vector3d currentLoc = get_vector_field(mesh, coordsField, node, dim);
  const stk::math::Vector3d thisNodeUnsnappedLoc = currentLoc - snapDisp;
  std::vector<stk::mesh::Entity> interpNodes;
  std::vector<double> interpWeights;
  fill_interplation_nodes_and_weights_at_location(mesh, cdfemSupport.get_active_part(), coordsField, node, 
    thisNodeUnsnappedLoc, interpNodes, interpWeights);

  stk::mesh::Entity otherNode = find_other_parent_node(node, interpNodes, interpWeights, lsField.isovar);

  std::vector<stk::mesh::EntityId> parentNodeIds{mesh.identifier(node), mesh.identifier(otherNode)};
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets(2);

  const std::array<double,2> dLocdLs = compute_d_loc_d_levelsets_for_edge(lsField.isovar, {{node, otherNode}});
  const stk::math::Vector3d otherNodeUnsnappedLoc = get_vector_field(mesh, coordsField, otherNode, dim) -
    get_vector_field(mesh, snapDisplacements, otherNode, dim);
  const stk::math::Vector3d dx = otherNodeUnsnappedLoc - thisNodeUnsnappedLoc;

  dCoordsdParentLevelSets[0] = (dLocdLs[0] * dx);
  dCoordsdParentLevelSets[1] = (dLocdLs[1] * dx);

  shapeSensitivities.emplace_back(mesh.identifier(node), parentNodeIds, dCoordsdParentLevelSets);
}

void append_sensitivities_for_background_nodes_on_interface(const stk::mesh::BulkData & mesh, const LS_Field & lsField, std::vector<LevelSetShapeSensitivity> & shapeSensitivities)
{
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsField.identifier);
  const stk::mesh::Selector interfaceNotChildSelector = (!cdfemSupport.get_child_node_part()) & interfaceSelector;
  const FieldRef coordsField = cdfemSupport.get_coords_field();

  std::vector<stk::mesh::EntityId> parentNodeIds(1);
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets(1);
  std::vector<stk::mesh::Entity> neighbors;

  const bool calcSensFromParents = (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE) &&
    cdfemSupport.get_resnap_method() == Resnap_Method::RESNAP_USING_INTERFACE_ON_PREVIOUS_SNAPPED_MESH &&
    cdfemSupport.get_cdfem_snap_displacements_field().valid();

  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part() & interfaceNotChildSelector ) )
  {
    for (auto node : *bucket)
    {
      const stk::math::Vector3d x0(field_data<double>(coordsField, node));
      if(calcSensFromParents)
      {
        calcuate_sensitivities_for_snapped_nodes_from_parents(mesh, cdfemSupport, lsField, node, shapeSensitivities);
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

std::vector<LevelSetShapeSensitivity> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & lsFields)
{
  STK_ThrowRequireMsg(1 == lsFields.size(), "Currently only one level set is supported.");
  std::vector<LevelSetShapeSensitivity> shapeSensitivities;
  append_sensitivities_for_child_nodes(mesh, lsFields[0].isovar, shapeSensitivities);
  append_sensitivities_for_background_nodes_on_interface(mesh, lsFields[0], shapeSensitivities);

  return shapeSensitivities;
}

}

