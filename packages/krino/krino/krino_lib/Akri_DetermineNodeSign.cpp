#include <Akri_DetermineNodeSign.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Edge.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Sign.hpp>
#include <Akri_Surface.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

using NodeToSignMap = std::map<stk::mesh::Entity, int8_t>;

static
void pack_shared_nodes_and_signs(const stk::mesh::BulkData & mesh,
    const NodeToSignMap & nodesToSign,
    stk::CommSparse &commSparse)
{
  std::vector<int> commSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto & [node, sign] : nodesToSign)
    {
      if (mesh.bucket(node).shared())
      {
        mesh.comm_shared_procs(node, commSharedProcs);
        for (int procId : commSharedProcs)
        {
          if (procId != commSparse.parallel_rank())
          {
            commSparse.send_buffer(procId).pack(mesh.identifier(node));
            commSparse.send_buffer(procId).pack(sign);
          }
        }
      }
    }
  });
}

static
void unpack_shared_nodes_and_signs(const stk::mesh::BulkData & mesh,
    std::vector<std::pair<stk::mesh::Entity,int8_t>> & receivedNodesAndSigns,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId elemId;
      commSparse.recv_buffer(procId).unpack(elemId);
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, elemId);
      STK_ThrowRequire(mesh.is_valid(node));
      int8_t sign = 0;
      commSparse.recv_buffer(procId).unpack(sign);
      receivedNodesAndSigns.push_back(std::make_pair(node, sign));
    }
  });
}

void insert_remote_node_and_sign(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const int8_t sign, NodeToSignMap & nodesToSign)
{
  STK_ThrowAssert(sign == -1 || sign == 1);
  auto iter = nodesToSign.lower_bound(node);
  if (iter == nodesToSign.end() || iter->first != node)
    nodesToSign.emplace_hint(iter, node, sign);
  else
    STK_ThrowRequireMsg(sign == iter->second, "Attempt to change sign of shared node " << int(iter->second) << " -> " << int(sign) << debug_entity(mesh, node));
}

void communicate_nodes_and_signs(const stk::mesh::BulkData & mesh, NodeToSignMap & nodesToSign)
{
  if (mesh.parallel_size() == 1)
    return;

  stk::CommSparse commSparse(mesh.parallel());

  pack_shared_nodes_and_signs(mesh, nodesToSign, commSparse);

  std::vector<std::pair<stk::mesh::Entity,int8_t>> receivedNodesAndSigns;
  unpack_shared_nodes_and_signs(mesh, receivedNodesAndSigns, commSparse);

  for (auto & [node, sign] : receivedNodesAndSigns)
    insert_remote_node_and_sign(mesh, node, sign, nodesToSign);
}

static bool can_surface_compute_node_sign_anywhere(const Surface & surface)
{
  return !surface.truncated_distance_may_have_wrong_sign();
}

static NodeToSignsMap initialize_node_signs(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const size_t numSurfaces)
{
  stk::mesh::Selector nodeSelector = elementSelector & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  const std::vector<int8_t> initVals(numSurfaces, -3);
  NodeToSignsMap nodesToSigns;
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, nodeSelector))
    for (const auto & node : *bucketPtr)
      nodesToSigns[node] = initVals;
  return nodesToSigns;
}

static void assign_initial_node_sign(const stk::mesh::BulkData & mesh, const unsigned surfIndex, const stk::mesh::Entity node, const int8_t sign, NodeToSignMap & nodesToSign, NodeToSignsMap & nodesToAllSigns)
{
  STK_ThrowAssert(sign == -1 || sign == 0 || sign == 1);
  int8_t & nodeSign = nodesToAllSigns.at(node)[surfIndex];

  if (nodeSign != 0)
  {
    const int8_t initialSign = 2*sign; // 0, -2, or 2
    STK_ThrowRequireMsg(nodeSign == -3 || nodeSign == initialSign, "Attempt to change initial sign of node " << debug_entity(mesh, node));
    nodeSign = initialSign;

    if (sign != 0)
      nodesToSign[node] = sign;
  }
}

static void assign_node_sign_and_append_if_new(const stk::mesh::BulkData & mesh, const unsigned surfIndex, const stk::mesh::Entity node, const int8_t sign, NodeToSignMap & nodesToSign, NodeToSignsMap & nodesToAllSigns)
{
  STK_ThrowAssert(sign == -1 || sign == 1);
  int8_t & nodeSign = nodesToAllSigns.at(node)[surfIndex];

  const bool hasInterfaceSign = nodeSign == 0 || nodeSign == -2 || nodeSign == 2;
  if (!hasInterfaceSign)
  {
    if (nodeSign == -3)
      nodesToSign[node] = sign;

    STK_ThrowRequireMsg(nodeSign == -3 || nodeSign == sign, "Attempt to change sign of node " << debug_entity(mesh, node));
    nodeSign = sign;
  }
}

static bool is_node_snapped_to_level_set(const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const stk::mesh::Entity node,
    const int lsIndex)
{
  const auto iter = nodesToCapturedDomains.find(node);
  if (iter == nodesToCapturedDomains.end())
    return false;
  const std::vector<int> & nodeSortedDomains = iter->second;
  return std::binary_search(nodeSortedDomains.begin(), nodeSortedDomains.end(), lsIndex);
}

static stk::mesh::Entity get_noninterface_node_if_all_other_nodes_of_element_are_on_interface(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity elem,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const unsigned surfIndex)
{
  stk::mesh::Entity nonInterfaceNodeOfElem;
  for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
  {
    if (!is_node_snapped_to_level_set(nodesToCapturedDomains, elemNode, surfIndex))
    {
      if (nonInterfaceNodeOfElem.is_local_offset_valid())
        return stk::mesh::Entity::InvalidEntity;
      nonInterfaceNodeOfElem = elemNode;
    }
  }
  return nonInterfaceNodeOfElem;
}

static void assign_initial_node_signs_from_snapping(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const unsigned surfIndex,
    const Surface & surface,
    NodeToSignMap & nodesToSign,
    NodeToSignsMap & nodesToAllSigns)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  for (auto & [node, snappedDomains] : nodesToCapturedDomains)
  {
    if ((mesh.bucket(node).owned() || mesh.bucket(node).shared()) &&
        binary_search(snappedDomains.begin(), snappedDomains.end(), surfIndex))
    {
      assign_initial_node_sign(mesh, surfIndex, node, 0, nodesToSign, nodesToAllSigns);
    }
  }

  const stk::mesh::Selector ownedElemSelector = elementSelector & mesh.mesh_meta_data().locally_owned_part();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedElemSelector))
  {
    for (const auto & elem : *bucketPtr)
    {
      const stk::mesh::Entity nonInterfaceNodeOfElem = get_noninterface_node_if_all_other_nodes_of_element_are_on_interface(mesh, elem, nodesToCapturedDomains, surfIndex);
      if (mesh.is_valid(nonInterfaceNodeOfElem))
      {
        const int nodeSign = sign(surface.point_signed_distance(get_vector_field(mesh, coordsField, nonInterfaceNodeOfElem, dim)));
        assign_initial_node_sign(mesh, surfIndex, nonInterfaceNodeOfElem, nodeSign, nodesToSign, nodesToAllSigns);
      }
    }
  }
}

static void assign_initial_node_signs_from_edge_intersections(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const unsigned surfIndex,
    const Surface & surface,
    NodeToSignMap & nodesToSign,
    NodeToSignsMap & nodesToAllSigns)
{
  const double largeEdgeTolToKeepCostLow = 0.1;
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  const stk::mesh::Selector ownedElemSelector = elementSelector & mesh.mesh_meta_data().locally_owned_part();
  const std::vector<Edge> edges = get_edges_of_selected_elements(mesh, ownedElemSelector);
  for (auto & edge : edges)
  {
    const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
    const auto [crossingSign, position] = surface.compute_intersection_with_segment(get_vector_field(mesh, coordsField, edgeNodes[0], dim), get_vector_field(mesh, coordsField, edgeNodes[1], dim), largeEdgeTolToKeepCostLow);
    if (crossingSign != 0)
    {
      assign_initial_node_sign(mesh, surfIndex, edgeNodes[0],-crossingSign, nodesToSign, nodesToAllSigns);
      assign_initial_node_sign(mesh, surfIndex, edgeNodes[1], crossingSign, nodesToSign, nodesToAllSigns);
    }
  }
}

static void assign_node_signs_from_surface_distance(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const unsigned surfIndex,
    const Surface & surface,
    NodeToSignsMap & nodesToAllSigns)
{
  NodeToSignMap nodesToSign;

  if (!nodesToCapturedDomains.empty())
    assign_initial_node_signs_from_snapping(mesh, coordsField, elementSelector, nodesToCapturedDomains, surfIndex, surface, nodesToSign, nodesToAllSigns);

  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  for (auto & [node, allSigns] : nodesToAllSigns)
    if (is_entity_selected(mesh, elementSelector, node))
      allSigns[surfIndex] = sign(surface.point_signed_distance(get_vector_field(mesh, coordsField, node, dim)));
}

static void get_initial_node_signs_from_snapped_nodes_and_edge_intersections(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const unsigned surfIndex,
    const Surface & surface,
    NodeToSignMap & nodesToSign,
    NodeToSignsMap & nodesToAllSigns)
{
  nodesToSign.clear();

  if (!nodesToCapturedDomains.empty())
    assign_initial_node_signs_from_snapping(mesh, coordsField, elementSelector, nodesToCapturedDomains, surfIndex, surface, nodesToSign, nodesToAllSigns);

  assign_initial_node_signs_from_edge_intersections(mesh, coordsField, elementSelector, surfIndex, surface, nodesToSign, nodesToAllSigns);
}

static stk::mesh::Entity first_selected_entity(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & entitySelector, const stk::mesh::EntityRank entityRank)
{
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, entitySelector))
    for (const auto & entity : *bucketPtr)
      return entity;
  return stk::mesh::Entity::InvalidEntity;
}

static int get_sign_of_first_node(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const Surface & surface)
{
  stk::mesh::Entity firstNode = first_selected_entity(mesh, elementSelector, stk::topology::NODE_RANK);

  int nodeSign = 0;
  if (mesh.is_valid(firstNode))
  {
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    nodeSign = sign(surface.point_signed_distance(get_vector_field(mesh, coordsField, firstNode, dim)));
  }

  int minNodeSign = nodeSign;
  stk::all_reduce_min(mesh.parallel(), &nodeSign, &minNodeSign, 1);
  int maxNodeSign = nodeSign;
  stk::all_reduce_max(mesh.parallel(), &nodeSign, &maxNodeSign, 1);

  STK_ThrowRequireMsg(!(minNodeSign == -1 && maxNodeSign == 1), "Disagreement on node sign in parallel.");
  if (minNodeSign == -1)
    return -1;
  STK_ThrowRequireMsg(maxNodeSign == 1, "Did not get a valid node sign in parallel.");
  return 1;
}

static void assign_all_nodes_the_sign_of_first_node(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const unsigned surfIndex,
    const Surface & surface,
    NodeToSignsMap & nodesToAllSigns)
{
  const int nodeSign = get_sign_of_first_node(mesh, coordsField, elementSelector, surface);

  for (auto & [node, allSigns] : nodesToAllSigns)
    if (is_entity_selected(mesh, elementSelector, node))
      allSigns[surfIndex] = nodeSign;
}

static void fill_neighbor_nodes_of_selected_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Selector & elementSelector, std::vector<stk::mesh::Entity> & neighbors)
{
  //NOTE: This is for simplex elements
  const unsigned numVertexNodes = mesh.mesh_meta_data().spatial_dimension() + 1;
  neighbors.clear();
  for (auto && elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    if (is_entity_selected(mesh, elementSelector, elem))
      for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.begin_nodes(elem)+numVertexNodes})
        if (elemNode != node && (mesh.bucket(node).owned() || mesh.bucket(node).shared()))
          neighbors.push_back(elemNode);
  stk::util::sort_and_unique(neighbors);
}

void assign_nodes_on_either_side_of_surface(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const unsigned surfIndex,
    NodeToSignMap & nodesToSign,
    NodeToSignsMap & nodesToAllSigns)
{
  std::vector<stk::mesh::Entity> nodeNbrs;
  const stk::mesh::Selector ownedElemSelector = elementSelector & mesh.mesh_meta_data().locally_owned_part();

  while (stk::is_true_on_any_proc(mesh.parallel(), !nodesToSign.empty()))
  {
    communicate_nodes_and_signs(mesh, nodesToSign);

    NodeToSignMap nextIterNodesToSign;

    for (auto & [node, sign] : nodesToSign)
    {
      fill_neighbor_nodes_of_selected_elements(mesh, node, ownedElemSelector, nodeNbrs);

      for (auto & nbrNode : nodeNbrs)
        assign_node_sign_and_append_if_new(mesh, surfIndex, nbrNode, sign, nextIterNodesToSign, nodesToAllSigns);
    }

    nodesToSign.swap(nextIterNodesToSign);
  }
}

NodeToSignsMap determine_node_signs(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Selector> & perSurfaceElementSelector,
    const std::vector<const Surface*> & surfaces,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  STK_ThrowAssert(perSurfaceElementSelector.size() == surfaces.size());
  NodeToSignsMap nodesToAllSigns = initialize_node_signs(mesh, elementSelector, surfaces.size());

  for (unsigned surfIndex=0; surfIndex<surfaces.size(); ++surfIndex)
  {
    stk::mesh::Selector elemSurfSelector = elementSelector & perSurfaceElementSelector[surfIndex];
    const Surface & surface = *surfaces[surfIndex];
    if (can_surface_compute_node_sign_anywhere(surface))
    {
      assign_node_signs_from_surface_distance(mesh, coordsField, elemSurfSelector, nodesToCapturedDomains, surfIndex, surface, nodesToAllSigns);
    }
    else
    {
      NodeToSignMap initialNodesToSign;
      get_initial_node_signs_from_snapped_nodes_and_edge_intersections(mesh, coordsField, elemSurfSelector, nodesToCapturedDomains, surfIndex, surface, initialNodesToSign, nodesToAllSigns);
      if (stk::is_true_on_all_procs(mesh.parallel(), initialNodesToSign.empty()))
        assign_all_nodes_the_sign_of_first_node(mesh, coordsField, elemSurfSelector, surfIndex, surface, nodesToAllSigns);
      else
        assign_nodes_on_either_side_of_surface(mesh, elemSurfSelector, surfIndex, initialNodesToSign, nodesToAllSigns);
    }
  }

  return nodesToAllSigns;
}


}

