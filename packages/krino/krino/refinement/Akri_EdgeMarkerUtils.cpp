#include <Akri_EdgeMarkerUtils.hpp>

#include <Akri_MeshHelpers.hpp>
#include <Akri_Refinement.hpp>

namespace krino {
namespace EdgeMarkerUtils {

static std::array<stk::mesh::EntityId,2> get_sorted_edge_node_ids(const stk::mesh::BulkData & mesh, const stk::mesh::Entity edgeNode0, const stk::mesh::Entity edgeNode1)
{
  const stk::mesh::EntityId edgeNode0Id = mesh.identifier(edgeNode0);
  const stk::mesh::EntityId edgeNode1Id = mesh.identifier(edgeNode1);
  return (edgeNode0Id<edgeNode1Id) ? std::array<stk::mesh::EntityId,2>{{edgeNode0Id, edgeNode1Id}} : std::array<stk::mesh::EntityId,2>{{edgeNode1Id, edgeNode0Id}};
}

bool node_is_parent_node(const StkMeshEntities & parentNodes, const stk::mesh::Entity node)
{
  return (std::find(parentNodes.begin(), parentNodes.end(), node) != parentNodes.end());
}

bool is_any_child_invalid_or_a_parent(const stk::mesh::BulkData & mesh, const Refinement & refinement, const std::vector<stk::mesh::Entity> & childElements)
{
  for (auto && childElem : childElements)
    if (!mesh.is_valid(childElem) || refinement.is_parent(childElem))
      return true;
  return false;
}

bool is_element_a_parent_but_not_a_grandparent(const stk::mesh::BulkData & mesh, const Refinement & refinement, const stk::mesh::Entity elem)
{
  if (refinement.is_parent(elem))
  {
    const std::vector<stk::mesh::Entity> childElems = refinement.get_children(elem);
    return !is_any_child_invalid_or_a_parent(mesh, refinement, childElems);
  }
  return false;
}

static void fill_refined_edge_nodes_for_marked_edges(const stk::mesh::BulkData & mesh,
    const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    bool & areAnyEdgesMarked)
{
  areAnyEdgesMarked = false;
  const unsigned elemNumEdges = elemTopology.num_edges();
  elemEdgeChildNodes.assign(elemNumEdges, stk::mesh::Entity());

  const stk::mesh::Entity * entityNodes = mesh.begin_nodes(elem);
  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, iEdge);
    const Edge edge = edge_from_edge_nodes(mesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
    const stk::mesh::Entity refinedEdgeNode = nodeRefiner.get_edge_child_node(edge);
    if (mesh.is_valid(refinedEdgeNode))
    {
      areAnyEdgesMarked = true;
      elemEdgeChildNodes[iEdge] = refinedEdgeNode;
    }
  }
}

static void keep_child_nodes_that_are_not_parent_nodes(const StkMeshEntities & parentNodes, std::vector<stk::mesh::Entity> & childNodesThatAreNotParentNodes)
{
  size_t numKeep = 0;
  for (size_t iChild=0; iChild<childNodesThatAreNotParentNodes.size(); ++iChild)
    if (!node_is_parent_node(parentNodes, childNodesThatAreNotParentNodes[iChild]))
      childNodesThatAreNotParentNodes[numKeep++] = childNodesThatAreNotParentNodes[iChild];
  childNodesThatAreNotParentNodes.resize(numKeep);
}

std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElems)
{
  std::vector<stk::mesh::Entity> childNodesThatAreNotParentNodes;
  if (!childElems.empty())
  {
    const StkMeshEntities parentNodes {mesh.begin_nodes(parentElem), mesh.end_nodes(parentElem)};
    childNodesThatAreNotParentNodes.reserve(parentNodes.size()*childElems.size());
    for (auto && childElem : childElems)
    {
      STK_ThrowAssertMsg(mesh.is_valid(childElem), "All child elements must by valid.  Should this element really be considered for adaptation?");
      StkMeshEntities childNodes{mesh.begin_nodes(childElem), mesh.end_nodes(childElem)};
      childNodesThatAreNotParentNodes.insert(childNodesThatAreNotParentNodes.end(), childNodes.begin(), childNodes.end());
    }
    stk::util::sort_and_unique(childNodesThatAreNotParentNodes);

    keep_child_nodes_that_are_not_parent_nodes(parentNodes, childNodesThatAreNotParentNodes);
  }
  return childNodesThatAreNotParentNodes;
}

std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const stk::mesh::Entity parentElem)
{
  const std::vector<stk::mesh::Entity> childElems = refinement.get_children(parentElem);
  return get_child_nodes_that_are_not_parent_nodes(mesh, parentElem, childElems);
}

static std::vector<std::array<stk::mesh::EntityId,2>>
get_parent_edge_node_ids_for_given_refined_edge_nodes(const Refinement & refinement, const std::vector<stk::mesh::Entity> & refinedEdgeNodes)
{
  std::vector<std::array<stk::mesh::EntityId,2>> refinedEdgeParentNodeIds;
  refinedEdgeParentNodeIds.reserve(refinedEdgeNodes.size());
  for (auto && refinedEdgeNode : refinedEdgeNodes)
  {
    const std::array<stk::mesh::EntityId,2> edgeParentNodeIds = refinement.get_edge_parent_node_ids(refinedEdgeNode);
    refinedEdgeParentNodeIds.push_back(edgeParentNodeIds);
  }
  return refinedEdgeParentNodeIds;
}

static void fill_in_existing_refined_edge_nodes_and_caseIds(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    ElementEdgeCaseIds & elementEdgeCaseIds)
{
  const unsigned elemNumEdges = elemTopology.num_edges();

  int preAdaptCaseId = 0;
  int postAdaptCaseId = 0;

  if (refinement.is_parent(elem))
  {
    const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_child_nodes_that_are_not_parent_nodes(mesh, refinement, elem);
    const std::vector<std::array<stk::mesh::EntityId,2>> refinedEdgeParentNodeIds = get_parent_edge_node_ids_for_given_refined_edge_nodes(refinement, refinedEdgeNodes);

    const stk::mesh::Entity * entityNodes = mesh.begin_nodes(elem);
    for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
    {
      const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, iEdge);
      const std::array<stk::mesh::EntityId,2> edgeNodeIds = get_sorted_edge_node_ids(mesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
      const auto iter = std::find(refinedEdgeParentNodeIds.begin(), refinedEdgeParentNodeIds.end(), edgeNodeIds);
      if (iter != refinedEdgeParentNodeIds.end())
      {
        preAdaptCaseId += 1<<iEdge;
        elemEdgeChildNodes[iEdge] = refinedEdgeNodes[std::distance(refinedEdgeParentNodeIds.begin(), iter)];
      }
    }
  }

  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
    if (elemEdgeChildNodes[iEdge] != stk::mesh::Entity::InvalidEntity)
      postAdaptCaseId += 1<<iEdge;

  elementEdgeCaseIds.set(preAdaptCaseId, postAdaptCaseId);
}

static void fill_post_unrefinement_edge_nodes_and_caseIds(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const std::vector<stk::mesh::Entity> & refinedEdgeNodes,
    const std::vector<bool> & areEdgeNodesBeingUnrefined,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    ElementEdgeCaseIds & elementEdgeCaseIds)
{
  const std::vector<std::array<stk::mesh::EntityId,2>> refinedEdgeParentNodeIds = get_parent_edge_node_ids_for_given_refined_edge_nodes(refinement, refinedEdgeNodes);

  int preAdaptCaseId = 0;
  int postAdaptCaseId = 0;
  const unsigned elemNumEdges = elemTopology.num_edges();
  elemEdgeChildNodes.assign(elemNumEdges, stk::mesh::Entity());

  const stk::mesh::Entity * entityNodes = mesh.begin_nodes(elem);
  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, iEdge);
    const std::array<stk::mesh::EntityId,2> edgeNodeIds = get_sorted_edge_node_ids(mesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
    const auto iter = std::find(refinedEdgeParentNodeIds.begin(), refinedEdgeParentNodeIds.end(), edgeNodeIds);
    if (iter != refinedEdgeParentNodeIds.end())
    {
      preAdaptCaseId += 1<<iEdge;
      const size_t index = std::distance(refinedEdgeParentNodeIds.begin(), iter);
      if (!areEdgeNodesBeingUnrefined[index])
      {
        postAdaptCaseId += 1<<iEdge;
        elemEdgeChildNodes[iEdge] = refinedEdgeNodes[index];
      }
    }
  }

  elementEdgeCaseIds.set(preAdaptCaseId, postAdaptCaseId);
}

static void fill_are_edge_nodes_being_unrefined(const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement,
    const std::vector<stk::mesh::Entity> & refinedEdgeNodes,
    std::vector<bool> & areEdgeNodesBeingUnrefined,
    bool & areAnyEdgeNodesBeingUnrefined)
{
  areAnyEdgeNodesBeingUnrefined = false;
  areEdgeNodesBeingUnrefined.assign(refinedEdgeNodes.size(), false);
  for (size_t iNode=0; iNode<refinedEdgeNodes.size(); ++iNode)
  {
    const auto iter = std::lower_bound(sortedOwnedOrSharedNodesToBeRemovedByUnrefinement.begin(), sortedOwnedOrSharedNodesToBeRemovedByUnrefinement.end(), refinedEdgeNodes[iNode]);
    if (iter != sortedOwnedOrSharedNodesToBeRemovedByUnrefinement.end() && *iter == refinedEdgeNodes[iNode])
    {
      areAnyEdgeNodesBeingUnrefined = true;
      areEdgeNodesBeingUnrefined[iNode] = true;
    }
  }
}

void fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const bool isElemCandidateForRefinement,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes)
{
  elementEdgeCaseIds.clear();

  if (isElemCandidateForRefinement)
  {
    bool areAnyEdgesMarked = false;
    fill_refined_edge_nodes_for_marked_edges(mesh, nodeRefiner, elem, elemTopology, elemEdgeChildNodes, areAnyEdgesMarked);
    if (areAnyEdgesMarked)
      fill_in_existing_refined_edge_nodes_and_caseIds(mesh, refinement, elem, elemTopology, elemEdgeChildNodes, elementEdgeCaseIds);
    else
      elemEdgeChildNodes.clear();
  }
}

void fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const bool isElemCandidateForUnrefinement,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes)
{
  elementEdgeCaseIds.clear();
  elemEdgeChildNodes.clear();
  if (isElemCandidateForUnrefinement)
  {
    const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_child_nodes_that_are_not_parent_nodes(mesh, refinement, elem);
    const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement = nodeRefiner.get_sorted_edge_nodes_that_will_be_removed_by_unrefinement();
    std::vector<bool> areEdgeNodesBeingUnrefined;
    bool areAnyEdgeNodesBeingUnrefined = false;
    fill_are_edge_nodes_being_unrefined(sortedOwnedOrSharedNodesToBeRemovedByUnrefinement, refinedEdgeNodes, areEdgeNodesBeingUnrefined, areAnyEdgeNodesBeingUnrefined);

    if (areAnyEdgeNodesBeingUnrefined)
      fill_post_unrefinement_edge_nodes_and_caseIds(mesh, refinement, elem, elemTopology, refinedEdgeNodes, areEdgeNodesBeingUnrefined, elemEdgeChildNodes, elementEdgeCaseIds);
  }
}

}
}


