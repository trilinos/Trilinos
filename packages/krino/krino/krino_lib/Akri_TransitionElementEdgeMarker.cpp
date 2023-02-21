/*
 * Akri_TransitionElementEdgeMarker.cpp
 *
 *  Created on: Oct 28, 2022
 *      Author: drnoble
 */
#include "Akri_TransitionElementEdgeMarker.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/RemoveIntersection.hpp>
#include "Akri_Edge.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_Refinement.hpp"

namespace krino {

UniformEdgeMarker::UniformEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement)
: myMesh(mesh), myRefinement(refinement)
{
}

void UniformEdgeMarker::locally_mark_edges_of_non_parent_elements(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & !myRefinement.parent_part();
  std::vector<Edge> edgesToRefineForElement;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_entity_edges(myMesh, elem, edgesToRefineForElement);
      for (auto && edgeToRefineForElement : edgesToRefineForElement)
        nodeRefiner.mark_edge_for_refinement(edgeToRefineForElement);
    }
  }
}

void UniformEdgeMarker::mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const
{
  locally_mark_edges_of_non_parent_elements(nodeRefiner);
  nodeRefiner.sync_shared_edges_from_other_procs_to_refine(myMesh);
}

bool UniformEdgeMarker::is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const
{
  return !myRefinement.is_parent(elem);
}

void UniformEdgeMarker::fill_element_refined_edge_nodes(const NodeRefiner & nodeRefiner, const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const unsigned elemNumEdges = elemTopology.num_edges();
  const stk::mesh::Entity * elemNodes = myMesh.begin_nodes(elem);
  std::array<stk::mesh::Entity,3> edgeNodes;
  elemEdgeChildNodes.resize(elemNumEdges);

  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
  {
    elemTopology.edge_nodes(elemNodes, iEdge, edgeNodes.data());
    const Edge edge = edge_from_edge_nodes(myMesh, edgeNodes[0], edgeNodes[1]);
    const stk::mesh::Entity refinedEdgeNode = nodeRefiner.get_edge_child_node(edge);
    elemEdgeChildNodes[iEdge] = refinedEdgeNode;
  }
}

void UniformEdgeMarker::fill_elements_modified_by_unrefinement(
    std::vector<stk::mesh::Entity> & parentElementsModifiedByUnrefinement,
    std::vector<stk::mesh::Entity> & childElementsToDeleteForUnrefinement) const
{
  parentElementsModifiedByUnrefinement.clear();
  childElementsToDeleteForUnrefinement.clear();
  // FIXME
}

ElementBasedEdgeMarker::ElementBasedEdgeMarker(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const std::string & elementMarkerFieldName)
    : myMesh(mesh), myRefinement(refinement)
{
  myElementMarkerField = static_cast<stk::mesh::Field<int> *>(mesh.mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, elementMarkerFieldName));
  ThrowRequireMsg(myElementMarkerField, "Element marker field \"" << elementMarkerFieldName << "\" not found.");
}

const std::string & ElementBasedEdgeMarker::get_marker_field_name() const
{
  ThrowAssert(myElementMarkerField);
  return myElementMarkerField->name();
}

const stk::mesh::Field<int> & ElementBasedEdgeMarker::get_marker_field() const
{
  ThrowAssert(myElementMarkerField);
  return *myElementMarkerField;
}

TransitionElementEdgeMarker::TransitionElementEdgeMarker(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const std::string & elementMarkerFieldName)
    : ElementBasedEdgeMarker(mesh, refinement, elementMarkerFieldName),
      myMesh(mesh),
      myRefinement(refinement)
{
}


static bool node_is_element_node(const StkMeshEntities & elementNodes, const stk::mesh::Entity node)
{
  return (std::find(elementNodes.begin(), elementNodes.end(), node) != elementNodes.end());
}

std::vector<stk::mesh::Entity> TransitionElementEdgeMarker::get_edge_nodes_of_transition_elements(
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & transitionElements) const
{
  std::vector<stk::mesh::Entity> transitionElementEdgeNodes;
  if (!transitionElements.empty())
  {
    const StkMeshEntities parentNodes {myMesh.begin_nodes(parentElem), myMesh.end_nodes(parentElem)};
    transitionElementEdgeNodes.reserve(2*transitionElements.size());
    for (auto && transitionElement : transitionElements)
      for (auto && elementNode : StkMeshEntities{myMesh.begin_nodes(transitionElement), myMesh.end_nodes(transitionElement)})
        if (!node_is_element_node(parentNodes, elementNode))
          transitionElementEdgeNodes.push_back(elementNode);
    stk::util::sort_and_unique(transitionElementEdgeNodes);
  }
  return transitionElementEdgeNodes;
}

bool TransitionElementEdgeMarker::is_transition(const stk::mesh::Entity elem) const
{
  // O(log(N)) search in get_parent().  If this becomes a hotspot, then maybe there should be a transition part.
  const stk::mesh::Bucket & bucket = myMesh.bucket(elem);
  if (!myRefinement.is_parent(bucket) && myRefinement.is_child(bucket))
  {
    auto parent = myRefinement.get_parent(elem);
    //If this is a leaf child but doesn't have its parent on this proc, it must be a
    //newly childless element that hasn't been returned to its originating proc yet.
    //Thus, it cannot be a transition element since a transition element couldn't have been
    //a parent in the first place
    if(!myMesh.is_valid(parent)) return false;
    return myRefinement.is_this_parent_element_partially_refined(myRefinement.get_parent(elem));
  }
  return false;
}

std::vector<Edge>
TransitionElementEdgeMarker::get_parent_edges_for_given_refined_edge_nodes(const std::vector<stk::mesh::Entity> & refinedEdgeNodes) const
{
  std::vector<Edge> edges;
  edges.reserve(refinedEdgeNodes.size());
  for (auto && refinedEdgeNode : refinedEdgeNodes)
  {
    const std::array<stk::mesh::Entity,2> edgeParentNodes = myRefinement.get_edge_parent_nodes(refinedEdgeNode);
    edges.push_back(edge_from_edge_nodes(myMesh, edgeParentNodes[0], edgeParentNodes[1]));
  }
  return edges;
}

void TransitionElementEdgeMarker::fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(
    const stk::mesh::Entity parentElem,
    const stk::topology & parentElemTopology,
    const unsigned parentElemNumEdges,
    const std::vector<stk::mesh::Entity> & childTransitionElements,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  elemEdgeChildNodes.clear();
  elemEdgeChildNodes.reserve(parentElemNumEdges);

  const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_edge_nodes_of_transition_elements(parentElem, childTransitionElements);
  const std::vector<Edge> refinedEdges = get_parent_edges_for_given_refined_edge_nodes(refinedEdgeNodes);

  const stk::mesh::Entity * entityNodes = myMesh.begin_nodes(parentElem);
  for (unsigned iEdge = 0; iEdge < parentElemNumEdges; ++iEdge)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(parentElemTopology, iEdge);
    const Edge edge = edge_from_edge_nodes(myMesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
    const auto iter = std::find(refinedEdges.begin(), refinedEdges.end(), edge);
    if (iter == refinedEdges.end())
      elemEdgeChildNodes.emplace_back(stk::mesh::Entity::InvalidEntity);
    else
      elemEdgeChildNodes.push_back(refinedEdgeNodes[std::distance(refinedEdges.begin(), iter)]);
  }
}

void TransitionElementEdgeMarker::fill_existing_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, const unsigned elemNumEdges, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  if (myRefinement.is_parent(elem))
  {
    ThrowAssert(myRefinement.is_this_parent_element_partially_refined(elem));
    const std::vector<stk::mesh::Entity> transitionElements = myRefinement.get_children(elem);
    fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(elem, elemTopology, elemNumEdges, transitionElements, elemEdgeChildNodes);
  }
  else
  {
    elemEdgeChildNodes.assign(elemNumEdges, stk::mesh::Entity());
  }
}

void TransitionElementEdgeMarker::fill_element_refined_edge_nodes(const NodeRefiner & nodeRefiner, const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const unsigned elemNumEdges = elemTopology.num_edges();

  fill_existing_element_refined_edge_nodes(elem, elemTopology, elemNumEdges, elemEdgeChildNodes);
  const stk::mesh::Entity * entityNodes = myMesh.begin_nodes(elem);

  bool haveNewNode = false;

  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
  {
    if (elemEdgeChildNodes[iEdge] == stk::mesh::Entity::InvalidEntity)
    {
      const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, iEdge);
      const Edge edge = edge_from_edge_nodes(myMesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
      const stk::mesh::Entity refinedEdgeNode = nodeRefiner.get_edge_child_node(edge);
      if (refinedEdgeNode != stk::mesh::Entity::InvalidEntity)
        haveNewNode = true;
      elemEdgeChildNodes[iEdge] = refinedEdgeNode;
    }
  }

  if (!haveNewNode) // all refined nodes already existed
    for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
      elemEdgeChildNodes[iEdge] = stk::mesh::Entity::InvalidEntity;
}

bool TransitionElementEdgeMarker::is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const
{
  if (myRefinement.is_parent(elem))
    return myRefinement.is_this_parent_element_partially_refined(elem);
  return !is_transition(elem);
}

static bool is_any_element_marked(FieldRef elementMarkerField, const std::vector<stk::mesh::Entity> & elems)
{
  for(const auto & elem : elems)
  {
    auto * elemMarker = field_data<int>(elementMarkerField, elem);
    if (elemMarker && *elemMarker == Refinement::REFINE)
      return true;
  }
  return false;
}

static bool is_element_marked_for_unrefinement(FieldRef elementMarkerField, const stk::mesh::Entity elem)
{
  auto * elemMarker = field_data<int>(elementMarkerField, elem);
  if (elemMarker && *elemMarker == Refinement::COARSEN)
    return true;
  return false;
}

void TransitionElementEdgeMarker::mark_unrefined_edges_of_partially_refined_parent_element(const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childTransitionElements,
    std::vector<Edge> & elemEdgesWorkspace,
    NodeRefiner & nodeRefiner,
    bool & wasAnyEdgeMarked) const
{
  const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_edge_nodes_of_transition_elements(parentElem, childTransitionElements);
  const std::vector<Edge> refinedEdges = get_parent_edges_for_given_refined_edge_nodes(refinedEdgeNodes);
  fill_entity_edges(myMesh, parentElem, elemEdgesWorkspace);
  for (auto && elemEdge : elemEdgesWorkspace)
    if (std::find(refinedEdges.begin(), refinedEdges.end(), elemEdge) == refinedEdges.end())
      wasAnyEdgeMarked |= nodeRefiner.mark_edge_for_refinement(elemEdge);
}

void TransitionElementEdgeMarker::locally_mark_edges_of_partially_refined_parent_elements_with_marked_children(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Field<int> & markerField = get_marker_field();
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & myRefinement.parent_part();
  std::vector<Edge> elemEdges;
  std::vector<stk::mesh::Entity> childTransitionElements;
  bool wasAnyEdgeMarked = false;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & parentElem : *bucketPtr)
    {
      if (myRefinement.is_this_parent_element_partially_refined(parentElem))
      {
        myRefinement.fill_children(parentElem, childTransitionElements);
        if (is_any_element_marked(markerField, childTransitionElements))
        {
          mark_unrefined_edges_of_partially_refined_parent_element(parentElem, childTransitionElements, elemEdges, nodeRefiner, wasAnyEdgeMarked);
        }
      }
    }
  }
}

static bool edge_is_parent_edge(const StkMeshEntities & parentNodes, const stk::mesh::Entity * edgeNodes)
{
  return node_is_element_node(parentNodes, edgeNodes[0]) &&
         node_is_element_node(parentNodes, edgeNodes[1]);
}

bool does_child_element_have_marked_edge_that_is_not_parent_edge(const stk::mesh::BulkData & mesh, const StkMeshEntities & parentNodes, const stk::mesh::Entity childElem, const NodeRefiner & nodeRefiner)
{
  const stk::topology childTopology = mesh.bucket(childElem).topology();
  const unsigned childNumEdges = childTopology.num_edges();
  const stk::mesh::Entity * childElemNodes = mesh.begin_nodes(childElem);
  std::array<stk::mesh::Entity,3> childEdgeNodes;

  for (unsigned iEdge = 0; iEdge < childNumEdges; ++iEdge)
  {
    childTopology.edge_nodes(childElemNodes, iEdge, childEdgeNodes.data());
    if (!edge_is_parent_edge(parentNodes, childEdgeNodes.data()))
    {
      const Edge edge = edge_from_edge_nodes(mesh, childEdgeNodes[0], childEdgeNodes[1]);
      if (nodeRefiner.is_edge_marked_for_refinement(edge))
        return true;
    }
  }
  return false;
}

bool does_any_child_element_have_marked_edge_that_is_not_parent_edge(const stk::mesh::BulkData & mesh, const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childElems, const NodeRefiner & nodeRefiner)
{
  const StkMeshEntities parentNodes {mesh.begin_nodes(parentElem), mesh.end_nodes(parentElem)};
  for(const auto & childElem : childElems)
    if(does_child_element_have_marked_edge_that_is_not_parent_edge(mesh, parentNodes, childElem, nodeRefiner))
      return true;
  return false;
}

void TransitionElementEdgeMarker::locally_mark_edges_of_partially_refined_parent_elements_to_satisfy_template(NodeRefiner & nodeRefiner, bool & wasAnyEdgeMarked) const
{
  const stk::mesh::Field<int> & markerField = get_marker_field();
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & myRefinement.parent_part();
  std::vector<Edge> elemEdges;
  std::vector<stk::mesh::Entity> childTransitionElements;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & parentElem : *bucketPtr)
    {
      if (myRefinement.is_this_parent_element_partially_refined(parentElem))
      {
        myRefinement.fill_children(parentElem, childTransitionElements);
        if (does_any_child_element_have_marked_edge_that_is_not_parent_edge(myMesh, parentElem, childTransitionElements, nodeRefiner))
          mark_unrefined_edges_of_partially_refined_parent_element(parentElem, childTransitionElements, elemEdges, nodeRefiner, wasAnyEdgeMarked);
      }
    }
  }
}

void TransitionElementEdgeMarker::locally_mark_edges_of_marked_non_transition_elements(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Field<int> & markerField = get_marker_field();
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & !myRefinement.parent_part();
  std::vector<Edge> edgesToRefineForElement;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const stk::mesh::Bucket & bucket = *bucketPtr;
    auto * elemMarker = field_data<int>(markerField, bucket);
    const unsigned markerFieldLength = stk::mesh::field_scalars_per_entity(markerField, bucket);

    for(size_t i=0; i<bucketPtr->size(); ++i)
    {
      const stk::mesh::Entity elem = bucket[i];

      if (elemMarker[i*markerFieldLength] == Refinement::REFINE && !is_transition(elem))
      {
        fill_entity_edges(myMesh, elem, edgesToRefineForElement);
        for (auto && edgeToRefineForElement : edgesToRefineForElement)
          nodeRefiner.mark_edge_for_refinement(edgeToRefineForElement);
      }
    }
  }
}

void TransitionElementEdgeMarker::mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const
{
//     Template Refinement algorithm
//     1. Mark all edges of marked, non-transition, leaf elements - one iter only
//     2. Mark all edges of partially refined parent element that are not already refined if any child element is marked - one iter only
//     3. Mark all edges of partially refined parent element that are not already refined if any edge of child element is marked that is not an edge of the parent - multiple iters

  nodeRefiner.clear_edges_to_refine();

  locally_mark_edges_of_marked_non_transition_elements(nodeRefiner); // stage 1
  locally_mark_edges_of_partially_refined_parent_elements_with_marked_children(nodeRefiner); // stage 2
  nodeRefiner.sync_shared_edges_from_other_procs_to_refine(myMesh);

  bool done = false;
  while (!done)
  {
    bool madeAnyChangeLocally = false;
    locally_mark_edges_of_partially_refined_parent_elements_to_satisfy_template(nodeRefiner, madeAnyChangeLocally); // stage 3
    done = stk::is_true_on_all_procs(myMesh.parallel(), !madeAnyChangeLocally);
    if (!done)
      nodeRefiner.sync_shared_edges_from_other_procs_to_refine(myMesh);
  }
}

bool TransitionElementEdgeMarker::are_all_children_leaves_and_marked_for_unrefinement(const std::vector<stk::mesh::Entity> & childElements) const
{
  const FieldRef markerField = get_marker_field();

  for (auto && childElem : childElements)
  {
    if (!myMesh.is_valid(childElem) || myRefinement.is_parent(childElem) || !is_element_marked_for_unrefinement(markerField, childElem))
      return false;
  }

  return true;
}

bool TransitionElementEdgeMarker::child_elements_are_all_leaves_and_are_transition_elements_or_marked_for_unrefinement(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childElements) const
{
  if (are_all_children_leaves_and_marked_for_unrefinement(childElements))
    return true;

  const bool isPartiallyRefined = myRefinement.get_num_children_when_fully_refined(parentElem) != childElements.size();
  if (isPartiallyRefined)
    return true;

  return false;
}

bool TransitionElementEdgeMarker::can_edge_node_be_unrefined_based_on_locally_owned_elements(const stk::mesh::Entity refinedNode,
    std::vector<stk::mesh::Entity> & workParentEdgeElements,
    std::vector<stk::mesh::Entity> & workChildElements) const
{
  // Fast, but only partial check.  If any element using this node is not a leaf child, then the node must be kept.
  for (auto && elem : StkMeshEntities{myMesh.begin_elements(refinedNode), myMesh.end_elements(refinedNode)})
    if (myRefinement.is_parent(elem))
      return false;

  const auto edgeNodeParents = myRefinement.get_edge_parent_nodes(refinedNode);
  stk::mesh::get_entities_through_relations(myMesh, {edgeNodeParents[0], edgeNodeParents[1]}, stk::topology::ELEMENT_RANK, workParentEdgeElements);

  for (auto && parentEdgeElem : workParentEdgeElements)
  {
    // Only examines locally owned elements so that any return value of false is sure to be correct.
    // Subsequent communication is then needed to make sure the final determination is parallel consistent.
    if (myMesh.bucket(parentEdgeElem).owned())
    {
      myRefinement.fill_children(parentEdgeElem, workChildElements);
      if (!child_elements_are_all_leaves_and_are_transition_elements_or_marked_for_unrefinement(parentEdgeElem, workChildElements))
        return false;
    }
  }

  return true;
}

static void communicate_to_get_sorted_edges_nodes_that_will_be_removed_by_unrefinement(const stk::mesh::BulkData & mesh,
  const std::vector<stk::mesh::Entity> & sharedEdgeNodesThatMustBeKept,
  std::vector<stk::mesh::Entity> & ownedOrSharedEdgeNodesThatMayBeUnrefined)
{
  stk::CommSparse commSparse(mesh.parallel());

  pack_entities_for_sharing_procs(mesh, sharedEdgeNodesThatMustBeKept, commSparse);
  std::vector<stk::mesh::Entity> sharedEdgeNodesThatMustBeKeptFromOtherProcs;
  unpack_shared_entities(mesh, sharedEdgeNodesThatMustBeKeptFromOtherProcs, commSparse);
  stk::util::remove_intersection_from_first(ownedOrSharedEdgeNodesThatMayBeUnrefined, sharedEdgeNodesThatMustBeKeptFromOtherProcs);
}

std::vector<stk::mesh::Entity> TransitionElementEdgeMarker::get_sorted_edge_nodes_that_will_be_removed_by_unrefinement() const
{
  std::vector<stk::mesh::Entity> ownedOrSharedEdgeNodesThatMayBeUnrefined;
  std::vector<stk::mesh::Entity> sharedEdgeNodesThatMustBeKept;
  std::vector<stk::mesh::Entity> workParentEdgeElements;
  std::vector<stk::mesh::Entity> workChildElements;

  const stk::mesh::Selector refinedEdgeNodeSelector = myRefinement.refined_edge_node_part() & (myMesh.mesh_meta_data().locally_owned_part() | myMesh.mesh_meta_data().globally_shared_part());

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::NODE_RANK, refinedEdgeNodeSelector))
  {
    for (auto && refinedNode : *bucketPtr)
    {
      if (can_edge_node_be_unrefined_based_on_locally_owned_elements(refinedNode, workParentEdgeElements, workChildElements))
      {
        ownedOrSharedEdgeNodesThatMayBeUnrefined.push_back(refinedNode);
      }
      else if (bucketPtr->shared())
      {
        sharedEdgeNodesThatMustBeKept.push_back(refinedNode);
      }
    }
  }

  communicate_to_get_sorted_edges_nodes_that_will_be_removed_by_unrefinement(myMesh, sharedEdgeNodesThatMustBeKept, ownedOrSharedEdgeNodesThatMayBeUnrefined);
  return ownedOrSharedEdgeNodesThatMayBeUnrefined;
}

static bool is_any_in_second_vec_in_first_sorted_vec(const std::vector<stk::mesh::Entity> & sortedFirstVec, const std::vector<stk::mesh::Entity> & secondVec)
{
  for (auto && entity : secondVec)
  {
    auto it = std::lower_bound(sortedFirstVec.begin(), sortedFirstVec.end(), entity);
    if (it != sortedFirstVec.end() && *it == entity)
      return true;
  }
  return false;
}

static bool is_any_child_invalid_or_a_parent(const stk::mesh::BulkData & mesh, const Refinement & refinement, const std::vector<stk::mesh::Entity> & childElements)
{
  for (auto && childElem : childElements)
    if (!mesh.is_valid(childElem) || refinement.is_parent(childElem))
      return true;
  return false;
}

bool TransitionElementEdgeMarker::is_parent_element_modified_by_unrefinement(const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElements,
    const std::vector<stk::mesh::Entity> & sortedEdgeNodesThatWillBeUnrefined,
    std::vector<stk::mesh::Entity> & workElemEdgeChildNodes) const
{
  if (is_any_child_invalid_or_a_parent(myMesh, myRefinement, childElements))
    return false;

  // If all edge nodes are going to unrefined, then this parent element will have no children after unrefinement
  // If only some of the edge nodes are going to unrefined, then the parent will still have children, but the current elements must still be deleted
  const stk::topology elemTopology = myMesh.bucket(parentElem).topology();
  const unsigned elemNumEdges = elemTopology.num_edges();
  fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(parentElem, elemTopology, elemNumEdges, childElements, workElemEdgeChildNodes);

  const bool isParentModifiedByUnrefinement = is_any_in_second_vec_in_first_sorted_vec(sortedEdgeNodesThatWillBeUnrefined, workElemEdgeChildNodes);

  return isParentModifiedByUnrefinement;
}

void TransitionElementEdgeMarker::fill_elements_modified_by_unrefinement(std::vector<stk::mesh::Entity> & parentElementsModifiedByUnrefinement,
    std::vector<stk::mesh::Entity> & childElementsToDeleteForUnrefinement) const
{
  const stk::mesh::Field<int> & markerField = get_marker_field();
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & myRefinement.parent_part();
  std::vector<stk::mesh::Entity> childElements;
  std::vector<stk::mesh::Entity> workElemEdgeChildNodes;

  parentElementsModifiedByUnrefinement.clear();
  childElementsToDeleteForUnrefinement.clear();

  const std::vector<stk::mesh::Entity> sortedEdgeNodesThatWillBeUnrefined = get_sorted_edge_nodes_that_will_be_removed_by_unrefinement();

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & parentElem : *bucketPtr)
    {
      myRefinement.fill_children(parentElem, childElements);
      if (is_parent_element_modified_by_unrefinement(parentElem, childElements, sortedEdgeNodesThatWillBeUnrefined, workElemEdgeChildNodes))
      {
        childElementsToDeleteForUnrefinement.insert(childElementsToDeleteForUnrefinement.end(), childElements.begin(), childElements.end());
        parentElementsModifiedByUnrefinement.push_back(parentElem);
      }
    }
  }
}

bool TransitionElementEdgeMarker::locally_have_elements_to_unrefine() const
{
  const stk::mesh::Field<int> & markerField = get_marker_field();
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & myRefinement.parent_part();
  std::vector<stk::mesh::Entity> childElements;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & parentElem : *bucketPtr)
    {
      myRefinement.fill_children(parentElem, childElements);
      const bool isPartiallyRefined = myRefinement.get_num_children_when_fully_refined(parentElem) != childElements.size();
      if (!isPartiallyRefined && are_all_children_leaves_and_marked_for_unrefinement(childElements))
        return true;
    }
  }

  return false;
}

}


