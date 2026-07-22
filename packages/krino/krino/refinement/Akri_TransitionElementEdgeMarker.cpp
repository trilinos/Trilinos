/*
 * Akri_TransitionElementEdgeMarker.cpp
 *
 *  Created on: Oct 28, 2022
 *      Author: drnoble
 */
#include "Akri_TransitionElementEdgeMarker.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/RemoveIntersection.hpp>
#include "Akri_Edge.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_Refinement.hpp"

namespace krino {

ElementBasedEdgeMarker::ElementBasedEdgeMarker(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const std::string & elementMarkerFieldName)
    : myMesh(mesh), myRefinement(refinement)
{
  myElementMarkerField = static_cast<stk::mesh::Field<int> *>(mesh.mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, elementMarkerFieldName));
  STK_ThrowRequireMsg(myElementMarkerField, "Element marker field \"" << elementMarkerFieldName << "\" not found.");
}

const std::string & ElementBasedEdgeMarker::get_marker_field_name() const
{
  STK_ThrowAssert(myElementMarkerField);
  return myElementMarkerField->name();
}

void ElementBasedEdgeMarker::set_marker_field(stk::mesh::Field<int> * field)
{
  STK_ThrowAssert(field);
  myElementMarkerField = field;
}

const stk::mesh::Field<int> & ElementBasedEdgeMarker::get_marker_field_and_sync_to_host() const
{
  STK_ThrowAssert(myElementMarkerField);
  myElementMarkerField->sync_to_host();
  return *myElementMarkerField;
}

TransitionElementEdgeMarker::TransitionElementEdgeMarker(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const std::string & elementMarkerFieldName)
    : ElementBasedEdgeMarker(mesh, refinement, elementMarkerFieldName)
{
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

void TransitionElementEdgeMarker::fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const bool isElemCandidateForRefinement = !myRefinement.is_parent(elem) || myRefinement.is_this_parent_element_partially_refined(elem);
  EdgeMarkerUtils::fill_refinement_caseIds_and_refined_edge_nodes_if_changed(myMesh, myRefinement, nodeRefiner, elem, isElemCandidateForRefinement, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
}

void TransitionElementEdgeMarker::fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const bool isElemCandidateForUnrefinement = myRefinement.is_parent(elem);
  EdgeMarkerUtils::fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(myMesh, myRefinement, nodeRefiner, elem, isElemCandidateForUnrefinement, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
}

void TransitionElementEdgeMarker::fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  if (doingRefinement)
    fill_refinement_caseIds_and_refined_edge_nodes_if_changed(nodeRefiner, elem, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
  else
    fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(nodeRefiner, elem, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
}

bool TransitionElementEdgeMarker::is_element_a_candidate_for_unrefinement(const stk::mesh::Entity elem) const
{
  return EdgeMarkerUtils::is_element_a_parent_but_not_a_grandparent(myMesh, myRefinement, elem);
}

bool TransitionElementEdgeMarker::is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const
{
  if (myRefinement.is_parent(elem))
    return myRefinement.is_this_parent_element_partially_refined(elem);
  return !is_transition(elem);
}

bool TransitionElementEdgeMarker::is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const
{
  if (doingRefinement)
    return is_element_a_candidate_for_refinement(elem);
  return is_element_a_candidate_for_unrefinement(elem);
}

static bool is_any_element_marked(FieldRef elementMarkerField, const std::vector<stk::mesh::Entity> & elems)
{
  for(const auto & elem : elems)
  {
    auto * elemMarker = field_data<int>(elementMarkerField, elem);
    if (elemMarker && *elemMarker == static_cast<int>(Refinement::RefinementMarker::REFINE))
      return true;
  }
  return false;
}

static bool is_element_marked_for_unrefinement(FieldRef elementMarkerField, const stk::mesh::Entity elem)
{
  auto * elemMarker = field_data<int>(elementMarkerField, elem);
  if (elemMarker && *elemMarker == static_cast<int>(Refinement::RefinementMarker::COARSEN))
    return true;
  return false;
}

void TransitionElementEdgeMarker::mark_unrefined_edges_of_partially_refined_parent_element(const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childTransitionElements,
    std::vector<Edge> & elemEdgesWorkspace,
    NodeRefiner & nodeRefiner,
    bool & wasAnyEdgeMarked) const
{
  const std::vector<stk::mesh::Entity> refinedEdgeNodes = EdgeMarkerUtils::get_child_nodes_that_are_not_parent_nodes(myMesh, parentElem, childTransitionElements);
  const std::vector<Edge> refinedEdges = get_parent_edges_for_given_refined_edge_nodes(refinedEdgeNodes);
  fill_entity_edges(myMesh, parentElem, elemEdgesWorkspace);
  for (auto && elemEdge : elemEdgesWorkspace)
    if (std::find(refinedEdges.begin(), refinedEdges.end(), elemEdge) == refinedEdges.end())
      wasAnyEdgeMarked |= nodeRefiner.mark_edge_for_refinement(elemEdge);
}

void TransitionElementEdgeMarker::locally_mark_edges_of_partially_refined_parent_elements_with_marked_children(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Field<int> & markerField = get_marker_field_and_sync_to_host();
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
  return EdgeMarkerUtils::node_is_parent_node(parentNodes, edgeNodes[0]) &&
         EdgeMarkerUtils::node_is_parent_node(parentNodes, edgeNodes[1]);
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
  const stk::mesh::Field<int> & markerField = get_marker_field_and_sync_to_host();
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
  const stk::mesh::Field<int> & markerField = get_marker_field_and_sync_to_host();
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

      if (elemMarker[i*markerFieldLength] == static_cast<int>(Refinement::RefinementMarker::REFINE) && !is_transition(elem))
      {
        fill_entity_edges(myMesh, elem, edgesToRefineForElement);
        for (auto && edgeToRefineForElement : edgesToRefineForElement)
          nodeRefiner.mark_edge_for_refinement(edgeToRefineForElement);
      }
    }
  }
}

void TransitionElementEdgeMarker::mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const
{
//     Template Refinement algorithm
//     1. Mark all edges of marked, non-transition, leaf elements - one iter only
//     2. Mark all edges of partially refined parent element that are not already refined if any child element is marked - one iter only
//     3. Mark all edges of partially refined parent element that are not already refined if any edge of child element is marked that is not an edge of the parent - multiple iters

  nodeRefiner.clear_entities_to_refine();

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

void TransitionElementEdgeMarker::mark_entities_to_be_unrefined(NodeRefiner & nodeRefiner) const
{
  nodeRefiner.set_sorted_edge_nodes_that_will_be_removed_by_unrefinement(find_sorted_edge_nodes_that_will_be_removed_by_unrefinement());
}

bool TransitionElementEdgeMarker::are_all_children_leaves_and_marked_for_unrefinement(const std::vector<stk::mesh::Entity> & childElements) const
{
  const FieldRef markerField = get_marker_field_and_sync_to_host();

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
  stk::mesh::get_entities_through_relations(myMesh, stk::mesh::EntityVector{edgeNodeParents[0], edgeNodeParents[1]}, stk::topology::ELEMENT_RANK, workParentEdgeElements);

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
  std::vector<stk::mesh::Entity> sharedEdgeNodesThatMustBeKeptFromOtherProcs = unpack_entities_from_other_procs(mesh, commSparse);
  stk::util::remove_intersection_from_first(ownedOrSharedEdgeNodesThatMayBeUnrefined, sharedEdgeNodesThatMustBeKeptFromOtherProcs);
}

std::vector<stk::mesh::Entity> TransitionElementEdgeMarker::find_sorted_edge_nodes_that_will_be_removed_by_unrefinement() const
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

bool TransitionElementEdgeMarker::is_parent_element_modified_by_unrefinement(const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElements,
    const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement) const
{
  if (EdgeMarkerUtils::is_any_child_invalid_or_a_parent(myMesh, myRefinement, childElements))
    return false;

  const std::vector<stk::mesh::Entity> refinedEdgeNodes = EdgeMarkerUtils::get_child_nodes_that_are_not_parent_nodes(myMesh, myRefinement, parentElem);

  const bool isParentModifiedByUnrefinement = is_any_in_second_vec_in_first_sorted_vec(sortedOwnedOrSharedNodesToBeRemovedByUnrefinement, refinedEdgeNodes);

  return isParentModifiedByUnrefinement;
}

std::vector<stk::mesh::Entity> TransitionElementEdgeMarker::get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & nodeRefiner) const
{
  std::vector<stk::mesh::Entity> parentElementsModifiedByUnrefinement;
  const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement = nodeRefiner.get_sorted_edge_nodes_that_will_be_removed_by_unrefinement();

  if (!sortedOwnedOrSharedNodesToBeRemovedByUnrefinement.empty())
  {
    const stk::mesh::Field<int> & markerField = get_marker_field_and_sync_to_host();
    const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(markerField) & myRefinement.parent_part();
    std::vector<stk::mesh::Entity> childElements;

    for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
    {
      for(const auto & parentElem : *bucketPtr)
      {
        myRefinement.fill_children(parentElem, childElements);
        if (is_parent_element_modified_by_unrefinement(parentElem, childElements, sortedOwnedOrSharedNodesToBeRemovedByUnrefinement))
          parentElementsModifiedByUnrefinement.push_back(parentElem);
      }
    }
  }

  return parentElementsModifiedByUnrefinement;
}

bool TransitionElementEdgeMarker::locally_have_elements_to_unrefine() const
{
  const stk::mesh::Field<int> & markerField = get_marker_field_and_sync_to_host();
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


