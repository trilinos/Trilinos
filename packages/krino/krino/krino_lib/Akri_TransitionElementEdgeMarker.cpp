/*
 * Akri_TransitionElementEdgeMarker.cpp
 *
 *  Created on: Oct 28, 2022
 *      Author: drnoble
 */
#include "Akri_TransitionElementEdgeMarker.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "Akri_Edge.hpp"
#include "Akri_FieldRef.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_Refinement.hpp"

namespace krino {

  ElementBasedEdgeMarker::ElementBasedEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement, const std::string & elementMarkerFieldName)
: myMesh(mesh), myRefinement(refinement)
{
  myElementMarkerField = mesh.mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, elementMarkerFieldName);
  ThrowRequireMsg(myElementMarkerField.valid(), "Element marker field \"" << elementMarkerFieldName << "\" not found.");
}

TransitionElementEdgeMarker::TransitionElementEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement, const std::string & elementMarkerFieldName)
: ElementBasedEdgeMarker(mesh, refinement, elementMarkerFieldName), myMesh(mesh), myRefinement(refinement)
{
}

std::vector<stk::mesh::Entity> TransitionElementEdgeMarker::get_edge_nodes_of_transition_elements(const std::vector<stk::mesh::Entity> & transitionElements) const
{
  std::vector<stk::mesh::Entity> transitionElementEdgeNodes;
  if (!transitionElements.empty())
  {
    transitionElementEdgeNodes.reserve(2*transitionElements.size());
    for (auto && transitionElement : transitionElements)
      for (auto && elementNode : StkMeshEntities{myMesh.begin_nodes(transitionElement), myMesh.end_nodes(transitionElement)})
        if (myRefinement.is_refined_edge_node(elementNode))
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
    return myRefinement.is_a_partially_refined_parent_element(myRefinement.get_parent(elem));
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

  const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_edge_nodes_of_transition_elements(childTransitionElements);
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
    ThrowAssert(myRefinement.is_a_partially_refined_parent_element(elem));
    const std::vector<stk::mesh::Entity> transitionElements = myRefinement.get_children(elem);
    fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(elem, elemTopology, elemNumEdges, transitionElements, elemEdgeChildNodes);
  }
  else
  {
    elemEdgeChildNodes.assign(elemNumEdges, stk::mesh::Entity());
  }
}

void TransitionElementEdgeMarker::fill_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const unsigned elemNumEdges = elemTopology.num_edges();

  fill_existing_element_refined_edge_nodes(elem, elemTopology, elemNumEdges, elemEdgeChildNodes);
  const stk::mesh::Entity * entityNodes = myMesh.begin_nodes(elem);

  for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
  {
    if (elemEdgeChildNodes[iEdge] == stk::mesh::Entity::InvalidEntity)
    {
      const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, iEdge);
      const Edge edge = edge_from_edge_nodes(myMesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]);
      elemEdgeChildNodes[iEdge] = myRefinement.get_edge_child_node(edge);
    }
  }
}

bool TransitionElementEdgeMarker::is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const
{
  if (myRefinement.is_parent(elem))
    return myRefinement.is_a_partially_refined_parent_element(elem);
  return !is_transition(elem);
}

void TransitionElementEdgeMarker::fill_edges_to_refine_for_marked_element(const stk::mesh::Entity elem, std::vector<Edge> & edgesToRefineForElement) const
{
//     Template Refinement algorithm
//     * If the marked element is a non-transition leaf element, mark all edges of the element
//     * If the marked element is a transition element, mark all edges of the parent of the transition element that are not already refined

  edgesToRefineForElement.clear();

  std::vector<Edge> elemEdges;
  std::vector<stk::mesh::Entity> childTransitionElements;

  stk::mesh::Entity parentElem;
  bool isTransition = false;
  if (myRefinement.is_child(elem))
  {
    parentElem = myRefinement.get_parent(elem);
    isTransition = myRefinement.is_a_partially_refined_parent_element(parentElem);
  }
  if (isTransition)
  {
    myRefinement.fill_children(parentElem, childTransitionElements);
    const std::vector<stk::mesh::Entity> refinedEdgeNodes = get_edge_nodes_of_transition_elements(childTransitionElements);
    const std::vector<Edge> refinedEdges = get_parent_edges_for_given_refined_edge_nodes(refinedEdgeNodes);

    fill_entity_edges(myMesh, parentElem, elemEdges);
    for (auto && elemEdge : elemEdges)
      if (std::find(refinedEdges.begin(), refinedEdges.end(), elemEdge) == refinedEdges.end())
        edgesToRefineForElement.push_back(elemEdge);
  }
  else
  {
    fill_entity_edges(myMesh, elem, elemEdges);
    for (auto && elemEdge : elemEdges)
      edgesToRefineForElement.push_back(elemEdge);
  }
}

void ElementBasedEdgeMarker::mark_local_edges_to_be_refined() const
{
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & stk::mesh::selectField(myElementMarkerField.field()) & !myRefinement.parent_part();
  std::vector<Edge> edgesToRefineForElement;

  for(const auto & bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & elem : *bucketPtr)
    {
      auto * elemMarker = field_data<int>(myElementMarkerField, elem);
      if (elemMarker && *elemMarker == Refinement::REFINE)
      {
        fill_edges_to_refine_for_marked_element(elem, edgesToRefineForElement);

        for (auto && edgeToRefineForElement : edgesToRefineForElement)
          myRefinement.mark_edge_for_refinement(edgeToRefineForElement);
      }
    }
  }
}


}


