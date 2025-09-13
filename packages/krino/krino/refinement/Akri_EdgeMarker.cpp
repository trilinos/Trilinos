#include <Akri_EdgeMarker.hpp>
#include <Akri_MeshHelpers.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "Akri_Edge.hpp"
#include "Akri_Refinement.hpp"

namespace krino {

static bool element_topology_has_refined_centroid_node(const stk::topology elemTopology)
{
  return (elemTopology == stk::topology::QUAD_4_2D ||
      elemTopology == stk::topology::HEXAHEDRON_8);
}

static bool element_topology_has_refined_face_nodes(const stk::topology elemTopology)
{
  return (elemTopology == stk::topology::HEXAHEDRON_8);
}


UniformRefinementEdgeMarker::UniformRefinementEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement)
: myMesh(mesh), myRefinement(refinement)
{
}

void UniformRefinementEdgeMarker::locally_mark_all_entities_of_non_parent_elements(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & !myRefinement.parent_part();
  std::vector<Edge> edgesToRefineForElement;
  std::vector<QuadFace> facesToRefineForElement;

  for(const auto * bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const bool elemHasRefinedCentroidNode = element_topology_has_refined_centroid_node(bucketPtr->topology());
    const bool elemHasRefinedFaceNodes = element_topology_has_refined_face_nodes(bucketPtr->topology());
    for(const auto & elem : *bucketPtr)
    {
      fill_entity_edges(myMesh, elem, edgesToRefineForElement);
      for (auto && edgeToRefineForElement : edgesToRefineForElement)
        nodeRefiner.mark_edge_for_refinement(edgeToRefineForElement);
      if (elemHasRefinedCentroidNode)
        nodeRefiner.mark_element_with_child_centroid_node_for_refinement(elem);
      if (elemHasRefinedFaceNodes)
      {
        fill_entity_quad_faces(myMesh, elem, facesToRefineForElement);
        for (auto && faceToRefineForElement : facesToRefineForElement)
          nodeRefiner.mark_quad_face_for_refinement(faceToRefineForElement);
      }
    }
  }
}

void UniformRefinementEdgeMarker::mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const
{
  nodeRefiner.clear_entities_to_refine();
  locally_mark_all_entities_of_non_parent_elements(nodeRefiner);
  nodeRefiner.sync_shared_edges_from_other_procs_to_refine(myMesh);
}

bool UniformRefinementEdgeMarker::is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const
{
  return doingRefinement && !myRefinement.is_parent(elem);
}

void UniformRefinementEdgeMarker::fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  elementEdgeCaseIds.clear();

  if (doingRefinement)
  {
    const unsigned elemNumEdges = elemTopology.num_edges();
    const stk::mesh::Entity * elemNodes = myMesh.begin_nodes(elem);
    std::array<stk::mesh::Entity,3> edgeNodes;
    elemEdgeChildNodes.resize(elemNumEdges);

    int postRefineCaseId = 0;
    for (unsigned iEdge = 0; iEdge < elemNumEdges; ++iEdge)
    {
      elemTopology.edge_nodes(elemNodes, iEdge, edgeNodes.data());
      const Edge edge = edge_from_edge_nodes(myMesh, edgeNodes[0], edgeNodes[1]);
      const stk::mesh::Entity refinedEdgeNode = nodeRefiner.get_edge_child_node(edge);
      elemEdgeChildNodes[iEdge] = refinedEdgeNode;
      if (elemEdgeChildNodes[iEdge] != stk::mesh::Entity::InvalidEntity)
        postRefineCaseId += 1<<iEdge;
    }

    const int preRefineCaseId = 0;
    elementEdgeCaseIds.set(preRefineCaseId, postRefineCaseId);
  }
}

EdgeBasedEdgeMarker::EdgeBasedEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement)
: myMesh(mesh), myRefinement(refinement)
{
}

bool EdgeBasedEdgeMarker::is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const
{
  if (doingRefinement)
      return !myRefinement.is_parent(elem);
  // Not supporting unrefinement yet
  return false; //return EdgeMarkerUtils::is_element_a_parent_but_not_a_grandparent(myMesh, myRefinement, elem);
}

void EdgeBasedEdgeMarker::locally_mark_selected_edges_of_non_parent_elements(NodeRefiner & nodeRefiner) const
{
  const stk::mesh::Selector selector = myMesh.mesh_meta_data().locally_owned_part() & !myRefinement.parent_part();
  std::vector<Edge> elemEdges;

  for(const auto * bucketPtr : myMesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_entity_edges(myMesh, elem, elemEdges);
      for (auto && edge : elemEdges)
        if (is_edge_marked_for_refinement(edge))
          nodeRefiner.mark_edge_for_refinement(edge);
    }
  }
}

void EdgeBasedEdgeMarker::mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const
{
  nodeRefiner.clear_entities_to_refine();
  locally_mark_selected_edges_of_non_parent_elements(nodeRefiner);
  nodeRefiner.sync_shared_edges_from_other_procs_to_refine(myMesh);
}

bool EdgeBasedEdgeMarker::locally_have_elements_to_unrefine() const
{
  // Not supporting unrefinement yet
  return false;
}

void EdgeBasedEdgeMarker::clear_marked_edges()
{
  myLocallyMarkedEdges.clear();
}

void EdgeBasedEdgeMarker::mark_edge_for_refinement(const Edge & edge)
{
  myLocallyMarkedEdges.insert(edge);
}

bool EdgeBasedEdgeMarker::is_edge_marked_for_refinement(const Edge & edge) const
{
  return myLocallyMarkedEdges.find(edge) != myLocallyMarkedEdges.end();
}

void EdgeBasedEdgeMarker::mark_entities_to_be_unrefined(NodeRefiner & /*nodeRefiner*/) const
{
  // Not supporting unrefinement yet
}

std::vector<stk::mesh::Entity> EdgeBasedEdgeMarker::get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & /*nodeRefiner*/) const
{
  // Not supporting unrefinement yet
  std::vector<stk::mesh::Entity> tmp;
  return tmp;
}

void EdgeBasedEdgeMarker::fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  const bool isElemCandidateForRefinement = !myRefinement.is_parent(elem);
  EdgeMarkerUtils::fill_refinement_caseIds_and_refined_edge_nodes_if_changed(myMesh, myRefinement, nodeRefiner, elem, isElemCandidateForRefinement, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
}

void EdgeBasedEdgeMarker::fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const
{
  if (doingRefinement)
    fill_refinement_caseIds_and_refined_edge_nodes_if_changed(nodeRefiner, elem, elemTopology, elementEdgeCaseIds, elemEdgeChildNodes);
  // Not supporting unrefinement yet
}

}

