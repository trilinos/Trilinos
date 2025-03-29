
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#include <unordered_map>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include "Akri_Edge.hpp"
#include "Akri_QuadFace.hpp"

namespace krino {

struct Edge;
class FieldRef;
class EdgeMarkerInterface;

class NodeRefiner
{
public:
  NodeRefiner(const bool force64Bit, const bool assert32Bit) : myForce64Bit(force64Bit), myAssert32Bit(assert32Bit) {}
  NodeRefiner ( const NodeRefiner & ) = delete;
  NodeRefiner & operator= ( const NodeRefiner & ) = delete;
  using EdgeToNodeMap = std::map<Edge,stk::mesh::Entity>;
  using RefinedEdgeMap = EdgeToNodeMap;
  using RefinedQuadFaceMap = std::map<QuadFace,stk::mesh::Entity>;
  using EdgeMidNodeMap = EdgeToNodeMap;

  void clear_entities_to_refine();
  void sync_shared_edges_from_other_procs_to_refine(const stk::mesh::BulkData & mesh);
  bool mark_edge_for_refinement(const Edge & edge);
  bool mark_edge_as_needing_midnode(const Edge & edge);
  bool is_edge_marked_for_refinement(const Edge & edge) const;
  bool mark_quad_face_for_refinement(const QuadFace & quadFace);
  void create_refined_edge_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedEdgeNodeParts, const bool useExistingEdgeMidNodes);
  void create_edge_midnodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedEdgeMidNodeParts);
  void create_refined_element_centroid_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedElemCentroidNodeParts);
  void assign_refined_quad_face_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedQuadFaceNodeParentIdsField) const;
  void create_refined_quad_face_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedQuadFaceNodeParts);
  size_t get_num_edges_to_refine() const;
  bool locally_have_edges_to_refine() const { return !myRefinedEdgesToChildNodes.empty(); }
  stk::mesh::Entity get_edge_child_node(const Edge edge) const;
  stk::mesh::Entity get_edge_midnode(const Edge edge) const;
  stk::mesh::Entity get_element_centroid_child_node(const stk::mesh::Entity elem) const;
  stk::mesh::Entity get_element_child_face_node(const QuadFace quadFace) const;
  std::vector<stk::mesh::Entity> get_element_child_face_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem) const;
  void prolong_refined_nodes_and_edge_midnodes(const stk::mesh::BulkData & mesh) const;
  void mark_element_with_child_centroid_node_for_refinement(const stk::mesh::Entity & elem);
  void set_sorted_edge_nodes_that_will_be_removed_by_unrefinement(const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedEdgeNodesToBeRemovedByUnrefinement)
    { mySortedOwnedOrSharedEdgeNodesToBeRemovedByUnrefinement = sortedOwnedOrSharedEdgeNodesToBeRemovedByUnrefinement; }
  const std::vector<stk::mesh::Entity> & get_sorted_edge_nodes_that_will_be_removed_by_unrefinement() const { return mySortedOwnedOrSharedEdgeNodesToBeRemovedByUnrefinement; }
  void assign_refined_edge_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedEdgeNodeParentIdsField) const;
private:
  RefinedEdgeMap myRefinedEdgesToChildNodes;
  EdgeMidNodeMap myRefinedEdgesToMidNodes;
  RefinedQuadFaceMap myRefinedQuadFacesToChildNodes;
  std::vector<std::pair<stk::mesh::Entity,stk::mesh::Entity>> myRefinedElementsToChildNodes;
  std::vector<stk::mesh::Entity> mySortedOwnedOrSharedEdgeNodesToBeRemovedByUnrefinement;
  bool myForce64Bit;
  bool myAssert32Bit;
};

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_ */
