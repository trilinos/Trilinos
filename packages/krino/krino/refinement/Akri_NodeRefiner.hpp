
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#include <unordered_map>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <Akri_Edge.hpp>
#include <Akri_QuadFace.hpp>

namespace krino {

class FieldRef;
class EdgeMarkerInterface;

//class ChildNodesAndParents
//{
//public:
//  using ChildNodeAnd2Parents = std::pair<stk::mesh::Entity, std::array<stk::mesh::Entity,2>>;
//  using ChildNodeAndNParents = std::pair<stk::mesh::Entity, std::vector<stk::mesh::Entity>>;
//  std::vector<ChildNodeAnd2Parents> childNodesWith2Parents;
//  std::vector<ChildNodeAndNParents> childNodesWithNParents;
//};

//class ChildNodeAndParentsIterator {
//public:
//  ChildNodeAndParentsIterator(const size_t index) : myIndex(index) {}
//  Node operator*() const;
//  ChildNodeAndParentsIterator& operator++() { myIndex++; return *this; }
//  ChildNodeAndParentsIterator operator++(int) { ChildNodeAndParentsIterator tmp = *this; ++(*this); return tmp; }
//
//  friend bool operator==(const ChildNodeAndParentsIterator& a, const ChildNodeAndParentsIterator& b) { return a.myIndex == b.myIndex; }
//  friend bool operator!=(const ChildNodeAndParentsIterator& a, const ChildNodeAndParentsIterator& b) { return a.myIndex != b.myIndex; }
//private:
//  size_t myIndex;
//};

//class ChildNodesAndParents
//{
//public:
//  ChildNodesAndParents()
//    : myNumChildNodesWith2Parents(0)
//  {
//    myIndexForChildNodesWithMoreThan2Parents.push_back(0);
//  }
//
//  size_t size() const { return myNumChildNodesWith2Parents + num_child_nodes_with_more_than_2_parents(); }
//  void fill_child_and_parents(const size_t childNodeIndex, stk::mesh::Entity & childNode, std::vector<stk::mesh::Entity> & parentNodes) const
//  {
//    STK_ThrowAssert(childNodeIndex < size());
//    size_t childNodeAndParentsIndex;
//    size_t numParents;
//    if (childNodeIndex < myNumChildNodesWith2Parents)
//    {
//      childNodeAndParentsIndex = childNodeIndex*3;
//      numParents = 2;
//    }
//    else
//    {
//      const size_t moreThan2ParentsIndex = childNodeIndex - myNumChildNodesWith2Parents;
//      childNodeAndParentsIndex = myIndexForChildNodesWithMoreThan2Parents[moreThan2ParentsIndex];
//      numParents = myIndexForChildNodesWithMoreThan2Parents[moreThan2ParentsIndex+1] - childNodeAndParentsIndex - 1;
//    }
//    childNode = myChildNodesAndParents[childNodeAndParentsIndex];
//    const auto parentsBegin = myChildNodesAndParents.begin() + childNodeAndParentsIndex+1;
//    const auto parentsEnd = parentsBegin + numParents;
//    parentNodes.assign(parentsBegin, parentsEnd);
//  }
//private:
//  size_t num_child_nodes_with_more_than_2_parents() const { return myIndexForChildNodesWithMoreThan2Parents.size()-1; }
//  std::vector<stk::mesh::Entity> myChildNodesAndParents;
//  size_t myNumChildNodesWith2Parents;
//  std::vector<size_t> myIndexForChildNodesWithMoreThan2Parents;
//};

using UnaryForChildNodeAndParents = std::function<void(const stk::mesh::Entity childNode, const stk::mesh::Entity * parentNodeBeg, const stk::mesh::Entity * parentNodeEnd)>;

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
  size_t num_refined_nodes_and_edge_midnodes() const;
  void for_each_refined_node_and_edge_midnode(const stk::mesh::BulkData & mesh, const UnaryForChildNodeAndParents & unary) const;
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
