
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_
#include <unordered_map>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include "Akri_Edge.hpp"

namespace krino {

class Edge;
class FieldRef;
class EdgeMarkerInterface;

class NodeRefiner
{
public:
  NodeRefiner(const bool force64Bit, const bool assert32Bit) : myForce64Bit(force64Bit), myAssert32Bit(assert32Bit) {}
  NodeRefiner ( const NodeRefiner & ) = delete;
  NodeRefiner & operator= ( const NodeRefiner & ) = delete;
  typedef std::unordered_map<Edge,stk::mesh::Entity> RefinedEdgeMap;

  void clear_edges_to_refine();
  void sync_shared_edges_from_other_procs_to_refine(const stk::mesh::BulkData & mesh);
  bool mark_edge_for_refinement(const Edge & edge);
  bool is_edge_marked_for_refinement(const Edge & edge) const;
  void create_refined_edge_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedEdgeNodeParts, FieldRef refinedEdgeNodeParentIdsField);
  size_t get_num_edges_to_refine() const;
  stk::mesh::Entity get_edge_child_node(const Edge edge) const;
  void prolong_refined_edge_nodes(const stk::mesh::BulkData & mesh) const;
  bool mark_already_refined_edge(const Edge & edge, const stk::mesh::Entity refinedEdgeNode);
private:
   void assign_refined_edge_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedEdgeNodeParentIdsField) const;
  RefinedEdgeMap myRefinedEdgesToChildNodes;
  bool myForce64Bit;
  bool myAssert32Bit;
};

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_ */
