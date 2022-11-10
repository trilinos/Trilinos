
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
  NodeRefiner() {}
  NodeRefiner ( const NodeRefiner & ) = delete;
  NodeRefiner & operator= ( const NodeRefiner & ) = delete;
  typedef std::unordered_map<Edge,stk::mesh::Entity> RefinedEdgeMap;

  void find_edges_to_refine(const stk::mesh::BulkData & mesh, const EdgeMarkerInterface & edgeMarker);
  void mark_edge_for_refinement(const Edge & edge);
  void create_refined_edge_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedEdgeNodeParts, FieldRef refinedEdgeNodeParentIdsField);
  size_t get_num_edges_to_refine() const;
  stk::mesh::Entity get_edge_child_node(const Edge edge) const;
  void prolong_refined_edge_nodes(const stk::mesh::BulkData & mesh) const;
private:
   void assign_refined_edge_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedEdgeNodeParentIdsField) const;
  RefinedEdgeMap myRefinedEdgesToChildNodes;
};

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODEREFINER_HPP_ */
