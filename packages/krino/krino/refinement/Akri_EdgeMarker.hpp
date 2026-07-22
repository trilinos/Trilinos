#ifndef KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKER_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKER_HPP_

#include <string>
#include <vector>
#include <set>

#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_EdgeMarkerUtils.hpp>

namespace krino {

class Refinement;
struct Edge;
class NodeRefiner;

class EdgeMarkerInterface
{
public:
  virtual ~EdgeMarkerInterface() {}
  virtual void mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const = 0;
  virtual void mark_entities_to_be_unrefined(NodeRefiner & nodeRefiner) const = 0;
  virtual bool is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const = 0;
  virtual void fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const = 0;
  virtual std::vector<stk::mesh::Entity> get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & nodeRefiner) const = 0;
  virtual bool locally_have_elements_to_unrefine() const = 0;
};

class EdgeBasedEdgeMarker : public EdgeMarkerInterface
{
public:
  EdgeBasedEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement);
  virtual ~EdgeBasedEdgeMarker() {}

  void clear_marked_edges();
  void mark_edge_for_refinement(const Edge & edge);

protected:
  virtual void mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual void mark_entities_to_be_unrefined(NodeRefiner & nodeRefiner) const override;
  virtual bool is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const override;
  virtual void fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual std::vector<stk::mesh::Entity> get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & nodeRefiner) const override;
  virtual bool locally_have_elements_to_unrefine() const override;
  bool is_edge_marked_for_refinement(const Edge & edge) const;
private:
  void locally_mark_selected_edges_of_non_parent_elements(NodeRefiner & nodeRefiner) const;
  void fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
  std::set<Edge> myLocallyMarkedEdges;
};

class UniformRefinementEdgeMarker : public EdgeMarkerInterface
{
public:
  UniformRefinementEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement);
  virtual ~UniformRefinementEdgeMarker() {}
protected:
  virtual void mark_entities_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual void mark_entities_to_be_unrefined(NodeRefiner & /*nodeRefiner*/) const override {}
  virtual bool is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const override;
  virtual void fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual std::vector<stk::mesh::Entity> get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & /*nodeRefiner*/) const override {std::vector<stk::mesh::Entity> tmp; return tmp;}
  virtual bool locally_have_elements_to_unrefine() const override { return false; }
private:
  void locally_mark_all_entities_of_non_parent_elements(NodeRefiner & nodeRefiner) const;
  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
};

}

#endif /* KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKER_HPP_ */
