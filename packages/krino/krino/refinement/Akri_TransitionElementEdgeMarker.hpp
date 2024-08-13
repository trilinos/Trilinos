/*
 * Akri_TransitionElementEdgeMarker.hpp
 *
 *  Created on: Oct 28, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_
#include <string>
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace krino {

class Refinement;
struct Edge;
class NodeRefiner;

class ElementEdgeCaseIds
{
public:
  void set(const int preAdaptCaseId, const int postAdaptCaseId) { myPreAdaptCaseId=preAdaptCaseId; myPostAdaptCaseId=postAdaptCaseId; }
  bool has_changed() const { return myPreAdaptCaseId != myPostAdaptCaseId; }
  void clear() { myPreAdaptCaseId=0; myPostAdaptCaseId=0; }
  int pre_adapt_case_id() const { STK_ThrowAssert(has_changed()); return myPreAdaptCaseId; }
  int post_adapt_case_id() const { STK_ThrowAssert(has_changed()); return myPostAdaptCaseId; }
private:
  int myPreAdaptCaseId{0};
  int myPostAdaptCaseId{0};
};

class EdgeMarkerInterface
{
public:
  virtual ~EdgeMarkerInterface() {}
  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const = 0;
  virtual void mark_edges_to_be_unrefined(NodeRefiner & nodeRefiner) const = 0;
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

class UniformEdgeMarker : public EdgeMarkerInterface
{
public:
  UniformEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement);
  virtual ~UniformEdgeMarker() {}
protected:
  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual void mark_edges_to_be_unrefined(NodeRefiner & nodeRefiner) const override {}
  virtual bool is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const override;
  virtual void fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual std::vector<stk::mesh::Entity> get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & nodeRefiner) const override {std::vector<stk::mesh::Entity> tmp; return tmp;};
  virtual bool locally_have_elements_to_unrefine() const override { return false; }
private:
  void locally_mark_edges_of_non_parent_elements(NodeRefiner & nodeRefiner) const;
  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
};

class ElementBasedEdgeMarker : public EdgeMarkerInterface
{
public:
  ElementBasedEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement, const std::string & elementMarkerFieldName);
  virtual ~ElementBasedEdgeMarker() {}
  const std::string & get_marker_field_name() const;
protected:
  const stk::mesh::Field<int> & get_marker_field_and_sync_to_host() const;
private:
  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
  stk::mesh::Field<int> * myElementMarkerField{nullptr};
};

class TransitionElementEdgeMarker : public ElementBasedEdgeMarker
{
public:
  TransitionElementEdgeMarker(const stk::mesh::BulkData & mesh,
      Refinement & refinement,
      const std::string & elementMarkerFieldName);
      
  TransitionElementEdgeMarker ( const TransitionElementEdgeMarker & ) = delete;
  TransitionElementEdgeMarker & operator= ( const TransitionElementEdgeMarker & ) = delete;
  virtual ~TransitionElementEdgeMarker() {}

  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual void mark_edges_to_be_unrefined(NodeRefiner & nodeRefiner) const override;
  virtual bool is_element_a_candidate_for_adaptation(const stk::mesh::Entity elem, const bool doingRefinement) const override;
  virtual void fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const bool doingRefinement,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual std::vector<stk::mesh::Entity> get_parent_elements_that_will_be_modified_by_unrefinement(const NodeRefiner & nodeRefiner) const override;
  virtual bool locally_have_elements_to_unrefine() const override;

  bool is_transition(const stk::mesh::Entity elem) const;

private:
  std::vector<stk::mesh::Entity> find_sorted_edge_nodes_that_will_be_removed_by_unrefinement() const;
  std::vector<stk::mesh::Entity> get_edge_nodes_of_transition_elements(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & transitionElements) const;
  std::vector<Edge> get_parent_edges_for_given_refined_edge_nodes(const std::vector<stk::mesh::Entity> & refinedEdgeNodes) const;
  std::vector<std::array<stk::mesh::EntityId,2>> get_parent_edge_node_ids_for_given_refined_edge_nodes(const std::vector<stk::mesh::Entity> & refinedEdgeNodes) const;
  void fill_are_edge_nodes_being_unrefined(const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement,
    const std::vector<stk::mesh::Entity> & refinedEdgeNodes,
    std::vector<bool> & areEdgeNodesBeingUnrefined,
    bool & areAnyEdgeNodesBeingUnrefined) const;
  void fill_post_unrefinement_edge_nodes_and_caseIds(const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    const std::vector<stk::mesh::Entity> & refinedEdgeNodes,
    const std::vector<bool> & areEdgeNodesBeingUnrefined,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    ElementEdgeCaseIds & elementEdgeCaseIds) const;
  std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childElems) const;
  std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::Entity parentElem) const;
  void fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(
    const stk::mesh::Entity parentElem,
    const stk::topology & parentElemTopology,
    const unsigned parentElemNumEdges,
    const std::vector<stk::mesh::Entity> & childTransitionElements,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void fill_refined_edge_nodes_for_marked_edges(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    bool & areAnyEdgesMarked) const;
  void fill_in_existing_refined_edge_nodes_and_caseIds(const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes,
    ElementEdgeCaseIds & elementEdgeCaseIds) const;
  void fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  bool is_element_a_candidate_for_unrefinement(const stk::mesh::Entity elem) const;
  bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const;
  void fill_existing_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, const unsigned elemNumEdges, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void mark_unrefined_edges_of_partially_refined_parent_element(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childTransitionElements, std::vector<Edge> & elemEdgesWorkspace, NodeRefiner & nodeRefiner, bool & wasEdgeMarked) const;
  void locally_mark_edges_of_partially_refined_parent_elements_with_marked_children(NodeRefiner & nodeRefiner) const;
  void locally_mark_edges_of_partially_refined_parent_elements_to_satisfy_template (NodeRefiner & nodeRefiner, bool & wasEdgeMarked) const;
  void locally_mark_edges_of_marked_non_transition_elements(NodeRefiner & nodeRefiner) const;
  bool are_all_children_leaves_and_marked_for_unrefinement(const std::vector<stk::mesh::Entity> & childElements) const;
  bool child_elements_are_all_leaves_and_are_transition_elements_or_marked_for_unrefinement(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childElements) const;
  bool can_edge_node_be_unrefined_based_on_locally_owned_elements(const stk::mesh::Entity refinedNode,
    std::vector<stk::mesh::Entity> & workParentEdgeElements,
    std::vector<stk::mesh::Entity> & workChildElements) const;

  bool is_parent_element_modified_by_unrefinement(const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElements,
    const std::vector<stk::mesh::Entity> & sortedOwnedOrSharedNodesToBeRemovedByUnrefinement) const;

  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_ */
