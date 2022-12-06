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
#include "Akri_FieldRef.hpp"

namespace krino {

class Refinement;
class Edge;
class NodeRefiner;

class EdgeMarkerInterface
{
public:
  virtual ~EdgeMarkerInterface() {}
  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const = 0;
  virtual bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const = 0;
  virtual void fill_element_refined_edge_nodes(const NodeRefiner & nodeRefiner, const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const = 0;
  virtual void fill_elements_modified_by_unrefinement(std::vector<stk::mesh::Entity> & elementsWithoutChildrenAfterUnrefinement, std::vector<stk::mesh::Entity> & childElementsToDeleteForUnrefinement) const = 0;
  virtual bool locally_have_elements_to_unrefine() const = 0;
};

class UniformEdgeMarker : public EdgeMarkerInterface
{
public:
  UniformEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement);
  virtual ~UniformEdgeMarker() {}
protected:
  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const override;
  virtual void fill_element_refined_edge_nodes(const NodeRefiner & nodeRefiner, const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual void fill_elements_modified_by_unrefinement(std::vector<stk::mesh::Entity> & elementsWithoutChildrenAfterUnrefinement, std::vector<stk::mesh::Entity> & childElementsToDeleteForUnrefinement) const override;
  virtual bool locally_have_elements_to_unrefine() const { return false; }
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
  FieldRef get_marker_field() const { return myElementMarkerField; }
private:
  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
  FieldRef myElementMarkerField;
};

class TransitionElementEdgeMarker : public ElementBasedEdgeMarker
{
public:
  TransitionElementEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement, const std::string & elementMarkerFieldName);
  TransitionElementEdgeMarker ( const TransitionElementEdgeMarker & ) = delete;
  TransitionElementEdgeMarker & operator= ( const TransitionElementEdgeMarker & ) = delete;
  virtual ~TransitionElementEdgeMarker() {}

  virtual void mark_edges_to_be_refined(NodeRefiner & nodeRefiner) const override;
  virtual bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const override;
  virtual void fill_element_refined_edge_nodes(const NodeRefiner & nodeRefiner, const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;
  virtual void fill_elements_modified_by_unrefinement(std::vector<stk::mesh::Entity> & elementsWithoutChildrenAfterUnrefinement, std::vector<stk::mesh::Entity> & childElementsToDeleteForUnrefinement) const override;
  virtual bool locally_have_elements_to_unrefine() const;

  bool is_transition(const stk::mesh::Entity elem) const;

private:
  std::vector<stk::mesh::Entity> get_edge_nodes_of_transition_elements(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & transitionElements) const;
  std::vector<Edge> get_parent_edges_for_given_refined_edge_nodes(const std::vector<stk::mesh::Entity> & refinedEdgeNodes) const;
  void fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(
      const stk::mesh::Entity parentElem,
      const stk::topology & parentElemTopology,
      const unsigned parentElemNumEdges,
      const std::vector<stk::mesh::Entity> & childTransitionElements,
      std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void fill_existing_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, const unsigned elemNumEdges, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void mark_unrefined_edges_of_partially_refined_parent_element(const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & childTransitionElements, std::vector<Edge> & elemEdgesWorkspace, NodeRefiner & nodeRefiner, bool & wasEdgeMarked) const;
  void locally_mark_edges_of_partially_refined_parent_elements_with_marked_children(NodeRefiner & nodeRefiner) const;
  void locally_mark_edges_of_partially_refined_parent_elements_to_satisfy_template (NodeRefiner & nodeRefiner, bool & wasEdgeMarked) const;
  void locally_mark_edges_of_marked_non_transition_elements(NodeRefiner & nodeRefiner) const;
  bool are_all_children_leaves_and_marked_for_unrefinement(const std::vector<stk::mesh::Entity> & childElements) const;

  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_ */
