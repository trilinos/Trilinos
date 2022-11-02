/*
 * Akri_TransitionElementEdgeMarker.hpp
 *
 *  Created on: Oct 28, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "Akri_FieldRef.hpp"

namespace krino {

class Refinement;
class Edge;

class EdgeMarkerInterface
{
public:
  virtual ~EdgeMarkerInterface() {}
  virtual void mark_local_edges_to_be_refined() const = 0;
  virtual bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const = 0;
  virtual void fill_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const = 0;
};

class ElementBasedEdgeMarker : public EdgeMarkerInterface
{
public:
  ElementBasedEdgeMarker(const stk::mesh::BulkData & mesh, Refinement & refinement, const std::string & elementMarkerFieldName);
  virtual ~ElementBasedEdgeMarker() {}
  virtual void mark_local_edges_to_be_refined() const override;
  virtual void fill_edges_to_refine_for_marked_element(const stk::mesh::Entity elem, std::vector<Edge> & elementEdgesToRefine) const = 0;

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

  virtual bool is_element_a_candidate_for_refinement(const stk::mesh::Entity elem) const override;
  virtual void fill_edges_to_refine_for_marked_element(const stk::mesh::Entity elem, std::vector<Edge> & elementEdgesToRefine) const override;
  virtual void fill_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const override;

  // public for unit testing
  bool is_transition(const stk::mesh::Entity elem) const;

private:
  std::vector<stk::mesh::Entity> get_edge_nodes_of_transition_elements(const std::vector<stk::mesh::Entity> & transitionElements) const;
  std::vector<Edge> get_parent_edges_for_given_refined_edge_nodes(const std::vector<stk::mesh::Entity> & refinedEdgeNodes) const;
  void fill_existing_element_refined_edge_nodes_for_partially_refined_parent_element(
      const stk::mesh::Entity parentElem,
      const stk::topology & parentElemTopology,
      const unsigned parentElemNumEdges,
      const std::vector<stk::mesh::Entity> & childTransitionElements,
      std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;
  void fill_existing_element_refined_edge_nodes(const stk::mesh::Entity elem, const stk::topology & elemTopology, const unsigned elemNumEdges, std::vector<stk::mesh::Entity> & elemEdgeChildNodes) const;

  const stk::mesh::BulkData & myMesh;
  Refinement & myRefinement;
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TRANSITIONELEMENTEDGEMARKER_HPP_ */
