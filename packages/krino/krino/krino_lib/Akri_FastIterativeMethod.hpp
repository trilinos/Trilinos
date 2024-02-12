// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_FASTITERATIVEMETHOD_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_FASTITERATIVEMETHOD_HPP_

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/diag/Timer.hpp>

namespace krino {

class ParallelErrorMessage;

class FastIterativeMethod {
public:
  FastIterativeMethod(const stk::mesh::BulkData & mesh,
      const stk::mesh::Selector & selector,
      const FieldRef& coordinates,
      const FieldRef& distance,
      const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> & get_interface_speed,
      stk::diag::Timer & parentTimer);

  void redistance();
  bool check_converged_solution() const;

private:
  std::set<stk::mesh::Entity> initialize(ParallelErrorMessage& err);
  std::set<stk::mesh::Entity> build_initial_working_set(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & initialNodes);
  void advance(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & initialNodes, std::set<stk::mesh::Entity> & workingSet);

  const stk::mesh::BulkData& mesh() const { return myMesh; }
  stk::mesh::Selector field_not_ghost_selector() const;
  stk::mesh::Selector field_selector() const;

  bool have_crossing(const stk::mesh::Entity & elem) const;
  void initialize_element(const stk::mesh::Entity & elem, const double speed);
  void fill_noninitial_node_neighbors(const std::set<stk::mesh::Entity> & initialNodes, stk::mesh::Entity node, std::vector<stk::mesh::Entity> & nodeNbrs) const;

  double updated_node_distance(ParallelErrorMessage& err, const stk::mesh::Entity & node) const;
  double element_signed_distance_for_node(ParallelErrorMessage& err, const stk::mesh::Entity & elem, const stk::mesh::Entity & node) const;

  void parallel_communicate_nodes(ParallelErrorMessage& err, std::set<stk::mesh::Entity> & nodes) const;
  void parallel_communicate_node_distances(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & nodes, std::vector<double> & distances) const;
  void parallel_communicate_nodes_and_nodal_distance(ParallelErrorMessage& err, std::set<stk::mesh::Entity> & nodes);

  std::vector<double> compute_updates_and_communicate(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & nodes) const;

  double update_triangle(const stk::mesh::Entity * elemNodes, int node_to_update, const double speed) const;
  double update_tetrahedron(const stk::mesh::Entity * elemNodes, int node_to_update, const double speed) const;

  const stk::mesh::BulkData & myMesh;
  stk::mesh::Selector mySelector;
  const FieldRef myCoordinates;
  const FieldRef myDistance;
  const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> my_get_interface_speed;
  mutable stk::diag::Timer myTimer;

  void check_error(const ParallelErrorMessage& err, const std::string & context) const;
};

} // krino


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_FASTITERATIVEMETHOD_HPP_ */
