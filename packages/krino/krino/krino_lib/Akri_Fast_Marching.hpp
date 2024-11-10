// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Fast_Marching_h
#define Akri_Fast_Marching_h

#include <array>
#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/diag/Timer.hpp>

namespace krino {

class SubElement;
class Mesh_Element;
class AuxMetaData;
class ParallelErrorMessage;
struct StkMeshEntities;

enum Enum_Fast_Marching_Node_Status{STATUS_UNUSED=0, STATUS_INITIAL, STATUS_ACCEPTED, STATUS_TRIAL, STATUS_FAR};

class Fast_Marching_Node
{
public:
  Fast_Marching_Node()
  : my_node(), my_status(STATUS_UNUSED), my_signed_dist(std::numeric_limits<double>::max()), my_on_neg_side(false) {}
  Fast_Marching_Node(stk::mesh::Entity in_node, Enum_Fast_Marching_Node_Status in_status, double in_dist, int in_sign, stk::math::Vector3d in_coords)
  : my_node(in_node), my_status(in_status), my_signed_dist(in_dist), my_on_neg_side(in_sign<0), my_coords(in_coords) {}

  int sign() const { return (my_on_neg_side ? (-1) : 1); }
  void set_sign(const int in_sign) { my_on_neg_side = (in_sign<0); }
  stk::mesh::Entity node() const { return my_node; }
  Enum_Fast_Marching_Node_Status status() const { return my_status; }
  const stk::math::Vector3d & coords() const { return my_coords; }
  double signed_dist() const { return my_signed_dist; }
  void set_signed_dist(double dist) { my_signed_dist = dist; }
  void set_status(Enum_Fast_Marching_Node_Status status) { my_status = status; }

private:
  stk::mesh::Entity my_node;
  Enum_Fast_Marching_Node_Status my_status;
  double my_signed_dist;
  bool my_on_neg_side;
  stk::math::Vector3d my_coords;
};

class Fast_Marching_Node_Less
{
public:
  Fast_Marching_Node_Less(const stk::mesh::BulkData & mesh) : mMesh(mesh) {}
  bool operator() (const Fast_Marching_Node *a, const Fast_Marching_Node *b) const
  {
    const double unsignedDistA = a->signed_dist()*a->sign();
    const double unsignedDistB = b->signed_dist()*b->sign();
    if (unsignedDistA < unsignedDistB) return true;
    if (unsignedDistB < unsignedDistA) return false;
    return mMesh.identifier(a->node()) < mMesh.identifier(b->node());
  }
private:
  const stk::mesh::BulkData & mMesh;
};

class Fast_Marching {
public:
  Fast_Marching(const stk::mesh::BulkData & mesh,
      const stk::mesh::Selector & activeElementSelector,
      const FieldRef& coordinates,
      const FieldRef& distance,
      const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> & get_interface_speed,
      stk::diag::Timer & parentTimer);

  void redistance();
  void redistance(const std::vector<stk::mesh::Entity> & initialNodes);
  void initialize(ParallelErrorMessage& err);
  void update_neighbors(Fast_Marching_Node & accepted_node, ParallelErrorMessage & err);
  void update_node(std::vector<Fast_Marching_Node *> & elem_nodes, int node_to_update, const double speed);

  bool have_crossing(const StkMeshEntities & elemNodes) const;
  bool have_crossing(const stk::mesh::Entity & elem) const;
  void initialize_element(const stk::mesh::Entity & elem, ParallelErrorMessage& err);
  double update_triangle(const std::array<Fast_Marching_Node *, 3> & elemNodes, const int node_to_update, const double speed);
  double update_tetrahedron(const std::array<Fast_Marching_Node *, 4> & elemNodes, const int node_to_update, const double speed);

  void add_trial_node(Fast_Marching_Node & add_trial_node);
  void update_trial_node(Fast_Marching_Node & add_trial_node, const double dist);
  Fast_Marching_Node * get_fm_node(stk::mesh::Entity node);

private:
  void initialize_algorithm();
  void initialize_nodes_of_crossed_elements();
  void initialize_given_nodes(const std::vector<stk::mesh::Entity> & initialNodes);
  void propagate_distance_to_convergence();

  void initialize_node(const stk::mesh::Entity node);
  double & get_node_distance(const stk::mesh::Entity node);
  double get_node_distance(const stk::mesh::Entity node) const;
  double get_element_interface_speed(ParallelErrorMessage& err, const stk::mesh::Entity elem) const;
  stk::mesh::Selector selected_with_field_not_ghost_selector() const;
  stk::mesh::Selector selected_with_field_selector() const;

  void update_neighbors_from_triangle_element(const Fast_Marching_Node & acceptedNode, const stk::mesh::Entity elem, const double elemSpeed);
  void update_neighbors_from_quadrilateral_element(const Fast_Marching_Node & acceptedNode, const stk::mesh::Entity elem, const double elemSpeed);
  void update_neighbors_from_tetrahedron_element(const Fast_Marching_Node & acceptedNode, const stk::mesh::Entity elem, const double elemSpeed);

  template<size_t NSIMPLICES, size_t NNODES> void update_neighbors_from_simplices(const Fast_Marching_Node & acceptedNode,
      const stk::mesh::Entity* elemNodes,
      const std::array<std::array<int, NNODES>,NSIMPLICES> & simplices,
      const double elemSpeed);
  template<size_t NNODES> void update_neighbors_from_simplex(const Fast_Marching_Node & acceptedNode, const std::array<stk::mesh::Entity, NNODES> & simplexNodes, const double elemSpeed);
  template <size_t NNODES> void update_trial_node_from_simplex(const std::array<Fast_Marching_Node *, NNODES> & simplexNodes, const int nodeToUpdate, const double speed);

  template<size_t NNODES> void initialize_from_simplex(const std::array<stk::mesh::Entity, NNODES> & simplexNodes,
      const double elemSpeed,
      const std::function<const stk::math::Vector3d &(stk::mesh::Entity)> & get_coordinates);
  template<size_t NSIMPLICES, size_t NNODES> void initialize_from_simplices(const stk::mesh::Entity* elemNodes, const std::array<std::array<int, NNODES>,NSIMPLICES> & simplices,
    const double elemSpeed,
    const std::function<const stk::math::Vector3d &(stk::mesh::Entity)> & get_coordinates);

  const stk::mesh::BulkData & myMesh;
  stk::mesh::Selector mySelector;
  const FieldRef myCoordinates;
  const FieldRef myDistance;
  const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> my_get_interface_speed;
  stk::diag::Timer myTimer;

  std::vector<Fast_Marching_Node> fm_nodes;
  Fast_Marching_Node_Less my_fm_node_less;
  std::set<Fast_Marching_Node*, Fast_Marching_Node_Less> trial_nodes; //set sorted by distance, then by id to break "ties"

  void check_error(const ParallelErrorMessage& err, const std::string & context) const;
};

} // namespace krino

#endif // Akri_Fast_Marching_h
