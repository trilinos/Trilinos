// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_SubElement_h
#define Akri_SubElement_h

#include <vector>
#include <set>
#include <map>
#include <memory>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Akri_CDMesh.hpp>
#include <Akri_Element.hpp>
#include <Akri_Facet.hpp>
#include <Akri_InterfaceID.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

class ProlongationNodeData;
class ProlongationPointData;
class SubElementNodeAncestry;
class Surface_Identifier;

bool node_on_negative_side_of_interface(const SubElementNode * node, const InterfaceID key);

class SubElement : public ElementObj {
public:
  SubElement( const stk::topology topo,
               const NodeVec & nodes,
               const std::vector<int> & side_ids,
               const Mesh_Element * owner);
  virtual ~SubElement() {}

  static double coordinate_relative_tol() { return 1.e-6; }
  static double parametric_distance( const SubElementNode * node0, const SubElementNode * node1 );

  void set_permutation();
  double relative_volume() const;
  double maximum_relative_angle() const;

  void decompose(CDMesh & mesh, const InterfaceID interface_key);
  virtual void decompose_edges(CDMesh & mesh, const InterfaceID interface_key);
  virtual void determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key);
  virtual void determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key);

  void handle_hanging_children(CDMesh & mesh, const InterfaceID & interface);
  virtual void fix_hanging_children(CDMesh & /*mesh*/, const InterfaceID & /*interface*/, const std::vector<int> & /*edges_with_children*/) {}

  virtual void build_quadratic_subelements(CDMesh & mesh);
  virtual void cut_face_interior_intersection_points(CDMesh & mesh, int level = 0);

  static void get_edge_position(const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, double & position );

  bool check_entity_nodes(const stk::mesh::BulkData & stkMesh) const;
  void get_owner_coord_transform(double * dOwnerdSub) const;
  const Mesh_Element & get_owner() const { return *my_owner; }

  int parent_side_id(const int iside) const;
  virtual void determine_decomposed_elem_phase(const std::vector<Surface_Identifier> & surfaceIDs) override final;
  std::vector<int> subelement_interface_signs(const InterfaceID interface, const int sign) const;
  const std::vector<int> & get_interface_signs() const { return myInterfaceSigns; }
  void initialize_interface_signs();
  void set_interface_signs(const std::vector<int> & interfaceSigns);
  void update_interface_signs(const InterfaceID interface_key, const int sign);
  unsigned num_sides() const { return my_parent_side_ids.size(); }
  void debug_subelements(const NodeVec & lnodes, const InterfaceID & interface, const int case_id) const;
  void debug() const;
  bool have_interface(const InterfaceID& interface) const;
  void determine_if_cut_by_interface(const InterfaceID & interface, bool & haveAnyCrossing, bool & haveRealCrossing) const;

protected:
  std::vector<int> get_edges_with_children(const InterfaceID & interface) const;
  bool determine_node_signs_on_edge( const CDMesh & mesh, const InterfaceID interface_key, const int i0, const int i1 );
  void set_node_signs_on_uncrossed_subelement( const InterfaceID interface );
  void process_edge( CDMesh & mesh, const InterfaceID interface_key, const int i0, const int i1 );
  void find_refined_edges(std::vector<unsigned> & refined_edges) const;
  int find_longest_bad_edge(std::vector<unsigned> & bad_edges) const;

protected:
  std::vector<int> my_parent_side_ids;
  const Mesh_Element * my_owner;
  std::vector<int> myInterfaceSigns;

private:
  //: Default constructor not allowed
  SubElement();
};

class SubElement_Tri_3 : public SubElement {
public:
  SubElement_Tri_3( const NodeVec & nodes,
                     const std::vector<int> & parent_side_ids,
                     const Mesh_Element * owner);
  virtual ~SubElement_Tri_3() {}

  virtual void decompose_edges(CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void fix_hanging_children(CDMesh & mesh, const InterfaceID & interface, const std::vector<int> & edges_with_children) override;
  virtual void build_quadratic_subelements(CDMesh & mesh) override;

  virtual void determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains) override;

protected:
  void perform_decomposition(CDMesh & mesh, const InterfaceID interface_key, const std::array<int,3> & node_signs);
  bool determine_diagonal_for_cut_triangle(const Simplex_Generation_Method & simplexMethod, const NodeVec & lnodes, const int i0, const int i1, const int i2, const int i3, const int i5);
  bool is_degenerate( NodeVec & lnodes,
                      const int i0, const int i1, const int i2 );
  void handle_tri( NodeVec & lnodes,
      const std::vector<int> & subInterfaceSigns,
      const int i0, const int i1, const int i2,
      const int s0, const int s1, const int s2,
      const bool is_interface0, const bool is_interface1, const bool is_interface2);
  void handle_quad( CDMesh & mesh,
      NodeVec & lnodes,
      const std::vector<int> & subInterfaceSigns,
      const int i0, const int i1, const int i2, const int i3,
      const int s0, const int s1, const int s2, const int s3,
      const bool is_interface0, const bool is_interface1, const bool is_interface2, const bool is_interface3,
      const bool face );
private:
  //: Default constructor not allowed
  SubElement_Tri_3();

};

class SubElement_Tri_6 : public SubElement {
public:
  SubElement_Tri_6( const NodeVec & nodes,
                     const std::vector<int> & side_ids,
                     const Mesh_Element * owner);
  virtual ~SubElement_Tri_6() {}

private:
  //: Default constructor not allowed
  SubElement_Tri_6();

};

class SubElement_Tet_4 : public SubElement {
public:
  SubElement_Tet_4( const NodeVec & nodes,
                     const std::vector<int> & side_ids,
                     const Mesh_Element * owner);
  virtual ~SubElement_Tet_4() {}

  virtual void decompose_edges(CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void fix_hanging_children(CDMesh & mesh, const InterfaceID & interface, const std::vector<int> & edges_with_children) override;
  virtual void build_quadratic_subelements(CDMesh & mesh) override;
  virtual void cut_face_interior_intersection_points(CDMesh & mesh, int level = 0) override;

  virtual void determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key) override;
  virtual void cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains) override;

protected:
  void perform_decomposition(CDMesh & mesh, const InterfaceID interface_key, const std::array<int,4> & node_signs);
  void determine_node_scores_on_face( const CDMesh & mesh, const InterfaceID interface, const int i0, const int i1, const int i2 );
  static double tet_volume(const std::array<stk::math::Vector3d,4> & nodes);

  bool is_degenerate( NodeVec & lnodes,
                      const int i0, const int i1, const int i2, const int i3 );
  void handle_tet( NodeVec & lnodes,
                   const std::vector<int> & subInterfaceSigns,
                   const int i0, const int i1, const int i2, const int i3,
                   const int s0, const int s1, const int s2, const int s3);
  void handle_pyramid( NodeVec & lnodes,
                       const std::vector<int> & subInterfaceSigns,
                       const int i0, const int i1, const int i2, const int i3, const int i4,
                       const int s0, const int s1, const int s2, const int s3, const int s4,
                       const bool face4 );
  void handle_wedge( CDMesh & mesh,
      NodeVec & lnodes,
      const std::vector<int> & subInterfaceSigns,
      const int i0, const int i1, const int i2, const int i3, const int i4, const int i5,
      const int s0, const int s1, const int s2, const int s3, const int s4,
      const bool face0, const bool face1, const bool face2 );
private:
  //: Default constructor not allowed
  SubElement_Tet_4();

  bool determine_diagonal_for_cut_triangular_face(const Simplex_Generation_Method & simplexMethod, const bool globalIDsAreParallelConsistent, const NodeVec & lnodes, const int i0, const int i1, const int i2, const int i3, const int i5);
  void cut_face_intersection_point_with_permutation(CDMesh & mesh, const std::array<int,4> & permuteNodes, const std::array<int,4> & permuteSides, const std::vector<double> & faceNodeWeights, const std::vector<int> & sortedDomains);
};

class SubElement_Tet_10 : public SubElement {
public:
  SubElement_Tet_10( const NodeVec & nodes,
                      const std::vector<int> & side_ids,
                      const Mesh_Element * owner);
  virtual ~SubElement_Tet_10() {}

private:
  //: Default constructor not allowed
  SubElement_Tet_10();

};

class SubElementNode {
public:
  SubElementNode(const Mesh_Element * in_owner)
  : my_entity(),
    my_entityId(0),
    my_is_prolonged_flag(false),
    my_cached_owner(in_owner)
  {}

  virtual ~SubElementNode() {}
  virtual bool is_mesh_node() const { return false; }
  virtual NodeVec get_parents() const { throw std::runtime_error("Incorrect usage of SubElementNode.  The type of node does not have parents."); }
  virtual size_t get_num_parents() const { throw std::runtime_error("Incorrect usage of SubElementNode.  The type of node does not have parents."); }
  virtual std::vector<double> get_parent_weights() const { throw std::runtime_error("Incorrect usage of SubElementNode.  The type of node does not have parents."); }
  virtual bool needs_to_be_ale_prolonged(const CDMesh & mesh) const = 0;
  virtual void prolongate_fields(const CDMesh & mesh) const = 0;

  void fill_parent_entity_pointers(std::vector<stk::mesh::Entity*> & parentEntities) const;

  void set_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity in_node_obj) const { my_entity = in_node_obj; my_entityId = mesh.identifier(my_entity); }
  void set_entity(stk::mesh::Entity nodeEntity, stk::mesh::EntityId nodeEntityId) const { my_entity = nodeEntity; my_entityId = nodeEntityId; }
  const Mesh_Element * get_active_element() const { return my_cached_owner; }
  bool is_prolonged() const { return my_is_prolonged_flag; }
  void set_prolonged_flag(const bool val) const { my_is_prolonged_flag = val; } // this is dangerous, know what you are doing
  bool on_common_edge(const SubElementNode * other) const;

  bool entity_is_valid(const stk::mesh::BulkData & mesh) const { return mesh.is_valid(my_entity); }
  bool check_entity(const stk::mesh::BulkData & mesh) const { return (my_entityId == 0) ? (my_entity == stk::mesh::Entity()) : (my_entityId == mesh.identifier(my_entity)); }
  stk::mesh::Entity & entity() const { return my_entity; }
  stk::mesh::EntityId entityId() const { return my_entityId; }
  void set_entityId_from_entity(const stk::mesh::BulkData & mesh) const { STK_ThrowAssert(mesh.is_valid(my_entity)); my_entityId = mesh.identifier(my_entity); }

  const stk::math::Vector3d & owner_coords( const Mesh_Element * in_owner ) const {
    if ( in_owner != my_cached_owner ) {
      // call derived types function for recomputing owner_coords
      my_cached_owner = in_owner;
      my_cached_owner_coords = compute_owner_coords( in_owner );
    }
    return my_cached_owner_coords;
  }
  bool get_node_on_interface() const { return my_node_sign == 0; }

  void set_node_domains(const std::vector<int> & nodeDomains) const { my_sorted_node_domains = nodeDomains; }
  const std::vector<int> & get_sorted_node_domains() const { return my_sorted_node_domains; }
  void insert_node_domains(const std::vector<int> & domainsToAdd) const;
  bool captures_intersection_point_domains(const std::vector<int> & intersectionPointDomains) const;
  bool captures_interface(const InterfaceID & interface) const;
  bool captures_intersection_point_domains(const InterfaceID & interface) const;
  static bool higher_priority_by_score_then_ancestry(const SubElementNode & a, const SubElementNode & b, const bool globalIDsAreParallelConsistent);
  static bool less_by_entity_id(const SubElementNode & a, const SubElementNode & b);
  static bool less_by_coordinates_then_by_entity_id(const SubElementNode & a, const SubElementNode & b);
  bool node_sign_is_set() const { return my_node_sign != -2; }
  void set_node_sign(const int sign) const { my_node_sign = (my_node_sign == 0 || my_node_sign == -sign) ? 0 : sign; }
  int get_node_sign() const { STK_ThrowRequire(node_sign_is_set()); return my_node_sign; }
  void clear_node_sign() const { my_node_sign = -2; }
  void clear_node_score() const { my_node_score = 1.; }
  bool node_score_is_set() const { return my_node_score != 1.; }
  void set_node_score(const double score) const { my_node_score = std::min(my_node_score, score); }
  double get_node_score() const { STK_ThrowRequire(node_sign_is_set()); return my_node_score; }

  // pure virtual function for recomputing local coordinates
  virtual stk::math::Vector3d compute_owner_coords( const Mesh_Element * in_owner ) const = 0;

  const stk::math::Vector3d & coordinates() const { return my_global_coords; }

  void add_child(const SubElementNode* child) const { my_children.push_back(child); }
  static const SubElementNode * common_child( const NodeVec & parents );
  bool have_child() const;
  bool have_child(const SubElementNode* child) const;

  // This is called from Mesh::find_prolongation_node to determine the set of candidate prolongation nodes.
  std::vector<unsigned> prolongation_node_fields(const CDMesh & mesh) const;

  void get_ancestors(NodeSet & ancestors) const;
  SubElementNodeAncestry get_ancestry() const;
  void build_stencil(std::map<const SubElementNode *, double> & stencil, const double self_weight = 1.0) const;
  void build_constraint_stencil(const FieldRef field, std::vector<stk::mesh::Entity> & entities, std::vector<double> & weights) const;

protected:
  void prolong_cdfem_displacements(const CDMesh & mesh,
      const ProlongationPointData * prolong_data,
      const bool zero_if_no_prolong_data = true) const;
  void prolong_zeroed_fields(const CDMesh & mesh, const ProlongationNodeData * nodeToExamineForExistingField) const;

  mutable stk::mesh::Entity my_entity;
  mutable stk::mesh::EntityId my_entityId;
  mutable bool my_is_prolonged_flag;
  stk::math::Vector3d my_global_coords;
  mutable signed char my_node_sign;
  mutable double my_node_score;

  mutable const Mesh_Element * my_cached_owner;
  mutable stk::math::Vector3d my_cached_owner_coords;

  mutable std::vector<const SubElementNode*> my_children;
  mutable std::vector<int> my_sorted_node_domains;

private:
  //: copy constructor not allowed
  SubElementNode(SubElementNode const & copy);
};

class SubElementChildNode : public SubElementNode {
public:
  SubElementChildNode( const Mesh_Element * owner,
      const NodeVec & parents,
      const std::vector<double> & weights );

  virtual ~SubElementChildNode() {}
  virtual bool needs_to_be_ale_prolonged(const CDMesh & mesh) const override;
  virtual void prolongate_fields(const CDMesh & mesh) const override;
  virtual stk::math::Vector3d compute_owner_coords( const Mesh_Element * in_owner ) const override;

  virtual NodeVec get_parents() const override { return my_parents; }
  virtual size_t get_num_parents() const override { return my_parents.size(); }
  virtual std::vector<double> get_parent_weights() const override { return my_weights; }

private:
  void prolong_ale_fields(const CDMesh & mesh, const ProlongationPointData * prolong_data) const;
  void prolong_interpolation_fields(const CDMesh & mesh) const;
  void prolong_edge_interpolation_fields(const FieldSet & edgeInterpFields) const;

  NodeVec my_parents;
  std::vector<double> my_weights;
};

class SubElementSteinerNode : public SubElementChildNode {
public:

  SubElementSteinerNode( const Mesh_Element * in_owner,
      const NodeVec & parents,
      const std::vector<double> & weights )
  : SubElementChildNode(in_owner, parents, weights) {}

  virtual ~SubElementSteinerNode() {}
  virtual bool needs_to_be_ale_prolonged(const CDMesh & /*mesh*/) const override { return false; }
  virtual void prolongate_fields(const CDMesh & mesh) const override;
  virtual stk::math::Vector3d compute_owner_coords( const Mesh_Element * /*owner*/ ) const override {
    throw std::runtime_error("Incorrect usage of SubElementSteinerNode.  The type of node only has one owner.");
  }
};

class SubElementEdgeNode : public SubElementChildNode {
public:
  SubElementEdgeNode( const Mesh_Element * owner,
      const double & position,
      const SubElementNode *parent1,
      const SubElementNode *parent2)
  : SubElementChildNode( owner, {parent1, parent2}, {1.-position, position}) {}

  virtual ~SubElementEdgeNode() {}

  double get_position() const { STK_ThrowAssert(get_num_parents() == 2); return get_parent_weights()[1]; }
  double get_position(const SubElementNode *parent1, [[maybe_unused]] const SubElementNode *parent2) const { STK_ThrowAssert(check_parents(parent1,parent2)); return ((parent1==get_parents()[0]) ? get_parent_weights()[1] : get_parent_weights()[0]); }

private:
  bool check_parents(const SubElementNode *parent1, const SubElementNode *parent2) const {
    return (get_num_parents() == 2 && ((parent1==get_parents()[0] && parent2==get_parents()[1]) || (parent2==get_parents()[0] && parent1==get_parents()[1])));
  }
};

  class SubElementMidSideNode : public SubElementNode {
public:
  SubElementMidSideNode() = delete;
  SubElementMidSideNode( const Mesh_Element * owner,
      const SubElementNode *parent1,
      const SubElementNode *parent2);
  SubElementMidSideNode( const Mesh_Element * owner,
      const SubElementNode *parent1,
      const SubElementNode *parent2,
      stk::mesh::Entity meshNode,
      stk::mesh::EntityId meshNodeId);

  virtual ~SubElementMidSideNode() {}
  virtual bool needs_to_be_ale_prolonged(const CDMesh & /*mesh*/) const override { return false; }
  virtual void prolongate_fields(const CDMesh & mesh) const override;
  virtual stk::math::Vector3d compute_owner_coords( const Mesh_Element * in_owner ) const override {
    return 0.5 * my_parent1->owner_coords(in_owner) + 0.5 * my_parent2->owner_coords(in_owner);
  }
  virtual bool is_mesh_node() const override { return my_is_mesh_node; }

  virtual NodeVec get_parents() const override { return NodeVec{my_parent1,my_parent2}; }
  virtual size_t get_num_parents() const override { return 2; }
  virtual std::vector<double> get_parent_weights() const override { return std::vector<double>{0.5,0.5}; }

protected:
  bool my_is_mesh_node;
  const SubElementNode *my_parent1;
  const SubElementNode *my_parent2;

private:
  void prolong_ale_fields(const CDMesh & mesh) const;
  void prolong_interpolation_fields(const CDMesh & mesh) const;
  void prolong_edge_interpolation_fields(const FieldSet & edgeInterpFields) const;
  bool is_mesh_node_that_needs_to_be_prolonged(const CDMesh & mesh) const;
};


class SubElementMeshNode : public SubElementNode {
public:

  SubElementMeshNode( const Mesh_Element * owner,
      stk::mesh::Entity nodeEntity,
      stk::mesh::EntityId nodeEntityId,
      const stk::math::Vector3d & owner_coords,
      const stk::math::Vector3d & global_coords );

  virtual ~SubElementMeshNode() {}

  virtual bool needs_to_be_ale_prolonged(const CDMesh & mesh) const override;
  virtual void prolongate_fields(const CDMesh & mesh) const override;
  virtual stk::math::Vector3d compute_owner_coords( const Mesh_Element * owner ) const override { return owner->get_node_parametric_coords(this); }
  virtual bool is_mesh_node() const override { return true; }
protected:

private:
  void prolong_ale_fields(const CDMesh & mesh,
      const ProlongationPointData * prolong_data,
      const ProlongationNodeData * old_node) const;
  void prolong_interpolation_fields(const CDMesh & mesh, const ProlongationNodeData * old_node) const;
  //: Default constructor not allowed
  SubElementMeshNode();
};

} // namespace krino

#endif // Akri_SubElement_h
