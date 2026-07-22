// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Element_h
#define Akri_Element_h

#include <Akri_TypeDefs.hpp>
#include <vector>
#include <set>
#include <map>

#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_MasterElement.hpp>
#include <stk_math/StkVector.hpp>


namespace krino {

class SubElement;
class SubElementNode;
class ProlongationElementData;
class InterfaceGeometry;
class IntersectionPoint;
struct ElementIntersection;
class Surface_Identifier;

class ElementObj {
public:

  ElementObj(const stk::mesh::BulkData & stkMesh, stk::mesh::Entity elemEntity);
  ElementObj(const stk::topology elem_topology, const NodeVec & nodes);

  virtual ~ElementObj(); // Definition must in implementation file because SubElement is incomplete

  static bool is_less(const ElementObj * elem0,const ElementObj * elem1) { return (elem0->entityId() < elem1->entityId()); }

  static bool will_cutting_quad_from_0to2_cut_largest_angle(const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3);
  static bool determine_diagonal_for_internal_quad_of_cut_tet_from_owner_nodes(const SubElementNode * pn0, const SubElementNode * pn1, const SubElementNode * pn2, const SubElementNode * pn3);
  static bool determine_diagonal_for_internal_quad_of_cut_tet_from_edge_nodes(const Simplex_Generation_Method simplexMethod, const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3,
      const bool face0, const bool face1, const bool face2, const bool face3);

  int spatial_dim() const { return my_master_elem.topology_dimension(); }
  const NodeVec & get_nodes() const { return my_nodes; }
  unsigned num_nodes() const { return my_nodes.size(); }
  const PhaseTag & get_phase() const { return my_phase; }
  void set_phase(const PhaseTag & phase) { my_phase = phase; }
  stk::mesh::EntityId entityId() const { return my_entityId; }
  bool check_entity(const stk::mesh::BulkData & mesh) const { return (my_entityId == 0) ? (my_entity == stk::mesh::Entity()) : (mesh.is_valid(my_entity) && my_entityId == mesh.identifier(my_entity)); }
  stk::mesh::Entity entity() const { return my_entity; }
  void set_entity( const stk::mesh::BulkData & mesh, stk::mesh::Entity meshEntity ) const { my_entity = meshEntity; my_entityId = mesh.identifier(my_entity); }

  const MasterElement & master_elem() const { return my_master_elem; }
  stk::topology topology() const { return my_master_elem.get_topology(); }

  void fill_node_owner_coords(const Mesh_Element * owner, std::vector<stk::math::Vector3d> & coords) const;
  stk::math::Vector3d compute_local_coords_from_owner_coordinates(const Mesh_Element * owner, const stk::math::Vector3d & coords) const;

  static void integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const MasterElement & me,
    const MasterElement & mesh_me);
  static void integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const MasterElement & me) { integration_weights(intg_weights, numCoordDims, mesh_coords, me, me); }
  static void integration_locations(
    std::vector<stk::math::Vector3d> & intg_pt_locations,
    const MasterElement & me);
  static void gather_nodal_field(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef field, std::vector<double> & result, unsigned dim = 1, unsigned npe = 0);
  static double volume(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coords_field);
  static PhaseTag update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const InterfaceID interface_key, const int sign);
  static PhaseTag update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const std::vector<InterfaceID> & interfaces, const std::vector<int> & interfaceSigns);

  const MasterElement* get_evaluation_master_element(const FieldRef field) const;
  void evaluate_prolongation_field(const CDMesh & mesh, const FieldRef field, const unsigned field_length, const stk::math::Vector3d & p_coords, double * result) const;

  void add_subelement(std::unique_ptr<SubElement> subelem);
  void get_subelements( std::vector<const SubElement *> & subelems ) const;
  void get_subelements( std::vector<SubElement *> & subelems );
  bool have_subelements() const { return !my_subelements.empty(); }
  unsigned num_subelements() const { return my_subelements.size(); }
  bool have_refined_edges() const;

  virtual void determine_decomposed_elem_phase(const std::vector<Surface_Identifier> & /*surfaceIDs*/) { STK_ThrowRequire(false); }
  virtual void cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains);
  virtual void cut_face_intersection_point(const int iFace, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains);
  bool have_edges_with_children() const;
  bool captures_intersection_point_domains(const std::vector<int> & intersectionPointDomains) const;
  bool face_captures_intersection_point_domains(const int iFace, const std::vector<int> & intersectionPointDomains) const;

  ProlongationElementData * get_prolongation_data() const { return my_prolongation_data; }
  void set_prolongation_data( ProlongationElementData * data ) const { my_prolongation_data = data; }

  void prolongate_fields(const CDMesh & mesh) const;

protected:
  const MasterElement& my_master_elem;
  std::vector<const SubElementNode *> my_nodes;
  std::vector<std::unique_ptr<SubElement>> my_subelements;
  PhaseTag my_phase;
  mutable stk::mesh::Entity my_entity;
  mutable stk::mesh::EntityId my_entityId;
  mutable ProlongationElementData * my_prolongation_data;

private:
  //: Default constructor not allowed
  ElementObj();
};

class Mesh_Element : public ElementObj {
public:

  Mesh_Element(CDMesh & mesh, stk::mesh::Entity elemEntity);

  virtual ~Mesh_Element();

  static bool is_supported_topology(stk::topology elem_topology) { return determine_subelement_topology(elem_topology).first != stk::topology::INVALID_TOPOLOGY; }
  static std::pair<stk::topology, unsigned> determine_subelement_topology(stk::topology elem_topology);

  stk::math::Vector3d get_node_parametric_coords( const int lnn ) const;
  stk::math::Vector3d get_node_parametric_coords( const SubElementNode * node ) const;

  const MasterElement & coord_master_elem() const { return my_master_elem; }

  stk::topology coord_topology() const { return my_master_elem.get_topology(); }

  const stk::topology & subelement_topology() const { return my_subelement_topology; }
  int subelement_order() const { return my_subelement_order; }

  std::string visualize(const CDMesh & mesh) const;
  double interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const;
  std::pair<int, double> interface_edge_crossing_sign_and_position(const InterfaceID interface, const SubElementNode * node1, const SubElementNode * node2) const;
  void fill_face_interior_intersections(const NodeVec & faceNodes, std::vector<ElementIntersection> & faceIntersectionPoints) const;
  double interface_crossing_position(const InterfaceID interface, const SubElementNode * node1, const SubElementNode * node2) const;
  std::function<bool(const std::array<unsigned,4> &)> get_diagonal_picker() const;

  const ElementCutter * get_cutter() const { return myCutter.get(); }
  ElementCutter * get_cutter() { return myCutter.get(); }

  bool is_prolonged() const;
  stk::math::Vector3d coordinates( const stk::math::Vector3d & p_coords ) const;

  PhaseTag determine_phase_at_location(const stk::math::Vector3d & location) const;
  bool have_interface() const { return my_have_interface; }
  bool have_interface(const InterfaceID interface) const { return std::binary_search(myCuttingInterfaces.begin(), myCuttingInterfaces.end(), interface); }
  int get_num_interfaces() const { return myCuttingInterfaces.size(); }
  int get_interface_index(const InterfaceID interface) const;
  int get_interface_sign_for_uncrossed_subelement(const InterfaceID interface, const std::vector<stk::math::Vector3d> & elemNodeCoords) const;
  const std::vector<InterfaceID> & get_sorted_cutting_interfaces() const { return myCuttingInterfaces; }
  virtual void determine_decomposed_elem_phase(const std::vector<Surface_Identifier> & surfaceIDs) override;
  void set_have_interface() { my_have_interface = true; }

  bool triangulate(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry); //return value indicates if any changes were made
  void create_cutter(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry);
  void create_base_subelement();
  void determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key);
  void determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key);
  void decompose(CDMesh & mesh, const InterfaceID interface_key);
  void handle_hanging_children(CDMesh & mesh, const InterfaceID & interface);
  void build_quadratic_subelements(CDMesh & mesh);
  std::vector<int> get_interface_signs_based_on_crossings(const NodeVec & nodes) const;
  void cut_interior_intersection_points(CDMesh & mesh);

  void find_child_coordinates_at_owner_coordinates(const stk::math::Vector3d & ownerCoordinates, const ElementObj *& child, stk::math::Vector3d & child_p_coords) const;

  bool is_single_coincident() const;

private:
  Mesh_Element() = delete;
  stk::math::Vector3d get_intersection_point_parametric_coordinates(const IntersectionPoint & intersectionPoint) const;
  std::vector<double> gather_nodal_coordinates() const;

  int my_subelement_order;
  stk::topology my_subelement_topology;
  bool my_have_interface;

  std::unique_ptr<ElementCutter> myCutter;
  std::vector<InterfaceID> myCuttingInterfaces;
};

} // namespace krino

#endif // Akri_Element_h
