// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Element_h
#define Akri_Element_h

#include <array>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_ParametricElement.hpp>
#include <Akri_TypeDefs.hpp>
#include <stk_math/StkVector.hpp>

#include <functional>

namespace krino {

class CubitJournal;
class SubElement;
class SubElementNode;
class ProlongationElementData;
class InterfaceGeometry;
class ElementCutter;
class IntersectionPoint;
struct ElementIntersection;
class Surface_Identifier;
class ParametricElementLibrary;

class ElementObj {
public:

  ElementObj(ParametricElementLibrary & parametricElemLibrary, const stk::mesh::BulkData & stkMesh, stk::mesh::Entity elemEntity);
  ElementObj(ParametricElementLibrary & parametricElemLibrary, const stk::topology elem_topology, const NodeVec & nodes);

  virtual ~ElementObj(); // Definition must in implementation file because SubElement is incomplete

  static bool is_less(const ElementObj * elem0,const ElementObj * elem1) { return (elem0->entityId() < elem1->entityId()); }

  static bool will_cutting_quad_from_0to2_cut_largest_angle(const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3);
  static bool determine_diagonal_for_internal_quad_of_cut_tet_from_owner_nodes(const SubElementNode * pn0, const SubElementNode * pn1, const SubElementNode * pn2, const SubElementNode * pn3);
  static bool determine_diagonal_for_internal_quad_of_cut_tet_from_edge_nodes(const Simplex_Generation_Method simplexMethod, const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3,
      const bool face0, const bool face1, const bool face2, const bool face3);

  int spatial_dim() const { return my_parametric_elem.topology_dimension(); }
  const NodeVec & get_nodes() const { return my_nodes; }
  unsigned num_nodes() const { return my_nodes.size(); }
  const PhaseTag & get_phase() const { return *my_phase; }
  // The phase must be registered (CDMesh::get_registered_phase), so that it outlives this element.
  void set_registered_phase(const PhaseTag & registeredPhase) { my_phase = &registeredPhase; }
  stk::mesh::EntityId entityId() const { return my_entityId; }
  bool check_entity(const stk::mesh::BulkData & mesh) const { return (my_entityId == 0) ? (my_entity == stk::mesh::Entity()) : (mesh.is_valid(my_entity) && my_entityId == mesh.identifier(my_entity)); }
  stk::mesh::Entity entity() const { return my_entity; }
  void set_entity( const stk::mesh::BulkData & mesh, stk::mesh::Entity meshEntity ) const { my_entity = meshEntity; my_entityId = mesh.identifier(my_entity); }

  const ParametricElement & parametric_elem() const { return my_parametric_elem; }
  ParametricElementLibrary & get_parametric_elem_library() const { return myParametricElemLib; }
  stk::topology topology() const { return my_parametric_elem.get_topology(); }

  void fill_node_owner_coords(const Mesh_Element * owner, std::vector<stk::math::Vector3d> & coords) const;
  stk::math::Vector3d compute_local_coords_from_owner_coordinates(const Mesh_Element * owner, const stk::math::Vector3d & coords) const;

  static void integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const ParametricElement & me,
    const ParametricElement & mesh_me);
  static void integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const ParametricElement & me) { integration_weights(intg_weights, numCoordDims, mesh_coords, me, me); }
  static void integration_locations(
    std::vector<stk::math::Vector3d> & intg_pt_locations,
    const ParametricElement & me);
  static void gather_nodal_field(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef field, std::vector<double> & result, unsigned dim = 1, unsigned npe = 0);
  static double volume(ParametricElementLibrary & parametricElemLibrary, const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coords_field);
  static PhaseTag update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const InterfaceID interface_key, const int sign);
  static PhaseTag update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const std::vector<InterfaceID> & interfaces, const std::vector<int> & interfaceSigns);

  const ParametricElement* get_evaluation_parametric_element(const FieldRef field) const;

  void add_subelement(std::unique_ptr<SubElement> subelem);
  // Replace leaves below this by new leaves that fill the same region, as by a swap: the new leaves
  // are children of the lowest subelement (or this) with all the old leaves below it, and a
  // subelement left with no children is removed, so that it does not become a leaf.
  void replace_leaf_subelements(const std::vector<const SubElement *> & oldLeaves,
      std::vector<std::unique_ptr<SubElement>> newLeaves);
  void get_subelements( std::vector<const SubElement *> & subelems ) const;
  void get_subelements( std::vector<SubElement *> & subelems );
  bool have_subelements() const { return !my_subelements.empty(); }
  unsigned num_subelements() const { return my_subelements.size(); }
  bool have_refined_edges() const;

  virtual void determine_decomposed_elem_phase(const CDMesh & /*mesh*/, const std::vector<Surface_Identifier> & /*surfaceIDs*/) { STK_ThrowRequire(false); }
  virtual void cut_interior_intersection_point(CDMesh & mesh, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains);
  virtual void cut_face_intersection_point(const int iFace, const stk::math::Vector3d & pCoords, const std::vector<int> & sortedDomains);
  bool have_edges_with_children() const;
  bool captures_intersection_point_domains(const std::vector<int> & intersectionPointDomains) const;
  bool face_captures_intersection_point_domains(const int iFace, const std::vector<int> & intersectionPointDomains) const;

  ProlongationElementData * get_prolongation_data() const { return my_prolongation_data; }
  void set_prolongation_data( ProlongationElementData * data ) const { my_prolongation_data = data; }

  void prolongate_fields(const CDMesh & mesh) const;

protected:
  ParametricElementLibrary & myParametricElemLib;
  const ParametricElement& my_parametric_elem;
  std::vector<const SubElementNode *> my_nodes;
  std::vector<std::unique_ptr<SubElement>> my_subelements;
  const PhaseTag * my_phase;
  mutable stk::mesh::Entity my_entity;
  mutable stk::mesh::EntityId my_entityId;
  mutable ProlongationElementData * my_prolongation_data;

private:
  //: Default constructor not allowed
  ElementObj();
  bool fill_path_to_leaf(const SubElement * leaf, std::vector<ElementObj *> & path);
  void remove_leaf_subelements(const std::set<const SubElement *> & leaves);
};

class Mesh_Element : public ElementObj {
public:

  Mesh_Element(CDMesh & mesh, stk::mesh::Entity elemEntity);

  virtual ~Mesh_Element();

  static bool is_supported_topology(stk::topology elem_topology) { return determine_subelement_topology(elem_topology).first != stk::topology::INVALID_TOPOLOGY; }
  static std::pair<stk::topology, unsigned> determine_subelement_topology(stk::topology elem_topology);

  stk::math::Vector3d get_node_parametric_coords( const int lnn ) const;
  stk::math::Vector3d get_node_parametric_coords( const SubElementNode * node ) const;

  const ParametricElement & coord_parametric_elem() const { return my_parametric_elem; }

  stk::topology coord_topology() const { return my_parametric_elem.get_topology(); }

  const stk::topology & subelement_topology() const { return my_subelement_topology; }
  int subelement_order() const { return my_subelement_order; }

  std::string visualize(const CDMesh & mesh) const;
  // For debugging: adds the parent (in the group element_<id>_parent) and its leaves to the journal,
  // the leaves in a group per phase (element_<id>_<phase name>), or with a mean ratio below
  // sliverMeanRatio in element_<id>_slivers.  Each leaf's mean ratio (MeanRatioQualityMetric) and
  // phase are in a comment.
  void add_to_journal(CubitJournal & journal,
      const std::function<std::string(const PhaseTag &)> & phaseName,
      const double sliverMeanRatio = 0.) const;
  double interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const;
  std::pair<int, double> interface_edge_crossing_sign_and_position(const InterfaceID interface, const SubElementNode * node1, const SubElementNode * node2) const;
  // The positions along the edge from node1 of the crossings of the interface, if there are some
  // that the signs at its ends do not show (see
  // ElementCutter::fill_multiple_edge_crossing_positions).
  void fill_multiple_edge_crossing_positions(const InterfaceID interface,
      const SubElementNode * node1,
      const SubElementNode * node2,
      std::vector<double> & positions) const;
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
  virtual void determine_decomposed_elem_phase(const CDMesh & mesh, const std::vector<Surface_Identifier> & surfaceIDs) override;
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
  void cut_interior_intersection_points(ParametricElementLibrary & parametricElemLibrary, CDMesh & mesh);
  // Cut the element's subelements at the crossings of the edges that the signs at their ends do
  // not show, so that the sign-based decomposition sees only edges it can cut.
  void cut_multiple_edge_crossings(CDMesh & mesh);
  // The linear element of the vertices: the element itself if it is linear.  A quadratic element is
  // cut as this element (its leaves are linear until build_quadratic_subelements).
  stk::topology vertex_topology() const { return coord_topology().base(); }
  NodeVec vertex_nodes() const { return NodeVec(my_nodes.begin(), my_nodes.begin() + vertex_topology().num_nodes()); }
  std::vector<stk::mesh::Entity> vertex_entities(const stk::mesh::BulkData & mesh) const;
  // Whether the element is cut at its intersection points by rank (insert_intersection_point): a
  // tet, or a tri in 2D, linear or quadratic, cut on its vertices.
  bool is_cut_by_rank() const
  {
    return vertex_topology() == stk::topology::TETRAHEDRON_4 ||
        vertex_topology() == stk::topology::TRIANGLE_3_2D;
  }
  // Cut the element (is_cut_by_rank) at an intersection point in the stencil of its nodes, already
  // reduced to the entity of the element it is on: the point gets a node on the entity of the leaf
  // it is in (inside it, or on a face or edge of it), or the leaf node it is at, with the point's
  // domains added, and every leaf of the element with that entity is split at the node.
  void insert_intersection_point(CDMesh & mesh, const IntersectionPoint & reducedPoint);
  // Split every leaf of the element that has all of entityNodes (all the nodes of a leaf, a face of
  // a tet leaf, or an edge) at node, on the entity they span.
  void split_leaves_at_node(const NodeVec & entityNodes, const SubElementNode * node);
  // The element's sharp-feature corners, each reduced by the edge tolerance onto the entity of the
  // element it is near, in canonical form.
  std::vector<IntersectionPoint> get_corner_intersection_points(const CDMesh & mesh) const;
  // The sharp-edge crossings of the element's leaf faces that no node of the leaf face captures
  // (get_uncaptured_sharp_edge_crossings_of_tet_face), each reduced by the edge tolerance,
  // measured on this element, onto the entity of the element it is near, in canonical form, without
  // those within the edge tolerance of a node of a leaf whose domains capture them.
  std::vector<IntersectionPoint> get_edge_intersection_points(const CDMesh & mesh) const;
  // The crossings of the element's leaf edges that the signs at their ends do not show, as points
  // with the interface's domains, reduced only to roundoff (so each stays on its leaf edge), in
  // canonical form.  Within the edge tolerance, times the extent of the mesh nodes the edge's nodes
  // come from, two neighboring crossings are close crossings, where the edge crosses the surface and
  // back, and are dropped.  A crossing near an end is at that end: that end node, unless it already
  // has the interface's domains.  Near is within the edge tolerance for a surface with sharp
  // features, and within 0.33, times the same extent, for one without.  An edge whose ends both
  // have the interface's domains (a chord of it) has none.  The crossings and the tolerance depend
  // only on the edge, so every element with the edge finds the same points.
  std::vector<IntersectionPoint> get_multiple_edge_crossing_points(const CDMesh & mesh) const;
  // Insert the points in canonical order (insert_intersection_point), without repeats: a point is
  // a repeat of one before it on the same entity whose domains contain its own, and whose weights
  // are within repeatTol of its own.
  void insert_intersection_points(
      CDMesh & mesh, std::vector<IntersectionPoint> points, const double repeatTol);
  // Make the segments of sharp edges inside the tet edges of its leaves, by swaps of its leaves
  // (SharpEdgeRecoveryBySwaps), with the crossings of their faces that the faces do not capture
  // (get_uncaptured_sharp_edge_crossings_of_tet_face).
  void recover_sharp_edges_in_tet_by_swaps(const CDMesh & mesh);
  // The sharp-edge crossings of a face of the tet's leaves that the face does not capture: no node
  // of the face of higher rank has the crossing's domains, and the crossing is not within the edge
  // tolerance, in the face's weights, of a node of the face with them, or of an edge of the face
  // whose two nodes have them (a chord of the sharp edge).
  std::vector<ElementIntersection>
  get_uncaptured_sharp_edge_crossings_of_tet_face(const NodeVec & faceNodes, const double edgeTol) const;

  // Set each leaf's interface signs from the crossings in it (get_interface_signs_based_on_crossings).
  void reset_leaf_interface_signs();

  void find_child_coordinates_at_owner_coordinates(const stk::math::Vector3d & ownerCoordinates, const ElementObj *& child, stk::math::Vector3d & child_p_coords) const;

  bool is_single_coincident() const;

private:
  Mesh_Element() = delete;
  stk::math::Vector3d get_intersection_point_parametric_coordinates(const IntersectionPoint & intersectionPoint) const;
  std::vector<double> gather_nodal_coordinates() const;
  // The nodes of each leaf, or of the vertices if the element has no subelements (it is its own leaf).
  std::vector<NodeVec> leaf_node_lists() const;
  // The sharp-edge crossings of the unique faces of the leaves that the faces do not capture
  // (get_uncaptured_sharp_edge_crossings_of_tet_face), as points in the stencil of the vertices.
  std::vector<IntersectionPoint> sharp_edge_crossing_points_of_leaf_faces(const std::vector<NodeVec> & leavesNodes,
      const std::vector<stk::mesh::Entity> & elementNodes,
      const double edgeTol) const;
  // The leaf nodes with domains, each as a point in the stencil of the vertices, with its domains,
  // reduced to roundoff, in canonical form.
  std::vector<IntersectionPoint> points_of_leaf_nodes_with_domains(const stk::mesh::BulkData & stkMesh,
      const std::vector<NodeVec> & leavesNodes,
      const std::vector<stk::mesh::Entity> & elementNodes) const;
  // The multiple crossings of the leaf edge to cut (multiple_crossing_positions_to_insert), of each
  // interface it is not a chord of, appended as points in the stencil of the vertices.  The end
  // tolerance is that of a surface with sharp features only if they are cut (sharpFeaturesAreCut).
  void append_multiple_crossing_points_of_leaf_edge(const SubElementNode * node0,
      const SubElementNode * node1,
      const std::vector<stk::mesh::Entity> & elementNodes,
      const double edgeTol,
      const bool sharpFeaturesAreCut,
      std::vector<IntersectionPoint> & crossings) const;
  void replace_tet_leaves(const std::vector<const SubElement *> & oldLeaves,
      const std::vector<std::array<const SubElementNode *, 4>> & tets);

  int my_subelement_order;
  stk::topology my_subelement_topology;
  bool my_have_interface;

  std::unique_ptr<ElementCutter> myCutter;
  std::vector<InterfaceID> myCuttingInterfaces;
};

} // namespace krino

#endif // Akri_Element_h
