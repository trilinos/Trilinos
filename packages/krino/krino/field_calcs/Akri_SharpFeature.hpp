#ifndef KRINO_KRINO_FIELD_CALCS_AKRI_SHARPFEATURE_HPP_
#define KRINO_KRINO_FIELD_CALCS_AKRI_SHARPFEATURE_HPP_
#include <array>
#include <map>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_FieldRef.hpp>

namespace krino {

struct Edge;

class SharpFeatureConstraint
{
public:
  bool is_pinned() const { return myConstrainedEdgeNeighbors[0] == invalid_entity() && myConstrainedEdgeNeighbors[1] == invalid_entity(); }
  bool is_constrained_on_edge() const { return myConstrainedEdgeNeighbors[0] != invalid_entity() && myConstrainedEdgeNeighbors[1] != invalid_entity(); }
  const std::array<stk::mesh::Entity,2> & get_sharp_edge_nodes() const { STK_ThrowAssert(is_constrained_on_edge()); return myConstrainedEdgeNeighbors; }
  static SharpFeatureConstraint edge_constraint(const stk::mesh::Entity entity0, const stk::mesh::Entity entity1) { return SharpFeatureConstraint{entity0, entity1}; }
  static SharpFeatureConstraint pinned_constraint() { return SharpFeatureConstraint(invalid_entity(),invalid_entity()); }
private:
  static stk::mesh::Entity invalid_entity() { static const stk::mesh::Entity invalidEntity; return invalidEntity; }
  SharpFeatureConstraint(const stk::mesh::Entity entity0, const stk::mesh::Entity entity1) : myConstrainedEdgeNeighbors{entity0, entity1} {}
  std::array<stk::mesh::Entity,2> myConstrainedEdgeNeighbors;
};

class SharpFeatureInfo
{
public:
  // The coordinates of a node, less an optional displacement.
  struct NodeCoordinates
  {
    FieldRef coordsField;
    FieldRef displacements;
    stk::math::Vector3d operator()(const stk::mesh::Entity node, const unsigned dim) const;
  };

  // The sharp edges are those between two sides in sideSelector, each with an element in
  // elementSelector, at a dihedral angle sharper than the feature angle.  Constraints are found for
  // the locally owned nodes, and also the shared ones if includeSharedNodes (which needs the aura).
  // If coordsDisplacements is valid, it is subtracted from the coordinates (for example to use the
  // unsnapped coordinates).
  void find_sharp_features(const stk::mesh::BulkData & mesh,
      const FieldRef coordsField,
      const stk::mesh::Selector & elementSelector,
      const stk::mesh::Selector & sideSelector,
      const double cosFeatureAngle,
      const bool includeSharedNodes = false,
      const FieldRef coordsDisplacements = FieldRef());
  const SharpFeatureConstraint * get_constraint(const stk::mesh::Entity node) const;
  const std::map<stk::mesh::Entity,SharpFeatureConstraint> & get_constraints() const { return myNodeToConstrainedNeighbors; }
private:
  void find_sharp_features_2D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const stk::mesh::Selector & elementSelector,
      const stk::mesh::Selector & sideSelector,
      const double cosFeatureAngle,
      const bool includeSharedNodes);
  void find_sharp_features_3D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const stk::mesh::Selector & elementSelector,
      const stk::mesh::Selector & sideSelector,
      const double cosFeatureAngle,
      const bool includeSharedNodes);
  static bool edge_has_sharp_feature_3D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const stk::mesh::Selector & elementSelector,
      const stk::mesh::Selector & sideSelector,
      const double cosFeatureAngle,
      const Edge edge);
  static bool node_has_sharp_feature_2D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const stk::mesh::Selector & elementSelector,
      const stk::mesh::Selector & sideSelector,
      const double cosFeatureAngle,
      const stk::mesh::Entity node);
  static bool angle_is_sharp_between_any_two_sides_3D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const double cosFeatureAngle,
      const std::array<stk::mesh::Entity, 2> & edgeNodes,
      const std::vector<stk::mesh::Entity> & sidesOfEdge);
  static bool angle_is_sharp_between_any_two_sides_2D(const stk::mesh::BulkData & mesh,
      const NodeCoordinates & coords,
      const double cosFeatureAngle,
      const stk::mesh::Entity node,
      const std::vector<stk::mesh::Entity> & sidesOfEdge);
  std::map<stk::mesh::Entity,SharpFeatureConstraint> myNodeToConstrainedNeighbors;
};

bool is_intersection_point_node_compatible_for_snapping_based_on_sharp_features(const SharpFeatureInfo & sharpFeatureInfo, const stk::mesh::Entity intPtNode, const std::vector<stk::mesh::Entity> & intPtNodes);

// Whether two tri sides sharing the edge meet at a dihedral angle sharper than the feature angle
// (its cosine greater than cosFeatureAngle), as find_sharp_features decides.
bool is_sharp_between_tri_sides(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::array<stk::mesh::Entity, 2> & edgeNodes,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1,
    const double cosFeatureAngle);

// Whether two line sides sharing the node meet at an angle sharper than the feature angle (2D).
bool is_sharp_between_line_sides(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Entity node,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1,
    const double cosFeatureAngle);

}

#endif /* KRINO_KRINO_FIELD_CALCS_AKRI_SHARPFEATURE_HPP_ */
