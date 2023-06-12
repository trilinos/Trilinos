#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SHARPFEATURE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SHARPFEATURE_HPP_
#include <map>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "Akri_FieldRef.hpp"

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
  void find_sharp_features(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const double cosFeatureAngle);
  const SharpFeatureConstraint * get_constraint(const stk::mesh::Entity node) const;
private:
  void find_sharp_features_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle);
  void find_sharp_features_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle);
  static bool edge_has_sharp_feature_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle, const Edge edge);
  static bool node_has_sharp_feature_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle, const stk::mesh::Entity node );
  static bool angle_is_sharp_between_any_two_sides_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double cosFeatureAngle, const std::array<stk::mesh::Entity,2> & edgeNodes, const std::vector<stk::mesh::Entity> & sidesOfEdge);
  static bool angle_is_sharp_between_any_two_sides_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double cosFeatureAngle, const stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & sidesOfEdge);
  std::map<stk::mesh::Entity,SharpFeatureConstraint> myNodeToConstrainedNeighbors;
};

bool is_intersection_point_node_compatible_for_snapping_based_on_sharp_features(const SharpFeatureInfo & sharpFeatureInfo, const stk::mesh::Entity intPtNode, const std::vector<stk::mesh::Entity> & intPtNodes);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SHARPFEATURE_HPP_ */
