#ifndef KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_LEVELSETSHAPESENSITIVITYTYPES_HPP_
#define KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_LEVELSETSHAPESENSITIVITYTYPES_HPP_

// The value types of the level set shape sensitivities, without the rest of krino, so that users
// of the results (such as PlatoKrinoInterface's consumers) need only STK and this header.

#include <array>
#include <cstddef>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino
{

// How the shape sensitivities, the derivatives of the interface node coordinates with respect to
// the nodal level set values, are computed.  Each name gives the direction the interface node moves
// and the gradient that sets its speed.
//
// EDGE_DIRECTION_USING_EDGE_GRADIENT
//   The node moves along the mesh edge it lies on: dx/dphi = (x1 - x0) / (phi0 - phi1).  Snapped
//   nodes use the edge they were snapped along, with the unsnapped coordinates and the original
//   level set.  This is the exact derivative of the discrete interface, so it agrees with finite
//   differences of any objective computed on the interface facets, with one exception: without
//   snapping, a background node on the interface (zero level set) has no cut edge.  Each of its
//   sides moves along an edge of its own parent element, so no single velocity is exact for every
//   objective.  The node moves along the surface normal N (the sum of its sides' area vectors) at
//   the speed that makes the enclosed volume's derivative exact (see
//   Akri_InterfaceInverseGradientArea.hpp): dx/dphi = -sum_s(A_s / |grad phi_s|) N / |N|^2.  On
//   the domain boundary, that velocity follows the interface's contact with the boundary, keeping
//   its normal speed: within the wall that contains the contact, or along the line where the two
//   walls that contain it meet (see compute_zero_levelset_boundary_node_constraint; walls are
//   separated by side parts and sharp edges).  That makes the area exact for a planar interface on
//   flat walls, including at edges and corners; away from the boundary and for curved interfaces it
//   is not.  With snapping, predictable snap selection moves such a node along one allowed edge, so
//   it is exact like any other snapped node.  At nodes with zero level set these are the
//   derivatives for increasing level set.
//
// SURFACE_NORMAL_USING_EDGE_GRADIENT
//   The edge direction velocity projected onto N, so the node moves normal to the discrete
//   interface at the edge form's normal speed.  N is the direction of the derivative of the
//   enclosed volume, so volume sensitivities equal the edge direction ones; area and other
//   objectives are not exact.  Background nodes on the interface with no cut edge move along N
//   (the boundary restriction above does not change their normal speed).
enum class ShapeSensitivityMethod
{
  EDGE_DIRECTION_USING_EDGE_GRADIENT,
  SURFACE_NORMAL_USING_EDGE_GRADIENT
};

// dCoordsdParentLevelSets[i] = d(interface node coordinates) / d(level set at parentNodeIds[i]),
// for 2 parents (one level set) to 4 (multiple level sets).  Includes shared interface nodes.
struct LevelSetShapeSensitivity
{
  LevelSetShapeSensitivity(const stk::mesh::EntityId inInterfaceNodeId,
      const std::vector<stk::mesh::EntityId> & inParentNodeIds,
      const std::vector<stk::math::Vector3d> & inDCoordsdParentLevelSets)
      : interfaceNodeId(inInterfaceNodeId),
        parentNodeIds(inParentNodeIds),
        dCoordsdParentLevelSets(inDCoordsdParentLevelSets)
  {
  }
  stk::mesh::EntityId interfaceNodeId;
  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;
};

// The interface facets (sides) of one level set's negative phase, oriented out of it, within its
// decomposition's interface restriction if any, with the shape sensitivities of their nodes:
// facetSensitivityIndices[f][i] indexes sensitivities, which are sorted by interface node id.
// NNODES is the spatial dimension.
template <size_t NNODES> struct InterfaceFacetsWithSensitivities
{
  std::vector<LevelSetShapeSensitivity> sensitivities;
  std::vector<std::array<size_t, NNODES>> facetSensitivityIndices;

  size_t num_facets() const { return facetSensitivityIndices.size(); }
  std::array<const LevelSetShapeSensitivity *, NNODES> node_sensitivities(const size_t facet) const;
  std::array<stk::mesh::Entity, NNODES>
  nodes(const stk::mesh::BulkData & mesh, const size_t facet) const;
};

} // namespace krino

#endif /* KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_LEVELSETSHAPESENSITIVITYTYPES_HPP_ */
