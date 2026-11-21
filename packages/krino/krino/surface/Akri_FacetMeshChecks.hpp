// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETMESHCHECKS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETMESHCHECKS_HPP_

#include <Akri_BoundingBox.hpp>
#include <stk_math/StkVector.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace krino
{

class FacetMesh;

// Checks of a FacetMesh's owned facets for remove_coincident_facets and check_watertight, and its
// sharp features.  Each has
// two forms: by vertex ids, exact when facets that meet share their vertices (merged vertices), and
// by location, for unconnected facets (a soup), which does not need that.  All are
// collective.  Removing facets drops the ghosts (see FacetMesh).

// Among owned facets with the same vertices, on any processors, those of opposite orientations
// cancel in pairs and the remaining ones of the same orientation collapse to one, as
// remove_coincident_facets does.  Returns the global number removed.
size_t remove_coincident_facets_by_vertex_ids(FacetMesh & mesh);
// Among owned facets at the same location, on any processors, as remove_coincident_facets finds
// them, those removed are the ones remove_coincident_facets would remove from all of them together,
// whatever processors own them.  Returns the global number removed.
size_t remove_coincident_facets_by_location(FacetMesh & mesh);

// The edges (vertices in 2D) of owned facets not shared by exactly two owned facets on any
// processors, matched by vertex ids, each reported by one processor as its two endpoints,
// restricted as is_reportable_non_watertight_edge does.
std::vector<std::array<stk::math::Vector3d, 2>> find_non_watertight_edges_by_vertex_ids(
    const FacetMesh & mesh, const std::optional<BoundingBox> & restrictBbox);
// The same, with edges matched by location, as find_non_watertight_edges does for disconnected
// facets.
std::vector<std::array<stk::math::Vector3d, 2>> find_non_watertight_edges_by_location(
    const FacetMesh & mesh, const std::optional<BoundingBox> & restrictBbox);

// The sharp features of a surface: in 3D its feature edges and corners, in 2D its corners.
struct FacetMeshSharpFeatures
{
  struct Edge
  {
    std::array<std::uint64_t, 2> vertexIds; // sorted
    std::array<stk::math::Vector3d, 2> coords;
  };
  struct Corner
  {
    std::uint64_t vertexId;
    stk::math::Vector3d coords;
  };
  std::vector<Edge> edges;
  std::vector<Corner> corners;
};

// Collective.  The sharp features of the owned facets on all processors, matched by vertex ids (so
// it needs shared vertices), each reported by one processor.  The angle between two facets, and the
// feature angle whose cosine is cosFeatureAngle, are as in Akri_SharpFeatureCalcs.
// In 3D an edge is a feature edge where its two facets meet at an angle sharper than the feature
// angle, where more than two facets meet (a non-manifold edge), or, if
// areOpenBoundaryEdgesFeatures, where only one facet has it.  A vertex is a corner where one
// feature edge ends, where more than two meet, or where two meet at an angle sharper than the
// feature angle. In 2D a vertex is a corner where its two segments meet at an angle sharper than
// the feature angle, where more than two meet, or, if areOpenBoundaryEdgesFeatures, where only one
// has it.
FacetMeshSharpFeatures find_sharp_features_by_vertex_ids(const FacetMesh & mesh,
    const double cosFeatureAngle,
    const bool areOpenBoundaryEdgesFeatures = false);

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_FACETMESHCHECKS_HPP_ */
