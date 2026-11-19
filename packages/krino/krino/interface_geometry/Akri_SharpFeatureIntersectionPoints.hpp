// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_SHARPFEATUREINTERSECTIONPOINTS_HPP_
#define KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_SHARPFEATUREINTERSECTIONPOINTS_HPP_

#include <vector>

#include <Akri_Intersection_Points.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace krino {

// The face and interior points of sharp features, by rank, in which a point within the snap
// tolerance of a lower-dimensional entity is moved onto it.

// Drop the nodes of the stencil of an intersection point whose weights are below tol (keeping the
// node with the largest weight if all are), and renormalize the weights of the rest.  This moves a
// point that is within tol of a lower-dimensional entity onto it: an interior point near a face
// becomes a face point, a face point near an edge an edge point, and a point near a node is left
// with that node alone.
void reduce_intersection_point_stencil(
    std::vector<stk::mesh::Entity> & nodes, std::vector<double> & weights, const double tol);

// The weights of the vertices of a simplex (a tet in 3D, a triangle in 2D) at these parametric
// coordinates, reduced by tol as by reduce_intersection_point_stencil: the ordinals of the vertices
// kept, and their weights.  So the vertices kept are those of the entity of the simplex (itself, a
// face, an edge, or a vertex) that the point is within tol of.
void fill_reduced_simplex_vertex_weights(const stk::math::Vector3d & paramCoords,
    const unsigned dim,
    const double tol,
    std::vector<int> & vertexOrdinals,
    std::vector<double> & weights);

// The codimension of the feature that a point with these sorted domains is on.  For one surface
// (all domains the same) it is the number of domains: {i} a crossing of the smooth surface, {i,i} a
// sharp edge, {i,i,i} a corner.  For level sets per phase (distinct domains) it is one less: {i,j}
// an interface, {i,j,k} a triple line, {i,j,k,l} a quadruple point.
int sharp_feature_codimension(const std::vector<int> & sortedDomains);

// Whether a node with these sorted domains captures a point with these: the node's domains contain
// the point's, counting repeats.  So a corner {i,i,i} captures the crossings {i,i} of its sharp
// edges, and a node on the smooth surface {i} does not.
bool node_domains_capture_point_domains(
    const std::vector<int> & sortedNodeDomains, const std::vector<int> & sortedPointDomains);
inline bool is_corner_intersection_point(const IntersectionPoint & intersectionPoint)
{
  return 3 == sharp_feature_codimension(intersectionPoint.get_sorted_domains());
}
inline bool is_sharp_edge_intersection_point(const IntersectionPoint & intersectionPoint)
{
  return 2 == sharp_feature_codimension(intersectionPoint.get_sorted_domains());
}

// Whether the points have the same domains and the same stencil nodes, in any order, with the
// weights of each node within tol.
bool is_same_intersection_point(
    const IntersectionPoint & a, const IntersectionPoint & b, const double tol);

// Whether point is at capturing, as the same point or one of its lower multiplicity: the same
// stencil nodes, in any order, with the weights of each node within tol, and domains that
// capturing's contain (as a multiset).  So a corner of a surface, with domains {0,0,0}, captures
// the crossings of its edges
// ({0,0}) as near it.
bool is_captured_by_intersection_point(
    const IntersectionPoint & capturing, const IntersectionPoint & point, const double tol);

// The point with its stencil reduced by tol and its stencil nodes sorted by global id.
IntersectionPoint reduced_intersection_point_in_canonical_form(
    const stk::mesh::BulkData & mesh, const IntersectionPoint & intersectionPoint, const double tol);

// Sort points in canonical form into an order that depends only on the model coordinates of their
// stencil nodes, their weights and domains, the same on every build and for any parallel
// decomposition: interior points, then face, edge, and single-node (captured) points; each by its
// stencil nodes' coordinates (in single precision), then by more domains first (so a point comes
// before the ones it can capture), then domains, then weights.  Global ids, which refinement assigns
// differently for different decompositions and builds, only break exact ties.
void sort_intersection_points_in_canonical_order(
    const stk::mesh::BulkData & mesh, std::vector<IntersectionPoint> & intersectionPoints);

// Remove the points that an earlier one captures (is_captured_by_intersection_point).  The points
// must be in canonical order, so that a point comes before the ones it can capture.
void remove_repeated_intersection_points(std::vector<IntersectionPoint> & intersectionPoints, const double tol);

// The intersection points of one codimension (3: corners; 2: sharp-edge crossings; 1: crossings
// of a surface), each reduced by tol to the entity it is within tol of, in canonical order, without
// repeats of the same point (to roundoff): distinct points within tol of each other are all kept.
// The points must be in the stencil of the parent element, so that tol is measured on it.
std::vector<IntersectionPoint>
reduced_intersection_points_of_codimension(const stk::mesh::BulkData & mesh,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const int codimension,
    const double tol);
}

#endif /* KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_SHARPFEATUREINTERSECTIONPOINTS_HPP_ */
