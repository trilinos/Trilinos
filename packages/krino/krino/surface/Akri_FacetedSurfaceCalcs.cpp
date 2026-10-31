/*
 * Akri_FacetedSurfaceCalcs.cpp
 *
 *  Created on: May 22, 2023
 *      Author: drnoble
 */
#include <Akri_BoundingBoxDistance.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>

#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Sign.hpp>
#include <Akri_QualityMetric.hpp>
#include <krino/mesh_utils/Akri_AllReduce.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace krino {

static double compute_bbox_padding(const std::vector<BoundingBox> & procFacetBboxes, const BoundingBox & localQueryBbox, const double truncationLength)
{
  double bboxPadding = std::numeric_limits<double>::max();
  for ( auto && procFacetBBox : procFacetBboxes )
  {
    if (procFacetBBox.valid())
    {
      const double upperBnd = std::sqrt(max_possible_closest_squared_distance_between_contained_points(procFacetBBox, localQueryBbox));
      bboxPadding = std::min(bboxPadding,upperBnd);
    }
  }
  if (std::numeric_limits<double>::max() == bboxPadding)
  {
    bboxPadding = truncationLength; // only should happen for no facets anywhere
  }
  if (truncationLength > 0.0)
  {
    bboxPadding = std::min(bboxPadding, truncationLength);
  }
  return bboxPadding;
}

std::vector<BoundingBox> fill_processor_bounding_boxes(const BoundingBox & localFacetBbox, const BoundingBox & localQueryBbox, const double truncationLength)
{
  std::vector<BoundingBox> procFacetBboxes;
  BoundingBox::gather_bboxes( localFacetBbox, procFacetBboxes );

  const double bboxPadding = compute_bbox_padding(procFacetBboxes, localQueryBbox, truncationLength);

  BoundingBox localBBox = localQueryBbox;
  localBBox.pad(bboxPadding);

  std::vector<BoundingBox> procPaddedQueryBboxes;
  BoundingBox::gather_bboxes( localBBox, procPaddedQueryBboxes );

  return procPaddedQueryBboxes;
}

template<class FACET>
static double compute_point_distance_squared(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  double minSqrDist = std::numeric_limits<double>::max();
  for ( auto&& facet : nearestFacets )
  {
    const double sqrDist = facet->point_distance_squared(x);
    if (sqrDist < minSqrDist)
      minSqrDist = sqrDist;
  }
  return minSqrDist;
}



template <class FACET>
double point_distance_given_nearest_facets(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance)
{
  if (nearestFacets.empty())
  {
    return far_field_value;
  }

  double dist = 0.0;
  if (compute_signed_distance)
  {
    dist = compute_point_to_facets_distance_by_average_normal(x, nearestFacets);
    if (0.0 != narrow_band_size && std::abs(dist) > narrow_band_size)
    {
      dist = far_field_value;
    }
  }
  else
  {
    const double minSqrDist = compute_point_distance_squared(x, nearestFacets);
    if (0.0 != narrow_band_size && minSqrDist > narrow_band_size*narrow_band_size)
    {
      dist = far_field_value;
    }
    else
    {
      dist = std::sqrt(minSqrDist);
    }
  }

  return dist;
}

template<class FACET>
std::vector<FacetDistanceQuery<FACET>> build_distance_queries(const stk::math::Vector3d &x, const std::vector<const FACET*> & facets)
{
  std::vector<FacetDistanceQuery<FACET>> facetDistQueries;
  facetDistQueries.reserve(facets.size());
  for ( auto&& facet : facets )
  {
    if (facet->degenerate()) continue; // Skip zero-sized facets
    facetDistQueries.emplace_back(*facet, x);
  }
  return facetDistQueries;
}

template<class FACET>
bool wins_tie_breaker(const FACET & a, const FACET & b)
{
  return is_less_than_in_x_then_y_then_z(a.centroid(), b.centroid());
}

template<class FACET>
bool is_closer(const double sqrDistA, const FACET * a, const double sqrDistB, const FACET * b)
{
  // The tie breaker in the test below is to guarantee that the resulting closest facet is independent of the order
  // in which the facets are stored/traversed.
  return sqrDistA < sqrDistB ||
    (sqrDistA == sqrDistB && is_less_than_in_x_then_y_then_z(a->centroid(), b->centroid()));
}

template<class FACET>
bool is_closer(const FacetDistanceQuery<FACET> & a, const FacetDistanceQuery<FACET> & b)
{
  return is_closer(a.distance_squared(), &a.facet(), b.distance_squared(), &b.facet());
}

template<class FACET>
stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  double minSqrDist = std::numeric_limits<double>::max();
  const FACET* closestFacet = nullptr;
  stk::math::Vector3d closestFacetPt;
  stk::math::Vector3d facetPt;
  for ( auto&& facet : nearestFacets )
  {
    facet->closest_point(x, facetPt);
    const double sqrDist = (x-facetPt).length_squared();
    if (is_closer(sqrDist, facet, minSqrDist, closestFacet))
    {
      closestFacet = facet;
      minSqrDist = sqrDist;
      closestFacetPt = facetPt;
    }
  }
  return closestFacetPt;
}

template<class FACET>
const FacetDistanceQuery<FACET> & find_closest_facet(const std::vector<FacetDistanceQuery<FACET>> & facetDistQueries)
{
  // The tie breaker in the test below is to guarantee that the resulting closest facet is independent of the order
  // in which the queries are stored.
  unsigned closestIndex = 0;
  for ( unsigned index=1; index<facetDistQueries.size(); ++index )
    if ( is_closer(facetDistQueries[index], facetDistQueries[closestIndex]) )
      closestIndex = index;
  return facetDistQueries[closestIndex];
}

template<class FACET>
stk::math::Vector3d compute_pseudo_normal(const std::vector<FacetDistanceQuery<FACET>> & facet_queries, const FacetDistanceQuery<FACET> & closestQuery)
{
  constexpr double tol = 1.e-6;
  stk::math::Vector3d pseudo_normal = stk::math::Vector3d::ZERO;
  stk::math::Vector3d average_normal = stk::math::Vector3d::ZERO;

  // NOTE: Because the facets are not sorted here, it is possible that different results can be obtained in the last
  // digits on procs with different facet order.  Rather than sort the facets, this is managed by using a very
  // small tolerance to check if the query point is on the nearest facet (sqrTolSoCloseThatItIsOnClosestFacet below).
  // It is expected that this should make it unnecessary for the pseudo_normal come out to exactly the same value
  // on all procs.

  const stk::math::Vector3d nearest_closest_point = closestQuery.closest_point();
  const double nearest_size2 = closestQuery.facet().mean_squared_edge_length();
  unsigned close_count = 0;
  for ( auto&& query : facet_queries )
  {
    const stk::math::Vector3d closest_point = query.closest_point();
    const double dist2_from_nearest = (closest_point-nearest_closest_point).length_squared();
    if (dist2_from_nearest < tol*tol*nearest_size2)
    {
      ++close_count;
      const FACET & facet = query.facet();

      average_normal += facet.facet_normal();

      if constexpr (3 == FACET::DIM)
      {
        const stk::math::Vector3d closest_pt_wts = query.closest_point_weights();
        const int closest_node = (closest_pt_wts[0] > closest_pt_wts[1]) ? ((closest_pt_wts[0] > closest_pt_wts[2]) ? 0 : 2) : ((closest_pt_wts[1] > closest_pt_wts[2]) ? 1 : 2);

        const int n0 = closest_node;
        const int n1 = (n0<2) ? (n0+1) : 0;
        const int n2 = (n1<2) ? (n1+1) : 0;
        const stk::math::Vector3d edge0 = facet.facet_vertex(n1) - facet.facet_vertex(n0);
        const stk::math::Vector3d edge1 = facet.facet_vertex(n2) - facet.facet_vertex(n0);
        const double facet_angle = std::acos(Dot(edge0, edge1)/(edge0.length()*edge1.length()));

        pseudo_normal += facet.facet_normal()*facet_angle;
      }
    }
  }
  STK_ThrowRequireMsg(close_count>0,"Issue with tolerance in compute_pseudo_normal.  No facet found within tolerance of closest point.");

  return (3 == FACET::DIM && close_count > 2) ? pseudo_normal : average_normal;
}

static bool is_closest_point_on_facet_edge_or_vertex(const stk::math::Vector3d & closestFacetWeights, const int dim)
{
  // If the closest_point weights are all larger than this value, then the closest point
  // is considered to be on the face of the closest facet rather than on the edges of the facet.
  constexpr double edgeTol = 1.e-6;

  for (int d=0; d<dim; ++d)
    if (closestFacetWeights[d] < edgeTol)
      return true;
  return false;
}

template<class FACET>
double compute_point_to_facets_distance_by_average_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & facets)
{
  const std::vector<FacetDistanceQuery<FACET>> facetDistQueries = build_distance_queries(x, facets);
  STK_ThrowRequireMsg(!facetDistQueries.empty(), "All facets are degenerate in compute_point_to_facets_distance_by_average_normal.");

  const FacetDistanceQuery<FACET> & closestQuery = find_closest_facet(facetDistQueries);
  constexpr double sqrTolSoCloseThatItIsOnClosestFacet = 1.e-18;

  if ( closestQuery.distance_squared()/closestQuery.facet().mean_squared_edge_length() < sqrTolSoCloseThatItIsOnClosestFacet )
  {
    return 0.0;
  }

  // If closest point is on the facet of the closest facet then only that facet is considered in the distance calculation.
  // Otherwise, all of the facets are considered to compute an average normal in order to compute the distance.

  if (!is_closest_point_on_facet_edge_or_vertex(closestQuery.closest_point_weights(), FACET::DIM))
  {
    return closestQuery.signed_distance(x);
  }

  const double min_sqr_dist = closestQuery.distance_squared();
  const stk::math::Vector3d pseudo_normal = compute_pseudo_normal(facetDistQueries, closestQuery);

  if (pseudo_normal.length_squared() == 0.0)
  {
    krinolog << "Warning:  Cannot determine the average facet normal for computing the level set distance at point " << x
        << ".  This can happen when faceted facet include coincident facets with opposite normals.  Arbitrarily setting the distance to be positive." << stk::diag::dendl;
    return std::sqrt(min_sqr_dist);
  }
  else
  {
    if (Dot(pseudo_normal, x-closestQuery.closest_point()) > 0)
    {
      return std::sqrt(min_sqr_dist);
    }
    else
    {
      return -std::sqrt(min_sqr_dist);
    }
  }
}

template<class FACET>
stk::math::Vector3d compute_closest_point_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  const std::vector<FacetDistanceQuery<FACET>> facetDistQueries = build_distance_queries(x, nearestFacets);
  STK_ThrowRequireMsg(!facetDistQueries.empty(), "All facets are degenerate in compute_point_to_facets_distance_by_average_normal.");

  const FacetDistanceQuery<FACET> & closestQuery = find_closest_facet(facetDistQueries);

  if (!is_closest_point_on_facet_edge_or_vertex(closestQuery.closest_point_weights(), FACET::DIM))
  {
    return closestQuery.facet().facet_normal();
  }

  const stk::math::Vector3d pseudoNormal = compute_pseudo_normal(facetDistQueries, closestQuery);

  if (pseudoNormal.length_squared() == 0.0)
  {
    krinolog << "Warning:  Cannot determine the average facet normal for computing the level set distance at point " << x
        << ".  This can happen when faceted facet include coincident facets with opposite normals.  The normal may be wrongly oriented." << stk::diag::dendl;
    return closestQuery.facet().facet_normal();
  }

  constexpr double sqrTolSoCloseThatItIsOnSurface = 1.e-12;
  if ( closestQuery.distance_squared()/closestQuery.facet().mean_squared_edge_length() < sqrTolSoCloseThatItIsOnSurface )
    return pseudoNormal.unit_vector();

  stk::math::Vector3d unitDiff = x-closestQuery.closest_point();
  unitDiff.unitize();

  if (Dot(pseudoNormal, x-closestQuery.closest_point()) > 0)
    return unitDiff;
  else
    return -unitDiff;
}

bool is_projection_of_point_inside_enlarged_triangle(const stk::math::Vector3d & triPt0, const stk::math::Vector3d & triPt1, const stk::math::Vector3d & triPt2, const stk::math::Vector3d& p)
{
  constexpr double expand {1e-10};
  const stk::math::Vector3d centroid = 1./3.*(triPt0+triPt1+triPt2);
  const stk::math::Vector3d p0 = triPt0 + expand*(triPt0-centroid);
  const stk::math::Vector3d p1 = triPt1 + expand*(triPt1-centroid);
  const stk::math::Vector3d p2 = triPt2 + expand*(triPt2-centroid);
  return Facet3d::Calc::is_projection_of_point_inside_triangle(p0, p1, p2, p);
}

bool is_projection_of_point_inside_enlarged_segment(const stk::math::Vector3d & segPt0, const stk::math::Vector3d & segPt1, const stk::math::Vector3d& p)
{
  constexpr double expand {1e-10};
  const stk::math::Vector3d p0 = segPt0 + expand*(segPt0-segPt1);
  const stk::math::Vector3d p1 = segPt1 + expand*(segPt1-segPt0);

  return Facet2d::Calc::is_projection_of_point_inside_segment(p0, p1, p);
}

bool is_projection_of_point_inside_enlarged_facet(const Facet3d & facet, const stk::math::Vector3d& p)
{
  return is_projection_of_point_inside_enlarged_triangle(facet.facet_vertex(0), facet.facet_vertex(1), facet.facet_vertex(2), p);
}

bool is_projection_of_point_inside_enlarged_facet(const Facet2d & facet, const stk::math::Vector3d& p)
{
  return is_projection_of_point_inside_enlarged_segment(facet.facet_vertex(0), facet.facet_vertex(1), p);
}

template<class FACET>
std::pair<int, double> compute_facet_edge_intersection(const FACET & facet,
  const stk::math::Vector3d& edgePt0,
  const stk::math::Vector3d& edgePt1)
{
  const double dist0 = facet.facet_plane_signed_distance(edgePt0);
  const double dist1 = facet.facet_plane_signed_distance(edgePt1);

  if (sign_change(dist0, dist1))
  {
    const double loc = dist0 / (dist0-dist1);
    const stk::math::Vector3d ptLoc = (1.-loc)*edgePt0 + loc*edgePt1;
    if (is_projection_of_point_inside_enlarged_facet(facet, ptLoc))
      return {sign(dist1), loc};
  }
  return {0, -1.};
}

template<class FACET>
double compute_intersection_between_surface_facets_and_edge(const std::vector<const FACET*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1)
{
  if (candidates.empty())
    return -1.;

  bool haveCrossing = false;
  double intersectionLoc = -1.;

  for (const FACET * surfFacet : candidates)
  {
    const auto [facetCrossingSign, facetIntersectionLoc] = compute_facet_edge_intersection(*surfFacet, edgePt0, edgePt1);
    if (facetCrossingSign != 0)
    {
      if (!haveCrossing || std::abs(facetIntersectionLoc-0.5) < std::abs(intersectionLoc-0.5)) // pick intersection closest to middle of edge
        intersectionLoc = facetIntersectionLoc;
      haveCrossing = true;
    }
  }

  return intersectionLoc;
}

template<class FACET>
stk::math::Vector3d compute_pseudo_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  const std::vector<FacetDistanceQuery<FACET>> facetDistQueries = build_distance_queries(x, nearestFacets);
  STK_ThrowRequireMsg(!facetDistQueries.empty(), "All facets are degenerate in compute_pseudo_normal.");

  const FacetDistanceQuery<FACET> & closestQuery = find_closest_facet(facetDistQueries);

  return compute_pseudo_normal(facetDistQueries, closestQuery);
}

template<class FACET>
static std::vector<long long> canonical_quantized_facet_key(const FACET & facet)
{
  // Winding-independent key: quantize each vertex to a facet-size-scaled grid, then
  // sort the vertices, so a facet and its reverse-wound twin produce the same key.
  constexpr double relTol = 1.e-6;
  const double scale = std::sqrt(facet.mean_squared_edge_length());
  const double cell = relTol * scale;

  std::vector<std::array<long long,3>> quantizedVerts;
  quantizedVerts.reserve(FACET::DIM);
  for (int v=0; v<FACET::DIM; ++v)
  {
    const stk::math::Vector3d & p = facet.facet_vertex(v);
    std::array<long long,3> q;
    for (int d=0; d<3; ++d)
      q[d] = (cell > 0.) ? std::llround(p[d]/cell) : 0ll;
    quantizedVerts.push_back(q);
  }
  std::sort(quantizedVerts.begin(), quantizedVerts.end());

  std::vector<long long> key;
  key.reserve(3*FACET::DIM);
  for (auto && q : quantizedVerts)
    for (int d=0; d<3; ++d)
      key.push_back(q[d]);
  return key;
}

template<class FACET>
size_t remove_coincident_facets(std::vector<FACET> & facets)
{
  const size_t numFacets = facets.size();
  if (numFacets < 2)
    return 0;

  // Group facets that occupy the same location (regardless of winding).
  std::map<std::vector<long long>, std::vector<size_t>> facetGroups;
  for (size_t i=0; i<numFacets; ++i)
    facetGroups[canonical_quantized_facet_key(facets[i])].push_back(i);

  // Every facet in a coincident group shares the same vertices, so each normal is either aligned
  // with or opposed to the group's reference normal.  Opposing orientations cancel pairwise (both
  // facets of each pair are degenerate for the average-normal distance); the surviving same-normal
  // duplicates then collapse to a single representative.  Net effect per group: keep one facet when
  // the two orientations are unbalanced, or none when they are balanced (everything cancels).
  std::vector<bool> shouldRemove(numFacets, false);
  size_t numToRemove = 0;
  for (auto && entry : facetGroups)
  {
    const std::vector<size_t> & members = entry.second;
    if (members.size() < 2)
      continue;

    const stk::math::Vector3d refNormal = facets[members[0]].facet_normal();
    std::vector<size_t> sameSide;      // normal aligned with refNormal (always includes members[0])
    std::vector<size_t> oppositeSide;  // normal opposes refNormal
    for (size_t idx : members)
    {
      if (Dot(facets[idx].facet_normal(), refNormal) < 0.)
        oppositeSide.push_back(idx);
      else
        sameSide.push_back(idx);
    }

    // Mark the whole group for removal, then keep a single majority-side representative unless the
    // orientations are balanced (in which case the group cancels entirely).
    for (size_t idx : members)
      shouldRemove[idx] = true;
    numToRemove += members.size();

    if (sameSide.size() != oppositeSide.size())
    {
      const size_t keepIdx = (sameSide.size() > oppositeSide.size()) ? sameSide[0] : oppositeSide[0];
      shouldRemove[keepIdx] = false;
      --numToRemove;
    }
  }

  if (numToRemove == 0)
    return 0;

  std::vector<FACET> keptFacets;
  keptFacets.reserve(numFacets - numToRemove);
  for (size_t i=0; i<numFacets; ++i)
    if (!shouldRemove[i])
      keptFacets.push_back(facets[i]);
  facets.swap(keptFacets);

  return numToRemove;
}

namespace {

std::array<long long,3> quantize_point(const stk::math::Vector3d & p, const double cell)
{
  std::array<long long,3> q;
  for (int d=0; d<3; ++d)
    q[d] = (cell > 0.) ? std::llround(p[d]/cell) : 0ll;
  return q;
}

// Enumerates the edges of a facet as (winding/order-independent quantized key, representative
// endpoint coordinates).  3D facets contribute their three triangle edges; 2D facets (segments)
// contribute their two endpoints, since a watertight closed curve requires each endpoint to be
// shared by exactly two segments.
template<class FACET>
std::vector<std::pair<std::vector<long long>, std::array<stk::math::Vector3d,2>>>
enumerate_facet_edges(const FACET & facet, const double cell)
{
  std::vector<std::pair<std::vector<long long>, std::array<stk::math::Vector3d,2>>> edges;
  if constexpr (3 == FACET::DIM)
  {
    static constexpr int edgeVerts[3][2] = {{0,1},{1,2},{2,0}};
    edges.reserve(3);
    for (auto && ev : edgeVerts)
    {
      const stk::math::Vector3d & p0 = facet.facet_vertex(ev[0]);
      const stk::math::Vector3d & p1 = facet.facet_vertex(ev[1]);
      std::array<long long,3> q0 = quantize_point(p0, cell);
      std::array<long long,3> q1 = quantize_point(p1, cell);
      if (q1 < q0) // unoriented: sort the two quantized endpoints so both windings match
        std::swap(q0, q1);
      edges.push_back({{q0[0],q0[1],q0[2],q1[0],q1[1],q1[2]}, {p0, p1}});
    }
  }
  else
  {
    edges.reserve(2);
    for (int v=0; v<2; ++v)
    {
      const stk::math::Vector3d & p = facet.facet_vertex(v);
      const std::array<long long,3> q = quantize_point(p, cell);
      edges.push_back({{q[0],q[1],q[2]}, {p, p}});
    }
  }
  return edges;
}

// Deterministic owner processor for an edge key (pure function of key contents, so every processor
// that produces this key agrees on where its tally is accumulated).  FNV-1a hash.
size_t hash_edge_key(const std::vector<long long> & key)
{
  size_t h = 1469598103934665603ull;
  for (long long k : key)
  {
    h ^= static_cast<size_t>(k);
    h *= 1099511628211ull;
  }
  return h;
}

}

template<class FACET>
double compute_global_edge_quantization_cell_size(const stk::ParallelMachine comm, const std::vector<FACET> & localFacets)
{
  constexpr double relTol = 1.e-6;
  double globalMinMeanSqEdge = std::numeric_limits<double>::max();
  for (auto && facet : localFacets)
    globalMinMeanSqEdge = std::min(globalMinMeanSqEdge, facet.mean_squared_edge_length());
  krino::all_reduce_min(comm, globalMinMeanSqEdge);
  if (globalMinMeanSqEdge == std::numeric_limits<double>::max() || globalMinMeanSqEdge <= 0.)
    return 0.;
  return relTol * std::sqrt(globalMinMeanSqEdge);
}

template<class FACET>
std::vector<std::array<stk::math::Vector3d,2>>
find_non_watertight_edges(const std::vector<FACET> & localFacets,
    const stk::ParallelMachine comm,
    const double globalEdgeCellSize,
    const std::optional<BoundingBox> & restrictBbox)
{
  struct EdgeData { int count = 0; std::array<stk::math::Vector3d,2> coords; };

  // Tally proc-local edge incidence counts keyed by the winding/order-independent quantized key.
  std::map<std::vector<long long>, EdgeData> localEdges;
  for (auto && facet : localFacets)
    for (auto && edge : enumerate_facet_edges(facet, globalEdgeCellSize))
    {
      EdgeData & data = localEdges[edge.first];
      if (data.count == 0)
        data.coords = edge.second;
      ++data.count;
    }

  auto edge_is_within_restriction = [&restrictBbox](const std::array<stk::math::Vector3d,2> & coords)
  {
    if (!restrictBbox.has_value())
      return true;
    return restrictBbox->contains(coords[0]) && restrictBbox->contains(coords[1]);
  };

  // An edge lying along a face of the restriction box is a boundary produced by intersecting the
  // surface with the meshed domain, not a crack in the surface, so it must not be flagged.  An edge
  // is "on the box" when both endpoints share the box's min or max plane in some dimension.  The box
  // is stored in single precision, so the on-plane tolerance admits the float round-off of the plane
  // coordinate (relative to its magnitude), floored by the edge quantization cell size.  Dimensions
  // in which the box has no extent (e.g. the z direction of a 2D domain) are not faces and so are
  // skipped -- otherwise every edge would trivially lie on that degenerate plane.
  auto edge_is_on_bounding_box = [&restrictBbox, globalEdgeCellSize](const std::array<stk::math::Vector3d,2> & coords)
  {
    if (!restrictBbox.has_value())
      return false;
    const BoundingBox & box = *restrictBbox;
    auto plane_tol = [globalEdgeCellSize](const double plane)
    {
      return std::max(globalEdgeCellSize, 8. * std::numeric_limits<float>::epsilon() * std::abs(plane));
    };
    auto both_on_plane = [&coords, &plane_tol](const int d, const double plane)
    {
      const double tol = plane_tol(plane);
      return std::abs(coords[0][d] - plane) <= tol && std::abs(coords[1][d] - plane) <= tol;
    };
    for (int d = 0; d < 3; ++d)
    {
      const double lo = static_cast<double>(box.get_min()[d]);
      const double hi = static_cast<double>(box.get_max()[d]);
      if (hi - lo <= plane_tol(hi) + plane_tol(lo))
        continue; // degenerate/flat dimension (e.g. the z of a 2D domain): not a face
      if (both_on_plane(d, lo) || both_on_plane(d, hi))
        return true;
    }
    return false;
  };

  auto edge_is_reportable = [&edge_is_within_restriction, &edge_is_on_bounding_box](const std::array<stk::math::Vector3d,2> & coords)
  {
    return edge_is_within_restriction(coords) && !edge_is_on_bounding_box(coords);
  };

  std::vector<std::array<stk::math::Vector3d,2>> badEdges;

  const int numProcs = stk::parallel_machine_size(comm);
  if (numProcs == 1)
  {
    for (auto && entry : localEdges)
      if (entry.second.count != 2 && edge_is_reportable(entry.second.coords))
        badEdges.push_back(entry.second.coords);
    return badEdges;
  }

  // Send each local edge tally to the deterministic owner of its key.  Two-step CommSparse
  // (size then communicate), mirroring parallel_distribute_facets in Akri_Faceted_Surface.cpp.
  stk::CommSparse commSparse(comm);
  for (int commStep = 0; commStep < 2; ++commStep)
  {
    for (auto && entry : localEdges)
    {
      const std::vector<long long> & key = entry.first;
      const int owner = static_cast<int>(hash_edge_key(key) % static_cast<size_t>(numProcs));
      stk::CommBuffer & b = commSparse.send_buffer(owner);
      const size_t keyLen = key.size();
      b.pack(keyLen);
      for (long long k : key)
        b.pack(k);
      b.pack(entry.second.count);
      for (int i=0; i<2; ++i)
        for (int d=0; d<3; ++d)
          b.pack(entry.second.coords[i][d]);
    }
    if (commStep == 0)
      commSparse.allocate_buffers();
    else
      commSparse.communicate();
  }

  // Owner sums per-key contributions from all processors and decides watertightness.
  std::map<std::vector<long long>, EdgeData> ownedEdges;
  for (int p = 0; p < numProcs; ++p)
  {
    stk::CommBuffer & b = commSparse.recv_buffer(p);
    while (b.remaining())
    {
      size_t keyLen = 0;
      b.unpack(keyLen);
      std::vector<long long> key(keyLen);
      for (size_t i=0; i<keyLen; ++i)
        b.unpack(key[i]);
      int count = 0;
      b.unpack(count);
      std::array<stk::math::Vector3d,2> coords;
      for (int i=0; i<2; ++i)
        for (int d=0; d<3; ++d)
          b.unpack(coords[i][d]);
      EdgeData & data = ownedEdges[key];
      if (data.count == 0)
        data.coords = coords;
      data.count += count;
    }
  }

  for (auto && entry : ownedEdges)
    if (entry.second.count != 2 && edge_is_reportable(entry.second.coords))
      badEdges.push_back(entry.second.coords);

  return badEdges;
}

// Explicit template instantiation.
//
// Every free function in this file is instantiated for each concrete facet type.
// The list is generated by the macro below so that adding a calc function is a
// single new line inside the macro (automatically covering all facet types) and
// adding a facet type is a single new macro invocation.  This keeps the four
// per-type blocks in lock-step -- there is no way to instantiate a function for
// some types but forget others.  The generated code is identical to writing each
// instantiation by hand, so compute/storage characteristics are unchanged.

#define KRINO_INSTANTIATE_FACETED_SURFACE_CALCS(FACET) \
  template stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets); \
  template stk::math::Vector3d compute_closest_point_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets); \
  template double point_distance_given_nearest_facets<FACET>(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance); \
  template double compute_point_to_facets_distance_by_average_normal<FACET>(const stk::math::Vector3d &x, const std::vector<const FACET*> & facets); \
  template stk::math::Vector3d compute_pseudo_normal<FACET>(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets); \
  template std::pair<int, double> compute_facet_edge_intersection<FACET>(const FACET & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1); \
  template double compute_intersection_between_surface_facets_and_edge<FACET>(const std::vector<const FACET*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1); \
  template size_t remove_coincident_facets<FACET>(std::vector<FACET> & facets); \
  template double compute_global_edge_quantization_cell_size<FACET>(const stk::ParallelMachine comm, const std::vector<FACET> & localFacets); \
  template std::vector<std::array<stk::math::Vector3d,2>> find_non_watertight_edges<FACET>(const std::vector<FACET> & localFacets, const stk::ParallelMachine comm, const double globalEdgeCellSize, const std::optional<BoundingBox> & restrictBbox);

KRINO_INSTANTIATE_FACETED_SURFACE_CALCS(Facet2d)
KRINO_INSTANTIATE_FACETED_SURFACE_CALCS(Facet3d)
KRINO_INSTANTIATE_FACETED_SURFACE_CALCS(FacetWithVelocity2d)
KRINO_INSTANTIATE_FACETED_SURFACE_CALCS(FacetWithVelocity3d)

#undef KRINO_INSTANTIATE_FACETED_SURFACE_CALCS

}
