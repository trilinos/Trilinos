/*
 * Akri_FacetedSurfaceCalcs.cpp
 *
 *  Created on: May 22, 2023
 *      Author: drnoble
 */
#include <Akri_BoundingBoxDistance.hpp>
#include "Akri_FacetedSurfaceCalcs.hpp"

#include <stk_util/util/ReportHandler.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Sign.hpp>

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

template<class FACET>
stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets)
{
  double minSqrDist = std::numeric_limits<double>::max();
  stk::math::Vector3d closestPt;
  stk::math::Vector3d facetClosestPt;
  for ( auto&& facet : nearestFacets )
  {
    facet->closest_point(x, facetClosestPt);
    const double sqrDist = (x-facetClosestPt).length_squared();
    if (sqrDist < minSqrDist)
    {
      minSqrDist = sqrDist;
      closestPt = facetClosestPt;
    }
  }
  return closestPt;
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
unsigned find_index_of_closest_facet(const std::vector<FacetDistanceQuery<FACET>> & facetDistQueries)
{
  unsigned closest = 0;
  for ( unsigned index=0; index<facetDistQueries.size(); ++index )
    if ( facetDistQueries[index].distance_squared() < facetDistQueries[closest].distance_squared() )
      closest = index;
  return closest;
}

template<class FACET>
stk::math::Vector3d compute_pseudo_normal(const std::vector<FacetDistanceQuery<FACET>> & facet_queries, const unsigned nearest)
{
  const double tol = 1.e-6;
  stk::math::Vector3d pseudo_normal = stk::math::Vector3d::ZERO;
  stk::math::Vector3d average_normal = stk::math::Vector3d::ZERO;

  const stk::math::Vector3d nearest_closest_point = facet_queries[nearest].closest_point();
  const double nearest_size2 = facet_queries[nearest].facet().mean_squared_edge_length();
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

template<class FACET>
double compute_point_to_facets_distance_by_average_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & facets)
{

  // If the closest_point weights are all larger than this value, then the closest point
  // is considered to be on the face of the closest facet rather than on the edges of the facet, and
  // therefore only the closest facet is considered in the distance calculation.  Otherwise, all of the
  // facets are considered to compute an average normal in order to compute the distance.
  const double edge_tol = 1.e-6;

  const std::vector<FacetDistanceQuery<FACET>> facetDistQueries = build_distance_queries(x, facets);
  STK_ThrowRequireMsg(!facetDistQueries.empty(), "All facets are degenerate in compute_point_to_facets_distance_by_average_normal.");

  const unsigned nearest = find_index_of_closest_facet(facetDistQueries);

  if ( facetDistQueries[nearest].distance_squared() == 0. )
  {
    return 0.0;
  }

  const stk::math::Vector3d closest_pt_wts = facetDistQueries[nearest].closest_point_weights();
  const int dim = FACET::DIM;

  bool closest_point_on_edge = false;
  for (int d=0; d<dim; ++d)
  {
    if (closest_pt_wts[d] < edge_tol)
    {
      closest_point_on_edge = true;
      break;
    }
  }

  if (!closest_point_on_edge)
  {
    return facetDistQueries[nearest].signed_distance(x);
  }

  const double min_sqr_dist = facetDistQueries[nearest].distance_squared();
  const stk::math::Vector3d pseudo_normal = compute_pseudo_normal(facetDistQueries, nearest);

  if (pseudo_normal.length_squared() == 0.0)
  {
    krinolog << "Warning:  Cannot determine the average facet normal for computing the level set distance at point " << x
        << ".  This can happen when faceted facet include coincident facets with opposite normals.  Arbitrarily setting the distance to be positive." << stk::diag::dendl;
    return std::sqrt(min_sqr_dist);
  }
  else
  {
    if (Dot(pseudo_normal, x-facetDistQueries[nearest].closest_point()) > 0)
    {
      return std::sqrt(min_sqr_dist);
    }
    else
    {
      return -std::sqrt(min_sqr_dist);
    }
  }
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

  const unsigned nearest = find_index_of_closest_facet(facetDistQueries);

  return compute_pseudo_normal(facetDistQueries, nearest);
}

// Explicit template instantiation

template stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const Facet2d*> & nearestFacets);
template stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const Facet3d*> & nearestFacets);
template double point_distance_given_nearest_facets<Facet2d>(const stk::math::Vector3d &x, const std::vector<const Facet2d*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);
template double point_distance_given_nearest_facets<Facet3d>(const stk::math::Vector3d &x, const std::vector<const Facet3d*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);
template double compute_point_to_facets_distance_by_average_normal<Facet2d>(const stk::math::Vector3d &x, const std::vector<const Facet2d*> & facets);
template double compute_point_to_facets_distance_by_average_normal<Facet3d>(const stk::math::Vector3d &x, const std::vector<const Facet3d*> & facets);
template stk::math::Vector3d compute_pseudo_normal<Facet2d>(const stk::math::Vector3d &x, const std::vector<const Facet2d*> & nearestFacets);
template stk::math::Vector3d compute_pseudo_normal<Facet3d>(const stk::math::Vector3d &x, const std::vector<const Facet3d*> & nearestFacets);
template std::pair<int, double> compute_facet_edge_intersection<Facet2d>(const Facet2d & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1);
template std::pair<int, double> compute_facet_edge_intersection<Facet3d>(const Facet3d & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1);
template double compute_intersection_between_surface_facets_and_edge<Facet2d>(const std::vector<const Facet2d*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1);
template double compute_intersection_between_surface_facets_and_edge<Facet3d>(const std::vector<const Facet3d*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1);

template stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity2d*> & nearestFacets);
template double point_distance_given_nearest_facets<FacetWithVelocity2d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity2d*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);
template double compute_point_to_facets_distance_by_average_normal<FacetWithVelocity2d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity2d*> & facets);
template stk::math::Vector3d compute_pseudo_normal<FacetWithVelocity2d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity2d*> & nearestFacets);
template std::pair<int, double> compute_facet_edge_intersection<FacetWithVelocity2d>(const FacetWithVelocity2d & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1);
template double compute_intersection_between_surface_facets_and_edge<FacetWithVelocity2d>(const std::vector<const FacetWithVelocity2d*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1);

template stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity3d*> & nearestFacets);
template double point_distance_given_nearest_facets<FacetWithVelocity3d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity3d*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);
template double compute_point_to_facets_distance_by_average_normal<FacetWithVelocity3d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity3d*> & facets);
template stk::math::Vector3d compute_pseudo_normal<FacetWithVelocity3d>(const stk::math::Vector3d &x, const std::vector<const FacetWithVelocity3d*> & nearestFacets);
template std::pair<int, double> compute_facet_edge_intersection<FacetWithVelocity3d>(const FacetWithVelocity3d & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1);
template double compute_intersection_between_surface_facets_and_edge<FacetWithVelocity3d>(const std::vector<const FacetWithVelocity3d*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1);

}
