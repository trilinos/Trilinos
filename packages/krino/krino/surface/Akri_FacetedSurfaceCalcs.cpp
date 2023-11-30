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

static double compute_point_distance_squared(const stk::math::Vector3d &x, const FacetVec & nearestFacets)
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

double
point_distance_given_nearest_facets(const stk::math::Vector3d &x, const FacetVec & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance)
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

double
compute_point_to_facets_distance_by_average_normal(const stk::math::Vector3d &x, const FacetVec & facets)
{

  // If the closest_point weights are all larger than this value, then the closest point
  // is considered to be on the face of the closest facet rather than on the edges of the facet, and
  // therefore only the closest facet is considered in the distance calculation.  Otherwise, all of the
  // facets are considered to compute an average normal in order to compute the distance.
  const double edge_tol = 1.e-6;

  std::vector<FacetDistanceQuery> facet_queries;
  facet_queries.reserve(facets.size());
  for ( auto&& facet : facets )
  {
    if (facet->degenerate()) continue; // Skip zero-sized facets
    facet_queries.emplace_back(*facet, x);
  }

  STK_ThrowRequireMsg(!facet_queries.empty(), "All facets are degenerate in compute_point_to_facets_distance_by_average_normal.");

  unsigned nearest = 0;
  for ( unsigned index=0; index<facet_queries.size(); ++index )
  {
    if ( facet_queries[index].distance_squared() < facet_queries[nearest].distance_squared() )
    {
      nearest = index;
    }
  }

  if ( facet_queries[nearest].distance_squared() == 0. )
  {
    return 0.0;
  }

  const stk::math::Vector3d closest_pt_wts = facet_queries[nearest].closest_point_weights();
  const int dim = dynamic_cast<krino::Facet3d *>(facets[nearest]) ? 3 : 2;

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
    return facet_queries[nearest].signed_distance(x);
  }

  const double min_sqr_dist = facet_queries[nearest].distance_squared();
  const stk::math::Vector3d pseudo_normal = compute_pseudo_normal(dim, facet_queries, nearest);

  if (pseudo_normal.length_squared() == 0.0)
  {
    krinolog << "Warning:  Cannot determine the average facet normal for computing the level set distance at point " << x
        << ".  This can happen when faceted facet include coincident facets with opposite normals.  Arbitrarily setting the distance to be positive." << stk::diag::dendl;
    return std::sqrt(min_sqr_dist);
  }
  else
  {
    if (Dot(pseudo_normal, x-facet_queries[nearest].closest_point()) > 0)
    {
      return std::sqrt(min_sqr_dist);
    }
    else
    {
      return -std::sqrt(min_sqr_dist);
    }
  }
}

stk::math::Vector3d
compute_pseudo_normal(const unsigned dim, const std::vector<FacetDistanceQuery> & facet_queries, const unsigned nearest)
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
      const Facet & facet = query.facet();

      average_normal += facet.facet_normal();

      if (3 == dim)
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

  return (3 == dim && close_count > 2) ? pseudo_normal : average_normal;
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

bool is_projection_of_point_inside_enlarged_facet(const Facet & facet, const stk::math::Vector3d& p)
{
  const int dim = (dynamic_cast<const krino::Facet3d *>(&facet)) ? 3 : 2;

  if (3 == dim)
    return is_projection_of_point_inside_enlarged_triangle(facet.facet_vertex(0), facet.facet_vertex(1), facet.facet_vertex(2), p);
  return is_projection_of_point_inside_enlarged_segment(facet.facet_vertex(0), facet.facet_vertex(1), p);
}

std::pair<int, double> compute_facet_edge_intersection(const Facet & facet,
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

std::pair<int, double> compute_intersection_between_and_surface_facets_and_edge(const std::vector<Facet*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1)
{
  if (candidates.empty())
    return {0, -1.};

  const double dist0 = compute_point_to_facets_distance_by_average_normal(edgePt0, candidates);
  const double dist1 = compute_point_to_facets_distance_by_average_normal(edgePt1, candidates);

  if (!sign_change(dist0, dist1))
    return {0, -1.};

  if (0. == dist0)
    return {-1, 0.};
  if (0. == dist1)
    return {1, 1.};

  bool haveCrossing = false;
  double intersectionLoc = -1.;

  for (const Facet * surfFacet : candidates)
  {
    const auto [facetCrossingSign, facetIntersectionLoc] = compute_facet_edge_intersection(*surfFacet, edgePt0, edgePt1);
    if (facetCrossingSign != 0)
    {
      if (!haveCrossing || std::abs(facetIntersectionLoc-0.5) < std::abs(intersectionLoc-0.5)) // pick intersection closest to middle of edge
        intersectionLoc = facetIntersectionLoc;
      haveCrossing = true;
    }
  }

  if (haveCrossing)
    return {sign(dist1), intersectionLoc};

  // Sign change, but no crossing. This could be because there is a small gap in the facets, or there could be
  // an unterminated surface far away.  Decide based on magnitude of the distance at the linear crossing location.

  const double linearCrossingLoc = dist0 / (dist0 - dist1);

  const stk::math::Vector3d linearCrossingPt = (1.-linearCrossingLoc) * edgePt0 + linearCrossingLoc * edgePt1;
  const double minSqrDist = compute_point_distance_squared(linearCrossingPt, candidates);
  const double edgeSqrLen = (edgePt1 - edgePt0).length_squared();

  constexpr double sqrTol = 1.e-4; // Large tol here is fine, right?
  if (minSqrDist < sqrTol*edgeSqrLen)
    return {sign(dist1), linearCrossingLoc};

  return {0, -1.};
}

}
