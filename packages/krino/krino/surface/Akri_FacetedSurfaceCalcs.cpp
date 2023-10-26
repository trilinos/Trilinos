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
    double min_sqr_dist = std::numeric_limits<double>::max();
    for ( auto&& facet : nearestFacets )
    {
      const double sqr_dist = facet->point_distance_squared(x);
      if (sqr_dist < min_sqr_dist)
      {
        min_sqr_dist = sqr_dist;
      }
    }
    if (0.0 != narrow_band_size && min_sqr_dist > narrow_band_size*narrow_band_size)
    {
      dist = far_field_value;
    }
    else
    {
      dist = std::sqrt(min_sqr_dist);
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

  if (!closest_point_on_edge) return facet_queries[nearest].signed_distance(x);

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

}


