/*
 * Akri_FacetedSurfaceCalcs.hpp
 *
 *  Created on: May 22, 2023
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_
#include <stk_math/StkVector.hpp>
#include <Akri_Facet.hpp>

namespace krino {
  std::vector<BoundingBox> fill_processor_bounding_boxes(const BoundingBox & localFacetBbox, const BoundingBox & localQueryBbox, const double truncationLength);
  double point_distance_given_nearest_facets(const stk::math::Vector3d &x, const FacetVec & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);
  double compute_point_to_facets_distance_by_average_normal(const stk::math::Vector3d &x, const FacetVec & facets);
  stk::math::Vector3d compute_pseudo_normal(const unsigned dim, const std::vector<FacetDistanceQuery> & facet_queries, const unsigned nearest);
}



#endif /* KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_ */
