/*
 * Akri_FacetedSurfaceCalcs.hpp
 *
 *  Created on: May 22, 2023
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_BoundingBox.hpp>

#include <array>
#include <optional>
#include <vector>

namespace krino {
  std::vector<BoundingBox> fill_processor_bounding_boxes(const BoundingBox & localFacetBbox, const BoundingBox & localQueryBbox, const double truncationLength);

  template<class FACET>
  stk::math::Vector3d compute_closest_point(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets);

  template<class FACET>
  stk::math::Vector3d compute_closest_point_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets);

  template<class FACET>
  double point_distance_given_nearest_facets(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance);

  template<class FACET>
  double compute_point_to_facets_distance_by_average_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & facets);

  template<class FACET>
  stk::math::Vector3d compute_pseudo_normal(const stk::math::Vector3d &x, const std::vector<const FACET*> & nearestFacets);

  template<class FACET>
  std::pair<int, double> compute_facet_edge_intersection(const FACET & facet, const stk::math::Vector3d& edgePt0, const stk::math::Vector3d& edgePt1);

  template<class FACET>
  double compute_intersection_between_surface_facets_and_edge(const std::vector<const FACET*> & candidates, const stk::math::Vector3d & edgePt0, const stk::math::Vector3d & edgePt1);

  // Removes redundant coincident facets that share the same vertices (to within a relative
  // tolerance scaled by facet size).  Within each coincident group: opposing-normal orientations
  // cancel pairwise (both facets of each pair are erased, since they otherwise produce an
  // indeterminate zero-length pseudo-normal), and the surviving same-normal duplicates collapse to
  // a single representative (extra copies double-count the facet normal in the average-normal
  // calculation).  A facet with no coincident partner is left untouched.  Net effect per group:
  // one facet remains when the orientations are unbalanced, none when they are balanced.  Returns
  // the number of facets removed.
  template<class FACET>
  size_t remove_coincident_facets(std::vector<FACET> & facets);

  // Computes a single global grid cell size for winding/order-independent edge matching in the
  // watertight check.  A shared edge between two differently-sized facets must quantize to the
  // same key, so the cell size is derived from a global (all-reduced) measure rather than each
  // facet's own size.  Returns relTol * sqrt(global min mean_squared_edge_length).  Collective.
  template<class FACET>
  double compute_global_edge_quantization_cell_size(const stk::ParallelMachine comm, const std::vector<FACET> & localFacets);

  // Detects non-watertight edges of a faceted surface: in 3D, triangle edges incident to a number
  // of facets other than 2; in 2D, segment endpoints incident to a number of segments other than 2.
  // A watertight surface has every edge/endpoint shared by exactly two facets; count 1 indicates a
  // boundary/crack and count > 2 a non-manifold junction.  Edges are matched across processors via a
  // deterministic owner assignment, so this is collective over comm and correct regardless of how
  // facets are distributed.  When restrictBbox is set, only edges with BOTH endpoints inside it are
  // reported, and edges lying on a (non-degenerate) face of the box are excluded: these are the open
  // boundary produced by clipping the surface to the meshed domain, not cracks in the surface.
  // Returns, on every processor, its owner-local list of offending edges (each as its two endpoint
  // coordinates); the union across processors is the full set (the caller reduces the count).
  template<class FACET>
  std::vector<std::array<stk::math::Vector3d,2>>
  find_non_watertight_edges(const std::vector<FACET> & localFacets,
      const stk::ParallelMachine comm,
      const double globalEdgeCellSize,
      const std::optional<BoundingBox> & restrictBbox);
}



#endif /* KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_ */
