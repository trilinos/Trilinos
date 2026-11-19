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
#include <Akri_Facet.hpp>

#include <array>
#include <functional>
#include <optional>
#include <vector>

namespace krino {
// The processor a facet read from a file is sent to: the first processor p >= 1 whose bounding box
// it intersects, else 0.
template <class FACET>
int find_destination_proc_for_facet(
    const std::vector<BoundingBox> & procBboxes, const FACET & facet)
{
  const int numProcs = procBboxes.size();
  for (int destProc = 1; destProc < numProcs; ++destProc)
    if (facet.does_intersect(procBboxes[destProc])) return destProc;
  return 0;
}

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

  // Distance queries for the non-degenerate facets, and the closest of them (ties broken by
  // centroid, so the result does not depend on the facet order).
  template <class FACET>
  std::vector<FacetDistanceQuery<FACET>>
  build_distance_queries(const stk::math::Vector3d & x, const std::vector<const FACET *> & facets);
  template <class FACET>
  const FacetDistanceQuery<FACET> &
  find_closest_facet(const std::vector<FacetDistanceQuery<FACET>> & facetDistQueries);

  // The pseudo normal at the closest point, from the facets whose closest points are within a small
  // tolerance of it (whether or not they are connected).
  template <class FACET>
  stk::math::Vector3d
  compute_pseudo_normal(const std::vector<FacetDistanceQuery<FACET>> & facetDistQueries,
      const FacetDistanceQuery<FACET> & closestQuery);

  // Whether a closest point with these facet weights is on an edge or vertex of the facet rather
  // than inside it.
  bool is_closest_point_on_facet_edge_or_vertex(
      const stk::math::Vector3d & closestFacetWeights, const int dim);

  template <class FACET>
  double signed_distance_given_pseudo_normal(const stk::math::Vector3d & x,
      const FacetDistanceQuery<FACET> & closestQuery,
      const stk::math::Vector3d & pseudoNormal);
  template <class FACET>
  stk::math::Vector3d closest_point_normal_given_pseudo_normal(const stk::math::Vector3d & x,
      const FacetDistanceQuery<FACET> & closestQuery,
      const stk::math::Vector3d & pseudoNormal);

  // The signed distance, and the closest point normal, given the closest facet.  When the closest
  // point is on an edge or vertex of that facet, the sign comes from the pseudo normal there, given
  // by compute_pseudo_normal (a callable returning it): the disconnected faceted surfaces estimate
  // it from the nearby facets, a connected one from the facets that share the edge or vertex.
  template <class FACET, class COMPUTE_PSEUDO_NORMAL>
  double signed_distance_given_closest_facet(const stk::math::Vector3d & x,
      const FacetDistanceQuery<FACET> & closestQuery,
      const COMPUTE_PSEUDO_NORMAL & compute_pseudo_normal)
  {
    constexpr double sqrTolSoCloseThatItIsOnClosestFacet = 1.e-18;
    if (closestQuery.distance_squared() / closestQuery.facet().mean_squared_edge_length() <
        sqrTolSoCloseThatItIsOnClosestFacet)
      return 0.0;
    if (!is_closest_point_on_facet_edge_or_vertex(closestQuery.closest_point_weights(), FACET::DIM))
      return closestQuery.signed_distance(x);
    return signed_distance_given_pseudo_normal(x, closestQuery, compute_pseudo_normal());
  }

  template <class FACET, class COMPUTE_PSEUDO_NORMAL>
  stk::math::Vector3d closest_point_normal_given_closest_facet(const stk::math::Vector3d & x,
      const FacetDistanceQuery<FACET> & closestQuery,
      const COMPUTE_PSEUDO_NORMAL & compute_pseudo_normal)
  {
    if (!is_closest_point_on_facet_edge_or_vertex(closestQuery.closest_point_weights(), FACET::DIM))
      return closestQuery.facet().facet_normal();
    return closest_point_normal_given_pseudo_normal(x, closestQuery, compute_pseudo_normal());
  }

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
  // The facets remove_coincident_facets removes.  Which facet of a group is kept depends only on
  // the coordinates of the group's vertices, not on the order of the facets.
  template <class FACET>
  std::vector<bool> find_coincident_facets_to_remove(const std::vector<FACET> & facets);
  // The location of a facet, the same for the facets coincident with it whatever their winding:
  // its vertices quantized to a grid scaled by its size, sorted, and padded with zeros to 3.
  template <class FACET> std::array<long long, 9> coincident_facet_key(const FACET & facet);
  // For a group of facets with the same coincident_facet_key, whether
  // find_coincident_facets_to_remove removes each.
  template <class FACET>
  std::vector<bool>
  find_coincident_group_facets_to_remove(const std::vector<const FACET *> & group);

  // Whether a non-watertight edge should be reported: with restrictBbox, only edges with both
  // endpoints inside it that do not lie on one of its (non-degenerate) faces, the open boundary
  // from clipping the surface to the meshed domain.  onBoxTol is the least tolerance for lying on a
  // face.
  bool is_reportable_non_watertight_edge(const std::array<stk::math::Vector3d, 2> & coords,
      const std::optional<BoundingBox> & restrictBbox,
      const double onBoxTol);

  // Collective.  Warns (once, on processor 0) if any processor found non-watertight edges, with a
  // sample of them.
  void report_non_watertight_edges(const stk::ParallelMachine comm,
      const std::vector<std::array<stk::math::Vector3d, 2>> & localBadEdges);

  // Computes a single global grid cell size for winding/order-independent edge matching in the
  // watertight check.  A shared edge between two differently-sized facets must quantize to the
  // same key, so the cell size is derived from a global (all-reduced) measure rather than each
  // facet's own size.  Returns relTol * sqrt(global min mean_squared_edge_length).  Collective.
  template<class FACET>
  double compute_global_edge_quantization_cell_size(const stk::ParallelMachine comm, const std::vector<FACET> & localFacets);
  // The same, with the numLocalFacets local facets given one at a time by local_facet, so that they
  // need not all be stored (as when made from a FacetMesh).
  template <class FACET>
  double compute_global_edge_quantization_cell_size(const stk::ParallelMachine comm,
      const size_t numLocalFacets,
      const std::function<FACET(size_t)> & local_facet);

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
  // The same, with the local facets given one at a time, as for
  // compute_global_edge_quantization_cell_size.
  template <class FACET>
  std::vector<std::array<stk::math::Vector3d, 2>>
  find_non_watertight_edges(const size_t numLocalFacets,
      const std::function<FACET(size_t)> & local_facet,
      const stk::ParallelMachine comm,
      const double globalEdgeCellSize,
      const std::optional<BoundingBox> & restrictBbox);
}



#endif /* KRINO_KRINO_SURFACE_AKRI_FACETEDSURFACECALCS_HPP_ */
