// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_MESHFACETEDSURFACE_HPP_
#define KRINO_KRINO_SURFACE_AKRI_MESHFACETEDSURFACE_HPP_

#include <Akri_Facet.hpp>
#include <Akri_FacetCoordinates.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetMeshChecks.hpp>
#include <Akri_SearchTree.hpp>
#include <Akri_Surface.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace krino
{

// A faceted surface (triangles in 3D, segments in 2D) whose facets either share vertices
// (connected) or each have their own copies of their vertices (an unconnected soup).  Each
// processor owns the facets build_local_facets gives it at the first prepare_to_compute (by default
// the ones it was constructed with), or at every prepare_to_compute for a subclass that rebuilds
// them.  Every prepare_to_compute then copies to each processor the remote facets within the narrow
// band of its query points, with their vertex ids, so that a vertex id names one local vertex
// shared by all the local facets that use it.
class MeshFacetedSurface
    : public SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign
{
public:
  // Collective.  The facets this processor owns, each given by the global ids of its vertices (3
  // per triangle in 3D, the first 2 of each entry per segment in 2D), and the coordinates of every
  // vertex those facets use.  A vertex id names the same vertex on every processor that uses it.
  // haveConnectedSurface says that facets meeting at a point share its vertex: the sign at a
  // closest point on an edge or vertex then comes from the facets sharing it.  Otherwise (for
  // example when coincident vertices have different ids) it comes from the facets near the closest
  // point.
  MeshFacetedSurface(const stk::ParallelMachine comm,
      const unsigned dim,
      const std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> & vertexIdsAndCoords,
      const std::vector<std::array<FacetVertexId, 3>> & facetVertexIds,
      const bool haveConnectedSurface = true);
  // As above, with vertex ids being indices (from 0) into vertices, which must then be the same on
  // every processor.
  MeshFacetedSurface(const stk::ParallelMachine comm,
      const unsigned dim,
      const std::vector<stk::math::Vector3d> & vertices,
      const std::vector<std::array<unsigned, 3>> & facetVertexIndices,
      const bool haveConnectedSurface = true);
  // As the first constructor, with the velocity of each vertex: vertexVelocities[i] is that of
  // vertexIdsAndCoords[i].
  MeshFacetedSurface(const stk::ParallelMachine comm,
      const unsigned dim,
      const std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> & vertexIdsAndCoords,
      const std::vector<std::array<FacetVertexId, 3>> & facetVertexIds,
      const std::vector<stk::math::Vector3d> & vertexVelocities,
      const bool haveConnectedSurface = true);
  // Collective.  The facets this processor owns, each given by the coordinates of its dim vertices
  // in order, with no shared vertices: the sign comes from the facets near the closest point, and
  // the facets are stored without vertex ids.
  MeshFacetedSurface(const stk::ParallelMachine comm,
      const unsigned dim,
      const std::vector<stk::math::Vector3d> & facetVertexCoords);
  virtual ~MeshFacetedSurface();

  virtual Surface_Type type() const override
  {
    return (myMesh.spatial_dimension() == 3) ? FACETED_SURFACE_3D : FACETED_SURFACE_2D;
  }
  virtual size_t storage_size() const override;
  virtual void prepare_to_compute(
      const double time, const BoundingBox & pointBbox, const double truncationLength) override;
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d & x,
      const double truncationLength,
      const double farFieldValue) const override;
  virtual std::pair<int, double> compute_intersection_with_segment(const stk::math::Vector3d & pt0,
      const stk::math::Vector3d & pt1,
      const double edgeCrossingTol) const override;
  // Every crossing of the segment with the facets near it, as where a segment crosses the surface
  // twice beside a sharp edge.
  virtual void fill_segment_crossings(const stk::math::Vector3d & pt0,
      const stk::math::Vector3d & pt1,
      const double edgeCrossingTol,
      std::vector<double> & crossingPositions) const override;
  virtual bool can_approximate_closest_point_normal() const override { return true; }
  virtual stk::math::Vector3d closest_point_normal(const stk::math::Vector3d & x) const override;
  // The facets searched by queries (owned and copied here by the last prepare_to_compute) whose
  // bounding boxes meet the region: the facets krino uses there.
  virtual void add_to_journal(CubitJournal & journal, const BoundingBox & region) const override;

  stk::math::Vector3d closest_point(const stk::math::Vector3d & x) const;
  // The velocity at the closest point on the surface, interpolated from its facet's vertices; zero
  // if there are no facets nearby.  Requires vertex velocities.
  stk::math::Vector3d velocity_at_closest_point(const stk::math::Vector3d & x) const;
  double point_unsigned_distance(const stk::math::Vector3d & x,
      const double truncationLength,
      const double farFieldValue) const;

  unsigned spatial_dimension() const { return myMesh.spatial_dimension(); }

  // Replace all the facets with these owned ones now, as the constructors give them: with no shared
  // vertices, or sharing vertices by id.  The vertex velocities, if any, are dropped, and so are
  // the ghosts until the next prepare_to_compute, which also removes coincident facets and checks
  // watertightness if asked.  Not for subclasses that build their own facets.
  void set_facets(const FacetCoordinates & facets);
  void set_facets(const FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords,
      const FacetMesh::FacetVertexIds & facetVertexIds,
      const bool haveConnectedSurface = true);
  // The coordinates of this processor's owned facets.
  FacetCoordinates owned_facet_coordinates() const;
  // Collective.  The global and local numbers of owned facets and their global areas (lengths in
  // 2D).
  std::string print_sizes() const;
  // The queries below are valid after the first prepare_to_compute, which builds the surface, or
  // after set_facets.
  const FacetMesh & facet_mesh() const { return myMesh; }
  size_t num_owned_facets() const { return myMesh.num_owned_facets(); }
  // Collective.  The number of distinct vertices of the facets over all processors.
  size_t global_num_vertices() const { return myMesh.global_num_vertices(); }
  // Facets copied to this processor by the last prepare_to_compute.
  size_t num_ghosted_facets() const { return myMesh.num_ghosted_facets(); }
  // Facets searched by queries after the last prepare_to_compute (owned and ghosted).
  size_t num_searched_facets() const { return myNumSearchedFacets; }

  // Applied to the owned facets each time they are built: coincident
  // facets are removed, then a warning is given for edges not shared by exactly two facets
  // (restricted to the box, if given).  A connected surface matches facets and edges by their
  // vertex ids, otherwise they are matched by location.
  void set_remove_coincident_facets(const bool flag) { myRemoveCoincidentFacets = flag; }
  void set_check_watertight(const bool flag) { myCheckWatertight = flag; }
  void set_watertight_bounding_box(const BoundingBox & bbox) { myWatertightBbox = bbox; }
  bool get_check_watertight() const { return myCheckWatertight; }
  const std::optional<BoundingBox> & get_watertight_bounding_box() const
  {
    return myWatertightBbox;
  }

  // Sharp features, off unless a feature angle (in degrees, as in Akri_SharpFeatureCalcs) is set:
  // the feature edges (3D) and corners of a connected surface, as find_sharp_features_by_vertex_ids
  // finds them each time the owned facets are built.  The surface then reports the feature edges
  // crossing a triangle and the corners in a tetrahedron (in a triangle in 2D), at which CDFEM cuts
  // the elements.  Open boundary edges are features only if asked for.
  void set_sharp_feature_angle(
      const double featureAngleInDegrees, const bool areOpenBoundaryEdgesFeatures = false);
  virtual bool has_edges() const override
  {
    return myCosSharpFeatureAngle.has_value() && spatial_dimension() == 3;
  }
  virtual bool has_corners() const override { return myCosSharpFeatureAngle.has_value(); }
  virtual void fill_triangle_intersection_parametric_coordinates(
      const std::array<stk::math::Vector3d, 3> & triNodes,
      std::vector<stk::math::Vector3d> & intParamCoords) const override;
  virtual void fill_tetrahedon_intersection_parametric_coordinates(
      const std::array<stk::math::Vector3d, 4> & tetNodes,
      std::vector<stk::math::Vector3d> & intParamCoords) const override;
  // The sharp features near this processor's query points, after prepare_to_compute.
  const FacetMeshSharpFeatures & local_sharp_features() const { return myLocalSharpFeatures; }

protected:
  using VertexIdsAndCoords = FacetMesh::VertexIdsAndCoords;
  using FacetVertexIds = FacetMesh::FacetVertexIds;

  stk::ParallelMachine parallel_comm() const { return myMesh.comm(); }

  // An empty surface whose facets come from build_local_facets, called by the first
  // prepare_to_compute, or by every one with rebuildOwnedFacetsEachPrepare (for facets that follow
  // a changing mesh).
  MeshFacetedSurface(const stk::ParallelMachine comm,
      const unsigned dim,
      const bool haveConnectedSurface,
      const FacetMeshOptions & meshOptions = FacetMeshOptions{},
      const bool rebuildOwnedFacetsEachPrepare = false);

  // Collective.  Called by prepare_to_compute (once, or every time if rebuilding), with every
  // processor's query box: the facets this processor will own and the coordinates of their
  // vertices.
  virtual void build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
      VertexIdsAndCoords & vertexIdsAndCoords,
      FacetVertexIds & facetVertexIds);

private:
  template <class FACET>
  std::vector<FACET> build_facets(const std::vector<unsigned> & facets) const;
  // The pseudo normal at the vertex or edge of closestFacet where the closest point (with these
  // weights) lies, from the facets that share it.
  template <class FACET>
  stk::math::Vector3d compute_connected_pseudo_normal(
      const unsigned closestFacet, const stk::math::Vector3d & closestPtWeights) const;
  template <class FACET>
  stk::math::Vector3d
  compute_pseudo_normal_at_closest_point(const std::vector<FacetDistanceQuery<FACET>> & queries,
      const FacetDistanceQuery<FACET> & closestQuery,
      const unsigned closestFacet) const;
  template <class FACET>
  double point_distance(const stk::math::Vector3d & x,
      const double truncationLength,
      const double farFieldValue,
      const bool computeSignedDistance) const;
  template <class FACET>
  stk::math::Vector3d
  closest_point_or_normal(const stk::math::Vector3d & x, const bool computeNormal) const;
  template <class FACET>
  stk::math::Vector3d velocity_at_closest_point_for_facets(const stk::math::Vector3d & x) const;
  template <class FACET>
  std::pair<int, double> compute_intersection_with_segment_for_facets(
      const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1) const;
  template <class FACET>
  void fill_segment_crossings_for_facets(const stk::math::Vector3d & pt0,
      const stk::math::Vector3d & pt1,
      const double edgeCrossingTol,
      std::vector<double> & crossingPositions) const;
  void build_search_tree(const BoundingBox & localPaddedQueryBbox);
  // Collective.  Gives each processor the sharp features within its box, and builds their search
  // trees.
  void distribute_sharp_features(const std::vector<BoundingBox> & procBboxes);
  void given_facets_were_set();
  std::vector<unsigned>
  find_nearest_facets(const stk::math::Vector3d & x, const double maxSearchRadius) const;

  FacetMesh myMesh;
  bool myHaveConnectedSurface;

  // Built by prepare_to_compute over the facets (owned and ghosted) within the local narrow band.
  std::unique_ptr<SearchTree<unsigned>> myTree;
  bool myRebuildOwnedFacetsEachPrepare{false};
  bool myHaveBuiltOwnedFacets{false};
  bool myHaveCheckedOwnedFacets{false};
  bool myRemoveCoincidentFacets{false};
  bool myCheckWatertight{false};
  std::optional<BoundingBox> myWatertightBbox;
  VertexIdsAndCoords myGivenVertexIdsAndCoords;
  FacetVertexIds myGivenFacetVertexIds;
  std::vector<stk::math::Vector3d> myGivenVertexVelocities;
  size_t myNumSearchedFacets{0};

  std::optional<double> myCosSharpFeatureAngle;
  bool myAreOpenBoundaryEdgesFeatures{false};
  FacetMeshSharpFeatures mySharpFeaturesFoundHere; // each found by one processor
  FacetMeshSharpFeatures myLocalSharpFeatures;
  std::unique_ptr<SearchTree<unsigned>> myFeatureEdgeTree;
  std::unique_ptr<SearchTree<unsigned>> myCornerTree;
};

// For each MeshFacetedSurface among the surfaces that checks watertightness but has no
// watertight bounding box, default its restriction box to bbox.  Other surfaces are ignored.
void default_watertight_bounding_box_from_domain(
    const std::vector<std::shared_ptr<Surface>> & surfaces, const BoundingBox & bbox);

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_MESHFACETEDSURFACE_HPP_ */
