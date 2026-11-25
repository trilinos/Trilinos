// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETMESH_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETMESH_HPP_

#include <Akri_BoundingBox.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <array>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

namespace krino
{

// The global id of a facet vertex: the same vertex has the same id on every processor.
using FacetVertexId = std::uint64_t;

// The kinds of vertices of a FacetMesh.  They must be the same on every processor.
struct FacetMeshOptions
{
  // Every vertex has a velocity, given with its coordinates and carried to the ghosts.
  bool withVertexVelocities = false;
  // Vertices are shared by id; otherwise each facet has its own, and there are no vertex ids.
  bool sharedVertices = true;
};

// Facets (triangles in 3D, segments in 2D) that share vertices: this processor's owned facets, then
// copies of nearby facets owned by other processors (ghosts).  A vertex id names one local vertex,
// used by every local facet (owned or ghosted) with that id.  Facets and vertices are referred to
// by local indices: the owned facets are 0 to num_owned_facets()-1, and the vertices of the owned
// facets are 0 to num_owned_facet_vertices()-1.  Changing the owned facets drops the ghosts.
// Without shared vertices, facet f's vertices are dim*f to dim*f+dim-1, and there are no vertex ids
// or vertex facets: vertex_id, vertex_facets, global_num_vertices and compute_vertex_owners throw.
class FacetMesh
{
public:
  using VertexIdsAndCoords = std::vector<std::pair<FacetVertexId, stk::math::Vector3d>>;
  using FacetVertexIds = std::vector<std::array<FacetVertexId, 3>>;

  FacetMesh(const stk::ParallelMachine comm,
      const unsigned dim,
      const FacetMeshOptions & options = FacetMeshOptions{});

  // Replaces all facets with these owned ones, each given by the ids of its vertices (the first dim
  // of each entry), with coordinates for every vertex they use.  With vertex velocities,
  // vertexVelocities[i] is the velocity of vertexIdsAndCoords[i]; otherwise it must be empty.
  void set_owned_facets(const VertexIdsAndCoords & vertexIdsAndCoords,
      const FacetVertexIds & facetVertexIds,
      const std::vector<stk::math::Vector3d> & vertexVelocities = {});
  // Replaces all facets with these owned ones, each given by the coordinates of its dim vertices in
  // order (and with vertex velocities, their velocities likewise).  With shared vertices, each
  // vertex gets a new id.
  void set_owned_facets_from_coordinates(const std::vector<stk::math::Vector3d> & facetVertexCoords,
      const std::vector<stk::math::Vector3d> & facetVertexVelocities = {});
  // Removes the owned facets flagged in isRemoved (indexed by owned facet), and the vertices only
  // they used.  The remaining facets keep their order.
  void remove_owned_facets(const std::vector<bool> & isRemoved);
  // Collective.  Replaces the ghosts with copies of the facets owned by other processors whose
  // bounding boxes intersect this processor's box in procBboxes (one box per processor).
  void update_ghost_facets(const std::vector<BoundingBox> & procBboxes);

  stk::ParallelMachine comm() const { return myComm; }
  unsigned spatial_dimension() const { return myDim; }
  size_t num_facets() const { return myFacetVertices.size(); }
  size_t num_owned_facets() const { return myNumOwnedFacets; }
  size_t num_ghosted_facets() const { return myFacetVertices.size() - myNumOwnedFacets; }
  bool is_owned(const unsigned facet) const { return facet < myNumOwnedFacets; }
  bool has_vertex_velocities() const { return myOptions.withVertexVelocities; }
  bool has_shared_vertices() const { return myOptions.sharedVertices; }
  // The local vertices (of owned and ghosted facets).
  size_t num_vertices() const { return myVertexCoords.size(); }
  // The vertices of the owned facets, which come first; the others are used only by ghosts.
  size_t num_owned_facet_vertices() const { return myNumOwnedFacetVertices; }
  // Collective.  The number of distinct vertices of the owned facets over all processors.
  size_t global_num_vertices() const;
  // Collective.  The owning processor of each vertex of the owned facets (0 to
  // num_owned_facet_vertices()-1): the lowest processor whose owned facets use it.
  std::vector<int> compute_vertex_owners() const;

  // The local vertices of the facet: 3 for a triangle, 2 for a segment.
  std::span<const unsigned> facet_vertices(const unsigned facet) const
  {
    return {myFacetVertices[facet].data(), myDim};
  }
  FacetVertexId vertex_id(const unsigned vertex) const
  {
    require_shared_vertices("vertex_id");
    return myVertexIds[vertex];
  }
  const stk::math::Vector3d & vertex_coordinates(const unsigned vertex) const
  {
    return myVertexCoords[vertex];
  }
  // The local facets (owned and ghosted) that use the vertex.
  std::span<const unsigned> vertex_facets(const unsigned vertex) const
  {
    require_shared_vertices("vertex_facets");
    return {myVertexFacets.data() + myVertexFacetOffsets[vertex],
        myVertexFacets.data() + myVertexFacetOffsets[vertex + 1]};
  }

  BoundingBox facet_bounding_box(const unsigned facet) const;
  stk::math::Vector3d facet_centroid(const unsigned facet) const;
  BoundingBox owned_facets_bounding_box() const;
  // The velocity of the vertex; requires has_vertex_velocities().
  const stk::math::Vector3d & vertex_velocity(const unsigned vertex) const
  {
    return myVertexVelocities[vertex];
  }

  // The facet as a Facet3d or Facet2d (a copy of its vertex coordinates), or as a
  // FacetWithVelocity3d or FacetWithVelocity2d (also a copy of its vertex velocities, which
  // requires has_vertex_velocities()).
  template <class FACET> FACET make_facet(const unsigned facet) const
  {
    if (FACET::DIM != myDim) throw_wrong_facet_dimension(FACET::DIM);
    const auto & v = myFacetVertices[facet];
    if constexpr (requires(const FACET & f) { f.get_velocity(); })
    {
      if (!myOptions.withVertexVelocities) throw_no_vertex_velocities();
      if constexpr (FACET::DIM == 3)
        return FACET(myVertexCoords[v[0]],
            myVertexCoords[v[1]],
            myVertexCoords[v[2]],
            myVertexVelocities[v[0]],
            myVertexVelocities[v[1]],
            myVertexVelocities[v[2]]);
      else
        return FACET(myVertexCoords[v[0]],
            myVertexCoords[v[1]],
            myVertexVelocities[v[0]],
            myVertexVelocities[v[1]]);
    }
    else if constexpr (FACET::DIM == 3)
      return FACET(myVertexCoords[v[0]], myVertexCoords[v[1]], myVertexCoords[v[2]]);
    else
      return FACET(myVertexCoords[v[0]], myVertexCoords[v[1]]);
  }

  // An estimate of the memory used.
  size_t storage_size() const;

private:
  unsigned add_vertex(const FacetVertexId id,
      const stk::math::Vector3d & coords,
      const stk::math::Vector3d & velocity);
  // A facet with its own vertices, appended after the current ones.
  void append_facet_with_own_vertices(
      const stk::math::Vector3d * coords, const stk::math::Vector3d * velocities);
  void build_vertex_facets();
  // The checks of the inline accessors; the throws are out of line, to keep the accessors small.
  void require_shared_vertices(const char * function) const
  {
    if (!myOptions.sharedVertices) throw_no_shared_vertices(function);
  }
  void throw_no_shared_vertices(const char * function) const;
  void throw_wrong_facet_dimension(const int facetDim) const;
  void throw_no_vertex_velocities() const;

  stk::ParallelMachine myComm;
  unsigned myDim;
  FacetMeshOptions myOptions;
  std::vector<FacetVertexId> myVertexIds; // empty without shared vertices
  std::vector<stk::math::Vector3d> myVertexCoords;
  std::vector<stk::math::Vector3d> myVertexVelocities; // empty without vertex velocities
  std::unordered_map<FacetVertexId, unsigned> myVertexIndexOfId;
  std::vector<std::array<unsigned, 3>> myFacetVertices;
  // myVertexFacets[myVertexFacetOffsets[v]] up to myVertexFacets[myVertexFacetOffsets[v+1]] are the
  // facets that use vertex v.
  std::vector<unsigned> myVertexFacetOffsets{0};
  std::vector<unsigned> myVertexFacets;
  size_t myNumOwnedFacets{0};
  size_t myNumOwnedFacetVertices{0}; // the vertices of the owned facets are the first ones
};

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_FACETMESH_HPP_ */
