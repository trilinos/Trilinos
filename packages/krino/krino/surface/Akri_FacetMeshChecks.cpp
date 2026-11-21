// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_FacetMeshChecks.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_SharpFeatureCalcs.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <algorithm>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

namespace krino
{

namespace
{

// The facet's vertex ids, sorted (with 0 after the 2 ids of a segment), and whether that is an even
// permutation of their order in the facet.
std::pair<std::array<FacetVertexId, 3>, bool>
sorted_vertex_ids(const FacetMesh & mesh, const unsigned facet)
{
  std::array<FacetVertexId, 3> ids{{0, 0, 0}};
  const unsigned dim = mesh.spatial_dimension();
  for (unsigned n = 0; n < dim; ++n)
    ids[n] = mesh.vertex_id(mesh.facet_vertices(facet)[n]);
  // Count the swaps of a bubble sort for the permutation parity.
  bool isEven = true;
  for (unsigned i = 0; i < dim; ++i)
    for (unsigned j = 0; j + 1 < dim - i; ++j)
      if (ids[j + 1] < ids[j])
      {
        std::swap(ids[j], ids[j + 1]);
        isEven = !isEven;
      }
  return {ids, isEven};
}

template <class KEY> size_t hash_key(const KEY & key)
{
  size_t h = 1469598103934665603ull; // FNV-1a
  for (auto value : key)
  {
    h ^= static_cast<size_t>(value);
    h *= 1099511628211ull;
  }
  return h;
}

template <class DATA> struct RendezvousFacet
{
  int proc;
  unsigned facet;
  DATA data;
};

// Each owned facet goes, with its key and data, to a rendezvous processor for its key, which gives
// each group of facets with the same key, ordered by (processor, facet) so the result does not
// depend on the order of arrival, to find_removed, and returns the ones it says to remove to their
// processors.  Returns the global number removed.
template <class KEY, class DATA, class KEY_AND_DATA, class FIND_REMOVED>
size_t remove_owned_facets_found_at_rendezvous(
    FacetMesh & mesh, const KEY_AND_DATA & key_and_data, const FIND_REMOVED & find_removed)
{
  const int numProcs = stk::parallel_machine_size(mesh.comm());
  stk::CommSparse toRendezvous(mesh.comm());
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for (unsigned facet = 0; facet < mesh.num_owned_facets(); ++facet)
        {
          const auto [key, data] = key_and_data(facet);
          stk::CommBuffer & buffer = toRendezvous.send_buffer(hash_key(key) % numProcs);
          buffer.pack(key);
          buffer.pack(data);
          buffer.pack(facet);
        }
      });

  std::map<KEY, std::vector<RendezvousFacet<DATA>>> groups;
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          KEY key;
          DATA data;
          unsigned facet = 0;
          buffer.unpack(key);
          buffer.unpack(data);
          buffer.unpack(facet);
          groups[key].push_back({proc, facet, data});
        }
      });

  stk::CommSparse fromRendezvous(mesh.comm());
  stk::pack_and_communicate(fromRendezvous,
      [&]()
      {
        for (auto & [key, members] : groups)
        {
          if (members.size() < 2) continue;
          std::sort(members.begin(),
              members.end(),
              [](const RendezvousFacet<DATA> & a, const RendezvousFacet<DATA> & b)
              { return std::tie(a.proc, a.facet) < std::tie(b.proc, b.facet); });
          const std::vector<bool> isRemoved = find_removed(members);
          for (size_t m = 0; m < members.size(); ++m)
            if (isRemoved[m]) fromRendezvous.send_buffer(members[m].proc).pack(members[m].facet);
        }
      });

  std::vector<bool> isRemoved(mesh.num_owned_facets(), false);
  size_t numRemoved = 0;
  stk::unpack_communications(fromRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = fromRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          unsigned facet = 0;
          buffer.unpack(facet);
          isRemoved[facet] = true;
          ++numRemoved;
        }
      });

  size_t globalNumRemoved = 0;
  stk::all_reduce_sum(mesh.comm(), &numRemoved, &globalNumRemoved, 1);
  if (globalNumRemoved > 0) mesh.remove_owned_facets(isRemoved);
  return globalNumRemoved;
}

// The facet with these vertex coordinates (3 per vertex, in order).
template <class FACET> FACET make_facet(const std::array<double, 9> & coords)
{
  auto vertex = [&](const unsigned v)
  { return stk::math::Vector3d(coords[3 * v], coords[3 * v + 1], coords[3 * v + 2]); };
  if constexpr (FACET::DIM == 3)
    return FACET(vertex(0), vertex(1), vertex(2));
  else
    return FACET(vertex(0), vertex(1));
}

template <class FACET> size_t remove_coincident_facets_by_location_for_dim(FacetMesh & mesh)
{
  // Each facet is sent with its vertices' coordinates, from which the rendezvous processor decides
  // as find_coincident_facets_to_remove does for all the facets together.
  auto key_and_coordinates = [&](const unsigned facet)
  {
    const FACET ownedFacet = mesh.make_facet<FACET>(facet);
    std::array<double, 9> coords{};
    for (int v = 0; v < FACET::DIM; ++v)
      for (int d = 0; d < 3; ++d)
        coords[3 * v + d] = ownedFacet.facet_vertex(v)[d];
    return std::pair{coincident_facet_key(ownedFacet), coords};
  };
  auto find_removed = [](const std::vector<RendezvousFacet<std::array<double, 9>>> & members)
  {
    std::vector<FACET> facets;
    facets.reserve(members.size());
    for (auto & member : members)
      facets.push_back(make_facet<FACET>(member.data));
    std::vector<const FACET *> group;
    for (auto & facet : facets)
      group.push_back(&facet);
    return find_coincident_group_facets_to_remove(group);
  };
  return remove_owned_facets_found_at_rendezvous<std::array<long long, 9>, std::array<double, 9>>(
      mesh, key_and_coordinates, find_removed);
}

template <class FACET>
std::vector<std::array<stk::math::Vector3d, 2>> find_non_watertight_edges_by_location_for_dim(
    const FacetMesh & mesh, const std::optional<BoundingBox> & restrictBbox)
{
  // The owned facets are made as they are used, rather than all stored.
  const std::function<FACET(size_t)> owned_facet = [&](const size_t facet)
  { return mesh.make_facet<FACET>(facet); };
  const double cellSize = compute_global_edge_quantization_cell_size<FACET>(
      mesh.comm(), mesh.num_owned_facets(), owned_facet);
  return find_non_watertight_edges<FACET>(
      mesh.num_owned_facets(), owned_facet, mesh.comm(), cellSize, restrictBbox);
}

} // namespace

size_t remove_coincident_facets_by_vertex_ids(FacetMesh & mesh)
{
  STK_ThrowRequireMsg(mesh.has_shared_vertices(),
      "remove_coincident_facets_by_vertex_ids needs shared vertices; use "
      "remove_coincident_facets_by_location.");
  // Each owned facet goes to a rendezvous processor for its sorted vertex ids, with whether they
  // are an even permutation of its vertices (its orientation).
  auto ids_and_orientation = [&](const unsigned facet) { return sorted_vertex_ids(mesh, facet); };
  auto find_removed = [](const std::vector<RendezvousFacet<bool>> & members)
  {
    // Keep one facet of the more common orientation, the first by (processor, facet), or none if
    // the orientations balance.
    const size_t numEven = std::count_if(
        members.begin(), members.end(), [](const RendezvousFacet<bool> & m) { return m.data; });
    const size_t numOdd = members.size() - numEven;
    std::vector<bool> isRemoved(members.size(), true);
    if (numEven != numOdd)
    {
      const auto kept = std::find_if(members.begin(),
          members.end(),
          [&](const RendezvousFacet<bool> & m) { return m.data == (numEven > numOdd); });
      isRemoved[kept - members.begin()] = false;
    }
    return isRemoved;
  };
  return remove_owned_facets_found_at_rendezvous<std::array<FacetVertexId, 3>, bool>(
      mesh, ids_and_orientation, find_removed);
}

std::vector<std::array<stk::math::Vector3d, 2>> find_non_watertight_edges_by_vertex_ids(
    const FacetMesh & mesh, const std::optional<BoundingBox> & restrictBbox)
{
  STK_ThrowRequireMsg(mesh.has_shared_vertices(),
      "find_non_watertight_edges_by_vertex_ids needs shared vertices; use "
      "find_non_watertight_edges_by_location.");
  // Each edge (sorted vertex ids, a vertex in 2D) of an owned facet goes to a rendezvous processor,
  // which counts its facets over all processors.
  const int numProcs = stk::parallel_machine_size(mesh.comm());
  const unsigned dim = mesh.spatial_dimension();
  auto for_each_owned_edge = [&](const auto & fn)
  {
    for (unsigned facet = 0; facet < mesh.num_owned_facets(); ++facet)
    {
      const auto & vertices = mesh.facet_vertices(facet);
      for (unsigned e = 0; e < dim; ++e)
      {
        // A 3D edge's two vertices, or a 2D segment's vertex (as both ends).
        std::array<unsigned, 2> edge{{vertices[e], vertices[(e + 1) % dim]}};
        if (dim == 2) edge[1] = edge[0];
        if (mesh.vertex_id(edge[1]) < mesh.vertex_id(edge[0])) std::swap(edge[0], edge[1]);
        fn(edge);
      }
    }
  };

  stk::CommSparse toRendezvous(mesh.comm());
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for_each_owned_edge(
            [&](const std::array<unsigned, 2> & edge)
            {
              const std::array<FacetVertexId, 3> key{
                  {mesh.vertex_id(edge[0]), mesh.vertex_id(edge[1]), 0}};
              stk::CommBuffer & buffer = toRendezvous.send_buffer(hash_key(key) % numProcs);
              buffer.pack(key[0]);
              buffer.pack(key[1]);
              buffer.pack(mesh.vertex_coordinates(edge[0]));
              buffer.pack(mesh.vertex_coordinates(edge[1]));
            });
      });

  struct EdgeCount
  {
    int count{0};
    std::array<stk::math::Vector3d, 2> coords;
  };
  std::map<std::pair<FacetVertexId, FacetVertexId>, EdgeCount> edgeCounts;
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          std::pair<FacetVertexId, FacetVertexId> key;
          std::array<stk::math::Vector3d, 2> coords;
          buffer.unpack(key.first);
          buffer.unpack(key.second);
          buffer.unpack(coords[0]);
          buffer.unpack(coords[1]);
          EdgeCount & edgeCount = edgeCounts[key];
          edgeCount.coords = coords;
          ++edgeCount.count;
        }
      });

  std::vector<std::array<stk::math::Vector3d, 2>> badEdges;
  for (auto & [key, edgeCount] : edgeCounts)
    if (edgeCount.count != 2 &&
        is_reportable_non_watertight_edge(edgeCount.coords, restrictBbox, 0.))
      badEdges.push_back(edgeCount.coords);
  return badEdges;
}

size_t remove_coincident_facets_by_location(FacetMesh & mesh)
{
  return (mesh.spatial_dimension() == 3) ?
      remove_coincident_facets_by_location_for_dim<Facet3d>(mesh) :
      remove_coincident_facets_by_location_for_dim<Facet2d>(mesh);
}

std::vector<std::array<stk::math::Vector3d, 2>> find_non_watertight_edges_by_location(
    const FacetMesh & mesh, const std::optional<BoundingBox> & restrictBbox)
{
  return (mesh.spatial_dimension() == 3) ?
      find_non_watertight_edges_by_location_for_dim<Facet3d>(mesh, restrictBbox) :
      find_non_watertight_edges_by_location_for_dim<Facet2d>(mesh, restrictBbox);
}

namespace
{

// A vertex and the other vertex of each of the segments (in 3D, the feature edges) it is in.
struct VertexAndOtherEnds
{
  stk::math::Vector3d coords;
  std::vector<stk::math::Vector3d> otherEnds;
};

// Collective.  Sends each (vertex id, coordinates, coordinates of the other end) from fill to the
// rendezvous processor for the vertex id, and returns those it receives, by vertex id.
std::map<FacetVertexId, VertexAndOtherEnds>
gather_vertex_other_ends_at_rendezvous(const stk::ParallelMachine comm,
    const std::function<void(const std::function<void(
            FacetVertexId, const stk::math::Vector3d &, const stk::math::Vector3d &)> &)> &
        for_each_vertex_and_other_end)
{
  const int numProcs = stk::parallel_machine_size(comm);
  stk::CommSparse toRendezvous(comm);
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for_each_vertex_and_other_end(
            [&](const FacetVertexId id,
                const stk::math::Vector3d & coords,
                const stk::math::Vector3d & otherEnd)
            {
              stk::CommBuffer & buffer =
                  toRendezvous.send_buffer(hash_key(std::array<FacetVertexId, 1>{{id}}) % numProcs);
              buffer.pack(id);
              buffer.pack(coords);
              buffer.pack(otherEnd);
            });
      });

  std::map<FacetVertexId, VertexAndOtherEnds> vertices;
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          FacetVertexId id = 0;
          stk::math::Vector3d coords, otherEnd;
          buffer.unpack(id);
          buffer.unpack(coords);
          buffer.unpack(otherEnd);
          VertexAndOtherEnds & vertex = vertices[id];
          vertex.coords = coords;
          vertex.otherEnds.push_back(otherEnd);
        }
      });
  return vertices;
}

bool is_sharp_between_segments(const VertexAndOtherEnds & vertex, const double cosFeatureAngle)
{
  return is_angle_sharper_than_feature_angle(
      cosine_of_angle_between_segments_at_vertex(
          vertex.coords, vertex.otherEnds[0], vertex.otherEnds[1]),
      cosFeatureAngle);
}

std::vector<FacetMeshSharpFeatures::Edge> find_feature_edges_3d(
    const FacetMesh & mesh, const double cosFeatureAngle, const bool areOpenBoundaryEdgesFeatures)
{
  // Each edge of an owned facet goes to a rendezvous processor for its sorted vertex ids, with the
  // facet's vertex off the edge.
  const int numProcs = stk::parallel_machine_size(mesh.comm());
  stk::CommSparse toRendezvous(mesh.comm());
  stk::pack_and_communicate(toRendezvous,
      [&]()
      {
        for (unsigned facet = 0; facet < mesh.num_owned_facets(); ++facet)
        {
          const auto facetVertices = mesh.facet_vertices(facet);
          for (unsigned e = 0; e < 3; ++e)
          {
            std::array<unsigned, 2> edge{{facetVertices[e], facetVertices[(e + 1) % 3]}};
            if (mesh.vertex_id(edge[1]) < mesh.vertex_id(edge[0])) std::swap(edge[0], edge[1]);
            const std::array<FacetVertexId, 2> key{
                {mesh.vertex_id(edge[0]), mesh.vertex_id(edge[1])}};
            stk::CommBuffer & buffer = toRendezvous.send_buffer(hash_key(key) % numProcs);
            buffer.pack(key);
            buffer.pack(mesh.vertex_coordinates(edge[0]));
            buffer.pack(mesh.vertex_coordinates(edge[1]));
            buffer.pack(mesh.vertex_coordinates(facetVertices[(e + 2) % 3]));
          }
        }
      });

  struct EdgeFacets
  {
    std::array<stk::math::Vector3d, 2> coords;
    std::vector<stk::math::Vector3d> verticesOffEdge;
  };
  std::map<std::array<FacetVertexId, 2>, EdgeFacets> edges;
  stk::unpack_communications(toRendezvous,
      [&](int proc)
      {
        stk::CommBuffer & buffer = toRendezvous.recv_buffer(proc);
        while (buffer.remaining())
        {
          std::array<FacetVertexId, 2> key;
          std::array<stk::math::Vector3d, 2> coords;
          stk::math::Vector3d vertexOffEdge;
          buffer.unpack(key);
          buffer.unpack(coords[0]);
          buffer.unpack(coords[1]);
          buffer.unpack(vertexOffEdge);
          EdgeFacets & edgeFacets = edges[key];
          edgeFacets.coords = coords;
          edgeFacets.verticesOffEdge.push_back(vertexOffEdge);
        }
      });

  std::vector<FacetMeshSharpFeatures::Edge> featureEdges;
  for (auto & [ids, edgeFacets] : edges)
  {
    const size_t numFacets = edgeFacets.verticesOffEdge.size();
    const bool isFeature = (numFacets == 1) ? areOpenBoundaryEdgesFeatures :
                                              (numFacets > 2) ||
            is_angle_sharper_than_feature_angle(
                cosine_of_angle_between_triangles_at_edge(edgeFacets.coords[0],
                    edgeFacets.coords[1],
                    edgeFacets.verticesOffEdge[0],
                    edgeFacets.verticesOffEdge[1]),
                cosFeatureAngle);
    if (isFeature) featureEdges.push_back({ids, edgeFacets.coords});
  }
  return featureEdges;
}

std::vector<FacetMeshSharpFeatures::Corner>
find_corners_of_feature_edges_3d(const stk::ParallelMachine comm,
    const std::vector<FacetMeshSharpFeatures::Edge> & featureEdges,
    const double cosFeatureAngle)
{
  const auto vertices = gather_vertex_other_ends_at_rendezvous(comm,
      [&](const auto & send)
      {
        for (auto & edge : featureEdges)
        {
          send(edge.vertexIds[0], edge.coords[0], edge.coords[1]);
          send(edge.vertexIds[1], edge.coords[1], edge.coords[0]);
        }
      });

  std::vector<FacetMeshSharpFeatures::Corner> corners;
  for (auto & [id, vertex] : vertices)
  {
    const size_t numFeatureEdges = vertex.otherEnds.size();
    if (numFeatureEdges != 2 || is_sharp_between_segments(vertex, cosFeatureAngle))
      corners.push_back({id, vertex.coords});
  }
  return corners;
}

std::vector<FacetMeshSharpFeatures::Corner> find_corners_2d(
    const FacetMesh & mesh, const double cosFeatureAngle, const bool areOpenBoundaryEdgesFeatures)
{
  const auto vertices = gather_vertex_other_ends_at_rendezvous(mesh.comm(),
      [&](const auto & send)
      {
        for (unsigned facet = 0; facet < mesh.num_owned_facets(); ++facet)
        {
          const auto segmentVertices = mesh.facet_vertices(facet);
          for (unsigned n = 0; n < 2; ++n)
            send(mesh.vertex_id(segmentVertices[n]),
                mesh.vertex_coordinates(segmentVertices[n]),
                mesh.vertex_coordinates(segmentVertices[1 - n]));
        }
      });

  std::vector<FacetMeshSharpFeatures::Corner> corners;
  for (auto & [id, vertex] : vertices)
  {
    const size_t numSegments = vertex.otherEnds.size();
    const bool isCorner = (numSegments == 1) ?
        areOpenBoundaryEdgesFeatures :
        (numSegments > 2) || is_sharp_between_segments(vertex, cosFeatureAngle);
    if (isCorner) corners.push_back({id, vertex.coords});
  }
  return corners;
}

} // namespace

FacetMeshSharpFeatures find_sharp_features_by_vertex_ids(
    const FacetMesh & mesh, const double cosFeatureAngle, const bool areOpenBoundaryEdgesFeatures)
{
  STK_ThrowRequireMsg(mesh.has_shared_vertices(),
      "find_sharp_features_by_vertex_ids needs shared vertices (a connected surface).");
  FacetMeshSharpFeatures features;
  if (mesh.spatial_dimension() == 3)
  {
    features.edges = find_feature_edges_3d(mesh, cosFeatureAngle, areOpenBoundaryEdgesFeatures);
    features.corners =
        find_corners_of_feature_edges_3d(mesh.comm(), features.edges, cosFeatureAngle);
  }
  else
  {
    features.corners = find_corners_2d(mesh, cosFeatureAngle, areOpenBoundaryEdgesFeatures);
  }
  return features;
}

} // namespace krino
