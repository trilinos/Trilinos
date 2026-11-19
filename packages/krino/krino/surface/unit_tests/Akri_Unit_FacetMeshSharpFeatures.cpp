// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_FacetMesh.hpp>
#include <Akri_FacetMeshChecks.hpp>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <set>
#include <vector>

namespace krino
{

namespace
{

double cos_of_degrees(const double degrees) { return std::cos(degrees * std::numbers::pi / 180.); }

// The facets with index i such that i % numProcs == rank, with ids (from 1) for the vertices.
FacetMesh build_mesh_with_facets_split_across_processors(const unsigned dim,
    const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & allFacets)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  std::set<unsigned> usedVertices;
  for (size_t i = 0; i < allFacets.size(); ++i)
  {
    if (int(i % numProcs) != me) continue;
    const auto & f = allFacets[i];
    facetVertexIds.push_back({{f[0] + 1u, f[1] + 1u, (dim == 3) ? f[2] + 1u : 0u}});
    for (unsigned n = 0; n < dim; ++n)
      usedVertices.insert(f[n]);
  }
  for (unsigned v : usedVertices)
    vertexIdsAndCoords.emplace_back(v + 1u, vertices[v]);
  FacetMesh mesh(MPI_COMM_WORLD, dim);
  mesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  return mesh;
}

// The feature edges found on all processors, by their (0-based) vertex indices, sorted.
std::vector<std::array<FacetVertexId, 2>>
global_feature_edges(const FacetMeshSharpFeatures & features)
{
  std::vector<std::array<FacetVertexId, 2>> local;
  for (auto & edge : features.edges)
    local.push_back({{edge.vertexIds[0] - 1, edge.vertexIds[1] - 1}});
  std::vector<std::array<FacetVertexId, 2>> global;
  stk::parallel_vector_concat(MPI_COMM_WORLD, local, global);
  std::sort(global.begin(), global.end());
  return global;
}

// The corners found on all processors, by their (0-based) vertex indices, sorted.
std::vector<FacetVertexId> global_corners(const FacetMeshSharpFeatures & features)
{
  std::vector<FacetVertexId> local;
  for (auto & corner : features.corners)
    local.push_back(corner.vertexId - 1);
  std::vector<FacetVertexId> global;
  stk::parallel_vector_concat(MPI_COMM_WORLD, local, global);
  std::sort(global.begin(), global.end());
  return global;
}

std::vector<std::array<FacetVertexId, 2>>
sorted_edges(std::vector<std::array<FacetVertexId, 2>> edges)
{
  for (auto & edge : edges)
    if (edge[1] < edge[0]) std::swap(edge[0], edge[1]);
  std::sort(edges.begin(), edges.end());
  return edges;
}

std::vector<FacetVertexId> all_indices(const size_t num)
{
  std::vector<FacetVertexId> indices(num);
  for (size_t i = 0; i < num; ++i)
    indices[i] = i;
  return indices;
}

// A closed prism over a counterclockwise footprint polygon, z in [0,1], with outward facets.  Its
// top and bottom are fans from footprint vertex fanApex, which sees the whole footprint.  Vertex i
// of the footprint is vertex i at the bottom and n+i at the top.
void fill_prism(const std::vector<std::array<double, 2>> & footprint,
    const unsigned fanApex,
    std::vector<stk::math::Vector3d> & vertices,
    std::vector<std::array<unsigned, 3>> & facets)
{
  const unsigned n = footprint.size();
  for (double z : {0., 1.})
    for (auto & xy : footprint)
      vertices.emplace_back(xy[0], xy[1], z);
  for (unsigned i = 1; i + 1 < n; ++i)
  {
    const unsigned a = (fanApex + i) % n;
    const unsigned b = (fanApex + i + 1) % n;
    facets.push_back({{n + fanApex, n + a, n + b}});
    facets.push_back({{fanApex, b, a}});
  }
  for (unsigned i = 0; i < n; ++i)
  {
    const unsigned j = (i + 1) % n;
    facets.push_back({{i, j, n + j}});
    facets.push_back({{i, n + j, n + i}});
  }
}

// The edges around the footprint at the bottom and top, and the vertical edges of a prism.
std::vector<std::array<FacetVertexId, 2>> prism_edges(const unsigned n)
{
  std::vector<std::array<FacetVertexId, 2>> edges;
  for (unsigned i = 0; i < n; ++i)
  {
    const unsigned j = (i + 1) % n;
    edges.push_back({{i, j}});
    edges.push_back({{n + i, n + j}});
    edges.push_back({{i, n + i}});
  }
  return sorted_edges(edges);
}

} // namespace

// A cube has its 12 edges as feature edges and its 8 vertices as corners at a feature angle of
// 135 degrees, and none at 60 degrees, sharper than its right angles.  The diagonals of its sides
// are flat.
TEST(FacetMeshSharpFeatures, cube_edgesAndCornersAtRightAngles)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> facets;
  fill_prism({{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}, 0, vertices, facets);
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(3, vertices, facets);

  const FacetMeshSharpFeatures features =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
  EXPECT_EQ(prism_edges(4), global_feature_edges(features));
  EXPECT_EQ(all_indices(8), global_corners(features));

  const FacetMeshSharpFeatures noFeatures =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(60.));
  EXPECT_TRUE(global_feature_edges(noFeatures).empty());
  EXPECT_TRUE(global_corners(noFeatures).empty());
}

// An L-shaped prism also has a concave (reentrant) vertical edge, which is as sharp as the convex
// ones.
TEST(FacetMeshSharpFeatures, lShapedPrism_convexAndConcaveEdges)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> facets;
  fill_prism({{0., 0.}, {2., 0.}, {2., 1.}, {1., 1.}, {1., 2.}, {0., 2.}}, 3, vertices, facets);
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(3, vertices, facets);

  const FacetMeshSharpFeatures features =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
  const auto featureEdges = global_feature_edges(features);
  EXPECT_EQ(prism_edges(6), featureEdges);
  EXPECT_TRUE(std::binary_search(featureEdges.begin(),
      featureEdges.end(),
      std::array<FacetVertexId, 2>{{3, 9}})); // the concave edge
  EXPECT_EQ(all_indices(12), global_corners(features));
}

// A faceted sphere has no facets meeting at less than 135 degrees.
TEST(FacetMeshSharpFeatures, sphere_noFeatures)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> facets;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, facets);
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(3, vertices, facets);

  const FacetMeshSharpFeatures features =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
  EXPECT_TRUE(global_feature_edges(features).empty());
  EXPECT_TRUE(global_corners(features).empty());
}

// The open boundary of a square plate is a feature only if asked for: then its 4 sides are feature
// edges and its vertices, where they meet at right angles, are corners.
TEST(FacetMeshSharpFeatures, openSquarePlate_boundaryEdgesOnlyIfAsked)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
  const std::vector<std::array<unsigned, 3>> facets{{{0, 1, 2}}, {{0, 2, 3}}};
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(3, vertices, facets);

  const FacetMeshSharpFeatures noFeatures =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
  EXPECT_TRUE(global_feature_edges(noFeatures).empty());
  EXPECT_TRUE(global_corners(noFeatures).empty());

  const FacetMeshSharpFeatures features =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.), true);
  EXPECT_EQ(sorted_edges({{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}), global_feature_edges(features));
  EXPECT_EQ(all_indices(4), global_corners(features));
}

// An edge of three facets is a feature edge, even where two of them are flat, and its ends are
// corners.
TEST(FacetMeshSharpFeatures, nonManifoldEdge_featureEdgeWithCornersAtItsEnds)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0.5, 1., 0.}, {0.5, -1., 0.}, {0.5, 0., 1.}};
  const std::vector<std::array<unsigned, 3>> facets{{{0, 1, 2}}, {{1, 0, 3}}, {{0, 1, 4}}};
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(3, vertices, facets);

  const FacetMeshSharpFeatures features =
      find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
  EXPECT_EQ(sorted_edges({{{0, 1}}}), global_feature_edges(features));
  EXPECT_EQ((std::vector<FacetVertexId>{0, 1}), global_corners(features));
}

// In 2D the corners of an L-shaped polygon, concave or convex, and of a square loop of segments
// with vertices at the middles of its sides, which are flat.
TEST(FacetMeshSharpFeatures, polygons2d_cornersAtSharpVertices)
{
  {
    const std::vector<stk::math::Vector3d> vertices{
        {0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}};
    std::vector<std::array<unsigned, 3>> segments;
    for (unsigned i = 0; i < 6; ++i)
      segments.push_back({{i, (i + 1) % 6, 0}});
    const FacetMesh mesh = build_mesh_with_facets_split_across_processors(2, vertices, segments);
    const FacetMeshSharpFeatures features =
        find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
    EXPECT_TRUE(features.edges.empty());
    EXPECT_EQ(all_indices(6), global_corners(features));
  }
  {
    const std::vector<stk::math::Vector3d> vertices{{0., 0., 0.},
        {0.5, 0., 0.},
        {1., 0., 0.},
        {1., 0.5, 0.},
        {1., 1., 0.},
        {0.5, 1., 0.},
        {0., 1., 0.},
        {0., 0.5, 0.}};
    std::vector<std::array<unsigned, 3>> segments;
    for (unsigned i = 0; i < 8; ++i)
      segments.push_back({{i, (i + 1) % 8, 0}});
    const FacetMesh mesh = build_mesh_with_facets_split_across_processors(2, vertices, segments);
    const FacetMeshSharpFeatures features =
        find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.));
    EXPECT_EQ((std::vector<FacetVertexId>{0, 2, 4, 6}), global_corners(features));
  }
}

// The ends of an open chain of segments are corners only if open boundaries are features.
TEST(FacetMeshSharpFeatures, openChain2d_endsAreCornersOnlyIfAsked)
{
  const std::vector<stk::math::Vector3d> vertices{{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}};
  const std::vector<std::array<unsigned, 3>> segments{{{0, 1, 0}}, {{1, 2, 0}}};
  const FacetMesh mesh = build_mesh_with_facets_split_across_processors(2, vertices, segments);
  EXPECT_TRUE(
      global_corners(find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.))).empty());
  EXPECT_EQ((std::vector<FacetVertexId>{0, 2}),
      global_corners(find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.), true)));
}

TEST(FacetMeshSharpFeatures, unconnectedFacets_throws)
{
  FacetMesh mesh(MPI_COMM_WORLD, 3, FacetMeshOptions{.sharedVertices = false});
  EXPECT_ANY_THROW(find_sharp_features_by_vertex_ids(mesh, cos_of_degrees(135.)));
}

} // namespace krino
