// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_BoundingBox.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetMeshChecks.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <optional>
#include <set>
#include <vector>

namespace krino
{

namespace
{

size_t global_sum(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// The facets with index i such that i % numProcs == rank, with ids (from 1) for the vertices.
void fill_facets_for_this_processor(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & allFacets,
    FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords,
    FacetMesh::FacetVertexIds & facetVertexIds)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::set<unsigned> usedVertices;
  for (size_t i = 0; i < allFacets.size(); ++i)
  {
    if (int(i % numProcs) != me) continue;
    const auto & f = allFacets[i];
    facetVertexIds.push_back({{f[0] + 1u, f[1] + 1u, f[2] + 1u}});
    usedVertices.insert(f.begin(), f.end());
  }
  for (unsigned v : usedVertices)
    vertexIdsAndCoords.emplace_back(v + 1u, vertices[v]);
}

BoundingBox facet_bounding_box(
    const std::vector<stk::math::Vector3d> & vertices, const std::array<unsigned, 3> & facet)
{
  BoundingBox bbox;
  for (unsigned v : facet)
    bbox.accommodate(vertices[v]);
  return bbox;
}

} // namespace

TEST(FacetMesh, twoTrianglesSharingAnEdge_storeEachVertexOnce)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets(
      {{1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {0., 1., 0.}}, {4, {1., 1., 0.}}},
      {{{1, 2, 3}}, {{2, 4, 3}}});

  EXPECT_EQ(2u, mesh.num_facets());
  EXPECT_EQ(2u, mesh.num_owned_facets());
  EXPECT_EQ(0u, mesh.num_ghosted_facets());
  EXPECT_EQ(4u, mesh.num_vertices());
  EXPECT_EQ(4u, mesh.global_num_vertices());

  // Vertices are numbered in the order the facets first use them.
  EXPECT_EQ(2u, mesh.vertex_id(1));
  EXPECT_EQ(stk::math::Vector3d(1., 0., 0.), mesh.vertex_coordinates(1));
  const auto sharedVertexFacets = mesh.vertex_facets(1);
  EXPECT_EQ((std::vector<unsigned>{0, 1}),
      std::vector<unsigned>(sharedVertexFacets.begin(), sharedVertexFacets.end()));
  EXPECT_EQ(1u, mesh.vertex_facets(0).size());

  EXPECT_EQ(stk::math::Vector3d(2. / 3., 2. / 3., 0.), mesh.facet_centroid(1));
  EXPECT_TRUE(mesh.facet_bounding_box(1).contains(stk::math::Vector3d(1., 1., 0.)));
}

TEST(FacetMesh, missingVertexCoordinates_throws)
{
  FacetMesh mesh(MPI_COMM_WORLD, 3);
  EXPECT_ANY_THROW(mesh.set_owned_facets({{1, {0., 0., 0.}}, {2, {1., 0., 0.}}}, {{{1, 2, 3}}}));
}

// Each processor owns every numProcs-th facet of a sphere and asks for the facets intersecting its
// own x-slab: it gets each of them once (owned or ghosted), with each vertex id once.
TEST(FacetMesh, sphereSplitAcrossProcessors_ghostsEveryRemoteFacetIntersectingTheLocalBoxOnce)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<BoundingBox> slabs;
  for (int p = 0; p < numProcs; ++p)
  {
    const double width = 2. * radius / numProcs;
    slabs.emplace_back(stk::math::Vector3d(-radius + p * width, -1., -1.),
        stk::math::Vector3d(-radius + (p + 1) * width, 1., 1.));
  }

  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  fill_facets_for_this_processor(vertices, allFacets, vertexIdsAndCoords, facetVertexIds);

  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  EXPECT_EQ(vertices.size(), mesh.global_num_vertices());

  size_t expectedNumGhosts = 0;
  for (size_t i = 0; i < allFacets.size(); ++i)
    if (int(i % numProcs) != me && slabs[me].intersects(facet_bounding_box(vertices, allFacets[i])))
      ++expectedNumGhosts;

  for (int pass = 0; pass < 2; ++pass) // ghosts are replaced, not accumulated
  {
    mesh.update_ghost_facets(slabs);

    EXPECT_EQ(facetVertexIds.size(), mesh.num_owned_facets());
    EXPECT_EQ(expectedNumGhosts, mesh.num_ghosted_facets());

    std::set<std::array<FacetVertexId, 3>> ghostIds;
    for (unsigned f = mesh.num_owned_facets(); f < mesh.num_facets(); ++f)
    {
      EXPECT_TRUE(slabs[me].intersects(mesh.facet_bounding_box(f)));
      std::array<FacetVertexId, 3> ids;
      for (unsigned n = 0; n < 3; ++n)
        ids[n] = mesh.vertex_id(mesh.facet_vertices(f)[n]);
      ghostIds.insert(ids);
    }
    EXPECT_EQ(mesh.num_ghosted_facets(), ghostIds.size());

    std::set<FacetVertexId> vertexIds;
    for (unsigned v = 0; v < mesh.num_vertices(); ++v)
    {
      vertexIds.insert(mesh.vertex_id(v));
      EXPECT_EQ(vertices[mesh.vertex_id(v) - 1], mesh.vertex_coordinates(v));
      for (unsigned f : mesh.vertex_facets(v))
      {
        const auto & fv = mesh.facet_vertices(f);
        EXPECT_NE(fv.end(), std::find(fv.begin(), fv.end(), v));
      }
    }
    EXPECT_EQ(mesh.num_vertices(), vertexIds.size());
    EXPECT_EQ(vertices.size(), mesh.global_num_vertices());

    // The owned facets' vertices come first; the others are used only by ghosts.
    for (unsigned f = 0; f < mesh.num_owned_facets(); ++f)
      for (unsigned v : mesh.facet_vertices(f))
        EXPECT_LT(v, mesh.num_owned_facet_vertices());
    for (unsigned v = mesh.num_owned_facet_vertices(); v < mesh.num_vertices(); ++v)
      for (unsigned f : mesh.vertex_facets(v))
        EXPECT_FALSE(mesh.is_owned(f));
  }
}

// Every processor holds the same two segments, which share vertex 2.
TEST(FacetMesh, segmentsInTwoDimensions_facetVerticesHaveTwoEntries)
{
  FacetMesh mesh(MPI_COMM_WORLD, 2);
  mesh.set_owned_facets(
      {{1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {1., 1., 0.}}}, {{{1, 2, 0}}, {{2, 3, 0}}});

  EXPECT_EQ(3u, mesh.num_vertices());
  EXPECT_EQ(3u, mesh.num_owned_facet_vertices());
  EXPECT_EQ(3u, mesh.global_num_vertices());
  const std::array<std::array<FacetVertexId, 2>, 2> expectedFacetIds{{{{1, 2}}, {{2, 3}}}};
  for (unsigned f = 0; f < 2; ++f)
  {
    const auto facetVertices = mesh.facet_vertices(f);
    ASSERT_EQ(2u, facetVertices.size());
    for (unsigned n = 0; n < 2; ++n)
      EXPECT_EQ(expectedFacetIds[f][n], mesh.vertex_id(facetVertices[n]));
  }
  const auto sharedVertexFacets = mesh.vertex_facets(mesh.facet_vertices(0)[1]);
  EXPECT_EQ((std::vector<unsigned>{0, 1}),
      std::vector<unsigned>(sharedVertexFacets.begin(), sharedVertexFacets.end()));
  EXPECT_EQ(stk::math::Vector3d(1., 0.5, 0.), mesh.facet_centroid(1));
}

TEST(FacetMesh, makeFacetOfOtherDimension_throws)
{
  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets({{1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {0., 1., 0.}}}, {{{1, 2, 3}}});
  EXPECT_NO_THROW(mesh.make_facet<Facet3d>(0));
  EXPECT_ANY_THROW(mesh.make_facet<Facet2d>(0));
}

namespace
{

std::vector<Facet3d> facets_from(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & connectivity)
{
  std::vector<Facet3d> facets;
  for (auto & f : connectivity)
    facets.emplace_back(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
  return facets;
}

FacetMesh
build_mesh_with_facets_split_across_processors(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & allFacets)
{
  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  fill_facets_for_this_processor(vertices, allFacets, vertexIdsAndCoords, facetVertexIds);
  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  return mesh;
}

} // namespace

// Coincident facets (the same vertices), wherever they are owned, are reduced as
// remove_coincident_facets reduces disconnected facets.
TEST(FacetMesh, coincidentFacetsOnAnyProcessors_removedAsForDisconnectedFacets)
{
  const std::vector<stk::math::Vector3d> vertices{{0., 0., 0.},
      {1., 0., 0.},
      {0., 1., 0.},
      {1., 1., 0.},
      {2., 0., 0.},
      {2., 1., 0.},
      {3., 0., 0.}};
  const std::vector<std::array<unsigned, 3>> allFacets{{{0, 1, 2}},
      {{1, 2, 0}}, // the same orientation (a rotation): one stays
      {{1, 3, 2}},
      {{2, 3, 1}}, // opposite orientations: both go
      {{1, 4, 3}},
      {{4, 3, 1}},
      {{3, 4, 1}},  // two of one orientation and one of the other: one stays
      {{4, 6, 5}}}; // no partner: stays

  std::vector<Facet3d> disconnected = facets_from(vertices, allFacets);
  const size_t expectedNumRemoved = remove_coincident_facets(disconnected);
  ASSERT_EQ(5u, expectedNumRemoved);

  FacetMesh mesh = build_mesh_with_facets_split_across_processors(vertices, allFacets);
  EXPECT_EQ(expectedNumRemoved, remove_coincident_facets_by_vertex_ids(mesh));
  EXPECT_EQ(disconnected.size(), global_sum(mesh.num_owned_facets()));
  // The kept facets use every vertex: 0-1-2, one of 1-3-4, and 4-5-6.
  EXPECT_EQ(7u, mesh.global_num_vertices());
}

TEST(FacetMesh, closedAndOpenSurfaces_openEdgesAsForDisconnectedFacets)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> sphereFacets;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, sphereFacets);

  std::vector<std::array<unsigned, 3>> hemisphereFacets;
  for (auto & f : sphereFacets)
    if (vertices[f[0]][2] + vertices[f[1]][2] + vertices[f[2]][2] > 0.)
      hemisphereFacets.push_back(f);

  const BoundingBox lowerHalf(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 0.1));
  for (auto * facets : {&sphereFacets, &hemisphereFacets})
  {
    const std::vector<Facet3d> disconnected = facets_from(vertices, *facets);
    const double cellSize = compute_global_edge_quantization_cell_size(MPI_COMM_SELF, disconnected);
    FacetMesh mesh = build_mesh_with_facets_split_across_processors(vertices, *facets);
    for (const std::optional<BoundingBox> & restrictBbox :
        {std::optional<BoundingBox>{}, std::optional<BoundingBox>{lowerHalf}})
    {
      const size_t expected =
          find_non_watertight_edges(disconnected, MPI_COMM_SELF, cellSize, restrictBbox).size();
      EXPECT_EQ(
          expected, global_sum(find_non_watertight_edges_by_vertex_ids(mesh, restrictBbox).size()));
    }
  }
  EXPECT_EQ(0u,
      global_sum(find_non_watertight_edges_by_vertex_ids(
          build_mesh_with_facets_split_across_processors(vertices, sphereFacets), {})
              .size()));
  EXPECT_LT(0u,
      global_sum(find_non_watertight_edges_by_vertex_ids(
          build_mesh_with_facets_split_across_processors(vertices, hemisphereFacets), {})
              .size()));
}

namespace
{

// Each facet with its own vertices (ids 3f+1..3f+3), as read without merging.
FacetMesh build_unmerged_mesh(const std::vector<Facet3d> & facets)
{
  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  for (size_t f = 0; f < facets.size(); ++f)
  {
    std::array<FacetVertexId, 3> ids;
    for (unsigned n = 0; n < 3; ++n)
    {
      ids[n] = 3 * f + n + 1;
      vertexIdsAndCoords.emplace_back(ids[n], facets[f].facet_vertex(n));
    }
    facetVertexIds.push_back(ids);
  }
  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  return mesh;
}

} // namespace

// Without shared vertices, coincident facets and open edges are found by location, as for
// disconnected facets; by vertex ids every facet would be unique and every edge open.
TEST(FacetMesh, unmergedFacets_coincidentFacetsAndOpenEdgesFoundByLocation)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {2., 0., 0.}};
  const std::vector<std::array<unsigned, 3>> connectivity{
      {{0, 1, 2}}, {{1, 2, 0}}, {{1, 3, 2}}, {{2, 3, 1}}, {{1, 4, 3}}};
  std::vector<Facet3d> facets = facets_from(vertices, connectivity);
  const double cellSize = compute_global_edge_quantization_cell_size(MPI_COMM_SELF, facets);
  const size_t expectedNumOpenEdges =
      find_non_watertight_edges(facets, MPI_COMM_SELF, cellSize, {}).size();

  // The facets owned by processor 0 only (the checks are over all processors).
  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  {
    FacetMesh mesh = build_unmerged_mesh(isProc0 ? facets : std::vector<Facet3d>{});
    EXPECT_EQ(
        expectedNumOpenEdges, global_sum(find_non_watertight_edges_by_location(mesh, {}).size()));
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
      EXPECT_EQ(3 * connectivity.size(), find_non_watertight_edges_by_vertex_ids(mesh, {}).size());
    }
  }

  FacetMesh mesh = build_unmerged_mesh(isProc0 ? facets : std::vector<Facet3d>{});
  const size_t expectedNumRemoved = remove_coincident_facets(facets);
  EXPECT_EQ(expectedNumRemoved, remove_coincident_facets_by_location(mesh));
  EXPECT_EQ(facets.size(), global_sum(mesh.num_owned_facets()));
}

TEST(FacetMesh, removeOwnedFacets_keepsTheOthersInOrderWithTheirVertices)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets({{1, {0., 0., 0.}},
                            {2, {1., 0., 0.}},
                            {3, {0., 1., 0.}},
                            {4, {1., 1., 0.}},
                            {5, {2., 0., 0.}}},
      {{{1, 2, 3}}, {{2, 4, 3}}, {{2, 5, 4}}});
  mesh.remove_owned_facets({true, false, false}); // vertex 1 is then unused

  EXPECT_EQ(2u, mesh.num_owned_facets());
  EXPECT_EQ(4u, mesh.num_vertices());
  EXPECT_EQ(4u, mesh.global_num_vertices());
  const std::array<FacetVertexId, 3> expectedFacetIds[2] = {{{2, 4, 3}}, {{2, 5, 4}}};
  for (unsigned f = 0; f < 2; ++f)
    for (unsigned n = 0; n < 3; ++n)
    {
      const unsigned v = mesh.facet_vertices(f)[n];
      EXPECT_EQ(expectedFacetIds[f][n], mesh.vertex_id(v));
    }
  for (unsigned v = 0; v < mesh.num_vertices(); ++v)
    EXPECT_NE(1u, mesh.vertex_id(v)); // the removed facet's own vertex
  // Vertex 2 (the first kept vertex) is used by both remaining facets.
  const auto vertex2Facets = mesh.vertex_facets(mesh.facet_vertices(0)[0]);
  EXPECT_EQ((std::vector<unsigned>{0, 1}),
      std::vector<unsigned>(vertex2Facets.begin(), vertex2Facets.end()));
  EXPECT_EQ(stk::math::Vector3d(2., 0., 0.), mesh.vertex_coordinates(mesh.facet_vertices(1)[1]));
}

namespace
{

// The velocity of the vertex with id: a function of the id, so that every processor can check it.
stk::math::Vector3d velocity_of_vertex(const FacetVertexId id)
{
  const double x = static_cast<double>(id);
  return stk::math::Vector3d(x, 2. * x, -x);
}

std::vector<stk::math::Vector3d>
velocities_of_vertices(const FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords)
{
  std::vector<stk::math::Vector3d> velocities;
  for (auto & idAndCoords : vertexIdsAndCoords)
    velocities.push_back(velocity_of_vertex(idAndCoords.first));
  return velocities;
}

void expect_every_vertex_has_the_velocity_of_its_id(const FacetMesh & mesh)
{
  for (unsigned v = 0; v < mesh.num_vertices(); ++v)
    EXPECT_EQ(velocity_of_vertex(mesh.vertex_id(v)), mesh.vertex_velocity(v)) << "vertex " << v;
}

const FacetMesh::VertexIdsAndCoords squareVertices{
    {1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {0., 1., 0.}}, {4, {1., 1., 0.}}, {5, {2., 0., 0.}}};

} // namespace

TEST(FacetMesh, twoTrianglesWithVelocities_facetsWithVelocityHaveTheirVerticesVelocities)
{
  FacetMesh mesh(MPI_COMM_WORLD, 3, {.withVertexVelocities = true});
  ASSERT_TRUE(mesh.has_vertex_velocities());
  mesh.set_owned_facets(
      squareVertices, {{{1, 2, 3}}, {{2, 4, 3}}}, velocities_of_vertices(squareVertices));

  expect_every_vertex_has_the_velocity_of_its_id(mesh);
  for (unsigned f = 0; f < mesh.num_facets(); ++f)
  {
    const FacetWithVelocity3d facetWithVelocity = mesh.make_facet<FacetWithVelocity3d>(f);
    const Facet3d facet = mesh.make_facet<Facet3d>(f);
    for (unsigned n = 0; n < 3; ++n)
    {
      const unsigned v = mesh.facet_vertices(f)[n];
      EXPECT_EQ(velocity_of_vertex(mesh.vertex_id(v)), facetWithVelocity.get_velocity()[n]);
      EXPECT_EQ(mesh.vertex_coordinates(v), facetWithVelocity.facet_vertex(n));
      EXPECT_EQ(mesh.vertex_coordinates(v), facet.facet_vertex(n));
    }
  }
}

TEST(FacetMesh, segmentsWithVelocities_facetWithVelocity2dHasTheirVerticesVelocities)
{
  FacetMesh mesh(MPI_COMM_WORLD, 2, {.withVertexVelocities = true});
  const FacetMesh::VertexIdsAndCoords vertices{
      {1, {0., 0., 0.}}, {2, {1., 0., 0.}}, {3, {1., 1., 0.}}};
  mesh.set_owned_facets(vertices, {{{1, 2, 0}}, {{2, 3, 0}}}, velocities_of_vertices(vertices));

  const FacetWithVelocity2d facet = mesh.make_facet<FacetWithVelocity2d>(1);
  EXPECT_EQ(velocity_of_vertex(2), facet.get_velocity()[0]);
  EXPECT_EQ(velocity_of_vertex(3), facet.get_velocity()[1]);
}

// As sphereSplitAcrossProcessors_ghostsEveryRemoteFacetIntersectingTheLocalBoxOnce, with
// velocities: the ghosts' vertices carry them.
TEST(FacetMesh, sphereWithVelocitiesSplitAcrossProcessors_ghostsCarryTheVelocities)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  std::vector<BoundingBox> slabs;
  for (int p = 0; p < numProcs; ++p)
  {
    const double width = 2. * radius / numProcs;
    slabs.emplace_back(stk::math::Vector3d(-radius + p * width, -1., -1.),
        stk::math::Vector3d(-radius + (p + 1) * width, 1., 1.));
  }

  FacetMesh::VertexIdsAndCoords vertexIdsAndCoords;
  FacetMesh::FacetVertexIds facetVertexIds;
  fill_facets_for_this_processor(vertices, allFacets, vertexIdsAndCoords, facetVertexIds);
  FacetMesh mesh(MPI_COMM_WORLD, 3, {.withVertexVelocities = true});
  mesh.set_owned_facets(
      vertexIdsAndCoords, facetVertexIds, velocities_of_vertices(vertexIdsAndCoords));

  for (int pass = 0; pass < 2; ++pass)
  {
    mesh.update_ghost_facets(slabs);
    expect_every_vertex_has_the_velocity_of_its_id(mesh);
  }
  if (numProcs > 1)
  {
    EXPECT_LT(0u, global_sum(mesh.num_ghosted_facets()));
  }
}

TEST(FacetMesh, removeOwnedFacetsWithVelocities_keepsTheVelocitiesOfTheRemainingVertices)
{
  FacetMesh mesh(MPI_COMM_WORLD, 3, {.withVertexVelocities = true});
  mesh.set_owned_facets(squareVertices,
      {{{1, 2, 3}}, {{2, 4, 3}}, {{2, 5, 4}}},
      velocities_of_vertices(squareVertices));
  mesh.remove_owned_facets({true, false, false});

  EXPECT_EQ(4u, mesh.num_vertices());
  expect_every_vertex_has_the_velocity_of_its_id(mesh);
}

TEST(FacetMesh, velocitiesMismatched_throw)
{
  const FacetMesh::FacetVertexIds facets{{{1, 2, 3}}};
  const std::vector<stk::math::Vector3d> velocities = velocities_of_vertices(squareVertices);

  FacetMesh withoutVelocities(MPI_COMM_WORLD, 3);
  EXPECT_ANY_THROW(withoutVelocities.set_owned_facets(squareVertices, facets, velocities));
  withoutVelocities.set_owned_facets(squareVertices, facets);
  EXPECT_ANY_THROW(withoutVelocities.make_facet<FacetWithVelocity3d>(0));

  FacetMesh withVelocities(MPI_COMM_WORLD, 3, {.withVertexVelocities = true});
  EXPECT_ANY_THROW(withVelocities.set_owned_facets(squareVertices, facets));
  EXPECT_ANY_THROW(withVelocities.set_owned_facets(squareVertices,
      facets,
      std::vector<stk::math::Vector3d>(velocities.begin(), velocities.end() - 1)));
}

namespace
{

const FacetMeshOptions unshared{.sharedVertices = false};

// The coordinates of each facet's vertices, in order.
std::vector<stk::math::Vector3d>
facet_vertex_coordinates(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facets)
{
  std::vector<stk::math::Vector3d> coords;
  for (auto & facet : facets)
    for (unsigned n = 0; n < 3; ++n)
      coords.push_back(vertices[facet[n]]);
  return coords;
}

// The velocity of a vertex at x: a function of its position, so that every processor can check it.
stk::math::Vector3d velocity_at(const stk::math::Vector3d & x)
{
  return stk::math::Vector3d(x[1], -x[0], 2. * x[2]);
}

std::array<stk::math::Vector3d, 3> facet_coordinates(const FacetMesh & mesh, const unsigned facet)
{
  std::array<stk::math::Vector3d, 3> coords;
  for (unsigned n = 0; n < 3; ++n)
    coords[n] = mesh.vertex_coordinates(mesh.facet_vertices(facet)[n]);
  return coords;
}

// Orders facets by their vertices' coordinates.
struct FacetCoordinatesLess
{
  bool operator()(const std::array<stk::math::Vector3d, 3> & a,
      const std::array<stk::math::Vector3d, 3> & b) const
  {
    for (unsigned n = 0; n < 3; ++n)
      for (unsigned d = 0; d < 3; ++d)
        if (a[n][d] != b[n][d]) return a[n][d] < b[n][d];
    return false;
  }
};
using FacetCoordinatesSet = std::set<std::array<stk::math::Vector3d, 3>, FacetCoordinatesLess>;

} // namespace

// Two triangles sharing an edge, without shared vertices: each has its own 3 vertices.
TEST(FacetMesh, unconnectedTriangles_eachFacetHasItsOwnVertices)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  const std::vector<std::array<unsigned, 3>> facets{{{0, 1, 2}}, {{1, 3, 2}}};

  FacetMesh mesh(MPI_COMM_WORLD, 3, unshared);
  EXPECT_FALSE(mesh.has_shared_vertices());
  mesh.set_owned_facets_from_coordinates(facet_vertex_coordinates(vertices, facets));
  EXPECT_EQ(2u, mesh.num_owned_facets());
  EXPECT_EQ(6u, mesh.num_vertices());
  EXPECT_EQ(6u, mesh.num_owned_facet_vertices());
  for (unsigned f = 0; f < 2; ++f)
  {
    const Facet3d facet = mesh.make_facet<Facet3d>(f);
    for (unsigned n = 0; n < 3; ++n)
    {
      EXPECT_EQ(3 * f + n, mesh.facet_vertices(f)[n]);
      EXPECT_EQ(vertices[facets[f][n]], mesh.vertex_coordinates(3 * f + n));
      EXPECT_EQ(vertices[facets[f][n]], facet.facet_vertex(n));
    }
  }
  EXPECT_ANY_THROW(mesh.vertex_id(0));
  EXPECT_ANY_THROW(mesh.vertex_facets(0));
  EXPECT_ANY_THROW(mesh.global_num_vertices());

  // Shared ids given to an unconnected mesh are not shared: each facet still has its own vertices.
  mesh.set_owned_facets(squareVertices, {{{1, 2, 3}}, {{2, 4, 3}}});
  EXPECT_EQ(6u, mesh.num_vertices());
  EXPECT_EQ(stk::math::Vector3d(1., 0., 0.), mesh.vertex_coordinates(mesh.facet_vertices(1)[0]));
}

// As sphereSplitAcrossProcessors_ghostsEveryRemoteFacetIntersectingTheLocalBoxOnce, without shared
// vertices and with velocities: each ghost, found by its coordinates, once, with its velocities.
TEST(FacetMesh, unconnectedSphereSplitAcrossProcessors_ghostsEveryRemoteFacetOnceWithVelocities)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> allFacets;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, allFacets);
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<BoundingBox> slabs;
  for (int p = 0; p < numProcs; ++p)
  {
    const double width = 2. * radius / numProcs;
    slabs.emplace_back(stk::math::Vector3d(-radius + p * width, -1., -1.),
        stk::math::Vector3d(-radius + (p + 1) * width, 1., 1.));
  }

  std::vector<std::array<unsigned, 3>> ownedFacets;
  FacetCoordinatesSet expectedGhosts;
  for (size_t i = 0; i < allFacets.size(); ++i)
  {
    if (int(i % numProcs) == me)
      ownedFacets.push_back(allFacets[i]);
    else if (slabs[me].intersects(facet_bounding_box(vertices, allFacets[i])))
      expectedGhosts.insert(
          {{vertices[allFacets[i][0]], vertices[allFacets[i][1]], vertices[allFacets[i][2]]}});
  }
  const std::vector<stk::math::Vector3d> coords = facet_vertex_coordinates(vertices, ownedFacets);
  std::vector<stk::math::Vector3d> velocities;
  for (auto & x : coords)
    velocities.push_back(velocity_at(x));

  FacetMesh mesh(MPI_COMM_WORLD, 3, {.withVertexVelocities = true, .sharedVertices = false});
  mesh.set_owned_facets_from_coordinates(coords, velocities);
  for (int pass = 0; pass < 2; ++pass)
  {
    mesh.update_ghost_facets(slabs);
    EXPECT_EQ(ownedFacets.size(), mesh.num_owned_facets());
    EXPECT_EQ(3 * mesh.num_facets(), mesh.num_vertices());
    FacetCoordinatesSet ghosts;
    for (unsigned f = mesh.num_owned_facets(); f < mesh.num_facets(); ++f)
      ghosts.insert(facet_coordinates(mesh, f));
    EXPECT_EQ(mesh.num_ghosted_facets(), ghosts.size());
    EXPECT_TRUE(expectedGhosts == ghosts);
    for (unsigned v = 0; v < mesh.num_vertices(); ++v)
      EXPECT_EQ(velocity_at(mesh.vertex_coordinates(v)), mesh.vertex_velocity(v));
  }
}

TEST(FacetMesh, unconnectedRemoveOwnedFacets_keepsTheOthersInOrder)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {2., 0., 0.}};
  const std::vector<std::array<unsigned, 3>> facets{{{0, 1, 2}}, {{1, 3, 2}}, {{1, 4, 3}}};
  FacetMesh mesh(MPI_COMM_WORLD, 3, unshared);
  mesh.set_owned_facets_from_coordinates(facet_vertex_coordinates(vertices, facets));
  mesh.remove_owned_facets({false, true, false});

  ASSERT_EQ(2u, mesh.num_owned_facets());
  EXPECT_EQ(6u, mesh.num_vertices());
  for (unsigned n = 0; n < 3; ++n)
  {
    EXPECT_EQ(vertices[facets[0][n]], facet_coordinates(mesh, 0)[n]);
    EXPECT_EQ(vertices[facets[2][n]], facet_coordinates(mesh, 1)[n]);
    EXPECT_EQ(3 + n, mesh.facet_vertices(1)[n]);
  }
}

// Without shared vertices, the checks by vertex ids refuse the mesh, and the checks by location
// work as for disconnected facets.
TEST(FacetMesh, unconnectedByVertexIdChecks_throwAndByLocationChecksWork)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {2., 0., 0.}};
  const std::vector<std::array<unsigned, 3>> connectivity{
      {{0, 1, 2}}, {{1, 2, 0}}, {{1, 3, 2}}, {{2, 3, 1}}, {{1, 4, 3}}};
  std::vector<Facet3d> facets = facets_from(vertices, connectivity);
  const double cellSize = compute_global_edge_quantization_cell_size(MPI_COMM_SELF, facets);
  const size_t expectedNumOpenEdges =
      find_non_watertight_edges(facets, MPI_COMM_SELF, cellSize, {}).size();
  // Every processor holds all the facets, so they are coincident with each other's copies too.
  std::vector<Facet3d> allCopies;
  for (int p = 0; p < stk::parallel_machine_size(MPI_COMM_WORLD); ++p)
    allCopies.insert(allCopies.end(), facets.begin(), facets.end());
  const size_t expectedNumRemoved = remove_coincident_facets(allCopies);

  FacetMesh mesh(MPI_COMM_WORLD, 3, unshared);
  mesh.set_owned_facets_from_coordinates(facet_vertex_coordinates(vertices, connectivity));
  EXPECT_ANY_THROW(remove_coincident_facets_by_vertex_ids(mesh));
  EXPECT_ANY_THROW(find_non_watertight_edges_by_vertex_ids(mesh, {}));

  EXPECT_EQ(expectedNumRemoved, remove_coincident_facets_by_location(mesh));
  EXPECT_EQ(allCopies.size(), global_sum(mesh.num_owned_facets()));
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    FacetMesh openMesh(MPI_COMM_WORLD, 3, unshared);
    openMesh.set_owned_facets_from_coordinates(facet_vertex_coordinates(vertices, connectivity));
    EXPECT_EQ(expectedNumOpenEdges, find_non_watertight_edges_by_location(openMesh, {}).size());
  }
}

namespace
{

// The coordinates of every processor's owned facets, sorted.
std::vector<std::array<double, 9>> sorted_global_facet_coordinates(const FacetMesh & mesh)
{
  std::vector<std::array<double, 9>> local;
  for (unsigned f = 0; f < mesh.num_owned_facets(); ++f)
  {
    std::array<double, 9> coords{};
    const std::array<stk::math::Vector3d, 3> vertices = facet_coordinates(mesh, f);
    for (unsigned v = 0; v < 3; ++v)
      for (unsigned d = 0; d < 3; ++d)
        coords[3 * v + d] = vertices[v][d];
    local.push_back(coords);
  }
  std::vector<std::array<double, 9>> global;
  stk::parallel_vector_concat(MPI_COMM_WORLD, local, global);
  std::sort(global.begin(), global.end());
  return global;
}

} // namespace

// Coincident facets are removed by location as remove_coincident_facets removes them from all the
// facets together, however the facets are distributed over the processors and in whatever order:
// the same facets remain, including which of two nearly coincident facets is kept.
TEST(FacetMesh, unconnectedCoincidentFacetsOnAnyProcessors_removedAsForAllTogether)
{
  const stk::math::Vector3d nudge(1.e-10, 1.e-10, 1.e-10);
  const std::vector<std::array<stk::math::Vector3d, 3>> allFacets{
      {{{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}},
      {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}}, // the same orientation: one of the two stays
      {{{1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}}},
      {{{1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}}}, // opposite orientations: both go
      {{stk::math::Vector3d(2., 0., 0.) + nudge, {2., 1., 0.}, {1., 0., 0.}}},
      {{{2., 0., 0.}, {2., 1., 0.}, {1., 0., 0.}}}, // nearly coincident: the unnudged one stays
      {{{3., 0., 0.}, {3., 1., 0.}, {2., 0., 0.}}},
      {{{3., 0., 0.}, {2., 0., 0.}, {3., 1., 0.}}},
      {{{3., 1., 0.}, {2., 0., 0.}, {3., 0., 0.}}},  // two of one orientation, one of the other
      {{{4., 0., 0.}, {5., 0., 0.}, {4., 1., 0.}}}}; // no partner: stays

  std::vector<Facet3d> expected;
  for (auto & f : allFacets)
    expected.emplace_back(f[0], f[1], f[2]);
  const size_t expectedNumRemoved = remove_coincident_facets(expected);
  ASSERT_EQ(6u, expectedNumRemoved);
  FacetMesh expectedMesh(MPI_COMM_WORLD, 3, unshared);
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    std::vector<stk::math::Vector3d> coords;
    for (auto & facet : expected)
      for (unsigned n = 0; n < 3; ++n)
        coords.push_back(facet.facet_vertex(n));
    expectedMesh.set_owned_facets_from_coordinates(coords);
  }

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const size_t numFacets = allFacets.size();
  // The processor of facet i, and whether each processor takes its facets in reverse order.
  const std::vector<std::pair<std::function<int(size_t)>, bool>> layouts{
      {[&](size_t i) { return int(i % numProcs); }, false},
      {[&](size_t) { return numProcs - 1; }, false},
      {[&](size_t i) { return int((numFacets - 1 - i) % numProcs); }, true}};
  for (auto & [proc_of, isReversed] : layouts)
  {
    std::vector<stk::math::Vector3d> coords;
    for (size_t k = 0; k < numFacets; ++k)
    {
      const size_t i = isReversed ? numFacets - 1 - k : k;
      if (proc_of(i) == me) coords.insert(coords.end(), allFacets[i].begin(), allFacets[i].end());
    }
    FacetMesh mesh(MPI_COMM_WORLD, 3, unshared);
    mesh.set_owned_facets_from_coordinates(coords);
    EXPECT_EQ(expectedNumRemoved, remove_coincident_facets_by_location(mesh));
    EXPECT_EQ(sorted_global_facet_coordinates(expectedMesh), sorted_global_facet_coordinates(mesh));
  }
}

// The same sphere facets, each with its own vertices: unconnected storage takes at most half the
// memory of shared storage with distinct ids.
TEST(FacetMesh, unconnectedStorage_atMostHalfOfSharedStorageWithDistinctIds)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> facets;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, facets);
  const std::vector<stk::math::Vector3d> coords = facet_vertex_coordinates(vertices, facets);

  FacetMesh shared(MPI_COMM_WORLD, 3);
  shared.set_owned_facets_from_coordinates(coords);
  FacetMesh unconnected(MPI_COMM_WORLD, 3, unshared);
  unconnected.set_owned_facets_from_coordinates(coords);

  EXPECT_EQ(shared.num_vertices(), unconnected.num_vertices());
  EXPECT_LE(2 * unconnected.storage_size(), shared.storage_size());
}

// A fan of triangles (1, 2+r, 3+r), one owned by each processor r: vertex 1 is used by every
// processor and vertex 2+r by processors r-1 and r, so the owner of vertex 2+k is max(0, k-1).
TEST(FacetMesh, computeVertexOwners_lowestProcessorUsingEachVertex)
{
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const FacetVertexId first = 2 + me;
  FacetMesh mesh(MPI_COMM_WORLD, 3);
  mesh.set_owned_facets(
      {{1, {0., 0., 0.}}, {first, {1., double(me), 0.}}, {first + 1, {1., me + 1., 0.}}},
      {{{1, first, first + 1}}});

  const std::vector<int> owners = mesh.compute_vertex_owners();
  ASSERT_EQ(3u, owners.size());
  for (unsigned v = 0; v < 3; ++v)
  {
    const FacetVertexId id = mesh.vertex_id(v);
    const int expectedOwner = (id == 1) ? 0 : std::max(0, int(id) - 3);
    EXPECT_EQ(expectedOwner, owners[v]) << "vertex " << id;
  }

  FacetMesh unconnected(MPI_COMM_WORLD, 3, {.sharedVertices = false});
  EXPECT_ANY_THROW(unconnected.compute_vertex_owners());
}

} // namespace krino
