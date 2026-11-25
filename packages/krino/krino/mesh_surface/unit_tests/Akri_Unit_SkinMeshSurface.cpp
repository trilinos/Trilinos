// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_BoundingBox.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_Unit_MultiBlockMesh.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <algorithm>
#include <array>
#include <memory>
#include <set>
#include <vector>

namespace krino
{

namespace
{

size_t global_count(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

std::set<stk::mesh::EntityId> facet_node_ids(const FacetMesh & facets, const unsigned facet)
{
  std::set<stk::mesh::EntityId> ids;
  for (auto v : facets.facet_vertices(facet))
    ids.insert(facets.vertex_id(v));
  return ids;
}

stk::math::Vector3d facet_normal_and_centroid(const FacetMesh & facets, const unsigned facet, stk::math::Vector3d & centroid)
{
  const auto verts = facets.facet_vertices(facet);
  const stk::math::Vector3d x0 = facets.vertex_coordinates(verts[0]);
  const stk::math::Vector3d x1 = facets.vertex_coordinates(verts[1]);
  const stk::math::Vector3d x2 = facets.vertex_coordinates(verts[2]);
  centroid = (x0 + x1 + x2) / 3.;
  return Cross(x1 - x0, x2 - x0);
}

bool is_subset(const std::set<stk::mesh::EntityId> & ids, const std::set<stk::mesh::EntityId> & of)
{
  return std::includes(of.begin(), of.end(), ids.begin(), ids.end());
}

// Each owned facet is on an element of the skinned blocks (its vertices, other than a quadrilateral's
// centroid, are nodes of that element) and faces out of it, or into it if sign = -1.
void expect_facets_face_out_of_their_elements(const MultiBlockMesh & mesh,
    const std::vector<unsigned> & blockIndices,
    const MeshSurface & skin,
    const int sign)
{
  const FacetMesh & facets = skin.facet_mesh();
  for (unsigned f = 0; f < facets.num_owned_facets(); ++f)
  {
    std::set<stk::mesh::EntityId> nodeIds = facet_node_ids(facets, f);
    std::erase_if(nodeIds, [](const stk::mesh::EntityId id) { return id >= (stk::mesh::EntityId(1) << 63); });
    stk::math::Vector3d centroid;
    const stk::math::Vector3d normal = facet_normal_and_centroid(facets, f, centroid);
    unsigned numElements = 0;
    for (auto b : blockIndices)
    {
      if (!is_subset(nodeIds, mesh.myElemNodeIds[b])) continue;
      ++numElements;
      EXPECT_LT(0., sign * Dot(normal, centroid - mesh.myElemCentroids[b])) << "facet at " << centroid;
    }
    EXPECT_EQ(1u, numElements) << "facet at " << centroid;
  }
}

size_t global_count_of_facets_on_face(const MeshSurface & skin, const std::set<stk::mesh::EntityId> & faceNodeIds)
{
  const FacetMesh & facets = skin.facet_mesh();
  size_t count = 0;
  for (unsigned f = 0; f < facets.num_owned_facets(); ++f)
  {
    std::set<stk::mesh::EntityId> nodeIds = facet_node_ids(facets, f);
    std::erase_if(nodeIds, [](const stk::mesh::EntityId id) { return id >= (stk::mesh::EntityId(1) << 63); });
    if (is_subset(nodeIds, faceNodeIds)) ++count;
  }
  return global_count(count);
}

// The owned facets on the quadrilateral face, each as its sorted vertex ids, from every processor: those
// whose vertices are its nodes or its centroid.
std::vector<std::vector<FacetVertexId>> global_facets_on_quadrilateral(const MultiBlockMesh & mesh,
    const MeshSurface & skin,
    const std::set<stk::mesh::EntityId> & faceNodeIds)
{
  stk::math::Vector3d faceCentroid = stk::math::Vector3d::ZERO;
  const auto coordsData = mesh.myCoords->data();
  for (auto id : faceNodeIds)
    faceCentroid += 0.25 * stk::math::Vector3d(coordsData.entity_values(mesh.myMesh->get_entity(stk::topology::NODE_RANK, id)), 3);
  const FacetMesh & facets = skin.facet_mesh();
  auto is_vertex_of_face = [&](const unsigned v)
  {
    return faceNodeIds.count(facets.vertex_id(v)) > 0 || (facets.vertex_coordinates(v) - faceCentroid).length() < 1.e-12;
  };
  std::vector<FacetVertexId> local;
  for (unsigned f = 0; f < facets.num_owned_facets(); ++f)
  {
    const auto verts = facets.facet_vertices(f);
    if (!std::all_of(verts.begin(), verts.end(), is_vertex_of_face)) continue;
    std::vector<FacetVertexId> ids;
    for (auto v : facets.facet_vertices(f))
      ids.push_back(facets.vertex_id(v));
    std::sort(ids.begin(), ids.end());
    local.insert(local.end(), ids.begin(), ids.end());
  }
  int numProcs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const int localCount = int(local.size());
  std::vector<int> counts(numProcs), offsets(numProcs, 0);
  MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  for (int p = 1; p < numProcs; ++p)
    offsets[p] = offsets[p - 1] + counts[p - 1];
  std::vector<FacetVertexId> all(offsets.back() + counts.back());
  MPI_Allgatherv(local.data(), localCount, MPI_UINT64_T, all.data(), counts.data(), offsets.data(), MPI_UINT64_T, MPI_COMM_WORLD);
  std::vector<std::vector<FacetVertexId>> result;
  for (size_t i = 0; i < all.size(); i += 3)
    result.push_back({all[i], all[i + 1], all[i + 2]});
  std::sort(result.begin(), result.end());
  return result;
}

} // namespace

// The skin of each tet's block is its 4 faces facing out of it, or into it with sign = -1, including the
// face between the blocks, which faces out of each block in its skin.
TEST(SkinMeshSurface, twoTetBlocks_skinOfEachBlockIsItsFacesFacingOutOfIt)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  const MultiBlockMesh mesh = two_tet_blocks();
  for (const unsigned b : {0u, 1u})
    for (const int sign : {1, -1})
    {
      const auto skin = mesh.build_skin({b}, sign);
      EXPECT_EQ(4u, global_count(skin->num_owned_facets()));
      EXPECT_EQ(4u, skin->global_num_vertices());
      expect_facets_face_out_of_their_elements(mesh, {b}, *skin, sign);
      EXPECT_EQ(1u, global_count_of_facets_on_face(*skin, {2, 3, 4}));
      EXPECT_GT(0., sign * skin->point_signed_distance(mesh.myElemCentroids[b]));
      EXPECT_LT(0., sign * skin->point_signed_distance(mesh.myElemCentroids[1 - b]));
    }
}

// The skin of both blocks is the union's 6 exposed faces, without the face between them.
TEST(SkinMeshSurface, twoTetBlocks_skinOfBothBlocksLeavesOutTheFaceBetweenThem)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  const MultiBlockMesh mesh = two_tet_blocks();
  const auto skin = mesh.build_skin({0, 1}, 1);
  EXPECT_EQ(6u, global_count(skin->num_owned_facets()));
  EXPECT_EQ(5u, skin->global_num_vertices());
  expect_facets_face_out_of_their_elements(mesh, {0, 1}, *skin, 1);
  EXPECT_EQ(0u, global_count_of_facets_on_face(*skin, {2, 3, 4}));
}

// Both blocks' skins split the quadrilateral between them into the same triangles, planar (2 about a
// diagonal) or not (4 about its centroid, with the same id in both), so their facets there match.
TEST(SkinMeshSurface, twoHexBlocks_quadrilateralBetweenThemSplitAlikeInBothSkins)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  for (const double sharedNodeShift : {0., 0.2})
  {
    const MultiBlockMesh mesh = two_hex_blocks(sharedNodeShift);
    const auto skin1 = mesh.build_skin({0}, 1);
    const auto skin2 = mesh.build_skin({1}, 1);
    expect_facets_face_out_of_their_elements(mesh, {0}, *skin1, 1);
    expect_facets_face_out_of_their_elements(mesh, {1}, *skin2, 1);
    const std::set<stk::mesh::EntityId> sharedFace{2, 3, 8, 9};
    const auto sharedFacets1 = global_facets_on_quadrilateral(mesh, *skin1, sharedFace);
    EXPECT_EQ((sharedNodeShift == 0.) ? 2u : 4u, sharedFacets1.size()) << "shift " << sharedNodeShift;
    EXPECT_EQ(sharedFacets1, global_facets_on_quadrilateral(mesh, *skin2, sharedFace)) << "shift " << sharedNodeShift;
  }
}

} // namespace krino

