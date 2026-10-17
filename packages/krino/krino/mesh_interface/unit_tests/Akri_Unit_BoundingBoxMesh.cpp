// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBoxMesh.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <vector>

namespace krino {

namespace {

// A unit cube [0,1]^3 (or unit square in 2D) meshed at mesh_size == 1/N produces
// an N x N (x N) lattice of cells, i.e. (N+1)^dim nodes, independent of decomposition.
std::vector<size_t> global_entity_counts(const stk::mesh::BulkData & mesh)
{
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);
  return counts;
}

void expect_coordinates_within_domain(const stk::mesh::BulkData & mesh,
    const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const unsigned dim = meta.spatial_dimension();
  const stk::mesh::FieldBase & coordsField = *meta.coordinate_field();
  const double tol = 1.e-9;
  for (const auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (const auto node : *bucketPtr)
    {
      const double * coords = static_cast<const double*>(stk::mesh::field_data(coordsField, node));
      for (unsigned i = 0; i < dim; ++i)
      {
        EXPECT_GE(coords[i], minCorner[i] - tol);
        EXPECT_LE(coords[i], maxCorner[i] + tol);
      }
    }
  }
}

} // namespace

TEST(BoundingBoxMesh, Hex8_UnitCube_ProducesExpectedLattice)
{
  const size_t N = 4;
  const stk::math::Vector3d minCorner(0., 0., 0.);
  const stk::math::Vector3d maxCorner(1., 1., 1.);

  BoundingBoxMesh bboxMesh(stk::topology::HEXAHEDRON_8, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  bboxMesh.set_domain(BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), 1.0/N);
  bboxMesh.populate_mesh();

  EXPECT_EQ(3u, bboxMesh.meta_data().spatial_dimension());
  EXPECT_TRUE(bboxMesh.has_flat_boundaries());

  const auto counts = global_entity_counts(bboxMesh.bulk_data());
  EXPECT_EQ((N+1)*(N+1)*(N+1), counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(N*N*N,             counts[stk::topology::ELEM_RANK]);

  expect_coordinates_within_domain(bboxMesh.bulk_data(), minCorner, maxCorner);
}

TEST(BoundingBoxMesh, Tet4_UnitCube_ProducesExpectedNodeLattice)
{
  const size_t N = 3;
  const stk::math::Vector3d minCorner(0., 0., 0.);
  const stk::math::Vector3d maxCorner(1., 1., 1.);

  BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  bboxMesh.set_domain(BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), 1.0/N);
  bboxMesh.populate_mesh();

  EXPECT_EQ(3u, bboxMesh.meta_data().spatial_dimension());

  // The default structure for tets is BCC (not flat-walled), which augments the
  // (N+1)^3 cell-corner lattice with body-centered nodes, so assert on the
  // corner-lattice lower bound and on a non-zero element count.
  const auto counts = global_entity_counts(bboxMesh.bulk_data());
  EXPECT_GE(counts[stk::topology::NODE_RANK], (N+1)*(N+1)*(N+1));
  EXPECT_GT(counts[stk::topology::ELEM_RANK], 0u);
}

TEST(BoundingBoxMesh, Quad4_UnitSquare_Is2DWithExpectedLattice)
{
  const size_t N = 5;
  const stk::math::Vector3d minCorner(0., 0., 0.);
  const stk::math::Vector3d maxCorner(1., 1., 0.);

  BoundingBoxMesh bboxMesh(stk::topology::QUADRILATERAL_4_2D, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  bboxMesh.set_domain(BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), 1.0/N);
  bboxMesh.populate_mesh();

  EXPECT_EQ(2u, bboxMesh.meta_data().spatial_dimension());
  EXPECT_TRUE(bboxMesh.has_flat_boundaries());

  const auto counts = global_entity_counts(bboxMesh.bulk_data());
  EXPECT_EQ((N+1)*(N+1), counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(N*N,         counts[stk::topology::ELEM_RANK]);

  expect_coordinates_within_domain(bboxMesh.bulk_data(), minCorner, maxCorner);
}

TEST(BoundingBoxMesh, StructureTypeControlsFlatBoundaries)
{
  BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD

  bboxMesh.set_mesh_structure_type(CUBIC_BOUNDING_BOX_MESH);
  EXPECT_TRUE(bboxMesh.has_flat_boundaries());

  bboxMesh.set_mesh_structure_type(BCC_BOUNDING_BOX_MESH);
  EXPECT_FALSE(bboxMesh.has_flat_boundaries());

  bboxMesh.set_mesh_structure_type(FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  EXPECT_TRUE(bboxMesh.has_flat_boundaries());
}

// The recursive-coordinate-bisection layout gives each rank a compact rectangular
// brick of cells rather than a 1D slab. On 4+ ranks of a cube, the bisection splits
// two distinct axes, so a rank's locally-owned nodes span strictly less than the full
// domain in at least two dimensions. A 1D slab layout would span the full domain in
// two of the three dimensions and fail this check.
TEST(BoundingBoxMesh, Hex8_RcbProducesCompactSubdomains)
{
  const int pSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (pSize < 4)
    return; // decomposition shape is only distinguishable with several ranks

  const size_t N = 8;
  const stk::math::Vector3d minCorner(0., 0., 0.);
  const stk::math::Vector3d maxCorner(1., 1., 1.);

  BoundingBoxMesh bboxMesh(stk::topology::HEXAHEDRON_8, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  bboxMesh.set_domain(BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), 1.0/N);
  bboxMesh.populate_mesh();

  const stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const stk::mesh::FieldBase & coordsField = *meta.coordinate_field();

  stk::math::Vector3d localMin(1.e30, 1.e30, 1.e30);
  stk::math::Vector3d localMax(-1.e30, -1.e30, -1.e30);
  const stk::mesh::Selector ownedSelector = meta.locally_owned_part();
  for (const auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedSelector))
  {
    for (const auto node : *bucketPtr)
    {
      const double * coords = static_cast<const double*>(stk::mesh::field_data(coordsField, node));
      for (unsigned i = 0; i < 3; ++i)
      {
        localMin[i] = std::min(localMin[i], coords[i]);
        localMax[i] = std::max(localMax[i], coords[i]);
      }
    }
  }

  const double fullSpan = 1.0;
  const double tol = 1.e-6;
  int numConfinedDims = 0;
  for (unsigned i = 0; i < 3; ++i)
    if ((localMax[i] - localMin[i]) < fullSpan - tol)
      ++numConfinedDims;

  EXPECT_GE(numConfinedDims, 2)
      << "rank " << stk::parallel_machine_rank(MPI_COMM_WORLD)
      << " owns a subdomain spanning the full domain in >1 dimension; "
      << "expected a compact RCB brick";
}

// Exercises MeshInterface polymorphism: drive a BoundingBoxMesh through a base
// MeshInterface& and confirm meta/bulk accessors resolve to the populated mesh.
TEST(MeshInterface, BoundingBoxMeshThroughBasePointer)
{
  const size_t N = 2;
  BoundingBoxMesh bboxMesh(stk::topology::HEXAHEDRON_8, MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  bboxMesh.set_domain(BoundingBoxMesh::BoundingBoxType(stk::math::Vector3d(0.,0.,0.), stk::math::Vector3d(1.,1.,1.)), 1.0/N);

  MeshInterface & meshInterface = bboxMesh;
  meshInterface.populate_mesh();

  EXPECT_EQ(3u, meshInterface.meta_data().spatial_dimension());
  EXPECT_EQ(&bboxMesh.bulk_data(), &meshInterface.bulk_data());

  const auto counts = global_entity_counts(meshInterface.bulk_data());
  EXPECT_EQ((N+1)*(N+1)*(N+1), counts[stk::topology::NODE_RANK]);
}

} // namespace krino
