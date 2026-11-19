// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBox.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetFileSurface.hpp>
#include <Akri_FacetFileReader.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_OutputUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <Ionit_Initializer.h>
#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

// The EXO round-trip test opens an exodusII database, which requires the Ioss
// database types (exodus) to be registered. The main app does this at startup,
// but the unit-test main does not, so register them here.
Ioss::Init::Initializer initialize_ioss;

namespace
{

// A small RAII helper that writes text content to a file on rank 0 and removes
// it on destruction, so the surface-reader tests leave no artifacts behind.
class ScopedTextFile
{
public:
  ScopedTextFile(const std::string & filename, const std::string & contents)
  : myFilename(filename)
  {
    if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
    {
      std::ofstream out(myFilename);
      out << contents;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  ~ScopedTextFile()
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
      std::remove(myFilename.c_str());
  }

  const std::string & name() const { return myFilename; }

private:
  std::string myFilename;
};

// A single, non-degenerate, finite-normal triangle in the z=0 plane, written in
// ASCII STL.  The header line starts with "solid", which is how the STL reader
// selects the ASCII path (the binary path is covered by the sphere tests, which
// use write_stl).
std::string ascii_stl_single_triangle()
{
  return
    "solid unit_test\n"
    "facet normal 0 0 1\n"
    "  outer loop\n"
    "    vertex 0 0 0\n"
    "    vertex 1 0 0\n"
    "    vertex 0 1 0\n"
    "  endloop\n"
    "endfacet\n"
    "endsolid unit_test\n";
}

// The same triangle expressed as a minimal ASCII PLY.
std::string ascii_ply_single_triangle()
{
  return
    "ply\n"
    "format ascii 1.0\n"
    "element vertex 3\n"
    "property float x\n"
    "property float y\n"
    "property float z\n"
    "element face 1\n"
    "property list uchar int vertex_index\n"
    "end_header\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "3 0 1 2\n";
}

// The same triangle in the FAC grammar: point count, then <id X Y Z> per point,
// then facet count, then <id p1 p2 p3> per facet.
std::string fac_single_triangle()
{
  return
    "3\n"
    "0 0 0 0\n"
    "1 1 0 0\n"
    "2 0 1 0\n"
    "1\n"
    "0 0 1 2\n";
}

}

namespace krino
{

//--------------------------------------------------------------------------
// compute_facet_file_read_batch_size is pure arithmetic that partitions a global facet count into
// read batches.
//--------------------------------------------------------------------------
TEST(FacetFileReading, batchSize_partitionCoversAllFacetsWithoutNegativeRemainder)
{
  // The facet counts are summed over the processors, so the global total scales with the number of
  // ranks.  Running this at 1-4 procs exercises that reduction as well as the arithmetic.
  const int numRanks = stk::parallel_machine_size(MPI_COMM_WORLD);

  for (const size_t localNumFacets : {0, 1, 100, 2048, 5000, 20000})
  {
    size_t batchSize = 0;
    size_t numBatches = 0;
    compute_facet_file_read_batch_size(MPI_COMM_WORLD, localNumFacets, batchSize, numBatches);

    const size_t globalNumFacets = localNumFacets * numRanks;

    EXPECT_GE(batchSize, 1u) << "batch size must be positive for localNumFacets=" << localNumFacets;
    EXPECT_GE(numBatches, 1u) << "there must be at least one batch for localNumFacets="
                              << localNumFacets;

    // The read loops consume the facets as min(batchSize, remaining) per batch, so the batches must
    // cover every facet and the final remainder must never be negative.
    EXPECT_GE(numBatches * batchSize, globalNumFacets)
        << "batches must cover all facets for localNumFacets=" << localNumFacets;
    EXPECT_LE((numBatches - 1) * batchSize, globalNumFacets)
        << "last-batch remainder must be non-negative for localNumFacets=" << localNumFacets;
  }
}

size_t global_num_owned_surface_facets(const MeshFacetedSurface & surface)
{
  const size_t local = surface.num_owned_facets();
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// A single triangle with vertices (0,0,0), (1,0,0), (0,1,0) read from the file: its bounding box
// and facet count.
void expect_single_triangle(const std::string & format, const std::string & filename)
{
  FacetFileSurface surface(MPI_COMM_WORLD, format, filename, 1, stk::math::Vector3d(1., 1., 1.));
  BoundingBox bbox;
  surface.insert_into(bbox);

  EXPECT_TRUE(bbox.valid()) << format;
  EXPECT_NEAR(bbox.get_min()[0], 0.0, 1.e-6) << format;
  EXPECT_NEAR(bbox.get_min()[1], 0.0, 1.e-6) << format;
  EXPECT_NEAR(bbox.get_min()[2], 0.0, 1.e-6) << format;
  EXPECT_NEAR(bbox.get_max()[0], 1.0, 1.e-6) << format;
  EXPECT_NEAR(bbox.get_max()[1], 1.0, 1.e-6) << format;
  EXPECT_NEAR(bbox.get_max()[2], 0.0, 1.e-6) << format;

  BoundingBox queryBbox = bbox;
  queryBbox.scale(2.0);
  surface.prepare_to_compute(0., queryBbox, 0.);
  EXPECT_EQ(1u, global_num_owned_surface_facets(surface))
      << "expected exactly one facet from the " << format << " file";
}

TEST(FacetFileReading, asciiStl_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  ScopedTextFile file("unit_test_ascii.stl", ascii_stl_single_triangle());
  expect_single_triangle("STL", file.name());
}

TEST(FacetFileReading, ply_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  ScopedTextFile file("unit_test.ply", ascii_ply_single_triangle());
  expect_single_triangle("PLY", file.name());
}

TEST(FacetFileReading, fac_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  ScopedTextFile file("unit_test.fac", fac_single_triangle());
  expect_single_triangle("FAC", file.name());
}

// EXO files of 3-node (trishell3) elements are what krino's write_facets produces, so the write and
// read paths pin each other.  The file is written on processor 0 alone so it is not split.
TEST(FacetFileReading, exoRoundTrip_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  const std::string baseName = "unit_test_exo_surface";
  const std::string writtenFile = create_file_name(baseName, 0);
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    FacetMesh facets(MPI_COMM_SELF, 3, {.sharedVertices = false});
    facets.set_owned_facets_from_coordinates({stk::math::Vector3d(0., 0., 0.),
        stk::math::Vector3d(1., 0., 0.),
        stk::math::Vector3d(0., 1., 0.)});
    write_facets(facets, baseName, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  ASSERT_TRUE(does_file_exist(writtenFile)) << "write_facets did not produce " << writtenFile;

  expect_single_triangle("EXO", writtenFile);

  MPI_Barrier(MPI_COMM_WORLD);
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD)) std::remove(writtenFile.c_str());
}

// Reading with sign=-1 flips the facet orientation, which negates the signed distance.
TEST(FacetFileReading, oppositeSign_negatesSignedDistance)
{
  ScopedTextFile file("unit_test_sign.stl", ascii_stl_single_triangle());

  FacetFileSurface positive(MPI_COMM_WORLD, "STL", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  FacetFileSurface negative(
      MPI_COMM_WORLD, "STL", file.name(), -1, stk::math::Vector3d(1., 1., 1.));

  BoundingBox bbox;
  positive.insert_into(bbox);
  bbox.scale(3.0);
  positive.prepare_to_compute(0., bbox, 0.);
  negative.prepare_to_compute(0., bbox, 0.);

  // A point clearly off one side of the triangle's plane.
  const stk::math::Vector3d queryPt(0.25, 0.25, 0.5);
  const double posDist = positive.point_signed_distance(queryPt);
  const double negDist = negative.point_signed_distance(queryPt);

  EXPECT_GT(std::abs(posDist), 1.e-6);
  EXPECT_NEAR(posDist, -negDist, 1.e-6) << "flipping the sign should negate the signed distance";
}

// A nonexistent file and a malformed PLY header throw on every processor.
TEST(FacetFileReading, missingFile_throws)
{
  FacetFileSurface surface(
      MPI_COMM_WORLD, "STL", "this_file_does_not_exist.stl", 1, stk::math::Vector3d(1., 1., 1.));
  BoundingBox bbox;
  EXPECT_ANY_THROW(surface.insert_into(bbox));
}

// FAC points must be listed with ids 0, 1, ... in order (facets refer to them by id).
TEST(FacetFileReading, facPointsOutOfOrder_throws)
{
  ScopedTextFile file("unit_test_out_of_order.fac", "3\n0 0 0 0\n2 0 1 0\n1 1 0 0\n1\n0 0 1 2\n");

  FacetFileSurface surface(MPI_COMM_WORLD, "FAC", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  BoundingBox bbox;
  EXPECT_ANY_THROW(surface.insert_into(bbox));
}

TEST(FacetFileReading, malformedPlyHeader_throws)
{
  // Missing the leading "ply" magic token that the header requires.
  ScopedTextFile file("unit_test_bad.ply", "not_a_ply_file\nelement vertex 0\nend_header\n");

  FacetFileSurface surface(MPI_COMM_WORLD, "PLY", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  BoundingBox bbox;
  EXPECT_ANY_THROW(surface.insert_into(bbox));
}

}
