// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBox.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_OutputUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/parallel/Parallel.hpp>
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
// the ASCII STL grammar that STLSurface::read_ascii_facets expects. The header
// line starts with "solid", which is how STLSurface::read_header selects the
// ASCII path (the binary path is already covered by the sphere test, which uses
// write_stl -> a "binary solid ..." header).
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

// The same triangle expressed as a minimal ASCII PLY, matching the tokens that
// PLYSurface::read_header / read_points / read_facets parse.
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

// The same triangle in the FAC grammar that FACSurface::read_points /
// read_facets expect: point count, then <id X Y Z> per point, then facet
// count, then <id p1 p2 p3> per facet.
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
// Parallel_Facet_File_Reader::get_batch_size is pure arithmetic that
// partitions a global facet count into read batches. It is used by every
// file reader but was previously untested.
//--------------------------------------------------------------------------
TEST(ParallelFacetFileReader, getBatchSize_partitionCoversAllFacetsWithoutNegativeRemainder)
{
  // Note: get_batch_size does an all_reduce_sum over the local counts, so the
  // global total scales with the number of ranks. Running this at 1-4 procs
  // exercises that reduction as well as the arithmetic.
  const int numRanks = stk::parallel_machine_size(MPI_COMM_WORLD);

  for (const int localNumFacets : {0, 1, 100, 2048, 5000, 20000})
  {
    int batchSize = -1;
    int numBatches = -1;
    Parallel_Facet_File_Reader::get_batch_size(localNumFacets, batchSize, numBatches);

    const int globalNumFacets = localNumFacets * numRanks;

    EXPECT_GE(batchSize, 1) << "batch size must be positive for localNumFacets=" << localNumFacets;
    EXPECT_GE(numBatches, 1) << "there must be at least one batch for localNumFacets=" << localNumFacets;

    // The read loops consume the facets as min(batchSize, remaining) per batch,
    // so the batches must cover every facet and the final remainder must never
    // be negative (batch*batchSize <= globalNumFacets for the last batch).
    EXPECT_GE(numBatches * batchSize, globalNumFacets)
        << "batches must cover all facets for localNumFacets=" << localNumFacets;
    EXPECT_LE((numBatches - 1) * batchSize, globalNumFacets)
        << "last-batch remainder must be non-negative for localNumFacets=" << localNumFacets;
  }
}

//--------------------------------------------------------------------------
// ASCII STL reading. The existing sphere test covers the binary path only;
// this pins the ASCII branch of STLSurface::read_header / read_ascii_facets.
//--------------------------------------------------------------------------
TEST(STLSurface, asciiFile_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return; // file readers distribute facets in parallel; assert counts in serial

  ScopedTextFile file("unit_test_ascii.stl", ascii_stl_single_triangle());

  STLSurface surface(file.name());
  const BoundingBox bbox = surface.get_bounding_box();

  EXPECT_TRUE(bbox.valid());
  const auto & lo = bbox.get_min();
  const auto & hi = bbox.get_max();
  EXPECT_NEAR(lo[0], 0.0, 1.e-6);
  EXPECT_NEAR(lo[1], 0.0, 1.e-6);
  EXPECT_NEAR(lo[2], 0.0, 1.e-6);
  EXPECT_NEAR(hi[0], 1.0, 1.e-6);
  EXPECT_NEAR(hi[1], 1.0, 1.e-6);
  EXPECT_NEAR(hi[2], 0.0, 1.e-6);

  BoundingBox queryBbox = bbox;
  queryBbox.scale(2.0);
  surface.prepare_to_compute(0., queryBbox, 0.);
  EXPECT_EQ(1u, surface.size()) << "expected exactly one facet from the ASCII STL";
}

//--------------------------------------------------------------------------
// PLY reading. No prior coverage for PLYSurface.
//--------------------------------------------------------------------------
TEST(PLYSurface, asciiFile_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  ScopedTextFile file("unit_test.ply", ascii_ply_single_triangle());

  PLYSurface surface(sierra::Diag::sierraTimer(), file.name(), 1, stk::math::Vector3d(1.,1.,1.));
  const BoundingBox bbox = surface.get_bounding_box();

  EXPECT_TRUE(bbox.valid());
  EXPECT_NEAR(bbox.get_min()[0], 0.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[0], 1.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[1], 1.0, 1.e-6);

  BoundingBox queryBbox = bbox;
  queryBbox.scale(2.0);
  surface.prepare_to_compute(0., queryBbox, 0.);
  EXPECT_EQ(1u, surface.size()) << "expected exactly one facet from the PLY file";
}

//--------------------------------------------------------------------------
// Sign convention: reading with sign=-1 flips the facet orientation, which
// negates the signed distance relative to sign=+1.
//--------------------------------------------------------------------------
TEST(STLSurface, oppositeSign_negatesSignedDistance)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  ScopedTextFile file("unit_test_sign.stl", ascii_stl_single_triangle());

  STLSurface positive(sierra::Diag::sierraTimer(), file.name(), 1, stk::math::Vector3d(1.,1.,1.));
  STLSurface negative(sierra::Diag::sierraTimer(), file.name(), -1, stk::math::Vector3d(1.,1.,1.));

  BoundingBox bbox = positive.get_bounding_box();
  bbox.scale(3.0);
  positive.prepare_to_compute(0., bbox, 0.);
  negative.prepare_to_compute(0., bbox, 0.);

  // A point clearly off one side of the triangle's plane.
  const stk::math::Vector3d queryPt(0.25, 0.25, 0.5);
  const double posDist = positive.point_distance(queryPt, 0., 0., true);
  const double negDist = negative.point_distance(queryPt, 0., 0., true);

  EXPECT_GT(std::abs(posDist), 1.e-6);
  EXPECT_NEAR(posDist, -negDist, 1.e-6) << "flipping the sign should negate the signed distance";
}

//--------------------------------------------------------------------------
// Error handling: a nonexistent file and a malformed PLY header should throw.
//--------------------------------------------------------------------------
TEST(STLSurface, missingFile_throws)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  STLSurface surface("this_file_does_not_exist.stl");
  EXPECT_ANY_THROW(surface.get_bounding_box());
}

TEST(PLYSurface, malformedHeader_throws)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  // Missing the leading "ply" magic token that read_header requires.
  ScopedTextFile file("unit_test_bad.ply", "not_a_ply_file\nelement vertex 0\nend_header\n");

  PLYSurface surface(sierra::Diag::sierraTimer(), file.name(), 1, stk::math::Vector3d(1.,1.,1.));
  EXPECT_ANY_THROW(surface.get_bounding_box());
}

//--------------------------------------------------------------------------
// FAC reading. No prior coverage for FACSurface.
//--------------------------------------------------------------------------
TEST(FACSurface, asciiFile_singleTriangle_getExpectedFacetCountAndBoundingBox)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  ScopedTextFile file("unit_test.fac", fac_single_triangle());

  FACSurface surface(sierra::Diag::sierraTimer(), file.name(), 1, stk::math::Vector3d(1.,1.,1.));
  const BoundingBox bbox = surface.get_bounding_box();

  EXPECT_TRUE(bbox.valid());
  EXPECT_NEAR(bbox.get_min()[0], 0.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[0], 1.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[1], 1.0, 1.e-6);

  BoundingBox queryBbox = bbox;
  queryBbox.scale(2.0);
  surface.prepare_to_compute(0., queryBbox, 0.);
  EXPECT_EQ(1u, surface.size()) << "expected exactly one facet from the FAC file";
}

//--------------------------------------------------------------------------
// EXO reading. EXOSurface reads a 3D ExodusII file of 3-node (trishell3)
// elements -- exactly what krino's own write_facets produces. Rather than
// check in a binary .exo artifact, this test writes a one-triangle surface to
// exodus with write_facets, then reads it back with EXOSurface, so the write
// and read paths pin each other.
//--------------------------------------------------------------------------
TEST(EXOSurface, roundTripSingleTriangle_getExpectedFacetCountAndBoundingBox)
{
  if (1 != stk::parallel_machine_size(MPI_COMM_WORLD))
    return;

  // A bare 3D mesh is all write_facets needs -- it only reads the parallel
  // machine and spatial dimension off the BulkData to pick the 3d facet path.
  std::shared_ptr<stk::mesh::BulkData> mesh =
      stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();

  std::unique_ptr<FacetedSurfaceBase> facets = FacetedSurfaceBase::build(3);
  facets->emplace_back_3d(stk::math::Vector3d(0.,0.,0.),
                          stk::math::Vector3d(1.,0.,0.),
                          stk::math::Vector3d(0.,1.,0.));

  // write_facets transforms the base name (adds an index and .exo extension),
  // so the actual file on disk is what create_file_name reports.
  const std::string baseName = "unit_test_exo_surface";
  write_facets(*mesh, *facets, baseName, 0);
  const std::string writtenFile = create_file_name(baseName, 0);
  ASSERT_TRUE(does_file_exist(writtenFile)) << "write_facets did not produce " << writtenFile;

  EXOSurface surface(sierra::Diag::sierraTimer(), writtenFile, 1, stk::math::Vector3d(1.,1.,1.));
  const BoundingBox bbox = surface.get_bounding_box();

  EXPECT_TRUE(bbox.valid());
  EXPECT_NEAR(bbox.get_min()[0], 0.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[0], 1.0, 1.e-6);
  EXPECT_NEAR(bbox.get_max()[1], 1.0, 1.e-6);

  BoundingBox queryBbox = bbox;
  queryBbox.scale(2.0);
  surface.prepare_to_compute(0., queryBbox, 0.);
  EXPECT_EQ(1u, surface.size()) << "expected exactly one facet from the EXO round-trip";

  std::remove(writtenFile.c_str());
}

}
