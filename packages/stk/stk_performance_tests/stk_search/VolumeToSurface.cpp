// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_unit_test_utils/MeshUtilsForBoundingVolumes.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <ArborX.hpp>

namespace {

template<typename BoxVectorType>
void run_volume_to_surface_test(const std::string& meshFileName,
                                stk::search::SearchMethod searchMethod,
                                bool enforceSearchResultSymmetry = true)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    BoxVectorType elemBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK, elemBoxes);

    BoxVectorType sideBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::FACE_RANK, sideBoxes);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      SearchResults searchResults;
      stk::search::coarse_search(elemBoxes, sideBoxes, searchMethod, comm, searchResults, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

template<typename BoxIdentProcType>
void run_volume_to_surface_test_with_views(const std::string& meshFileName,
                                           stk::search::SearchMethod searchMethod,
                                           bool enforceSearchResultSymmetry = true)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<BoxIdentProcType *, ExecSpace> elemBoxes = createBoundingBoxesForEntities<BoxIdentProcType>(*bulkPtr,
                                                                                              stk::topology::ELEM_RANK);

    Kokkos::View<BoxIdentProcType *, ExecSpace> sideBoxes = createBoundingBoxesForEntities<BoxIdentProcType>(*bulkPtr,
                                                                                              stk::topology::FACE_RANK);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      Kokkos::View<IdentProcIntersection*, ExecSpace> searchResults;
      stk::search::coarse_search(elemBoxes, sideBoxes, searchMethod, comm, searchResults,
                                 ExecSpace{}, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_KDTREE)
{
  run_volume_to_surface_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_MORTON_LBVH)
{
  run_volume_to_surface_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_ARBORX)
{
  run_volume_to_surface_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_with_views_MORTON_LBVH)
{
  run_volume_to_surface_test_with_views<FloatBoxIdentProc>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_with_views_ARBORX)
{
  run_volume_to_surface_test_with_views<FloatBoxIdentProc>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}


TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_KDTREE)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test<FloatBoxIdentProcVector>(meshFileName, stk::search::KDTREE);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_MORTON_LBVH)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test<FloatBoxIdentProcVector>(meshFileName, stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_ARBORX)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test<FloatBoxIdentProcVector>(meshFileName, stk::search::ARBORX);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_with_views_MORTON_LBVH)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_with_views<FloatBoxIdentProc>(meshFileName, stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_with_views_ARBORX)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_with_views<FloatBoxIdentProc>(meshFileName, stk::search::ARBORX);
}

template<typename BoxVectorType>
void run_volume_to_surface_test_local(const std::string& meshFileName,
                                      stk::search::SearchMethod searchMethod)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    BoxVectorType elemBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK, elemBoxes);

    BoxVectorType sideBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::FACE_RANK, sideBoxes);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      LocalSearchResults searchResults;
      stk::search::local_coarse_search(elemBoxes, sideBoxes, searchMethod, searchResults);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

template<typename BoxIdentType>
void run_volume_to_surface_test_local_with_views(const std::string& meshFileName,
                                                 stk::search::SearchMethod searchMethod)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<BoxIdentType *, ExecSpace> elemBoxes = createBoundingBoxesForEntities<BoxIdentType>(*bulkPtr,
                                                                                              stk::topology::ELEM_RANK);

    Kokkos::View<BoxIdentType *, ExecSpace> sideBoxes = createBoundingBoxesForEntities<BoxIdentType>(*bulkPtr,
                                                                                              stk::topology::FACE_RANK);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
      stk::search::local_coarse_search(elemBoxes, sideBoxes, searchMethod, searchResults, ExecSpace{});
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_local_KDTREE)
{
  run_volume_to_surface_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_local_MORTON_LBVH)
{
  run_volume_to_surface_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_local_ARBORX)
{
  run_volume_to_surface_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_local_with_views_MORTON_LBVH)
{
  run_volume_to_surface_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_local_with_views_ARBORX)
{
  run_volume_to_surface_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}


TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_local_KDTREE)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_local<FloatBoxIdentVector>(meshFileName, stk::search::KDTREE);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_local_MORTON_LBVH)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_local<FloatBoxIdentVector>(meshFileName, stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_local_ARBORX)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_local<FloatBoxIdentVector>(meshFileName, stk::search::ARBORX);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_local_with_views_MORTON_LBVH)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_local_with_views<FloatBoxIdent>(meshFileName, stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_local_with_views_ARBORX)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_volume_to_surface_test_local_with_views<FloatBoxIdent>(meshFileName, stk::search::ARBORX);
}

using ExecSpace = stk::ngp::ExecSpace;
using MemSpace = stk::ngp::ExecSpace::memory_space;

inline Kokkos::View<ArborX::Box *, MemSpace>
createArborXBoundingBoxesForEntities(const stk::mesh::BulkData &bulk,
                                     stk::mesh::EntityRank rank)
{
  stk::mesh::EntityVector entities;
  const bool sortById = true;
  stk::mesh::get_entities(bulk, rank, bulk.mesh_meta_data().locally_owned_part(), entities, sortById);

  size_t numberBoundingBoxes = entities.size();
  Kokkos::View<ArborX::Box *, MemSpace> boundingBoxes("ArborX BBs", numberBoundingBoxes);
  auto boundingBoxesHost = Kokkos::create_mirror_view(boundingBoxes);

  stk::mesh::FieldBase const * coords = bulk.mesh_meta_data().coordinate_field();

  std::vector<double> boxCoordinates(6);

  for (size_t i = 0; i < entities.size(); ++i) {
    unsigned num_nodes = bulk.num_nodes(entities[i]);
    std::vector<double> coordinates(3*num_nodes,0);
    const stk::mesh::Entity* nodes = bulk.begin_nodes(entities[i]);
    for (unsigned j = 0; j < num_nodes; ++j) {
      double* data = static_cast<double*>(stk::mesh::field_data(*coords, nodes[j]));
      coordinates[3*j] = data[0];
      coordinates[3*j+1] = data[1];
      coordinates[3*j+2] = data[2];
    }
    findBoundingBoxCoordinates(coordinates, boxCoordinates);
    ArborX::Point min_point(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2]);
    ArborX::Point max_point(boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]);
    boundingBoxesHost(i) = {min_point, max_point};
  }

  Kokkos::deep_copy(boundingBoxes, boundingBoxesHost);

  return boundingBoxes;
}

void arborx_coarse_search(Kokkos::View<ArborX::Box *, MemSpace> elemBoxes,
                          Kokkos::View<ArborX::Box *, MemSpace> sideBoxes)
{
  ExecSpace execSpace;

  ArborX::BVH<MemSpace> bvh{execSpace, elemBoxes};

  const int numQueries = sideBoxes.extent(0);
  Kokkos::View<ArborX::Intersects<ArborX::Box> *, MemSpace> queries(Kokkos::ViewAllocateWithoutInitializing("queries"), numQueries);

  Kokkos::parallel_for("setup_queries", Kokkos::RangePolicy<ExecSpace>(0, numQueries),
                       KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(sideBoxes(i)); });
  Kokkos::fence();

  Kokkos::View<int *, MemSpace> indices("Example::indices", 0);
  Kokkos::View<int *, MemSpace> offsets("Example::offsets", 0);

  bvh.query(execSpace, queries, indices, offsets);
}

void run_search_test_arborx(const std::string& meshFileName)
{
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<ArborX::Box *, MemSpace> elemBoxes = createArborXBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK);
    Kokkos::View<ArborX::Box *, MemSpace> sideBoxes = createArborXBoundingBoxesForEntities(*bulkPtr, stk::topology::FACE_RANK);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      arborx_coarse_search(elemBoxes, sideBoxes);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_rawARBORX) {
  run_search_test_arborx("generated:40x80x20|sideset:xXyYzZ");
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_rawARBORX) {

  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_search_test_arborx(meshFileName);
}

void distributed_arborx_coarse_search(Kokkos::View<ArborX::Box *, MemSpace> elemBoxes,
                                      Kokkos::View<ArborX::Box *, MemSpace> sideBoxes,
                                      MPI_Comm comm)
{

  ExecSpace execSpace;

  ArborX::DistributedTree<MemSpace> tree(comm, execSpace, elemBoxes);

  const int numQueries = sideBoxes.extent(0);
  Kokkos::View<ArborX::Intersects<ArborX::Box> *, MemSpace> queries(Kokkos::ViewAllocateWithoutInitializing("queries"), numQueries);

  Kokkos::parallel_for("setup_queries", Kokkos::RangePolicy<ExecSpace>(0, numQueries),
                       KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(sideBoxes(i)); });
  Kokkos::fence();

  Kokkos::View<ArborX::PairIndexRank *, MemSpace> values("indicesAndRanks", 0);
  Kokkos::View<int *, MemSpace> offsets("offsets", 0);

  tree.query(execSpace, queries, values, offsets);
}

void run_search_test_distributed_arborx(const std::string& meshFileName)
{
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<ArborX::Box *, MemSpace> elemBoxes = createArborXBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK);
    Kokkos::View<ArborX::Box *, MemSpace> sideBoxes = createArborXBoundingBoxesForEntities(*bulkPtr, stk::topology::FACE_RANK);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      distributed_arborx_coarse_search(elemBoxes, sideBoxes, MPI_COMM_WORLD);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToSurface, generatedMesh_floatBox_distributed_rawARBORX) {
  run_search_test_distributed_arborx("generated:40x80x20|sideset:xXyYzZ");
}

TEST(StkSearch_VolumeToSurface, casaMesh_floatBox_distributed_rawARBORX) {

  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") GTEST_SKIP();

  run_search_test_distributed_arborx(meshFileName);
}

} // namespace

