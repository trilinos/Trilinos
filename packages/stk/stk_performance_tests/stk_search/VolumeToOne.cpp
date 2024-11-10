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
#include "stk_search/SearchMethod.hpp"
#include <Kokkos_Core.hpp>
#include "stk_util/ngp/NgpSpaces.hpp"

namespace {

template<typename BoxVectorType>
void run_volume_to_one_test(const std::string& meshFileName,
                            stk::search::SearchMethod searchMethod,
                            bool enforceSearchResultSymmetry = true)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  const int pRank = stk::parallel_machine_rank(comm);
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    BoxVectorType elemBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK, elemBoxes);

    BoxVectorType supersetBoxVec { elemBoxes[0] };
    auto & [supersetBox, supersetIdent] = supersetBoxVec[0];
    supersetIdent.set_id(pRank);
    supersetIdent.set_proc(pRank);

    for (const auto & [box, ident] : elemBoxes) {
      stk::search::add_to_box(supersetBox, box);
    }

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      SearchResults searchResults;
      stk::search::coarse_search(elemBoxes, supersetBoxVec, searchMethod, comm, searchResults, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

template<typename BoxIdentProcType>
void run_volume_to_one_test_with_views(const std::string& meshFileName,
                                       stk::search::SearchMethod searchMethod,
                                       bool enforceSearchResultSymmetry = true)
{

  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  const int pRank = stk::parallel_machine_rank(comm);
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<BoxIdentProcType *, ExecSpace> elemBoxes = createBoundingBoxesForEntities<BoxIdentProcType>(*bulkPtr, stk::topology::ELEM_RANK);

    Kokkos::View<BoxIdentProcType *, ExecSpace> supersetBoxes("Range Boxes", 1);
    supersetBoxes(0) = {elemBoxes[0].box, IdentProc(pRank, pRank)};

    for (unsigned i = 0; i != elemBoxes.extent(0); ++i) {
      stk::search::add_to_box(supersetBoxes(0).box, elemBoxes(i).box);
    }


    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      Kokkos::View<IdentProcIntersection*, ExecSpace> searchResults;
      stk::search::coarse_search(elemBoxes, supersetBoxes, searchMethod, comm, searchResults,
                                 ExecSpace{}, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_KDTREE)
{
  run_volume_to_one_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_MORTON_LBVH)
{
  run_volume_to_one_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_ARBORX)
{
  run_volume_to_one_test<FloatBoxIdentProcVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_with_views_MORTON_LBVH)
{
  run_volume_to_one_test_with_views<FloatBoxIdentProc>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_with_views_ARBORX)
{
  run_volume_to_one_test_with_views<FloatBoxIdentProc>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

template<typename BoxVectorType>
void run_volume_to_one_test_local(const std::string& meshFileName,
                                  stk::search::SearchMethod searchMethod)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    BoxVectorType elemBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK, elemBoxes);

    BoxVectorType supersetBoxVec { elemBoxes[0] };
    auto & [supersetBox, supersetIdent] = supersetBoxVec[0];

    for (const auto & [box, ident] : elemBoxes) {
      stk::search::add_to_box(supersetBox, box);
    }

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      LocalSearchResults searchResults;
      stk::search::local_coarse_search(elemBoxes, supersetBoxVec, searchMethod, searchResults);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

template<typename BoxIdentType>
void run_volume_to_one_test_local_with_views(const std::string& meshFileName,
                                             stk::search::SearchMethod searchMethod)
{

  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<BoxIdentType *, ExecSpace> elemBoxes = createBoundingBoxesForEntities<BoxIdentType>(*bulkPtr, stk::topology::ELEM_RANK);

    Kokkos::View<BoxIdentType *, ExecSpace> supersetBoxes("Range Boxes", 1);
    supersetBoxes(0) = {elemBoxes[0].box, stk::parallel_machine_rank(comm)};

    for (unsigned i = 0; i != elemBoxes.extent(0); ++i) {
      stk::search::add_to_box(supersetBoxes(0).box, elemBoxes(i).box);
    }

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
      stk::search::local_coarse_search(elemBoxes, supersetBoxes, searchMethod, searchResults, ExecSpace{});
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

template<typename BoxIdentType>
void run_one_to_volume_test_local_with_views(const std::string& meshFileName,
                                             stk::search::SearchMethod searchMethod)
{

  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1000;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    Kokkos::View<BoxIdentType *, ExecSpace> elemBoxes = createBoundingBoxesForEntities<BoxIdentType>(*bulkPtr, stk::topology::ELEM_RANK);

    Kokkos::View<BoxIdentType *, ExecSpace> supersetBoxes("Range Boxes", 1);
    supersetBoxes(0) = {elemBoxes[0].box, stk::parallel_machine_rank(comm)};

    for (unsigned i = 0; i != elemBoxes.extent(0); ++i) {
      stk::search::add_to_box(supersetBoxes(0).box, elemBoxes(i).box);
    }

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
      stk::search::local_coarse_search(supersetBoxes, elemBoxes, searchMethod, searchResults, ExecSpace{});
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_local_KDTREE)
{
  run_volume_to_one_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_local_MORTON_LBVH)
{
  run_volume_to_one_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_local_ARBORX)
{
  run_volume_to_one_test_local<FloatBoxIdentVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_local_with_views_MORTON_LBVH)
{
  run_volume_to_one_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

TEST(StkSearch_VolumeToOne, generatedMesh_floatBox_local_with_views_ARBORX)
{
  run_volume_to_one_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_OneToVolume, generatedMesh_floatBox_local_with_views_ARBORX)
{
  run_one_to_volume_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::ARBORX);
}

TEST(StkSearch_OneToVolume, generatedMesh_floatBox_local_with_views_MORTON_LBVH)
{
  run_one_to_volume_test_local_with_views<FloatBoxIdent>("generated:40x80x20|sideset:xXyYzZ", stk::search::MORTON_LBVH);
}

} // namespace

