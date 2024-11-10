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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_unit_test_utils/MeshUtilsForBoundingVolumes.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include "stk_search/SearchMethod.hpp"

namespace {

std::string get_parallel_filename(const std::string& baseName, MPI_Comm comm)
{
  const int nProcs = stk::parallel_machine_size(comm);
  const int pRank = stk::parallel_machine_rank(comm);

  return baseName + "_" + std::to_string(nProcs) + "." + std::to_string(pRank);
}

template <typename BoxVectorType>
void sequentially_number_all_boxes(BoxVectorType & boxVector, MPI_Comm comm)
{
  const int nProcs = stk::parallel_machine_size(comm);
  const int pRank = stk::parallel_machine_rank(comm);
  std::vector<unsigned long> allNumBoxes(nProcs, 0);
  const unsigned long localNumBoxes = boxVector.size();

  MPI_Allgather(&localNumBoxes, 1, MPI_UNSIGNED_LONG, allNumBoxes.data(), 1, MPI_UNSIGNED_LONG, comm);

  unsigned long id = 0;
  for (int proc = 0; proc < pRank; ++proc) {
    id += allNumBoxes[proc];
  }

  for (auto & [box, ident] : boxVector) {
    ident.set_proc(pRank);
    ident.set_id(id++);
  }
}

template <typename BoxVectorType>
BoxVectorType read_boxes_from_file(const std::string& baseName, MPI_Comm comm)
{
  const std::string fileName = get_parallel_filename(baseName, comm);
  std::ifstream infile(fileName);
  STK_ThrowRequireMsg(infile.good(), "Unable to open file " + fileName);

  BoxVectorType boxVector;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    std::istringstream iss(line);
    boxVector.emplace_back();
    auto & [box, ident] = boxVector.back();
    iss >> box;
  }

  sequentially_number_all_boxes(boxVector, comm);

  return boxVector;
}

template<typename BoxVectorType>
void run_imported_surface_to_surface_test(const std::string& boxFileBaseName,
                                          const int numIterations,
                                          stk::search::SearchMethod searchMethod,
                                          bool enforceSearchResultSymmetry = true)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxes = read_boxes_from_file<BoxVectorType>(boxFileBaseName + ".txt_dice", comm);
  BoxVectorType toolBoxes = read_boxes_from_file<BoxVectorType>(boxFileBaseName + ".txt_tool", comm);

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      SearchResults searchResults;
      stk::search::coarse_search(diceBoxes, toolBoxes, searchMethod, comm, searchResults, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

template<typename BoxIdentProcType>
void run_imported_surface_to_surface_test_with_views(const std::string& boxFileBaseName,
                                                     const int numIterations,
                                                     stk::search::SearchMethod searchMethod,
                                                     bool enforceSearchResultSymmetry = true)
{
  using BoxType = typename BoxIdentProcType::box_type;
  using IdentProcType = typename BoxIdentProcType::ident_proc_type;
  using BoxVectorType = typename std::vector<std::pair<BoxType, IdentProcType>>;
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxesVector = read_boxes_from_file<BoxVectorType>(boxFileBaseName + ".txt_dice", comm);
  BoxVectorType toolBoxesVector = read_boxes_from_file<BoxVectorType>(boxFileBaseName + ".txt_tool", comm);

  Kokkos::View<BoxIdentProcType *, ExecSpace> diceBoxes("diceBoxes", diceBoxesVector.size());
  Kokkos::View<BoxIdentProcType *, ExecSpace> toolBoxes("toolBoxes", toolBoxesVector.size());
  auto diceBoxesHost = Kokkos::create_mirror_view(diceBoxes);
  auto toolBoxesHost = Kokkos::create_mirror_view(toolBoxes);

  for (unsigned i = 0; i < diceBoxesVector.size(); i++) {
    auto boxIdentProcPair = diceBoxesVector[i];
    BoxIdentProcType domainBoxIdentProc{boxIdentProcPair.first, boxIdentProcPair.second};
    diceBoxesHost(i) = domainBoxIdentProc;
  }

  for (unsigned i = 0; i < toolBoxesVector.size(); i++) {
    auto boxIdentProcPair = toolBoxesVector[i];
    BoxIdentProcType rangeBoxIdentProc{boxIdentProcPair.first, boxIdentProcPair.second};
    toolBoxesHost(i) = rangeBoxIdentProc;
  }

  Kokkos::deep_copy(diceBoxes, diceBoxesHost);
  Kokkos::deep_copy(toolBoxes, toolBoxesHost);

  for (unsigned run = 0; run < NUM_RUNS; ++run) {
    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      Kokkos::View<IdentProcIntersection*, ExecSpace> searchResults;
      stk::search::coarse_search(diceBoxes, toolBoxes, searchMethod, comm, searchResults,
                                 ExecSpace{}, enforceSearchResultSymmetry);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, jenga_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxIdentProcVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test_with_views<FloatBoxIdentProc>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

template <typename BoxVectorType>
BoxVectorType read_local_boxes_from_file_and_number(const std::string& baseName)
{
  const std::string fileName = baseName + "_1.0";
  std::ifstream infile(fileName);
  STK_ThrowRequireMsg(infile.good(), "Unable to open file " + fileName);

  BoxVectorType boxVector;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    std::istringstream iss(line);
    boxVector.emplace_back();
    auto & [box, ident] = boxVector.back();
    iss >> box;
  }

  unsigned long id = 0;
  for (auto & boxIdent : boxVector) {
    boxIdent.second = id;
    id++;
  }

  return boxVector;
}

template <typename BoxVectorType>
BoxVectorType read_local_boxes_from_file_and_number(const std::string& baseName, MPI_Comm comm)
{
  const std::string fileName = get_parallel_filename(baseName, comm);
  std::ifstream infile(fileName);
  STK_ThrowRequireMsg(infile.good(), "Unable to open file " + fileName);

  BoxVectorType boxVector;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    std::istringstream iss(line);
    boxVector.emplace_back();
    auto & [box, ident] = boxVector.back();
    iss >> box;
  }

  unsigned long id = 0;
  for (auto & boxIdent : boxVector) {
    boxIdent.second = id;
    id++;
  } 

  return boxVector;
}

template<typename BoxVectorType>
void run_imported_surface_to_surface_test_local(const std::string& boxFileBaseName,
                                                const int numIterations,
                                                stk::search::SearchMethod searchMethod)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxes = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_dice");
  BoxVectorType toolBoxes = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_tool");

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      LocalSearchResults searchResults;
      stk::search::local_coarse_search(diceBoxes, toolBoxes, searchMethod, searchResults);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

template<typename BoxIdentType>
void run_imported_surface_to_surface_test_local_with_views(const std::string& boxFileBaseName,
                                                           const int numIterations,
                                                           stk::search::SearchMethod searchMethod)
{
  using BoxType = typename BoxIdentType::box_type;
  using IdentType = typename BoxIdentType::second_type;
  using BoxVectorType = typename std::vector<std::pair<BoxType, IdentType>>;
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_dice");
  BoxVectorType toolBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_tool");

  Kokkos::View<BoxIdentType *, ExecSpace> diceBoxes("diceBoxes", diceBoxesVector.size());
  Kokkos::View<BoxIdentType *, ExecSpace> toolBoxes("toolBoxes", toolBoxesVector.size());
  auto diceBoxesHost = Kokkos::create_mirror_view(diceBoxes);
  auto toolBoxesHost = Kokkos::create_mirror_view(toolBoxes);

  for (unsigned i = 0; i < diceBoxesVector.size(); i++) {
    auto boxIdentPair = diceBoxesVector[i];
    BoxIdentType domainBoxIdent{boxIdentPair.first, boxIdentPair.second};
    diceBoxesHost(i) = domainBoxIdent;
  }

  for (unsigned i = 0; i < toolBoxesVector.size(); i++) {
    auto boxIdentPair = toolBoxesVector[i];
    BoxIdentType rangeBoxIdent{boxIdentPair.first, boxIdentPair.second};
    toolBoxesHost(i) = rangeBoxIdent;
  }

  Kokkos::deep_copy(diceBoxes, diceBoxesHost);
  Kokkos::deep_copy(toolBoxes, toolBoxesHost);

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
      stk::search::local_coarse_search(diceBoxes, toolBoxes, searchMethod, searchResults, ExecSpace{});
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

template<typename BoxIdentType>
void run_imported_surface_to_surface_test_pll_local_with_views(const std::string& boxFileBaseName,
                                                           const int numIterations,
                                                           stk::search::SearchMethod searchMethod)
{
  using BoxType = typename BoxIdentType::box_type;
  using IdentType = typename BoxIdentType::second_type;
  using BoxVectorType = typename std::vector<std::pair<BoxType, IdentType>>;
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_dice", comm);
  BoxVectorType toolBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_tool", comm);

  Kokkos::View<BoxIdentType *, ExecSpace> diceBoxes("diceBoxes", diceBoxesVector.size());
  Kokkos::View<BoxIdentType *, ExecSpace> toolBoxes("diceBoxes", toolBoxesVector.size());
  auto diceBoxesHost = Kokkos::create_mirror_view(diceBoxes);
  auto toolBoxesHost = Kokkos::create_mirror_view(toolBoxes);

  for (unsigned i = 0; i < diceBoxesVector.size(); i++) {
    auto boxIdentPair = diceBoxesVector[i];
    BoxIdentType domainBoxIdent{boxIdentPair.first, boxIdentPair.second};
    diceBoxesHost(i) = domainBoxIdent;
  }

  for (unsigned i = 0; i < toolBoxesVector.size(); i++) {
    auto boxIdentPair = toolBoxesVector[i];
    BoxIdentType rangeBoxIdent{boxIdentPair.first, boxIdentPair.second};
    toolBoxesHost(i) = rangeBoxIdent;
  }

  Kokkos::deep_copy(diceBoxes, diceBoxesHost);
  Kokkos::deep_copy(toolBoxes, toolBoxesHost);

  std::shared_ptr<stk::search::SearchData> searchData;
  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
      constexpr bool sortResults = false;
      searchData = stk::search::local_coarse_search(diceBoxes, toolBoxes, searchMethod, searchResults, ExecSpace{}, sortResults, searchData);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

using MemSpace = stk::ngp::ExecSpace::memory_space;

template<typename BoxIdentType>
void run_imported_surface_to_surface_test_local_with_views_rawArborX(const std::string& boxFileBaseName,
                                                           const int numIterations,
                                                           stk::search::SearchMethod searchMethod)
{
  using BoxType = typename BoxIdentType::box_type;
  using IdentType = typename BoxIdentType::second_type;
  using BoxVectorType = typename std::vector<std::pair<BoxType, IdentType>>;
  using ExecSpace = Kokkos::DefaultExecutionSpace;

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  BoxVectorType diceBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_dice", comm);
  BoxVectorType toolBoxesVector = read_local_boxes_from_file_and_number<BoxVectorType>(boxFileBaseName + ".txt_tool", comm);

  Kokkos::View<ArborX::Box *, MemSpace> diceBoxes("diceBoxes", diceBoxesVector.size());
  Kokkos::View<ArborX::Box *, MemSpace> toolBoxes("diceBoxes", toolBoxesVector.size());
  auto diceBoxesHost = Kokkos::create_mirror_view(diceBoxes);
  auto toolBoxesHost = Kokkos::create_mirror_view(toolBoxes);

  for (unsigned i = 0; i < diceBoxesVector.size(); i++) {
    auto boxIdentPair = diceBoxesVector[i];
    auto stkMinCorner = boxIdentPair.first.min_corner();
    auto stkMaxCorner = boxIdentPair.first.max_corner();
    ArborX::Point min_point(stkMinCorner[0], stkMinCorner[1], stkMinCorner[2]);
    ArborX::Point max_point(stkMaxCorner[0], stkMaxCorner[1], stkMaxCorner[2]);
    diceBoxesHost(i) = {min_point, max_point};
  }

  for (unsigned i = 0; i < toolBoxesVector.size(); i++) {
    auto boxIdentPair = toolBoxesVector[i];
    auto stkMinCorner = boxIdentPair.first.min_corner();
    auto stkMaxCorner = boxIdentPair.first.max_corner();
    ArborX::Point min_point(stkMinCorner[0], stkMinCorner[1], stkMinCorner[2]);
    ArborX::Point max_point(stkMaxCorner[0], stkMaxCorner[1], stkMaxCorner[2]);
    toolBoxesHost(i) = {min_point, max_point};
  }

  Kokkos::deep_copy(diceBoxes, diceBoxesHost);
  Kokkos::deep_copy(toolBoxes, toolBoxesHost);

  for (unsigned run = 0; run < NUM_RUNS; ++run) {
    batchTimer.start_batch_timer();
    ExecSpace execSpace{};
    for (int i = 0; i < numIterations; ++i) {
      Kokkos::Profiling::pushRegion("Raw ArborX");
      Kokkos::View<int *, MemSpace> indices("ArborX::indices", 0);
      Kokkos::View<int *, MemSpace> offsets("ArborX::offsets", 0);

      ArborX::BVH<MemSpace> bvh{execSpace, toolBoxes};
      const int numQueries = diceBoxes.extent(0);
      Kokkos::View<ArborX::Intersects<ArborX::Box> *, MemSpace> queries(Kokkos::ViewAllocateWithoutInitializing("queries"), numQueries);

      Kokkos::parallel_for("setup_queries", Kokkos::RangePolicy<ExecSpace>(0, numQueries),
                           KOKKOS_LAMBDA(int i) { queries(i) = ArborX::intersects(diceBoxes(i)); });
      Kokkos::fence();
      bvh.query(execSpace, queries, indices, offsets);
      Kokkos::Profiling::popRegion();
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numIterations);
}

TEST(StkSearch_SurfaceToSurface, a001_intent_strong_link_floatBox_local_with_views_rawARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 40) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test_local_with_views_rawArborX<FloatBoxIdent>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, a001_intent_strong_link_floatBox_local_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 40) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test_pll_local_with_views<FloatBoxIdent>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, a001_intent_strong_link_floatBox_local_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 40) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test_pll_local_with_views<FloatBoxIdent>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_local_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_local<FloatBoxIdentVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_local_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_local<FloatBoxIdentVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_local_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_local<FloatBoxIdentVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_local_with_views_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_local_with_views<FloatBoxIdent>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_local_with_views_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test_local_with_views<FloatBoxIdent>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

} // namespace

