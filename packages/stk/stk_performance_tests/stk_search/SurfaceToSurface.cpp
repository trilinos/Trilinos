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
    SearchResults searchResults;

    batchTimer.start_batch_timer();
    for (int i = 0; i < numIterations; ++i) {
      stk::search::coarse_search(diceBoxes, toolBoxes, searchMethod, comm, searchResults, enforceSearchResultSymmetry);
//      std::cout << "Num intersections = " << searchResults.size() << std::endl;
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
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, ecsl_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  const int numIterations = 4;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, newtonCradle_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 16) GTEST_SKIP();

  const int numIterations = 2000;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, b61NoseCrush_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, coneCrush_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 200;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, jenga_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, jenga_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 32) GTEST_SKIP();

  const int numIterations = 500;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}


TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_KDTREE)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::KDTREE);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_MORTON_LBVH)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::MORTON_LBVH);
}

TEST(StkSearch_SurfaceToSurface, tractorTrailerCrash_floatBox_ARBORX)
{
  std::string boxFileBaseName = stk::unit_test_util::get_option("-m", "none-specified");
  if (boxFileBaseName == "none-specified") GTEST_SKIP();
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 128) GTEST_SKIP();

  const int numIterations = 20;
  run_imported_surface_to_surface_test<FloatBoxVector>(boxFileBaseName, numIterations, stk::search::ARBORX);
}

} // namespace

