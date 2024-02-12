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
#include <stk_search/CoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_unit_test_utils/MeshUtilsForBoundingVolumes.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace {

template<typename BoxVectorType>
void run_search_test(const std::string& meshFileName,
                     stk::search::SearchMethod searchMethod,
                     bool communicateResultInfo = true)
{
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 100;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

    stk::io::fill_mesh_with_auto_decomp(meshFileName, *bulkPtr);

    BoxVectorType elemBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::ELEM_RANK, elemBoxes);
    BoxVectorType sideBoxes;
    createBoundingBoxesForEntities(*bulkPtr, stk::topology::FACE_RANK, sideBoxes);

    batchTimer.start_batch_timer();

    for(unsigned i=0; i<NUM_ITERS; ++i) {
      SearchResults searchResults;
      stk::search::coarse_search(elemBoxes, sideBoxes, stk::search::KDTREE, MPI_COMM_WORLD, searchResults);
    }

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkSearch, hex_mesh_KDTREE_float)
{
  run_search_test<FloatBoxVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch, hex_mesh_KDTREE_double)
{
  run_search_test<StkBoxVector>("generated:40x80x20|sideset:xXyYzZ", stk::search::KDTREE);
}

TEST(StkSearch, perf_mesh_KDTREE_float)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") {
    std::cout<<"Skipping test, no mesh specified via '-m' command-line flag."<<std::endl;
    return;
  }

  std::cout<<"running search test with '"<<meshFileName<<"'."<<std::endl;

  run_search_test<FloatBoxVector>(meshFileName, stk::search::KDTREE);
}

TEST(StkSearch, perf_mesh_KDTREE_double)
{
  std::string meshFileName = stk::unit_test_util::get_option("-m", "none-specified");
  if (meshFileName == "none-specified") {
    std::cout<<"Skipping test, no mesh specified via '-m' command-line flag."<<std::endl;
    return;
  }

  std::cout<<"running search test with '"<<meshFileName<<"'."<<std::endl;

  run_search_test<StkBoxVector>(meshFileName, stk::search::KDTREE);
}

} // namespace

