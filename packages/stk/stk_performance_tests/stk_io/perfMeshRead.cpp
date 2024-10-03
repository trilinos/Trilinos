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
#include <stk_util/stk_config.h>
#include <stk_unit_test_utils/getOption.h>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace stk_perf_io_mesh_read
{

TEST(StkIo, meshRead_hex_noAura)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 10;
  int ELEMS_PER_DIM = stk::unit_test_util::get_command_line_option("--ne", 80);
  std::string elems = std::to_string(ELEMS_PER_DIM);
  std::string meshSpec = "generated:"+elems+"x"+elems+"x"+elems;

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    std::cout << "Using mesh-spec: " << meshSpec << std::endl;
  }

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {

    stk::parallel_machine_barrier(MPI_COMM_WORLD);

    batchTimer.start_batch_timer();

    for(unsigned i=0; i<NUM_ITERS; ++i) {
      std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                          .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                          .create();
      stk::io::fill_mesh(meshSpec, *bulkPtr);
    }

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(StkIo, meshRead_hex_shells_sidesets_aura)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 10;
  int ELEMS_PER_DIM = stk::unit_test_util::get_command_line_option("--ne", 80);
  std::string elems = std::to_string(ELEMS_PER_DIM);
  std::string meshSpec = "generated:"+elems+"x"+elems+"x"+elems+"|shell:xyzXYZ|sideset:xyzXYZ";

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    std::cout << "Using mesh-spec: " << meshSpec << std::endl;
  }

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {

    stk::parallel_machine_barrier(MPI_COMM_WORLD);

    batchTimer.start_batch_timer();

    for(unsigned i=0; i<NUM_ITERS; ++i) {
      std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                          .create();
      stk::io::fill_mesh(meshSpec, *bulkPtr);
    }

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
