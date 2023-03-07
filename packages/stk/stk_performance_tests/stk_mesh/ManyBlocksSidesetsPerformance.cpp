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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/ExodusTranslator.hpp"
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectUtils.hpp>
namespace stk_perf_many_blocks
{

class ManyBlocksSidesets : public stk::unit_test_util::simple_fields::MeshFixture
{
public:
  ManyBlocksSidesets()
    : batchTimer(get_comm())
  { }

protected:
  std::string get_mesh_spec(unsigned numElemsPerDim)
  {
    std::ostringstream os;
    os<<"generated:"<<numElemsPerDim<<"x"<<numElemsPerDim<<"x"<<numElemsPerDim
      <<"|shell:xXyYzZ|sideset:xXyYzZ";
    return os.str();
  }

  std::string get_mesh_spec_without_sidesets(unsigned numElemsPerDim, unsigned numBlocks)
  {
    std::ostringstream os;
    os<<"generated:"<<numBlocks<<"x"<<numElemsPerDim<<"x"<<numElemsPerDim;
    return os.str();
  }

  void setup_multi_block_mesh(unsigned numElemsPerDim, unsigned numBlocks)
  {
    stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA;
    setup_empty_mesh(auraOption);
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    stk::performance_tests::setup_sidesets_between_blocks(get_meta());
    setup_mesh(get_mesh_spec(numElemsPerDim), auraOption);
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
    stk::performance_tests::fill_sidesets_between_blocks(get_bulk());
  }

  void setup_multi_block_mesh_without_sidesets(unsigned numElemsPerDim, unsigned numBlocks)
  {
    stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA;
    setup_empty_mesh(auraOption);
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    setup_mesh(get_mesh_spec_without_sidesets(numElemsPerDim, numBlocks), auraOption);
    stk::performance_tests::move_elements_to_other_contiguous_blocks(get_bulk(), numBlocks);
  }

  void empty_mod_cycle()
  {
    get_bulk().modification_begin();
    get_bulk().modification_end();
  }

  void output_mesh(stk::mesh::BulkData & bulk, const std::string & fileName)
  {
    std::string writeOutput = stk::unit_test_util::simple_fields::get_option("--output", "off");
    if (writeOutput == "on") {
      stk::io::write_mesh(fileName, bulk);
    }
  }

  void output_mesh(stk::mesh::BulkData & bulk)
  {
    const std::string fileName = std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".g";
    output_mesh(bulk, fileName);
  }

  void print_memory_stats(const std::string& preamble, std::ostream &stream = std::cout)
  {
    size_t maxHwm = 0, minHwm = 0, avgHwm = 0;
    stk::get_memory_high_water_mark_across_processors(get_comm(), maxHwm, minHwm, avgHwm);

    if (get_parallel_rank() == 0) {
      std::ostringstream os;
      const double bytesInMegabyte = 1024*1024;
      os << preamble << " "
         << std::setw(6) << std::fixed << std::setprecision(1)
         << "Max HWM: " <<double(maxHwm)/double(bytesInMegabyte) << "MB"
         <<", Min HWM: "<<double(minHwm)/double(bytesInMegabyte) << "MB"
         <<", Avg HWM: "<<double(avgHwm)/double(bytesInMegabyte) << "MB" <<std::endl;

      stream << os.str();
    }
  }

  void print_stats(const std::string& preamble, std::ostream &stream = std::cout)
  {
    print_memory_stats(preamble, stream);

    unsigned localNumFaces = stk::mesh::count_entities(get_bulk(), stk::topology::FACE_RANK, get_meta().locally_owned_part());
    unsigned globalNumFaces;
    stk::all_reduce_sum(get_comm(), &localNumFaces, &globalNumFaces, 1);

    if (get_parallel_rank() == 0) {
      stream << "Global number of face created: " << globalNumFaces << "\n" << std::endl;
    }
  }

  stk::unit_test_util::BatchTimer batchTimer;
};

TEST_F(ManyBlocksSidesets, timing)
{
  if (get_parallel_size() > 10) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 4000;
  const int ELEMS_PER_DIM = 100;
  const int NUM_BLOCKS = 100;


  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_multi_block_mesh(ELEMS_PER_DIM, NUM_BLOCKS);
    batchTimer.start_batch_timer();
    for (int i = 0; i < NUM_ITERS; i++) {
      empty_mod_cycle();
    }
    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(ManyBlocksSidesets, disconnect_blocks_face_creation)
{
  if (get_parallel_size() > 16) return;

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 1;
  int ELEMS_PER_DIM = stk::unit_test_util::simple_fields::get_command_line_option("--ne", 400);
  int NUM_BLOCKS = stk::unit_test_util::simple_fields::get_command_line_option("--nb", 10);
  bool verbose = stk::unit_test_util::simple_fields::has_option("--v");

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    setup_multi_block_mesh_without_sidesets(ELEMS_PER_DIM, NUM_BLOCKS);

    stk::tools::disconnect_all_blocks(get_bulk());

    stk::performance_tests::setup_sidesets_for_blocks(get_meta());

    stk::mesh::PartVector elemBlocks;
    stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlocks);

    print_memory_stats("Before face creation:");

    batchTimer.start_batch_timer();

    for(stk::mesh::Part* block : elemBlocks) {
      stk::mesh::Part& blockPart = *block;
      unsigned partId = blockPart.id();

      std::string sidesetName = "surface_" + std::to_string(partId);

      stk::mesh::Part* surface = get_meta().get_part(sidesetName);
      EXPECT_TRUE(nullptr != surface) << "Could not find " << sidesetName;
      stk::mesh::create_exposed_block_boundary_sides(get_bulk(), *block, stk::mesh::PartVector{surface});

      if(verbose) {
        print_stats("After face creation for " + sidesetName + ":");
      }
    }

    print_stats("After face creation :");

    batchTimer.stop_batch_timer();
    reset_mesh();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}


}
