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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_performance_tests/stk_mesh/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>
namespace stk_perf_many_blocks
{

class ManyBlocksSidesets : public stk::unit_test_util::simple_fields::MeshFixture
{
public:
  ManyBlocksSidesets()
    : timer(get_comm())
  { }

protected:
  std::string get_mesh_spec(unsigned numElemsPerDim)
  {
    std::ostringstream os;
    os<<"generated:"<<numElemsPerDim<<"x"<<numElemsPerDim<<"x"<<numElemsPerDim
      <<"|shell:xXyYzZ|sideset:xXyYzZ";
    return os.str();
  }

  void setup_multi_block_mesh(unsigned numElemsPerDim, unsigned numBlocks)
  {
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    stk::performance_tests::setup_sidesets_between_blocks(get_meta());
    setup_mesh(get_mesh_spec(numElemsPerDim), stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
    stk::performance_tests::fill_sidesets_between_blocks(get_bulk());
  }

  void empty_mod_cycle()
  {
    get_bulk().modification_begin();
    get_bulk().modification_end();
  }

  stk::performance_tests::Timer timer;
};

TEST_F(ManyBlocksSidesets, Timing)
{
  if (get_parallel_size() > 10) return;

  const int NUM_RUNS = 100;
  const int ELEMS_PER_DIM = 100;
  const int NUM_BLOCKS = 100;

  setup_multi_block_mesh(ELEMS_PER_DIM, NUM_BLOCKS);

  timer.start_timing();
  for (int run=0; run<NUM_RUNS; run++) {
    empty_mod_cycle();
  }

  timer.update_timing();
  timer.print_timing(NUM_RUNS);
}

}
