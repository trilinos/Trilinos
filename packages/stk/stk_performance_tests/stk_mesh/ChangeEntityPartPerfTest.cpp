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

#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace
{

class ChangePartsTest : public stk::unit_test_util::MeshFixture
{
public:
  ChangePartsTest() : stk::unit_test_util::MeshFixture(),
    batchTimer(get_comm()),
    elementsOnBlock1(true)
  {
  }

  void setup_mesh_with_many_blocks_many_elements_in_one_block()
  {
    std::string meshDesc = "generated:" + std::to_string(numElemPerDim) + "x"
                                        + std::to_string(numElemPerDim) + "x"
                                        + std::to_string(numElemPerDim);
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    setup_multiple_blocks();

    stk::io::fill_mesh(meshDesc, get_bulk());

    setup_block_part_vectors();
  }

  void setup_block_part_vectors()
  {
    stk::mesh::Part* block1Part = get_meta().get_part("block_1");
    stk::mesh::Part* block2Part = &get_meta().declare_part("block_2", stk::topology::ELEM_RANK);

    block1PartVector.push_back(block1Part);
    block2PartVector.push_back(block2Part);
  }

  void setup_multiple_blocks()
  {
    std::string block1Name("block_1");
    stk::mesh::Part& block1Part = get_meta().declare_part_with_topology(block1Name, stk::topology::HEX_8);
    get_meta().set_part_id(block1Part, 1);
    parts.push_back(&block1Part);

    for (unsigned i = 1; i < numBlocks; i++) {
      unsigned partId = i+1;
      std::string blockName = "block_" + std::to_string(partId);
      stk::mesh::Part& part = get_meta().declare_part_with_topology(blockName, stk::topology::HEX_8);
      get_meta().set_part_id(part, partId);
      parts.push_back(&part);
    }
  }

  void elements_from_block1_to_other_blocks()
  {
    if(numBlocks < std::pow(numElemPerDim, 3)) { return; }

    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, entities);
    stk::mesh::Part* removePart = parts[0];

    for(unsigned i = 1; i < entities.size(); i++) {
      stk::mesh::Part* addPart = parts[i];

      get_bulk().modification_begin();
      get_bulk().change_entity_parts(entities[i], stk::mesh::PartVector{addPart}, stk::mesh::PartVector{removePart});
      get_bulk().modification_end();
    }
  }

  void elements_back_to_block1()
  {
    if(numBlocks < std::pow(numElemPerDim, 3)) { return; }

    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, entities);

    for(unsigned i = 1; i < entities.size(); i++) {
      stk::mesh::Part* addPart = parts[0];
      stk::mesh::Part* removePart = parts[i];

      get_bulk().modification_begin();
      get_bulk().change_entity_parts(entities[i], stk::mesh::PartVector{addPart}, stk::mesh::PartVector{removePart});
      get_bulk().modification_end();
    }
  }

  void elements_back_to_block1_using_selector()
  {
    if(numBlocks < std::pow(numElemPerDim, 3)) { return; }

    stk::mesh::Selector selector = stk::mesh::selectUnion(parts);
    stk::mesh::Part* addPart = parts[0];

    get_bulk().batch_change_entity_parts(selector, stk::topology::ELEM_RANK, stk::mesh::PartVector{addPart}, stk::mesh::PartVector{});
  }

  void prepare_move_elements_between_block1_and_block2(stk::mesh::PartVector& addParts, stk::mesh::PartVector& removeParts)
  {
    if(elementsOnBlock1) {
      addParts = block2PartVector;
      removeParts = block1PartVector;
    } else {
      addParts = block1PartVector;
      removeParts = block2PartVector;
    }
    elementsOnBlock1 = !elementsOnBlock1;
  }

protected:
  unsigned numBlocks;
  unsigned numElemPerDim;
  stk::mesh::PartVector parts;
  stk::mesh::PartVector block1PartVector;
  stk::mesh::PartVector block2PartVector;
  stk::unit_test_util::BatchTimer batchTimer;

private:
  bool elementsOnBlock1;
};

TEST_F(ChangePartsTest, changeEntityPartsUsingEntityVectorSimplePerfTest)
{
  if(get_parallel_size() > 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  numElemPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  numBlocks = 1;

  batchTimer.initialize_batch_timer();

  setup_mesh_with_many_blocks_many_elements_in_one_block();

  stk::mesh::PartVector addParts;
  stk::mesh::PartVector removeParts;
  stk::mesh::EntityVector elements;
  get_bulk().get_entities(stk::topology::ELEM_RANK, stk::mesh::Selector(*block1PartVector[0]), elements);

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for(unsigned i = 0; i < NUM_ITERS; i++) {
      prepare_move_elements_between_block1_and_block2(addParts, removeParts);
      get_bulk().batch_change_entity_parts(elements, addParts, removeParts);
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(ChangePartsTest, changeEntityPartsUsingSelectorSimplePerfTest)
{
  if(get_parallel_size() > 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 200);
  numElemPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  numBlocks = 1;

  batchTimer.initialize_batch_timer();

  setup_mesh_with_many_blocks_many_elements_in_one_block();

  stk::mesh::PartVector addParts;
  stk::mesh::PartVector removeParts;

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for(unsigned i = 0; i < NUM_ITERS; i++) {
      prepare_move_elements_between_block1_and_block2(addParts, removeParts);
      stk::mesh::Selector elemSelector(*removeParts[0]);
      get_bulk().batch_change_entity_parts(elemSelector, stk::topology::ELEM_RANK, addParts, removeParts);
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(ChangePartsTest, cacheRemovalImpactChangeEntityPartsWithEntityVector)
{
  if(get_parallel_size() > 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 500000);
  numElemPerDim = stk::unit_test_util::get_command_line_option("-e", 80);
  numBlocks = stk::unit_test_util::get_command_line_option("-b", 125);

  batchTimer.initialize_batch_timer();

  setup_mesh_with_many_blocks_many_elements_in_one_block();

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for(unsigned i = 0; i < NUM_ITERS; i++) {
      elements_from_block1_to_other_blocks();
      elements_back_to_block1();

      for(unsigned block = 0; block < numBlocks; block++) {
        get_bulk().get_buckets(stk::topology::ELEM_RANK, stk::mesh::Selector(*parts[block]));
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(ChangePartsTest, cacheRemovalImpactChangeEntityPartsWithSelector)
{
  if(get_parallel_size() > 1) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 500000);
  numElemPerDim = stk::unit_test_util::get_command_line_option("-e", 80);
  numBlocks = stk::unit_test_util::get_command_line_option("-b", 125);

  batchTimer.initialize_batch_timer();

  setup_mesh_with_many_blocks_many_elements_in_one_block();

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for(unsigned i = 0; i < NUM_ITERS; i++) {
      elements_from_block1_to_other_blocks();
      elements_back_to_block1_using_selector();

      for(unsigned block = 0; block < numBlocks; block++) {
        get_bulk().get_buckets(stk::topology::ELEM_RANK, stk::mesh::Selector(*parts[block]));
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
