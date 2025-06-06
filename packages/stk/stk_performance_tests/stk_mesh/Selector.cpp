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
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include "stk_unit_test_utils/stk_mesh_fixtures/SelectorFixture.hpp"
#include <stdexcept>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace {

using stk::mesh::fixtures::VariableSelectorFixture;

}

TEST(selector_timings, selector_timings)
{
  // A scaling test for selector and bucket operations

  // If we are running with STL in debug mode we shrink the problem
  // down in order to keep things running in a reasonable amount of
  // time.
#ifdef _GLIBCXX_DEBUG
  size_t N = 500;
#else
  size_t N = 50000;
#endif

  const unsigned NUM_RUNS = 5;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    VariableSelectorFixture fix(N);

    size_t total_buckets_grabbed = 0;

    batchTimer.start_batch_timer();

    for (size_t n = 1 ; n<N; n*=2) {
      // Selector creation
      stk::mesh::Selector selectUnion;
      for (size_t part_i = 0 ; part_i<n ; ++part_i) {
        selectUnion |= *fix.m_declared_part_vector[part_i];
      }

      // Selector usage:
      stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK;
      stk::mesh::BucketVector const& buckets_out =  fix.m_BulkData.get_buckets(entity_rank, selectUnion);
      total_buckets_grabbed += buckets_out.size();
    }

    std::cout << "total_buckets_grabbed: " << total_buckets_grabbed << std::endl;

    batchTimer.stop_batch_timer();
  }
  const int numIterations = 1;
  batchTimer.print_batch_timing(numIterations);
}

void create_elem_block_parts(int numParts, stk::mesh::MetaData& meta)
{
  const int firstBlockAfterBlock1 = 2;
  for(int i=0; i<numParts; ++i) {
    stk::mesh::Part& part = meta.declare_part(std::string("elemBlock")+std::to_string(i+firstBlockAfterBlock1), stk::topology::ELEM_RANK);
    meta.set_part_id(part, i+firstBlockAfterBlock1);
  }
}

void move_elems_to_each_part(stk::mesh::BulkData& bulk, int elemsPerPart)
{
  stk::mesh::PartVector blocks = stk::mesh::ExodusTranslator(bulk).get_element_block_parts();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Part* block1 = meta.get_part("block_1");
  STK_ThrowRequire(block1 != nullptr);
  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, meta.locally_owned_part(), elems);

  bulk.modification_begin();
  
  unsigned elemCounter = 0;
  unsigned blockCounter = 0;
  stk::mesh::PartVector addParts = {blocks[blockCounter++]};
  stk::mesh::PartVector rmParts = {block1};

  while(elemCounter < elems.size()) {
    if (addParts[0] != block1) {
      for(int e=0; e<elemsPerPart; ++e) {
        if (elemCounter >= elems.size()) {
          break;
        }
        bulk.change_entity_parts(elems[elemCounter++], addParts, rmParts);
      }
    }
    if (blockCounter >= blocks.size()) {
      break;
    }
    addParts[0] = blocks[blockCounter++];
  }

  bulk.modification_end();
}

TEST(selector_timings, get_buckets)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();

  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  const int numParts = 5000;
  create_elem_block_parts(numParts, meta);

  stk::io::fill_mesh("generated:100x100x100", *bulkPtr);
  const int numElems = 100*100*100;
  const int elemsPerPart = numElems/numParts;
  move_elems_to_each_part(*bulkPtr, elemsPerPart);

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  batchTimer.start_batch_timer();

  stk::mesh::PartVector blocks = stk::mesh::ExodusTranslator(*bulkPtr).get_element_block_parts();
  for(stk::mesh::Part* partPtr : blocks) {
    stk::mesh::Selector selector(*partPtr);
    selector.get_buckets(stk::topology::NODE_RANK);
  }

  batchTimer.stop_batch_timer();

  const int numIterations = 1;
  batchTimer.print_batch_timing(numIterations);
}

TEST(selector_timings, selectPart)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();

  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  const int numParts = 5000;
  create_elem_block_parts(numParts, meta);

  stk::io::fill_mesh("generated:100x100x100", *bulkPtr);
  const int numElems = 100*100*100;
  const int elemsPerPart = numElems/numParts;
  move_elems_to_each_part(*bulkPtr, elemsPerPart);

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  stk::mesh::PartVector blocks = stk::mesh::ExodusTranslator(*bulkPtr).get_element_block_parts();
  stk::mesh::Selector selector = stk::mesh::selectUnion(blocks);
  stk::mesh::PartVector blockSubsets;
  unsigned counter = 1;
  for(stk::mesh::Part* blockPtr : blocks) {
    stk::mesh::Part& blkSubset = meta.declare_part("sub"+std::to_string(counter), stk::topology::ELEM_RANK);
    ++counter;
    blockSubsets.push_back(&blkSubset);
    meta.declare_part_subset(*blockPtr, blkSubset);
  }

  batchTimer.start_batch_timer();

  for(int n=0; n<100; ++n) {
    for(stk::mesh::Part* partPtr : blockSubsets) {
      STK_ThrowRequire((selector(*partPtr)));
    }
  }

  batchTimer.stop_batch_timer();

  const int numIterations = 1;
  batchTimer.print_batch_timing(numIterations);
}

TEST(Verify, selectorAlgorithmicComplexity)
{
  //
  //  Verify expected algorithmic complexity is obtained by selector classes.
  //  Only test serial, in parallel this gets a bit ill posed as some procs have empty bucket selections
  //  and the complexity is not meaningful.
  //

  stk::mesh::fixtures::SelectorFixture fix;
  fix.m_meta_data.commit();
  fix.m_bulk_data.modification_begin();
  fix.generate_mesh();


  if (fix.get_BulkData().parallel_size() > 1) return;

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Part & partC = fix.m_partC;
  stk::mesh::Part & partD = fix.m_partD;

  stk::mesh::Selector selectABC = partA & partB & partC;
  stk::mesh::Selector selectA   = partA;
  stk::mesh::Selector selectB   = partB;
  stk::mesh::Selector selectC   = partC;
  stk::mesh::Selector selectD   = partD;

  unsigned numCopiesBase = 1000000;
  std::vector<unsigned> trialLen;

  trialLen.push_back(128);
  trialLen.push_back(256);
  trialLen.push_back(512);
  trialLen.push_back(1024);
  trialLen.push_back(2048);
  trialLen.push_back(4096);

  std::vector<double> elapsedTime(trialLen.size());

  //
  //  Test complexity of construction of large binary selector as would be made of selectUnion of a large number of parts
  //
  for(unsigned j=0; j<trialLen.size(); ++j) {

    unsigned numCopies = (numCopiesBase/2)/trialLen[j];

    double startTime = stk::cpu_time();
    for(unsigned i=0; i<numCopies; ++i) {
      stk::mesh::Selector selectBig ;
      for(unsigned k=0; k< trialLen[j] ; ++k) {
        selectBig |= selectC & selectD;
      }
    }
    elapsedTime[j] = (stk::cpu_time()-startTime)/(double)numCopies;
    std::cout<<"Selector Constructor Complexity Check, Length: "<<trialLen[j]
               <<" time: "<<elapsedTime[j]<<" Expecting linear complexity"<<std::endl;
  }
  //  Expected factor = N, which implies the |= algorithm is O(N) where N is the number of times or is called
  //  Due to cache effects and random runtime variability though keep a pretty loose tolerane around this expectation
  double expectedFactor = ((double)trialLen.back())/trialLen[0];
  double factor         = elapsedTime.back()/elapsedTime[0];
  double relative_factor = std::abs(expectedFactor-factor)/expectedFactor;

  std::cout << "relative_factor= " << relative_factor << std::endl;

  double gold_relative_factor = 0.40;
  EXPECT_TRUE( relative_factor < gold_relative_factor );

  std::cout<<"  Speedup factors: "<<expectedFactor<<" vs. "<<factor<<std::endl;

  //
  //  Test complexity of traversal of the large selector.  Some buckets will be found immediately and some will require traversal of
  //  the entire selector to exclude thus this operation should be O(N) in selector size.
  const stk::mesh::BucketVector& buckets = fix.get_BulkData().get_buckets(stk::topology::NODE_RANK, fix.get_MetaData().locally_owned_part());

  for(unsigned j=0; j<trialLen.size(); ++j) {
    unsigned count=0;

    unsigned numCopies = numCopiesBase/trialLen[j];

    stk::mesh::Selector selectCD = selectC | selectD;
    stk::mesh::Selector selectBig = selectCD;
    for(int i=0; i< ((int)trialLen[j]-1) ; ++i) {
      selectBig &= selectCD;
    }
    double startTime = stk::cpu_time();
    for(unsigned i=0; i<numCopies; ++i) {
      for(unsigned ibucket=0, n=buckets.size(); ibucket<n; ++ibucket) {
        if(selectBig(*(buckets[ibucket]))) {
          count++;
        }
      }
    }

    EXPECT_TRUE(count == 2 * numCopies);

    elapsedTime[j] = (stk::cpu_time()-startTime)/(double)numCopies;
    std::cout<<"Selector Traversal Complexity Check, Length: "<<trialLen[j]
               <<" time: "<<elapsedTime[j]<<" Expecting linear complexity"<<std::endl;
  }

  expectedFactor = ((double)trialLen.back())/trialLen[0];
  factor         = elapsedTime.back()/elapsedTime[0];

  std::cout<<"  Speedup factors: "<<expectedFactor<<" vs. "<<factor<<std::endl;

  relative_factor = std::abs(expectedFactor-factor)/expectedFactor;
  std::cout <<" relative_factor= " << relative_factor << "\n";
  EXPECT_TRUE( relative_factor < 0.60);
}
