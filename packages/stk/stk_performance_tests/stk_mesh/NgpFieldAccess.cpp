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
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_util/stk_config.h>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace ngp_field_perf_test
{

class NgpFieldAccess : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldAccess()
    : batchTimer(get_comm())
  { }

protected:
  void declare_centroid_field()
  {
    centroid = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*centroid, get_meta().universal_part(), 3, nullptr);

    hostCentroid = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "hostCentroid");
    stk::mesh::put_field_on_mesh(*hostCentroid, get_meta().universal_part(), 3, nullptr);
  }

  void declare_centroid_partial_mesh(unsigned numBlocks)
  {
    centroid = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "centroid");
    for(unsigned i = 1; i <= numBlocks; i++) {
      const std::string partName = "block_" + std::to_string(i);
      stk::mesh::Part& part = get_meta().declare_part(partName, stk::topology::ELEM_RANK);
      stk::mesh::put_field_on_mesh(*centroid, part, 3, nullptr);
    }
  }

  void setup_multi_block_mesh(unsigned numElemsPerDim, unsigned numBlocks)
  {
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
  }

  void verify_averaged_centroids_are_center_of_mesh(int elemsPerDim, stk::mesh::Field<double>& centroidField)
  {
    std::vector<double> average = stk::performance_tests::get_centroid_average_from_host(get_bulk(), centroidField, stk::mesh::Selector(get_meta().universal_part()));
    double meshCenter = elemsPerDim/2.0;
    for(size_t dim = 0; dim < 3; dim++) {
      EXPECT_DOUBLE_EQ(meshCenter, average[dim]);
    }
  }

  void verify_averaged_centroids_are_center_of_mesh(int elemsPerDim, const stk::mesh::Selector& selector)
  {
    std::vector<double> hostAverage = stk::performance_tests::get_centroid_average_from_host(get_bulk(), *centroid, selector);
    std::vector<double> deviceAverage = stk::performance_tests::get_centroid_average_from_device(get_bulk(), *centroid, selector);
    for(size_t dim = 0; dim < 3; dim++) {
      EXPECT_DOUBLE_EQ(hostAverage[dim], deviceAverage[dim]);
    }
  }

  void compare_and_verify_average_centroids(int elemsPerDim, const stk::mesh::Selector& selector)
  {
    std::vector<double> centroidAverage = stk::performance_tests::get_centroid_average_from_host(get_bulk(), *centroid, selector);
    std::vector<double> hostCentroidAverage = stk::performance_tests::get_centroid_average_from_host(get_bulk(), *hostCentroid, selector);
    for(size_t dim = 0; dim < 3; dim++) {
      EXPECT_DOUBLE_EQ(centroidAverage[dim], hostCentroidAverage[dim]);
    }
  }

  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::Field<double> *centroid;
  stk::mesh::Field<double> *hostCentroid;
};

TEST_F(NgpFieldAccess, Centroid)
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 100;
  const int ELEMS_PER_DIM = 120;

  batchTimer.initialize_batch_timer();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_centroid_field();
  stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for (int i = 0; i <NUM_ITERS; i++) {
      stk::performance_tests::calculate_centroid_using_coord_field<stk::mesh::NgpField<double>>(get_bulk(), *centroid);
    }
    verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM, *centroid);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAccess, HostCentroid)
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 100;
  const int ELEMS_PER_DIM = 120;

  batchTimer.initialize_batch_timer();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_centroid_field();
  stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), get_bulk());

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (int i = 0; i < NUM_ITERS; i++) {
      stk::performance_tests::calculate_centroid_using_host_coord_fields(get_bulk(), *hostCentroid);
    }
    verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM, *hostCentroid);

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAccess, CentroidMultiBlock)
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 100;
  const int ELEMS_PER_DIM = 100;
  const int NUM_BLOCKS = 100;

  batchTimer.initialize_batch_timer();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_centroid_field();
  setup_multi_block_mesh(ELEMS_PER_DIM, NUM_BLOCKS);

  stk::mesh::PartVector elemBlockParts;
  stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlockParts);

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    for (int i = 0; i < NUM_ITERS; i++) {
      for (const stk::mesh::Part* blockPart : elemBlockParts) {
        stk::performance_tests::calculate_centroid_using_coord_field<stk::mesh::NgpField<double>>(
              get_bulk(), *blockPart, *centroid);
      }
    }

    verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM, get_meta().universal_part());
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAccess, CentroidPartialBlock)
{
  if (get_parallel_size() != 1) return;

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 250;
  const int ELEMS_PER_DIM = 100;
  const int NUM_BLOCKS = 100;
  int BLOCKS = stk::unit_test_util::get_command_line_option<int>("-n", 50);
  BLOCKS = std::max(BLOCKS, 1);
  BLOCKS = std::min(BLOCKS, NUM_BLOCKS);

  batchTimer.initialize_batch_timer();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_centroid_partial_mesh(BLOCKS);
  setup_multi_block_mesh(ELEMS_PER_DIM, NUM_BLOCKS);

  stk::mesh::PartVector elemBlockParts;
  stk::mesh::fill_element_block_parts(get_meta(), stk::topology::HEX_8, elemBlockParts);

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
  
    stk::mesh::get_updated_ngp_mesh(get_bulk());
    for (int i = 0; i < NUM_ITERS; i++) {
      for (const stk::mesh::Part* blockPart : elemBlockParts) {
        stk::performance_tests::calculate_centroid_using_coord_field<stk::mesh::NgpField<double>>(
              get_bulk(), *blockPart, *centroid);
        verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM, *blockPart);
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
