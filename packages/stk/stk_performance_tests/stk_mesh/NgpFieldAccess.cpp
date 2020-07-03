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
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/stk_config.h>
#include <stk_performance_tests/stk_mesh/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace ngp_field_perf_test
{

class NgpFieldAccess : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldAccess()
    : timer(get_comm())
  { }

protected:
  void declare_centroid_field()
  {
    centroid = &get_meta().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field_on_mesh(*centroid, get_meta().universal_part(), 3,
                                 (stk::mesh::FieldTraits<stk::mesh::Field<double, stk::mesh::Cartesian3d> >::data_type*) nullptr);
  }

  void setup_multi_block_mesh(unsigned numElemsPerDim, unsigned numBlocks)
  {
    stk::performance_tests::setup_multiple_blocks(get_meta(), numBlocks);
    setup_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim, numBlocks);
  }

  stk::mesh::Selector block_selector(unsigned blockId)
  {
    std::string blockName = "block_" + std::to_string(blockId);
    stk::mesh::Selector selector(*get_meta().get_part(blockName));
    return selector;
  }

  void verify_averaged_centroids_are_center_of_mesh(int elemsPerDim)
  {
    std::vector<double> average = stk::performance_tests::get_centroid_average(get_bulk(), *centroid);
    double meshCenter = elemsPerDim/2.0;
    for(size_t dim = 0; dim < 3; dim++) {
      EXPECT_EQ(meshCenter, average[dim]);
    }
  }

  stk::performance_tests::Timer timer;
  stk::mesh::Field<double, stk::mesh::Cartesian3d> *centroid;
};

TEST_F(NgpFieldAccess, Centroid)
{
  if (get_parallel_size() != 1) return;

  const int NUM_RUNS = 400;
  const int ELEMS_PER_DIM = 100;

  declare_centroid_field();
  setup_mesh(stk::unit_test_util::get_mesh_spec(ELEMS_PER_DIM), stk::mesh::BulkData::NO_AUTO_AURA);

  timer.start_timing();
  for (int run=0; run<NUM_RUNS; run++) {
    stk::performance_tests::calculate_centroid_using_coord_field<stk::mesh::NgpField<double>>(get_bulk(), *centroid);
    verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM);
  }
  timer.update_timing();
  timer.print_timing(NUM_RUNS);
}

TEST_F(NgpFieldAccess, CentroidMultiBlock)
{
  if (get_parallel_size() != 1) return;

  const int NUM_RUNS = 5;
  const int ELEMS_PER_DIM = 100;
  const int NUM_BLOCKS = 100;

  declare_centroid_field();
  setup_multi_block_mesh(ELEMS_PER_DIM, NUM_BLOCKS);

  timer.start_timing();
  for (int run=0; run<NUM_RUNS; run++) {
    for (int blockId=1; blockId<=NUM_BLOCKS; blockId++) {
      stk::performance_tests::calculate_centroid_using_coord_field<stk::mesh::NgpField<double>>(
            get_bulk(), block_selector(blockId), *centroid);
    }

    verify_averaged_centroids_are_center_of_mesh(ELEMS_PER_DIM);
  }
  timer.update_timing();
  timer.print_timing(NUM_RUNS);
}

}
