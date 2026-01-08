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
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_util/util/FieldDataAllocator.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace {

class FieldStateRotation : public stk::unit_test_util::MeshFixture
{
public:
  FieldStateRotation()
    : batchTimer(get_comm())
  {}

protected:
  void declare_multistate_fields(int numFields, int numStates)
  {
    for (int i = 0; i < numFields; ++i) {
      stk::mesh::FieldBase& newField = get_meta().declare_field<double>(stk::topology::ELEM_RANK,
                                                                        "multistate_"+std::to_string(i), numStates);
      stk::mesh::put_field_on_mesh(newField, get_meta().universal_part(), nullptr);
      m_multistateFields.push_back(&newField);
    }
  }

  void fill_multi_block_mesh(unsigned numElemsPerDim)
  {
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
  }

  void initialize_all_device_data()
  {
    for (stk::mesh::FieldBase* field : m_multistateFields) {
      const int numStates = field->number_of_states();
      for (int state = 0; state < numStates; ++state) {
        stk::mesh::FieldBase* stateField = field->field_state(static_cast<stk::mesh::FieldState>(state));
        stateField->data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Trigger device data creation
      }

    }
  }

  template <typename RotateFunctor>
  void run_host_multiple_block_test(int numIters, int numFields, int numStates, const RotateFunctor& rotateFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;
    const int NUM_BLOCKS = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::setup_multiple_blocks(get_meta(), NUM_BLOCKS);
    declare_multistate_fields(numFields, numStates);
    fill_multi_block_mesh(ELEMS_PER_DIM);

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        rotateFunctor(m_multistateFields);
      }
      batchTimer.stop_batch_timer();
    }

    batchTimer.print_batch_timing(numIters);
  }

  template <typename RotateFunctor>
  void run_host_and_device_multiple_block_test(int numIters, int numFields, int numStates,
                                               const RotateFunctor& rotateFunctor)
  {
    const unsigned NUM_RUNS = 5;
    const int ELEMS_PER_DIM = 100;
    const int NUM_BLOCKS = 100;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::performance_tests::setup_multiple_blocks(get_meta(), NUM_BLOCKS);
    declare_multistate_fields(numFields, numStates);
    fill_multi_block_mesh(ELEMS_PER_DIM);
    initialize_all_device_data();

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        rotateFunctor(m_multistateFields);
      }
      batchTimer.stop_batch_timer();
    }

    batchTimer.print_batch_timing(numIters);
  }
  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::FieldVector m_multistateFields;
};


auto host_rotate_fields = [](const stk::mesh::FieldVector& multistateFields)
{
  stk::mesh::BulkData& bulk = multistateFields.front()->get_mesh();

  for (stk::mesh::FieldBase* field : multistateFields) {
    bulk.update_field_data_states(field);
  }
};

auto host_and_device_through_sync_rotate_fields = [](const stk::mesh::FieldVector& multistateFields)
{
  stk::mesh::BulkData& bulk = multistateFields.front()->get_mesh();
  const bool rotateDeviceData = false;

  for (stk::mesh::FieldBase* field : multistateFields) {
    const int numStates = field->number_of_states();
    for (int state = 0; state < numStates; ++state) {
      stk::mesh::FieldBase* stateField = field->field_state(static_cast<stk::mesh::FieldState>(state));
      stateField->synchronize<stk::mesh::ReadWrite>();  // Modify on host
    }

    bulk.update_field_data_states(field, rotateDeviceData);

    for (int state = 0; state < numStates; ++state) {
      stk::mesh::FieldBase* stateField = field->field_state(static_cast<stk::mesh::FieldState>(state));
      stateField->synchronize<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Sync to device
    }
  }
};

auto host_and_device_through_swap_rotate_fields = [](const stk::mesh::FieldVector& multistateFields)
{
  stk::mesh::BulkData& bulk = multistateFields.front()->get_mesh();
  const bool rotateDeviceData = true;

  for (stk::mesh::FieldBase* field : multistateFields) {
    bulk.update_field_data_states(field, rotateDeviceData);
  }
};

//------------------------------------------------------------------------------
TEST_F(FieldStateRotation, host_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 10000;
  const int numFields = 100;
  const int numStates = 2;

  run_host_multiple_block_test(numIters, numFields, numStates, host_rotate_fields);
}

TEST_F(FieldStateRotation, host_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 10000;
  const int numFields = 100;
  const int numStates = 3;

  run_host_multiple_block_test(numIters, numFields, numStates, host_rotate_fields);
}


//------------------------------------------------------------------------------
TEST_F(FieldStateRotation, host_and_deviceThroughSync_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 10;
  const int numFields = 100;
  const int numStates = 2;

  run_host_and_device_multiple_block_test(numIters, numFields, numStates, host_and_device_through_sync_rotate_fields);
}

TEST_F(FieldStateRotation, host_and_deviceThroughSync_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 10;
  const int numFields = 100;
  const int numStates = 3;

  run_host_and_device_multiple_block_test(numIters, numFields, numStates, host_and_device_through_sync_rotate_fields);
}


//------------------------------------------------------------------------------
TEST_F(FieldStateRotation, host_and_deviceThroughSwap_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 1000;
  const int numFields = 100;
  const int numStates = 2;

  run_host_and_device_multiple_block_test(numIters, numFields, numStates, host_and_device_through_swap_rotate_fields);
}

TEST_F(FieldStateRotation, host_and_deviceThroughSwap_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numIters = 1000;
  const int numFields = 100;
  const int numStates = 3;

  run_host_and_device_multiple_block_test(numIters, numFields, numStates, host_and_device_through_swap_rotate_fields);
}

}
