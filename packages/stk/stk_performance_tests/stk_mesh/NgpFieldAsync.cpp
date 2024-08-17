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
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpExecutionSpace.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>
#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <string>
#include <time.h>
#include <chrono>
#include <thread>

#define SPEEDUP_DELTA 1.0

class NgpFieldAsyncTest : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldAsyncTest()
  : stk::unit_test_util::MeshFixture(),
    m_numBlocks(1),
    m_numElemsPerDim(100),
    m_numElements(std::pow(m_numElemsPerDim, 3)),
    m_numComponents(3),
    m_increment(10),
    m_defaultLaunchBlockingEnvVarSet(false),
    m_defaultLaunchBlockingEnvVarValue(0)
  {
    set_launch_blocking_env_var();
  }

  ~NgpFieldAsyncTest()
  {
    revert_launch_blocking_env_var();
  }

  void setup_simple_mesh_with_fields(unsigned numElemsPerDim)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    m_numElemsPerDim = numElemsPerDim;
    setup_fields();
    setup_mesh_with_many_blocks_many_elements();
  }

  void setup_multi_block_mesh_with_field_per_block(unsigned numElemsPerDim, unsigned numBlocks, unsigned numFields)
  {
    m_numElemsPerDim = numElemsPerDim;
    m_numBlocks = numBlocks;

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    setup_fields_on_all_blocks(numFields);
    setup_mesh_with_many_blocks_many_elements();
  }

  void setup_fields()
  {
    std::vector<int> init;
    setup_field_component_data(init);

    auto field1 = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "intField1", 1);
    auto field2 = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "intField2", 1);
    auto field3 = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "intField3", 1);

    stk::mesh::put_field_on_mesh(*field1, get_meta().universal_part(), m_numComponents, init.data());
    stk::mesh::put_field_on_mesh(*field2, get_meta().universal_part(), m_numComponents, init.data());
    stk::mesh::put_field_on_mesh(*field3, get_meta().universal_part(), m_numComponents, init.data());
  }

  void setup_fields_on_all_blocks(unsigned numFields)
  {
    unsigned numStates = 1;
    std::vector<int> init;
    setup_field_component_data(init);

    for(unsigned i = 1; i <= m_numBlocks; i++) {
      std::string blockName = "block_" + std::to_string(i);

      stk::mesh::Part& part = get_meta().declare_part_with_topology(blockName, stk::topology::HEX_8);
      get_meta().set_part_id(part, i);
      EXPECT_NE(&part, nullptr);

      for(unsigned j = 1; j <= numFields; j++) {
        std::string fieldName = "intField" + std::to_string(j);
        stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::ELEM_RANK, fieldName, numStates);
        stk::mesh::put_field_on_mesh(field, part, m_numComponents, init.data());
      }
    }
  }

  void setup_mesh_with_many_blocks_many_elements(stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA)
  {
    std::string meshDesc = "generated:" + std::to_string(m_numElemsPerDim) + "x"
                                        + std::to_string(m_numElemsPerDim) + "x"
                                        + std::to_string(m_numElemsPerDim);
    stk::performance_tests::setup_multiple_blocks(get_meta(), m_numBlocks);
    setup_mesh(meshDesc, auraOption);
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), m_numElemsPerDim);
  }

  void pass_time_on_device(const stk::ngp::ExecSpace& space, unsigned iterationSpent = 100)
  {
    typedef typename stk::ngp::DeviceTeamPolicy::member_type TeamHandleType;
    const auto& teamPolicy = stk::ngp::DeviceTeamPolicy(space, 10, Kokkos::AUTO);

    Kokkos::parallel_for("run_with_team_policy", teamPolicy,
                        KOKKOS_LAMBDA(const TeamHandleType & team) {
                          for(unsigned j = 0; j < iterationSpent; ++j) {
                            clock_t start = clock();
                            clock_t now;
                            for (;;) {
                              now = clock();
                              clock_t cycles = now > start ? now - start : now + (0xffffffff - start);
                              if (cycles >= 1e6) {
                                break;
                              }
                            }
                          }
                        }
                        );
  }
  
  void pass_time_on_host(unsigned sleepMilliseconds)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(sleepMilliseconds));
  }

  template<typename Func>
  void set_fields_values_on_host(stk::mesh::FieldVector& fields, Func&& setValue)
  {
    for(auto field : fields) {
      stk::mesh::EntityVector elems;
      stk::mesh::get_selected_entities(stk::mesh::Selector(*field), get_bulk().buckets(stk::topology::ELEM_RANK), elems);

      for(auto elem : elems) {
        int* data = reinterpret_cast<int*>(stk::mesh::field_data(*field, elem));
        unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, elem);
        for(unsigned j = 0; j < numComponents; j++) {
          setValue(data, j);
        }
      }
    }
  }
  
  void update_fields_values_on_host(stk::mesh::FieldVector& fields)
  {
    auto updateValueFunc = [this](int* data, unsigned component)
                           {
                             data[component] += m_increment;
                           };

    set_fields_values_on_host(fields, updateValueFunc);
  }

  void reset_fields_values_on_host(stk::mesh::FieldVector& fields)
  {
    auto updateValueFunc = [](int* data, unsigned component)
                           {
                             data[component] = component;
                           };

    set_fields_values_on_host(fields, updateValueFunc);
  }

  void set_fields_values_on_device(stk::mesh::FieldVector& fields, unsigned value)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

    for(auto field : fields) {
      stk::mesh::NgpField<int> ngpField = stk::mesh::get_updated_ngp_field<int>(*field);

      stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, stk::mesh::Selector(*field),
                                     KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
                                       const int numScalarsPerEntity = ngpField.get_num_components_per_entity(entityIndex);
                                     
                                       for (int component = 0; component < numScalarsPerEntity; component++) {
                                         ngpField(entityIndex, component) = value;
                                       }
                                     });
    }
  }

  void verify_values_on_device(stk::mesh::FieldVector& fields)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    unsigned numComponents = m_numComponents;
    unsigned increments = m_increment;

    for(auto field : fields) {
      auto ngpField = stk::mesh::get_updated_ngp_field<int>(*field);

      stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, stk::mesh::Selector(*field),
                                    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                    {
                                      for(unsigned i = 0; i < numComponents; i++) {
                                        int expected = i + increments;
                                        int fieldValue = ngpField(elem, i);
                                        NGP_EXPECT_EQ(expected, fieldValue);
                                      }
                                    });
    }

    Kokkos::fence();
  }

  void verify_values_on_host(stk::mesh::FieldVector& fields, unsigned expectedValue)
  {
    for(auto field : fields) {
      auto ngpField = stk::mesh::get_updated_ngp_field<int>(*field);

      stk::mesh::EntityVector elems;
      stk::mesh::get_selected_entities(stk::mesh::Selector(*field), get_bulk().buckets(stk::topology::ELEM_RANK), elems);

      for(auto elem : elems) {
        int* data = reinterpret_cast<int*>(stk::mesh::field_data(*field, elem));
        unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, elem);
        for(unsigned j = 0; j < numComponents; j++) {
          EXPECT_EQ((int)expectedValue, data[j]);
        }
      }
    }
  }

private:
  void setup_field_component_data(std::vector<int>& init)
  {
    for(unsigned i = 0; i < m_numComponents; i++) {
      init.push_back(i);
    }
  }

  void revert_launch_blocking_env_var()
  {
    if(m_defaultLaunchBlockingEnvVarSet) {
      setenv(m_launchBlockingEnvVar.c_str(), std::to_string(m_defaultLaunchBlockingEnvVarValue).c_str(), 1);
    } else {
      unsetenv(m_launchBlockingEnvVar.c_str());
    }
  }

  void set_launch_blocking_env_var()
  {
    char* varValue = std::getenv(m_launchBlockingEnvVar.c_str());
    m_defaultLaunchBlockingEnvVarSet = (varValue != nullptr);

    if(m_defaultLaunchBlockingEnvVarSet) {
      m_defaultLaunchBlockingEnvVarValue = stk::get_env_var_as_int(m_launchBlockingEnvVar, 0);
      setenv(m_launchBlockingEnvVar.c_str(), "0", 1);
    }
  }

  unsigned m_numBlocks;
  unsigned m_numElemsPerDim;
  unsigned m_numElements;
  unsigned m_numComponents;
  unsigned m_increment;

  bool m_defaultLaunchBlockingEnvVarSet;
  int m_defaultLaunchBlockingEnvVarValue;
  const std::string m_launchBlockingEnvVar = "CUDA_LAUNCH_BLOCKING";
};

TEST_F(NgpFieldAsyncTest, SyncToDeviceAsyncTiming)
{
  if(get_parallel_size() != 1) return;

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 3);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 100);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();

  setup_simple_mesh_with_fields(numElemsPerDim);

  stk::mesh::FieldBase* intField1 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField1");
  stk::mesh::FieldBase* intField2 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField2");
  stk::mesh::FieldBase* intField3 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField3");
  stk::mesh::NgpField<int>& ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(*intField1);
  stk::mesh::NgpField<int>& ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(*intField2);
  stk::mesh::NgpField<int>& ngpIntField3 = stk::mesh::get_updated_ngp_field<int>(*intField3);
  stk::mesh::FieldVector fields{intField1, intField2, intField3};
  std::vector<stk::mesh::NgpField<int>*> ngpFields = {&ngpIntField1, &ngpIntField2, &ngpIntField3};

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();
      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_host();
        ngpField->sync_to_device();
        pass_time_on_device(defaultExecSpace, waitIteration);
        ngpField->fence();
      }

      verify_values_on_device(fields);
    }
    batchTimer.stop_batch_timer();
  }
  
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer2.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {
      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
      for(unsigned s = 0; s < numStreams; s++) {
        auto space = stk::mesh::get_execution_space_with_stream();
        spaces.push_back(space);
      }

      for(unsigned f = 0; f < ngpFields.size(); f++) {
        auto ngpField = ngpFields[f];
        auto space = spaces[f % spaces.size()];
        ngpField->modify_on_host();
        ngpField->sync_to_device(space);
        pass_time_on_device(space, waitIteration);
      }

      stk::mesh::ngp_field_fence(get_meta());
      verify_values_on_device(fields);
    }
    batchTimer2.stop_batch_timer();
  }

    double blockingSyncTime = batchTimer.get_min_batch_time();
    double nonBlockingSyncTime = batchTimer2.get_min_batch_time();
    double speedup = blockingSyncTime / nonBlockingSyncTime;

    EXPECT_GE(speedup, 1.0);
    EXPECT_LE(speedup, numStreams + SPEEDUP_DELTA);

  batchTimer2.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAsyncTest, SyncToHostAsyncTiming)
{
  if(get_parallel_size() != 1) return;

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 3);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 100);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();

  setup_simple_mesh_with_fields(numElemsPerDim);

  stk::mesh::FieldBase* intField1 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField1");
  stk::mesh::FieldBase* intField2 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField2");
  stk::mesh::FieldBase* intField3 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField3");
  stk::mesh::NgpField<int>& ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(*intField1);
  stk::mesh::NgpField<int>& ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(*intField2);
  stk::mesh::NgpField<int>& ngpIntField3 = stk::mesh::get_updated_ngp_field<int>(*intField3);
  stk::mesh::FieldVector fields{intField1, intField2, intField3};
  std::vector<stk::mesh::NgpField<int>*> ngpFields = {&ngpIntField1, &ngpIntField2, &ngpIntField3};

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      unsigned initialValue = 0;
      unsigned setValue = i+1;

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();
      set_fields_values_on_device(fields, initialValue);
      set_fields_values_on_device(fields, setValue);

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_device();
        ngpField->sync_to_host();
        pass_time_on_device(defaultExecSpace, waitIteration);
        ngpField->fence();
      }

      verify_values_on_host(fields, setValue);
    }
    batchTimer.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer2.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      unsigned initialValue = 0;
      unsigned setValue = i+1;

      set_fields_values_on_device(fields, initialValue);
      set_fields_values_on_device(fields, setValue);

      std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
      for(unsigned s = 0; s < numStreams; s++) {
        auto space = stk::mesh::get_execution_space_with_stream();
        spaces.push_back(space);
      }

      for(unsigned f = 0; f < ngpFields.size(); f++) {
        auto ngpField = ngpFields[f];
        auto space = spaces[f % spaces.size()];
        ngpField->modify_on_device();
        ngpField->sync_to_host(space);
        pass_time_on_device(space, waitIteration);
      }

      stk::mesh::ngp_field_fence(get_meta());
      verify_values_on_host(fields, setValue);
    }
    batchTimer2.stop_batch_timer();
  }

    double blockingSyncTime = batchTimer.get_min_batch_time();
    double nonBlockingSyncTime = batchTimer2.get_min_batch_time();
    double speedup = blockingSyncTime / nonBlockingSyncTime;

    EXPECT_GE(speedup, 1.0);
    EXPECT_LE(speedup, numStreams + SPEEDUP_DELTA);

  batchTimer2.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAsyncTest, SyncAsyncTiming)
{
  if(get_parallel_size() != 1) return;

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 3);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 100);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();
  
  setup_simple_mesh_with_fields(numElemsPerDim);

  stk::mesh::FieldBase* intField1 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField1");
  stk::mesh::FieldBase* intField2 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField2");
  stk::mesh::FieldBase* intField3 = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField3");
  stk::mesh::NgpField<int>& ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(*intField1);
  stk::mesh::NgpField<int>& ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(*intField2);
  stk::mesh::NgpField<int>& ngpIntField3 = stk::mesh::get_updated_ngp_field<int>(*intField3);
  stk::mesh::FieldVector fields{intField1, intField2, intField3};
  std::vector<stk::mesh::NgpField<int>*> ngpFields = {&ngpIntField1, &ngpIntField2, &ngpIntField3};

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();
      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_host();
        ngpField->sync_to_device();
        pass_time_on_device(defaultExecSpace, waitIteration);
        ngpField->modify_on_device();
        ngpField->sync_to_host();
        ngpField->fence();
      }

    }
    batchTimer.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer2.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {
      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
      for(unsigned s = 0; s < numStreams; s++) {
        auto space = stk::mesh::get_execution_space_with_stream();
        spaces.push_back(space);
      }

      for(unsigned f = 0; f < ngpFields.size(); f++) {
        auto ngpField = ngpFields[f];
        auto space = spaces[f % spaces.size()];

        ngpField->modify_on_host();
        ngpField->sync_to_device(space);
        pass_time_on_device(space, waitIteration);
        ngpField->modify_on_device();
        ngpField->sync_to_host(space);
      }

      stk::mesh::ngp_field_fence(get_meta());
    }
    batchTimer2.stop_batch_timer();
  }

    double blockingSyncTime = batchTimer.get_min_batch_time();
    double nonBlockingSyncTime = batchTimer2.get_min_batch_time();
    double speedup = blockingSyncTime / nonBlockingSyncTime;

    EXPECT_GE(speedup, 1.0);
    EXPECT_LE(speedup, numStreams + SPEEDUP_DELTA);

    batchTimer2.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAsyncTest, PartialSyncToDeviceAsyncTiming)
{
  if(get_parallel_size() != 1) return;

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 3);
  unsigned numFields = stk::unit_test_util::get_command_line_option("-f", 3);
  unsigned numBlocks = stk::unit_test_util::get_command_line_option("-b", 3);
  unsigned numBlocksToSync = stk::unit_test_util::get_command_line_option("-c", 1);
  EXPECT_TRUE(numBlocksToSync <= numBlocks && numBlocksToSync >= 1);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 100);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();

  setup_multi_block_mesh_with_field_per_block(numElemsPerDim, numBlocks, numFields);

  stk::mesh::FieldVector fields;
  std::vector<stk::mesh::NgpField<int>*> ngpFields;

  for(unsigned i = 1; i <= numFields; i++) {
    stk::mesh::FieldBase* intField = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField" + std::to_string(i));
    EXPECT_NE(nullptr, intField);
    fields.push_back(intField);

    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(*intField);
    ngpFields.push_back(&ngpIntField);
  }

  stk::mesh::Selector selector;
  for(unsigned i = 1; i <= numBlocksToSync; i++) {
    selector = selector | stk::mesh::Selector(*get_meta().get_part("block_" + std::to_string(i)));
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();
      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_host(selector);
        ngpField->sync_to_device();
        pass_time_on_device(defaultExecSpace, waitIteration);
        ngpField->fence();
      }

    }
    batchTimer.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer2.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      reset_fields_values_on_host(fields);
      update_fields_values_on_host(fields);

      std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
      for(unsigned s = 0; s < numStreams; s++) {
        auto space = stk::mesh::get_execution_space_with_stream();
        spaces.push_back(space);
      }

      for(unsigned f = 0; f < ngpFields.size(); f++) {
        auto ngpField = ngpFields[f];
        auto space = spaces[f % spaces.size()];

        ngpField->modify_on_host(selector);
        ngpField->sync_to_device(space);
        pass_time_on_device(space, waitIteration);
      }

      stk::mesh::ngp_field_fence(get_meta());
    }
    batchTimer2.stop_batch_timer();
  }

    double blockingSyncTime = batchTimer.get_min_batch_time();
    double nonBlockingSyncTime = batchTimer2.get_min_batch_time();
    double speedup = blockingSyncTime / nonBlockingSyncTime;

    EXPECT_GE(speedup, 1.0);
    EXPECT_LE(speedup, numStreams + SPEEDUP_DELTA);

  batchTimer2.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAsyncTest, PartialSyncToHostAsyncTiming)
{
  if(get_parallel_size() != 1) { GTEST_SKIP(); }

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 3);
  unsigned numFields = stk::unit_test_util::get_command_line_option("-f", 3);
  unsigned numBlocks = stk::unit_test_util::get_command_line_option("-b", 3);
  unsigned numBlocksToSync = stk::unit_test_util::get_command_line_option("-c", 1);
  EXPECT_TRUE(numBlocksToSync <= numBlocks && numBlocksToSync >= 1);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 100);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();

  setup_multi_block_mesh_with_field_per_block(numElemsPerDim, numBlocks, numFields);

  stk::mesh::FieldVector fields;
  std::vector<stk::mesh::NgpField<int>*> ngpFields;

  for(unsigned i = 1; i <= numFields; i++) {
    stk::mesh::FieldBase* intField = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField" + std::to_string(i));
    EXPECT_NE(nullptr, intField);
    fields.push_back(intField);

    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(*intField);
    ngpFields.push_back(&ngpIntField);
  }

  stk::mesh::Selector selector;
  for(unsigned i = 1; i <= numBlocksToSync; i++) {
    selector = selector | stk::mesh::Selector(*get_meta().get_part("block_" + std::to_string(i)));
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      unsigned initialValue = 0;
      unsigned setValue = i+1;

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();
      set_fields_values_on_device(fields, initialValue);
      set_fields_values_on_device(fields, setValue);

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_device(selector);
        ngpField->sync_to_host();
        pass_time_on_device(defaultExecSpace, waitIteration);
        ngpField->fence();
      }

    }
    batchTimer.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer2.start_batch_timer();

    for (unsigned i = 0; i < NUM_ITERS; i++) {

      unsigned initialValue = 0;
      unsigned setValue = i+1;

      set_fields_values_on_device(fields, initialValue);
      set_fields_values_on_device(fields, setValue);

      std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
      for(unsigned s = 0; s < numStreams; s++) {
        auto space = stk::mesh::get_execution_space_with_stream();
        spaces.push_back(space);
      }

      for(unsigned f = 0; f < ngpFields.size(); f++) {
        auto ngpField = ngpFields[f];
        auto space = spaces[f % spaces.size()];

        ngpField->modify_on_device(selector);
        ngpField->sync_to_host(space);
        pass_time_on_device(space, waitIteration);
      }

      stk::mesh::ngp_field_fence(get_meta());
    }
    batchTimer2.stop_batch_timer();
  }

  double blockingSyncTime = batchTimer.get_min_batch_time();
  double nonBlockingSyncTime = batchTimer2.get_min_batch_time();
  double speedup = blockingSyncTime / nonBlockingSyncTime;

  EXPECT_GE(speedup, 1.0);
  EXPECT_LE(speedup, numStreams + SPEEDUP_DELTA);

  batchTimer2.print_batch_timing(NUM_ITERS);
}

TEST_F(NgpFieldAsyncTest, AsyncDeepCopyTiming)
{
  if(get_parallel_size() != 1) { GTEST_SKIP(); }

  unsigned NUM_RUNS = 5;
  unsigned NUM_ITERS = stk::unit_test_util::get_command_line_option("-r", 50);
  unsigned numStreams = stk::unit_test_util::get_command_line_option("-s", 10);
  unsigned numFields = stk::unit_test_util::get_command_line_option("-f", 10);
  unsigned numBlocks = stk::unit_test_util::get_command_line_option("-b", 1);
  unsigned numElemsPerDim = stk::unit_test_util::get_command_line_option("-e", 100);
  unsigned sleepTime = stk::unit_test_util::get_command_line_option("-m", 50);
  unsigned waitIteration = stk::unit_test_util::get_command_line_option("-p", 20);
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer2(MPI_COMM_WORLD);
  stk::unit_test_util::BatchTimer batchTimer3(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  batchTimer2.initialize_batch_timer();
  batchTimer3.initialize_batch_timer();

  setup_multi_block_mesh_with_field_per_block(numElemsPerDim, numBlocks, numFields);

  stk::mesh::FieldVector fields;
  std::vector<stk::mesh::NgpField<int>*> ngpFields;

  for(unsigned i = 1; i <= numFields; i++) {
    stk::mesh::FieldBase* intField = get_meta().get_field(stk::topology::ELEMENT_RANK, "intField" + std::to_string(i));
    EXPECT_NE(nullptr, intField);
    fields.push_back(intField);

    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(*intField);
    ngpFields.push_back(&ngpIntField);
  }

  std::vector<stk::mesh::ExecSpaceWrapper<stk::ngp::ExecSpace>> spaces;
  for(unsigned i = 0; i < numStreams; i++) {
    auto space = stk::mesh::get_execution_space_with_stream();
    spaces.push_back(space);
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    reset_fields_values_on_host(fields);
    update_fields_values_on_host(fields);

    batchTimer.start_batch_timer();

    for (unsigned iter = 0; iter < NUM_ITERS; iter++) {

      Kokkos::fence();

      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_host();
        ngpField->sync_to_device();
        pass_time_on_device(defaultExecSpace, waitIteration);
      }
      pass_time_on_host(sleepTime);

    }
    batchTimer.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    reset_fields_values_on_host(fields);
    update_fields_values_on_host(fields);

    batchTimer2.start_batch_timer();

    for (unsigned iter = 0; iter < NUM_ITERS; iter++) {
    
      auto defaultExecSpace = Kokkos::DefaultExecutionSpace();

      for(auto ngpField : ngpFields) {
        ngpField->modify_on_host();
        ngpField->sync_to_device(defaultExecSpace);
        pass_time_on_device(defaultExecSpace, waitIteration);
      }
      pass_time_on_host(sleepTime);

      stk::mesh::ngp_field_fence(get_meta());
    }
    batchTimer2.stop_batch_timer();
  }

  for (unsigned j = 0; j < NUM_RUNS; j++) {

    reset_fields_values_on_host(fields);
    update_fields_values_on_host(fields);

    batchTimer3.start_batch_timer();

    for (unsigned iter = 0; iter < NUM_ITERS; iter++) {
    
      for(unsigned i = 0; i < ngpFields.size(); i++) {
        auto ngpField = ngpFields[i];
        auto space = spaces[i % spaces.size()];
        ngpField->modify_on_host();
        ngpField->sync_to_device(space);
        pass_time_on_device(space, waitIteration);
      }
      pass_time_on_host(sleepTime);

      stk::mesh::ngp_field_fence(get_meta());
    }
    batchTimer3.stop_batch_timer();
  }

  EXPECT_LE(batchTimer2.get_min_batch_time(), batchTimer.get_min_batch_time());
  EXPECT_LE(batchTimer3.get_min_batch_time(), batchTimer2.get_min_batch_time());

  batchTimer3.print_batch_timing(NUM_ITERS);
}
