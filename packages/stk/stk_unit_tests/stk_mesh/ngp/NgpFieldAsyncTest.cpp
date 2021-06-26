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
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpExecutionSpace.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/NgpSpaces.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/util/StkNgpVector.hpp>
#include "NgpUnitTestUtils.hpp"
#include <Kokkos_Core.hpp>
#include <string>
#include <cstdlib>

#ifndef KOKKOS_ENABLE_CUDA
#define TEST_ONLY_ON_CUDA(testname) DISABLED_##testname
#else
#define TEST_ONLY_ON_CUDA(testname) testname
#endif

class NgpAsyncDeepCopyFixture : public stk::unit_test_util::MeshFixture
{
public:
  NgpAsyncDeepCopyFixture()
  : m_numComponents(3),
    m_bucketCapacity(5),
    m_numBlocks(3),
    m_numFields(m_numBlocks),
    m_numStreams(3),
    m_multiplier(5),
    m_defaultLaunchBlockingEnvVarSet(false),
    m_defaultLaunchBlockingEnvVarValue(0)
  {
    set_launch_blocking_env_var();
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, m_bucketCapacity);
  }

  ~NgpAsyncDeepCopyFixture()
  {
    revert_launch_blocking_env_var();
  }

  std::vector<stk::mesh::ExecSpaceWrapper<>> get_execution_spaces_with_streams(unsigned numStreams)
  {
    std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces;

    for(unsigned i = 0; i < numStreams; i++) {
      execSpaces.push_back(stk::mesh::get_execution_space_with_stream());
    }

    return execSpaces;
  }

  void setup_multi_block_mesh_with_field_per_block()
  {
    std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(m_numBlocks);
    std::vector<double> coordinates = stk::unit_test_util::get_many_block_coordinates(m_numBlocks);

    setup_field_per_block();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
    construct_ngp_fields();
  }

  void setup_multi_block_mesh_with_fields_on_all_blocks()
  {
    std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(m_numBlocks);
    std::vector<double> coordinates = stk::unit_test_util::get_many_block_coordinates(m_numBlocks);

    setup_fields_on_all_blocks();
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
    construct_ngp_fields();
  }

  void construct_ngp_fields()
  {
    for(auto field : m_fields) {
      stk::mesh::get_updated_ngp_field<int>(*field);
    }
  }

  void setup_field_per_block()
  {
    unsigned numStates = 1;
    std::vector<int> init;
    setup_field_component_data(init);

    for(unsigned i = 1; i <= m_numBlocks; i++) {
      std::string blockName = "block_" + std::to_string(i);
      std::string fieldName = "field_" + std::to_string(i);

      stk::mesh::Part& part = get_meta().declare_part_with_topology(blockName, stk::topology::HEX_8);
      get_meta().set_part_id(part, i);
      EXPECT_NE(&part, nullptr);

      stk::mesh::Field<int>& field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, fieldName, numStates);
      m_fields.push_back(&field);
      stk::mesh::put_field_on_mesh(field, part, m_numComponents, init.data());
    }
  }

  void setup_fields_on_all_blocks()
  {
    unsigned numStates = 1;
    std::vector<int> init;
    setup_field_component_data(init);

    for(unsigned i = 1; i <= m_numBlocks; i++) {
      std::string blockName = "block_" + std::to_string(i);

      stk::mesh::Part& part = get_meta().declare_part_with_topology(blockName, stk::topology::HEX_8);
      get_meta().set_part_id(part, i);
      EXPECT_NE(&part, nullptr);

      for(unsigned j = 1; j <= m_numFields; j++) {
        std::string fieldName = "field_on_all_blocks_" + std::to_string(j);
        stk::mesh::Field<int>& field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, fieldName, numStates);
        stk::mesh::put_field_on_mesh(field, part, m_numComponents, init.data());
        m_fields.push_back(&field);
      }
    }
  }

  std::vector<stk::mesh::Field<int>*> get_fields()
  {
    return m_fields;
  } 

  stk::mesh::PartVector get_parts()
  {
    auto allParts = get_meta().get_parts();
    stk::mesh::PartVector elemParts;

    for(auto part : allParts) {
      if(part->primary_entity_rank() == stk::topology::ELEM_RANK && part->id() != stk::mesh::Part::INVALID_ID) {
        elemParts.push_back(part);
      }
    }

    return elemParts;
  }

  void setup_test(unsigned numBlocks, unsigned numStreams, unsigned multiplier)
  {
    m_numBlocks = numBlocks;
    m_numStreams = numStreams;
    m_multiplier = multiplier;
    m_numFields = m_numBlocks;
  }

  void setup_field_data_on_host()
  {
    for(auto field : m_fields) {
      set_element_field_data(*field, stk::mesh::Selector(*field), m_multiplier);
    }
  }

  void setup_selected_field_data_on_host(stk::mesh::Selector& selector)
  {
    for(auto field : m_fields) {
      set_element_field_data(*field, selector, m_multiplier);
    }
  }

  stk::mesh::Selector get_block_selector_for_partial_sync()
  {
    unsigned numPartialSyncBlocks = stk::unit_test_util::get_command_line_option("-b", m_numBlocks/2);

    numPartialSyncBlocks = std::min(numPartialSyncBlocks, m_numBlocks);

    stk::mesh::PartVector syncParts;
    for(unsigned i = 1; i <= numPartialSyncBlocks; i++) {
      std::string partName = "block_" + std::to_string(i);

      stk::mesh::Part* part = get_meta().get_part(partName);
      EXPECT_NE(part, nullptr);

      syncParts.push_back(part);
    }

    stk::mesh::Selector selector = stk::mesh::selectUnion(syncParts);
    return selector;
  }

  void setup_field_data_on_device()
  {
    for(auto field : m_fields) {
      auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
      set_element_field_data_on_device(ngpMesh, *field, stk::mesh::Selector(*field), m_multiplier);
    }
  }

  void setup_selected_field_data_on_device(stk::mesh::Selector& selector)
  {
    for(auto field : m_fields) {
      auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
      set_element_field_data_on_device(ngpMesh, *field, selector, m_multiplier);
    }
  }

  void add_parts_to_all_blocks(stk::mesh::PartVector& addParts)
  {
    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, entities);

    get_bulk().batch_change_entity_parts(entities, addParts, {});
  }

  void change_parts_on_selected_blocks(stk::mesh::PartVector& addParts, stk::mesh::PartVector& removeParts,
                                       stk::mesh::Selector& blockSelector)
  {
    stk::mesh::EntityVector entities;
    stk::mesh::get_selected_entities(blockSelector, get_bulk().buckets(stk::topology::ELEM_RANK), entities);

    get_bulk().batch_change_entity_parts(entities, addParts, removeParts);
  }

  void sync_fields_to_host()
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
      ngpField.modify_on_device();
      ngpField.sync_to_host();
    }
  }

  void sync_fields_to_host_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
      ngpField.modify_on_device();
      ngpField.sync_to_host(execSpaces[i % execSpaces.size()]);
    }
  }

  void sync_fields_to_host_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces,
                                 stk::mesh::Selector& selector)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
      ngpField.modify_on_device(selector);
      ngpField.sync_to_host(execSpaces[i % execSpaces.size()]);
    }
  }

  void sync_fields_to_device_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
      ngpField.modify_on_host();
      ngpField.sync_to_device(execSpaces[i % execSpaces.size()]);
    }
  }

  void sync_fields_to_device_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces,
                                   stk::mesh::Selector& selector)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
      ngpField.modify_on_host(selector);
      ngpField.sync_to_device(execSpaces[i % execSpaces.size()]);
    }
  }

  template<typename FieldType>
  void compare_device_data_to_init_data(FieldType& field)
  {
    unsigned numComponents = m_numComponents;
    unsigned multiplier = m_multiplier;
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, stk::mesh::Selector(*field),
                                  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                  {
                                    auto entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, elem);

                                    for(unsigned i = 0; i < numComponents; i++) {
                                      int expected = ngpMesh.identifier(entity) * multiplier + i;
                                      int fieldValue = ngpField(elem, i);
                                      NGP_EXPECT_EQ(expected, fieldValue);
                                    }
                                  });
    Kokkos::fence();
  }

  template<typename FieldType>
  void test_partial_copy_to_device_result(FieldType field, stk::mesh::Selector& selector)
  {
    unsigned numComponents = m_numComponents;
    unsigned multiplier = m_multiplier;
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                  {
                                    auto entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, elem);

                                    for(unsigned i = 0; i < numComponents; i++) {
                                      int expected = ngpMesh.identifier(entity) * multiplier + i;
                                      int fieldValue = ngpField(elem, i);
                                      NGP_EXPECT_EQ(expected, fieldValue);
                                    }
                                  });
    Kokkos::fence();

    stk::NgpVector<int> ngpVector;
    std::vector<int> init;
    setup_field_component_data(init);

    for(auto val : init) {
      ngpVector.push_back(val);
    }
    ngpVector.copy_host_to_device();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, !selector,
                                  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                  {
                                    for(unsigned i = 0; i < numComponents; i++) {
                                      int expected = ngpVector.device_get(i);
                                      int fieldValue = ngpField(elem, i);
                                      NGP_EXPECT_EQ(expected, fieldValue);
                                    }
                                  });
    Kokkos::fence();
  }

  template<typename FieldType>
  void compare_host_data_to_modified_data(FieldType field, unsigned scale = 1)
  {
    stk::mesh::EntityVector elems;
    stk::mesh::get_selected_entities(stk::mesh::Selector(*field), get_bulk().buckets(stk::topology::ELEM_RANK), elems);
    unsigned multiplier = m_multiplier * scale;

    for(auto elem : elems) {
      int* data = reinterpret_cast<int*>(stk::mesh::field_data(*field, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, elem);
      for(unsigned j = 0; j < numComponents; j++) {
        unsigned expectedValue = get_bulk().identifier(elem) * multiplier + j;
        EXPECT_EQ((int)expectedValue, data[j]);
      }
    }
  }

  template<typename FieldType>
  void test_partial_copy_to_host_result(FieldType field, stk::mesh::Selector& selector)
  {
    stk::mesh::EntityVector elems;
    stk::mesh::get_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK), elems);

    check_result_on_host_expect_multiplied_data(field, elems, m_multiplier);

    stk::mesh::Selector otherSelector = stk::mesh::Selector(*field) - selector;

    stk::mesh::EntityVector notSelectedElems;
    stk::mesh::get_selected_entities(otherSelector, get_bulk().buckets(stk::topology::ELEM_RANK), notSelectedElems);

    check_result_on_host_expect_init_data(field, notSelectedElems);
  }

  void set_element_field_data(stk::mesh::Field<int>& stkIntField, stk::mesh::Selector selector, unsigned multiplier)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK), elements);

    for(stk::mesh::Entity elem : elements) {
      int* data = reinterpret_cast<int*>(stk::mesh::field_data(stkIntField, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkIntField, elem);
      for(unsigned j = 0; j < numComponents; j++) {
        data[j] = get_bulk().identifier(elem) * multiplier + j;
      }
    }
  }

  void set_element_field_data_on_device(stk::mesh::NgpMesh& ngpMesh, stk::mesh::Field<int>& stkIntField,
                                        const stk::mesh::Selector& selector, unsigned multiplier)
  {
    stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
                                    const int numScalarsPerEntity = ngpField.get_num_components_per_entity(entityIndex);
                                    for (int component=0; component<numScalarsPerEntity; component++) {
                                      stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, entityIndex);
                                      ngpField(entityIndex, component) = ngpMesh.identifier(entity) * multiplier + component;
                                    }
                                  });
  }

private:
  template<typename FieldType, typename Func>
  void check_result_on_host(FieldType field, stk::mesh::EntityVector elems, Func&& testValues)
  {
    for(auto elem : elems) {
      int* data = reinterpret_cast<int*>(stk::mesh::field_data(*field, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, elem);
      for(unsigned j = 0; j < numComponents; j++) {
        testValues(data, elem, j);
      }
    }
  }

  template<typename FieldType>
  void check_result_on_host_expect_init_data(FieldType field, stk::mesh::EntityVector& elems)
  {
    auto expectInitData = [](int* data, stk::mesh::Entity entity, unsigned component)
                          {
                            int expectedValue = component;
                            EXPECT_EQ(data[component], expectedValue);
                          };

    check_result_on_host(field, elems, expectInitData);
  }

  template<typename FieldType>
  void check_result_on_host_expect_multiplied_data(FieldType field, stk::mesh::EntityVector& elems, unsigned multiplier)
  {
    auto expectInitData = [this, multiplier](int* data, stk::mesh::Entity elem, unsigned component)
                          {
                            int expectedValue = get_bulk().identifier(elem) * multiplier + component;
                            EXPECT_EQ(data[component], expectedValue);
                          };

    check_result_on_host(field, elems, expectInitData);
  }

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

  unsigned m_numComponents;
  unsigned m_bucketCapacity;
  unsigned m_numBlocks;
  unsigned m_numFields;
  unsigned m_numStreams;
  unsigned m_multiplier;

  bool m_defaultLaunchBlockingEnvVarSet;
  int m_defaultLaunchBlockingEnvVarValue;

  std::vector<stk::mesh::Field<int>*> m_fields;
  const std::string m_launchBlockingEnvVar = "CUDA_LAUNCH_BLOCKING";
};

TEST_F(NgpAsyncDeepCopyFixture, TwoStreamsAsyncSyncToDevice)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 2;
  unsigned numStreams = 2;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_device_data_to_init_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TwoStreamsAsyncSyncToHostFenceAll)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 2;
  unsigned numStreams = 2;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_device();

  sync_fields_to_host_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TwoStreamsAsyncSyncToHostFenceEach)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 2;
  unsigned numStreams = 2;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_device();

  sync_fields_to_host_async(execSpaces);

  for(auto field : get_fields()) {
    auto& ngpField = stk::mesh::get_updated_ngp_field<int>(*field);
    ngpField.fence();
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TwoStreamsAsyncSync)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);
  sync_fields_to_host_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncSyncUsingSameStreams)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 1;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);
  sync_fields_to_host_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_device_data_to_init_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncSyncToDeviceThenMeshMod)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 2;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  stk::mesh::Part& newPart = get_meta().declare_part("testPart", stk::topology::ELEMENT_RANK);
  stk::mesh::PartVector addParts(1, &newPart);

  add_parts_to_all_blocks(addParts);

  for(auto field : get_fields()) {
    compare_device_data_to_init_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncCopyFollowedBySyncCopy)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  sync_fields_to_host();

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncCopyFollowedByGetUpdatedNgpField)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_host();

  sync_fields_to_device_async(execSpaces);

  for(auto field : get_fields()) {
    stk::mesh::get_updated_ngp_field<int>(*field);
  }
  sync_fields_to_host();

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncSyncToHostFollowedByMeshMod)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_device();

  sync_fields_to_host_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  stk::mesh::Part& newPart = get_meta().declare_part("testPart", stk::topology::ELEMENT_RANK);
  stk::mesh::PartVector addParts(1, &newPart);

  add_parts_to_all_blocks(addParts);

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncSyncToHostFollowedByDataModOnHostThenGetUpdatedNgpField)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 1;
  unsigned numStreams = 1;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);

  setup_multi_block_mesh_with_field_per_block();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);

  setup_field_data_on_device();

  sync_fields_to_host_async(execSpaces);

  stk::mesh::ngp_field_fence(get_meta());

  unsigned scale = 3;

  for(auto field : get_fields()) {
    set_element_field_data(*field, stk::mesh::Selector(*field), multiplier*scale);
  }

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field, scale);
  }

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field, scale);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TEST_ONLY_ON_CUDA(ThreeStreamsAsyncPartialSyncToDevice))
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  auto selector = get_block_selector_for_partial_sync();

  stk::mesh::Selector allBlockSelector = get_meta().universal_part();
  setup_selected_field_data_on_host(allBlockSelector);

  sync_fields_to_device_async(execSpaces, selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_device_result(field, selector);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TEST_ONLY_ON_CUDA(FourStreamsAsyncPartialSyncToDevice_MeshModAfterPartialSync))
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 4;
  unsigned numStreams = 4;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  const stk::mesh::PartVector& parts = get_parts();
  auto block1Selector = stk::mesh::Selector(*parts[0]);
  auto block2Selector = stk::mesh::Selector(*parts[1]);
  auto block3Selector = stk::mesh::Selector(*parts[2]);
  auto block4Selector = stk::mesh::Selector(*parts[3]);
  auto allBlockSelector = stk::mesh::Selector(get_meta().universal_part());

  setup_selected_field_data_on_host(allBlockSelector);

  sync_fields_to_device_async(execSpaces, block1Selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_device_result(field, block1Selector);
  }

  stk::mesh::PartVector addParts(1, parts[1]);
  stk::mesh::PartVector removeParts(1, parts[3]);

  stk::mesh::Selector block1And2Selector = block1Selector | block2Selector;
  change_parts_on_selected_blocks(addParts, removeParts, block4Selector);

  for(auto field : get_fields()) {
    stk::mesh::get_updated_ngp_field<int>(*field);
  }

  for(auto field : get_fields()) {
    test_partial_copy_to_device_result(field, block1And2Selector);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TEST_ONLY_ON_CUDA(ThreeStreamsAsyncPartialSyncToHost))
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 3;
  unsigned numStreams = 3;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  auto selector = get_block_selector_for_partial_sync();

  stk::mesh::Selector allBlockSelector = get_meta().universal_part();
  setup_selected_field_data_on_device(allBlockSelector);

  sync_fields_to_host_async(execSpaces, selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_host_result(field, selector);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TEST_ONLY_ON_CUDA(FourStreamsAsyncPartialSyncToHost_MeshModAfterPartialSync))
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 4;
  unsigned numStreams = 4;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  const stk::mesh::PartVector& parts = get_parts();
  auto block1Selector = stk::mesh::Selector(*parts[0]);
  auto block2Selector = stk::mesh::Selector(*parts[1]);
  auto block3Selector = stk::mesh::Selector(*parts[2]);
  auto block4Selector = stk::mesh::Selector(*parts[3]);
  auto allBlockSelector = stk::mesh::Selector(get_meta().universal_part());

  setup_selected_field_data_on_device(allBlockSelector);

  sync_fields_to_host_async(execSpaces, block1Selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_host_result(field, block1Selector);
  }

  stk::mesh::PartVector addParts(1, parts[1]);
  stk::mesh::PartVector removeParts(1, parts[3]);

  stk::mesh::Selector block1And2Selector = block1Selector | block2Selector;
  change_parts_on_selected_blocks(addParts, removeParts, block4Selector);

  for(auto field : get_fields()) {
    stk::mesh::get_updated_ngp_field<int>(*field);
  }

  for(auto field : get_fields()) {
    test_partial_copy_to_host_result(field, block1Selector);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, TEST_ONLY_ON_CUDA(FourStreamsAsyncPartialSyncToDeviceThenPartialSyncToHost))
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 4;
  unsigned numStreams = 4;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  const stk::mesh::PartVector& parts = get_parts();
  auto block1Selector = stk::mesh::Selector(*parts[0]);
  auto block2Selector = stk::mesh::Selector(*parts[1]);
  auto block3Selector = stk::mesh::Selector(*parts[2]);
  auto block4Selector = stk::mesh::Selector(*parts[3]);
  auto allBlockSelector = stk::mesh::Selector(get_meta().universal_part());

  setup_selected_field_data_on_host(allBlockSelector);

  auto block1And3Selector = block1Selector | block3Selector;
  sync_fields_to_device_async(execSpaces, block1And3Selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_device_result(field, block1And3Selector);
  }

  auto block2And4Selector = block2Selector | block4Selector;
  sync_fields_to_host_async(execSpaces, block2And4Selector);

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    test_partial_copy_to_host_result(field, block1And3Selector);
  }
}

TEST_F(NgpAsyncDeepCopyFixture, AsyncGetUpdatedNgpField)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 2;
  unsigned numStreams = 2;
  unsigned multiplier = 5;

  setup_test(numBlocks, numStreams, multiplier);
  setup_multi_block_mesh_with_fields_on_all_blocks();
  std::vector<stk::mesh::ExecSpaceWrapper<>> execSpaces = get_execution_spaces_with_streams(numStreams);
  const stk::mesh::PartVector& parts = get_parts();
  auto block1Selector = stk::mesh::Selector(*parts[0]);
  auto block2Selector = stk::mesh::Selector(*parts[1]);
  auto allBlockSelector = stk::mesh::Selector(get_meta().universal_part());

  stk::mesh::PartVector addParts(1, parts[0]);
  stk::mesh::PartVector removeParts(1, parts[1]);

  change_parts_on_selected_blocks(addParts, removeParts, block2Selector);

  setup_field_data_on_host();

  auto fields = get_fields();
  for(unsigned i = 0; i < fields.size(); i++) {
    auto& ngpField = stk::mesh::get_updated_ngp_field_async<int>(*fields[i], execSpaces[i % execSpaces.size()]);
    ngpField.modify_on_host();
    ngpField.sync_to_device(execSpaces[i % execSpaces.size()]);
  }

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : fields) {
    test_partial_copy_to_device_result(field, allBlockSelector);
  }
}