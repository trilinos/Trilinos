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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpExecutionSpace.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/util/StkNgpVector.hpp>
#include "ngp/NgpUnitTestUtils.hpp"
#include <Kokkos_Core.hpp>
#include <string>
#include <cstdlib>

namespace {

#ifndef KOKKOS_ENABLE_CUDA
#define TEST_ONLY_ON_CUDA(testname) DISABLED_##testname
#else
#define TEST_ONLY_ON_CUDA(testname) testname
#endif

class UnitTestFieldAsync : public stk::unit_test_util::MeshFixture
{
public:
  UnitTestFieldAsync()
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
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, m_bucketCapacity, m_bucketCapacity);
  }

  ~UnitTestFieldAsync()
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
    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    construct_device_data();
  }

  void setup_multi_block_mesh_with_fields_on_all_blocks()
  {
    std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(m_numBlocks);
    std::vector<double> coordinates = stk::unit_test_util::get_many_block_coordinates(m_numBlocks);

    setup_fields_on_all_blocks();
    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    construct_device_data();
  }

  void construct_device_data()
  {
    for(auto field : m_fields) {
      field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
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

      stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::ELEM_RANK, fieldName, numStates);
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
        stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::ELEM_RANK, fieldName, numStates);
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
      field->synchronize<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    }
  }

  void sync_fields_to_host_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      field->synchronize<stk::mesh::ReadOnly, stk::ngp::HostSpace>(execSpaces[i % execSpaces.size()].get_execution_space());
    }
  }

  void sync_fields_to_device_async(std::vector<stk::mesh::ExecSpaceWrapper<>>& execSpaces)
  {
    for(unsigned i = 0; i < m_fields.size(); i++) {
      auto field = m_fields[i];
      field->synchronize<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(execSpaces[i % execSpaces.size()].get_execution_space());
    }
  }

  template<typename FieldType>
  void compare_device_data_to_init_data(FieldType& field)
  {
    int numComponents = m_numComponents;
    int multiplier = m_multiplier;
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto fieldData = field->template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, stk::mesh::Selector(*field),
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                   {
                                     auto entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, elem);
                                     auto fieldValues = fieldData.entity_values(elem);

                                     STK_NGP_ThrowRequire(numComponents == fieldValues.num_components());
                                     for (stk::mesh::ComponentIdx component : fieldValues.components()) {
                                       int expected = ngpMesh.identifier(entity) * multiplier + component;
                                       NGP_EXPECT_EQ(expected, fieldValues(component));
                                     }
                                   });
    Kokkos::fence();
  }

  template<typename FieldType>
  void test_partial_copy_to_device_result(FieldType field, stk::mesh::Selector& selector)
  {
    int multiplier = m_multiplier;
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto fieldData = field->template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                   {
                                     auto entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, elem);
                                     auto fieldValues = fieldData.entity_values(elem);

                                     for (stk::mesh::ComponentIdx component : fieldValues.components()) {
                                       int expected = ngpMesh.identifier(entity) * multiplier + component;
                                       NGP_EXPECT_EQ(expected, fieldValues(component));
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
                                     auto fieldValues = fieldData.entity_values(elem);

                                     for (stk::mesh::ComponentIdx component : fieldValues.components()) {
                                       int expected = ngpVector.device_get(component);
                                       NGP_EXPECT_EQ(expected, fieldValues(component));
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
    auto fieldData = field->template data<>();

    for(auto elem : elems) {
      auto fieldValues = fieldData.entity_values(elem);
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        int expectedValue = get_bulk().identifier(elem) * multiplier + component;
        EXPECT_EQ(expectedValue, fieldValues(component));
      }
    }
  }

  template<typename FieldType>
  void test_partial_copy_to_host_result(FieldType field, const stk::mesh::Selector& selector)
  {
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, selector, elems);

    check_result_on_host_expect_multiplied_data(field, elems, m_multiplier);

    stk::mesh::Selector otherSelector = stk::mesh::Selector(*field) - selector;

    stk::mesh::EntityVector notSelectedElems;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, otherSelector, notSelectedElems);

    check_result_on_host_expect_init_data(field, notSelectedElems);
  }

  void set_element_field_data(stk::mesh::Field<int>& stkIntField, const stk::mesh::Selector& selector, unsigned multiplier)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, selector, elements);
    auto fieldData = stkIntField.data<stk::mesh::ReadWrite>();

    for(stk::mesh::Entity elem : elements) {
      auto fieldValues = fieldData.entity_values(elem);
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        fieldValues(component) = get_bulk().identifier(elem) * multiplier + component;
      }
    }
  }

  void set_element_field_data_on_device(stk::mesh::NgpMesh& ngpMesh, stk::mesh::Field<int>& stkIntField,
                                        const stk::mesh::Selector& selector, unsigned multiplier)
  {
    auto fieldData = stkIntField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
                                     auto fieldValues = fieldData.entity_values(entityIndex);
                                     for (stk::mesh::ComponentIdx component : fieldValues.components()) {
                                       stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::ELEM_RANK, entityIndex);
                                       fieldValues(component) = ngpMesh.identifier(entity) * multiplier + component;
                                     }
                                   });
  }

private:
  template<typename FieldType, typename Func>
  void check_result_on_host(FieldType field, stk::mesh::EntityVector elems, Func&& testValues)
  {
    auto fieldData = field->data();
    for(auto elem : elems) {
      auto fieldValues = fieldData.entity_values(elem);
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        testValues(fieldValues(component), elem, component);
      }
    }
  }

  template<typename FieldType>
  void check_result_on_host_expect_init_data(FieldType field, stk::mesh::EntityVector& elems)
  {
    auto expectInitData = [](const int& value, stk::mesh::Entity /*entity*/, int component)
    {
      int expectedValue = component;
      EXPECT_EQ(value, expectedValue);
    };

    check_result_on_host(field, elems, expectInitData);
  }

  template<typename FieldType>
  void check_result_on_host_expect_multiplied_data(FieldType field, stk::mesh::EntityVector& elems, int multiplier)
  {
    auto expectInitData = [this, multiplier](const int& value, stk::mesh::Entity elem, int component)
    {
      int expectedValue = get_bulk().identifier(elem) * multiplier + component;
      EXPECT_EQ(value, expectedValue);
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

NGP_TEST_F(UnitTestFieldAsync, TwoStreamsAsyncSyncToDevice)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

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

NGP_TEST_F(UnitTestFieldAsync, TwoStreamsAsyncSyncToHostFenceAll)
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

NGP_TEST_F(UnitTestFieldAsync, TwoStreamsAsyncSyncToHostFenceEach)
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
    field->fence();
    compare_host_data_to_modified_data(field);
  }
}

NGP_TEST_F(UnitTestFieldAsync, TwoStreamsAsyncSync)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncSyncUsingSameStreams)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncSyncToDeviceThenMeshMod)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncCopyFollowedBySyncCopy)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncCopyFollowedByReacquiringFieldData)
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
    field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  }
  sync_fields_to_host();

  for(auto field : get_fields()) {
    compare_host_data_to_modified_data(field);
  }
}

NGP_TEST_F(UnitTestFieldAsync, AsyncSyncToHostFollowedByMeshMod)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncSyncToHostFollowedByDataModOnHostThenGetDeviceData)
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

NGP_TEST_F(UnitTestFieldAsync, AsyncAcquireFieldData)
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

  for(unsigned i = 0; i < get_fields().size(); i++) {
    auto field = get_fields()[i];
    field->data<stk::mesh::ReadWrite>();  // Modify on host
  }

  for(unsigned i = 0; i < get_fields().size(); i++) {
    auto field = get_fields()[i];
    auto execSpace = execSpaces[i % execSpaces.size()].get_execution_space();
    field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(execSpace);
  }

  stk::mesh::ngp_field_fence(get_meta());

  for(auto field : get_fields()) {
    compare_device_data_to_init_data(field);
  }
}

}
