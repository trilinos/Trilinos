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
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/DestroyRelations.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/NgpFieldBLAS.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include "NgpUnitTestUtils.hpp"
#include "NgpFieldTestUtils.hpp"
#include <Kokkos_Core.hpp>
#include <string>
#include <cstdlib>

namespace ngp_field_test {

class NgpFieldFixture : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  stk::mesh::Field<T> & create_field(stk::topology::rank_t rank, const std::string & name, unsigned numComponent = 1)
  {
    unsigned numStates = 1;
    const std::vector<T> init(numComponent, 1);
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), numComponent, init.data());
    return field;
  }

  void setup_one_field_one_element_mesh()
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stk::mesh::Field<int>& stkField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
    stk::mesh::Part& block = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);

    const int init1 = 1;
    stk::mesh::put_field_on_mesh(stkField, block, 1, &init1);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,6,block_1\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_two_field_two_element_mesh()
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stk::mesh::Field<int>& stkField1 = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "variableLengthField");
    stk::mesh::Field<int>& stkField2 = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "potentiallyOverwrittenField");
    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);

    const int init1[2] = {1, 1};
    stk::mesh::put_field_on_mesh(stkField1, block1, 1, init1);
    stk::mesh::put_field_on_mesh(stkField1, block2, 2, init1);

    const int init2 = 1;
    stk::mesh::put_field_on_mesh(stkField2, block1, 1, &init2);
    stk::mesh::put_field_on_mesh(stkField2, block2, 1, &init2);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,6,block_1\n"
                                 "0,2,SHELL_QUAD_4,2,3,4,5,block_2\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_two_variable_fields_two_element_mesh()
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);

    stk::mesh::Field<int>& stkField1 = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "variableLengthField");
    stk::mesh::Field<int>& stkField2 = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "partiallyAllocatedField");
    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);

    const int init1[] = { 1,  2,  3,
                         11, 12, 13,
                         21, 22, 23,
                         31, 32, 33,
                         41, 42, 43};
    stk::mesh::put_field_on_mesh(stkField1, block1, 3, 4, init1);
    stk::mesh::put_field_on_mesh(stkField1, block2, 3, 5, init1);

    const int init2[] = {1, 2};
    stk::mesh::put_field_on_mesh(stkField2, block1, 2, init2);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,6,block_1\n"
                                 "0,2,SHELL_QUAD_4,2,3,4,5,block_2\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  template<typename T>
  void setup_two_fields_five_hex_three_block_mesh(const int numComponent1, const int numComponent2, const int numStates = 1)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stk::mesh::Field<T>& stkField1 = get_meta().declare_field<T>(stk::topology::ELEM_RANK, "variableLengthField1", numStates);
    stk::mesh::Field<T>& stkField2 = get_meta().declare_field<T>(stk::topology::ELEM_RANK, "variableLengthField2", numStates);

    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);

    const std::vector<T> init1(numComponent1, -1);
    stk::mesh::put_field_on_mesh(stkField1, block1, numComponent1, init1.data());

    const std::vector<T> init2(numComponent2, -2);
    stk::mesh::put_field_on_mesh(stkField1, block2, numComponent2, init2.data());

    const std::vector<T> init3(numComponent1, -1);
    stk::mesh::put_field_on_mesh(stkField2, block1, numComponent1, init3.data());

    const std::vector<T> init4(numComponent2, -2);
    stk::mesh::put_field_on_mesh(stkField2, block2, numComponent2, init4.data());

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                                 "0,3,HEX_8,9,13,14,15,16,17,18,19,block_2\n"
                                 "0,4,HEX_8,9,20,21,22,23,24,25,26,block_2\n"
                                 "0,5,HEX_8,9,27,28,29,30,31,32,33,block_3";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  template<typename T>
  void setup_three_fields_five_hex_three_block_mesh(const int numComponent1, const int numComponent2, const int numStates = 1)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stk::mesh::Field<T>& stkField1 = get_meta().declare_field<T>(stk::topology::NODE_RANK, "variableLengthField1", numStates);
    stk::mesh::Field<T>& stkField2 = get_meta().declare_field<T>(stk::topology::NODE_RANK, "variableLengthField2", numStates);
    stk::mesh::Field<T>& stkField3 = get_meta().declare_field<T>(stk::topology::NODE_RANK, "variableLengthField3", numStates);

    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);

    const std::vector<T> init1(numComponent1, -1);
    stk::mesh::put_field_on_mesh(stkField1, block1, numComponent1, init1.data());

    const std::vector<T> init2(numComponent2, -2);
    stk::mesh::put_field_on_mesh(stkField1, block2, numComponent2, init2.data());

    const std::vector<T> init3(numComponent1, -1);
    stk::mesh::put_field_on_mesh(stkField2, block1, numComponent1, init3.data());

    const std::vector<T> init4(numComponent2, -2);
    stk::mesh::put_field_on_mesh(stkField2, block2, numComponent2, init4.data());

    const std::vector<T> init5(numComponent1, 0);
    stk::mesh::put_field_on_mesh(stkField3, block1, numComponent1, init5.data());

    const std::vector<T> init6(numComponent2, 0);
    stk::mesh::put_field_on_mesh(stkField3, block2, numComponent2, init6.data());

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                                 "0,3,HEX_8,9,13,14,15,16,17,18,19,block_2\n"
                                 "0,4,HEX_8,9,20,21,22,23,24,25,26,block_2\n"
                                 "0,5,HEX_8,9,27,28,29,30,31,32,33,block_3";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_two_element_mesh_field_on_each_element(stk::mesh::Field<int>& intField, stk::mesh::Field<double>& doubleField)
  {
    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);

    int initIntVal = 1;
    double initDoubleVal = 2.0;
    stk::mesh::put_field_on_mesh(intField, block1, 1, &initIntVal);
    stk::mesh::put_field_on_mesh(doubleField, block2, 1, &initDoubleVal);

    const std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,4,3,block_1\n"
                                 "0,2,SHELL_QUAD_4,2,5,6,4,block_2\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void add_3rd_element_to_2hex_3block_mesh()
  {
    get_bulk().modification_begin();
    stk::mesh::PartVector parts {get_meta().get_part("block_3")};
    stk::mesh::EntityIdVector nodeIds(8);
    std::iota(nodeIds.begin(), nodeIds.end(), 9);
    stk::mesh::declare_element(get_bulk(), parts, 3, nodeIds);
    get_bulk().modification_end();
  }

  template<typename T>
  void copy_fields_on_device(stk::mesh::NgpMesh& ngpMesh, stk::mesh::FieldBase* inputField,
                             stk::mesh::FieldBase* outputField, const stk::mesh::Selector& selector)
  {
    stk::mesh::NgpField<T> inputNgpField = stk::mesh::get_updated_ngp_field<T>(*inputField);
    stk::mesh::NgpField<T> outputNgpField = stk::mesh::get_updated_ngp_field<T>(*outputField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
                                     const int numScalarsPerEntity = inputNgpField.get_num_components_per_entity(entityIndex);

                                     for (int component=0; component<numScalarsPerEntity; component++) {
                                       outputNgpField(entityIndex, component) = inputNgpField(entityIndex, component);
                                     }
                                   });

    outputNgpField.modify_on_device();
  }

  template<typename T>
  void check_field_data_on_host(stk::mesh::Field<T>& stkField, const stk::mesh::Selector& selector, unsigned multiplier)
  {
    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), rank, selector, elements);

    for(stk::mesh::Entity element : elements) {
      T* data = stk::mesh::field_data(stkField, element);
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, element);
      for(unsigned j = 0; j < numComponents; j++) {
        int expectedVal = get_bulk().identifier(element) * multiplier + j;
        EXPECT_EQ(data[j], expectedVal);
      }
    }
  }

  template<typename T>
  void check_field_data_on_host(stk::mesh::Field<T>& stkField, unsigned multiplier)
  {
    check_field_data_on_host(stkField, stkField, multiplier);
  }

  template<typename T, typename FieldData, typename FieldDataMirror>
  void get_field_data_from_device(stk::mesh::EntityVector& elements, stk::mesh::NgpField<T>& ngpField,
                                  FieldData& deviceData, FieldDataMirror& hostData)
  {
    unsigned numElems = elements.size();
    stk::NgpVector<stk::mesh::Entity> ngpElements(numElems);

    for(unsigned i = 0; i < numElems; i++) {
      ngpElements[i] = elements[i];
    }
    ngpElements.copy_host_to_device();

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, numElems),
                         KOKKOS_LAMBDA(const int& elemIdx) {
                           stk::mesh::Entity elem = ngpElements.device_get(elemIdx);
                           auto meshIndex = ngpMesh.fast_mesh_index(elem);
                           unsigned numScalarsPerEntity = ngpField.get_num_components_per_entity(meshIndex);

                           for(unsigned i = 0; i < numScalarsPerEntity; i++) {
                             deviceData(elemIdx,i) = ngpField(meshIndex,i);
                           }
                         }
                         );
    Kokkos::deep_copy(hostData, deviceData);
  }

  template<typename T, typename FieldDataMirror, typename Func>
  void verify_field_data_on_device(const stk::mesh::EntityVector& elements, const stk::mesh::Field<T>& stkField,
                                   const FieldDataMirror& hostData, Func&& checkFunc)
  {
    for(unsigned i = 0; i < elements.size(); i++) {
      T* data = stk::mesh::field_data(stkField, elements[i]);
      unsigned numComponents = stk::mesh::field_scalars_per_entity(stkField, elements[i]);
      for(unsigned j = 0; j < numComponents; j++) {
        checkFunc(hostData(i,j), data[j]);
      }
    }
  }

  template<typename T, typename Func>
  void check_field_data_equality_on_device(stk::mesh::EntityVector& elements,
                                           stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField,
                                           Func&& checkFunc)
  {
    using FieldData = Kokkos::View<T**, Kokkos::LayoutRight, stk::ngp::MemSpace>;

    unsigned numElems = elements.size();
    unsigned numPerEntity = stkField.max_size();
    FieldData deviceData = FieldData("deviceData", numElems, numPerEntity);
    typename FieldData::HostMirror hostData = Kokkos::create_mirror_view(deviceData);

    get_field_data_from_device<T, FieldData, typename FieldData::HostMirror>(elements, ngpField, deviceData, hostData);
    verify_field_data_on_device<T, typename FieldData::HostMirror>(elements, stkField, hostData, checkFunc);
  }

  template<typename T>
  void check_mismatched_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField, stk::mesh::Selector& syncSelector)
  {
    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    fieldSelector &= syncSelector;
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(fieldSelector, get_bulk().buckets(rank), elements);
    auto checkFunc = [](T hostData, T stkData)
    {
      EXPECT_NE(hostData, stkData);
    };
    check_field_data_equality_on_device<T>(elements, ngpField, stkField, checkFunc);
  }

  template<typename T>
  void check_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField, stk::mesh::Selector& syncSelector)
  {
    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::Selector fieldSelector(stkField);
    fieldSelector &= syncSelector;
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(fieldSelector, get_bulk().buckets(rank), elements);
    auto checkFunc = [](T hostData, T stkData)
    {
      EXPECT_EQ(hostData, stkData);
    };
    check_field_data_equality_on_device<T>(elements, ngpField, stkField, checkFunc);
  }

  template<typename T>
  void check_field_data_on_device(stk::mesh::NgpField<T>& ngpField, stk::mesh::Field<T>& stkField)
  {
    stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), rank, elements);
    auto checkFunc = [](T hostData, T stkData)
    {
      EXPECT_EQ(hostData, stkData);
    };
    check_field_data_equality_on_device<T>(elements, ngpField, stkField, checkFunc);
  }

  void set_element_field_data(stk::mesh::FieldBase* field)
  {
    stk::mesh::EntityVector elements;
    const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
    stk::mesh::get_selected_entities(get_meta().universal_part(), buckets, elements);

    for(stk::mesh::Entity elem : elements) {
      int* data = static_cast<int*>(stk::mesh::field_data(*field, elem));
      unsigned numComponents = stk::mesh::field_scalars_per_entity(*field, get_bulk().bucket(elem));
      for(unsigned i = 0; i < numComponents; i++) {
        data[i] = get_bulk().identifier(elem) * 10 + i;
      }
    }
  }
};

class OptimizedNgpFieldFixture : public NgpFieldFixture
{

public:
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

  void setup_3hex_3block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    ngp_unit_test_utils::setup_mesh_3hex_3block(get_bulk(), bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_3hex_2block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    ngp_unit_test_utils::setup_mesh_3hex_2block(get_bulk(), bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void setup_2hex_3block_mesh_with_field(unsigned bucketCapacity, stk::mesh::Field<int>& stkIntField)
  {
    get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);
    ngp_unit_test_utils::setup_mesh_2hex_2block(get_bulk(), bucketCapacity);
    set_element_field_data(stkIntField, get_meta().universal_part(), 10u);
  }

  void run_add_and_delete_bucket3(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket3(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_and_delete_bucket2(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket2(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_and_delete_bucket(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_2hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_add_and_delete_bucket(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }


  void run_delete_bucket_in_middle(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    modify_mesh_delete_bucket_in_middle(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_add_new_element(unsigned numComponents, unsigned bucketCapacity)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_element(stkIntField, ngpIntField, bucketCapacity);
  }

  void run_add_new_element_and_modify_from_device(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    {
      stk::mesh::get_updated_ngp_mesh(get_bulk());
      stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

      modify_and_test_add_element(stkIntField, ngpIntField, bucketCapacity);
    }

    {
      stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
      stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

      check_field_data_on_host<int>(stkIntField, 10u);

      ngpIntField.sync_to_device();
      ngpIntField.clear_sync_state();
      ngpIntField.modify_on_device();

      stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEMENT_RANK, get_meta().universal_part(),
                                     KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
                                       ngpIntField(elem, 0) *= 2;
                                     }
                                     );

      ngpIntField.sync_to_host();

      check_field_data_on_host<int>(stkIntField, 20u);
    }
  }

  void run_add_bucket_in_middle_copy(unsigned numComponents)
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void run_add_bucket_in_middle_external(unsigned numComponents)
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField(get_bulk(), stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void run_add_bucket_in_middle_internal(unsigned numComponents)
  {
    const unsigned bucketCapacity = 1;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_3block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    modify_and_test_add_bucket_in_middle(stkIntField, ngpIntField);
  }

  void modify_and_test_add_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    modify_mesh_add_bucket_in_middle(stkIntField, ngpIntField);

    ngpIntField.modify_on_host();
    stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void modify_and_test_add_element(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField, unsigned bucketCapacity)
  {
    modify_mesh_add_element(stkIntField, ngpIntField, bucketCapacity);

    ngpIntField.modify_on_host();
    stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_change_bucket_content_by_mesh_modification(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int> ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    modify_mesh_change_bucket_content(stkIntField, ngpIntField);

    ngpIntField.clear_sync_state();
    ngpIntField.modify_on_host();
    ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void run_change_bucket_content_by_user(unsigned numComponents)
  {
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
    stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

    setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);
    stk::mesh::NgpField<int>& ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    set_element_field_data(stkIntField, get_meta().universal_part(), 20u);

    ngpIntField.modify_on_host();
    ngpIntField.sync_to_device();

    check_field_data_on_device<int>(ngpIntField, stkIntField);
  }

  void modify_mesh_add_and_delete_bucket(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_3")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 1), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {1}}, {"block_2", {2}}});
  }

  void fill_nodes(const stk::mesh::Entity element, unsigned numNodes, stk::mesh::EntityVector& nodes)
  {
    const stk::mesh::Entity* elementNodes = get_bulk().begin_nodes(element);

    for(unsigned i = 0; i < numNodes; i++) {
      nodes.push_back(elementNodes[i]);
    }
  }

  stk::mesh::Entity add_element_and_place_in_block(const std::string& newBlockName)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
    stk::mesh::EntityVector nodes;

    unsigned numNodes = get_bulk().num_nodes(element);

    fill_nodes(element, numNodes, nodes);

    stk::mesh::Part* part = get_meta().get_part(newBlockName);
    stk::mesh::PartVector parts = {part};
    parts.push_back(&get_meta().get_topology_root_part(stk::topology::HEX_8));
    stk::mesh::Entity newElement = get_bulk().declare_element(4, parts);

    for(unsigned i = 0; i < numNodes; i++) {
      get_bulk().declare_relation(newElement, nodes[i], i);
    }

    return newElement;
  }

  void replace_element_and_place_in_block(const std::string& newBlockName)
  {
    add_element_and_place_in_block(newBlockName);

    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
    stk::mesh::destroy_relations(get_bulk(), element, stk::topology::NODE_RANK);
    get_bulk().destroy_entity(element);
  }

  void modify_host_bucket_value_with_selector(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField,
                                              stk::mesh::Selector selector, int newMultiplier)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector);
    set_element_field_data(stkIntField, selector, newMultiplier);
    ngpMesh.update_mesh();
  }

  void modify_device_bucket_value_with_selector(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField,
                                                stk::mesh::Selector selector, int newMultiplier)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    check_field_data_on_device<int>(ngpIntField, stkIntField, selector);
    set_element_field_data_on_device(ngpMesh, stkIntField, selector, newMultiplier);
    // ngpMesh.update_mesh();
  }

  void modify_mesh_add_and_delete_bucket3(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    replace_element_and_place_in_block("block_3");
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {4}}, {"block_2", {2}}});
  }

  void modify_mesh_add_and_delete_bucket2(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_2")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 1), addParts, removeParts);
    addParts[0] = get_meta().get_part("block_3");
    removeParts[0] = get_meta().get_part("block_2");
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_3", {2}}, {"block_2", {1}}});
  }

  void modify_mesh_delete_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_2")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1,2}}, {"block_3", {3}}});
  }

  void modify_mesh_add_bucket_in_middle(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts{get_meta().get_part("block_1")};
    stk::mesh::PartVector removeParts{get_meta().get_part("block_3")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 3), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_1", {3}}, {"block_2", {2}}});
  }

  void modify_mesh_add_element(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField, unsigned bucketCapacity)
  {
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);

    get_bulk().modification_begin();
    stk::mesh::Entity newElement = add_element_and_place_in_block("block_3");
    get_bulk().modification_end();

    int* data = stk::mesh::field_data(stkIntField, newElement);
    *data = get_bulk().identifier(newElement) * 10u;
    ngpMesh.update_mesh();

    if(bucketCapacity == 1) {
      ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3}}, {"block_3", {4}}});
    }
    else if(bucketCapacity == 2) {
      ngp_unit_test_utils::check_bucket_layout(get_bulk(), {{"block_1", {1}}, {"block_2", {2}}, {"block_3", {3,4}}});
    }
  }

  void modify_mesh_change_bucket_content(stk::mesh::Field<int>& stkIntField, stk::mesh::NgpField<int>& ngpIntField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), { {"block_1", {1, 2}}, {"block_3", {3}}});
    check_field_data_on_device<int>(ngpIntField, stkIntField);
    get_bulk().modification_begin();
    stk::mesh::PartVector addParts {get_meta().get_part("block_3")};
    stk::mesh::PartVector removeParts {get_meta().get_part("block_1")};
    get_bulk().change_entity_parts(get_bulk().get_entity(stk::topology::ELEM_RANK, 2), addParts, removeParts);
    get_bulk().modification_end();
    ngpMesh.update_mesh();
    ngp_unit_test_utils::check_bucket_layout(get_bulk(), { {"block_1", {1}}, {"block_3", {2, 3}}});
  }
};

void move_data_between_fields_on_host(const stk::mesh::BulkData & bulk,
                                      const stk::mesh::Field<int>& source,
                                      stk::mesh::Field<int>& dest)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEMENT_RANK);
  stk::mesh::NgpField<int>& ngpSource = stk::mesh::get_updated_ngp_field<int>(source);
  ngpSource.sync_to_host();

  for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
  {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];

    int* sourceData = static_cast<int*>(stk::mesh::field_data(source, bucket));
    int* destData   = static_cast<int*>(stk::mesh::field_data(dest, bucket));
    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
    {
      *destData = *sourceData;
    }
  }

  stk::mesh::NgpField<int>& ngpDest = stk::mesh::get_updated_ngp_field<int>(dest);
  ngpDest.modify_on_host();
}

struct CheckFieldValues {
  CheckFieldValues(stk::mesh::NgpField<int>& _ngpField, unsigned _numScalarsPerEntity, int _expectedFieldValue)
    : ngpField(_ngpField), numScalarsPerEntity(_numScalarsPerEntity), expectedFieldValue(_expectedFieldValue)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entity) const
  {
    for (unsigned component = 0; component < numScalarsPerEntity; component++) {
      STK_NGP_ThrowRequire(ngpField(entity, component) == expectedFieldValue);
    }
  }

private:
  stk::mesh::NgpField<int> ngpField;
  unsigned numScalarsPerEntity;
  int expectedFieldValue;
};

void test_field_values_on_device_without_initial_sync(stk::mesh::BulkData& bulk,
                                                      const stk::mesh::Field<int>& stkField,
                                                      const stk::mesh::Part& part,
                                                      const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField = stk::mesh::get_updated_ngp_field<int>(stkField);

  stk::mesh::Selector selection = bulk.mesh_meta_data().locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, *buckets[0]);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  CheckFieldValues checkFieldValues(ngpField, numScalarsPerEntity, expectedFieldValue);
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selection, checkFieldValues);
}

void test_field_values_on_device(stk::mesh::BulkData& bulk,
                                 const stk::mesh::Field<int>& stkField,
                                 const stk::mesh::Part& part,
                                 const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField = stk::mesh::get_updated_ngp_field<int>(stkField);
  ngpField.sync_to_device();

  test_field_values_on_device_without_initial_sync(bulk, stkField, part, expectedFieldValue);
}

void test_field_values_on_device(stk::mesh::BulkData &bulk,
                                 const stk::mesh::Field<int> & stkField,
                                 const int expectedFieldValue)
{
  test_field_values_on_device(bulk, stkField, bulk.mesh_meta_data().locally_owned_part(), expectedFieldValue);
}

void test_field_values_on_host_without_initial_sync(const stk::mesh::BulkData& bulk,
                                                    const stk::mesh::Field<int>& stkField,
                                                    const stk::mesh::Part& part,
                                                    const int expectedFieldValue)
{
  stk::mesh::Selector selection = bulk.mesh_meta_data().locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  for (size_t iBucket=0; iBucket<buckets.size(); iBucket++) {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];
    const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, bucket);

    int* fieldData = reinterpret_cast<int*>(stk::mesh::field_data(stkField, bucket));
    for (size_t iEntity=0; iEntity<bucket.size(); iEntity++) {
      for (unsigned component=0; component<numScalarsPerEntity; component++) {
        EXPECT_EQ(expectedFieldValue, fieldData[component]);
      }
    }
  }
}

void test_field_values_on_host(const stk::mesh::BulkData& bulk,
                               const stk::mesh::Field<int>& stkField,
                               const stk::mesh::Part& part,
                               const int expectedFieldValue)
{
  stkField.sync_to_host();

  test_field_values_on_host_without_initial_sync(bulk, stkField, part, expectedFieldValue);
}

void test_field_values_on_host(const stk::mesh::BulkData& bulk,
                               const stk::mesh::Field<int>& stkField,
                               const int expectedFieldValue)
{
  test_field_values_on_host(bulk, stkField, bulk.mesh_meta_data().locally_owned_part(), expectedFieldValue);
}

template <typename T>
void initialize_ngp_field(stk::mesh::Field<T> & stkField)
{
  stk::mesh::get_updated_ngp_field<T>(stkField);
}

template <typename T>
void multiply_field_data_on_host(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & field, int multiplier)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(field.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(field, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      fieldData[iEntity] *= multiplier;
    }
  }
}

template <typename T>
void modify_field_on_host(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & field, int multiplier)
{
  multiply_field_data_on_host(bulk, field, multiplier);
  field.modify_on_host();
}

void modify_mesh_part_membership(stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  stk::mesh::PartVector addParts(1, &part);
  stk::mesh::EntityVector allElements;

  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, allElements);

  bulk.modification_begin();
  bulk.change_entity_parts(allElements[0], addParts);
  bulk.modification_end();
}

template <typename T>
void sync_field_to_device(stk::mesh::Field<T> & field)
{
  field.sync_to_device();
}

template <typename T>
void multiply_field_data_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, int multiplier)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   const int numScalarsPerEntity = ngpField.get_num_components_per_entity(entity);
                                   for (int component=0; component<numScalarsPerEntity; component++) {
                                     ngpField(entity, component) *= multiplier;
                                   }
                                 });
}

template <typename T>
void modify_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, int multiplier)
{
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  multiply_field_data_on_device(bulk, stkField, multiplier);

  stkField.modify_on_device();
}

template <typename T>
void modify_field_on_device(stk::mesh::BulkData& bulk, stk::mesh::Field<T>& stkField, stk::mesh::Part& part, int value)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  stk::mesh::Selector selection = meta.locally_owned_part() & part;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selection);
  const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(stkField, *buckets[0]);
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selection,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   for (unsigned component=0; component<numScalarsPerEntity; component++) {
                                     ngpField(entity, component) = value;
                                   }
                                 });
  stkField.modify_on_device();
}

template <typename T>
void sync_field_to_host(stk::mesh::Field<T> & stkField)
{
  stkField.sync_to_host();
}

template <typename T>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         const stk::mesh::Field<T>& stkField,
                         int expectedValue)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      EXPECT_EQ(fieldData[iEntity], expectedValue);
    }
  }
}

TEST_F(NgpFieldFixture, noFieldDataTest)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  unsigned numBlocks = 1;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(numBlocks);
  std::vector<double> coordinates = stk::unit_test_util::get_many_block_coordinates(numBlocks);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "", 1);
  stk::unit_test_util::setup_text_mesh(get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_field<int>(field));
}

TEST_F(NgpFieldFixture, ModifyOnHostFlagClearedOnInitialNgpFieldConstruction)
{
  if (get_parallel_size() != 1) return;

  setup_one_field_one_element_mesh();

  stk::mesh::Field<int>& field1 = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "field1");
  field1.modify_on_host();

  auto ngpfield = stk::mesh::get_updated_ngp_field<int>(field1);
  EXPECT_FALSE(field1.need_sync_to_device());
}

TEST_F(NgpFieldFixture, InvalidModifyFlagCondition)
{
  if (get_parallel_size() != 1) return;

  setup_one_field_one_element_mesh();

  stk::mesh::Field<int>& field1 = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "field1");

  auto ngpfield = stk::mesh::get_updated_ngp_field<int>(field1);
  EXPECT_FALSE(field1.need_sync_to_device());

  field1.modify_on_host();
  EXPECT_THROW(ngpfield.modify_on_device(), std::logic_error);
}

TEST_F(NgpFieldFixture, PersistentModifyOnDeviceFlag)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  setup_one_field_one_element_mesh();

  stk::mesh::Field<int>& field1 = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "field1");
  EXPECT_FALSE(field1.need_sync_to_host());
  field1.sync_to_device();
  field1.modify_on_device();

  auto ngpfield = stk::mesh::get_updated_ngp_field<int>(field1);
  EXPECT_TRUE(field1.need_sync_to_host());

  multiply_field_data_on_device(get_bulk(), field1, 2);
  ngpfield.sync_to_host();

  test_field_values_on_host_without_initial_sync(get_bulk(), field1, get_meta().universal_part(), 2);
}

TEST_F(NgpFieldFixture, noOverwriteInVariableLengthFields)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  setup_two_field_two_element_mesh();

  stk::mesh::Field<int>& variableLengthField = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "variableLengthField");
  stk::mesh::Field<int>& potentiallyOverwrittenField = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "potentiallyOverwrittenField");

  modify_field_on_device(get_bulk(), variableLengthField, 2);
  modify_field_on_host(get_bulk(), potentiallyOverwrittenField, 3);
  variableLengthField.sync_to_host();

  test_field_values_on_host(get_bulk(), potentiallyOverwrittenField, 3);
}

TEST_F(NgpFieldFixture, FieldCopyVariableLengthField)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numComponent1 = 8;
  const int numComponent2 = 3;
  setup_two_fields_five_hex_three_block_mesh<int>(numComponent1, numComponent2);

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto stkField1 = get_meta().get_field<int>(stk::topology::ELEM_RANK, "variableLengthField1");
  auto stkField2 = get_meta().get_field<int>(stk::topology::ELEM_RANK, "variableLengthField2");
  auto elem2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
  auto elem4 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 4);

  int* inData2 = stk::mesh::field_data(*stkField1, elem2);
  int* inData4 = stk::mesh::field_data(*stkField1, elem4);

  for(int c = 0; c < numComponent1; ++c) {
    inData2[c] = 2 + c;
  }
  for(int c = 0; c < numComponent2; ++c) {
    inData4[c] = numComponent1 + 3 + c;
  }

  copy_fields_on_device<int>(ngpMesh, stkField1, stkField2, get_meta().universal_part());
  stkField2->sync_to_host();

  int* outData2 = stk::mesh::field_data(*stkField2, elem2);
  int* outData4 = stk::mesh::field_data(*stkField2, elem4);

  for(int c = 0; c < numComponent1; ++c) {
    EXPECT_EQ(2 + c, outData2[c]);
  }
  for(int c = 0; c < numComponent2; ++c) {
    EXPECT_EQ(numComponent1 + 3 + c, outData4[c]);
  }
}

#ifdef STK_USE_DEVICE_MESH
TEST_F(NgpFieldFixture, DeviceField_set_all_after_modified_on_host)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numComponent1 = 8;
  const int numComponent2 = 3;
  setup_two_fields_five_hex_three_block_mesh<double>(numComponent1, numComponent2);

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto stkField1 = get_meta().get_field<double>(stk::topology::ELEM_RANK, "variableLengthField1");

  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField1->need_sync_to_device());

  stk::mesh::NgpField<double> ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
 
  ngpField1.modify_on_host();
  EXPECT_TRUE(stkField1->need_sync_to_device());

  ngpField1.set_all(ngpMesh, 97.9);

  EXPECT_TRUE(ngpField1.need_sync_to_host());
}

TEST_F(NgpFieldFixture, blas_field_copy_device_to_device)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const int numComponent1 = 8;
  const int numComponent2 = 3;
  setup_two_fields_five_hex_three_block_mesh<double>(numComponent1, numComponent2);

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto stkField1 = get_meta().get_field<double>(stk::topology::ELEM_RANK, "variableLengthField1");
  auto stkField2 = get_meta().get_field<double>(stk::topology::ELEM_RANK, "variableLengthField2");

  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
  EXPECT_FALSE(stkField1->need_sync_to_device());
  EXPECT_FALSE(stkField2->need_sync_to_device());

 
  const double myConstantValue = 97.9;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());

#ifdef STK_ENABLE_GPU
  stk::mesh::NgpField<double>& ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
  EXPECT_TRUE(ngpField1.need_sync_to_host());
#endif

  stk::mesh::field_copy(*stkField1, *stkField2);

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#endif

  stk::mesh::Selector selector(*stkField2);
  stk::mesh::NgpField<double>& ngpField2 = stk::mesh::get_updated_ngp_field<double>(*stkField2);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField2, selector, myConstantValue);
}

#endif

void check_expected_num_elements(const stk::mesh::BulkData & bulk, unsigned numElements)
{
  std::vector<size_t> counts;
  const stk::mesh::Part & locallyOwnedPart = bulk.mesh_meta_data().locally_owned_part();
  stk::mesh::count_entities(locallyOwnedPart, bulk, counts);
  ASSERT_EQ(counts[stk::topology::ELEM_RANK], numElements);
}

class NumScalarsPerEntity
{
public:
  NumScalarsPerEntity() = default;
  NumScalarsPerEntity(unsigned scalar0, unsigned scalar1)
  {
    m_numScalars[0] = scalar0;
    m_numScalars[1] = scalar1;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned& operator[](unsigned entityIndex)
  {
    return m_numScalars[entityIndex];
  }

  KOKKOS_INLINE_FUNCTION
  const unsigned& operator[](const unsigned entityIndex) const
  {
    return m_numScalars[entityIndex];
  }

private:
  unsigned m_numScalars[2];
};

void fill_gold_num_scalars_per_entity(const stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & variableLengthField,
                                      unsigned numElements, NumScalarsPerEntity& goldNumScalarsPerEntity)
{
  stk::mesh::EntityVector elements;
  stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), elements);
  for (const stk::mesh::Entity & element : elements) {
    const unsigned index = bulk.identifier(element) - 1;
    ASSERT_LT(index, numElements);
    goldNumScalarsPerEntity[index] = stk::mesh::field_scalars_per_entity(variableLengthField, element);
  }
}

struct CheckNumScalarsPerEntity {
  CheckNumScalarsPerEntity(stk::mesh::NgpMesh& _ngpMesh,
                           stk::mesh::NgpField<int>& _ngpVariableLengthField,
                           NumScalarsPerEntity _goldNumScalarsPerEntity)
    : ngpMesh(_ngpMesh),
      ngpVariableLengthField(_ngpVariableLengthField),
      goldNumScalarsPerEntity(_goldNumScalarsPerEntity)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& index) const
  {
    const stk::mesh::Entity element = ngpMesh.get_entity(stk::topology::ELEM_RANK, index);
    const unsigned goldIndex = ngpMesh.identifier(element) - 1;
    NGP_EXPECT_EQ(ngpVariableLengthField.get_num_components_per_entity(index), goldNumScalarsPerEntity[goldIndex]);
  }

private:
  stk::mesh::NgpMesh ngpMesh;
  stk::mesh::NgpField<int> ngpVariableLengthField;
  NumScalarsPerEntity goldNumScalarsPerEntity;
};

struct CheckExtentsPerEntity {
  CheckExtentsPerEntity(stk::mesh::NgpMesh& _ngpMesh,
                        stk::mesh::NgpField<int>& _ngpVariableLengthField,
                        NumScalarsPerEntity _goldVariableLengthExtent0PerEntity,
                        NumScalarsPerEntity _goldVariableLengthExtent1PerEntity,
                        stk::mesh::NgpField<int>& _ngpPartiallyAllocatedField,
                        NumScalarsPerEntity _goldPartiallyAllocatedExtent0PerEntity,
                        NumScalarsPerEntity _goldPartiallyAllocatedExtent1PerEntity)
    : ngpMesh(_ngpMesh),
      ngpVariableLengthField(_ngpVariableLengthField),
      goldVariableLengthExtent0PerEntity(_goldVariableLengthExtent0PerEntity),
      goldVariableLengthExtent1PerEntity(_goldVariableLengthExtent1PerEntity),
      ngpPartiallyAllocatedField(_ngpPartiallyAllocatedField),
      goldPartiallyAllocatedExtent0PerEntity(_goldPartiallyAllocatedExtent0PerEntity),
      goldPartiallyAllocatedExtent1PerEntity(_goldPartiallyAllocatedExtent1PerEntity)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& index) const
  {
    const stk::mesh::Entity element = ngpMesh.get_entity(stk::topology::ELEM_RANK, index);
    const unsigned goldIndex = ngpMesh.identifier(element) - 1;

    NGP_EXPECT_EQ(ngpVariableLengthField.get_extent0_per_entity(index), goldVariableLengthExtent0PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpVariableLengthField.get_extent1_per_entity(index), goldVariableLengthExtent1PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpVariableLengthField.get_extent_per_entity(index, 0), goldVariableLengthExtent0PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpVariableLengthField.get_extent_per_entity(index, 1), goldVariableLengthExtent1PerEntity[goldIndex]);

    NGP_EXPECT_EQ(ngpPartiallyAllocatedField.get_extent0_per_entity(index), goldPartiallyAllocatedExtent0PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpPartiallyAllocatedField.get_extent1_per_entity(index), goldPartiallyAllocatedExtent1PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpPartiallyAllocatedField.get_extent_per_entity(index, 0), goldPartiallyAllocatedExtent0PerEntity[goldIndex]);
    NGP_EXPECT_EQ(ngpPartiallyAllocatedField.get_extent_per_entity(index, 1), goldPartiallyAllocatedExtent1PerEntity[goldIndex]);
  }

private:
  stk::mesh::NgpMesh ngpMesh;
  stk::mesh::NgpField<int> ngpVariableLengthField;
  NumScalarsPerEntity goldVariableLengthExtent0PerEntity;
  NumScalarsPerEntity goldVariableLengthExtent1PerEntity;
  stk::mesh::NgpField<int> ngpPartiallyAllocatedField;
  NumScalarsPerEntity goldPartiallyAllocatedExtent0PerEntity;
  NumScalarsPerEntity goldPartiallyAllocatedExtent1PerEntity;
};

void test_num_scalars_per_entity(stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & variableLengthField)
{
  const unsigned numElements = 2;
  check_expected_num_elements(bulk, numElements);

  NumScalarsPerEntity goldNumScalarsPerEntity;
  fill_gold_num_scalars_per_entity(bulk, variableLengthField, numElements, goldNumScalarsPerEntity);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int> ngpVariableLengthField = stk::mesh::get_updated_ngp_field<int>(variableLengthField);

  CheckNumScalarsPerEntity checkNumScalarsPerEntity(ngpMesh, ngpVariableLengthField, goldNumScalarsPerEntity);
  stk::mesh::for_each_entity_run(
        ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), checkNumScalarsPerEntity);
}

void test_extents_per_entity(stk::mesh::BulkData & bulk,
                             const stk::mesh::FieldBase & variableLengthField,
                             const stk::mesh::FieldBase & partiallyAllocatedField)
{
  const unsigned numElements = 2;
  check_expected_num_elements(bulk, numElements);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int> ngpVariableLengthField = stk::mesh::get_updated_ngp_field<int>(variableLengthField);
  stk::mesh::NgpField<int> ngpPartiallyAllocatedField = stk::mesh::get_updated_ngp_field<int>(partiallyAllocatedField);

  NumScalarsPerEntity goldVariableLengthExtent0PerEntity(3, 3);
  NumScalarsPerEntity goldVariableLengthExtent1PerEntity(4, 5);
  NumScalarsPerEntity goldPartiallyAllocatedExtent0PerEntity(2, 0);
  NumScalarsPerEntity goldPartiallyAllocatedExtent1PerEntity(1, 0);

  CheckExtentsPerEntity checkExtentsPerEntity(ngpMesh,
                                              ngpVariableLengthField, goldVariableLengthExtent0PerEntity, goldVariableLengthExtent1PerEntity,
                                              ngpPartiallyAllocatedField, goldPartiallyAllocatedExtent0PerEntity, goldPartiallyAllocatedExtent1PerEntity);

  stk::mesh::for_each_entity_run(
        ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), checkExtentsPerEntity);
}

TEST_F(NgpFieldFixture, NumScalarsPerEntityOnDevice)
{
  if (get_parallel_size() != 1) return;

  setup_two_field_two_element_mesh();

  stk::mesh::Field<int>& variableLengthField = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "variableLengthField");

  test_num_scalars_per_entity(get_bulk(), variableLengthField);
}

TEST_F(NgpFieldFixture, ExtentsPerEntityOnDevice)
{
  if (get_parallel_size() != 1) return;

  setup_two_variable_fields_two_element_mesh();

  stk::mesh::Field<int>& variableLengthField = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "variableLengthField");
  stk::mesh::Field<int>& partiallyAllocatedField = *get_meta().get_field<int>(stk::topology::ELEM_RANK, "partiallyAllocatedField");

  test_extents_per_entity(get_bulk(), variableLengthField, partiallyAllocatedField);
}

TEST_F(NgpFieldFixture, TestAriaAlgorithm)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Field<int> & f1 = create_field<int>(stk::topology::ELEM_RANK, "f1");
  stk::mesh::Field<int> & f2 = create_field<int>(stk::topology::ELEM_RANK, "f2");

  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpMesh ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

  int multiplier = 2;
  modify_field_on_device(get_bulk(), f1, multiplier);

  move_data_between_fields_on_host(get_bulk(), f1, f2);

  test_field_values_on_device(get_bulk(), f1, multiplier);
  test_field_values_on_device(get_bulk(), f2, multiplier);

  test_field_values_on_host(get_bulk(), f1, multiplier);
  test_field_values_on_host(get_bulk(), f2, multiplier);
}

TEST_F(NgpFieldFixture, GetNgpField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int> & ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  EXPECT_EQ(ngpIntField.get_ordinal(), stkIntField.mesh_meta_data_ordinal());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::DeviceField<int>&>::value));
#else
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::HostField<int>&>::value));
#endif
}

TEST_F(NgpFieldFixture, GetMultipleNgpFields)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField1 = create_field<int>(stk::topology::ELEM_RANK, "intField1");
  stk::mesh::Field<int> & stkIntField2 = create_field<int>(stk::topology::ELEM_RANK, "intField2");
  stk::mesh::Field<double> & stkDoubleField1 = create_field<double>(stk::topology::ELEM_RANK, "doubleField1");
  stk::mesh::Field<double> & stkDoubleField2 = create_field<double>(stk::topology::ELEM_RANK, "doubleField2");

  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int> & ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(stkIntField1);
  stk::mesh::NgpField<int> & ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(stkIntField2);
  stk::mesh::NgpField<double> & ngpDoubleField1 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField1);
  stk::mesh::NgpField<double> & ngpDoubleField2 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField2);

  EXPECT_EQ(ngpIntField1.get_ordinal(), stkIntField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpIntField2.get_ordinal(), stkIntField2.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField1.get_ordinal(), stkDoubleField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField2.get_ordinal(), stkDoubleField2.mesh_meta_data_ordinal());
}

TEST_F(NgpFieldFixture, ModifyAndSync)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 2;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 3;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_device(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);

  sync_field_to_host(stkIntField);
}

TEST_F(NgpFieldFixture, ConsistentNeedToSyncForAllCopies)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int> copyNgpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  copyNgpField.sync_to_device();
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  copyNgpField.sync_to_host();
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

  size_t expectedSyncsToDevice = 2;
  size_t expectedSyncsToHost = 1;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  EXPECT_EQ(stkIntField.num_syncs_to_device(), copyNgpField.num_syncs_to_device());
  EXPECT_EQ(stkIntField.num_syncs_to_host(), copyNgpField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, ConsistentNeedToSyncForAllCopiesNoModifyCall)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int> copyNgpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 1;
  check_field_on_host(get_bulk(), stkIntField, multiplier);
  copyNgpField.sync_to_device();
  copyNgpField.sync_to_host();
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  size_t expectedSyncsToDevice = 1;
  size_t expectedSyncsToHost = 0;

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  EXPECT_EQ(stkIntField.num_syncs_to_device(), copyNgpField.num_syncs_to_device());
  EXPECT_EQ(stkIntField.num_syncs_to_host(), copyNgpField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, RequireSyncBetweenModifyOnHostAndDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  stk::io::fill_mesh("generated:1x1x2", get_bulk());

  initialize_ngp_field(stkIntField);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
  EXPECT_THROW(ngpField.modify_on_device(), std::logic_error);
#endif

  stkIntField.clear_sync_state();
  ngpField.modify_on_device();
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
  EXPECT_THROW(ngpField.modify_on_host(), std::logic_error);
#endif
}

TEST_F(NgpFieldFixture, ClearSyncStateAfterModifyOnDevice)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  ngpField.modify_on_device();

  EXPECT_TRUE(ngpField.need_sync_to_host());

  stkIntField.clear_sync_state();
  EXPECT_FALSE(ngpField.need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearHostSyncState)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  ngpField.modify_on_host();

  EXPECT_TRUE(ngpField.need_sync_to_device());

  ngpField.clear_host_sync_state();

  EXPECT_FALSE(ngpField.need_sync_to_device());
}

TEST_F(NgpFieldFixture, ClearDeviceSyncState)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  ngpField.modify_on_device();

  EXPECT_TRUE(ngpField.need_sync_to_host());

  ngpField.clear_device_sync_state();

  EXPECT_FALSE(ngpField.need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearHostSyncState_doesntClearDeviceMod)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  ngpField.modify_on_device();

  EXPECT_TRUE(ngpField.need_sync_to_host());

  ngpField.clear_host_sync_state();

  EXPECT_TRUE(ngpField.need_sync_to_host());
}

TEST_F(NgpFieldFixture, ClearDeviceSyncState_doesntClearHostMod)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(stkIntField);

  ngpField.modify_on_host();

  EXPECT_TRUE(ngpField.need_sync_to_device());

  ngpField.clear_device_sync_state();

  EXPECT_TRUE(ngpField.need_sync_to_device());
}

TEST_F(NgpFieldFixture, updateBucketPtrView)
{
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
  stk::mesh::Field<int>& field = create_field<int>(stk::topology::ELEM_RANK, "field");
  get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);
  ngp_unit_test_utils::setup_mesh_2hex_2block(get_bulk(), bucketCapacity);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);

  add_3rd_element_to_2hex_3block_mesh();

  ngpField.update_bucket_pointer_view();
}

TEST_F(NgpFieldFixture, LateFieldUsage)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  initialize_ngp_field(stkIntField);

  get_meta().enable_late_fields();
  stk::mesh::Field<int> & stkLateIntField = create_field<int>(stk::topology::ELEM_RANK, "lateIntField");

  initialize_ngp_field(stkLateIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  modify_field_on_host(get_bulk(), stkLateIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkLateIntField, multiplier);

  sync_field_to_device(stkIntField);
  sync_field_to_device(stkLateIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);
  modify_field_on_device(get_bulk(), stkLateIntField, multiplier);

  sync_field_to_host(stkIntField);
  sync_field_to_host(stkLateIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);
  check_field_on_host(get_bulk(), stkLateIntField, multiplier*multiplier);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_1|block_3|        |block_1|block_3|block_3|
//   -------------------------        -------------------------
//
TEST_F(OptimizedNgpFieldFixture, ChangeBucketContentsUsingLayoutModificationWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_change_bucket_content_by_mesh_modification(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, ChangeBucketContentsUsingLayoutModificationWithThreeComponents)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 3;

  run_change_bucket_content_by_mesh_modification(numComponents);
}

void construct_ngp_field(stk::mesh::Field<int>& stkIntField)
{
  EXPECT_NO_THROW(stk::mesh::get_updated_ngp_field<int>(stkIntField));
}

TEST_F(OptimizedNgpFieldFixture, CreateConsecutiveNgpFields)
{
  if (get_parallel_size() != 1) return;

  const unsigned bucketCapacity = 2;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

  unsigned numComponents = 1;

  stk::mesh::Field<int>& stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField", numComponents);

  setup_3hex_2block_mesh_with_field(bucketCapacity, stkIntField);

  construct_ngp_field(stkIntField);
  construct_ngp_field(stkIntField);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_2|block_1|
//   -------------------------        -------------------------
//

TEST_F(OptimizedNgpFieldFixture, AddBucketInMiddleWithSingleComponent_FieldReference)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_bucket_in_middle_internal(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, AddBucketInMiddleWithSingleComponent_FieldCopy)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_bucket_in_middle_copy(numComponents);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |  3,4  |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//

TEST_F(OptimizedNgpFieldFixture, AddNewElementWithBucketCapacity2)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;
  unsigned bucketCapacity = 2;

  run_add_new_element(numComponents, bucketCapacity);
}

TEST_F(OptimizedNgpFieldFixture, AddNewElementAndModifyFromDevice)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_new_element_and_modify_from_device(numComponents);
}

TEST_F(OptimizedNgpFieldFixture, AddNewElementWithBucketCapacity1)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;
  unsigned bucketCapacity = 1;

  run_add_new_element(numComponents, bucketCapacity);
}

//   -------------------------        -------------------------
//   |       |       |       |        |       |       |       |
//   |   1   |   2   |   3   |  ===>  |   1   |   2   |   3   |
//   |block_1|block_2|block_3|        |block_1|block_1|block_3|
//   -------------------------        -------------------------
//
TEST_F(OptimizedNgpFieldFixture, DeleteBucketInMiddleWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_delete_bucket_in_middle(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   1   |   2   |
//   |block_1|block_2|        |block_3|block_1|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   1   |   2   |
//   |block_1|block_2|        |block_2|block_3|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent2)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket2(numComponents);
}

//   -----------------        -----------------
//   |       |       |        |       |       |
//   |   1   |   2   |  ===>  |   2   |   3   |
//   |block_1|block_2|        |block_2|block_3|
//   -----------------        -----------------
//

TEST_F(OptimizedNgpFieldFixture, AddAndDeleteBucketWithSingleComponent3)
{
  if (get_parallel_size() != 1) return;

  unsigned numComponents = 1;

  run_add_and_delete_bucket3(numComponents);
}

class ModifyBySelectorFixture : public OptimizedNgpFieldFixture
{
public:
  ModifyBySelectorFixture()
    : entityIdMultiplier(10),
      m_numComponents(2),
      m_bucketCapacity(1)
  { }

  stk::mesh::Field<int>& setup_field_on_mesh_two_buckets_two_components()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, m_bucketCapacity, m_bucketCapacity);
    stk::mesh::Field<int>& stkField = create_field<int>(stk::topology::ELEM_RANK, "intField", m_numComponents);
    ngp_unit_test_utils::setup_mesh_2hex_2block(get_bulk(), m_bucketCapacity);
    return stkField;
  }

  stk::mesh::Field<int>& setup_partial_field_on_mesh_three_buckets_two_components(const std::vector<std::string>& partNames)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, m_bucketCapacity, m_bucketCapacity);
    unsigned numStates = 1;
    const std::vector<int> init(m_numComponents, 1);
    stk::mesh::Field<int>& stkField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "intField", numStates);

    for (const std::string& name : partNames) {
      stk::mesh::Part& part = get_meta().declare_part_with_topology(name, stk::topology::HEX_8);
      stk::mesh::put_field_on_mesh(stkField, part, m_numComponents, init.data());
    }

    ngp_unit_test_utils::setup_mesh_3hex_3block(get_bulk(), m_bucketCapacity);
    return stkField;
  }

  void fill_field_on_device(stk::mesh::Field<int>& stkField)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    set_element_field_data_on_device(ngpMesh, stkField, get_meta().universal_part(), entityIdMultiplier);
  }

protected:
  unsigned entityIdMultiplier;

private:
  unsigned m_numComponents;
  unsigned m_bucketCapacity;
};

TEST_F(ModifyBySelectorFixture, hostToDevice_dontSpecifySelector_byReference)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Field<int>& stkField = setup_field_on_mesh_two_buckets_two_components();
  set_element_field_data(stkField, get_meta().universal_part(), entityIdMultiplier);

  stk::mesh::NgpField<int>& ngpFieldByRef = stk::mesh::get_updated_ngp_field<int>(stkField);
  ngpFieldByRef.modify_on_host();
  ngpFieldByRef.sync_to_device();

  check_field_data_on_device<int>(ngpFieldByRef, stkField);
}

TEST_F(ModifyBySelectorFixture, hostToDevice_partialField_byReference)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Field<int>& stkField = setup_partial_field_on_mesh_three_buckets_two_components({"block_2", "block_3"});
  set_element_field_data(stkField, get_meta().universal_part(), entityIdMultiplier);

  stk::mesh::NgpField<int>& ngpFieldByRef = stk::mesh::get_updated_ngp_field<int>(stkField);
  ngpFieldByRef.modify_on_host(get_meta().universal_part());
  ngpFieldByRef.sync_to_device();

  check_field_data_on_device<int>(ngpFieldByRef, stkField);
}

TEST(DeviceField, checkSizeof)
{
  size_t expectedNumBytes = 384;
  std::cout << "sizeof(stk::mesh::DeviceField<double>): " << sizeof(stk::mesh::DeviceField<double>) << std::endl;
  EXPECT_TRUE(sizeof(stk::mesh::DeviceField<double>) <= expectedNumBytes);
}

TEST(DeviceBucket, checkSizeof)
{
#ifndef STK_HIDE_DEPRECATED_CODE  // Delete after 2024/06/26
  size_t expectedNumBytes = 176;
#else
  size_t expectedNumBytes = 152;  // Value after removing DeviceBucket::m_hostEntities
#endif
  std::cout << "sizeof(stk::mesh::DeviceBucket): " << sizeof(stk::mesh::DeviceBucket) << std::endl;
  EXPECT_TRUE(sizeof(stk::mesh::DeviceBucket) <= expectedNumBytes);
}


enum PartIds : int {
  part_1  = 1,
  part_2  = 2,
  part_3  = 3,
  part_4  = 4,
  part_5  = 5,
  part_6  = 6,
  part_7  = 7,
  part_8  = 8,
  part_9  = 9,
  part_10 = 10,
  part_11 = 11,
  part_12 = 12,
  part_13 = 13,
  part_14 = 14,
  part_15 = 15
};

struct NodeIdPartId {
  int     nodeId;
  PartIds partId;
};

class NgpFieldUpdate : public ::ngp_testing::Test
{
public:
  NgpFieldUpdate()
    : m_field(nullptr)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    m_meta = builder.create_meta_data();
  }

  void create_mesh_with_parts_and_nodes(std::vector<PartIds> partIds, std::vector<NodeIdPartId> nodes)
  {
    for (PartIds partId : partIds) {
      stk::mesh::Part & part = m_meta->declare_part_with_topology("part_" + std::to_string(partId),
                                                                   stk::topology::NODE);
      m_parts[partId] = &part;
    }

    m_field = &m_meta->declare_field<int>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*m_field, m_meta->universal_part(), nullptr);

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    m_bulk = builder.create(m_meta);

    m_bulk->modification_begin();
    for (NodeIdPartId node : nodes) {
      const stk::mesh::Entity newNode = m_bulk->declare_node(node.nodeId,
                                                             stk::mesh::PartVector{m_parts.at(node.partId)});
      *stk::mesh::field_data(*m_field, newNode) = node.nodeId;
    }
    m_bulk->modification_end();

    stk::mesh::get_updated_ngp_mesh(*m_bulk);
    stk::mesh::get_updated_ngp_field<int>(*m_field);
  }

  void add_node(NodeIdPartId newNodeInfo)
  {
    m_bulk->modification_begin();
    const stk::mesh::Entity node = m_bulk->declare_node(newNodeInfo.nodeId,
                                                        stk::mesh::PartVector{m_parts.at(newNodeInfo.partId)});
    *stk::mesh::field_data(*m_field, node) = newNodeInfo.nodeId;
    m_bulk->modification_end();
  }

  void remove_node(int removedNodeId)
  {
    m_bulk->modification_begin();
    const stk::mesh::Entity node = m_bulk->get_entity(stk::topology::NODE_RANK, removedNodeId);
    m_bulk->destroy_entity(node);
    m_bulk->modification_end();
  }

  void check_field_values()
  {
    const stk::mesh::BucketVector& buckets = m_bulk->buckets(stk::topology::NODE_RANK);
    unsigned nodeIdx = 0;
    for (const stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const int* fieldValue = stk::mesh::field_data(*m_field, node);
        const int expectedValue = m_bulk->identifier(node);
        EXPECT_EQ(*fieldValue, expectedValue);
      }
    }

    const unsigned numNodes = stk::mesh::count_entities(*m_bulk, stk::topology::NODE_RANK, m_meta->universal_part());
    Kokkos::View<int*> deviceValues("deviceValues", numNodes);
    Kokkos::View<int*>::HostMirror hostValuesFromDevice = Kokkos::create_mirror_view(deviceValues);

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulk);
    stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(*m_field);

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(size_t /*index*/) {
        unsigned nodeGlobalIndex = 0;
        const unsigned numBuckets = ngpMesh.num_buckets(stk::topology::NODE_RANK);
        for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
          const stk::mesh::NgpMesh::BucketType & deviceBucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, bucketId);
          const unsigned numNodesInBucket = deviceBucket.size();
          for (unsigned nodeOrdinal = 0; nodeOrdinal < numNodesInBucket; ++nodeOrdinal) {
            const stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(deviceBucket[nodeOrdinal]);
            deviceValues[nodeGlobalIndex++] = ngpField.get(nodeIndex, 0);
          }
        }
      });

    Kokkos::deep_copy(hostValuesFromDevice, deviceValues);

    nodeIdx = 0;
    for (const stk::mesh::Bucket * bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const int expectedValue = m_bulk->identifier(node);
        EXPECT_EQ(hostValuesFromDevice[nodeIdx++], expectedValue);
      }
    }
  }

protected:
  std::shared_ptr<stk::mesh::BulkData> m_bulk;
  std::shared_ptr<stk::mesh::MetaData> m_meta;
  std::map<PartIds, stk::mesh::Part*> m_parts;
  stk::mesh::Field<int> * m_field;
};

//   o = Node that persists through modification
//   * = New node
//   x = Node that will be deleted

// |   2   |   3   |       |   1   |   2   |   3   |
// |   o   |   o   |  ==>  |   *   |   o   |   o   |
// |part_2 |part_3 |       |part_1 |part_2 |part_3 |
//
TEST_F(NgpFieldUpdate, AddBucketAtBeginning_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{2, part_2}, {3, part_3}});
  add_node({1, part_1});

  check_field_values();
}

// |   1   |   3   |       |   1   |   2   |   3   |
// |   o   |   o   |  ==>  |   o   |   *   |   o   |
// |part_1 |part_3 |       |part_1 |part_2 |part_3 |
//
TEST_F(NgpFieldUpdate, AddBucketInMiddle_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{1, part_1}, {3, part_3}});
  add_node({2, part_2});

  check_field_values();
}

// |   1   |   2   |       |   1   |   2   |   3   |
// |   o   |   o   |  ==>  |   o   |   o   |   *   |
// |part_1 |part_2 |       |part_1 |part_2 |part_3 |
//
TEST_F(NgpFieldUpdate, AddBucketAtEnd_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{1, part_1}, {2, part_2}});
  add_node({3, part_3});

  check_field_values();
}

// |   1   |   2   |   3   |       |   2   |   3   |
// |   x   |   o   |   o   |  ==>  |   o   |   o   |
// |part_1 |part_2 |part_3 |       |part_2 |part_3 |
//
TEST_F(NgpFieldUpdate, RemoveBucketAtFront_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{1, part_1}, {2, part_2}, {3, part_3}});
  remove_node(1);

  check_field_values();
}

// |   1   |   2   |   3   |       |   1   |   3   |
// |   o   |   x   |   o   |  ==>  |   o   |   o   |
// |part_1 |part_2 |part_3 |       |part_1 |part_3 |
//
TEST_F(NgpFieldUpdate, RemoveBucketInMiddle_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{1, part_1}, {2, part_2}, {3, part_3}});
  remove_node(2);

  check_field_values();
}

// |   1   |   2   |   3   |       |   1   |   2   |
// |   o   |   o   |   x   |  ==>  |   o   |   o   |
// |part_1 |part_2 |part_3 |       |part_1 |part_2 |
//
TEST_F(NgpFieldUpdate, RemoveBucketAtEnd_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3},
                                   {{1, part_1}, {2, part_2}, {3, part_3}});
  remove_node(3);

  check_field_values();
}

// |   1   |   3   |   4   |       | 1   2 |   3   |   4   |   5   |
// |   o   |   o   |   o   |  ==>  | o   * |   o   |   o   |   *   |
// |part_1 |part_3 |part_4 |       |part_1 |part_3 |part_4 |part_5 |
//
TEST_F(NgpFieldUpdate, ModifyBucketAtBeginning_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_3, part_4, part_5},
                                   {{1, part_1}, {3, part_3}, {4, part_4}});
  add_node({2, part_1});
  add_node({5, part_5});

  check_field_values();
}

// |   1   |   2   |   4   |       |   1   | 2   3 |   4   |   5   |
// |   o   |   o   |   o   |  ==>  |   o   | o   * |   o   |   *   |
// |part_1 |part_2 |part_4 |       |part_1 |part_2 |part_4 |part_5 |
//
TEST_F(NgpFieldUpdate, ModifyBucketInMiddle_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_4, part_5},
                                   {{1, part_1}, {2, part_2}, {4, part_4}});
  add_node({3, part_2});
  add_node({5, part_5});

  check_field_values();
}

// |   1   |   3   |   4   |       |   1   |   2   |   3   | 4   5 |
// |   o   |   o   |   o   |  ==>  |   o   |   *   |   o   | o   * |
// |part_1 |part_3 |part_4 |       |part_1 |part_2 |part_3 |part_4 |
//
TEST_F(NgpFieldUpdate, ModifyBucketAtEnd_WhileReallocating)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4},
                                   {{1, part_1}, {3, part_3}, {4, part_4}});
  add_node({2, part_2});
  add_node({5, part_4});

  check_field_values();
}

// |   1   |   2   |   4   |       |   1   |   3   |   4   |
// |   o   |   x   |   o   |  ==>  |   o   |   *   |   o   |
// |part_1 |part_2 |part_4 |       |part_1 |part_3 |part_4 |
//
TEST_F(NgpFieldUpdate, AddAndRemoveBucketSameLocation)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4},
                                   {{1, part_1}, {2, part_2}, {4, part_4}});
  remove_node(2);
  add_node({3, part_3});

  check_field_values();
}

// |   1   |   2   |   3   |   4   |   6   |
// |   x   |   x   |   o   |   o   |   o   |  ==>
// |part_1 |part_2 |part_3 |part_4 |part_6 |
//                 ____/           ____/
//            ____/           ____/
//       ____/           ____/
//      /               /
//     V               V
// |   3   | 4   5 |   6   |   7   |   8   |
// |   o   | o   * |   o   |   *   |   *   |
// |part_3 |part_4 |part_6 |part_7 |part_8 |
//
TEST_F(NgpFieldUpdate, OverlapBackwardMovedRangesDueToModifiedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_6, part_7, part_8},
                                   {{1, part_1}, {2, part_2}, {3, part_3}, {4, part_4}, {6, part_6}});
  remove_node(1);
  remove_node(2);
  add_node({5, part_4});
  add_node({7, part_7});
  add_node({8, part_8});

  check_field_values();
}

// |   1   |   2   |   3   |   4   |   5   |
// |   o   |   x   |   o   |   x   |   o   |  ==>
// |part_1 |part_2 |part_3 |part_4 |part_5 |
//     |             __/           ____/
//     |           _/         ____/
//     |         _/      ____/
//     |        /       /
//     V       V       V
// |   1   |   3   |   5   |   6   |   7   |
// |   o   |   o   |   o   |   *   |   *   |
// |part_1 |part_3 |part_5 |part_6 |part_7 |
//
TEST_F(NgpFieldUpdate, OverlapBackwardMovedRangesDueToDeletedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7},
                                   {{1, part_1}, {2, part_2}, {3, part_3}, {4, part_4}, {5, part_5}});
  remove_node(2);
  remove_node(4);
  add_node({6, part_6});
  add_node({7, part_7});

  check_field_values();
}

// |   1   |   2   |   3   |   4   |   5   |   6   |   8   |   9   |  11   |
// |   x   |   x   |   x   |   x   |   o   |   o   |   o   |   o   |   o   |  ==>
// |part_1 |part_2 |part_3 |part_4 |part_5 |part_6 |part_8 |part_9 |part_11|
//                           __________/     __________/      _________/
//                 _________/      _________/       _________/
//       _________/      _________/       _________/
//      /               /                /
//     V               V                V
// |   5   | 6   7 |   8   | 9  10 |  11   |  12   |  13   |  14   |  15   |
// |   o   | o   * |   o   | o   * |   o   |   *   |   *   |   *   |   *   |
// |part_5 |part_6 |part_8 |part_9 |part_11|part_12|part_13|part_14|part_15|
//
TEST_F(NgpFieldUpdate, DoubleOverlapBackwardMovedRangesDueToModifiedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_8, part_9, part_11,
                                    part_12, part_13, part_14, part_15},
                                   {{1, part_1}, {2, part_2}, {3, part_3}, {4, part_4}, {5, part_5}, {6, part_6},
                                    {8, part_8}, {9, part_9}, {11, part_11}});
  remove_node(1);
  remove_node(2);
  remove_node(3);
  remove_node(4);
  add_node({7, part_6});
  add_node({10, part_9});
  add_node({12, part_12});
  add_node({13, part_13});
  add_node({14, part_14});
  add_node({15, part_15});

  check_field_values();
}

// |   1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   |
// |   x   |   x   |   o   |   o   |   x   |   o   |   x   |   o   |  ==>
// |part_1 |part_2 |part_3 |part_4 |part_5 |part_6 |part_7 |part_8 |
//                 ____/   ____/        _______/     __________/
//            ____/   ____/     _______/   _________/
//       ____/   ____/   ______/ _________/
//      /       /       /       /
//     V       V       V       V
// |   3   |   4   |   6   |   8   |   9   |  10   |  11   |  12   |
// |   o   |   o   |   o   |   o   |   *   |   *   |   *   |   *   |
// |part_3 |part_4 |part_6 |part_8 |part_9 |part_10|part_11|part_12|
//
TEST_F(NgpFieldUpdate, DoubleOverlapBackwardMovedRangesDueToDeletedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7, part_8, part_9,
                                    part_10, part_11, part_12},
                                   {{1, part_1}, {2, part_2}, {3, part_3}, {4, part_4}, {5, part_5}, {6, part_6},
                                    {7, part_7}, {8, part_8}});
  remove_node(1);
  remove_node(2);
  remove_node(5);
  remove_node(7);
  add_node({9, part_9});
  add_node({10, part_10});
  add_node({11, part_11});
  add_node({12, part_12});

  check_field_values();
}

// |   3   |   4   |   6   |   7   |   8   |
// |   o   |   o   |   o   |   x   |   x   |  ==>
// |part_3 |part_4 |part_6 |part_7 |part_8 |
//     \____           \____
//          \____           \____
//               \____           \____
//                    \               \    .
//                     V               V
// |   1   |   2   |   3   | 4   5 |   6   |
// |   *   |   *   |   o   | o   * |   o   |
// |part_1 |part_2 |part_3 |part_4 |part_6 |
//
TEST_F(NgpFieldUpdate, OverlapForwardMovedRangesDueToModifiedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_6, part_7, part_8},
                                   {{3, part_3}, {4, part_4}, {6, part_6}, {7, part_7}, {8, part_8}});
  remove_node(7);
  remove_node(8);
  add_node({1, part_1});
  add_node({2, part_2});
  add_node({5, part_4});

  check_field_values();
}

// |   3   |   4   |   5   |   6   |   7   |
// |   o   |   x   |   o   |   x   |   o   |  ==>
// |part_3 |part_4 |part_5 |part_6 |part_7 |
//     \____           \__             |
//          \____         \_           |
//               \____      \_         |
//                    \       \        |
//                     V       V       V
// |   1   |   2   |   3   |   5   |   7   |
// |   *   |   *   |   o   |   o   |   o   |
// |part_1 |part_2 |part_3 |part_5 |part_7 |
//
TEST_F(NgpFieldUpdate, OverlapForwardMovedRangesDueToDeletedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7},
                                   {{3, part_3}, {4, part_4}, {5, part_5}, {6, part_6}, {7, part_7}});
  remove_node(4);
  remove_node(6);
  add_node({1, part_1});
  add_node({2, part_2});

  check_field_values();
}

// |   5   |   6   |   8   |   9   |  11   |  12   |  13   |  14   |  15   |
// |   o   |   o   |   o   |   o   |   o   |   x   |   x   |   x   |   x   |  ==>
// |part_5 |part_6 |part_8 |part_9 |part_11|part_12|part_13|part_14|part_15|
//     \__________     \__________     \__________
//                \_________      \_________      \_________
//                          \_________      \_________      \_________
//                                    \               \               \    .
//                                     V               V               V
// |   1   |   2   |   3   |   4   |   5   | 6   7 |   8   | 9  10 |  11   |
// |   *   |   *   |   *   |   *   |   o   | o   * |   o   | o   * |   o   |
// |part_1 |part_2 |part_3 |part_4 |part_5 |part_6 |part_8 |part_9 |part_11|
//
TEST_F(NgpFieldUpdate, DoubleOverlapForwardMovedRangesDueToModifiedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_8, part_9, part_11,
                                    part_12, part_13, part_14, part_15},
                                   {{5, part_5}, {6, part_6}, {8, part_8}, {9, part_9}, {11, part_11}, {12, part_12},
                                    {13, part_13}, {14, part_14}, {15, part_15}});
  remove_node(12);
  remove_node(13);
  remove_node(14);
  remove_node(15);
  add_node({1, part_1});
  add_node({2, part_2});
  add_node({3, part_3});
  add_node({4, part_4});
  add_node({7, part_6});
  add_node({10, part_9});

  check_field_values();
}

// |   5   |   6   |   7   |   8   |   9   |  10   |  11   |  12   |
// |   o   |   x   |   o   |   x   |   o   |   o   |   x   |   x   |  ==>
// |part_5 |part_6 |part_7 |part_8 |part_9 |part_10|part_11|part_12|
//     \__________     \_______        \____   \____
//                \_________   \_______     \____   \____
//                          \_________ \______   \____   \____
//                                    \       \       \       \    .
//                                     V       V       V       V
// |   1   |   2   |   3   |   4   |   5   |   7   |   9   |  10   |
// |   *   |   *   |   *   |   *   |   o   |   o   |   o   |   o   |
// |part_1 |part_2 |part_3 |part_4 |part_5 |part_7 |part_9 |part_10|
//
TEST_F(NgpFieldUpdate, DoubleOverlapForwardMovedRangesDueToDeletedBucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7, part_8, part_9,
                                    part_10, part_11, part_12},
                                   {{5, part_5}, {6, part_6}, {7, part_7}, {8, part_8}, {9, part_9}, {10, part_10},
                                    {11, part_11}, {12, part_12}});
  remove_node(6);
  remove_node(8);
  remove_node(11);
  remove_node(12);
  add_node({1, part_1});
  add_node({2, part_2});
  add_node({3, part_3});
  add_node({4, part_4});

  check_field_values();
}

// |   2   |   3   |   4   |   5   |   8   |   9   |
// |   o   |   x   |   x   |   o   |   o   |   x   |  ==>
// |part_2 |part_3 |part_4 |part_5 |part_8 |part_9 |
//     \__                   __/       \__
//        \_               _/             \_
//          \_           _/                 \_
//            \         /                     \    .
//             V       V                       V
// |   1   |   2   |   5   |   6   |   7   |   8   |
// |   *   |   o   |   o   |   *   |   *   |   o   |
// |part_1 |part_2 |part_5 |part_6 |part_7 |part_8 |
//
TEST_F(NgpFieldUpdate, MoveForwardBackwardForward)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7, part_8, part_9},
                                   {{2, part_2}, {3, part_3}, {4, part_4}, {5, part_5}, {8, part_8}, {9, part_9}});
  remove_node(3);
  remove_node(4);
  remove_node(9);
  add_node({1, part_1});
  add_node({6, part_6});
  add_node({7, part_7});

  check_field_values();
}

// |   1   |   2   |   5   |   6   |   7   |   8   |
// |   x   |   o   |   o   |   x   |   x   |   o   |  ==>
// |part_1 |part_2 |part_5 |part_6 |part_7 |part_8 |
//           __/       \__                   __/
//         _/             \_               _/
//       _/                 \_           _/
//      /                     \         /
//     V                       V       V
// |   2   |   3   |   4   |   5   |   8   |   9   |
// |   o   |   *   |   *   |   o   |   o   |   *   |
// |part_2 |part_3 |part_4 |part_5 |part_8 |part_9 |
//
TEST_F(NgpFieldUpdate, MoveBackwardForwardBackward)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  create_mesh_with_parts_and_nodes({part_1, part_2, part_3, part_4, part_5, part_6, part_7, part_8, part_9},
                                   {{1, part_1}, {2, part_2}, {5, part_5}, {6, part_6}, {7, part_7}, {8, part_8}});
  remove_node(1);
  remove_node(6);
  remove_node(7);
  add_node({3, part_3});
  add_node({4, part_4});
  add_node({9, part_9});

  check_field_values();
}

}
