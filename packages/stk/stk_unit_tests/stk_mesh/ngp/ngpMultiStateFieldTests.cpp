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

#include <stk_util/stk_config.h>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/TestHexFixture.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/NgpFieldBLAS.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>

namespace {

class ClassWithNgpField
{
public:
  KOKKOS_FUNCTION
  ClassWithNgpField(const stk::mesh::NgpField<double>& ngpField)
    : m_ngpField(ngpField)
  {}

  KOKKOS_FUNCTION
  unsigned num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const
  { return m_ngpField.get_num_components_per_entity(entity); }

  KOKKOS_FUNCTION
  double access_field_data(const stk::mesh::FastMeshIndex& entity, unsigned component) const
  { return m_ngpField(entity, component); }

private:
  stk::mesh::NgpField<double> m_ngpField;
};

#ifdef STK_ENABLE_GPU
#define MY_LAMBDA KOKKOS_LAMBDA
#else
#define MY_LAMBDA [&]
#endif

ClassWithNgpField* create_class_on_device(const stk::mesh::NgpField<double>& ngpField)
{
  ClassWithNgpField* devicePtr = static_cast<ClassWithNgpField*>(
        Kokkos::kokkos_malloc<stk::ngp::MemSpace>("device class memory", sizeof(ClassWithNgpField)));
  Kokkos::parallel_for("construct class on device", stk::ngp::DeviceRangePolicy(0, 1),
                       MY_LAMBDA(const int /*i*/) { new (devicePtr) ClassWithNgpField(ngpField); }
                       );
  Kokkos::fence();
  return devicePtr;
}

void delete_class_on_device(ClassWithNgpField* devicePtr)
{
  Kokkos::parallel_for("device_destruct", stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const int /*i*/) { devicePtr->~ClassWithNgpField(); }
                       );
  Kokkos::fence();
  Kokkos::kokkos_free<stk::ngp::MemSpace>(static_cast<void*>(devicePtr));
}

class LegacyNgpMultiStateFieldTest : public stk::mesh::fixtures::TestHexFixture
{
public:

  template <typename T>
  stk::mesh::Field<T> & create_multistate_field(stk::topology::rank_t rank, const std::string & name, unsigned numStates)
  {
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
    return field;
  }

  void setup_multistate_field(int numStates)
  {
    stk::mesh::FieldBase& field = create_multistate_field<double>(stk::topology::NODE_RANK, "myField", numStates);
    for (int state = 0; state < numStates; ++state) {
      m_fieldStates.push_back(field.field_state(static_cast<stk::mesh::FieldState>(state)));
    }

    get_meta().declare_part("test_part");
  }

  stk::mesh::Field<double>& get_field(stk::mesh::FieldState state) {
    return *static_cast<stk::mesh::Field<double>*>(m_fieldStates[state]);
  }

  template <typename ValueType>
  struct CheckValueUsingNgpField {
    CheckValueUsingNgpField(const stk::mesh::NgpField<ValueType>& _ngpField, ValueType _expectedValue)
      : ngpField(_ngpField), expectedValue(_expectedValue)
    {
    }

    KOKKOS_FUNCTION
    void operator()(const stk::mesh::FastMeshIndex& entity) const
    {
      unsigned numComponents = ngpField.get_num_components_per_entity(entity);
      for (unsigned component = 0; component < numComponents; ++component) {
        NGP_EXPECT_EQ(expectedValue, ngpField(entity, component));
      }
    }

  private:
    stk::mesh::NgpField<ValueType> ngpField;
    ValueType expectedValue;
  };

  template<typename T>
  void check_field_data_value_on_device(const stk::mesh::NgpMesh& ngpMesh,
                                        const stk::mesh::NgpField<T>& ngpField,
                                        T expectedValue)
  {
    stk::mesh::Selector owned = ngpMesh.get_bulk_on_host().mesh_meta_data().locally_owned_part();
    CheckValueUsingNgpField<T> checkValueUsingNgpField(ngpField, expectedValue);
    stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), owned, checkValueUsingNgpField);
  }

  template<typename T>
  void check_field_data_value_on_host(const stk::mesh::BulkData& bulk,
                                      const stk::mesh::Field<T>& hostField,
                                      T expectedValue)
  {
    const stk::mesh::BucketVector& buckets = bulk.get_buckets(hostField.entity_rank(), hostField);

    auto hostFieldData = hostField.data();
    for (const stk::mesh::Bucket* bucket : buckets) {
      for (const stk::mesh::Entity& entity : *bucket) {
        auto value = hostFieldData.entity_values(entity);
        EXPECT_EQ(value(0_comp), expectedValue);
      }
    }
  }

  template <typename ValueType>
  struct CheckValueUsingClass {
    CheckValueUsingClass(const ClassWithNgpField* _deviceClassPointer, ValueType _expectedValue)
      : deviceClassPointer(_deviceClassPointer), expectedValue(_expectedValue)
    {
    }

    KOKKOS_FUNCTION
    void operator()(const stk::mesh::FastMeshIndex& entity) const
    {
      unsigned numComponents = deviceClassPointer->num_components_per_entity(entity);
      for (unsigned component = 0; component < numComponents; ++component) {
        NGP_EXPECT_EQ(expectedValue, deviceClassPointer->access_field_data(entity, component));
      }
    }

  private:
    const ClassWithNgpField* deviceClassPointer;
    ValueType expectedValue;
  };

  template<typename T>
  void check_field_data_value_on_device(const stk::mesh::NgpMesh& ngpMesh,
                                        stk::mesh::EntityRank rank,
                                        const ClassWithNgpField* deviceClassPointer,
                                        T expectedValue)
  {
    stk::mesh::Selector owned = ngpMesh.get_bulk_on_host().mesh_meta_data().locally_owned_part();
    CheckValueUsingClass<T> checkValueUsingClass(deviceClassPointer, expectedValue);
    stk::mesh::for_each_entity_run(ngpMesh, rank, owned, checkValueUsingClass);
  }

  void perform_field_state_rotation()
  {
    stk::mesh::sync_to_host_and_mark_modified(get_meta());
    get_bulk().update_field_data_states();
  }

  void put_all_nodes_in_part(stk::mesh::BulkData& bulk, stk::mesh::Part& part)
  {
    stk::mesh::EntityVector allNodes;
    bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), allNodes);

    bulk.batch_change_entity_parts(allNodes, {&part}, {});
  }

  stk::mesh::FieldVector m_fieldStates;
};


NGP_TEST_F(LegacyNgpMultiStateFieldTest, rotateOnDevice_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 11.1;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueOld);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, rotateOnDevice_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 3;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  stk::mesh::field_fill(valueNp1, get_field(stk::mesh::StateNP1));
  stk::mesh::field_fill(valueN,   get_field(stk::mesh::StateN));
  stk::mesh::field_fill(valueNm1, get_field(stk::mesh::StateNM1));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNp1 = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNP1));
  stk::mesh::NgpField<double>& ngpFieldN   = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateN));
  stk::mesh::NgpField<double>& ngpFieldNm1 = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNM1));

  check_field_data_value_on_device(ngpMesh, ngpFieldNp1, valueNp1);
  check_field_data_value_on_device(ngpMesh, ngpFieldN,   valueN);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm1, valueNm1);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, ngpFieldNp1, valueNm1);
  check_field_data_value_on_device(ngpMesh, ngpFieldN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm1, valueN);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, rotateOnDevice_fourStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 4;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  const double valueNm2 = 44.4;
  stk::mesh::field_fill(valueNp1, get_field(stk::mesh::StateNP1));
  stk::mesh::field_fill(valueN,   get_field(stk::mesh::StateN));
  stk::mesh::field_fill(valueNm1, get_field(stk::mesh::StateNM1));
  stk::mesh::field_fill(valueNm2, get_field(stk::mesh::StateNM2));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNp1 = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNP1));
  stk::mesh::NgpField<double>& ngpFieldN   = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateN));
  stk::mesh::NgpField<double>& ngpFieldNm1 = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNM1));
  stk::mesh::NgpField<double>& ngpFieldNm2 = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNM2));

  check_field_data_value_on_device(ngpMesh, ngpFieldNp1, valueNp1);
  check_field_data_value_on_device(ngpMesh, ngpFieldN,   valueN);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm1, valueNm1);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm2, valueNm2);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, ngpFieldNp1, valueNm2);
  check_field_data_value_on_device(ngpMesh, ngpFieldN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm1, valueN);
  check_field_data_value_on_device(ngpMesh, ngpFieldNm2, valueNm1);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_rotateDeviceStates_syncStatesUnchanged)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueOld);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  EXPECT_FALSE(ngpFieldNew.need_sync_to_host());
  EXPECT_FALSE(ngpFieldOld.need_sync_to_host());
  EXPECT_FALSE(ngpFieldNew.need_sync_to_device());
  EXPECT_FALSE(ngpFieldOld.need_sync_to_device());

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_copyHasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  stk::mesh::NgpField<double> copyOfNgpFieldNew(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueOld);

  perform_field_state_rotation();

  ngpFieldNew.sync_to_device();
  ngpFieldOld.sync_to_device();

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);

  check_field_data_value_on_device(ngpMesh, copyOfNgpFieldNew, valueOld);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_copyHasCorrectDataAfterDeviceStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  stk::mesh::NgpField<double> copyOfNgpFieldNew(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueOld);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);

  check_field_data_value_on_device(ngpMesh, copyOfNgpFieldNew, valueOld);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, persistentDeviceField_hasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  ClassWithNgpField* persistentDeviceClass = create_class_on_device(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  perform_field_state_rotation();

  ngpFieldNew.sync_to_device();
  ngpFieldOld.sync_to_device();

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, persistentDeviceField_hasCorrectDataAfterDeviceStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  ClassWithNgpField* persistentDeviceClass = create_class_on_device(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  ngpFieldNew.sync_to_device();
  ngpFieldOld.sync_to_device();

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, persistentSyncToDeviceCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device());

  ngpFieldNew.modify_on_host();
  ngpFieldNew.sync_to_device();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device()+1);

  perform_field_state_rotation();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device()+1, ngpFieldOld.num_syncs_to_device());
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, persistentSyncToHostCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device());

  ngpFieldOld.modify_on_device();
  ngpFieldOld.sync_to_host();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_host()+1, ngpFieldOld.num_syncs_to_host());

  perform_field_state_rotation();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_host(), ngpFieldOld.num_syncs_to_host()+1);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, persistentModifyOnHostAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  ngpFieldNew.clear_sync_state();
  ngpFieldOld.clear_sync_state();

  ngpFieldNew.modify_on_host();

  EXPECT_TRUE(ngpFieldNew.need_sync_to_device());
  EXPECT_FALSE(ngpFieldOld.need_sync_to_device());

  get_bulk().update_field_data_states();

  EXPECT_FALSE(ngpFieldNew.need_sync_to_device());
  EXPECT_TRUE(ngpFieldOld.need_sync_to_device());
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, createStateOldNgpFieldAfterRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));
  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);

  get_bulk().update_field_data_states();
  ngpFieldNew.modify_on_host();
  ngpFieldNew.sync_to_device();

  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));

  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
}

NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_rotateAllStates_withMeshModification)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));

  stk::mesh::Part& testPart = *get_meta().get_part("test_part");
  stk::mesh::Field<double>& fieldNew = *static_cast<stk::mesh::Field<double>*>(get_meta().get_field(stk::topology::NODE_RANK,
                                                                                                    "myField"));
  stk::mesh::Field<double>& fieldOld = fieldNew.field_of_state(stk::mesh::StateOld);

  stk::mesh::NgpMesh* ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>* ngpFieldOld = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  stk::mesh::NgpField<double>* ngpFieldNew = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));

  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldOld, valueOld);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldNew, valueNew);

  put_all_nodes_in_part(get_bulk(), testPart);
  ngpFieldOld = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  ngpFieldNew = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldOld, valueNew);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldNew, valueOld);
}

// This is a meaningless test without separate device storage, because it explicitly exercises a fundamentally-broken
// swap function that only partially rotates the data.  The HostMesh version of this function does nothing.
// If in a unified-memory build, then leaving the other memory space unrotated and relying on syncs to fill in
// the correct value will not work as expected.  It's not that the other memory space has an unrotated copy of
// the data, but that the other memory space has a pointer back into the wrong state in the opposite memory
// space.
#if defined(STK_USE_DEVICE_MESH) && !defined(STK_UNIFIED_MEMORY)
NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_useAwfulSwapFunction_withSyncsBackToHost)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));

  stk::mesh::Field<double>& fieldNew = *static_cast<stk::mesh::Field<double>*>(get_meta().get_field(stk::topology::NODE_RANK,
                                                                                                    "myField"));
  stk::mesh::Field<double>& fieldOld = fieldNew.field_of_state(stk::mesh::StateOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>* ngpFieldOld = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  stk::mesh::NgpField<double>* ngpFieldNew = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));

  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, *ngpFieldOld, valueOld);
  check_field_data_value_on_device(ngpMesh, *ngpFieldNew, valueNew);

  ngpFieldNew->swap(*ngpFieldOld);  /// This only swaps device data and not pointers back to host data

  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);  // Unrotated values
  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, *ngpFieldOld, valueNew);  // Rotated values
  check_field_data_value_on_device(ngpMesh, *ngpFieldNew, valueOld);

  ngpFieldOld->modify_on_device();
  ngpFieldNew->modify_on_device();
  ngpFieldOld->sync_to_host();  // Uses unrotated pointers to unrotated host data
  ngpFieldNew->sync_to_host();

  check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);  // Overwritten values with correct rotated state
  check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, *ngpFieldOld, valueNew);  // Rotated values
  check_field_data_value_on_device(ngpMesh, *ngpFieldNew, valueOld);
}
#endif

NGP_TEST_F(LegacyNgpMultiStateFieldTest, multistateField_updateUnmodifiedBucket_afterMeshModification)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, get_field(stk::mesh::StateOld));
  stk::mesh::field_fill(valueNew, get_field(stk::mesh::StateNew));

  stk::mesh::Part& testPart = *get_meta().get_part("test_part");
  stk::mesh::Field<double>& fieldNew = *static_cast<stk::mesh::Field<double>*>(get_meta().get_field(stk::topology::NODE_RANK,
                                                                                                    "myField"));
  stk::mesh::Field<double>& fieldOld = fieldNew.field_of_state(stk::mesh::StateOld);

  stk::mesh::NgpMesh* ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>* ngpFieldOld = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateOld));
  stk::mesh::NgpField<double>* ngpFieldNew = &stk::mesh::get_updated_ngp_field<double>(get_field(stk::mesh::StateNew));

  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldOld, valueOld);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldNew, valueNew);

  // New node goes in new Bucket, leaving others alone
  get_bulk().modification_begin();
  {
    const stk::mesh::Entity newNode = get_bulk().declare_node(1000, stk::mesh::PartVector{&testPart});
    auto fieldOldData = fieldOld.data<stk::mesh::ReadWrite>();
    auto fieldNewData = fieldNew.data<stk::mesh::ReadWrite>();
    auto fieldOldValue = fieldOldData.entity_values(newNode);
    auto fieldNewValue = fieldNewData.entity_values(newNode);
    fieldOldValue(0_comp) = valueOld;
    fieldNewValue(0_comp) = valueNew;
  }
  get_bulk().modification_end();

  fieldNew.rotate_multistate_data();  // Don't rotate device data; DeviceFields are updated as part of rotating host pointers

  fieldOld.modify_on_host();
  fieldNew.modify_on_host();
  fieldOld.sync_to_device();
  fieldNew.sync_to_device();

  ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);  // Rotated properly on host
  check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);
  check_field_data_value_on_device(*ngpMesh, *ngpFieldOld, valueNew);  // Rotated data synced properly to device
  check_field_data_value_on_device(*ngpMesh, *ngpFieldNew, valueOld);
}



class ClassWithDeviceFieldData
{
public:
  KOKKOS_FUNCTION
  ClassWithDeviceFieldData(const stk::mesh::ConstFieldData<double, stk::ngp::DeviceSpace>& deviceFieldData)
    : m_deviceFieldData(deviceFieldData)
  {}

  KOKKOS_FUNCTION
  unsigned num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const
  { return m_deviceFieldData.entity_values(entity).num_components(); }

  KOKKOS_FUNCTION
  double access_field_data(const stk::mesh::FastMeshIndex& entity, stk::mesh::ComponentIdx component) const
  { return m_deviceFieldData.entity_values(entity)(component); }

private:
  stk::mesh::ConstFieldData<double, stk::ngp::DeviceSpace> m_deviceFieldData;
};


ClassWithDeviceFieldData* create_class_on_device(
    const stk::mesh::ConstFieldData<double, stk::ngp::DeviceSpace>& deviceFieldData)
{
  ClassWithDeviceFieldData* devicePtr = static_cast<ClassWithDeviceFieldData*>(
        Kokkos::kokkos_malloc<stk::ngp::MemSpace>("device class memory", sizeof(ClassWithDeviceFieldData)));
  Kokkos::parallel_for("construct class on device", stk::ngp::DeviceRangePolicy(0, 1),
                         MY_LAMBDA(const int /*i*/) { new (devicePtr) ClassWithDeviceFieldData(deviceFieldData); }
                       );
  Kokkos::fence();
  return devicePtr;
}


void delete_class_on_device(ClassWithDeviceFieldData* devicePtr)
{
  Kokkos::parallel_for("device_destruct", stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const int /*i*/) { devicePtr->~ClassWithDeviceFieldData(); }
                       );
  Kokkos::fence();
  Kokkos::kokkos_free<stk::ngp::MemSpace>(static_cast<void*>(devicePtr));
}

class NgpMultiStateFieldTest : public stk::mesh::fixtures::TestHexFixture
{
public:

  template <typename T>
  stk::mesh::Field<T> & create_multistate_field(stk::topology::rank_t rank, const std::string & name, unsigned numStates)
  {
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
    return field;
  }

  void setup_multistate_field(int numStates)
  {
    stk::mesh::FieldBase& field = create_multistate_field<double>(stk::topology::NODE_RANK, "myField", numStates);
    for (int state = 0; state < numStates; ++state) {
      m_fieldStates.push_back(field.field_state(static_cast<stk::mesh::FieldState>(state)));
    }

    get_meta().declare_part("test_part");
  }

  stk::mesh::Field<double>& get_field(stk::mesh::FieldState state) {
    return *static_cast<stk::mesh::Field<double>*>(m_fieldStates[state]);
  }

  template <typename T>
  struct CheckValueUsingDeviceFieldData {
    CheckValueUsingDeviceFieldData(const stk::mesh::ConstFieldData<T, stk::ngp::DeviceSpace>& deviceFieldData,
                                   T expectedValue)
      : m_deviceFieldData(deviceFieldData),
        m_expectedValue(expectedValue)
    {
    }

    KOKKOS_FUNCTION
    void operator()(const stk::mesh::FastMeshIndex& entityIdx) const
    {
      auto entityValues = m_deviceFieldData.entity_values(entityIdx);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        NGP_EXPECT_EQ(entityValues(component), m_expectedValue);
      }
    }

  private:
    stk::mesh::ConstFieldData<T, stk::ngp::DeviceSpace> m_deviceFieldData;
    T m_expectedValue;
  };

  template<typename T>
  void check_field_data_value_on_device(const stk::mesh::NgpMesh& ngpMesh,
                                        const stk::mesh::ConstFieldData<T, stk::ngp::DeviceSpace>& deviceFieldData,
                                        T expectedValue)
  {
    stk::mesh::Selector owned = ngpMesh.get_bulk_on_host().mesh_meta_data().locally_owned_part();
    CheckValueUsingDeviceFieldData<T> checkValueUsingDeviceFieldData(deviceFieldData, expectedValue);
    stk::mesh::for_each_entity_run(ngpMesh, deviceFieldData.entity_rank(), owned, checkValueUsingDeviceFieldData);
  }

  template<typename T>
  void check_field_data_value_on_host(const stk::mesh::BulkData& bulk,
                                      const stk::mesh::Field<T>& hostField,
                                      T expectedValue)
  {
    const stk::mesh::BucketVector& buckets = bulk.get_buckets(hostField.entity_rank(), hostField);

    auto hostFieldData = hostField.data();
    for (const stk::mesh::Bucket* bucket : buckets) {
      for (const stk::mesh::Entity& entity : *bucket) {
        auto value = hostFieldData.entity_values(entity);
        EXPECT_EQ(value(0_comp), expectedValue);
      }
    }
  }

  template <typename ValueType>
  struct CheckValueUsingClass {
    CheckValueUsingClass(const ClassWithDeviceFieldData* deviceClassPointer, ValueType expectedValue)
      : m_deviceClassPointer(deviceClassPointer),
        m_expectedValue(expectedValue)
    {
    }

    KOKKOS_FUNCTION
    void operator()(const stk::mesh::FastMeshIndex& entity) const
    {
      int numComponents = m_deviceClassPointer->num_components_per_entity(entity);
      for (stk::mesh::ComponentIdx component(0); component < numComponents; ++component) {
        NGP_EXPECT_EQ(m_expectedValue, m_deviceClassPointer->access_field_data(entity, component));
      }
    }

  private:
    const ClassWithDeviceFieldData* m_deviceClassPointer;
    ValueType m_expectedValue;
  };

  template<typename T>
  void check_field_data_value_on_device(const stk::mesh::NgpMesh& ngpMesh,
                                        stk::mesh::EntityRank rank,
                                        const ClassWithDeviceFieldData* deviceClassPointer,
                                        T expectedValue)
  {
    stk::mesh::Selector owned = ngpMesh.get_bulk_on_host().mesh_meta_data().locally_owned_part();
    CheckValueUsingClass<T> checkValueUsingClass(deviceClassPointer, expectedValue);
    stk::mesh::for_each_entity_run(ngpMesh, rank, owned, checkValueUsingClass);
  }

  void put_all_nodes_in_part(stk::mesh::BulkData& bulk, stk::mesh::Part& part)
  {
    stk::mesh::EntityVector allNodes;
    bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), allNodes);

    bulk.batch_change_entity_parts(allNodes, {&part}, {});
  }

  stk::mesh::FieldVector m_fieldStates;
};




NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::Field<double>& fieldNew = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateNew));
  stk::mesh::Field<double>& fieldOld = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateOld));

  const double valueNew = 11.1;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);
  }

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);  // Confirm host rotated
  check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);  // Confirm device rotated
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);
  }

  // Trigger sync to device, to make sure data pointers to different states are in correct location
  fieldNew.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  fieldOld.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  fieldNew.synchronize<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  fieldOld.synchronize<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);  // Confirm device rotated
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);
  }

}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 3;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNp1 = get_field(stk::mesh::StateNP1);
  stk::mesh::FieldBase& fieldN   = get_field(stk::mesh::StateN);
  stk::mesh::FieldBase& fieldNm1 = get_field(stk::mesh::StateNM1);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  stk::mesh::field_fill(valueNp1, fieldNp1);
  stk::mesh::field_fill(valueN,   fieldN);
  stk::mesh::field_fill(valueNm1, fieldNm1);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNp1 = fieldNp1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataN   = fieldN.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm1 = fieldNm1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueN);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueNm1);

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNm1);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueN);
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_fourStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 4;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNp1 = get_field(stk::mesh::StateNP1);
  stk::mesh::FieldBase& fieldN   = get_field(stk::mesh::StateN);
  stk::mesh::FieldBase& fieldNm1 = get_field(stk::mesh::StateNM1);
  stk::mesh::FieldBase& fieldNm2 = get_field(stk::mesh::StateNM2);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  const double valueNm2 = 44.4;
  stk::mesh::field_fill(valueNp1, fieldNp1);
  stk::mesh::field_fill(valueN,   fieldN);
  stk::mesh::field_fill(valueNm1, fieldNm1);
  stk::mesh::field_fill(valueNm2, fieldNm2);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNp1 = fieldNp1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataN   = fieldN.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm1 = fieldNm1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm2 = fieldNm2.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueN);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueNm1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm2, valueNm2);

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNm2);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueN);
  check_field_data_value_on_device(ngpMesh, fieldDataNm2, valueNm1);
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_notAllDeviceStatesExist_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 11.1;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNew, valueNew);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#else
  EXPECT_FALSE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNew.has_device_data());
  EXPECT_TRUE(fieldOld.has_device_data());
#else
  EXPECT_FALSE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#endif

  auto fieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNew, valueOld);
  check_field_data_value_on_device(ngpMesh, fieldDataOld, valueNew);
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_notAllDeviceStatesExist_threeStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 3;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNp1 = get_field(stk::mesh::StateNP1);
  stk::mesh::FieldBase& fieldN   = get_field(stk::mesh::StateN);
  stk::mesh::FieldBase& fieldNm1 = get_field(stk::mesh::StateNM1);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  stk::mesh::field_fill(valueNp1, fieldNp1);
  stk::mesh::field_fill(valueN,   fieldN);
  stk::mesh::field_fill(valueNm1, fieldNm1);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNp1 = fieldNp1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNp1);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
#else
  EXPECT_FALSE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNp1.has_device_data());
  EXPECT_TRUE(fieldN.has_device_data());
  EXPECT_TRUE(fieldNm1.has_device_data());
#else
  EXPECT_FALSE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
#endif

  auto fieldDataN   = fieldN.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm1 = fieldNm1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNm1);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueN);
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_notAllDeviceStatesExist_fourStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 4;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNp1 = get_field(stk::mesh::StateNP1);
  stk::mesh::FieldBase& fieldN   = get_field(stk::mesh::StateN);
  stk::mesh::FieldBase& fieldNm1 = get_field(stk::mesh::StateNM1);
  stk::mesh::FieldBase& fieldNm2 = get_field(stk::mesh::StateNM2);

  const double valueNp1 = 11.1;
  const double valueN   = 22.2;
  const double valueNm1 = 33.3;
  const double valueNm2 = 44.4;
  stk::mesh::field_fill(valueNp1, fieldNp1);
  stk::mesh::field_fill(valueN,   fieldN);
  stk::mesh::field_fill(valueNm1, fieldNm1);
  stk::mesh::field_fill(valueNm2, fieldNm2);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNp1 = fieldNp1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm2 = fieldNm2.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm2, valueNm2);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
  EXPECT_TRUE(fieldNm2.has_device_data());
#else
  EXPECT_FALSE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
  EXPECT_FALSE(fieldNm2.has_device_data());
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNp1.has_device_data());
  EXPECT_TRUE(fieldN.has_device_data());
  EXPECT_TRUE(fieldNm1.has_device_data());
  EXPECT_TRUE(fieldNm2.has_device_data());
#else
  EXPECT_FALSE(fieldNp1.has_device_data());
  EXPECT_FALSE(fieldN.has_device_data());
  EXPECT_FALSE(fieldNm1.has_device_data());
  EXPECT_FALSE(fieldNm2.has_device_data());
#endif

  auto fieldDataN   = fieldN.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto fieldDataNm1 = fieldNm1.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNp1, valueNm2);
  check_field_data_value_on_device(ngpMesh, fieldDataN,   valueNp1);
  check_field_data_value_on_device(ngpMesh, fieldDataNm1, valueN);
  check_field_data_value_on_device(ngpMesh, fieldDataNm2, valueNm1);
}

NGP_TEST_F(NgpMultiStateFieldTest, syncBasedRotation_notAllDeviceStatesExist_twoStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 11.1;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto fieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, fieldDataNew, valueNew);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#else
  EXPECT_FALSE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#endif

  stk::mesh::sync_to_host_and_mark_modified(get_meta());
  get_bulk().update_field_data_states();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_TRUE(fieldNew.has_device_data());  // No change.  Nothing auto-registered when not rotating on device.
  EXPECT_FALSE(fieldOld.has_device_data());
#else
  EXPECT_FALSE(fieldNew.has_device_data());
  EXPECT_FALSE(fieldOld.has_device_data());
#endif

  {
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);  // All rotated data synced properly to device
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);
  }
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateDeviceStates_syncStatesUnchanged)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  EXPECT_FALSE(fieldNew.need_sync_to_host());
  EXPECT_FALSE(fieldOld.need_sync_to_host());
#if defined(STK_USE_DEVICE_MESH)
  EXPECT_FALSE(fieldNew.need_sync_to_device());
  EXPECT_FALSE(fieldOld.need_sync_to_device());
#else
  EXPECT_TRUE(fieldNew.need_sync_to_device());  // No easy way to clear initial flag state in host build
  EXPECT_TRUE(fieldOld.need_sync_to_device());
#endif

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);
}

NGP_TEST_F(NgpMultiStateFieldTest, copyHasCorrectDataAfterSyncBasedStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto copyOfFieldDataNew(deviceFieldDataNew);

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);

  stk::mesh::sync_to_host_and_mark_modified(get_meta());
  get_bulk().update_field_data_states();

  deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Reacquire (and sync) device FieldData
  deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);

  check_field_data_value_on_device(ngpMesh, copyOfFieldDataNew, valueOld);
}

NGP_TEST_F(NgpMultiStateFieldTest, copyHasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto copyOfFieldDataNew(deviceFieldDataNew);

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);
  check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);

  check_field_data_value_on_device(ngpMesh, copyOfFieldDataNew, valueOld);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentDeviceFieldData_hasCorrectDataAfterSyncBasedStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  ClassWithDeviceFieldData* persistentDeviceClass = create_class_on_device(deviceFieldDataNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  stk::mesh::sync_to_host_and_mark_modified(get_meta());
  get_bulk().update_field_data_states();

  deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Reacquire (and sync) device FieldData
  deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentDeviceFieldData_hasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, fieldNew);
  stk::mesh::field_fill(valueOld, fieldOld);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  ClassWithDeviceFieldData* persistentDeviceClass = create_class_on_device(deviceFieldDataNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToDeviceCountAfterSyncBasedStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());
  }

  fieldNew.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device()+1);
#else
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());  // No device data, so no automatic syncs
#endif

  stk::mesh::sync_to_host_and_mark_modified(get_meta());
  get_bulk().update_field_data_states();  // Must also re-sync to device for this rotation mode, but don't bother here

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_device()+1, fieldOld.num_syncs_to_device());  // Sync count follows states
#else
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());  // No device data, so no automatic syncs
#endif
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToDeviceCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());
  }

  fieldNew.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device()+1);
#else
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());  // No device data, so no automatic syncs
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_device()+1, fieldOld.num_syncs_to_device());  // Sync count follows states
#else
  EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());  // No device data, so no automatic syncs
#endif
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToHostCountAfterSyncBasedStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

    EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());
  }

  fieldOld.synchronize<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_host()+1, fieldOld.num_syncs_to_host());
#else
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host());  // No device data, so no automatic syncs
#endif

  get_bulk().update_field_data_states();  // Must also re-sync to device for this rotation mode, but don't bother here

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host()+1);  // Sync count follows states
#else
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host());  // No device data, so no automatic syncs
#endif
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToHostCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

    EXPECT_EQ(fieldNew.num_syncs_to_device(), fieldOld.num_syncs_to_device());
  }

  fieldOld.synchronize<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_host()+1, fieldOld.num_syncs_to_host());
#else
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host());  // No device data, so no automatic syncs
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host()+1);  // Sync count follows states
#else
  EXPECT_EQ(fieldNew.num_syncs_to_host(), fieldOld.num_syncs_to_host());  // No device data, so no automatic syncs
#endif
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentModifyOnHostAfterSyncBasedStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    EXPECT_FALSE(fieldNew.need_sync_to_device());
    EXPECT_FALSE(fieldOld.need_sync_to_device());
  }

  fieldNew.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  EXPECT_TRUE(fieldNew.need_sync_to_device());
  EXPECT_FALSE(fieldOld.need_sync_to_device());

  get_bulk().update_field_data_states();  // Must also re-sync to device for this rotation mode, but don't bother here

  EXPECT_FALSE(fieldNew.need_sync_to_device());
  EXPECT_TRUE(fieldOld.need_sync_to_device());
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentModifyOnHostAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::FieldBase& fieldNew = get_field(stk::mesh::StateNew);
  stk::mesh::FieldBase& fieldOld = get_field(stk::mesh::StateOld);

  {
    auto deviceFieldDataNew = fieldNew.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataOld = fieldOld.data<double, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    EXPECT_FALSE(fieldNew.need_sync_to_device());
    EXPECT_FALSE(fieldOld.need_sync_to_device());
  }

  fieldNew.synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  EXPECT_TRUE(fieldNew.need_sync_to_device());
  EXPECT_FALSE(fieldOld.need_sync_to_device());

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  EXPECT_FALSE(fieldNew.need_sync_to_device());
  EXPECT_TRUE(fieldOld.need_sync_to_device());
}

NGP_TEST_F(NgpMultiStateFieldTest, rotateOnDevice_withMeshModification)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::Field<double>& fieldOld = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateOld));
  stk::mesh::Field<double>& fieldNew = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateNew));

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, fieldOld);
  stk::mesh::field_fill(valueNew, fieldNew);

  stk::mesh::Part& testPart = *get_meta().get_part("test_part");

  stk::mesh::NgpMesh* ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
  check_field_data_value_on_device(*ngpMesh, deviceFieldDataOld, valueOld);
  check_field_data_value_on_device(*ngpMesh, deviceFieldDataNew, valueNew);

  put_all_nodes_in_part(get_bulk(), testPart);
  deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  ngpMesh = &stk::mesh::get_updated_ngp_mesh(get_bulk());
  check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);
  check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);
  check_field_data_value_on_device(*ngpMesh, deviceFieldDataOld, valueNew);
  check_field_data_value_on_device(*ngpMesh, deviceFieldDataNew, valueOld);
}

NGP_TEST_F(NgpMultiStateFieldTest, updateUnmodifiedBucketWithSyncBasedRotation_afterMeshModification)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::Field<double>& fieldOld = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateOld));
  stk::mesh::Field<double>& fieldNew = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateNew));

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, fieldOld);
  stk::mesh::field_fill(valueNew, fieldNew);

  stk::mesh::Part& testPart = *get_meta().get_part("test_part");

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
    check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
  }

  // New node goes in new Bucket, leaving others alone
  get_bulk().modification_begin();
  {
    const stk::mesh::Entity newNode = get_bulk().declare_node(1000, stk::mesh::PartVector{&testPart});
    auto hostFieldDataOld = fieldOld.data<stk::mesh::ReadWrite>();
    auto hostFieldDataNew = fieldNew.data<stk::mesh::ReadWrite>();
    auto fieldOldValue = hostFieldDataOld.entity_values(newNode);
    auto fieldNewValue = hostFieldDataNew.entity_values(newNode);
    fieldOldValue(0_comp) = valueOld;
    fieldNewValue(0_comp) = valueNew;
  }
  get_bulk().modification_end();

  stk::mesh::sync_to_host_and_mark_modified(get_meta());
  get_bulk().update_field_data_states();

  {
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    check_field_data_value_on_host(get_bulk(), fieldOld, valueNew);  // Rotated properly on host
    check_field_data_value_on_host(get_bulk(), fieldNew, valueOld);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueNew);  // All rotated data synced properly to device
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueOld);
  }
}

template <typename Space = stk::ngp::HostSpace>
class TestConstFieldDataBytes : public stk::mesh::ConstFieldDataBytes<Space>
{
public:
  using stk::mesh::ConstFieldDataBytes<Space>::ConstFieldDataBytes;

  int get_field_meta_data_mod_count() const {
    return this->m_localFieldMetaDataModCount;
  }
};

NGP_TEST_F(NgpMultiStateFieldTest, deviceRotationUpdatesAllStates)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int numStates = 2;
  setup_multistate_field(numStates);
  setup_mesh(2, 2, 2);

  stk::mesh::Field<double>& fieldOld = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateOld));
  stk::mesh::Field<double>& fieldNew = static_cast<stk::mesh::Field<double>&>(get_field(stk::mesh::StateNew));

  const double valueOld = 2.0;
  const double valueNew = 4.0;
  stk::mesh::field_fill(valueOld, fieldOld);
  stk::mesh::field_fill(valueNew, fieldNew);

  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    auto deviceFieldDataOld = fieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto deviceFieldDataNew = fieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    check_field_data_value_on_host(get_bulk(), fieldOld, valueOld);
    check_field_data_value_on_host(get_bulk(), fieldNew, valueNew);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataOld, valueOld);
    check_field_data_value_on_device(ngpMesh, deviceFieldDataNew, valueNew);
  }

  get_bulk().modification_begin();
  const stk::mesh::Entity newNode = get_bulk().declare_node(1000);
  fieldOld.data<stk::mesh::ReadWrite>().entity_values(newNode)() = valueOld;
  fieldNew.data<stk::mesh::ReadWrite>().entity_values(newNode)() = valueNew;
  get_bulk().modification_end();

  fieldNew.synchronize<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Trigger an update of just this state

  auto dataBytesOld = fieldOld.data_bytes<const std::byte, stk::ngp::DeviceSpace>();
  auto dataBytesNew = fieldNew.data_bytes<const std::byte, stk::ngp::DeviceSpace>();
  auto* testDataBytesOld = static_cast<TestConstFieldDataBytes<stk::ngp::DeviceSpace>*>(&dataBytesOld);
  auto* testDataBytesNew = static_cast<TestConstFieldDataBytes<stk::ngp::DeviceSpace>*>(&dataBytesNew);

  const int fieldMetaDataModCountOld = testDataBytesOld->get_field_meta_data_mod_count();
  const int fieldMetaDataModCountNew = testDataBytesNew->get_field_meta_data_mod_count();

#if defined(STK_USE_DEVICE_MESH)
  EXPECT_GT(fieldMetaDataModCountNew, fieldMetaDataModCountOld);  // State NEW updated, state OLD not
#else
  EXPECT_EQ(fieldMetaDataModCountNew, fieldMetaDataModCountOld);  // Both updated when no real device data
#endif

  const bool alsoRotateOnDevice = true;
  get_bulk().update_field_data_states(alsoRotateOnDevice);

  dataBytesOld = fieldOld.data_bytes<const std::byte, stk::ngp::DeviceSpace>();
  dataBytesNew = fieldNew.data_bytes<const std::byte, stk::ngp::DeviceSpace>();
  testDataBytesOld = static_cast<TestConstFieldDataBytes<stk::ngp::DeviceSpace>*>(&dataBytesOld);
  testDataBytesNew = static_cast<TestConstFieldDataBytes<stk::ngp::DeviceSpace>*>(&dataBytesNew);

  const int rotatedFieldMetaDataModCountOld = testDataBytesOld->get_field_meta_data_mod_count();
  const int rotatedFieldMetaDataModCountNew = testDataBytesNew->get_field_meta_data_mod_count();

  EXPECT_EQ(rotatedFieldMetaDataModCountNew, rotatedFieldMetaDataModCountOld);  // Both states updated to same level
}

}
