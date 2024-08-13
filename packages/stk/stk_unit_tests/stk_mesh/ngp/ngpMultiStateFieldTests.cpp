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
                       MY_LAMBDA(const int i) { new (devicePtr) ClassWithNgpField(ngpField); }
                       );
  Kokkos::fence();
  return devicePtr;
}

void delete_class_on_device(ClassWithNgpField* devicePtr)
{
  Kokkos::parallel_for("device_destruct", stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const int i) { devicePtr->~ClassWithNgpField(); }
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

  void setup_multistate_field()
  {
    const unsigned numStates = 2;
    m_fieldNew = &create_multistate_field<double>(stk::topology::NODE_RANK, "myField", numStates);
    EXPECT_EQ(stk::mesh::StateNew, m_fieldNew->state());

    m_fieldOld = &m_fieldNew->field_of_state(stk::mesh::StateOld);
    EXPECT_EQ(stk::mesh::StateOld, m_fieldOld->state());
  }

  stk::mesh::Field<double>& get_field_new() { return *m_fieldNew; }
  stk::mesh::Field<double>& get_field_old() { return *m_fieldOld; }

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

  stk::mesh::Field<double>* m_fieldNew;
  stk::mesh::Field<double>* m_fieldOld;
};

NGP_TEST_F(NgpMultiStateFieldTest, multistateField_rotateDeviceStates_syncStatesUnchanged)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field_new(), stk::ngp::HostExecSpace());
  stk::mesh::field_fill(valueOld, get_field_old(), stk::ngp::HostExecSpace());

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());

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

NGP_TEST_F(NgpMultiStateFieldTest, multistateField_copyHasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field_new(), stk::ngp::HostExecSpace());
  stk::mesh::field_fill(valueOld, get_field_old(), stk::ngp::HostExecSpace());

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());
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

NGP_TEST_F(NgpMultiStateFieldTest, multistateField_copyHasWrongDataAfterDeviceStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field_new(), stk::ngp::HostExecSpace());
  stk::mesh::field_fill(valueOld, get_field_old(), stk::ngp::HostExecSpace());

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());
  stk::mesh::NgpField<double> copyOfNgpFieldNew(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueNew);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueOld);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
  check_field_data_value_on_device(ngpMesh, ngpFieldOld, valueNew);

#ifdef STK_USE_DEVICE_MESH
  const double valueNewBecauseCopyOfDeviceFieldGotDisconnectedByDeviceRotation = valueNew;
  check_field_data_value_on_device(ngpMesh, copyOfNgpFieldNew, valueNewBecauseCopyOfDeviceFieldGotDisconnectedByDeviceRotation);
#else
  check_field_data_value_on_device(ngpMesh, copyOfNgpFieldNew, valueOld);
#endif
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentDeviceField_hasCorrectDataAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field_new(), stk::ngp::HostExecSpace());
  stk::mesh::field_fill(valueOld, get_field_old(), stk::ngp::HostExecSpace());

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());

  ClassWithNgpField* persistentDeviceClass = create_class_on_device(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  perform_field_state_rotation();

  ngpFieldNew.sync_to_device();
  ngpFieldOld.sync_to_device();

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentDeviceField_hasWrongDataAfterDeviceStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  const double valueNew = 44.4;
  const double valueOld = 22.2;
  stk::mesh::field_fill(valueNew, get_field_new(), stk::ngp::HostExecSpace());
  stk::mesh::field_fill(valueOld, get_field_old(), stk::ngp::HostExecSpace());

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());

  ClassWithNgpField* persistentDeviceClass = create_class_on_device(ngpFieldNew);

  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNew);

  const bool rotateDeviceNgpFieldStates = true;
  get_bulk().update_field_data_states(rotateDeviceNgpFieldStates);

  ngpFieldNew.sync_to_device();
  ngpFieldOld.sync_to_device();

  check_field_data_value_on_device(ngpMesh, ngpFieldNew, valueOld);
#ifdef STK_USE_DEVICE_MESH
  const double valueNewBecausePersistentDeviceFieldGotDisconnectedByDeviceRotation = valueNew;
  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueNewBecausePersistentDeviceFieldGotDisconnectedByDeviceRotation);
#else
  check_field_data_value_on_device(ngpMesh, stk::topology::NODE_RANK, persistentDeviceClass, valueOld);
#endif

  delete_class_on_device(persistentDeviceClass);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToDeviceCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device());

  ngpFieldNew.modify_on_host();
  ngpFieldNew.sync_to_device();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device()+1);

  perform_field_state_rotation();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device()+1, ngpFieldOld.num_syncs_to_device());
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentSyncToHostCountAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());

  EXPECT_EQ(ngpFieldNew.num_syncs_to_device(), ngpFieldOld.num_syncs_to_device());

  ngpFieldOld.modify_on_device();
  ngpFieldOld.sync_to_host();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_host()+1, ngpFieldOld.num_syncs_to_host());

  perform_field_state_rotation();

  EXPECT_EQ(ngpFieldNew.num_syncs_to_host(), ngpFieldOld.num_syncs_to_host()+1);
}

NGP_TEST_F(NgpMultiStateFieldTest, persistentModifyOnHostAfterStateRotation)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  setup_multistate_field();
  setup_mesh(2, 2, 2);

  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(get_field_new());
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(get_field_old());
  ngpFieldNew.clear_sync_state();
  ngpFieldOld.clear_sync_state();

  ngpFieldNew.modify_on_host();

  EXPECT_TRUE(ngpFieldNew.need_sync_to_device());
  EXPECT_FALSE(ngpFieldOld.need_sync_to_device());

  get_bulk().update_field_data_states();

  EXPECT_FALSE(ngpFieldNew.need_sync_to_device());
  EXPECT_TRUE(ngpFieldOld.need_sync_to_device());
}
