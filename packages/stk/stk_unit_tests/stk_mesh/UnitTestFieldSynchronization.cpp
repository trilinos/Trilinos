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
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include "ngp/NgpUnitTestUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/ConstFieldData.hpp>
#include <stk_mesh/base/FieldDataBase.hpp>
#include <stk_mesh/base/FieldIndexTypes.hpp>

namespace {

class FieldDataSynchronization : public stk::unit_test_util::MeshFixture
{
public:
  FieldDataSynchronization()
    : m_field(nullptr)
  {}

  void build_two_element_mesh_with_nodal_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    m_field = &get_meta().declare_field<int>(stk::topology::NODE_RANK, "IntField", 1);
    stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), nullptr);

    stk::mesh::fixtures::HexFixture::fill_mesh(2, 1, 1, get_bulk());

    m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Trigger creation of default-initialized device data object
  }

  void set_field_values_on_host(stk::mesh::FieldData<int, stk::ngp::HostSpace>& fieldData, int scaleFactor)
  {
    stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::NODE_RANK);
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const stk::mesh::EntityId nodeId = get_bulk().identifier(node);
        stk::mesh::EntityValues fieldValues = fieldData.entity_values(node);
        fieldValues() = nodeId * scaleFactor;
      }
    }
  }

  void check_field_values_on_host(stk::mesh::ConstFieldData<int, stk::ngp::HostSpace>& fieldData, int scaleFactor)
  {
    stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::NODE_RANK);
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const stk::mesh::EntityId nodeId = get_bulk().identifier(node);
        stk::mesh::EntityValues fieldValues = fieldData.entity_values(node);
        EXPECT_EQ(fieldValues(), static_cast<int>(nodeId * scaleFactor));
      }
    }
  }

  void set_field_values_on_device(stk::mesh::FieldData<int, stk::ngp::DeviceSpace>& fieldData, int scaleFactor)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, get_meta().universal_part(),
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node) {
        const stk::mesh::EntityId nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, node));
        auto fieldValues = fieldData.entity_values(node);
        fieldValues() = nodeId * scaleFactor;
      }
    );
  }

  void check_field_values_on_device(stk::mesh::ConstFieldData<int, stk::ngp::DeviceSpace>& fieldData, int scaleFactor)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, get_meta().universal_part(),
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node) {
        const stk::mesh::EntityId nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, node));
        auto fieldValues = fieldData.entity_values(node);
        NGP_EXPECT_EQ(fieldValues(), static_cast<int>(nodeId * scaleFactor));
      }
    );
  }

  void check_field_values_on_device(stk::mesh::NgpField<int>& ngpField, int scaleFactor)
  {
    stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, get_meta().universal_part(),
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node) {
        const stk::mesh::EntityId nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, node));
        const int& fieldValues = ngpField(node, 0);
        NGP_EXPECT_EQ(fieldValues, static_cast<int>(nodeId * scaleFactor));
      }
    );
  }

  void check_field_values_on_host(stk::mesh::HostField<int>& hostField, int scaleFactor)
  {
    stk::mesh::HostMesh hostMesh(get_bulk());

    stk::mesh::for_each_entity_run(hostMesh, stk::topology::NODE_RANK, get_meta().universal_part(),
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node) {
        const stk::mesh::EntityId nodeId = hostMesh.identifier(hostMesh.get_entity(stk::topology::NODE_RANK, node));
        const int& fieldValues = hostField(node, 0);
        NGP_EXPECT_EQ(fieldValues, static_cast<int>(nodeId * scaleFactor));
      }, stk::ngp::HostExecSpace()
    );
  }

protected:
  stk::mesh::Field<int>* m_field;
};


//==============================================================================
// Normal usage

TEST_F(FieldDataSynchronization, simpleHostUsage)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  stk::mesh::FieldData<int> hostFieldData = m_field->data<stk::mesh::ReadWrite>();
  set_field_values_on_host(hostFieldData, 1);
  check_field_values_on_host(hostFieldData, 1);
}

TEST_F(FieldDataSynchronization, simpleHostUsage_twoOverlappingFieldData)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  stk::mesh::FieldData<int> hostFieldData = m_field->data<stk::mesh::ReadWrite>();
  set_field_values_on_host(hostFieldData, 1);

  stk::mesh::ConstFieldData<int> constHostFieldData = m_field->data();
  check_field_values_on_host(constHostFieldData, 1);
}

TEST_F(FieldDataSynchronization, mixedApiUsage)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  set_field_values_on_host(hostFieldData, 1);

  auto ngpField = stk::mesh::get_updated_ngp_field<int>(*m_field);
  ngpField.modify_on_host();
  ngpField.sync_to_device();
  check_field_values_on_device(ngpField, 1);
}

TEST_F(FieldDataSynchronization, mixedApiUsage_HostField_syncToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  set_field_values_on_host(hostFieldData, 1);

  auto ngpField = stk::mesh::get_updated_ngp_field<int>(*m_field);
  ngpField.modify_on_host();
  stk::mesh::HostField<int> hostField(*m_field);
  hostField.sync_to_device();
  check_field_values_on_device(ngpField, 1);
}

TEST_F(FieldDataSynchronization, mixedApiUsage_HostField_syncToHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  set_field_values_on_device(deviceFieldData, 1);

  auto ngpField = stk::mesh::get_updated_ngp_field<int>(*m_field);
  ngpField.modify_on_device();
  stk::mesh::HostField<int> hostField(*m_field);
  hostField.sync_to_host();
  check_field_values_on_host(hostField, 1);
}

TEST_F(FieldDataSynchronization, writeOnHost_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, readWriteOnDevice_readOnlyOnHost_syncedToHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, overwriteAllOnHost_overwriteAllOnDevice_overwriteAllOnHost_notSynced)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto deviceFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>();
#if defined(STK_USE_DEVICE_MESH) && !defined(STK_UNIFIED_MEMORY)
    check_field_values_on_device(deviceFieldData, 0);  // Host values not synced to device
#else
    check_field_values_on_device(deviceFieldData, 1);  // Host values already implicitly on device
#endif

    set_field_values_on_device(deviceFieldData, 2);
  }
  {
    auto hostFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>();
#if defined(STK_USE_DEVICE_MESH) && !defined(STK_UNIFIED_MEMORY)
    check_field_values_on_host(hostFieldData, 1);  // Device values not synced to host
#else
    check_field_values_on_host(hostFieldData, 2);  // Host values already implicitly on device
#endif
    set_field_values_on_host(hostFieldData, 3);
  }
}


TEST_F(FieldDataSynchronization, readWriteOnHostWithoutDestroying_readWriteOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnHostWithoutDestroying_readOnlyOnDevice_noError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnHostWithoutDestroying_overwriteAllOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readWriteOnHostWithoutDestroying_readOnlyOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnHostWithoutDestroying_readWriteOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readWriteOnHostWithoutDestroying_overwriteAllOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnHostWithoutDestroying_readWriteOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnHostWithoutDestroying_overwriteAllOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnHostWithoutDestroying_readOnlyOnDevice_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto hostFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}


TEST_F(FieldDataSynchronization, readWriteOnDeviceWithoutDestroying_readWriteOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnDeviceWithoutDestroying_readOnlyOnHost_noError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnDeviceWithoutDestroying_overwriteAllOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readWriteOnDeviceWithoutDestroying_readOnlyOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnDeviceWithoutDestroying_readWriteOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readWriteOnDeviceWithoutDestroying_overwriteAllOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnDeviceWithoutDestroying_readWriteOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, readOnlyOnDeviceWithoutDestroying_overwriteAllOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::OverwriteAll, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, overwriteAllOnDeviceWithoutDestroying_readOnlyOnHost_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  [[maybe_unused]] auto deviceFieldData = m_field->data<stk::mesh::OverwriteAll, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}


struct FakeHostKernel {
  FakeHostKernel() = default;
  ~FakeHostKernel() = default;
  stk::mesh::FieldData<int, stk::ngp::HostSpace> m_fieldData;
};

using HostKernelContainer = std::vector<FakeHostKernel>;

TEST_F(FieldDataSynchronization, persistentHostCopy_triggerSynchronizeManually)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  HostKernelContainer hostKernelContainer(1);
  hostKernelContainer[0].m_fieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();

  {
    set_field_values_on_host(hostKernelContainer[0].m_fieldData, 1);
  }

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
#if defined(STK_USE_DEVICE_MESH) && !defined(STK_UNIFIED_MEMORY)
    check_field_values_on_device(deviceFieldData, 0);  // No sync to device; initial values
#else
    check_field_values_on_device(deviceFieldData, 1);  // Value implicitly synced to device
#endif
  }

  {
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();  // Flag so that synced to device next time
    set_field_values_on_host(hostKernelContainer[0].m_fieldData, 2);
  }

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(deviceFieldData, 2);
  }
}


struct FakeDeviceKernel {
  KOKKOS_DEFAULTED_FUNCTION FakeDeviceKernel() = default;
  KOKKOS_DEFAULTED_FUNCTION ~FakeDeviceKernel() = default;
  stk::mesh::FieldData<int, stk::ngp::DeviceSpace> m_fieldData;
};

using DeviceKernelContainer = Kokkos::View<FakeDeviceKernel*, stk::ngp::UVMMemSpace>;

TEST_F(FieldDataSynchronization, persistentDeviceCopy_triggerSynchronizeManually)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  DeviceKernelContainer deviceKernelContainer(Kokkos::view_alloc(Kokkos::WithoutInitializing, "DeviceKernelView"), 1);

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    new (&deviceKernelContainer[0].m_fieldData) stk::mesh::FieldData<int, stk::ngp::DeviceSpace>(deviceFieldData);
  }

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

  {
#if defined(STK_USE_DEVICE_MESH) && !defined(STK_UNIFIED_MEMORY)
    check_field_values_on_device(deviceKernelContainer[0].m_fieldData, 0);  // No sync to device; initial values
#else
    check_field_values_on_device(deviceKernelContainer[0].m_fieldData, 1);  // Value implicitly synced to device
#endif
  }

  {
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    check_field_values_on_device(deviceKernelContainer[0].m_fieldData, 1);
  }

  deviceKernelContainer[0].m_fieldData.~FieldData();
}


//==============================================================================
// Copy construction on host

TEST_F(FieldDataSynchronization, writeOnHostFromCopyConstruction_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    auto copyHostFieldData = hostFieldData;
    set_field_values_on_host(copyHostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromCopyConstruction_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    auto copyHostFieldData = hostFieldData;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(copyHostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostFromCopyConstruction_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  {
    auto copyHostFieldData = hostFieldData;
    set_field_values_on_host(copyHostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
  EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromCopyConstruction_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  {
    auto copyHostFieldData = hostFieldData;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(copyHostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}


//==============================================================================
// Move construction on host

TEST_F(FieldDataSynchronization, writeOnHostFromMoveConstruction_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    // C++17 copy elision is used for the return value of the data() method, so move construction must
    // be triggered manually.
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    auto copyHostFieldData = std::move(hostFieldData);
    set_field_values_on_host(copyHostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromMoveConstruction_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    auto copyHostFieldData = std::move(hostFieldData);
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(copyHostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostFromMoveConstruction_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  auto copyHostFieldData = std::move(hostFieldData);
  {
    set_field_values_on_host(copyHostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
  EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromMoveConstruction_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  auto copyHostFieldData = std::move(hostFieldData);
  {
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(copyHostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}


//==============================================================================
// Copy assignment on host

TEST_F(FieldDataSynchronization, writeOnHostFromCopyAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    auto hostFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromCopyAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    auto hostFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostAddingUnsynchronizedFromCopyAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    auto hostFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostRemovingUnsynchronizedFromCopyAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    auto hostFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostFromCopyAssignment_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  {
    auto hostFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromCopyAssignment_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  {
    auto hostFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnHostAddingUnsynchronizedFromCopyAssignment_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  {
    auto hostFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnHostRemovingUnsynchronizedFromCopyAssignment_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  {
    auto hostFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = hostFieldData2;
    set_field_values_on_host(hostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}


//==============================================================================
// Move assignment on host

TEST_F(FieldDataSynchronization, writeOnHostFromMoveAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromMoveAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostAddingUnsynchronizedFromMoveAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostRemovingUnsynchronizedFromMoveAssignment_readOnDevice_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }
  {
    auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(constDeviceFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnHostFromMoveAssignment_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  {
    hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnHostFromMoveAssignment_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  {
    hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnHostAddingUnsynchronizedFromMoveAssignment_readOnDeviceWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
  {
    hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
  auto constDeviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_values_on_device(constDeviceFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnHostRemovingUnsynchronizedFromMoveAssignment_readOnDeviceWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto hostFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::HostSpace>();
  {
    hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>()));
#endif
}


//==============================================================================
// Copy construction on device

TEST_F(FieldDataSynchronization, writeOnDeviceFromCopyConstruction_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto copyDeviceFieldData = deviceFieldData;
    set_field_values_on_device(copyDeviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromCopyConstruction_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    auto copyDeviceFieldData = deviceFieldData;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(copyDeviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceFromCopyConstruction_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  {
    auto copyDeviceFieldData = deviceFieldData;
    set_field_values_on_device(copyDeviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
  EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromCopyConstruction_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  {
    auto copyDeviceFieldData = deviceFieldData;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(copyDeviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}


//==============================================================================
// Move construction on device

TEST_F(FieldDataSynchronization, writeOnDeviceFromMoveConstruction_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    // C++17 copy elision is used for the return value of the data() method, so move construction must
    // be triggered manually.
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto copyDeviceFieldData = std::move(deviceFieldData);
    set_field_values_on_device(copyDeviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromMoveConstruction_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    auto copyDeviceFieldData = std::move(deviceFieldData);
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(copyDeviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceFromMoveConstruction_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  auto copyDeviceFieldData = std::move(deviceFieldData);
  {
    set_field_values_on_device(copyDeviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
  EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromMoveConstruction_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  auto copyDeviceFieldData = std::move(deviceFieldData);
  {
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(copyDeviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}


//==============================================================================
// Copy assignment on device

TEST_F(FieldDataSynchronization, writeOnDeviceFromCopyAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto deviceFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromCopyAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    auto deviceFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceAddingUnsynchronizedFromCopyAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    auto deviceFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceRemovingUnsynchronizedFromCopyAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    auto deviceFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceFromCopyAssignment_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  {
    auto deviceFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromCopyAssignment_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  {
    auto deviceFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnDeviceAddingUnsynchronizedFromCopyAssignment_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  {
    auto deviceFieldData2 = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnDeviceRemovingUnsynchronizedFromCopyAssignment_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  {
    auto deviceFieldData2 = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = deviceFieldData2;
    set_field_values_on_device(deviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}


//==============================================================================
// Move assignment on device

TEST_F(FieldDataSynchronization, writeOnDeviceFromMoveAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromMoveAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceAddingUnsynchronizedFromMoveAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceRemovingUnsynchronizedFromMoveAssignment_readOnHost_syncedToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }
  {
    auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
    check_field_values_on_host(constHostFieldData, 1);
  }
}

TEST_F(FieldDataSynchronization, writeOnDeviceFromMoveAssignment_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  {
    deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, writeUnsynchronizedOnDeviceFromMoveAssignment_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  {
    deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    m_field->synchronize<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnDeviceAddingUnsynchronizedFromMoveAssignment_readOnHostWithoutDestroying_noThrow)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  {
    deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }

  EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
  auto constHostFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  check_field_values_on_host(constHostFieldData, 1);
}

TEST_F(FieldDataSynchronization, writeOnDeviceRemovingUnsynchronizedFromMoveAssignment_readOnHostWithoutDestroying_throwsSyncError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  auto deviceFieldData = m_field->data<stk::mesh::Unsynchronized, stk::ngp::DeviceSpace>();
  {
    deviceFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    set_field_values_on_device(deviceFieldData, 1);
  }

#ifdef STK_USE_DEVICE_MESH
    EXPECT_ANY_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#else
    EXPECT_NO_THROW((m_field->data<stk::mesh::ReadOnly, stk::ngp::HostSpace>()));
#endif
}

TEST_F(FieldDataSynchronization, interleavedOldAndNewDeviceAccess_properlySyncsToDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_two_element_mesh_with_nodal_field();

  {
    // Initial host values implicitly synced to device during construction
    auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  }

  {
    // Set values on host and automatically mark as modified
    auto hostFieldData = m_field->data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();
    set_field_values_on_host(hostFieldData, 1);
  }

  {
    // This needs to grab a handle to device data without messing with the modify/sync state
    auto ngpField = stk::mesh::get_updated_ngp_field<int>(*m_field);
  }

  {
    // Access data on device, which should automatically sync the modified host values
    auto deviceFieldData = m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_values_on_device(deviceFieldData, 1);
  }
}

}
