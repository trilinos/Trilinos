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

#ifndef STK_MESH_DEVICE_FIELD_HPP
#define STK_MESH_DEVICE_FIELD_HPP

#include "stk_util/stk_config.h"
#include "Kokkos_Core.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/FieldData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"
#include "stk_mesh/base/EntityFieldData.hpp"

namespace stk::mesh {

constexpr unsigned INVALID_ORDINAL = 9999999;

namespace impl {
KOKKOS_FORCEINLINE_FUNCTION
constexpr bool is_on_host() {
  KOKKOS_IF_ON_HOST(return true;);
  KOKKOS_IF_ON_DEVICE(return false;);
}
}

template <typename T, typename NgpMemSpace>
class DeviceField : public NgpFieldBase
{
public:
  using ExecSpace = stk::ngp::ExecSpace;
  using MemSpace = NgpMemSpace;
  using value_type = T;

  KOKKOS_FUNCTION
  DeviceField()
    : NgpFieldBase(),
      m_hostField(nullptr)
  {}

  DeviceField(const FieldBase& stkField, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      m_hostField(&stkField),
      m_deviceFieldData(stkField.data<T, Unsynchronized, stk::ngp::DeviceSpace>())
  {
    STK_ThrowRequireMsg(isFromGetUpdatedNgpField, "NgpField must be obtained from get_updated_ngp_field()");
  }

  KOKKOS_DEFAULTED_FUNCTION DeviceField(const DeviceField&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField(DeviceField&&) = default;
  KOKKOS_FUNCTION ~DeviceField() {}
  KOKKOS_DEFAULTED_FUNCTION DeviceField& operator=(const DeviceField&) = default;
  KOKKOS_DEFAULTED_FUNCTION DeviceField& operator=(DeviceField&&) = default;

  void update_field(const ExecSpace& newExecSpace) override
  {
    set_execution_space(newExecSpace);
    update_field();
  }

  void update_field(ExecSpace&& newExecSpace) override
  {
    set_execution_space(std::forward<ExecSpace>(newExecSpace));
    update_field();
  }

  size_t num_syncs_to_host() const override { return m_hostField->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return m_hostField->num_syncs_to_device(); }

  void fence() override {
    get_execution_space().fence();
    reset_execution_space();
  }

  void modify_on_host() override
  {
    set_modify_on_host();
  }

  void modify_on_host(const Selector& /*selector*/) override
  {
    modify_on_host();
  }

  void modify_on_device() override
  {
    set_modify_on_device();
  }

  void modify_on_device(const Selector& /*selector*/) override
  {
    modify_on_device();
  }

  void clear_sync_state() override
  {
    m_hostField->clear_sync_state();
  }

  void clear_host_sync_state() override
  {
    m_hostField->clear_host_sync_state();
  }

  void clear_device_sync_state() override
  {
    m_hostField->clear_device_sync_state();
  }

  void sync_to_host() override
  {
    sync_to_host(stk::ngp::ExecSpace());
  }

  void sync_to_host(const ExecSpace& execSpace) override
  {
    ProfilingBlock prof("DeviceField::sync_to_host() for " + m_hostField->name());
    if (m_hostField->need_sync_to_host()) {
      set_execution_space(execSpace);
      if (not m_hostField->has_unified_device_storage()) {
        m_deviceFieldData.sync_to_host(get_execution_space(), m_hostField->host_data_layout());
      }
      else {
        execSpace.fence();
      }
      m_hostField->increment_num_syncs_to_host();
      m_hostField->clear_device_sync_state();
    }
  }

  void sync_to_host(ExecSpace&& execSpace) override
  {
    sync_to_host(execSpace);
  }

  void sync_to_device() override
  {
    sync_to_device(stk::ngp::ExecSpace());
  }

  void sync_to_device(const ExecSpace& execSpace) override
  {
    ProfilingBlock prof("DeviceField::sync_to_device() for " + m_hostField->name());
    if (m_hostField->need_sync_to_device()) {
      set_execution_space(execSpace);
      if (m_deviceFieldData.needs_update()) {
        m_deviceFieldData.update(get_execution_space(), m_hostField->host_data_layout(), need_sync_to_device());
      }
      if (not m_hostField->has_unified_device_storage()) {
        m_deviceFieldData.sync_to_device(get_execution_space(), m_hostField->host_data_layout());
      }
      else {
        execSpace.fence();
      }
      m_hostField->increment_num_syncs_to_device();
      m_hostField->clear_host_sync_state();
    }
  }

  void sync_to_device(ExecSpace&& execSpace) override
  {
    sync_to_device(execSpace);
  }

  size_t synchronized_count() const override {
    return m_deviceFieldData.field_data_synchronized_count();
  }

  KOKKOS_FUNCTION
  unsigned get_component_stride(const FastMeshIndex& entityIndex) const
  {
    if constexpr (impl::is_on_host()) {
      return m_deviceFieldData.mesh().buckets(m_deviceFieldData.entity_rank())[entityIndex.bucket_id]->capacity();
    }
    else {
      return m_deviceFieldData.m_deviceFieldMetaData(entityIndex.bucket_id).m_bucketCapacity;
    }
    return 0;  // Keep Nvidia compiler happy about always having return value
  }

  KOKKOS_FUNCTION
  unsigned get_component_stride(unsigned bucketId) const
  {
    if constexpr (impl::is_on_host()) {
      return m_deviceFieldData.mesh().buckets(m_deviceFieldData.entity_rank())[bucketId]->capacity();
    }
    else {
      return m_deviceFieldData.m_deviceFieldMetaData(bucketId).m_bucketCapacity;
    }
    return 0;  // Keep Nvidia compiler happy about always having return value
  }

  KOKKOS_FUNCTION
  unsigned get_num_components_per_entity(const FastMeshIndex& entityIndex) const {
    const DeviceFieldMetaData& fieldMetaData = m_deviceFieldData.m_deviceFieldMetaData(entityIndex.bucket_id);
    return fieldMetaData.m_numComponentsPerEntity*fieldMetaData.m_numCopiesPerEntity;
  }

  KOKKOS_FUNCTION
  unsigned get_extent0_per_entity(const FastMeshIndex& entityIndex) const {
    return m_deviceFieldData.m_deviceFieldMetaData(entityIndex.bucket_id).m_numComponentsPerEntity;
  }

  KOKKOS_FUNCTION
  unsigned get_extent1_per_entity(const FastMeshIndex& entityIndex) const {
    return m_deviceFieldData.m_deviceFieldMetaData(entityIndex.bucket_id).m_numCopiesPerEntity;
  }

  KOKKOS_FUNCTION
  unsigned get_extent_per_entity(const FastMeshIndex& entityIndex, unsigned dimension) const {
    const unsigned bucketId = entityIndex.bucket_id;
    if (dimension == 0) {
      return m_deviceFieldData.m_deviceFieldMetaData(bucketId).m_numComponentsPerEntity;
    }
    else if (dimension == 1) {
      return m_deviceFieldData.m_deviceFieldMetaData(entityIndex.bucket_id).m_numCopiesPerEntity;
    }
    else {
      const unsigned numComponents = m_deviceFieldData.m_deviceFieldMetaData(bucketId).m_numComponentsPerEntity;
      return (numComponents != 0) ? 1 : 0;
    }
  }

  // NgpMesh no longer needed to call this, but keeping the API the same because this will eventually
  // be deprecated in place of a totally new Field data access API.
  template <typename Mesh> KOKKOS_FUNCTION
  T& get(const Mesh& /*ngpMesh*/, Entity entity, int component) const
  {
    FastMeshIndex fastMeshIndex = m_deviceFieldData.m_deviceFastMeshIndices[entity.local_offset()];
    return get(fastMeshIndex, component);
  }

  KOKKOS_FUNCTION
  T& get(FastMeshIndex index, int component) const
  {
    const DeviceFieldMetaData& fieldMetaData = m_deviceFieldData.m_deviceFieldMetaData(index.bucket_id);
    return *(reinterpret_cast<T*>(fieldMetaData.m_data) + component*fieldMetaData.m_bucketCapacity + index.bucket_ord);
  }

  KOKKOS_FUNCTION
  T& operator()(const FastMeshIndex& index, int component) const
  {
    const DeviceFieldMetaData& fieldMetaData = m_deviceFieldData.m_deviceFieldMetaData(index.bucket_id);
    return *(reinterpret_cast<T*>(fieldMetaData.m_data) + component*fieldMetaData.m_bucketCapacity + index.bucket_ord);
  }

  KOKKOS_FUNCTION
  EntityFieldData<T> operator()(const FastMeshIndex& index) const
  {
    const DeviceFieldMetaData& fieldMetaData = m_deviceFieldData.m_deviceFieldMetaData(index.bucket_id);
    T* entityPtr = reinterpret_cast<T*>(fieldMetaData.m_data) + index.bucket_ord;
    return EntityFieldData<T>(entityPtr, fieldMetaData.m_numComponentsPerEntity, fieldMetaData.m_bucketCapacity);
  }

  template <typename Mesh>
  void set_all(const Mesh& /*ngpMesh*/, const T& value)
  {
    clear_sync_state();
    impl::fill_field_with_value<T>(get_execution_space(), m_deviceFieldData.m_deviceFieldMetaData, value);
    modify_on_device();
  }

  KOKKOS_FUNCTION
  EntityRank get_rank() const { return m_deviceFieldData.entity_rank(); }

  KOKKOS_FUNCTION
  unsigned get_ordinal() const { return m_deviceFieldData.field_ordinal(); }

  const std::string& name() const { return m_hostField->name(); }

  const BulkData& get_bulk() const { return m_deviceFieldData.mesh(); }

  FieldState state() const { return m_hostField->state(); }

  const FieldBase* get_field_base() const { return m_hostField; }

  void swap(DeviceField<T, NgpMemSpace>& otherDeviceField)
  {
    m_deviceFieldData.incomplete_swap_field_data(otherDeviceField.m_deviceFieldData);
  }

  bool need_sync_to_host() const override
  {
    return m_hostField->need_sync_to_host();
  }

  bool need_sync_to_device() const override
  {
    return m_hostField->need_sync_to_device();
  }

private:
  ExecSpace& get_execution_space() const { return m_hostField->get_execution_space(); }

  void set_execution_space(const ExecSpace& executionSpace)
  {
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, NgpMemSpace>::accessible);
    m_hostField->set_execution_space(executionSpace);
  }

  void set_execution_space(ExecSpace&& executionSpace)
  {
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, NgpMemSpace>::accessible);
    m_hostField->set_execution_space(std::forward<ExecSpace>(executionSpace));
  }

  void reset_execution_space() { m_hostField->reset_execution_space(); }

  void set_modify_on_host() { m_hostField->modify_on_host(); }

  void set_modify_on_device() { m_hostField->modify_on_device(); }

  void clear_sync_state_flags() { m_hostField->clear_sync_state(); }

  void update_field()
  {
    m_deviceFieldData = m_hostField->data<T, Unsynchronized, stk::ngp::DeviceSpace>();  // Reacquire current copy
  }

  const FieldBase* m_hostField;
  FieldData<T, stk::ngp::DeviceSpace> m_deviceFieldData;
};

}

#endif

