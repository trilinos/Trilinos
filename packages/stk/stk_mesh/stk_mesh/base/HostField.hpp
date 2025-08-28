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

#ifndef STK_MESH_HOST_FIELD_HPP
#define STK_MESH_HOST_FIELD_HPP

#include "stk_util/stk_config.h"
#include "Kokkos_Core.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/FieldData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/EntityFieldData.hpp"
#include "stk_mesh/baseImpl/NgpFieldAux.hpp"

namespace stk {
namespace mesh {

template<typename T, typename NgpMemSpace>
class HostField : public NgpFieldBase
{
 public:
  using ExecSpace = stk::ngp::ExecSpace;
  using MemSpace = NgpMemSpace;
  using value_type = T;

  HostField()
    : NgpFieldBase(),
      m_hostField(nullptr)
  {}

  HostField(const stk::mesh::BulkData& /*bulk*/, const stk::mesh::FieldBase& f,
            [[maybe_unused]] bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      m_hostField(&f),
      m_hostFieldData(m_hostField->data<T, Unsynchronized, stk::ngp::HostMemSpace>())
  {}

  HostField(const stk::mesh::FieldBase& f, [[maybe_unused]] bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      m_hostField(&f),
      m_hostFieldData(m_hostField->data<T, Unsynchronized, stk::ngp::HostMemSpace>())
  {}

  HostField(const HostField&) = default;
  HostField(HostField&&) = default;
  HostField& operator=(const HostField&) = default;
  HostField& operator=(HostField&&) = default;

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
  void fence() override {}

  unsigned get_component_stride(const FastMeshIndex& /*entityIndex*/) const { return 1; }

  unsigned get_component_stride(unsigned /*bucketId*/) const { return 1; }

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after March 5, 2025
  STK_DEPRECATED_MSG("The component stride on device is now a function of the Bucket you are accessing.  Please use "
                     "one of the other overloads of get_component_stride().")
  unsigned get_component_stride() const { return 1; }
#endif

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_scalars_per_entity(*m_hostField, entity.bucket_id);
  }

  unsigned get_extent0_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_extent0_per_entity(*m_hostField, entity.bucket_id);
  }

  unsigned get_extent1_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_extent1_per_entity(*m_hostField, entity.bucket_id);
  }

  unsigned get_extent_per_entity(const stk::mesh::FastMeshIndex& entity, unsigned dimension) const {
    return stk::mesh::field_extent_per_entity(*m_hostField, dimension, entity.bucket_id);
  }

  T& get(const HostMesh& /*ngpMesh*/, stk::mesh::Entity entity, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*m_hostField, entity));
    STK_ThrowAssert(data);
    return data[component];
  }

  T& get(stk::mesh::FastMeshIndex entity, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*m_hostField, entity.bucket_id, entity.bucket_ord));
    STK_ThrowAssert(data);
    return data[component];
  }

  T& operator()(const stk::mesh::FastMeshIndex& index, int component) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*m_hostField, index.bucket_id, index.bucket_ord));
    STK_ThrowAssert(data);
    return data[component];
  }

  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*m_hostField, index.bucket_id, index.bucket_ord));
    unsigned numScalars = stk::mesh::field_scalars_per_entity(*m_hostField, index.bucket_id);
    STK_ThrowAssert(data);
    return EntityFieldData<T>(data, numScalars, 1);
  }

  void set_all(const HostMesh& ngpMesh, const T& value)
  {
    stk::mesh::for_each_entity_run(ngpMesh, m_hostField->entity_rank(), *m_hostField, [&](const FastMeshIndex& entity) {
      T* fieldPtr = static_cast<T*>(stk::mesh::field_data(*m_hostField, entity.bucket_id, entity.bucket_ord));
      int numScalars = stk::mesh::field_scalars_per_entity(*m_hostField, entity.bucket_id);
      for (int i = 0; i < numScalars; i++) {
        fieldPtr[i] = value;
      }
    });
  }

  void modify_on_host() override { m_hostField->modify_on_host(); }

  void modify_on_device() override { m_hostField->modify_on_device(); }

  void clear_sync_state() override
  {
    m_hostField->clear_sync_state();
  }

  void modify_on_host(const Selector& selector) override { m_hostField->modify_on_host(selector); }

  void modify_on_device(const Selector& selector) override { m_hostField->modify_on_device(selector); }

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
    m_hostField->sync_to_host();
    if constexpr (!std::is_same_v<Kokkos::DefaultExecutionSpace,Kokkos::Serial>) {
      Kokkos::fence();
    }
  }

  void sync_to_host(const ExecSpace& execSpace) override
  {
    m_hostField->sync_to_host(execSpace);
  }

  void sync_to_host(ExecSpace&& execSpace) override
  {
    m_hostField->sync_to_host(execSpace);
  }

  void sync_to_device() override 
  {
    m_hostField->sync_to_device();
    Kokkos::fence();
  }

  void sync_to_device(const ExecSpace& execSpace) override
  {
    m_hostField->sync_to_device(execSpace);
  }

  void sync_to_device(ExecSpace&& execSpace) override
  {
    sync_to_device(execSpace);
  }

  size_t synchronized_count() const override { return m_hostFieldData.field_data_synchronized_count(); }

  FieldState state() const { return m_hostField->state(); }

  stk::mesh::EntityRank get_rank() const { return m_hostFieldData.m_rank; }

  unsigned get_ordinal() const { return m_hostFieldData.m_ordinal; }

  const std::string& name() const { return m_hostField->name(); }

  const BulkData& get_bulk() const { return m_hostFieldData.mesh(); }

  const FieldBase* get_field_base() const { return m_hostField; }

  bool need_sync_to_host() const override { return m_hostField->need_sync_to_host(); }

  bool need_sync_to_device() const override { return m_hostField->need_sync_to_device(); }

 private:
  ExecSpace& get_execution_space() const { return m_hostField->get_execution_space(); }

  void set_execution_space(const ExecSpace& executionSpace) { m_hostField->set_execution_space(executionSpace); }

  void set_execution_space(ExecSpace&& executionSpace)
  {
    m_hostField->set_execution_space(std::forward<ExecSpace>(executionSpace));
  }

  void reset_execution_space() { m_hostField->reset_execution_space(); }

  void update_field()
  {
    if (m_hostFieldData.needs_update()) {
      m_hostFieldData.update(get_execution_space(), m_hostField->host_data_layout());
      m_hostField->increment_num_syncs_to_device();
    }
  }

  void set_modify_on_device() {
    m_hostField->modify_on_device();
  }

  void copy_device_to_host()
  {
    m_hostField->increment_num_syncs_to_host();
  }

  const stk::mesh::FieldBase* m_hostField;
  FieldData<T, stk::ngp::HostMemSpace> m_hostFieldData;
};

}
}

#endif

