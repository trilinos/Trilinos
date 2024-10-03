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
#include "Kokkos_DualView.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpProfilingBlock.hpp"
#include "stk_mesh/base/FieldSyncDebugging.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/EntityFieldData.hpp"
#include "stk_mesh/baseImpl/NgpFieldAux.hpp"
#include "stk_mesh/base/NgpFieldSyncDebugger.hpp"

namespace stk {
namespace mesh {

template<typename T, template <typename> class NgpDebugger>
class HostField : public NgpFieldBase
{
 public:
  using ExecSpace = stk::ngp::ExecSpace;
  using value_type = T;
  using StkDebugger = typename NgpDebugger<T>::StkFieldSyncDebuggerType;

  HostField()
    : NgpFieldBase(),
      hostBulk(nullptr),
      field(nullptr),
      synchronizedCount(0)
  {
  }

  HostField(const stk::mesh::BulkData& b, const stk::mesh::FieldBase& f, bool isFromGetUpdatedNgpField = false)
    : NgpFieldBase(),
      hostBulk(&b),
      field(&f),
      synchronizedCount(0)
  {
    field->template make_field_sync_debugger<StkDebugger>();
  }

  HostField(const HostField<T, NgpDebugger>&) = default;
  HostField(HostField<T, NgpDebugger>&&) = default;
  HostField<T, NgpDebugger>& operator=(const HostField<T, NgpDebugger>&) = default;
  HostField<T, NgpDebugger>& operator=(HostField<T, NgpDebugger>&&) = default;

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

  void set_field_states(HostField<T, NgpDebugger>* fields[]) {}

  size_t num_syncs_to_host() const override { return field->num_syncs_to_host(); }
  size_t num_syncs_to_device() const override { return field->num_syncs_to_device(); }
  void fence() override {}

  unsigned get_component_stride() const { return 1; }

  unsigned get_num_components_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
  }

  unsigned get_extent0_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_extent0_per_entity(*field, entity.bucket_id);
  }

  unsigned get_extent1_per_entity(const stk::mesh::FastMeshIndex& entity) const {
    return stk::mesh::field_extent1_per_entity(*field, entity.bucket_id);
  }

  unsigned get_extent_per_entity(const stk::mesh::FastMeshIndex& entity, unsigned dimension) const {
    return stk::mesh::field_extent_per_entity(*field, dimension, entity.bucket_id);
  }

  unsigned debug_get_bucket_offset(unsigned bucketOrdinal) const override {
    return bucketOrdinal;
  }

  T& get(const HostMesh& ngpMesh, stk::mesh::Entity entity, int component,
         const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity));
    STK_ThrowAssert(data);
    return data[component];
  }

  T& get(stk::mesh::FastMeshIndex entity, int component,
         const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, entity.bucket_id, entity.bucket_ord));
    STK_ThrowAssert(data);
    return data[component];
  }

  T& operator()(const stk::mesh::FastMeshIndex& index, int component,
                const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    STK_ThrowAssert(data);
    return data[component];
  }

  EntityFieldData<T> operator()(const stk::mesh::FastMeshIndex& index,
                                const char * fileName = HOST_DEBUG_FILE_NAME, int lineNumber = HOST_DEBUG_LINE_NUMBER) const
  {
    T *data = static_cast<T *>(stk::mesh::field_data(*field, index.bucket_id, index.bucket_ord));
    unsigned numScalars = stk::mesh::field_scalars_per_entity(*field, index.bucket_id);
    STK_ThrowAssert(data);
    return EntityFieldData<T>(data, numScalars);
  }

  void set_all(const HostMesh& ngpMesh, const T& value)
  {
    stk::mesh::for_each_entity_run(ngpMesh, field->entity_rank(), *field, [&](const FastMeshIndex& entity) {
      T* fieldPtr = static_cast<T*>(stk::mesh::field_data(*field, entity.bucket_id, entity.bucket_ord));
      int numScalars = stk::mesh::field_scalars_per_entity(*field, entity.bucket_id);
      for (int i = 0; i < numScalars; i++) {
        fieldPtr[i] = value;
      }
    });
  }

  void modify_on_host() override { field->modify_on_host(); }

  void modify_on_device() override { field->modify_on_device(); }

  void clear_sync_state() override
  {
    field->clear_sync_state();
  }

  void modify_on_host(const Selector& selector) override { field->modify_on_host(selector); }

  void modify_on_device(const Selector& selector) override { field->modify_on_device(selector); }

  void clear_host_sync_state() override
  {
    field->clear_host_sync_state();
  }

  void clear_device_sync_state() override
  {
    field->clear_device_sync_state();
  }

  void sync_to_host() override
  {
    sync_to_host(Kokkos::DefaultExecutionSpace());
    Kokkos::fence();
  }

  void sync_to_host(const ExecSpace& execSpace) override
  {
    if (need_sync_to_host()) {
      copy_device_to_host();
      clear_device_sync_state();
    }
  }

  void sync_to_host(ExecSpace&& execSpace) override
  {
    if (need_sync_to_host()) {
      copy_device_to_host();
      clear_device_sync_state();
    }
  }

  void sync_to_device() override 
  {
    sync_to_device(Kokkos::DefaultExecutionSpace());
    Kokkos::fence();
  }

  void sync_to_device(const ExecSpace& execSpace) override
  {
    if (need_sync_to_device()) {
      if (hostBulk->synchronized_count() != synchronizedCount) {
        update_field();
      }
      copy_host_to_device();
      clear_host_sync_state();
    }
  }

  void sync_to_device(ExecSpace&& execSpace) override
  {
    if (need_sync_to_device()) {
      if (hostBulk->synchronized_count() != synchronizedCount) {
        update_field();
      }
      copy_host_to_device();
      clear_host_sync_state();
    }
  }

  size_t synchronized_count() const override { return synchronizedCount; }

  FieldState state() const { return field->state(); }

  void update_bucket_pointer_view() override { }

  void swap_field_views(NgpFieldBase *other) override { }
  void swap(HostField<T> &other) { }

  stk::mesh::EntityRank get_rank() const { return field ? field->entity_rank() : stk::topology::INVALID_RANK; }

  unsigned get_ordinal() const { return field->mesh_meta_data_ordinal(); }

  const FieldBase* get_field_base() const { return field; }

  void debug_initialize_debug_views() override {}
  void debug_modification_begin() override {}
  void debug_modification_end(size_t) override {}
  void debug_detect_device_field_modification() override {}

  bool need_sync_to_host() const override { return field->need_sync_to_host(); }

  bool need_sync_to_device() const override { return field->need_sync_to_device(); }

  void notify_sync_debugger_clear_sync_state() override {}
  void notify_sync_debugger_clear_host_sync_state() override {}
  void notify_sync_debugger_clear_device_sync_state() override {}

 private:
  ExecSpace& get_execution_space() const { return field->get_execution_space(); }

  void set_execution_space(const ExecSpace& executionSpace) { field->set_execution_space(executionSpace); }

  void set_execution_space(ExecSpace&& executionSpace)
  {
    field->set_execution_space(std::forward<ExecSpace>(executionSpace));
  }

  void reset_execution_space() { field->reset_execution_space(); }

  void update_field()
  {
    field->increment_num_syncs_to_device();
    synchronizedCount = hostBulk->synchronized_count();
  }

  void set_modify_on_device() {
    field->modify_on_device();
  }

  void copy_host_to_device()
  {
    field->increment_num_syncs_to_device();
  }

  void copy_device_to_host()
  {
    field->increment_num_syncs_to_host();
  }

  const BulkData* hostBulk;
  const stk::mesh::FieldBase * field;

  size_t synchronizedCount;
};

}
}

#endif

