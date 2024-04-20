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

#ifndef NGPFIELDSYNCDEBUGGER_HPP
#define NGPFIELDSYNCDEBUGGER_HPP

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/DeviceMesh.hpp"
#include "stk_mesh/base/NgpFieldBase.hpp"
#include "stk_mesh/base/FieldSyncDebugging.hpp"
#include "stk_mesh/base/StkFieldSyncDebugger.hpp"
#include "stk_mesh/base/NgpForEachEntity.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include <memory>

namespace stk {
namespace mesh {

//==============================================================================
template <typename T>
class EmptyNgpFieldSyncDebugger
{
public:
  using StkFieldSyncDebuggerType = EmptyStkFieldSyncDebugger;

  KOKKOS_FUNCTION EmptyNgpFieldSyncDebugger(StkFieldSyncDebuggerType *) {}
  KOKKOS_DEFAULTED_FUNCTION ~EmptyNgpFieldSyncDebugger() = default;

  KOKKOS_DEFAULTED_FUNCTION EmptyNgpFieldSyncDebugger(const EmptyNgpFieldSyncDebugger &) = default;
  KOKKOS_DEFAULTED_FUNCTION EmptyNgpFieldSyncDebugger& operator=(const EmptyNgpFieldSyncDebugger &) = default;

  template <typename NgpField>
  inline void initialize_debug_views(NgpField*) {}

  template <typename NgpField>
  inline void update_field(NgpField*) {}

  template <typename NgpField>
  inline void clear_sync_state(NgpField*) {}

  template <typename NgpField>
  inline void clear_host_sync_state(NgpField*) {}

  template <typename NgpField>
  inline void clear_device_sync_state(NgpField*) {}

  template <typename NgpField>
  inline void sync_to_host(NgpField*) {}

  template <typename NgpField>
  inline void sync_to_device(NgpField*) {}

  template <typename ViewType>
  inline void initialize_view(ViewType&) {}

  inline void set_any_potential_device_field_modification(bool) {}

  inline unsigned get_bucket_offset(unsigned) const { return 0; }

  template <typename NgpField>
  inline void update_debug_storage(NgpField*, size_t) {}

  template <typename NgpField>
  inline void detect_device_field_modification(NgpField*) {}

  template <typename NgpField>
  KOKKOS_INLINE_FUNCTION
  void device_stale_access_check(NgpField *, const stk::mesh::FastMeshIndex &, int, const char *, int) const {}

  template <typename NgpField>
  KOKKOS_INLINE_FUNCTION
  void device_stale_access_check(NgpField *, const stk::mesh::FastMeshIndex &, const char *, int) const {}

  template <typename NgpField>
  void device_stale_access_check(NgpField *, const char *, int) const {}
};

//==============================================================================
template <typename T>
class NgpFieldSyncDebugger
{
public:
  using StkFieldSyncDebuggerType = StkFieldSyncDebugger;

  KOKKOS_FUNCTION NgpFieldSyncDebugger(StkFieldSyncDebuggerType * stkFieldSyncDebugger)
    : m_stkFieldSyncDebugger(stkFieldSyncDebugger)
  {}

  KOKKOS_DEFAULTED_FUNCTION ~NgpFieldSyncDebugger() = default;

  KOKKOS_DEFAULTED_FUNCTION NgpFieldSyncDebugger(const NgpFieldSyncDebugger &) = default;
  KOKKOS_DEFAULTED_FUNCTION NgpFieldSyncDebugger& operator=(const NgpFieldSyncDebugger &) = default;

  template <typename NgpField>
  void initialize_debug_views(NgpField* ngpField)
  {
    const stk::mesh::FieldBase & stkField = *ngpField->hostField;
    const stk::mesh::BulkData & bulk = *ngpField->hostBulk;
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    const stk::mesh::BucketVector& allBuckets = bulk.buckets(stkField.entity_rank());
    const stk::mesh::BucketVector& buckets = bulk.get_buckets(stkField.entity_rank(), selector);
    unsigned numPerEntity = stkField.max_size();

    construct_debug_views(ngpField, allBuckets, buckets, numPerEntity);

    std::strcpy(fieldName.data(), stkField.name().data());
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;

    Kokkos::deep_copy(lastFieldValue, ngpField->deviceData);

    reset_last_modification_state(ngpField, LastModLocation::HOST_OR_DEVICE);

    hostSynchronizedCount() = bulk.synchronized_count();
    isFieldLayoutConsistent = true;
  }

  template <typename NgpField>
  void update_field(NgpField* ngpField)
  {
    if (lastFieldValue.extent(0) != 0) {
      lostDeviceFieldData() = false;
      anyPotentialDeviceFieldModification() = false;

      Kokkos::deep_copy(lastFieldValue, ngpField->deviceData);

      reset_last_modification_state(ngpField, LastModLocation::HOST_OR_DEVICE);
      hostSynchronizedCount() = ngpField->hostBulk->synchronized_count();
    }

    isFieldLayoutConsistent = true;
  }

  template <typename NgpField>
  void clear_sync_state(NgpField* ngpField)
  {
    clear_host_sync_state(ngpField);
    clear_device_sync_state(ngpField);
  }

  template <typename NgpField>
  void clear_host_sync_state(NgpField* ngpField)
  {
    set_last_modification_state_bit(LastModLocation::DEVICE);
  }

  template <typename NgpField>
  void clear_device_sync_state(NgpField* ngpField)
  {
    if (ngpField->hostBulk->synchronized_count() != ngpField->synchronizedCount) {
      ngpField->update_field();
    }
    Kokkos::deep_copy(lastFieldValue, ngpField->deviceData);
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;
    set_last_modification_state_bit(LastModLocation::HOST);
  }

  template <typename NgpField>
  void sync_to_host(NgpField* ngpField)
  {
    set_last_device_field_value(ngpField, get_modified_selector(ngpField));
    reset_last_modification_state(ngpField, LastModLocation::HOST_OR_DEVICE);
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;
  }

  template <typename NgpField>
  void sync_to_device(NgpField* ngpField)
  {
    Kokkos::deep_copy(lastFieldValue, ngpField->deviceData);
    reset_last_modification_state(ngpField, LastModLocation::HOST_OR_DEVICE);
    lostDeviceFieldData() = false;
    anyPotentialDeviceFieldModification() = false;
  }

  template <typename NgpField>
  KOKKOS_FUNCTION
  void device_stale_access_check(NgpField* ngpField, const stk::mesh::FastMeshIndex& index, int component,
                                 const char* fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod(ngpField->synchronizedCount)) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    if (data_is_stale_on_device(index, component)) {
      print_stale_data_warning(ngpField, index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
    }
  }

  template <typename NgpField>
  KOKKOS_FUNCTION
  void device_stale_access_check(NgpField* ngpField, const stk::mesh::FastMeshIndex& index,
                                 const char* fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod(ngpField->synchronizedCount)) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    for (unsigned component = 0; component < ngpField->get_num_components_per_entity(index); ++component) {
      if (data_is_stale_on_device(index, component)) {
        print_stale_data_warning(ngpField, index.bucket_id, index.bucket_ord, component, fileName, lineNumber);
      }
    }
  }

  template <typename NgpField>
  void device_stale_access_check(NgpField* ngpField,
                                 const char* fileName, int lineNumber) const
  {
    anyPotentialDeviceFieldModification() = true;

    if (field_not_updated_after_mesh_mod(ngpField->synchronizedCount)) {
      print_unupdated_field_warning(fileName, lineNumber);
      return;
    }

    for (size_t i = 0; i < lastFieldModLocation.extent(0); ++i) {
      for (size_t j = 0; j < lastFieldModLocation.extent(1); ++j) {
        for (size_t k = 0; k < lastFieldModLocation.extent(2); ++k) {
          if (!(lastFieldModLocation(i,j,k) & LastModLocation::DEVICE)) {
            print_stale_data_warning_without_field_values(i, ORDER_INDICES(j, k), fileName, lineNumber);
          }
        }
      }
    }
  }

  template <typename ViewType>
  void initialize_view(ViewType& view)
  {
    Kokkos::deep_copy(view, 0);
  }

  void set_any_potential_device_field_modification(bool value)
  {
    anyPotentialDeviceFieldModification() = value;
  }

  template <typename NgpField>
  void detect_device_field_modification(NgpField* ngpField) {
    if (!anyPotentialDeviceFieldModification()) {
      return;
    }

    if (field_not_updated_after_mesh_mod(ngpField->synchronizedCount)) {
      return;
    }

    const stk::mesh::BulkData & bulk = *ngpField->hostBulk;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    stk::mesh::Selector fieldSelector(*(ngpField->hostField));

    UnsignedViewType & localDeviceNumComponentsPerEntity = ngpField->deviceFieldBucketsNumComponentsPerEntity;
    FieldDataDeviceViewType<T> & localDeviceData = ngpField->deviceData;
    FieldDataDeviceViewType<T> & localLastFieldValue = lastFieldValue;
    LastFieldModLocationType & localLastFieldModLocation = lastFieldModLocation;
    ScalarUvmType<bool> & localLostDeviceFieldData = lostDeviceFieldData;
    UnsignedViewType & localDebugDeviceSelectedBucketOffset = debugDeviceSelectedBucketOffset;
    const bool inModCycle = bulk.in_modifiable_state();

    stk::mesh::for_each_entity_run(ngpMesh, ngpField->rank, fieldSelector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& index)
    {
      const unsigned offsetBucketId = localDebugDeviceSelectedBucketOffset(index.bucket_id);
      const unsigned numComponents = localDeviceNumComponentsPerEntity(offsetBucketId);
      for (unsigned component = 0; component < numComponents; ++component) {
        if (localDeviceData(offsetBucketId, ORDER_INDICES(index.bucket_ord, component)) !=
            localLastFieldValue(offsetBucketId, ORDER_INDICES(index.bucket_ord, component)))
        {
          localLastFieldModLocation(offsetBucketId, ORDER_INDICES(index.bucket_ord, component)) = LastModLocation::DEVICE;
          if (inModCycle) {
            localLostDeviceFieldData() = true;
          }
        }
      }
    });

    anyPotentialDeviceFieldModification() = false;
    Kokkos::fence();
  }

  template <typename NgpField>
  void update_debug_storage(NgpField* ngpField, size_t hostSyncCount) {
    hostSynchronizedCount() = hostSyncCount;

    const stk::mesh::FieldBase & stkField = *ngpField->hostField;
    const stk::mesh::BulkData & bulk = *ngpField->hostBulk;
    stk::mesh::Selector selector = stk::mesh::selectField(stkField);
    const stk::mesh::BucketVector & allBuckets = bulk.buckets(stkField.entity_rank());
    const stk::mesh::BucketVector & buckets = bulk.get_buckets(stkField.entity_rank(), selector);
    const unsigned numPerEntity = stkField.max_size();

    ngpField->construct_new_index_view(allBuckets);
    construct_debug_views(ngpField, allBuckets, buckets, numPerEntity);

    isFieldLayoutConsistent = false;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned get_bucket_offset(unsigned bucketOrdinal) const {
    return debugHostSelectedBucketOffset(bucketOrdinal);
  }

  template <typename NgpField>
  void set_last_device_field_value(NgpField* ngpField, const stk::mesh::Selector& modifiedSelector)
  {
    const stk::mesh::BulkData & bulk = *ngpField->hostBulk;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    UnsignedViewType & localDeviceNumComponentsPerEntity = ngpField->deviceFieldBucketsNumComponentsPerEntity;
    FieldDataDeviceViewType<T> & localDeviceData = ngpField->deviceData;
    FieldDataDeviceViewType<T> & localLastFieldValue = lastFieldValue;
    UnsignedViewType & localDebugDeviceSelectedBucketOffset = debugDeviceSelectedBucketOffset;

    stk::mesh::for_each_entity_run(ngpMesh, ngpField->rank, modifiedSelector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& index)
    {
      const unsigned offsetBucketId = localDebugDeviceSelectedBucketOffset(index.bucket_id);
      const unsigned numComponents = localDeviceNumComponentsPerEntity(offsetBucketId);
      for (unsigned component = 0; component < numComponents; ++component) {
        localLastFieldValue(offsetBucketId, ORDER_INDICES(index.bucket_ord, component)) =
                localDeviceData(offsetBucketId, ORDER_INDICES(index.bucket_ord, component));
      }
    });
  }

private:
  template <typename NgpField>
  void construct_debug_views(NgpField* ngpField, const BucketVector & allBuckets, const BucketVector & buckets, unsigned numPerEntity) {
    const stk::mesh::FieldBase & stkField = *ngpField->hostField;

    if (buckets.size() != 0) {
      lastFieldValue = FieldDataDeviceViewType<T>(stkField.name()+"_lastValue", buckets.size(),
                                                  ORDER_INDICES(ngpField->bucketCapacity, numPerEntity));
      lastFieldModLocation = LastFieldModLocationType(stkField.name()+"_lastModLocation", buckets.size(),
                                                      ORDER_INDICES(ngpField->bucketCapacity, numPerEntity));
    }

    if (fieldName.extent(0) == 0u) {
      fieldName = FieldNameType(stkField.name()+"_name", stkField.name().size()+1);
      hostSynchronizedCount = ScalarUvmType<size_t>(stkField.name()+"_hostSyncCount");
      lostDeviceFieldData = ScalarUvmType<bool>(stkField.name()+"_anyDeviceFieldMod");
      anyPotentialDeviceFieldModification = ScalarUvmType<bool>(stkField.name()+"_anyPotentialDeviceFieldMod");
    }

    debugHostSelectedBucketOffset = ngpField->newHostSelectedBucketOffset;
    debugDeviceSelectedBucketOffset = ngpField->newDeviceSelectedBucketOffset;

    m_stkFieldSyncDebugger->set_last_modification_view(lastFieldModLocation);
    m_stkFieldSyncDebugger->set_lost_device_field_data_view(lostDeviceFieldData);
    m_stkFieldSyncDebugger->set_bucket_offset_view(debugHostSelectedBucketOffset);
    m_stkFieldSyncDebugger->mark_data_initialized();
  }

  template <typename NgpField>
  void reset_last_modification_state(NgpField* ngpField, LastModLocation value) {
    Kokkos::deep_copy(lastFieldModLocation, value);
  }

  void set_last_modification_state_bit(LastModLocation value) {
    for (size_t i = 0; i < lastFieldModLocation.extent(0); ++i) {
      for (size_t j = 0; j < lastFieldModLocation.extent(1); ++j) {
        for (size_t k = 0; k < lastFieldModLocation.extent(2); ++k) {
          lastFieldModLocation(i, j, k) = static_cast<LastModLocation>(lastFieldModLocation(i, j, k) | value);
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  bool field_not_updated_after_mesh_mod(size_t synchronizedCount) const {
    return hostSynchronizedCount() != synchronizedCount;
  }

  KOKKOS_INLINE_FUNCTION
  bool data_is_stale_on_device(const stk::mesh::FastMeshIndex & index, int component) const {
    return data_is_stale_on_device(index.bucket_id, index.bucket_ord, component);
  }

  KOKKOS_INLINE_FUNCTION
  bool data_is_stale_on_device(int bucketId, int bucketOrdinal, int component) const {
    return !(lastFieldModLocation(debugDeviceSelectedBucketOffset(bucketId),
                                  ORDER_INDICES(bucketOrdinal, component)) & LastModLocation::DEVICE);
  }

  void print_outside_selector_warning(const Selector& badSelector) const
  {
    std::cout << "*** WARNING: Marked field " << fieldName.data() << " modified with selector " << badSelector
      << " that includes buckets outside the subset of the mesh that the field is defined on." << std::endl;
  }

  KOKKOS_INLINE_FUNCTION
  void print_unupdated_field_warning(const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing un-updated Field %s on Device after mesh modification\n", fieldName.data());
    }
    else {
      printf("%s:%i *** WARNING: Accessing un-updated Field %s on Device after mesh modification\n",
             fileName, lineNumber, fieldName.data());
    }
  }

  template <typename NgpField>
  KOKKOS_INLINE_FUNCTION
  void print_stale_data_warning(NgpField* ngpField, unsigned bucketId, unsigned bucketOrdinal, int component,
                                const char * fileName, int lineNumber) const
  {
    const double currentValue = ngpField->deviceData(debugDeviceSelectedBucketOffset(bucketId),
                                                     ORDER_INDICES(bucketOrdinal, component));
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing stale data on Device for Field %s[%i]=%f\n",
             fieldName.data(), component, currentValue);
    }
    else {
      printf("%s:%i *** WARNING: Accessing stale data on Device for Field %s[%i]=%f\n",
             fileName, lineNumber, fieldName.data(), component, currentValue);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void print_stale_data_warning_without_field_values(unsigned bucketId, unsigned bucketOrdinal, int component,
                                                     const char * fileName, int lineNumber) const
  {
    if (lineNumber == -1) {
      printf("*** WARNING: Accessing stale data on Device for Field %s[%i]\n", fieldName.data(), component);
    }
    else {
      printf("%s:%i *** WARNING: Accessing stale data on Device for Field %s[%i]\n",
             fileName, lineNumber, fieldName.data(), component);
    }
  }

  template <typename NgpField>
  stk::mesh::Selector get_modified_selector(NgpField* ngpField)
  {
    stk::mesh::Selector modifiedSelector(*(ngpField->hostField));
    return modifiedSelector;
  }


  using FieldNameType = Kokkos::View<char*, stk::ngp::UVMMemSpace>;
  template <typename U> using HostArrayType = Kokkos::View<U*, Kokkos::HostSpace>;

  StkFieldSyncDebuggerType* m_stkFieldSyncDebugger;
  FieldNameType fieldName;
  bool isFieldLayoutConsistent;
  ScalarUvmType<size_t> hostSynchronizedCount;
  ScalarUvmType<bool> lostDeviceFieldData;
  ScalarUvmType<bool> anyPotentialDeviceFieldModification;
  LastFieldModLocationType lastFieldModLocation;
  FieldDataDeviceViewType<T> lastFieldValue;
  typename UnsignedViewType::HostMirror debugHostSelectedBucketOffset;
  UnsignedViewType debugDeviceSelectedBucketOffset;
};

//==============================================================================

}
}

#endif // NGPFIELDSYNCDEBUGGER_HPP
