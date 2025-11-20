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

#ifndef STK_MESH_DEVICEFIELDDATAMANAGER_HPP
#define STK_MESH_DEVICEFIELDDATAMANAGER_HPP

#include "stk_mesh/base/DeviceFieldDataManagerBase.hpp"
#include "stk_util/stk_config.h"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/AdjustForAlignment.hpp"
#include "stk_util/util/FieldDataAllocator.hpp"
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/FieldDataBytes.hpp"
#include "stk_mesh/base/FindRestriction.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include <vector>
#include <algorithm>

namespace stk::mesh {

// All device-side allocations are automatically aligned on 128 or 256-byte boundaries,
// we don't have SIMD to worry about, and only alignment of the fundamental datatype
// matters for fast access.  Just use the biggest fundamental size for our
// intra-Bucket alignment requirements.
constexpr int DeviceFieldAlignmentSize = alignof(std::max_align_t);

template <typename Space>
class DeviceFieldDataManager : public DeviceFieldDataManagerBase
{
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;

  using DeviceFieldMetaDataCollectionType = Kokkos::View<DeviceFieldMetaDataArrayType<mem_space>*, stk::ngp::HostExecSpace>;
  using HostFieldMetaDataCollectionType = Kokkos::View<HostFieldMetaDataArrayType<mem_space>*, stk::ngp::HostExecSpace>;

  using DeviceBucketsModifiedCollectionType = Kokkos::View<int**, Kokkos::LayoutRight, mem_space>;
  using HostBucketsModifiedCollectionType = typename DeviceBucketsModifiedCollectionType::host_mirror_type;

  using AllocationType = FieldDataAllocator<std::byte>::DeviceAllocationType;
  using BucketRawDataArrayType = Kokkos::View<AllocationType*, stk::ngp::HostExecSpace>;
  using BucketCapacityType = std::vector<size_t>;
  using FieldArrayType = std::array<std::vector<FieldBase*>, stk::topology::NUM_RANKS>;

  struct BucketShift {
    BucketShift(int _oldIndex, int _newIndex)
      : oldIndex(_oldIndex),
        newIndex(_newIndex)
    {}
    int oldIndex;
    int newIndex;
  };

public:
  DeviceFieldDataManager(const BulkData& bulk)
    : DeviceFieldDataManagerBase(),
      m_bulk(bulk),
      m_synchronizedCount(0),
      m_totalNumFields(0)
  {}

  virtual ~DeviceFieldDataManager() override = default;

  // Reorganize the device-side Field data allocations to match changes on the host side.
  // No effort is made to initialize any storage.  Sections of Fields in Buckets that have
  // either had their contents or their allocation modifed will be synced in-full to device
  // during an update step.
  virtual bool update_all_bucket_allocations() override;
  virtual void update_host_bucket_pointers(Ordinal fieldOrdinal) override;
  virtual void swap_field_data(Ordinal fieldOrdinal1, Ordinal fieldOrdinal2) override;
  virtual void clear_bucket_is_modified(Ordinal fieldOrdinal) override;
  virtual std::any get_device_bucket_is_modified(Ordinal fieldOrdinal, int& fieldIndex) override;
  virtual size_t get_num_bytes_allocated_on_field(const FieldBase& field) const override;
  virtual bool has_unified_device_storage(Ordinal fieldOrdinal) const override;

  virtual void set_device_field_meta_data(FieldDataBase& fieldDataBase) override;

private:
  void allocate_bucket(EntityRank rank, const FieldVector& fields, const PartVector& parts,
                       int bucketId, int size, int capacity);
  void resize_field_arrays(int oldNumAllFields, int newNumAllFields);
  void resize_field_meta_data_arrays(const FieldVector& fields, int oldNumBuckets, int newNumBuckets);
  void resize_bucket_arrays(EntityRank rank, int oldNumBuckets, int newNumBuckets);
  void resize_bucket_modified_array(EntityRank rank, int newNumFields, int newNumBuckets);
  void reorder_buckets(EntityRank rank);
  void shift_bucket_data(EntityRank rank, const std::vector<BucketShift>& bucketShiftList);
  void shift_field_and_bucket_data(EntityRank rank, const FieldVector& fieldsOfRank,
                                   const std::vector<BucketShift>& bucketShiftList);
  void update_field_meta_data(const FieldVector& fields, int bucketId, int bucketSize);
  void fill_field_meta_data_pointers_from_offsets(int bucketId, const FieldVector& fields,
                                                  const AllocationType& bucketRawData,
                                                  HostFieldMetaDataCollectionType& hostFieldMetaData);
  std::byte* get_host_bucket_pointer_for_device(const stk::mesh::FieldBase& field, int bucketId) const;

  FieldDataAllocator<std::byte> m_fieldDataAllocator;
  const BulkData& m_bulk;
  int m_synchronizedCount;
  int m_totalNumFields;
  BucketRawDataArrayType m_bucketRawData[stk::topology::NUM_RANKS];
  DeviceFieldMetaDataCollectionType m_deviceFieldMetaData;
  HostFieldMetaDataCollectionType m_hostFieldMetaData;
  BucketCapacityType m_bucketCapacity[stk::topology::NUM_RANKS];
  DeviceBucketsModifiedCollectionType m_deviceBucketIsModified[stk::topology::NUM_RANKS];
  HostBucketsModifiedCollectionType m_hostBucketIsModified[stk::topology::NUM_RANKS];
};


template <typename Space>
std::byte* DeviceFieldDataManager<Space>::get_host_bucket_pointer_for_device(const stk::mesh::FieldBase& field,
                                                                             int bucketId) const
{
  return m_fieldDataAllocator.get_host_pointer_for_device(field.get_meta_data_for_field()[bucketId].m_data);
}

template <typename Space>
bool DeviceFieldDataManager<Space>::update_all_bucket_allocations()
{
  if (static_cast<int>(m_bulk.synchronized_count()) == m_synchronizedCount) {
    return false;  // Nothing to do
  }

  const MetaData& meta = m_bulk.mesh_meta_data();
  const FieldVector allFields = meta.get_fields();
  const int oldNumAllFields = m_totalNumFields;
  const int newNumAllFields = allFields.size();
  FieldArrayType oldFields {};
  FieldArrayType newFields {};

  for (int i = 0; i < newNumAllFields; ++i) {
    FieldBase* field = allFields[i];
    const int fieldRank = field->entity_rank();
    if (i < oldNumAllFields) {
      oldFields[fieldRank].push_back(field);
    }
    else {
      newFields[fieldRank].push_back(field);
    }
  }

  if (newNumAllFields != oldNumAllFields) {
    resize_field_arrays(oldNumAllFields, newNumAllFields);
  }

  for (EntityRank rank = stk::topology::NODE_RANK; rank < stk::topology::NUM_RANKS; ++rank) {
    const BucketVector& bucketsOfRank = m_bulk.buckets(rank);
    const FieldVector& allFieldsOfRank = meta.get_fields(rank);
    const int oldNumBuckets = m_bucketCapacity[rank].size();
    const int newNumBuckets = bucketsOfRank.size();
    const bool rankHasNewFields = not newFields[rank].empty();

    if (rankHasNewFields) {
      resize_field_meta_data_arrays(newFields[rank], 0, newNumBuckets);
    }

    if (newNumBuckets > oldNumBuckets) {
      resize_bucket_arrays(rank, oldNumBuckets, newNumBuckets);
      resize_field_meta_data_arrays(oldFields[rank], oldNumBuckets, newNumBuckets);
    }

    const int clippedNewNumBuckets = std::max(newNumBuckets, oldNumBuckets);
    if (rankHasNewFields || (clippedNewNumBuckets > oldNumBuckets)) {  // Only allow growth here
      resize_bucket_modified_array(rank, allFieldsOfRank.size(), clippedNewNumBuckets);
    }

    const BucketCapacityType oldBucketCapacity = m_bucketCapacity[rank];
    if (not rankHasNewFields) {  // Otherwise, all Buckets will be re-allocated and filled in again below
      reorder_buckets(rank);

    }

    if (newNumBuckets < oldNumBuckets) {  // Allow shrinkage now that things are reorganized
      resize_bucket_arrays(rank, oldNumBuckets, newNumBuckets);
      resize_field_meta_data_arrays(allFieldsOfRank, oldNumBuckets, newNumBuckets);
      resize_bucket_modified_array(rank, allFieldsOfRank.size(), newNumBuckets);
    }

    for (Bucket* bucket : bucketsOfRank) {
      const PartVector& supersetParts = bucket->supersets();
      const unsigned oldBucketId = bucket->ngp_field_bucket_id();
      const unsigned newBucketId = bucket->bucket_id();

      const bool isNewBucket = (oldBucketId == INVALID_BUCKET_ID);
      const bool isGrownBucket = isNewBucket || (bucket->capacity() > oldBucketCapacity[oldBucketId]);
      const bool bucketIsModified = bucket->ngp_field_bucket_is_modified();

      if (isNewBucket || isGrownBucket || rankHasNewFields) {
        allocate_bucket(rank, allFieldsOfRank, supersetParts, newBucketId, bucket->size(), bucket->capacity());
      }
      else {
        update_field_meta_data(allFieldsOfRank, newBucketId, bucket->size());
      }

      auto& hostBucketIsModified = m_hostBucketIsModified[rank];
      for (unsigned fieldIndex = 0; fieldIndex < allFieldsOfRank.size(); ++fieldIndex) {
        hostBucketIsModified(fieldIndex, newBucketId) += bucketIsModified;
      }
      bucket->set_ngp_field_bucket_id(newBucketId);
    }

    Kokkos::deep_copy(m_deviceBucketIsModified[rank], m_hostBucketIsModified[rank]);
  }

  for (int i = 0; i < newNumAllFields; ++i) {
    Kokkos::deep_copy(m_deviceFieldMetaData[i], m_hostFieldMetaData[i]);
  }

  m_synchronizedCount = m_bulk.synchronized_count();
  return true;
}

template <typename Space>
void DeviceFieldDataManager<Space>::update_host_bucket_pointers(Ordinal fieldOrdinal)
{
  const MetaData& meta = m_bulk.mesh_meta_data();
  const FieldBase& fieldBase = *meta.get_fields()[fieldOrdinal];
  const EntityRank rank = fieldBase.entity_rank();

  const BucketVector& bucketsOfRank = m_bulk.buckets(rank);

  HostFieldMetaDataArrayType<mem_space>& hostFieldMetaDataArray = m_hostFieldMetaData[fieldOrdinal];
  for (int bucketId = 0; bucketId < static_cast<int>(bucketsOfRank.size()); ++bucketId) {
    if (fieldBase.has_unified_device_storage()) {
      hostFieldMetaDataArray[bucketId].m_data = get_host_bucket_pointer_for_device(fieldBase, bucketId);
      hostFieldMetaDataArray[bucketId].m_hostData = hostFieldMetaDataArray[bucketId].m_data;
    }
    else {
      hostFieldMetaDataArray[bucketId].m_hostData = get_host_bucket_pointer_for_device(fieldBase, bucketId);
    }
  }
  Kokkos::deep_copy(m_deviceFieldMetaData[fieldOrdinal], hostFieldMetaDataArray);
}

template <typename Space>
void DeviceFieldDataManager<Space>::swap_field_data(Ordinal fieldOrdinal1, Ordinal fieldOrdinal2)
{
  std::swap(m_hostFieldMetaData[fieldOrdinal1], m_hostFieldMetaData[fieldOrdinal2]);
  std::swap(m_deviceFieldMetaData[fieldOrdinal1], m_deviceFieldMetaData[fieldOrdinal2]);
}

template <typename Space>
void DeviceFieldDataManager<Space>::clear_bucket_is_modified(Ordinal fieldOrdinal)
{
  const MetaData& meta = m_bulk.mesh_meta_data();
  const FieldBase& fieldBase = *meta.get_fields()[fieldOrdinal];
  const unsigned fieldRankedOrdinal = fieldBase.field_ranked_ordinal();
  const EntityRank rank = fieldBase.entity_rank();

  auto& hostBucketIsModified = m_hostBucketIsModified[rank];
  for (int bucketId = 0; bucketId < static_cast<int>(hostBucketIsModified.extent(1)); ++bucketId) {
    hostBucketIsModified(fieldRankedOrdinal, bucketId) = 0;
  }
}

template <typename Space>
size_t DeviceFieldDataManager<Space>::get_num_bytes_allocated_on_field(const FieldBase& field) const
{
  if (field.has_unified_device_storage()) {
    return 0;
  }
  else {
    const int fieldOrdinal = field.mesh_meta_data_ordinal();
    const HostFieldMetaDataArrayType<mem_space>& hostFieldMetaDataArray = m_hostFieldMetaData[fieldOrdinal];
    size_t numBytes = 0;
    for (int bucketId = 0; bucketId < static_cast<int>(hostFieldMetaDataArray.extent(0)); ++bucketId) {
      const DeviceFieldMetaData& deviceFieldMetaData = hostFieldMetaDataArray[bucketId];
      if (deviceFieldMetaData.m_data != nullptr) {
        const size_t bytesPerScalar = field.data_traits().alignment_of;
        const size_t bytesPerEntity = bytesPerScalar * deviceFieldMetaData.m_numComponentsPerEntity *
                                      deviceFieldMetaData.m_numCopiesPerEntity;
        const int bytesThisBucket = stk::adjust_up_to_alignment_boundary(bytesPerEntity *
                                                                         deviceFieldMetaData.m_bucketCapacity,
                                                                         DeviceFieldAlignmentSize);
        numBytes += bytesThisBucket;
      }
    }
    return numBytes;
  }
}

template <typename Space>
bool DeviceFieldDataManager<Space>::has_unified_device_storage(Ordinal fieldOrdinal) const
{
  const FieldBase& fieldBase = *m_bulk.mesh_meta_data().get_fields()[fieldOrdinal];
  return fieldBase.has_unified_device_storage();
}

template <typename Space>
void
DeviceFieldDataManager<Space>::set_device_field_meta_data(FieldDataBase& fieldDataBase)
{
  FieldDataBytes<Space>* fieldDataBytes = dynamic_cast<FieldDataBytes<Space>*>(&fieldDataBase);
  STK_ThrowRequireMsg(fieldDataBytes != nullptr,
                      "All device Fields must live in the same memory space.  Found a Field with a different memory "
                      "space than that of the DeviceFieldDataManager");

  const Ordinal fieldOrdinal = fieldDataBytes->field_ordinal();

  fieldDataBytes->m_deviceFieldMetaData = m_deviceFieldMetaData[fieldOrdinal];
}

template <typename Space>
std::any
DeviceFieldDataManager<Space>::get_device_bucket_is_modified(Ordinal fieldOrdinal, int& fieldRankedOrdinal)
{
  const MetaData& meta = m_bulk.mesh_meta_data();
  const FieldBase& fieldBase = *meta.get_fields()[fieldOrdinal];
  const EntityRank rank = fieldBase.entity_rank();

  fieldRankedOrdinal = fieldBase.field_ranked_ordinal();

  return std::any(m_deviceBucketIsModified[rank]);
}

struct DeviceFieldLayoutData
{
  int numBytesPerEntity {};
  int numComponents {};
  int numCopies {};
};

inline DeviceFieldLayoutData get_device_field_layout_data(const FieldBase& field, EntityRank rank,
                                                          const PartVector& parts)
{
  const FieldBase::Restriction& restriction = find_restriction(field, rank, parts);

  DeviceFieldLayoutData layout;

  if (restriction.num_scalars_per_entity() > 0) {
    const int bytesPerScalar = field.data_traits().alignment_of;
    layout.numBytesPerEntity = bytesPerScalar * restriction.num_scalars_per_entity();
    layout.numComponents = restriction.dimension();
    layout.numCopies = restriction.num_scalars_per_entity() / restriction.dimension();
  }

  return layout;
}

template <typename Space>
void DeviceFieldDataManager<Space>::allocate_bucket(EntityRank rank, const FieldVector& fields,
                                                    const PartVector& parts, int bucketId, int size, int capacity)
{
  int totalFieldBytesThisBucket = 0;

  FieldVector separateFields;
  separateFields.reserve(fields.size());

  for (FieldBase* field : fields) {
    const int fieldOrdinal = field->mesh_meta_data_ordinal();
    DeviceFieldMetaData& fieldMetaData = m_hostFieldMetaData[fieldOrdinal][bucketId];
    const DeviceFieldLayoutData layout = get_device_field_layout_data(*field, rank, parts);

    if (layout.numComponents > 0) {
      if (field->has_unified_device_storage()) {
        // Aim this Field at the pre-existing host allocation
        fieldMetaData.m_data = get_host_bucket_pointer_for_device(*field, bucketId);
        fieldMetaData.m_hostData = fieldMetaData.m_data;
      }
      else {
        separateFields.push_back(field);

        const int fieldBytesThisBucket = stk::adjust_up_to_alignment_boundary(static_cast<size_t>(layout.numBytesPerEntity)*
                                                                              capacity, DeviceFieldAlignmentSize);
        totalFieldBytesThisBucket += fieldBytesThisBucket;

        // Temporarily store chunk size in pointer variable; use it to set all pointers later
        fieldMetaData.m_data = reinterpret_cast<std::byte*>(fieldBytesThisBucket);
        fieldMetaData.m_hostData = get_host_bucket_pointer_for_device(*field, bucketId);
      }

      fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
      fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
      fieldMetaData.m_bucketSize = size;
      fieldMetaData.m_bucketCapacity = capacity;
    }
    else {
      fieldMetaData = DeviceFieldMetaData{};  // Field not on this bucket
      fieldMetaData.m_bucketSize = size;
      fieldMetaData.m_bucketCapacity = capacity;
    }
  }

  if (not separateFields.empty()) {
    auto bucketRawData = m_fieldDataAllocator.device_allocate(totalFieldBytesThisBucket);
    fill_field_meta_data_pointers_from_offsets(bucketId, separateFields, bucketRawData, m_hostFieldMetaData);
    m_bucketRawData[rank][bucketId] = bucketRawData;
  }
  else {
    m_bucketRawData[rank][bucketId] = AllocationType();
  }
  m_bucketCapacity[rank][bucketId] = capacity;
}

template <typename Space>
void DeviceFieldDataManager<Space>::resize_field_arrays(int oldNumAllFields, int newNumAllFields)
{
  if (oldNumAllFields == 0) {
    m_deviceFieldMetaData = DeviceFieldMetaDataCollectionType(Kokkos::view_alloc("deviceFieldMetaData_collection",
                                                                                 Kokkos::SequentialHostInit),
                                                              newNumAllFields);
    m_hostFieldMetaData = HostFieldMetaDataCollectionType(Kokkos::view_alloc("hostFieldMetaData_collection",
                                                                             Kokkos::SequentialHostInit),
                                                          newNumAllFields);
  }
  else {
    Kokkos::resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), m_deviceFieldMetaData, newNumAllFields);
    Kokkos::resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), m_hostFieldMetaData, newNumAllFields);
  }

  m_totalNumFields = newNumAllFields;
}

template <typename Space>
void DeviceFieldDataManager<Space>::resize_field_meta_data_arrays(const FieldVector& fields,
                                                                  int oldNumBuckets, int newNumBuckets)
{
  for (const FieldBase* field : fields) {
    const int fieldOrdinal = field->mesh_meta_data_ordinal();
    if (oldNumBuckets == 0) {
      m_deviceFieldMetaData[fieldOrdinal] = DeviceFieldMetaDataArrayType<mem_space>("deviceFieldMetaDataArray_" +
                                                                                    field->name(), newNumBuckets);
      m_hostFieldMetaData[fieldOrdinal] = Kokkos::create_mirror_view(m_deviceFieldMetaData[fieldOrdinal]);
    }
    else {
      Kokkos::resize(m_deviceFieldMetaData[fieldOrdinal], newNumBuckets);
      Kokkos::resize(m_hostFieldMetaData[fieldOrdinal], newNumBuckets);
    }
  }
}

template <typename Space>
void DeviceFieldDataManager<Space>::resize_bucket_arrays(EntityRank rank, int oldNumBuckets, int newNumBuckets)
{
  m_bucketCapacity[rank].resize(newNumBuckets);

  if (oldNumBuckets == 0) {
    m_bucketRawData[rank] = BucketRawDataArrayType("bucketRawDataArray_" + std::to_string(rank), newNumBuckets);
  }
  else {
    Kokkos::resize(m_bucketRawData[rank], newNumBuckets);
  }
}

template <typename Space>
void DeviceFieldDataManager<Space>::resize_bucket_modified_array(EntityRank rank, int newNumFields, int newNumBuckets)
{
  const int oldNumFields = m_hostBucketIsModified[rank].extent(0);
  const int oldNumBuckets = m_hostBucketIsModified[rank].extent(1);

  if (oldNumFields == 0) {
    m_deviceBucketIsModified[rank] = DeviceBucketsModifiedCollectionType("deviceBucketModified_" + std::to_string(rank),
                                                                         newNumFields, newNumBuckets);
    m_hostBucketIsModified[rank] = Kokkos::create_mirror_view(m_deviceBucketIsModified[rank]);
  }
  else {
    auto newDeviceBucketIsModified = DeviceBucketsModifiedCollectionType("deviceBucketModified_" + std::to_string(rank),
                                                                         newNumFields, newNumBuckets);
    auto newHostBucketIsModified = Kokkos::create_mirror_view(newDeviceBucketIsModified);

    const int minNumBuckets = std::min(oldNumBuckets, newNumBuckets);
    auto& oldHostBucketIsModified = m_hostBucketIsModified[rank];
    for (int fieldIdx = 0; fieldIdx < oldNumFields; ++fieldIdx) {
      for (int bucketIdx = 0; bucketIdx < minNumBuckets; ++bucketIdx) {
        newHostBucketIsModified(fieldIdx, bucketIdx) = oldHostBucketIsModified(fieldIdx, bucketIdx);
      }
    }

    m_deviceBucketIsModified[rank] = newDeviceBucketIsModified;
    m_hostBucketIsModified[rank] = newHostBucketIsModified;
  }
}

template <typename Space>
void DeviceFieldDataManager<Space>::reorder_buckets(EntityRank rank)
{
  const MetaData& meta = m_bulk.mesh_meta_data();
  const BucketVector& bucketsOfRank = m_bulk.buckets(rank);
  const FieldVector& fieldsOfRank = meta.get_fields(rank);

  std::vector<BucketShift> bucketShiftList;
  bool anyBucketMovement = false;

  for (const Bucket* bucket : bucketsOfRank) {
    const unsigned oldBucketId = bucket->ngp_field_bucket_id();
    const unsigned newBucketId = bucket->bucket_id();
    const bool isNewBucket = (oldBucketId == INVALID_BUCKET_ID);

    if (not isNewBucket) {  // Otherwise, values will be filled in later
      const bool bucketHasMoved = (oldBucketId != newBucketId);
      anyBucketMovement |= bucketHasMoved;

      bucketShiftList.emplace_back(oldBucketId, newBucketId);
    }
  }

  if (anyBucketMovement) {
    shift_bucket_data(rank, bucketShiftList);
    shift_field_and_bucket_data(rank, fieldsOfRank, bucketShiftList);
  }
}

template <typename Space>
void DeviceFieldDataManager<Space>::shift_bucket_data(EntityRank rank, const std::vector<BucketShift>& bucketShiftList)
{
  const int numBuckets = m_bucketCapacity[rank].size();

  const BucketCapacityType& oldBucketCapacity = m_bucketCapacity[rank];
  BucketCapacityType newBucketCapacity(numBuckets);

  const BucketRawDataArrayType& oldBucketRawData = m_bucketRawData[rank];
  BucketRawDataArrayType newBucketRawData("bucketRawDataArray_" + std::to_string(rank), numBuckets);

  for (const BucketShift& shift : bucketShiftList) {
    newBucketCapacity[shift.newIndex] = oldBucketCapacity[shift.oldIndex];
    newBucketRawData[shift.newIndex] = oldBucketRawData[shift.oldIndex];
  }

  m_bucketCapacity[rank].swap(newBucketCapacity);
  m_bucketRawData[rank] = newBucketRawData;
}

template <typename Space>
void DeviceFieldDataManager<Space>::shift_field_and_bucket_data(EntityRank rank, const FieldVector& fieldsOfRank,
                                                                const std::vector<BucketShift>& bucketShiftList)
{
  const HostBucketsModifiedCollectionType& oldBucketsModified = m_hostBucketIsModified[rank];
  HostBucketsModifiedCollectionType newBucketsModified("hostBucketModified_" + std::to_string(rank),
                                                       oldBucketsModified.extent(0), oldBucketsModified.extent(1));

  for (const FieldBase* field : fieldsOfRank) {
    const int fieldOrdinal = field->mesh_meta_data_ordinal();
    const int fieldRankedOrdinal = field->field_ranked_ordinal();

    const HostFieldMetaDataArrayType<mem_space>& oldHostFieldMetaData = m_hostFieldMetaData[fieldOrdinal];
    HostFieldMetaDataArrayType<mem_space> newHostFieldMetaData("hostFieldMetaData" + std::to_string(fieldOrdinal),
                                                               oldHostFieldMetaData.extent(0));

    for (const BucketShift& shift : bucketShiftList) {
      newHostFieldMetaData[shift.newIndex] = oldHostFieldMetaData[shift.oldIndex];
      newBucketsModified(fieldRankedOrdinal, shift.newIndex) = oldBucketsModified(fieldRankedOrdinal, shift.oldIndex);
    }

    m_hostFieldMetaData[fieldOrdinal] = newHostFieldMetaData;
  }

  m_hostBucketIsModified[rank] = newBucketsModified;
}

template <typename Space>
void DeviceFieldDataManager<Space>::update_field_meta_data(const FieldVector& fields, int bucketId, int bucketSize)
{
  for (const FieldBase* field : fields) {
    const int fieldOrdinal = field->mesh_meta_data_ordinal();
    m_hostFieldMetaData[fieldOrdinal][bucketId].m_bucketSize = bucketSize;

    if (m_hostFieldMetaData[fieldOrdinal][bucketId].m_data != nullptr) {
      if (field->has_unified_device_storage()) {
        m_hostFieldMetaData[fieldOrdinal][bucketId].m_data = get_host_bucket_pointer_for_device(*field, bucketId);
        m_hostFieldMetaData[fieldOrdinal][bucketId].m_hostData = m_hostFieldMetaData[fieldOrdinal][bucketId].m_data;
      }
      else {
        m_hostFieldMetaData[fieldOrdinal][bucketId].m_hostData = get_host_bucket_pointer_for_device(*field, bucketId);
      }
    }
  }
}

template <typename Space>
void DeviceFieldDataManager<Space>::fill_field_meta_data_pointers_from_offsets(
    int bucketId, const FieldVector& fields, const AllocationType& bucketRawData,
    HostFieldMetaDataCollectionType& hostFieldMetaData)
{
  std::byte* bucketPointer = bucketRawData.data();
  uintptr_t pointerOffset = 0;
  for (const FieldBase* field : fields) {
    const int fieldOrdinal = field->mesh_meta_data_ordinal();
    const uintptr_t chunkSize = reinterpret_cast<uintptr_t>(hostFieldMetaData[fieldOrdinal][bucketId].m_data);
    hostFieldMetaData[fieldOrdinal][bucketId].m_data = (chunkSize > 0) ? bucketPointer + pointerOffset : nullptr;
    pointerOffset += chunkSize;
  }
}

}

#endif // DEVICEFIELDDATAMANAGERBASE_HPP
