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

#include "stk_mesh/base/FieldDataManager.hpp"
#include <string.h>                     // for memcpy, memmove, memset
#include <algorithm>                    // for swap
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>  // for FieldMetaData, etc
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/FindRestriction.hpp>
#include "stk_mesh/base/Types.hpp"      // for EntityRank, PartVector
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire, etc
#include "stk_util/util/AdjustForAlignment.hpp"

namespace stk {
namespace mesh {

size_t get_field_bucket_size_in_bytes(const FieldMetaDataArrayType& fieldMetaDataArray, const unsigned bucketId,
                                      const std::byte* endOfField);

struct FieldLayoutData
{
  int numBytesPerEntity {};
  int numComponents {};
  int numCopies {};
};

FieldLayoutData get_field_layout_data(const FieldBase &field, const EntityRank rank, const PartVector& supersetParts)
{
  FieldLayoutData layout;

  const FieldBase::Restriction & restriction = find_and_check_restriction(field, rank, supersetParts);

  if (restriction.num_scalars_per_entity() > 0) {
    const unsigned typeStride = field.data_traits().stride_of;
    layout.numBytesPerEntity = typeStride * restriction.num_scalars_per_entity();
    layout.numComponents = restriction.dimension();
    layout.numCopies = restriction.num_scalars_per_entity() / restriction.dimension();
  }

  return layout;
}

void update_field_pointer(FieldMetaData& fieldMetaData, const size_t capacity, size_t &currentFieldOffset,
                          std::byte* allData, size_t alignmentIncrementBytes)
{
  if (fieldMetaData.m_bytesPerEntity > 0)
  {
    currentFieldOffset = stk::adjust_up_to_alignment_boundary(currentFieldOffset, alignmentIncrementBytes);
    fieldMetaData.m_data = allData + currentFieldOffset;
    currentFieldOffset += fieldMetaData.m_bytesPerEntity * capacity;
  }
}

void initialize_field(FieldMetaData& fieldMetaData, const std::byte* initVal, [[maybe_unused]] unsigned size,
                      [[maybe_unused]] unsigned capacity)
{
#ifdef STK_ASAN_IS_ON
  if (fieldMetaData.m_bytesPerEntity > 0) {
    ASAN_UNPOISON_MEMORY_REGION(fieldMetaData.m_data, size * fieldMetaData.m_bytesPerEntity);
    if (initVal != nullptr) {
      for (unsigned j = 0; j < size; ++j) {
        std::memcpy(fieldMetaData.m_data + j * fieldMetaData.m_bytesPerEntity, initVal,
                    fieldMetaData.m_bytesPerEntity);
      }
    }
    else {
      std::memset(fieldMetaData.m_data, 0, size * fieldMetaData.m_bytesPerEntity);
    }
  }
#else
  if (fieldMetaData.m_bytesPerEntity > 0) {
    if (initVal != nullptr) {
      for (unsigned j = 0; j < capacity; ++j) {
        std::memcpy(fieldMetaData.m_data + j * fieldMetaData.m_bytesPerEntity, initVal,
                    fieldMetaData.m_bytesPerEntity);
      }
    }
    else {
      std::memset(fieldMetaData.m_data, 0, static_cast<size_t>(capacity) * fieldMetaData.m_bytesPerEntity);
    }
  }
#endif
}

void setInitialValue(std::byte* dataLocation, const FieldBase& field, const int numBytesPerEntity)
{
  ASAN_UNPOISON_MEMORY_REGION(dataLocation, numBytesPerEntity);
  const std::byte* initVal = reinterpret_cast<const std::byte*>(field.get_initial_value());
  if (initVal != nullptr) {
    std::memcpy( dataLocation, initVal, numBytesPerEntity );
  }
  else {
    std::memset( dataLocation, 0, numBytesPerEntity );
  }
}

void resize_field_meta_data(FieldBase& field, int newSize)
{
  FieldMetaDataArrayType& fieldMetaDataArray = field.get_internal_field_meta_data();
  if (fieldMetaDataArray.extent(0) == 0u) {
    fieldMetaDataArray = FieldMetaDataArrayType("fieldMetaDataArray_" + field.name(), newSize);
  }
  else {
    Kokkos::resize(fieldMetaDataArray, newSize);
  }
}

FieldDataManager::FieldDataManager(unsigned alignmentIncrementBytes,
                                   std::unique_ptr<AllocatorAdaptorInterface> allocatorAdaptor)
  : m_fieldDataAllocator(std::move(allocatorAdaptor)),
    m_alignmentIncrementBytes(alignmentIncrementBytes)
{
  if (not m_fieldDataAllocator) {
    m_fieldDataAllocator = std::make_unique<AllocatorAdaptor<stk::impl::FieldDataAllocator<std::byte>>>();
  }
}

void
DefaultFieldDataManager::resize_bucket_arrays(const EntityRank rank, const std::vector<FieldBase*>& allFields,
                                              int newNumBuckets)
{
  if (m_fieldRawData.empty()) {
    STK_ThrowRequireMsg(!allFields.empty(),
                        "allocate_bucket_field_data ERROR, field-data-manager was constructed with 0 "
                        "entity-ranks, and there are no fields. Mesh has not been initialized correctly.");

    for (FieldBase* field : allFields) {
      if (field != nullptr) {
        m_fieldRawData.resize(field->get_mesh().mesh_meta_data().entity_rank_count());
        m_bucketCapacity.resize(field->get_mesh().mesh_meta_data().entity_rank_count());
        break;
      }
    }
  }

  if (rank >= static_cast<int>(m_fieldRawData.size())) {
    return;
  }

  m_fieldRawData[rank].resize(newNumBuckets);
  m_bucketCapacity[rank].resize(newNumBuckets);

  for (FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
      if (fieldMetaDataArray.extent(0) == 0u) {
        fieldMetaDataArray = FieldMetaDataArrayType("FieldMetaDataArray_" + field->name(), newNumBuckets);
      }
      else {
        Kokkos::resize(fieldMetaDataArray, newNumBuckets);
      }
    }
  }
}

void
DefaultFieldDataManager::allocate_bucket_ordinal_field_data(const EntityRank rank,
                                                            const std::vector<FieldBase *>& fields,
                                                            const PartVector& supersetParts,
                                                            unsigned bucketOrd,
                                                            unsigned size,
                                                            unsigned capacity)
{
  m_bucketCapacity[rank][bucketOrd] = capacity;

  if (m_numBytesAllocatedPerField.empty()) {
    m_numBytesAllocatedPerField.resize(fields.size(), 0);
  }

  size_t numFields = fields.size();
  // Sizing loop
  size_t totalFieldDataSize = 0;

  for (size_t i = 0; i < numFields; ++i) {
    FieldMetaData newFieldMetaData;

    FieldBase& field = *fields[i];
    if (field.entity_rank() == rank) {
      const FieldLayoutData layout = get_field_layout_data(field, rank, supersetParts);
      if (layout.numBytesPerEntity > 0) {
        newFieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
        newFieldMetaData.m_numComponentsPerEntity = layout.numComponents;
        newFieldMetaData.m_numCopiesPerEntity = layout.numCopies;
        newFieldMetaData.m_bucketSize = size;
        newFieldMetaData.m_bucketCapacity = capacity;
        size_t field_data_size_this_bucket =
            stk::adjust_up_to_alignment_boundary(static_cast<size_t>(layout.numBytesPerEntity)*capacity,
                                                 m_alignmentIncrementBytes);
        totalFieldDataSize += field_data_size_this_bucket;
        m_numBytesAllocatedPerField[i] += field_data_size_this_bucket;
      }
      fields[i]->get_internal_field_meta_data()[bucketOrd] = newFieldMetaData;
      fields[i]->update_cached_field_meta_data();
    }
  }

  // Allocate all field data for this bucket
  if (totalFieldDataSize > 0) {
    std::byte* allData = m_fieldDataAllocator->allocate(totalFieldDataSize);
    m_fieldRawData[rank][bucketOrd] = allData;

    // Set data ptrs in field meta datas
    size_t currentFieldOffset = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
      const FieldBase & field = *fields[i];
      if (field.entity_rank() == rank) {
        const std::byte* initVal = reinterpret_cast<const std::byte*>(field.get_initial_value());
        FieldMetaDataArrayType& fieldMetaDataArray = fields[i]->get_internal_field_meta_data();
        FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(fieldMetaDataArray[bucketOrd]);
        update_field_pointer(fieldMetaData, capacity, currentFieldOffset, allData, m_alignmentIncrementBytes);
        initialize_field(fieldMetaData, initVal, size, capacity);
      }
    }
  }
  else {
    m_fieldRawData[rank][bucketOrd] = nullptr;
  }
}

void
DefaultFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
                                                    const std::vector<FieldBase*>& allFields,
                                                    const PartVector& supersetParts,
                                                    unsigned size,
                                                    unsigned capacity)
{
  const unsigned newNumBuckets = (m_bucketCapacity.empty()) ? 1 : m_bucketCapacity[rank].size() + 1;
  resize_bucket_arrays(rank, allFields, newNumBuckets);

  const unsigned bucketOrd = newNumBuckets - 1;  // New Bucket at the end of the list
  allocate_bucket_ordinal_field_data(rank, allFields, supersetParts, bucketOrd, size, capacity);
}

void DefaultFieldDataManager::allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId,
                                                           const std::vector<FieldBase*>& allFields)
{
  for (FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const unsigned currentSize = field->get_internal_field_meta_data().extent(0);
      if (bucketId >= currentSize) {
        resize_field_meta_data(*field, currentSize+1);
        field->update_cached_field_meta_data();
      }
    }
  }
}

std::vector<BucketFieldSegment>
DefaultFieldDataManager::get_old_bucket_field_offsets(const EntityRank rank,
                                                      const unsigned bucketId,
                                                      const std::vector<FieldBase*>& allFields,
                                                      const unsigned capacity) const
{
  const std::byte* oldAllocationStart = m_fieldRawData[rank][bucketId];

  std::vector<BucketFieldSegment> oldOffsetForField;
  oldOffsetForField.reserve(allFields.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase* field : allFields) {
    const bool isFieldValid = (field != nullptr);
    const bool isTargetRank = isFieldValid && (field->entity_rank() == rank);

    if (isFieldValid && isTargetRank) {
      const FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
      const bool isBucketInRange = (bucketId < fieldMetaDataArray.extent(0));
      const bool hasAllocation = (isBucketInRange) ? (fieldMetaDataArray[bucketId].m_data != nullptr) : false;
      const size_t oldOffsetIntoBucket = (hasAllocation) ? fieldMetaDataArray[bucketId].m_data - oldAllocationStart : 0;
      const size_t oldBytesPerEntity = (hasAllocation) ? fieldMetaDataArray[bucketId].m_bytesPerEntity : 0;
      const size_t oldSizeThisBucket = stk::adjust_up_to_alignment_boundary(oldBytesPerEntity*capacity,
                                                                            m_alignmentIncrementBytes);

      oldOffsetForField.emplace_back(oldOffsetIntoBucket, oldSizeThisBucket);
      totalAllocationSize += oldSizeThisBucket;
    }
  }

  oldOffsetForField.emplace_back(totalAllocationSize, 0);  // Add a one-past-the-end entry for bookkeeping

  return oldOffsetForField;
}

std::vector<BucketFieldSegment>
DefaultFieldDataManager::get_new_bucket_field_offsets(const EntityRank rank,
                                                      const unsigned bucketId,
                                                      const std::vector<FieldBase*>& allFields,
                                                      const unsigned capacity) const
{
  std::vector<BucketFieldSegment> newOffsetForField;
  newOffsetForField.reserve(allFields.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const FieldMetaData & fieldMetaDataForBucket = field->get_internal_field_meta_data()[bucketId];
      size_t newSizeThisBucket =
          stk::adjust_up_to_alignment_boundary(static_cast<size_t>(fieldMetaDataForBucket.m_bytesPerEntity)*capacity,
                                               m_alignmentIncrementBytes);

      newOffsetForField.emplace_back(totalAllocationSize, newSizeThisBucket);
      totalAllocationSize += newSizeThisBucket;
    }
  }

  newOffsetForField.emplace_back(totalAllocationSize, 0);  // Add a one-past-the-end entry for bookkeeping

  return newOffsetForField;
}

void
DefaultFieldDataManager::update_field_meta_data(const EntityRank rank, const unsigned bucketId,
                                                const std::vector<FieldBase*> & allFields,
                                                const PartVector & supersetParts,
                                                unsigned bucketSize,
                                                unsigned bucketCapacity)
{
  for (const FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const FieldLayoutData layout = get_field_layout_data(*field, rank, supersetParts);
      FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(field->get_internal_field_meta_data()[bucketId]);
      fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
      fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
      fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
      fieldMetaData.m_bucketSize = bucketSize;
      fieldMetaData.m_bucketCapacity = bucketCapacity;
    }
  }
}

void
DefaultFieldDataManager::copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                                                unsigned bucketSize,
                                                                unsigned bucketId,
                                                                const std::vector<FieldBase*>& allFields,
                                                                const std::vector<BucketFieldSegment>& oldOffsetForField,
                                                                const std::vector<BucketFieldSegment>& newOffsetForField,
                                                                const std::byte* oldAllocationAllFields,
                                                                std::byte* newAllocationAllFields)
{
  size_t fieldIndex = 0;
  for (const FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const unsigned bytesPerEntity = field->get_internal_field_meta_data()[bucketId].m_bytesPerEntity;
      const bool oldHasAllocation = (oldOffsetForField[fieldIndex].size > 0);
      const unsigned oldNumBytesUsed = (oldHasAllocation) ? bucketSize * bytesPerEntity : 0;

      ASAN_UNPOISON_MEMORY_REGION(newAllocationAllFields + newOffsetForField[fieldIndex].offset, oldNumBytesUsed);
      std::memcpy(newAllocationAllFields + newOffsetForField[fieldIndex].offset,
                  oldAllocationAllFields + oldOffsetForField[fieldIndex].offset, oldNumBytesUsed);
      ++fieldIndex;
    }
  }
}

void
DefaultFieldDataManager::update_field_pointers_to_new_bucket(const EntityRank rank,
                                                             const unsigned bucketId,
                                                             const std::vector<FieldBase*>& allFields,
                                                             const size_t capacity,
                                                             std::byte* newAllocationAllFields)
{
  size_t currentFieldOffset = 0;
  for (size_t i = 0; i < allFields.size(); ++i) {
    const FieldBase& field = *allFields[i];
    if (field.entity_rank() == rank) {
      FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(field.get_internal_field_meta_data()[bucketId]);
      fieldMetaData.m_bucketCapacity = capacity;
      update_field_pointer(fieldMetaData, capacity, currentFieldOffset, newAllocationAllFields,
                           m_alignmentIncrementBytes);
    }
  }
}

void
DefaultFieldDataManager::initialize_new_field_values(FieldBase& newField, const EntityRank rank,
                                                     const unsigned bucketId, unsigned size, unsigned capacity)
{
  if (newField.entity_rank() == rank) {
    const std::byte* initVal = reinterpret_cast<const std::byte*>(newField.get_initial_value());
    FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(newField.get_internal_field_meta_data()[bucketId]);
    initialize_field(fieldMetaData, initVal, size, capacity);
  }
}

void
DefaultFieldDataManager::reallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
                                                      FieldBase & newField,
                                                      const std::vector<FieldBase *> & allFields,
                                                      const PartVector& supersetParts, unsigned bucketSize,
                                                      unsigned bucketCapacity)
{
  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, allFields,
                                                                                   bucketCapacity);
  const int oldBucketAllocationSize = oldOffsetForField.back().offset;

  allocate_new_field_meta_data(rank, bucketId, allFields);
  update_field_meta_data(rank, bucketId, allFields, supersetParts, bucketSize, bucketCapacity);

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, allFields,
                                                                                   bucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;

  if (newBucketAllocationSize > oldBucketAllocationSize) {
    const BucketFieldSegment & lastFieldSegment = newOffsetForField[newOffsetForField.size()-2];
    m_numBytesAllocatedPerField.back() += lastFieldSegment.size;

    std::byte* newAllocationAllFields = m_fieldDataAllocator->allocate(newBucketAllocationSize);
    const std::byte* oldAllocationAllFields = m_fieldRawData[rank][bucketId];

    copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields, oldOffsetForField, newOffsetForField,
                                           oldAllocationAllFields, newAllocationAllFields);

    m_fieldDataAllocator->deallocate(m_fieldRawData[rank][bucketId], oldBucketAllocationSize);
    m_fieldRawData[rank][bucketId] = newAllocationAllFields;
    m_bucketCapacity[rank][bucketId] = bucketCapacity;

    update_field_pointers_to_new_bucket(rank, bucketId, allFields, bucketCapacity, newAllocationAllFields);
    initialize_new_field_values(newField, rank, bucketId, bucketSize, bucketCapacity);
  }
}

void
DefaultFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
                                                      const size_t capacity,
                                                      const std::vector<FieldBase*>&  fields)
{
  if (fields.empty()) {
    return;
  }

  if (m_fieldRawData[rank][bucketId] != nullptr) {
    size_t bytes_to_delete = 0;
    for (unsigned int i = 0; i < fields.size(); ++i) {
      if (fields[i] == nullptr ||
          fields[i]->entity_rank() != rank ||
          fields[i]->get_internal_field_meta_data().extent(0) <= bucketId) {
        continue;
      }

      FieldMetaData& fieldMetaData = fields[i]->get_internal_field_meta_data()[bucketId];
      if (fieldMetaData.m_data != nullptr) {
        const size_t bytes_to_delete_this_field =
            stk::adjust_up_to_alignment_boundary(static_cast<size_t>(fieldMetaData.m_bytesPerEntity)*capacity,
                                                 m_alignmentIncrementBytes);
        m_numBytesAllocatedPerField[i] -= bytes_to_delete_this_field;
        bytes_to_delete += bytes_to_delete_this_field;
        fieldMetaData = FieldMetaData{};
      }
    }

    m_fieldDataAllocator->deallocate(m_fieldRawData[rank][bucketId], bytes_to_delete);
    m_fieldRawData[rank][bucketId] = nullptr;
    m_bucketCapacity[rank][bucketId] = 0;
  }
}

void
DefaultFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                                   const std::vector<unsigned>& reorderedBucketIds)
{
  std::vector<std::byte*> fieldRawData(reorderedBucketIds.size());
  std::vector<unsigned> bucketCapacity(reorderedBucketIds.size());
  for (unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m) {
    fieldRawData[m] = m_fieldRawData[rank][reorderedBucketIds[m]];
    bucketCapacity[m] = m_bucketCapacity[rank][reorderedBucketIds[m]];
  }
  m_fieldRawData[rank].swap(fieldRawData);
  m_bucketCapacity[rank].swap(bucketCapacity);

  for (size_t i = 0; i < fields.size(); ++i) {
    if (fields[i]->entity_rank() == rank) {
      FieldMetaDataArrayType newFieldMetaDataArray("FieldMetaDataArray_" + fields[i]->name(),
                                                   reorderedBucketIds.size());
      FieldMetaDataArrayType& oldFieldMetaDataArray = fields[i]->get_internal_field_meta_data();
      for (unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m) {
        newFieldMetaDataArray[m] = oldFieldMetaDataArray[reorderedBucketIds[m]];
      }
      std::swap(oldFieldMetaDataArray, newFieldMetaDataArray);
      fields[i]->update_cached_field_meta_data();
    }
  }
}

void
DefaultFieldDataManager::allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                             const std::vector<FieldBase*>& allFields)
{
  m_numBytesAllocatedPerField.resize(allFields.size(), 0);
  resize_bucket_arrays(rank, allFields, buckets.size());

  for (Bucket* bucket : buckets) {
    const PartVector& supersetParts = bucket->supersets();
    allocate_bucket_ordinal_field_data(rank, allFields, supersetParts, bucket->bucket_id(), bucket->size(),
                                       bucket->capacity());
  }
}

void
DefaultFieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                               FieldBase & currentField,
                                               const std::vector<FieldBase *> & allFields)
{
  m_numBytesAllocatedPerField.resize(allFields.size(), 0);
  for (size_t i = 0; i < buckets.size(); ++i) {
    const PartVector& supersetParts = buckets[i]->supersets();
    reallocate_bucket_field_data(rank, buckets[i]->bucket_id(), currentField, allFields, supersetParts,
                                 buckets[i]->size(), buckets[i]->capacity());
  }
}

void
DefaultFieldDataManager::remove_field_data_for_entity(EntityRank rank, unsigned bucketId,
                                                      unsigned /*bucketOrd*/, unsigned newBucketSize,
                                                      const std::vector<FieldBase *>& fields)
{
  for (size_t i = 0; i < fields.size(); ++i) {
    const FieldBase& field = *fields[i];
    if (field.entity_rank() == rank) {
      FieldMetaData& fieldMetaData = fields[i]->get_internal_field_meta_data()[bucketId];
      if (fieldMetaData.m_bytesPerEntity > 0) {
        fieldMetaData.m_bucketSize = newBucketSize;
      }
    }
  }
}

void
DefaultFieldDataManager::initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd,
                                                      unsigned newBucketSize, const std::vector<FieldBase *> &fields)
{
  // bucket of bucketId shrinks by one
  for (size_t i = 0; i < fields.size(); ++i) {
    const FieldBase& field = *fields[i];
    if (field.entity_rank() == rank) {
      FieldMetaData& fieldMetaData = fields[i]->get_internal_field_meta_data()[bucketId];
      const int numBytesPerEntity = fieldMetaData.m_bytesPerEntity;

      if (numBytesPerEntity > 0) {
        fieldMetaData.m_bucketSize = newBucketSize;
        setInitialValue(fieldMetaData.m_data + bucketOrd * numBytesPerEntity, field, numBytesPerEntity);
      }
    }
  }
}

void
DefaultFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                                   unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize)
{
  initialize_entity_field_data(dstRank, dstBucketId, dstBucketOrd, newBucketSize, allFields);
}

void
DefaultFieldDataManager::grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                                              unsigned bucketSize, unsigned bucketCapacity)
{
  const int oldBucketCapacity = m_bucketCapacity[rank][bucketId];
  m_bucketCapacity[rank][bucketId] = bucketCapacity;

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, allFields,
                                                                                   bucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;

  if (newBucketAllocationSize == 0) {
    return;
  }

  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, allFields,
                                                                                   oldBucketCapacity);
  const int oldBucketAllocationSize = oldOffsetForField.back().offset;

  unsigned i = 0;
  for (const stk::mesh::FieldBase * field : allFields) {
    if (field->entity_rank() == rank) {
      m_numBytesAllocatedPerField[field->mesh_meta_data_ordinal()] += newOffsetForField[i].size -
                                                                      oldOffsetForField[i].size;
      ++i;
    }
  }

  std::byte* newAllocationAllFields = m_fieldDataAllocator->allocate(newBucketAllocationSize);
  const std::byte* oldAllocationAllFields = m_fieldRawData[rank][bucketId];

  copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields, oldOffsetForField, newOffsetForField,
                                         oldAllocationAllFields, newAllocationAllFields);

  m_fieldDataAllocator->deallocate(m_fieldRawData[rank][bucketId], oldBucketAllocationSize);
  m_fieldRawData[rank][bucketId] = newAllocationAllFields;

  update_field_pointers_to_new_bucket(rank, bucketId, allFields, bucketCapacity, newAllocationAllFields);
}

void
DefaultFieldDataManager::reset_empty_field_data(EntityRank /*rank*/, unsigned bucketId, unsigned bucketSize,
                                                unsigned bucketCapacity, const FieldVector & fields)
{
  for (const FieldBase * field : fields) {
    const FieldMetaData & fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + bucketSize * fieldMetaData.m_bytesPerEntity,
                              (bucketCapacity - bucketSize) * fieldMetaData.m_bytesPerEntity);
  }
}

void reset_field_meta_data_pointers(const size_t bucketIndexBegin, const size_t bucketIndexEnd,
                                    FieldMetaDataArrayType& fieldMetaDataArray, std::byte* oldFieldData,
                                    std::byte* newFieldData)
{
  for (size_t j = bucketIndexBegin; j < bucketIndexEnd; ++j) {
    FieldMetaData &fieldMetaData = fieldMetaDataArray[j];
    if (fieldMetaData.m_bytesPerEntity > 0) {
      size_t sizeOfPreviousBuckets = fieldMetaData.m_data - oldFieldData;
      fieldMetaData.m_data = newFieldData + sizeOfPreviousBuckets;
    }
  }
}

//////////////////////////////////////////////////////////////
ContiguousFieldDataManager::~ContiguousFieldDataManager()
{
  for (size_t i = 0; i < m_fieldRawData.size(); ++i) {
    m_fieldDataAllocator->deallocate(m_fieldRawData[i], m_numBytesAllocatedPerField[i]);
    m_fieldRawData[i] = nullptr;
    m_numBytesAllocatedPerField[i] = 0;
    m_numBytesUsedPerField[i] = 0;
  }
}

void
ContiguousFieldDataManager::initialize_entity_field_data(EntityRank /*rank*/, unsigned /*bucketId*/,
                                                         unsigned /*bucketOrd*/, unsigned /*newBucketSize*/,
                                                         const std::vector<FieldBase *> & /*fields*/)
{
}

void
ContiguousFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
                                                       const std::vector<FieldBase *> & fields,
                                                       const PartVector& supersetParts,
                                                       unsigned size,
                                                       unsigned capacity)
{
  if (m_fieldRawData.empty()) {
    m_fieldRawData.resize(fields.size(), nullptr);
    m_numEntitiesInFieldForBucket.resize(fields.size());
    m_numBytesAllocatedPerField.resize(fields.size(), 0);
    m_numBytesUsedPerField.resize(fields.size(), 0);
  }

  for (size_t i = 0; i < fields.size(); i++) {
    if (fields[i]->entity_rank() == rank) {
      unsigned fieldOrdinal = fields[i]->mesh_meta_data_ordinal();
      const FieldLayoutData layout = get_field_layout_data(*fields[i], rank, supersetParts);
      FieldMetaData fieldMetaData;
      if (layout.numBytesPerEntity > 0) {
        fieldMetaData.m_data = m_fieldRawData[fieldOrdinal] + m_numBytesUsedPerField[fieldOrdinal];
        fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
        fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
        fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
        fieldMetaData.m_bucketSize = size;
        fieldMetaData.m_bucketCapacity = capacity;
      }
      const unsigned oldSize = fields[i]->get_internal_field_meta_data().extent(0);
      resize_field_meta_data(*fields[i], oldSize+1);
      fields[i]->update_cached_field_meta_data();
      fields[i]->get_internal_field_meta_data()[oldSize] = fieldMetaData;
      m_numEntitiesInFieldForBucket[fieldOrdinal].push_back(0);
    }
  }
}

void
ContiguousFieldDataManager::clear_bucket_field_data(const EntityRank rmRank, const unsigned rmBucketId,
                                                    const std::vector<FieldBase*>& allFields)
{
  for (size_t fieldIndex = 0; fieldIndex < allFields.size(); ++fieldIndex) {
    const FieldBase& field = *allFields[fieldIndex];
    unsigned fieldOrdinal = field.mesh_meta_data_ordinal();

    if (field.entity_rank() == rmRank && m_numEntitiesInFieldForBucket[fieldOrdinal][rmBucketId] > 0) {
      int numBytesPerEntity = field.get_internal_field_meta_data()[rmBucketId].m_bytesPerEntity;
      const std::byte* endOfField = m_fieldRawData[fieldOrdinal] + m_numBytesUsedPerField[fieldOrdinal];
      size_t sizeOfBucketToRemove = get_field_bucket_size_in_bytes(field.get_internal_field_meta_data(), rmBucketId,
                                                                   endOfField);

      if (numBytesPerEntity > 0) {
        std::byte* newFieldData = m_fieldRawData[fieldOrdinal];

        FieldMetaDataArrayType& fieldMetaDataArray =
            const_cast<FieldMetaDataArrayType&>(field.get_internal_field_meta_data());

        FieldMetaData &fieldMetaDataForModifiedBucket = fieldMetaDataArray[rmBucketId];
        size_t sizeOfBucketsToTheLeft = fieldMetaDataForModifiedBucket.m_data - m_fieldRawData[fieldOrdinal];
        size_t rightHalfSize = m_numBytesUsedPerField[fieldOrdinal] - sizeOfBucketsToTheLeft - sizeOfBucketToRemove;

        fieldMetaDataForModifiedBucket.m_data = nullptr;

        size_t numBucketsOfRank = fieldMetaDataArray.size();
        reset_field_meta_data_pointers(rmBucketId+1, numBucketsOfRank, fieldMetaDataArray,
                                       m_fieldRawData[fieldOrdinal], newFieldData-sizeOfBucketToRemove);

        ASAN_UNPOISON_MEMORY_REGION(newFieldData + sizeOfBucketsToTheLeft, sizeOfBucketToRemove + rightHalfSize);
        std::memmove(newFieldData + sizeOfBucketsToTheLeft,
                     m_fieldRawData[fieldOrdinal] + sizeOfBucketsToTheLeft + sizeOfBucketToRemove, rightHalfSize);

        m_numBytesUsedPerField[fieldOrdinal] -= sizeOfBucketToRemove;
        m_fieldRawData[fieldOrdinal] = newFieldData;
        m_numEntitiesInFieldForBucket[fieldOrdinal][rmBucketId] = 0;
      }
    }
  }
}


void
ContiguousFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
                                                         const size_t /*capacity*/,
                                                         const std::vector<FieldBase*>&  fields)
{
  if (fields.empty()) {
    return;
  }

  this->clear_bucket_field_data(rank, bucketId, fields);

  for (size_t fieldIndex = 0; fieldIndex < fields.size(); ++fieldIndex) {
    if (fields[fieldIndex]->entity_rank() == rank) {
      unsigned fieldOrdinal = fields[fieldIndex]->mesh_meta_data_ordinal();

      FieldMetaData& fieldMetaData = fields[fieldIndex]->get_internal_field_meta_data()[bucketId];
      fieldMetaData = FieldMetaData{};

      STK_ThrowRequireMsg(m_numEntitiesInFieldForBucket[fieldOrdinal][bucketId] == 0, "Bucket not empty!");
    }
  }
}

void
ContiguousFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                                      const std::vector<unsigned> & reorderedBucketIds)
{
  for (size_t fieldIndex = 0; fieldIndex < fields.size(); ++fieldIndex) {
    if (fields[fieldIndex]->entity_rank() == rank) {
      FieldMetaDataArrayType &oldMetaData =
          const_cast<FieldMetaDataArrayType&>(fields[fieldIndex]->get_internal_field_meta_data());
      unsigned fieldOrdinal = fields[fieldIndex]->mesh_meta_data_ordinal();
      const size_t newFieldSize = m_numBytesUsedPerField[fieldOrdinal] + m_extraCapacity;
      std::byte* newFieldData = m_fieldDataAllocator->allocate(newFieldSize);

      FieldMetaDataArrayType newMetaData("FieldMetaDataArray_" + fields[fieldIndex]->name(), reorderedBucketIds.size());
      std::vector<size_t> newNumEntitiesPerBucket(reorderedBucketIds.size(), 0);
      unsigned newOffset = 0;
      for (unsigned bucketIndex = 0, bucketEnd = reorderedBucketIds.size(); bucketIndex < bucketEnd; ++bucketIndex) {
        unsigned oldBucketIndex = reorderedBucketIds[bucketIndex];
        const std::byte* bucketStartPtr = oldMetaData[oldBucketIndex].m_data;

        if (oldMetaData[oldBucketIndex].m_bytesPerEntity > 0) {
          newNumEntitiesPerBucket[bucketIndex] = m_numEntitiesInFieldForBucket[fieldOrdinal][oldBucketIndex];

          const std::byte* endOfField = m_fieldRawData[fieldOrdinal]+m_numBytesUsedPerField[fieldOrdinal];
          unsigned bucketSize = get_field_bucket_size_in_bytes(oldMetaData, oldBucketIndex, endOfField);
          ASAN_UNPOISON_MEMORY_REGION(bucketStartPtr, bucketSize);
          ASAN_UNPOISON_MEMORY_REGION(newFieldData+newOffset, bucketSize);
          newMetaData[bucketIndex] = oldMetaData[oldBucketIndex];
          newMetaData[bucketIndex].m_data = &newFieldData[newOffset];

          std::memcpy(newFieldData+newOffset, bucketStartPtr, bucketSize);
          newOffset += bucketSize;
        }
      }

      m_fieldDataAllocator->deallocate(m_fieldRawData[fieldOrdinal], m_numBytesAllocatedPerField[fieldOrdinal]);
      STK_ThrowRequire(newOffset == m_numBytesUsedPerField[fieldOrdinal]);
      m_numBytesAllocatedPerField[fieldOrdinal] = newFieldSize;
      m_fieldRawData[fieldOrdinal] = newFieldData;
      m_numEntitiesInFieldForBucket[fieldOrdinal].swap(newNumEntitiesPerBucket);
      std::swap(oldMetaData, newMetaData);
      fields[fieldIndex]->update_cached_field_meta_data();
    }
  }
}

void
ContiguousFieldDataManager::remove_field_data_for_entity(EntityRank rmRank, unsigned rmBucketId,
                                                         unsigned /*rmBucketOrd*/, unsigned newBucketSize,
                                                         const std::vector<FieldBase *> &allFields)
{
  for(size_t fieldIndex = 0; fieldIndex < allFields.size(); ++fieldIndex) {
    const FieldBase& field = *allFields[fieldIndex];
    unsigned fieldOrdinal = field.mesh_meta_data_ordinal();
    if (field.entity_rank() == rmRank) {
      int numBytesPerEntity = field.get_internal_field_meta_data()[rmBucketId].m_bytesPerEntity;
      if (numBytesPerEntity > 0) {
        std::byte* newFieldData = m_fieldRawData[fieldOrdinal];
        const std::byte* endOfField = m_fieldRawData[fieldOrdinal]+m_numBytesUsedPerField[fieldOrdinal];

        const size_t currentBucketStorageUsed = numBytesPerEntity *
            m_numEntitiesInFieldForBucket[fieldOrdinal][rmBucketId];
        const size_t currentBucketAllocation = get_field_bucket_size_in_bytes(field.get_internal_field_meta_data(),
                                                                              rmBucketId, endOfField);
        const size_t newBucketStorageUsed = currentBucketStorageUsed - numBytesPerEntity;
        const size_t newBucketAllocation = stk::adjust_up_to_alignment_boundary(newBucketStorageUsed,
                                                                                m_alignmentIncrementBytes);
        const size_t allocationToRemove = currentBucketAllocation - newBucketAllocation;

        FieldMetaDataArrayType& fieldMetaDataArray =
            const_cast<FieldMetaDataArrayType&>(field.get_internal_field_meta_data());
        FieldMetaData &fieldMetaDataForModifiedBucket = fieldMetaDataArray[rmBucketId];
        size_t sizeOfBucketsToTheLeft = fieldMetaDataForModifiedBucket.m_data - m_fieldRawData[fieldOrdinal];
        size_t rightHalfSize = m_numBytesUsedPerField[fieldOrdinal] - sizeOfBucketsToTheLeft -
            currentBucketAllocation;

        fieldMetaDataForModifiedBucket.m_data = newFieldData + sizeOfBucketsToTheLeft;
        fieldMetaDataForModifiedBucket.m_bucketSize = newBucketSize;

        size_t numBucketsOfRank = fieldMetaDataArray.size();
        reset_field_meta_data_pointers(rmBucketId+1, numBucketsOfRank, fieldMetaDataArray,
                                       m_fieldRawData[fieldOrdinal], newFieldData-allocationToRemove);

        ASAN_UNPOISON_MEMORY_REGION(newFieldData + sizeOfBucketsToTheLeft + currentBucketStorageUsed,
                                    currentBucketAllocation - currentBucketStorageUsed);
        ASAN_UNPOISON_MEMORY_REGION(newFieldData + sizeOfBucketsToTheLeft + currentBucketAllocation,
                                    rightHalfSize);
        std::memmove(newFieldData + sizeOfBucketsToTheLeft + newBucketAllocation,
                     m_fieldRawData[fieldOrdinal] + sizeOfBucketsToTheLeft + currentBucketAllocation, rightHalfSize);
        m_numBytesUsedPerField[fieldOrdinal] -= allocationToRemove;
        m_fieldRawData[fieldOrdinal] = newFieldData;
        m_numEntitiesInFieldForBucket[fieldOrdinal][rmBucketId] -= 1;
      }
    }
  }
}

void
ContiguousFieldDataManager::allocate_field_data(EntityRank rank,
                                                const std::vector<Bucket*>& buckets,
                                                const std::vector<FieldBase *> & fields)
{
  m_fieldRawData.resize(fields.size(), nullptr);
  m_numEntitiesInFieldForBucket.resize(fields.size());
  for (size_t i=0; i<fields.size(); i++) {
    if (fields[i]->entity_rank() == rank) {
      unsigned fieldOrdinal = fields[i]->mesh_meta_data_ordinal();
      m_numEntitiesInFieldForBucket[fieldOrdinal].resize(buckets.size(),0);
    }
  }

  m_numBytesAllocatedPerField.resize(fields.size(), m_extraCapacity);
  m_numBytesUsedPerField.resize(fields.size(), 0);
  for (size_t fieldIndex = 0; fieldIndex != fields.size(); fieldIndex++) {
    FieldBase& field = *fields[fieldIndex];
    if (field.entity_rank() == rank) {
      unsigned fieldOrdinal = fields[fieldIndex]->mesh_meta_data_ordinal();
      for (size_t i = 0; i < buckets.size(); ++i) {
        const Bucket& bucket = *buckets[i];
        const PartVector& supersetParts = bucket.supersets();
        const FieldLayoutData layout = get_field_layout_data(field, rank, supersetParts);
        if (layout.numBytesPerEntity > 0) {
          m_numEntitiesInFieldForBucket[fieldOrdinal][i] = bucket.size();
          m_numBytesUsedPerField[fieldOrdinal] +=
              stk::adjust_up_to_alignment_boundary(layout.numBytesPerEntity * bucket.size(),
                                                   m_alignmentIncrementBytes);
        }
        FieldMetaData fieldMetaData;
        fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
        fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
        fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
        fieldMetaData.m_bucketSize = bucket.size();
        fieldMetaData.m_bucketCapacity = bucket.capacity();
        const unsigned oldSize = field.get_internal_field_meta_data().extent(0);
        resize_field_meta_data(field, oldSize+1);
        field.update_cached_field_meta_data();
        field.get_internal_field_meta_data()[oldSize] = fieldMetaData;
      }

      m_numBytesAllocatedPerField[fieldOrdinal] += m_numBytesUsedPerField[fieldOrdinal];
      m_fieldRawData[fieldOrdinal] = m_fieldDataAllocator->allocate(m_numBytesAllocatedPerField[fieldOrdinal]);

      size_t offset = 0;
      for (size_t i = 0; i < buckets.size(); ++i) {
        const std::byte* initVal = reinterpret_cast<const std::byte*>(field.get_initial_value());
        FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(field.get_internal_field_meta_data()[i]);
        update_field_pointer(fieldMetaData, buckets[i]->size(), offset, m_fieldRawData[fieldOrdinal],
                             m_alignmentIncrementBytes);
        initialize_field(fieldMetaData, initVal, buckets[i]->size(), buckets[i]->size());
      }
    }
  }
}

void
ContiguousFieldDataManager::allocate_new_field_meta_data(const EntityRank /*rank*/,
                                                         const std::vector<Bucket*> & buckets,
                                                         const std::vector<FieldBase*>& allFields)
{
  for (FieldBase* field : allFields) {
    for (stk::mesh::Bucket* bucket : buckets) {
      const unsigned oldSize = field->get_internal_field_meta_data().extent(0);
      if (bucket->bucket_id() >= oldSize) {
        resize_field_meta_data(*field, oldSize+1);
      }
    }
  }
}

std::vector<size_t>
ContiguousFieldDataManager::get_field_bucket_offsets(const std::vector<Bucket*> & buckets,
                                                     FieldBase & currentField) const
{
  std::vector<size_t> offsetForBucket;
  size_t currentBucketOffset = 0;
  for (unsigned int i = 0; i < buckets.size(); ++i) {
    const size_t bytesPerEntity = currentField.get_internal_field_meta_data()[i].m_bytesPerEntity;
    size_t bucketDataSizeThisField = stk::adjust_up_to_alignment_boundary(bytesPerEntity*buckets[i]->size(),
                                                                          m_alignmentIncrementBytes);

    offsetForBucket.push_back(currentBucketOffset);
    currentBucketOffset += bucketDataSizeThisField;
  }
  offsetForBucket.push_back(currentBucketOffset);
  return offsetForBucket;
}

void
ContiguousFieldDataManager::copy_bucket_data_from_old_to_new_field(const std::vector<size_t>& oldOffsetForBucket,
                                                                   const std::vector<size_t>& newOffsetForBucket,
                                                                   const std::byte* oldAllocationAllBuckets,
                                                                   std::byte* newAllocationAllBuckets)
{
  for (size_t bucketIndex = 0; bucketIndex < oldOffsetForBucket.size()-1; ++bucketIndex) {
    size_t oldBucketAllocationSize = oldOffsetForBucket[bucketIndex + 1] - oldOffsetForBucket[bucketIndex];
    ASAN_UNPOISON_MEMORY_REGION(oldAllocationAllBuckets + oldOffsetForBucket[bucketIndex], oldBucketAllocationSize);
    ASAN_UNPOISON_MEMORY_REGION(newAllocationAllBuckets + newOffsetForBucket[bucketIndex], oldBucketAllocationSize);
    std::memcpy(newAllocationAllBuckets + newOffsetForBucket[bucketIndex],
                oldAllocationAllBuckets + oldOffsetForBucket[bucketIndex], oldBucketAllocationSize);
  }
}

void
ContiguousFieldDataManager::update_bucket_storage_for_field(EntityRank rank,
                                                            const std::vector<Bucket*>& buckets,
                                                            FieldBase& currentField)
{
  const unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
  for (size_t i = 0; i < buckets.size(); ++i) {
    const Bucket& bucket = *buckets[i];
    const PartVector& supersetParts = bucket.supersets();
    const FieldLayoutData layout = get_field_layout_data(currentField, rank, supersetParts);
    if (layout.numBytesPerEntity > 0) {
      m_numEntitiesInFieldForBucket[fieldOrdinal][i] = bucket.size();
      m_numBytesUsedPerField[fieldOrdinal] +=
          stk::adjust_up_to_alignment_boundary(layout.numBytesPerEntity * bucket.size(),
                                               m_alignmentIncrementBytes);
    }
    FieldMetaData& fieldMetaData = currentField.get_internal_field_meta_data()[i];
    fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
    fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
    fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
    fieldMetaData.m_bucketSize = bucket.size();
    fieldMetaData.m_bucketCapacity = bucket.capacity();
  }
}

void
ContiguousFieldDataManager::update_bucket_pointers_to_new_field(const std::vector<Bucket*>& buckets,
                                                                FieldBase& currentField)
{
  const unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
  size_t offset = 0;
  for (size_t i = 0; i < buckets.size(); ++i) {
    FieldMetaData& fieldMetaData = currentField.get_internal_field_meta_data()[i];
    update_field_pointer(fieldMetaData, buckets[i]->size(), offset, m_fieldRawData[fieldOrdinal],
                         m_alignmentIncrementBytes);
  }
}

void
ContiguousFieldDataManager::initialize_new_bucket_values(const std::vector<Bucket*>& buckets,
                                                         const std::vector<size_t> & oldOffsetForBucket,
                                                         const std::vector<size_t> & newOffsetForBucket,
                                                         FieldBase& currentField)
{
  for (size_t i = 0; i < buckets.size(); ++i) {
    const size_t oldSize = oldOffsetForBucket[i+1] - oldOffsetForBucket[i];
    const size_t newSize = newOffsetForBucket[i+1] - newOffsetForBucket[i];
    const bool isNewBucketForField = (oldSize == 0u) && (newSize > 0u);
    if (isNewBucketForField) {
      const std::byte* initVal = reinterpret_cast<const std::byte*>(currentField.get_initial_value());
      FieldMetaData& fieldMetaData = currentField.get_internal_field_meta_data()[i];
      initialize_field(fieldMetaData, initVal, buckets[i]->size(), buckets[i]->size());
    }
  }
}

void
ContiguousFieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                                  FieldBase& currentField,
                                                  const std::vector<FieldBase*>& allFields)
{
  m_fieldRawData.resize(allFields.size(), nullptr);
  m_numBytesAllocatedPerField.resize(allFields.size(), 0);
  m_numBytesUsedPerField.resize(allFields.size(), 0);

  m_numEntitiesInFieldForBucket.resize(allFields.size());
  for (unsigned fieldOrdinal = 0; fieldOrdinal < allFields.size(); ++fieldOrdinal) {
    if (allFields[fieldOrdinal]->entity_rank() == rank) {
      auto& bucketListForField = m_numEntitiesInFieldForBucket[fieldOrdinal];
      bucketListForField.resize(buckets.size(), 0);
    }
  }

  if (currentField.entity_rank() == rank) {
    allocate_new_field_meta_data(rank, buckets, allFields);
    currentField.update_cached_field_meta_data();

    unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
    std::vector<size_t> oldOffsetForBucket = get_field_bucket_offsets(buckets, currentField);
    const size_t oldFieldBucketAllocationSize = m_numBytesAllocatedPerField[fieldOrdinal];

    update_bucket_storage_for_field(rank, buckets, currentField);

    std::vector<size_t> newOffsetForBucket = get_field_bucket_offsets(buckets, currentField);
    const size_t newFieldBucketAllocationSize = newOffsetForBucket.back() + m_extraCapacity;

    if (newFieldBucketAllocationSize > oldFieldBucketAllocationSize) {
      std::byte* newAllocationAllBuckets = m_fieldDataAllocator->allocate(newFieldBucketAllocationSize);
      const std::byte* oldAllocationAllBuckets = m_fieldRawData[fieldOrdinal];

      copy_bucket_data_from_old_to_new_field(oldOffsetForBucket, newOffsetForBucket, oldAllocationAllBuckets,
                                             newAllocationAllBuckets);

      m_fieldDataAllocator->deallocate(m_fieldRawData[fieldOrdinal], oldFieldBucketAllocationSize);
      m_fieldRawData[fieldOrdinal] = newAllocationAllBuckets;
      m_numBytesAllocatedPerField[fieldOrdinal] = newFieldBucketAllocationSize;

      update_bucket_pointers_to_new_field(buckets, currentField);
      initialize_new_bucket_values(buckets, oldOffsetForBucket, newOffsetForBucket, currentField);
    }
  }
}

size_t get_field_bucket_size_in_bytes(const FieldMetaDataArrayType& fieldMetaDataArray, const unsigned bucketId,
                                      const std::byte* endOfField)
{
  size_t sizeFieldThisBucketInBytes = endOfField - fieldMetaDataArray[bucketId].m_data;

  for (unsigned nextBucket = bucketId+1; nextBucket < fieldMetaDataArray.size(); ++nextBucket) {
    if ( fieldMetaDataArray[nextBucket].m_bytesPerEntity > 0 ) {
      sizeFieldThisBucketInBytes = fieldMetaDataArray[nextBucket].m_data - fieldMetaDataArray[bucketId].m_data;
      break;
    }
  }

  return sizeFieldThisBucketInBytes;
}

void
ContiguousFieldDataManager::grow_bucket_capacity(const FieldVector& allFields, EntityRank rank,
                                                 unsigned bucketId, unsigned /*bucketSize*/, unsigned bucketCapacity)
{
  for (size_t fieldIndex = 0; fieldIndex < allFields.size(); ++fieldIndex) {
    const FieldBase& field = *allFields[fieldIndex];
    if (field.entity_rank() == rank) {
      FieldMetaDataArrayType& fieldMetaDataArray =
          const_cast<FieldMetaDataArrayType&>(field.get_internal_field_meta_data());
      FieldMetaData &fieldMetaData = fieldMetaDataArray[bucketId];
      fieldMetaData.m_bucketCapacity = bucketCapacity;
    }
  }
}

void
ContiguousFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase*>& allFields,
                                                      EntityRank dstRank, unsigned dstBucketId,
                                                      unsigned dstBucketOrd, unsigned newBucketSize)
{
  for (size_t fieldIndex = 0; fieldIndex < allFields.size(); ++fieldIndex) {
    const FieldBase& field = *allFields[fieldIndex];
    unsigned fieldOrdinal = field.mesh_meta_data_ordinal();
    if (field.entity_rank() == dstRank) {
      int numBytesPerEntity = field.get_internal_field_meta_data()[dstBucketId].m_bytesPerEntity;
      if (numBytesPerEntity > 0) {
        const std::byte* endOfField = m_fieldRawData[fieldOrdinal]+m_numBytesUsedPerField[fieldOrdinal];
        const size_t currentBucketStorageUsed = numBytesPerEntity *
            m_numEntitiesInFieldForBucket[fieldOrdinal][dstBucketId];
        const size_t currentBucketAllocation = get_field_bucket_size_in_bytes(field.get_internal_field_meta_data(),
                                                                              dstBucketId, endOfField);
        const size_t newBucketStorageUsed = currentBucketStorageUsed + numBytesPerEntity;
        const size_t newBucketAllocation = stk::adjust_up_to_alignment_boundary(newBucketStorageUsed,
                                                                                m_alignmentIncrementBytes);
        const size_t extraAllocationNeeded = newBucketAllocation - currentBucketAllocation;
        const size_t newFieldSizeNeeded = m_numBytesUsedPerField[fieldOrdinal] + extraAllocationNeeded;

        bool requiresNewAllocation = false;
        size_t newFieldSize = m_numBytesAllocatedPerField[fieldOrdinal];

        // Only reallocate if we've outgrown the extra capacity
        if (newFieldSizeNeeded > m_numBytesAllocatedPerField[fieldOrdinal]) {
          requiresNewAllocation = true;
          newFieldSize = newFieldSizeNeeded + m_extraCapacity;
        }

        std::byte* newFieldData = m_fieldRawData[fieldOrdinal];
        FieldMetaDataArrayType& fieldMetaDataArray =
            const_cast<FieldMetaDataArrayType&>(field.get_internal_field_meta_data());

        if (requiresNewAllocation) {
          newFieldData = m_fieldDataAllocator->allocate(newFieldSize);
          reset_field_meta_data_pointers(0, dstBucketId, fieldMetaDataArray, m_fieldRawData[fieldOrdinal],
                                         newFieldData);
        }

        FieldMetaData &fieldMetaDataForModifiedBucket = fieldMetaDataArray[dstBucketId];
        size_t sizeOfBucketsToTheLeft = fieldMetaDataForModifiedBucket.m_data - m_fieldRawData[fieldOrdinal];
        size_t leftHalfSize = sizeOfBucketsToTheLeft + currentBucketStorageUsed;
        size_t rightHalfSize = m_numBytesUsedPerField[fieldOrdinal] - sizeOfBucketsToTheLeft -
            currentBucketAllocation;
        fieldMetaDataForModifiedBucket.m_data = newFieldData + sizeOfBucketsToTheLeft;
        fieldMetaDataForModifiedBucket.m_bucketSize = newBucketSize;

        size_t numBucketsOfRank = fieldMetaDataArray.extent(0);
        reset_field_meta_data_pointers(dstBucketId + 1, numBucketsOfRank, fieldMetaDataArray,
                                       m_fieldRawData[fieldOrdinal], newFieldData + extraAllocationNeeded);

        if (requiresNewAllocation) {
          // We lose some safety here because, ideally, we would only unpoison the used portions of the
          // allocation and skip over the SIMD padding, but that requires querying each Bucket's size.
          // During mesh construction, querying the Buckets would be done too early which would put the
          // BucketRegistrar in a bad state.
          ASAN_UNPOISON_MEMORY_REGION(m_fieldRawData[fieldOrdinal],
                                      leftHalfSize + currentBucketAllocation + rightHalfSize);
          ASAN_UNPOISON_MEMORY_REGION(newFieldData,
                                      leftHalfSize + newBucketAllocation + rightHalfSize);
          std::memcpy(newFieldData, m_fieldRawData[fieldOrdinal], leftHalfSize);
          std::memcpy(newFieldData + sizeOfBucketsToTheLeft + newBucketAllocation,
                      m_fieldRawData[fieldOrdinal] + sizeOfBucketsToTheLeft + currentBucketAllocation,
                      rightHalfSize);
          m_fieldDataAllocator->deallocate(m_fieldRawData[fieldOrdinal],
                                           m_numBytesAllocatedPerField[fieldOrdinal]);
          m_numBytesAllocatedPerField[fieldOrdinal] = newFieldSize;
        }
        else {
          ASAN_UNPOISON_MEMORY_REGION(m_fieldRawData[fieldOrdinal] + sizeOfBucketsToTheLeft,
                                      newBucketAllocation + rightHalfSize);
          std::memmove(newFieldData + sizeOfBucketsToTheLeft + newBucketAllocation,
                       m_fieldRawData[fieldOrdinal] + sizeOfBucketsToTheLeft + currentBucketAllocation,
                       rightHalfSize);
        }

        setInitialValue(newFieldData+leftHalfSize, field, numBytesPerEntity);

        m_numBytesUsedPerField[fieldOrdinal] += extraAllocationNeeded;
        m_fieldRawData[fieldOrdinal] = newFieldData;
        STK_ThrowRequire(dstBucketOrd == m_numEntitiesInFieldForBucket[fieldOrdinal][dstBucketId]);
        m_numEntitiesInFieldForBucket[fieldOrdinal][dstBucketId] += 1;
      }
    }
  }
}

void
ContiguousFieldDataManager::swap_fields(const int field1, const int field2)
{
  std::swap(m_fieldRawData[field1], m_fieldRawData[field2]);
  STK_ThrowRequire(m_numBytesAllocatedPerField[field1] == m_numBytesAllocatedPerField[field2]);
  STK_ThrowRequire(m_numBytesUsedPerField[field1] == m_numBytesUsedPerField[field2]);
  STK_ThrowRequire(m_numEntitiesInFieldForBucket[field1].size() == m_numEntitiesInFieldForBucket[field2].size());
}

void
ContiguousFieldDataManager::reset_empty_field_data(EntityRank /*rank*/, unsigned bucketId, unsigned bucketSize,
                                                   unsigned /*bucketCapacity*/, const FieldVector & fields)
{
  for (const FieldBase* field : fields) {
    const FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    const unsigned bucketStorageUsed = bucketSize * fieldMetaData.m_bytesPerEntity;
    const unsigned bucketStorageAllocated = stk::adjust_up_to_alignment_boundary(bucketStorageUsed,
                                                                                 m_alignmentIncrementBytes);
    ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + bucketStorageUsed, bucketStorageAllocated - bucketStorageUsed);
  }
}

}
}
