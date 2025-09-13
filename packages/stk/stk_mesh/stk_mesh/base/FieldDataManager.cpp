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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>  // for FieldMetaData, etc
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/FindRestriction.hpp>
#include "stk_mesh/base/Types.hpp"      // for EntityRank, PartVector
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire, etc
#include "stk_util/util/AdjustForAlignment.hpp"
#include <algorithm>                    // for swap

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

void update_field_pointers_to_new_bucket(const EntityRank rank,
                                         const unsigned bucketId,
                                         const std::vector<FieldBase*>& allFields,
                                         const size_t capacity,
                                         std::byte* newAllocationAllFields,
                                         unsigned alignmentIncrementBytes)
{
  size_t currentFieldOffset = 0;
  for (size_t i = 0; i < allFields.size(); ++i) {
    const FieldBase& field = *allFields[i];
    if (field.entity_rank() == rank) {
      FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(field.get_internal_field_meta_data()[bucketId]);
      fieldMetaData.m_bucketCapacity = capacity;
      update_field_pointer(fieldMetaData, capacity, currentFieldOffset, newAllocationAllFields,
                           alignmentIncrementBytes);
    }
  }
}

void initialize_field_on_bucket(const stk::mesh::FieldBase& field, int bucketId, const FieldMetaData& fieldMetaData,
                                unsigned size, unsigned capacity)
{
  // Poison it all, and then unpoison each byte before writing to it below
  ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data, capacity * fieldMetaData.m_bytesPerEntity);

  const std::byte* initVal = reinterpret_cast<const std::byte*>(field.get_initial_value());

  stk::mesh::bucket_bytes_execute<std::byte>(field, bucketId,
    [&](auto& bucketBytes) {
      if (bucketBytes.is_field_defined()) {
        if (initVal == nullptr) {
          std::vector<std::byte> zeroInit(size*bucketBytes.num_bytes());

          for (stk::mesh::EntityIdx entity : bucketBytes.entities()) {
            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              ASAN_UNPOISON_MEMORY_REGION(&bucketBytes(entity, byte), 1);
              bucketBytes(entity, byte) = zeroInit[byte];
            }
          }
        }
        else {
          for (stk::mesh::EntityIdx entity : bucketBytes.entities()) {
            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              ASAN_UNPOISON_MEMORY_REGION(&bucketBytes(entity, byte), 1);
              bucketBytes(entity, byte) = initVal[byte];
            }
          }
        }
      }
    }
  );
}

template <typename EntityBytesType>
void copy_init_into_entity(EntityBytesType& entityBytes, const std::byte* initVal)
{
  if (entityBytes.is_field_defined()) {
    if (initVal == nullptr) {
      std::vector<std::byte> zeroInit(entityBytes.num_bytes());

      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
        entityBytes(byte) = zeroInit[byte];
      }
    }
    else {
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
        entityBytes(byte) = initVal[byte];
      }
    }
  }
}

void initialize_field_on_entity(const stk::mesh::FieldBase& field, unsigned bucketId, unsigned bucketOrd)
{
  const std::byte* initVal = reinterpret_cast<const std::byte*>(field.get_initial_value());

  stk::mesh::entity_bytes_execute<std::byte>(field, stk::mesh::FastMeshIndex{bucketId, bucketOrd},
    [&](auto& entityBytes) {
      if (entityBytes.is_field_defined()) {
        if (initVal == nullptr) {
          std::vector<std::byte> zeroInit(entityBytes.num_bytes());

          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
            entityBytes(byte) = zeroInit[byte];
          }
        }
        else {
          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
            entityBytes(byte) = initVal[byte];
          }
        }
      }
    }
  );
}

void copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                            unsigned bucketSize,
                                            unsigned bucketId,
                                            const std::vector<FieldBase*>& allFields,
                                            const std::byte* oldAllocationAllFields,
                                            std::byte* newAllocationAllFields,
                                            const std::vector<BucketFieldSegment>& oldOffsetForField,
                                            const std::vector<BucketFieldSegment>& newOffsetForField,
                                            unsigned oldBucketCapacity,
                                            unsigned newBucketCapacity)
{
  size_t fieldIndex = 0;
  for (const FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const unsigned bytesPerEntity = field->get_internal_field_meta_data()[bucketId].m_bytesPerEntity;
      const bool oldHasAllocation = (oldOffsetForField[fieldIndex].size > 0);

      auto copy_bucket_to_bucket = [&](auto& oldBucketBytes, auto& newBucketBytes) {
        for (stk::mesh::EntityIdx entityIdx : oldBucketBytes.entities()) {
          for (stk::mesh::ByteIdx byte : oldBucketBytes.bytes()) {
            ASAN_UNPOISON_MEMORY_REGION(&newBucketBytes(entityIdx, byte), 1);
            newBucketBytes(entityIdx, byte) = oldBucketBytes(entityIdx, byte);
          }
        }
      };

      if (oldHasAllocation) {
        // We are in an in-between state where the internal FieldMetaData array for this Field is invalid,
        // and we are managing raw pointers to the old and new Bucket allocations.  Build up the BucketBytes
        // objects directly instead of going through the FieldBase::bytes() API.
        if (field->host_data_layout() == stk::mesh::Layout::Right) {
          BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right> oldBucketBytes(
                oldAllocationAllFields + oldOffsetForField[fieldIndex].offset, bytesPerEntity,
                field->data_traits().alignment_of, bucketSize);
          BucketBytes<std::byte, stk::ngp::HostMemSpace, Layout::Right> newBucketBytes(
                newAllocationAllFields + newOffsetForField[fieldIndex].offset, bytesPerEntity,
                field->data_traits().alignment_of, bucketSize);

          copy_bucket_to_bucket(oldBucketBytes, newBucketBytes);
        }
        else if (field->host_data_layout() == stk::mesh::Layout::Left) {
          BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left> oldBucketBytes(
                oldAllocationAllFields + oldOffsetForField[fieldIndex].offset,
                bytesPerEntity, field->data_traits().alignment_of, bucketSize, oldBucketCapacity);
          BucketBytes<std::byte, stk::ngp::HostMemSpace, Layout::Left> newBucketBytes(
                newAllocationAllFields + newOffsetForField[fieldIndex].offset,
                bytesPerEntity, field->data_traits().alignment_of, bucketSize, newBucketCapacity);

          copy_bucket_to_bucket(oldBucketBytes, newBucketBytes);
        }
        else {
          STK_ThrowErrorMsg("Unsupported Field host data layout: " << field->host_data_layout());
        }
      }

      ++fieldIndex;
    }
  }
}

void initialize_new_field_values(FieldBase& newField, const EntityRank rank, const unsigned bucketId,
                                 unsigned size, unsigned capacity)
{
  if (newField.entity_rank() == rank) {
    FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(newField.get_internal_field_meta_data()[bucketId]);
    initialize_field_on_bucket(newField, bucketId, fieldMetaData, size, capacity);
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

void update_field_meta_data(const EntityRank rank, const unsigned bucketId, const std::vector<FieldBase*> & allFields,
                            const PartVector & supersetParts, unsigned bucketSize, unsigned bucketCapacity)
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


FieldDataManager::FieldDataManager(const unsigned num_ranks,
                                   unsigned alignmentIncrementBytes,
                                   std::unique_ptr<AllocatorAdaptorInterface> allocatorAdaptor)
  : m_fieldDataAllocator(std::move(allocatorAdaptor)),
    m_alignmentIncrementBytes(alignmentIncrementBytes),
    m_fieldRawData(num_ranks),
    m_bucketCapacity(num_ranks),
    m_numBytesAllocatedPerField()
{
  if (not m_fieldDataAllocator) {
    m_fieldDataAllocator = std::make_unique<AllocatorAdaptor<stk::impl::FieldDataAllocator<std::byte>>>();
  }
}

void
FieldDataManager::resize_bucket_arrays(const EntityRank rank, const std::vector<FieldBase*>& allFields,
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
FieldDataManager::allocate_bucket_ordinal_field_data(const EntityRank rank,
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
        size_t fieldDataSizeThisBucket =
            stk::adjust_up_to_alignment_boundary(static_cast<size_t>(layout.numBytesPerEntity)*capacity,
                                                 m_alignmentIncrementBytes);
        totalFieldDataSize += fieldDataSizeThisBucket;
        m_numBytesAllocatedPerField[i] += fieldDataSizeThisBucket;
      }
      fields[i]->get_internal_field_meta_data()[bucketOrd] = newFieldMetaData;
      fields[i]->update_cached_field_meta_data();
    }
  }

  if (totalFieldDataSize > 0) {
    std::byte* allData = m_fieldDataAllocator->allocate(totalFieldDataSize);
    m_fieldRawData[rank][bucketOrd] = allData;

    size_t currentFieldOffset = 0;
    for (stk::mesh::FieldBase* field : fields) {
      if (field->entity_rank() == rank) {
        FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
        FieldMetaData& fieldMetaData = fieldMetaDataArray[bucketOrd];
        update_field_pointer(fieldMetaData, capacity, currentFieldOffset, allData, m_alignmentIncrementBytes);
        initialize_field_on_bucket(*field, bucketOrd, fieldMetaData, size, capacity);
      }
    }
  }
  else {
    m_fieldRawData[rank][bucketOrd] = nullptr;
  }
}

void
FieldDataManager::allocate_bucket_field_data(const EntityRank rank,
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

void
FieldDataManager::allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId,
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
FieldDataManager::get_old_bucket_field_offsets(const EntityRank rank,
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
FieldDataManager::get_new_bucket_field_offsets(const EntityRank rank,
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
FieldDataManager::reallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
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

    copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields,
                                           oldAllocationAllFields, newAllocationAllFields,
                                           oldOffsetForField, newOffsetForField,
                                           bucketCapacity, bucketCapacity);

    m_fieldDataAllocator->deallocate(m_fieldRawData[rank][bucketId], oldBucketAllocationSize);
    m_fieldRawData[rank][bucketId] = newAllocationAllFields;
    m_bucketCapacity[rank][bucketId] = bucketCapacity;

    update_field_pointers_to_new_bucket(rank, bucketId, allFields, bucketCapacity, newAllocationAllFields,
                                        m_alignmentIncrementBytes);
    initialize_new_field_values(newField, rank, bucketId, bucketSize, bucketCapacity);
  }
}

void
FieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
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
FieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
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
      FieldMetaDataArrayType newFieldMetaDataArray(Kokkos::view_alloc("FieldMetaDataArray_" + fields[i]->name(), Kokkos::WithoutInitializing),
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
FieldDataManager::allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
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
FieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
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
FieldDataManager::remove_field_data_for_entity(EntityRank rank, unsigned bucketId,
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
FieldDataManager::initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd,
                                               unsigned newBucketSize, const std::vector<FieldBase *> &fields)
{
  for (size_t i = 0; i < fields.size(); ++i) {
    const FieldBase& field = *fields[i];
    if (field.entity_rank() == rank) {
      FieldMetaData& fieldMetaData = fields[i]->get_internal_field_meta_data()[bucketId];
      const int numBytesPerEntity = fieldMetaData.m_bytesPerEntity;

      if (numBytesPerEntity > 0) {
        fieldMetaData.m_bucketSize = newBucketSize;
        initialize_field_on_entity(field, bucketId, bucketOrd);
      }
    }
  }
}

void
FieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                            unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize)
{
  initialize_entity_field_data(dstRank, dstBucketId, dstBucketOrd, newBucketSize, allFields);
}

void
FieldDataManager::grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                                       unsigned bucketSize, unsigned newBucketCapacity)
{
  const int oldBucketCapacity = m_bucketCapacity[rank][bucketId];
  m_bucketCapacity[rank][bucketId] = newBucketCapacity;

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, allFields,
                                                                                   newBucketCapacity);
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

  copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields,
                                         oldAllocationAllFields, newAllocationAllFields,
                                         oldOffsetForField, newOffsetForField,
                                         oldBucketCapacity, newBucketCapacity);

  m_fieldDataAllocator->deallocate(m_fieldRawData[rank][bucketId], oldBucketAllocationSize);
  m_fieldRawData[rank][bucketId] = newAllocationAllFields;

  update_field_pointers_to_new_bucket(rank, bucketId, allFields, newBucketCapacity, newAllocationAllFields,
                                      m_alignmentIncrementBytes);
}

void
FieldDataManager::reset_empty_field_data(EntityRank /*rank*/, unsigned bucketId, unsigned bucketSize,
                                         unsigned bucketCapacity, const FieldVector & fields)
{
  for (const FieldBase * field : fields) {
    const FieldMetaData & fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    if (field->host_data_layout() == Layout::Right) {
      ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + bucketSize * fieldMetaData.m_bytesPerEntity,
                                (bucketCapacity - bucketSize) * fieldMetaData.m_bytesPerEntity);
    }
    else if (field->host_data_layout() == Layout::Left) {
      const int numScalars = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity;
      const int bytesPerScalar = field->data_traits().alignment_of;
      for (int scalar = 0; scalar < numScalars; ++scalar) {
        ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + scalar*bucketCapacity*bytesPerScalar +
                                  bucketSize*bytesPerScalar,
                                  (bucketCapacity - bucketSize) * bytesPerScalar);

      }
    }
    else {
      STK_ThrowErrorMsg("Unsupported Field host data layout: " << field->host_data_layout());
    }
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

}
}
