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
#include <stk_mesh/base/FieldBase.hpp>  // for FieldMetaData, etc
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

void check_field_rank(const FieldBase* field, EntityRank expectedRank)
{
  STK_ThrowAssertMsg(field->entity_rank() == expectedRank, "Processing Field '" << field->name() << "' of rank " <<
                     field->entity_rank() << " when expecting only Fields with rank " << expectedRank << ".");
}

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
                          std::byte* allData, size_t alignmentPaddingSize)
{
  if (fieldMetaData.m_bytesPerEntity > 0)
  {
    currentFieldOffset = stk::adjust_up_to_alignment_boundary(currentFieldOffset, alignmentPaddingSize);
    fieldMetaData.m_data = allData + currentFieldOffset;
    currentFieldOffset += fieldMetaData.m_bytesPerEntity * capacity;
  }
}

void update_field_pointers_to_new_bucket(const EntityRank rank,
                                         const unsigned bucketId,
                                         const std::vector<FieldBase*>& fieldsOfRank,
                                         const size_t capacity,
                                         std::byte* newAllocationAllFields,
                                         unsigned alignmentPaddingSize)
{
  size_t currentFieldOffset = 0;
  for (FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    fieldMetaData.m_bucketCapacity = capacity;
    update_field_pointer(fieldMetaData, capacity, currentFieldOffset, newAllocationAllFields, alignmentPaddingSize);
  }
}

void initialize_field_on_bucket(const stk::mesh::FieldBase& field, int bucketId, const FieldMetaData& fieldMetaData,
                                unsigned size, unsigned capacity)
{
  // Poison it all, and then unpoison each byte before writing to it below
  ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data, capacity * fieldMetaData.m_bytesPerEntity);

  const auto& initVal = field.get_initial_value_bytes();

  stk::mesh::bucket_bytes_execute<std::byte>(field, bucketId,
    [&](auto& bucketBytes) {
      if (bucketBytes.is_field_defined()) {
        if (initVal.extent(0) == 0) {
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
              bucketBytes(entity, byte) = initVal(byte());
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
      constexpr std::byte zeroInit{0};

      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
        entityBytes(byte) = zeroInit;
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
  const Kokkos::View<std::byte*, stk::ngp::HostPinnedSpace>& initVal = field.get_initial_value_bytes();

  stk::mesh::entity_bytes_execute<std::byte>(field, stk::mesh::FastMeshIndex{bucketId, bucketOrd},
    [&](auto& entityBytes) {
      if (entityBytes.is_field_defined()) {
        if (initVal.extent(0) == 0) {
          constexpr std::byte zeroInit{0};

          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
            entityBytes(byte) = zeroInit;
          }
        }
        else {
          STK_ThrowRequireMsg(static_cast<int>(initVal.extent(0)) >= entityBytes.num_bytes(),
              "Field "<<field.name()<<"'s get_initial_value_bytes() returns view of size "
              <<initVal.extent(0)<<" but entityBytes.num_bytes() = "<<entityBytes.num_bytes());

          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            ASAN_UNPOISON_MEMORY_REGION(&entityBytes(byte), 1);
            entityBytes(byte) = initVal(byte());
          }
        }
      }
    }
  );
}

void copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                            unsigned bucketSize,
                                            unsigned bucketId,
                                            const std::vector<FieldBase*>& fieldsOfRank,
                                            const std::byte* oldAllocationAllFields,
                                            std::byte* newAllocationAllFields,
                                            const std::vector<BucketFieldSegment>& oldOffsetForField,
                                            const std::vector<BucketFieldSegment>& newOffsetForField,
                                            unsigned oldBucketCapacity,
                                            unsigned newBucketCapacity)
{
  size_t fieldIndex = 0;
  for (const FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
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
        BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Right> oldBucketBytes(
              oldAllocationAllFields + oldOffsetForField[fieldIndex].offset, bytesPerEntity,
              field->data_traits().alignment_of, bucketSize);
        BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Right> newBucketBytes(
              newAllocationAllFields + newOffsetForField[fieldIndex].offset, bytesPerEntity,
              field->data_traits().alignment_of, bucketSize);

        copy_bucket_to_bucket(oldBucketBytes, newBucketBytes);
      }
      else if (field->host_data_layout() == stk::mesh::Layout::Left) {
        BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Left> oldBucketBytes(
              oldAllocationAllFields + oldOffsetForField[fieldIndex].offset,
              bytesPerEntity, field->data_traits().alignment_of, bucketSize, oldBucketCapacity);
        BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Left> newBucketBytes(
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

void initialize_new_field_values(FieldBase& newField, const EntityRank rank, const unsigned bucketId,
                                 unsigned size, unsigned capacity)
{
  FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(newField.get_internal_field_meta_data()[bucketId]);
  initialize_field_on_bucket(newField, bucketId, fieldMetaData, size, capacity);
}

void resize_field_meta_data(FieldBase& field, int newSize)
{
  FieldMetaDataArrayType& fieldMetaDataArray = field.get_internal_field_meta_data();
  if (fieldMetaDataArray.capacity() == 0u) {
    fieldMetaDataArray = FieldMetaDataArrayType("fieldMetaDataArray_" + field.name(), newSize);
  }
  else {
    fieldMetaDataArray.resize_scale(newSize);
  }
}

void update_field_meta_data(const EntityRank rank, const unsigned bucketId, const std::vector<FieldBase*> & allFields,
                            const PartVector & supersetParts, unsigned bucketSize, unsigned bucketCapacity)
{
  for (const FieldBase* field : allFields) {
    check_field_rank(field, rank);
    const FieldLayoutData layout = get_field_layout_data(*field, rank, supersetParts);
    FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
    fieldMetaData.m_numComponentsPerEntity = layout.numComponents;
    fieldMetaData.m_numCopiesPerEntity = layout.numCopies;
    fieldMetaData.m_bucketSize = bucketSize;
    fieldMetaData.m_bucketCapacity = bucketCapacity;
  }
}

FieldDataManager::FieldDataManager(const unsigned numRanks,
                                   unsigned alignmentPaddingSize)
  : m_fieldDataAllocator(),
    m_alignmentPaddingSize(alignmentPaddingSize),
    m_fieldRawData(numRanks),
    m_bucketCapacity(numRanks),
    m_numBytesAllocatedPerField()
{
}

void
FieldDataManager::resize_bucket_arrays(const EntityRank rank, const std::vector<FieldBase*>& fieldsOfRank,
                                       int newNumBuckets)
{
  if (rank >= static_cast<int>(m_fieldRawData.size())) {
    m_fieldRawData.resize(rank+1);
    m_bucketCapacity.resize(rank+1);
  }

  m_fieldRawData[rank].resize(newNumBuckets);
  m_bucketCapacity[rank].resize(newNumBuckets);

  for (FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);

    FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
    if (fieldMetaDataArray.capacity() == 0u) {
      fieldMetaDataArray = FieldMetaDataArrayType("FieldMetaDataArray_" + field->name(), newNumBuckets);
    }
    else {
      fieldMetaDataArray.resize_scale(newNumBuckets);
    }
  }
}

void
FieldDataManager::allocate_bucket_ordinal_field_data(const EntityRank rank,
                                                     const std::vector<FieldBase *>& fieldsOfRank,
                                                     const PartVector& supersetParts,
                                                     unsigned totalNumFields,
                                                     unsigned bucketOrd,
                                                     unsigned size,
                                                     unsigned capacity)
{
  m_bucketCapacity[rank][bucketOrd] = capacity;
  m_numBytesAllocatedPerField.resize(totalNumFields, 0);
  size_t totalFieldDataSize = 0;

  for (FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);

    const unsigned fieldOrdinal = field->mesh_meta_data_ordinal();
    FieldMetaData newFieldMetaData;

    newFieldMetaData.m_bucketSize = size;
    newFieldMetaData.m_bucketCapacity = capacity;

    const FieldLayoutData layout = get_field_layout_data(*field, rank, supersetParts);
    if (layout.numBytesPerEntity > 0) {
      newFieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
      newFieldMetaData.m_numComponentsPerEntity = layout.numComponents;
      newFieldMetaData.m_numCopiesPerEntity = layout.numCopies;
      size_t fieldDataSizeThisBucket =
          stk::adjust_up_to_alignment_boundary(static_cast<size_t>(layout.numBytesPerEntity)*capacity,
                                               m_alignmentPaddingSize);
      totalFieldDataSize += fieldDataSizeThisBucket;
      m_numBytesAllocatedPerField[fieldOrdinal] += fieldDataSizeThisBucket;
    }
    field->get_internal_field_meta_data()[bucketOrd] = newFieldMetaData;
    field->update_cached_field_meta_data();
  }

  if (totalFieldDataSize > 0) {
    auto allData = m_fieldDataAllocator.host_allocate(totalFieldDataSize);
    m_fieldRawData[rank][bucketOrd] = allData;

    size_t currentFieldOffset = 0;
    for (FieldBase* field : fieldsOfRank) {
      FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
      FieldMetaData& fieldMetaData = fieldMetaDataArray[bucketOrd];
      update_field_pointer(fieldMetaData, capacity, currentFieldOffset, allData.data(), m_alignmentPaddingSize);
      initialize_field_on_bucket(*field, bucketOrd, fieldMetaData, size, capacity);
    }
  }
  else {
    m_fieldRawData[rank][bucketOrd] = AllocationType();
  }
}

void
FieldDataManager::allocate_bucket_field_data(const EntityRank rank,
                                             const std::vector<FieldBase*>& fieldsOfRank,
                                             const PartVector& supersetParts,
                                             unsigned totalNumFields,
                                             unsigned size,
                                             unsigned capacity)
{
  const unsigned newNumBuckets = (static_cast<int>(m_bucketCapacity.size()) <= rank) ? 1
                                                                                     : m_bucketCapacity[rank].size() + 1;
  resize_bucket_arrays(rank, fieldsOfRank, newNumBuckets);

  const unsigned bucketOrd = newNumBuckets - 1;  // New Bucket at the end of the list
  allocate_bucket_ordinal_field_data(rank, fieldsOfRank, supersetParts, totalNumFields, bucketOrd, size, capacity);
}

void
FieldDataManager::allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId,
                                               const std::vector<FieldBase*>& fieldsOfRank)
{
  for (FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    const unsigned currentSize = field->get_internal_field_meta_data().size();
    if (bucketId >= currentSize) {
      resize_field_meta_data(*field, currentSize+1);
      field->update_cached_field_meta_data();
    }
  }
}

std::vector<BucketFieldSegment>
FieldDataManager::get_old_bucket_field_offsets(const EntityRank rank,
                                               const unsigned bucketId,
                                               const std::vector<FieldBase*>& fieldsOfRank,
                                               const unsigned capacity) const
{
  const auto& oldAllocationStart = m_fieldRawData[rank][bucketId];

  std::vector<BucketFieldSegment> oldOffsetForField;
  oldOffsetForField.reserve(fieldsOfRank.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase* field : fieldsOfRank) {
    const bool isFieldValid = (field != nullptr);
    check_field_rank(field, rank);

    if (isFieldValid) {
      const FieldMetaDataArrayType& fieldMetaDataArray = field->get_internal_field_meta_data();
      const bool isBucketInRange = (bucketId < fieldMetaDataArray.size());
      const bool hasAllocation = (isBucketInRange) ? (fieldMetaDataArray[bucketId].m_data != nullptr) : false;
      const size_t oldOffsetIntoBucket = (hasAllocation) ? fieldMetaDataArray[bucketId].m_data -
                                                           oldAllocationStart.data()
                                                         : 0u;
      const size_t oldBytesPerEntity = (hasAllocation) ? fieldMetaDataArray[bucketId].m_bytesPerEntity : 0;
      const size_t oldSizeThisBucket = stk::adjust_up_to_alignment_boundary(oldBytesPerEntity*capacity,
                                                                            m_alignmentPaddingSize);

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
                                               const std::vector<FieldBase*>& fieldsOfRank,
                                               const unsigned capacity) const
{
  std::vector<BucketFieldSegment> newOffsetForField;
  newOffsetForField.reserve(fieldsOfRank.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    const FieldMetaData & fieldMetaDataForBucket = field->get_internal_field_meta_data()[bucketId];
    size_t newSizeThisBucket =
        stk::adjust_up_to_alignment_boundary(static_cast<size_t>(fieldMetaDataForBucket.m_bytesPerEntity)*capacity,
                                             m_alignmentPaddingSize);

    newOffsetForField.emplace_back(totalAllocationSize, newSizeThisBucket);
    totalAllocationSize += newSizeThisBucket;
  }

  newOffsetForField.emplace_back(totalAllocationSize, 0);  // Add a one-past-the-end entry for bookkeeping

  return newOffsetForField;
}

void
FieldDataManager::reallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId,
                                               FieldBase & targetField,
                                               const std::vector<FieldBase*> & fieldsOfRank,
                                               const PartVector& supersetParts, unsigned bucketSize,
                                               unsigned bucketCapacity)
{
  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, fieldsOfRank,
                                                                                   bucketCapacity);
  const int oldBucketAllocationSize = oldOffsetForField.back().offset;

  allocate_new_field_meta_data(rank, bucketId, fieldsOfRank);
  update_field_meta_data(rank, bucketId, fieldsOfRank, supersetParts, bucketSize, bucketCapacity);

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, fieldsOfRank,
                                                                                   bucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;
  const unsigned numberOfStates = targetField.number_of_states();
  const unsigned fieldOrdinal = targetField.mesh_meta_data_ordinal();
  const unsigned fieldRankedOrdinal = targetField.field_ranked_ordinal();

  // What is the code path for an early field into a late part?  Why does this not increase the size?
  if (newBucketAllocationSize > oldBucketAllocationSize) {
    for (unsigned state = 0; state < numberOfStates; ++state) {
      const BucketFieldSegment & lastFieldSegment = newOffsetForField[fieldRankedOrdinal+state];
      m_numBytesAllocatedPerField[fieldOrdinal+state] += lastFieldSegment.size;
    }

    auto newAllocationAllFields = m_fieldDataAllocator.host_allocate(newBucketAllocationSize);
    const auto& oldAllocationAllFields = m_fieldRawData[rank][bucketId];

    copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, fieldsOfRank,
                                           oldAllocationAllFields.data(), newAllocationAllFields.data(),
                                           oldOffsetForField, newOffsetForField,
                                           bucketCapacity, bucketCapacity);

    m_fieldRawData[rank][bucketId] = newAllocationAllFields;
    m_bucketCapacity[rank][bucketId] = bucketCapacity;

    update_field_pointers_to_new_bucket(rank, bucketId, fieldsOfRank, bucketCapacity, newAllocationAllFields.data(),
                                        m_alignmentPaddingSize);
    for (unsigned state = 0; state < numberOfStates; ++state) {
      FieldBase* fieldOfState = targetField.field_state(static_cast<FieldState>(state));
      initialize_new_field_values(*fieldOfState, rank, bucketId, bucketSize, bucketCapacity);
    }
  }
}

void
FieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                               const std::vector<FieldBase*>& fieldsOfRank)
{
  if (fieldsOfRank.empty()) {
    return;
  }

  if (m_fieldRawData[rank][bucketId].is_allocated()) {
    for (FieldBase* field : fieldsOfRank) {
      if (field == nullptr || field->get_internal_field_meta_data().size() <= bucketId) {
        continue;
      }

      FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
      if (fieldMetaData.m_data != nullptr) {
        const size_t bytes_to_delete_this_field =
            stk::adjust_up_to_alignment_boundary(static_cast<size_t>(fieldMetaData.m_bytesPerEntity)*capacity,
                                                 m_alignmentPaddingSize);
        const unsigned fieldOrdinal = field->mesh_meta_data_ordinal();
        m_numBytesAllocatedPerField[fieldOrdinal] -= bytes_to_delete_this_field;
        fieldMetaData = FieldMetaData{};
      }
    }

    m_fieldRawData[rank][bucketId] = AllocationType();
    m_bucketCapacity[rank][bucketId] = 0;
  }
}

void
FieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fieldsOfRank,
                                            const std::vector<unsigned>& reorderedBucketIds)
{
  std::vector<AllocationType> fieldRawData(reorderedBucketIds.size());
  std::vector<unsigned> bucketCapacity(reorderedBucketIds.size());
  for (unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m) {
    fieldRawData[m] = m_fieldRawData[rank][reorderedBucketIds[m]];
    bucketCapacity[m] = m_bucketCapacity[rank][reorderedBucketIds[m]];
  }
  m_fieldRawData[rank].swap(fieldRawData);
  m_bucketCapacity[rank].swap(bucketCapacity);

  for (FieldBase* field : fieldsOfRank) {
    FieldMetaDataArrayType newFieldMetaDataArray("FieldMetaDataArray_" + field->name());
    newFieldMetaDataArray.reserve(reorderedBucketIds.size());
    FieldMetaDataArrayType& oldFieldMetaDataArray = field->get_internal_field_meta_data();
    for (unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m) {
      newFieldMetaDataArray.push_back(oldFieldMetaDataArray[reorderedBucketIds[m]]);
    }
    std::swap(oldFieldMetaDataArray, newFieldMetaDataArray);
    field->update_cached_field_meta_data();
  }
}

void
FieldDataManager::allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                      const std::vector<FieldBase*>& fieldsOfRank, unsigned totalNumFields)
{
  m_numBytesAllocatedPerField.resize(totalNumFields, 0);
  resize_bucket_arrays(rank, fieldsOfRank, buckets.size());

  for (Bucket* bucket : buckets) {
    const PartVector& supersetParts = bucket->supersets();
    allocate_bucket_ordinal_field_data(rank, fieldsOfRank, supersetParts, totalNumFields, bucket->bucket_id(),
                                       bucket->size(), bucket->capacity());
  }
}

void
FieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                        FieldBase& targetField, const std::vector<FieldBase*>& fieldsOfRank,
                                        unsigned totalNumFields)
{
  m_numBytesAllocatedPerField.resize(totalNumFields, 0);
  for (Bucket* bucket : buckets) {
    const PartVector& supersetParts = bucket->supersets();
    reallocate_bucket_field_data(rank, bucket->bucket_id(), targetField, fieldsOfRank, supersetParts,
                                 bucket->size(), bucket->capacity());
  }
}

void
FieldDataManager::remove_field_data_for_entity(EntityRank rank, unsigned bucketId,
                                               unsigned /*bucketOrd*/, unsigned newBucketSize,
                                               const std::vector<FieldBase*>& fieldsOfRank)
{
  for (const FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    field->get_internal_field_meta_data()[bucketId].m_bucketSize = newBucketSize;
  }
}

void
FieldDataManager::initialize_entity_field_data(const std::vector<FieldBase*>& fieldsOfRank, EntityRank rank,
                                               unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize)
{
  for (FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
    fieldMetaData.m_bucketSize = newBucketSize;

    const int numBytesPerEntity = fieldMetaData.m_bytesPerEntity;
    if (numBytesPerEntity > 0) {
      initialize_field_on_entity(*field, bucketId, bucketOrd);
    }
  }
}

void
FieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *>& fieldsOfRank, EntityRank rank,
                                            unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize)
{
  for (const FieldBase* field : fieldsOfRank) {
    check_field_rank(field, rank);
    field->get_internal_field_meta_data()[bucketId].m_bucketSize = newBucketSize;
  }
}

void
FieldDataManager::grow_bucket_capacity(const FieldVector& fieldsOfRank, EntityRank rank, unsigned bucketId,
                                       unsigned bucketSize, unsigned newBucketCapacity)
{
  const int oldBucketCapacity = m_bucketCapacity[rank][bucketId];
  m_bucketCapacity[rank][bucketId] = newBucketCapacity;

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, fieldsOfRank,
                                                                                   newBucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;

  if (newBucketAllocationSize == 0) {
    for (const stk::mesh::FieldBase* field : fieldsOfRank) {
      check_field_rank(field, rank);
      FieldMetaData& fieldMetaData = field->get_internal_field_meta_data()[bucketId];
      fieldMetaData.m_bucketCapacity = newBucketCapacity;
    }

    return;
  }

  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, fieldsOfRank,
                                                                                   oldBucketCapacity);

  unsigned i = 0;
  for (const stk::mesh::FieldBase* field : fieldsOfRank) {
    m_numBytesAllocatedPerField[field->mesh_meta_data_ordinal()] += newOffsetForField[i].size - oldOffsetForField[i].size;
    ++i;
  }

  auto newAllocationAllFields = m_fieldDataAllocator.host_allocate(newBucketAllocationSize);
  const auto& oldAllocationAllFields = m_fieldRawData[rank][bucketId];

  copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, fieldsOfRank,
                                         oldAllocationAllFields.data(), newAllocationAllFields.data(),
                                         oldOffsetForField, newOffsetForField,
                                         oldBucketCapacity, newBucketCapacity);

  m_fieldRawData[rank][bucketId] = newAllocationAllFields;

  update_field_pointers_to_new_bucket(rank, bucketId, fieldsOfRank, newBucketCapacity, newAllocationAllFields.data(),
                                      m_alignmentPaddingSize);
}

void
FieldDataManager::reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                         unsigned bucketCapacity, const FieldVector& fieldsOfRank)
{
  for (const FieldBase * field : fieldsOfRank) {
    check_field_rank(field, rank);
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
