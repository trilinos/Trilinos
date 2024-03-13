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

size_t getFieldBucketSizeInBytes(const FieldMetaDataVector& field_meta_data_vector, const unsigned bucket_id, const unsigned char *end_of_field);

struct FieldLayoutData
{
  int numBytesPerEntity;
  int firstDimension;
};

FieldLayoutData getFieldLayoutData(const FieldBase &field, const EntityRank rank, const PartVector& superset_parts)
{
  FieldLayoutData layout{0, 0};

  const FieldBase::Restriction & restriction = find_and_check_restriction(field, rank, superset_parts);

  if (restriction.num_scalars_per_entity() > 0) {
    const unsigned type_stride = field.data_traits().stride_of;
    layout.numBytesPerEntity = type_stride * restriction.num_scalars_per_entity();
    layout.firstDimension = restriction.dimension();
  }

  return layout;
}

void updateFieldPointer(FieldMetaData& field_meta_data, const size_t capacity, size_t &current_field_offset, unsigned char* all_data, size_t alignment_increment_bytes)
{
    if (field_meta_data.m_bytesPerEntity > 0)
    {
        current_field_offset = stk::adjust_up_to_alignment_boundary(current_field_offset, alignment_increment_bytes);
        field_meta_data.m_data = all_data + current_field_offset;
        current_field_offset += field_meta_data.m_bytesPerEntity * capacity;
    }
}

void initializeField(FieldMetaData& field_meta_data, const unsigned char* init_val, [[maybe_unused]] unsigned size,
                     [[maybe_unused]] unsigned capacity)
{
#ifdef STK_ASAN_IS_ON
  if (field_meta_data.m_bytesPerEntity > 0) {
    ASAN_UNPOISON_MEMORY_REGION(field_meta_data.m_data, size * field_meta_data.m_bytesPerEntity);
    if (init_val != nullptr) {
      for (unsigned j = 0; j < size; ++j) {
        std::memcpy(field_meta_data.m_data + j * field_meta_data.m_bytesPerEntity, init_val,
                    field_meta_data.m_bytesPerEntity);
      }
    }
    else {
      std::memset(field_meta_data.m_data, 0, size * field_meta_data.m_bytesPerEntity);
    }
  }
#else
  if (field_meta_data.m_bytesPerEntity > 0) {
    if (init_val != nullptr) {
      for (unsigned j = 0; j < capacity; ++j) {
        std::memcpy(field_meta_data.m_data + j * field_meta_data.m_bytesPerEntity, init_val,
                    field_meta_data.m_bytesPerEntity);
      }
    }
    else {
      std::memset(field_meta_data.m_data, 0, capacity * field_meta_data.m_bytesPerEntity);
    }
  }
#endif
}

void setInitialValue(unsigned char* data_location, const FieldBase& field, const int numBytesPerEntity)
{
  ASAN_UNPOISON_MEMORY_REGION(data_location, numBytesPerEntity);
  const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
  if (init_val != nullptr) {
    std::memcpy( data_location, init_val, numBytesPerEntity );
  }
  else {
    std::memset( data_location, 0, numBytesPerEntity );
  }
}

FieldDataManager::FieldDataManager(unsigned alignmentIncrementBytes,
                                   std::unique_ptr<AllocatorAdaptorInterface> allocatorAdaptor)
  : m_fieldDataAllocator(std::move(allocatorAdaptor)),
    alignment_increment_bytes(alignmentIncrementBytes)
{
  if (not m_fieldDataAllocator) {
    m_fieldDataAllocator = std::make_unique<AllocatorAdaptor<stk::impl::FieldDataAllocator<unsigned char>>>();
  }
}

void DefaultFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
                                                         const std::vector<FieldBase *> & fields,
                                                         const PartVector& superset_parts,
                                                         unsigned size,
                                                         unsigned capacity)
{
  if (m_field_raw_data.empty()) {
    STK_ThrowRequireMsg(!fields.empty(),"allocate_bucket_field_data ERROR, field-data-manager was constructed with 0 entity-ranks, and there are no fields. Mesh has not been initialized correctly.");

    for (size_t i=0; i<fields.size(); ++i) {
      if (fields[i] != nullptr) {
        m_field_raw_data.resize(fields[i]->get_mesh().mesh_meta_data().entity_rank_count());
        m_bucketCapacity.resize(fields[i]->get_mesh().mesh_meta_data().entity_rank_count());
        break;
      }
    }
  }

  m_bucketCapacity[rank].push_back(capacity);

  if (m_num_bytes_allocated_per_field.empty()) {
    m_num_bytes_allocated_per_field.resize(fields.size(), 0);
  }

  size_t num_fields = fields.size();
  // Sizing loop
  size_t total_field_data_size = 0;

  for (size_t i = 0; i < num_fields; ++i) {
    FieldMetaData field_meta_data;

    const FieldBase & field = *fields[i];
    if (field.entity_rank() == rank) {
      const FieldLayoutData layout = getFieldLayoutData(field, rank, superset_parts);
      if (layout.numBytesPerEntity > 0) {
        field_meta_data.m_bytesPerEntity = layout.numBytesPerEntity;
        field_meta_data.m_firstDimension = layout.firstDimension;
        size_t field_data_size_this_bucket = stk::adjust_up_to_alignment_boundary(layout.numBytesPerEntity*capacity,
                                                                                  alignment_increment_bytes);
        total_field_data_size += field_data_size_this_bucket;
        m_num_bytes_allocated_per_field[i] += field_data_size_this_bucket;
      }
      fields[i]->get_meta_data_for_field().push_back(field_meta_data);
    }
  }

  // Allocate all field data for this bucket
  if (total_field_data_size > 0) {
    unsigned char* all_data = m_fieldDataAllocator->allocate(total_field_data_size);
    m_field_raw_data[rank].push_back(all_data);

    // Set data ptrs in field meta datas
    size_t current_field_offset = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
      const FieldBase & field = *fields[i];
      if (field.entity_rank() == rank) {
        const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
        FieldMetaData& field_meta_data = const_cast<FieldMetaData&>(field.get_meta_data_for_field().back());
        updateFieldPointer(field_meta_data, capacity, current_field_offset, all_data, alignment_increment_bytes);
        initializeField(field_meta_data, init_val, size, capacity);
      }
    }
  }
  else {
    m_field_raw_data[rank].push_back(nullptr);
  }
}

void DefaultFieldDataManager::allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId, const std::vector<FieldBase*>& allFields)
{
    for (FieldBase* field : allFields) {
        if (field->entity_rank() == rank) {
            if (bucketId >= field->get_meta_data_for_field().size()) {
                FieldMetaData fieldMetaData;
                field->get_meta_data_for_field().push_back(fieldMetaData);
            }
        }
    }
}

std::vector<BucketFieldSegment> DefaultFieldDataManager::get_old_bucket_field_offsets(const EntityRank rank,
                                                                                      const unsigned bucketId,
                                                                                      const std::vector<FieldBase*>& allFields,
                                                                                      const unsigned capacity) const
{
  const unsigned char* oldAllocationStart = m_field_raw_data[rank][bucketId];

  std::vector<BucketFieldSegment> oldOffsetForField;
  oldOffsetForField.reserve(allFields.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase * field : allFields) {
    const bool isFieldValid = (field != nullptr);
    const bool isTargetRank = (field->entity_rank() == rank);

    if (isFieldValid && isTargetRank) {
      const std::vector<FieldMetaData> & fieldMetaData = field->get_meta_data_for_field();
      const bool isBucketInRange = (bucketId < fieldMetaData.size());
      const bool hasAllocation = (isBucketInRange) ? (fieldMetaData[bucketId].m_data != nullptr) : false;
      const int oldOffsetIntoBucket = (hasAllocation) ? fieldMetaData[bucketId].m_data - oldAllocationStart : 0;
      const int oldBytesPerEntity = (hasAllocation) ? fieldMetaData[bucketId].m_bytesPerEntity : 0;
      const int oldSizeThisBucket = stk::adjust_up_to_alignment_boundary(oldBytesPerEntity*capacity, alignment_increment_bytes);

      oldOffsetForField.emplace_back(oldOffsetIntoBucket, oldSizeThisBucket);
      totalAllocationSize += oldSizeThisBucket;
    }
  }

  oldOffsetForField.emplace_back(totalAllocationSize, 0);  // Add a one-past-the-end entry for bookkeeping

  return oldOffsetForField;
}

std::vector<BucketFieldSegment> DefaultFieldDataManager::get_new_bucket_field_offsets(const EntityRank rank,
                                                                                      const unsigned bucketId,
                                                                                      const std::vector<FieldBase*>& allFields,
                                                                                      const unsigned capacity) const
{
  std::vector<BucketFieldSegment> newOffsetForField;
  newOffsetForField.reserve(allFields.size()+1);
  int totalAllocationSize = 0;

  for (const FieldBase * field : allFields) {
    if (field->entity_rank() == rank) {
      const FieldMetaData & fieldMetaDataForBucket = field->get_meta_data_for_field()[bucketId];
      size_t newSizeThisBucket = stk::adjust_up_to_alignment_boundary(fieldMetaDataForBucket.m_bytesPerEntity * capacity,
                                                                      alignment_increment_bytes);

      newOffsetForField.emplace_back(totalAllocationSize, newSizeThisBucket);
      totalAllocationSize += newSizeThisBucket;
    }
  }

  newOffsetForField.emplace_back(totalAllocationSize, 0);  // Add a one-past-the-end entry for bookkeeping

  return newOffsetForField;
}

void DefaultFieldDataManager::update_field_meta_data(const EntityRank rank, const unsigned bucketId,
                                                   const std::vector<FieldBase*> & allFields,
                                                   const PartVector & supersetParts)
{
  for (const FieldBase * field : allFields) {
    if (field->entity_rank() == rank) {
      const FieldLayoutData layout = getFieldLayoutData(*field, rank, supersetParts);
      FieldMetaData & fieldMetaData = const_cast<FieldMetaData&>(field->get_meta_data_for_field()[bucketId]);
      fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
      fieldMetaData.m_firstDimension = layout.firstDimension;
    }
  }
}

void DefaultFieldDataManager::copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                                                     unsigned bucketSize,
                                                                     unsigned bucketId,
                                                                     const std::vector<FieldBase*>& allFields,
                                                                     const std::vector<BucketFieldSegment>& oldOffsetForField,
                                                                     const std::vector<BucketFieldSegment>& newOffsetForField,
                                                                     const unsigned char* oldAllocationAllFields,
                                                                     unsigned char* newAllocationAllFields)
{
  size_t fieldIndex = 0;
  for (const FieldBase* field : allFields) {
    if (field->entity_rank() == rank) {
      const unsigned bytesPerEntity = field->get_meta_data_for_field()[bucketId].m_bytesPerEntity;
      const bool oldHasAllocation = (oldOffsetForField[fieldIndex].size > 0);
      const unsigned oldNumBytesUsed = (oldHasAllocation) ? bucketSize * bytesPerEntity : 0;

      ASAN_UNPOISON_MEMORY_REGION(newAllocationAllFields + newOffsetForField[fieldIndex].offset, oldNumBytesUsed);
      std::memcpy(newAllocationAllFields + newOffsetForField[fieldIndex].offset,
                  oldAllocationAllFields + oldOffsetForField[fieldIndex].offset, oldNumBytesUsed);
      ++fieldIndex;
    }
  }
}

void DefaultFieldDataManager::update_field_pointers_to_new_bucket(const EntityRank rank,
                                                                  const unsigned bucketId,
                                                                  const std::vector<FieldBase*>& allFields,
                                                                  const size_t capacity,
                                                                  unsigned char* newAllocationAllFields)
{
    size_t currentFieldOffset = 0;
    for (size_t i = 0; i < allFields.size(); ++i) {
        const FieldBase& field = *allFields[i];
        if (field.entity_rank() == rank) {
            FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(field.get_meta_data_for_field()[bucketId]);
            updateFieldPointer(fieldMetaData, capacity, currentFieldOffset, newAllocationAllFields, alignment_increment_bytes);
        }
    }
}

void DefaultFieldDataManager::initialize_new_field_values(FieldBase& newField, const EntityRank rank, const unsigned bucketId,
                                                          unsigned size, unsigned capacity)
{
    if (newField.entity_rank() == rank) {
        const unsigned char* initVal = reinterpret_cast<const unsigned char*>(newField.get_initial_value());
        FieldMetaData& fieldMetaData = const_cast<FieldMetaData&>(newField.get_meta_data_for_field()[bucketId]);
        initializeField(fieldMetaData, initVal, size, capacity);
    }
}

void DefaultFieldDataManager::reallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, FieldBase & newField,
                                                           const std::vector<FieldBase *> & allFields, const PartVector& supersetParts,
                                                           unsigned bucketSize, unsigned bucketCapacity)
{
  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, allFields, bucketCapacity);
  const int oldBucketAllocationSize = oldOffsetForField.back().offset;

  allocate_new_field_meta_data(rank, bucketId, allFields);
  update_field_meta_data(rank, bucketId, allFields, supersetParts);

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, allFields, bucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;

  if (newBucketAllocationSize > oldBucketAllocationSize) {
    const BucketFieldSegment & lastFieldSegment = newOffsetForField[newOffsetForField.size()-2];
    m_num_bytes_allocated_per_field.back() += lastFieldSegment.size;

    unsigned char* newAllocationAllFields = m_fieldDataAllocator->allocate(newBucketAllocationSize);
    const unsigned char* oldAllocationAllFields = m_field_raw_data[rank][bucketId];

    copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields, oldOffsetForField, newOffsetForField,
                                           oldAllocationAllFields, newAllocationAllFields);

    m_fieldDataAllocator->deallocate(m_field_raw_data[rank][bucketId], oldBucketAllocationSize);
    m_field_raw_data[rank][bucketId] = newAllocationAllFields;
    m_bucketCapacity[rank][bucketId] = bucketCapacity;

    update_field_pointers_to_new_bucket(rank, bucketId, allFields, bucketCapacity, newAllocationAllFields);
    initialize_new_field_values(newField, rank, bucketId, bucketSize, bucketCapacity);
  }
}

void DefaultFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
        const std::vector<FieldBase*>&  fields)
{
    if(fields.empty())
        return;

    if(m_field_raw_data[rank][bucket_id] != nullptr)
    {
        size_t bytes_to_delete = 0;
        for(unsigned int i = 0; i < fields.size(); ++i)
        {
            if(fields[i] == nullptr ||
               fields[i]->entity_rank() != rank ||
               fields[i]->get_meta_data_for_field().size() <= bucket_id)
                continue;
            FieldMetaData& field_data = fields[i]->get_meta_data_for_field()[bucket_id];
            if(field_data.m_data != nullptr)
            {
                const size_t bytes_to_delete_this_field = stk::adjust_up_to_alignment_boundary(field_data.m_bytesPerEntity*capacity,
                                                                                               alignment_increment_bytes);
                m_num_bytes_allocated_per_field[i] -= bytes_to_delete_this_field;
                bytes_to_delete += bytes_to_delete_this_field;
                field_data.m_bytesPerEntity = 0;
                field_data.m_firstDimension = 0;
                field_data.m_data = nullptr;
            }
        }
        m_fieldDataAllocator->deallocate(m_field_raw_data[rank][bucket_id], bytes_to_delete);
        m_field_raw_data[rank][bucket_id] = nullptr;
        m_bucketCapacity[rank][bucket_id] = 0;
    }
}

void DefaultFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                                        const std::vector<unsigned>& reorderedBucketIds)
{
    std::vector<unsigned char*> field_raw_data(reorderedBucketIds.size());
    std::vector<unsigned> bucketCapacity(reorderedBucketIds.size());
    for(unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m)
    {
        field_raw_data[m] = m_field_raw_data[rank][reorderedBucketIds[m]];
        bucketCapacity[m] = m_bucketCapacity[rank][reorderedBucketIds[m]];
    }
    m_field_raw_data[rank].swap(field_raw_data);
    m_bucketCapacity[rank].swap(bucketCapacity);

    for(size_t i = 0; i < fields.size(); ++i)
    {
        if (fields[i]->entity_rank() == rank)
        {
            FieldMetaDataVector new_field_meta_data_vector(reorderedBucketIds.size());
            for(unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m)
            {
                new_field_meta_data_vector[m] = fields[i]->get_meta_data_for_field()[reorderedBucketIds[m]];
            }
            new_field_meta_data_vector.swap(fields[i]->get_meta_data_for_field());
        }
    }
}

void DefaultFieldDataManager::allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & fields)
{
    m_num_bytes_allocated_per_field.resize(fields.size(), 0);
    for(size_t i=0; i<buckets.size(); ++i)
    {
        const PartVector& superset_parts = buckets[i]->supersets();
        allocate_bucket_field_data(rank, fields, superset_parts, buckets[i]->size(), buckets[i]->capacity());
    }
}

void DefaultFieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField, const std::vector<FieldBase *> & allFields)
{
    m_num_bytes_allocated_per_field.resize(allFields.size(), 0);
    for(size_t i=0; i<buckets.size(); ++i)
    {
        const PartVector& superset_parts = buckets[i]->supersets();
        reallocate_bucket_field_data(rank, buckets[i]->bucket_id(), currentField, allFields, superset_parts,
                                     buckets[i]->size(), buckets[i]->capacity());
    }
}

void DefaultFieldDataManager::remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields)
{

}
void DefaultFieldDataManager::initialize_entity_field_data(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields)
{
    // bucket of bucket_id shrinks by one
    for(size_t i = 0; i < fields.size(); ++i)
    {
        const FieldBase & field = *fields[i];
        if (field.entity_rank() == rank)
        {
            const FieldMetaData& field_meta_data = fields[i]->get_meta_data_for_field()[bucket_id];
            const int num_bytes_per_entity = field_meta_data.m_bytesPerEntity;

            if(num_bytes_per_entity > 0)
            {
                setInitialValue(field_meta_data.m_data + bucket_ord * num_bytes_per_entity, field, num_bytes_per_entity);
            }
        }
    }
}

void DefaultFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dst_rank,
                                                        unsigned dst_bucket_id, unsigned dst_bucket_ord)
{
    initialize_entity_field_data(dst_rank, dst_bucket_id, dst_bucket_ord, allFields);
}

void DefaultFieldDataManager::grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                                                   unsigned bucketSize, unsigned bucketCapacity)
{
  const int oldBucketCapacity = m_bucketCapacity[rank][bucketId];
  m_bucketCapacity[rank][bucketId] = bucketCapacity;

  std::vector<BucketFieldSegment> newOffsetForField = get_new_bucket_field_offsets(rank, bucketId, allFields, bucketCapacity);
  const int newBucketAllocationSize = newOffsetForField.back().offset;

  if (newBucketAllocationSize == 0) {
    return;
  }

  std::vector<BucketFieldSegment> oldOffsetForField = get_old_bucket_field_offsets(rank, bucketId, allFields, oldBucketCapacity);
  const int oldBucketAllocationSize = oldOffsetForField.back().offset;

  unsigned i = 0;
  for (const stk::mesh::FieldBase * field : allFields) {
    if (field->entity_rank() == rank) {
      m_num_bytes_allocated_per_field[field->mesh_meta_data_ordinal()] += newOffsetForField[i].size - oldOffsetForField[i].size;
      ++i;
    }
  }

  unsigned char* newAllocationAllFields = m_fieldDataAllocator->allocate(newBucketAllocationSize);
  const unsigned char* oldAllocationAllFields = m_field_raw_data[rank][bucketId];

  copy_field_data_from_old_to_new_bucket(rank, bucketSize, bucketId, allFields, oldOffsetForField, newOffsetForField,
                                         oldAllocationAllFields, newAllocationAllFields);

  m_fieldDataAllocator->deallocate(m_field_raw_data[rank][bucketId], oldBucketAllocationSize);
  m_field_raw_data[rank][bucketId] = newAllocationAllFields;

  update_field_pointers_to_new_bucket(rank, bucketId, allFields, bucketCapacity, newAllocationAllFields);
}

void
DefaultFieldDataManager::reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                                unsigned bucketCapacity, const FieldVector & fields)
{
  for (const FieldBase * field : fields) {
    const FieldMetaData & fieldMetaData = field->get_meta_data_for_field()[bucketId];
    ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + bucketSize * fieldMetaData.m_bytesPerEntity,
                              (bucketCapacity - bucketSize) * fieldMetaData.m_bytesPerEntity);
  }
}

void resetFieldMetaDataPointers(const size_t bucket_index_begin, const size_t bucket_index_end,
                                FieldMetaDataVector& field_meta_data_vector, unsigned char* old_field_data,
                                unsigned char* new_field_data)
{
    for (size_t j=bucket_index_begin;j<bucket_index_end;j++)
    {
        FieldMetaData &field_meta_data = field_meta_data_vector[j];
        if ( field_meta_data.m_bytesPerEntity > 0 )
        {
            size_t sizeOfPreviousBuckets = field_meta_data.m_data - old_field_data;
            field_meta_data.m_data = new_field_data + sizeOfPreviousBuckets;
        }
    }
}

//////////////////////////////////////////////////////////////
ContiguousFieldDataManager::~ContiguousFieldDataManager()
{
    for (size_t i=0;i<m_field_raw_data.size();i++)
    {
        m_fieldDataAllocator->deallocate(m_field_raw_data[i], m_num_bytes_allocated_per_field[i]);
        m_field_raw_data[i] = nullptr;
        m_num_bytes_allocated_per_field[i] = 0;
        m_num_bytes_used_per_field[i] = 0;
    }
}

void ContiguousFieldDataManager::initialize_entity_field_data(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields)
{

}
void ContiguousFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
                                                            const std::vector<FieldBase *> & fields,
                                                            const PartVector& superset_parts,
                                                            unsigned size,
                                                            unsigned capacity)
{
  if (m_field_raw_data.empty()) {
    m_field_raw_data.resize(fields.size(), nullptr);
    m_num_entities_in_field_for_bucket.resize(fields.size());
    m_num_bytes_allocated_per_field.resize(fields.size(), 0);
    m_num_bytes_used_per_field.resize(fields.size(), 0);
  }

  for (size_t i = 0; i < fields.size(); i++) {
    if (fields[i]->entity_rank() == rank) {
      unsigned field_ordinal = fields[i]->mesh_meta_data_ordinal();
      const FieldLayoutData layout = getFieldLayoutData(*fields[i], rank, superset_parts);
      FieldMetaData field_meta_data;
      if (layout.numBytesPerEntity > 0) {
        field_meta_data.m_bytesPerEntity = layout.numBytesPerEntity;
        field_meta_data.m_firstDimension = layout.firstDimension;
        field_meta_data.m_data = m_field_raw_data[field_ordinal] + m_num_bytes_used_per_field[field_ordinal];
      }
      fields[i]->get_meta_data_for_field().push_back(field_meta_data);
      m_num_entities_in_field_for_bucket[field_ordinal].push_back(0);
    }
  }
}

void ContiguousFieldDataManager::clear_bucket_field_data(const EntityRank rm_rank, const unsigned rm_bucket_id, const std::vector<FieldBase*>  &allFields)
{
    for (size_t field_index = 0; field_index < allFields.size(); field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();

        if (field.entity_rank() == rm_rank && m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] > 0)
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[rm_bucket_id].m_bytesPerEntity;
            const unsigned char* endOfField = m_field_raw_data[field_ordinal] + m_num_bytes_used_per_field[field_ordinal];
            size_t sizeOfBucketToRemove = getFieldBucketSizeInBytes(field.get_meta_data_for_field(), rm_bucket_id, endOfField);

            if (numBytesPerEntity > 0)
            {
                unsigned char* new_field_data = m_field_raw_data[field_ordinal];

                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());

                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[rm_bucket_id];
                size_t sizeOfBucketsToTheLeft = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t rightHalfSize = m_num_bytes_used_per_field[field_ordinal] - sizeOfBucketsToTheLeft - sizeOfBucketToRemove;

                field_meta_data_for_modified_bucket.m_data = nullptr;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(rm_bucket_id+1, numBucketsOfRank, field_meta_data_vector,
                                           m_field_raw_data[field_ordinal], new_field_data-sizeOfBucketToRemove);

                ASAN_UNPOISON_MEMORY_REGION(new_field_data + sizeOfBucketsToTheLeft, sizeOfBucketToRemove + rightHalfSize);
                std::memmove(new_field_data + sizeOfBucketsToTheLeft,
                             m_field_raw_data[field_ordinal] + sizeOfBucketsToTheLeft + sizeOfBucketToRemove, rightHalfSize);

                m_num_bytes_used_per_field[field_ordinal] -= sizeOfBucketToRemove;
                m_field_raw_data[field_ordinal] = new_field_data;
                m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] = 0;
            }
        }
    }
}


void ContiguousFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
        const std::vector<FieldBase*>&  fields)
{
    if(fields.empty())
        return;

    this->clear_bucket_field_data(rank, bucket_id, fields);
    
    for(size_t field_index = 0; field_index < fields.size(); field_index++)
    {
        if(fields[field_index]->entity_rank() == rank)
        {
            unsigned field_ordinal = fields[field_index]->mesh_meta_data_ordinal();

            FieldMetaData& field_data = fields[field_index]->get_meta_data_for_field()[bucket_id];
            field_data.m_data = nullptr;
            field_data.m_bytesPerEntity = 0;
            field_data.m_firstDimension = 0;

            STK_ThrowRequireMsg(m_num_entities_in_field_for_bucket[field_ordinal][bucket_id] == 0, "Bucket not empty!");
        }
    }
}

void ContiguousFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                                           const std::vector<unsigned> & reorderedBucketIds)
{
    for (size_t field_index = 0; field_index < fields.size(); ++field_index)
    {
        if (fields[field_index]->entity_rank() == rank)
        {
            FieldMetaDataVector &oldMetaData = const_cast<FieldMetaDataVector &>(fields[field_index]->get_meta_data_for_field());
            unsigned field_ordinal = fields[field_index]->mesh_meta_data_ordinal();
            const size_t newFieldSize = m_num_bytes_used_per_field[field_ordinal] + m_extra_capacity;
            unsigned char* new_field_data = m_fieldDataAllocator->allocate(newFieldSize);

            FieldMetaData fieldMeta;
            FieldMetaDataVector newMetaData(reorderedBucketIds.size(), fieldMeta);
            std::vector<size_t> newNumEntitiesPerBucket(reorderedBucketIds.size(), 0);
            unsigned new_offset = 0;
            for (unsigned bucket_index = 0, bucket_end = reorderedBucketIds.size(); bucket_index < bucket_end; ++bucket_index)
            {
                unsigned oldBucketIndex = reorderedBucketIds[bucket_index];
                const unsigned char* bucket_start_ptr = oldMetaData[oldBucketIndex].m_data;

                if (oldMetaData[oldBucketIndex].m_bytesPerEntity > 0)
                {
                    newNumEntitiesPerBucket[bucket_index] = m_num_entities_in_field_for_bucket[field_ordinal][oldBucketIndex];

                    const unsigned char* end_of_field = m_field_raw_data[field_ordinal]+m_num_bytes_used_per_field[field_ordinal];
                    unsigned bucket_size = getFieldBucketSizeInBytes(oldMetaData, oldBucketIndex, end_of_field);
                    ASAN_UNPOISON_MEMORY_REGION(bucket_start_ptr, bucket_size);
                    ASAN_UNPOISON_MEMORY_REGION(new_field_data+new_offset, bucket_size);
                    newMetaData[bucket_index].m_data = &new_field_data[new_offset];
                    newMetaData[bucket_index].m_bytesPerEntity = oldMetaData[oldBucketIndex].m_bytesPerEntity;
                    newMetaData[bucket_index].m_firstDimension = oldMetaData[oldBucketIndex].m_firstDimension;

                    std::memcpy(new_field_data+new_offset, bucket_start_ptr, bucket_size);
                    new_offset += bucket_size;
                }
            }

            m_fieldDataAllocator->deallocate(m_field_raw_data[field_ordinal], m_num_bytes_allocated_per_field[field_ordinal]);
            STK_ThrowRequire(new_offset == m_num_bytes_used_per_field[field_ordinal]);
            m_num_bytes_allocated_per_field[field_ordinal] = newFieldSize;
            m_field_raw_data[field_ordinal] = new_field_data;
            m_num_entities_in_field_for_bucket[field_ordinal].swap(newNumEntitiesPerBucket);
            newMetaData.swap(oldMetaData);
        }
    }
}

void ContiguousFieldDataManager::remove_field_data_for_entity(EntityRank rm_rank, unsigned rm_bucket_id,
                                                              unsigned rm_bucket_ord, const std::vector<FieldBase *> &allFields)
{
    for(size_t field_index = 0; field_index < allFields.size(); field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();
        if(field.entity_rank() == rm_rank)
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[rm_bucket_id].m_bytesPerEntity;
            if(numBytesPerEntity > 0)
            {
                unsigned char* new_field_data = m_field_raw_data[field_ordinal];
                const unsigned char* endOfField = m_field_raw_data[field_ordinal]+m_num_bytes_used_per_field[field_ordinal];

                const size_t currentBucketStorageUsed = numBytesPerEntity * m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id];
                const size_t currentBucketAllocation = getFieldBucketSizeInBytes(field.get_meta_data_for_field(), rm_bucket_id, endOfField);
                const size_t newBucketStorageUsed = currentBucketStorageUsed - numBytesPerEntity;
                const size_t newBucketAllocation = stk::adjust_up_to_alignment_boundary(newBucketStorageUsed, alignment_increment_bytes);
                const size_t allocationToRemove = currentBucketAllocation - newBucketAllocation;

                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());
                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[rm_bucket_id];
                size_t sizeOfBucketsToTheLeft = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t rightHalfSize = m_num_bytes_used_per_field[field_ordinal] - sizeOfBucketsToTheLeft - currentBucketAllocation;

                field_meta_data_for_modified_bucket.m_data = new_field_data + sizeOfBucketsToTheLeft;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(rm_bucket_id+1, numBucketsOfRank, field_meta_data_vector,
                                           m_field_raw_data[field_ordinal], new_field_data-allocationToRemove);

                std::memmove(new_field_data + sizeOfBucketsToTheLeft + newBucketAllocation,
                             m_field_raw_data[field_ordinal] + sizeOfBucketsToTheLeft + currentBucketAllocation, rightHalfSize);
                m_num_bytes_used_per_field[field_ordinal] -= allocationToRemove;
                m_field_raw_data[field_ordinal] = new_field_data;
                m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] -= 1;
            }
        }
    }
}

void ContiguousFieldDataManager::allocate_field_data(EntityRank rank,
                                                     const std::vector<Bucket*>& buckets,
                                                     const std::vector<FieldBase *> & fields)
{
  m_field_raw_data.resize(fields.size(), nullptr);
  m_num_entities_in_field_for_bucket.resize(fields.size());
  for (size_t i=0; i<fields.size(); i++) {
    if (fields[i]->entity_rank() == rank) {
      unsigned field_ordinal = fields[i]->mesh_meta_data_ordinal();
      m_num_entities_in_field_for_bucket[field_ordinal].resize(buckets.size(),0);
    }
  }

  m_num_bytes_allocated_per_field.resize(fields.size(), m_extra_capacity);
  m_num_bytes_used_per_field.resize(fields.size(), 0);
  for (size_t fieldIndex = 0; fieldIndex != fields.size(); fieldIndex++) {
    const FieldBase& field = *fields[fieldIndex];
    if (field.entity_rank() == rank) {
      unsigned field_ordinal = fields[fieldIndex]->mesh_meta_data_ordinal();
      for (size_t i = 0; i < buckets.size(); ++i) {
        const PartVector& superset_parts = buckets[i]->supersets();
        const FieldLayoutData layout = getFieldLayoutData(field, rank, superset_parts);
        if (layout.numBytesPerEntity > 0) {
          m_num_entities_in_field_for_bucket[field_ordinal][i] = buckets[i]->size();
          m_num_bytes_used_per_field[field_ordinal] += stk::adjust_up_to_alignment_boundary(layout.numBytesPerEntity * buckets[i]->size(),
                                                                                            alignment_increment_bytes);
        }
        FieldMetaData fieldMetaData;
        fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
        fieldMetaData.m_firstDimension = layout.firstDimension;
        const_cast<FieldBase&>(field).get_meta_data_for_field().push_back(fieldMetaData);
      }

      m_num_bytes_allocated_per_field[field_ordinal] += m_num_bytes_used_per_field[field_ordinal];
      m_field_raw_data[field_ordinal] = m_fieldDataAllocator->allocate(m_num_bytes_allocated_per_field[field_ordinal]);

      size_t offset = 0;
      for(size_t i = 0; i < buckets.size(); ++i)
      {
        const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
        FieldMetaData& field_meta_data = const_cast<FieldMetaData&>(field.get_meta_data_for_field()[i]);
        updateFieldPointer(field_meta_data, buckets[i]->size(), offset, m_field_raw_data[field_ordinal], alignment_increment_bytes);
        initializeField(field_meta_data, init_val, buckets[i]->size(), buckets[i]->size());
      }
    }
  }
}

void ContiguousFieldDataManager::allocate_new_field_meta_data(const EntityRank rank, const std::vector<Bucket*> & buckets, const std::vector<FieldBase*>& allFields)
{
    for (FieldBase* field : allFields) {
        for (stk::mesh::Bucket * bucket : buckets) {
            if (bucket->bucket_id() >= field->get_meta_data_for_field().size()) {
                FieldMetaData fieldMetaData;
                field->get_meta_data_for_field().push_back(fieldMetaData);
            }
        }
    }
}

std::vector<size_t> ContiguousFieldDataManager::get_field_bucket_offsets(const std::vector<Bucket*> & buckets,
                                                                         FieldBase & currentField) const
{
    std::vector<size_t> offsetForBucket;
    size_t currentBucketOffset = 0;
    for (unsigned int i = 0; i < buckets.size(); ++i) {
        const size_t bytesPerEntity = currentField.get_meta_data_for_field()[i].m_bytesPerEntity;
        size_t bucketDataSizeThisField = stk::adjust_up_to_alignment_boundary(bytesPerEntity*buckets[i]->size(), alignment_increment_bytes);

        offsetForBucket.push_back(currentBucketOffset);
        currentBucketOffset += bucketDataSizeThisField;
    }
    offsetForBucket.push_back(currentBucketOffset);
    return offsetForBucket;
}

void ContiguousFieldDataManager::copy_bucket_data_from_old_to_new_field(const std::vector<size_t>& oldOffsetForBucket,
                                                                        const std::vector<size_t>& newOffsetForBucket,
                                                                        const unsigned char* oldAllocationAllBuckets,
                                                                        unsigned char* newAllocationAllBuckets)
{
    for (size_t bucketIndex = 0; bucketIndex < oldOffsetForBucket.size()-1; ++bucketIndex) {
        size_t oldBucketAllocationSize = oldOffsetForBucket[bucketIndex + 1] - oldOffsetForBucket[bucketIndex];
        ASAN_UNPOISON_MEMORY_REGION(oldAllocationAllBuckets + oldOffsetForBucket[bucketIndex], oldBucketAllocationSize);
        ASAN_UNPOISON_MEMORY_REGION(newAllocationAllBuckets + newOffsetForBucket[bucketIndex], oldBucketAllocationSize);
        std::memcpy(newAllocationAllBuckets + newOffsetForBucket[bucketIndex],
                    oldAllocationAllBuckets + oldOffsetForBucket[bucketIndex], oldBucketAllocationSize);
    }
}

void ContiguousFieldDataManager::update_bucket_storage_for_field(EntityRank rank,
                                                                 const std::vector<Bucket*>& buckets,
                                                                 FieldBase& currentField)
{
  const unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
  for (size_t i = 0; i < buckets.size(); ++i) {
    const PartVector& supersetParts = buckets[i]->supersets();
    const FieldLayoutData layout = getFieldLayoutData(currentField, rank, supersetParts);
    if (layout.numBytesPerEntity > 0) {
      m_num_entities_in_field_for_bucket[fieldOrdinal][i] = buckets[i]->size();
      m_num_bytes_used_per_field[fieldOrdinal] += stk::adjust_up_to_alignment_boundary(layout.numBytesPerEntity * buckets[i]->size(),
                                                                                       alignment_increment_bytes);
    }
    FieldMetaData& fieldMetaData = currentField.get_meta_data_for_field()[i];
    fieldMetaData.m_bytesPerEntity = layout.numBytesPerEntity;
    fieldMetaData.m_firstDimension = layout.firstDimension;
  }
}

void ContiguousFieldDataManager::update_bucket_pointers_to_new_field(const std::vector<Bucket*>& buckets, FieldBase& currentField)
{
    const unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
    size_t offset = 0;
    for (size_t i = 0; i < buckets.size(); ++i) {
        FieldMetaData& fieldMetaData = currentField.get_meta_data_for_field()[i];
        updateFieldPointer(fieldMetaData, buckets[i]->size(), offset, m_field_raw_data[fieldOrdinal], alignment_increment_bytes);
    }
}

void ContiguousFieldDataManager::initialize_new_bucket_values(const std::vector<Bucket*>& buckets,
                                                              const std::vector<size_t> & oldOffsetForBucket,
                                                              const std::vector<size_t> & newOffsetForBucket,
                                                              FieldBase& currentField)
{
    for (size_t i = 0; i < buckets.size(); ++i) {
        const size_t oldSize = oldOffsetForBucket[i+1] - oldOffsetForBucket[i];
        const size_t newSize = newOffsetForBucket[i+1] - newOffsetForBucket[i];
        const bool isNewBucketForField = (oldSize == 0u) && (newSize > 0u);
        if (isNewBucketForField) {
            const unsigned char* initVal = reinterpret_cast<const unsigned char*>(currentField.get_initial_value());
            FieldMetaData& fieldMetaData = currentField.get_meta_data_for_field()[i];
            initializeField(fieldMetaData, initVal, buckets[i]->size(), buckets[i]->size());
        }
    }
}

void ContiguousFieldDataManager::reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField, const std::vector<FieldBase *> & allFields)
{
    m_field_raw_data.resize(allFields.size(), nullptr);
    m_num_bytes_allocated_per_field.resize(allFields.size(), 0);
    m_num_bytes_used_per_field.resize(allFields.size(), 0);
    m_num_entities_in_field_for_bucket.resize(allFields.size(), std::vector<size_t>(buckets.size(), 0));

    if (currentField.entity_rank() == rank)
    {
        allocate_new_field_meta_data(rank, buckets, allFields);

        unsigned fieldOrdinal = currentField.mesh_meta_data_ordinal();
        std::vector<size_t> oldOffsetForBucket = get_field_bucket_offsets(buckets, currentField);
        const size_t oldFieldBucketAllocationSize = m_num_bytes_allocated_per_field[fieldOrdinal];

        update_bucket_storage_for_field(rank, buckets, currentField);

        std::vector<size_t> newOffsetForBucket = get_field_bucket_offsets(buckets, currentField);
        const size_t newFieldBucketAllocationSize = newOffsetForBucket.back() + m_extra_capacity;

        if (newFieldBucketAllocationSize > oldFieldBucketAllocationSize) {
            unsigned char* newAllocationAllBuckets = m_fieldDataAllocator->allocate(newFieldBucketAllocationSize);
            const unsigned char* oldAllocationAllBuckets = m_field_raw_data[fieldOrdinal];

            copy_bucket_data_from_old_to_new_field(oldOffsetForBucket, newOffsetForBucket, oldAllocationAllBuckets, newAllocationAllBuckets);

            m_fieldDataAllocator->deallocate(m_field_raw_data[fieldOrdinal], oldFieldBucketAllocationSize);
            m_field_raw_data[fieldOrdinal] = newAllocationAllBuckets;
            m_num_bytes_allocated_per_field[fieldOrdinal] = newFieldBucketAllocationSize;

            update_bucket_pointers_to_new_field(buckets, currentField);
            initialize_new_bucket_values(buckets, oldOffsetForBucket, newOffsetForBucket, currentField);
        }
    }
}

size_t getFieldBucketSizeInBytes(const FieldMetaDataVector& field_meta_data_vector, const unsigned bucket_id, const unsigned char *end_of_field)
{
    size_t sizeFieldThisBucketInBytes = end_of_field - field_meta_data_vector[bucket_id].m_data;

    for (unsigned nextBucket=bucket_id+1;nextBucket<field_meta_data_vector.size();nextBucket++)
    {
        if ( field_meta_data_vector[nextBucket].m_bytesPerEntity > 0 )
        {
            sizeFieldThisBucketInBytes = field_meta_data_vector[nextBucket].m_data - field_meta_data_vector[bucket_id].m_data;
            break;
        }
    }

    return sizeFieldThisBucketInBytes;
}

void ContiguousFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> & allFields,
                                                           EntityRank dst_rank, unsigned dst_bucket_id,
                                                           unsigned dst_bucket_ord)
{
    for (size_t field_index=0;field_index<allFields.size();field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();
        if ( field.entity_rank() == dst_rank )
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[dst_bucket_id].m_bytesPerEntity;
            if ( numBytesPerEntity > 0 )
            {
                const unsigned char* endOfField = m_field_raw_data[field_ordinal]+m_num_bytes_used_per_field[field_ordinal];
                const size_t currentBucketStorageUsed = numBytesPerEntity * m_num_entities_in_field_for_bucket[field_ordinal][dst_bucket_id];
                const size_t currentBucketAllocation = getFieldBucketSizeInBytes(field.get_meta_data_for_field(), dst_bucket_id, endOfField);
                const size_t newBucketStorageUsed = currentBucketStorageUsed + numBytesPerEntity;
                const size_t newBucketAllocation = stk::adjust_up_to_alignment_boundary(newBucketStorageUsed,
                                                                                        alignment_increment_bytes);
                const size_t extraAllocationNeeded = newBucketAllocation - currentBucketAllocation;
                const size_t newFieldSizeNeeded = m_num_bytes_used_per_field[field_ordinal] + extraAllocationNeeded + m_extra_capacity;
                const size_t newFieldSize = std::max(newFieldSizeNeeded, m_num_bytes_allocated_per_field[field_ordinal]);
                bool requiresNewAllocation = newFieldSize > m_num_bytes_allocated_per_field[field_ordinal];

                unsigned char* new_field_data = m_field_raw_data[field_ordinal];
                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());

                if (requiresNewAllocation)
                {
                    new_field_data = m_fieldDataAllocator->allocate(newFieldSize);
                    resetFieldMetaDataPointers(0, dst_bucket_id, field_meta_data_vector, m_field_raw_data[field_ordinal], new_field_data);
                }

                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[dst_bucket_id];
                size_t sizeOfBucketsToTheLeft = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t leftHalfSize = sizeOfBucketsToTheLeft + currentBucketStorageUsed;
                size_t rightHalfSize = m_num_bytes_used_per_field[field_ordinal] - sizeOfBucketsToTheLeft - currentBucketAllocation;
                field_meta_data_for_modified_bucket.m_data = new_field_data + sizeOfBucketsToTheLeft;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(dst_bucket_id + 1, numBucketsOfRank, field_meta_data_vector, m_field_raw_data[field_ordinal],
                                           new_field_data + extraAllocationNeeded);

                if (requiresNewAllocation)
                {
                    // We lose some safety here because, ideally, we would only unpoison the used portions of the
                    // allocation and skip over the SIMD padding, but that requires querying each Bucket's size.
                    // During mesh construction, querying the Buckets would be done too early which would put the
                    // BucketRegistrar in a bad state.
                    ASAN_UNPOISON_MEMORY_REGION(m_field_raw_data[field_ordinal],
                                                leftHalfSize + currentBucketAllocation + rightHalfSize);
                    ASAN_UNPOISON_MEMORY_REGION(new_field_data,
                                                leftHalfSize + newBucketAllocation + rightHalfSize);
                    std::memcpy(new_field_data, m_field_raw_data[field_ordinal], leftHalfSize);
                    std::memcpy(new_field_data + sizeOfBucketsToTheLeft + newBucketAllocation,
                                m_field_raw_data[field_ordinal] + sizeOfBucketsToTheLeft + currentBucketAllocation, rightHalfSize);
                    m_fieldDataAllocator->deallocate(m_field_raw_data[field_ordinal], m_num_bytes_allocated_per_field[field_ordinal]);
                    m_num_bytes_allocated_per_field[field_ordinal] = newFieldSize;
                }
                else
                {
                    ASAN_UNPOISON_MEMORY_REGION(m_field_raw_data[field_ordinal] + sizeOfBucketsToTheLeft,
                                                newBucketAllocation + rightHalfSize);
                    std::memmove(new_field_data + sizeOfBucketsToTheLeft + newBucketAllocation,
                                 m_field_raw_data[field_ordinal] + sizeOfBucketsToTheLeft + currentBucketAllocation, rightHalfSize);
                }

                setInitialValue(new_field_data+leftHalfSize, field, numBytesPerEntity);

                m_num_bytes_used_per_field[field_ordinal] += extraAllocationNeeded;
                m_field_raw_data[field_ordinal] = new_field_data;
                STK_ThrowRequire(dst_bucket_ord == m_num_entities_in_field_for_bucket[field_ordinal][dst_bucket_id]);
                m_num_entities_in_field_for_bucket[field_ordinal][dst_bucket_id] += 1;
            }
        }
    }
}

void ContiguousFieldDataManager::swap_fields(const int field1, const int field2)
{
    std::swap(m_field_raw_data[field1], m_field_raw_data[field2]);
    STK_ThrowRequire(m_num_bytes_allocated_per_field[field1] == m_num_bytes_allocated_per_field[field2]);
    STK_ThrowRequire(m_num_bytes_used_per_field[field1] == m_num_bytes_used_per_field[field2]);
    STK_ThrowRequire(m_num_entities_in_field_for_bucket[field1].size() == m_num_entities_in_field_for_bucket[field2].size());
}

void
ContiguousFieldDataManager::reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                                   unsigned bucketCapacity, const FieldVector & fields)
{
  for (const FieldBase * field : fields) {
    const FieldMetaData & fieldMetaData = field->get_meta_data_for_field()[bucketId];
    const unsigned bucketStorageUsed = bucketSize * fieldMetaData.m_bytesPerEntity;
    const unsigned bucketStorageAllocated = stk::adjust_up_to_alignment_boundary(bucketStorageUsed,
                                                                                 alignment_increment_bytes);
    ASAN_POISON_MEMORY_REGION(fieldMetaData.m_data + bucketStorageUsed, bucketStorageAllocated - bucketStorageUsed);
  }
}

}
}
