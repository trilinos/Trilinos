// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/FieldBase.hpp>  // for FieldMetaData, etc
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/FindRestriction.hpp>
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Types.hpp"      // for EntityRank, PartVector
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequire, etc

namespace stk {
namespace mesh {

size_t getFieldBucketSizeInBytes(FieldMetaDataVector& field_meta_data_vector, const unsigned bucket_id, const unsigned char *end_of_field);

int getNumBytesForField(const FieldBase &field, const EntityRank rank, const PartVector& superset_parts)
{
    int num_bytes_per_entity = 0;

    const FieldBase::Restriction & restriction =
            find_and_check_restriction(field, rank, superset_parts);

    if(restriction.num_scalars_per_entity() > 0)
    { // Exists

        const unsigned type_stride = field.data_traits().stride_of;
        const unsigned field_array_rank = field.field_array_rank();

        num_bytes_per_entity = type_stride *
                (field_array_rank ? restriction.num_scalars_per_entity() : 1);
    }
    return num_bytes_per_entity;
}

void initializeField(FieldMetaData& field_meta_data, const unsigned char* init_val, const size_t capacity, size_t &current_field_offset, unsigned char* all_data)
{
    if (field_meta_data.m_bytes_per_entity > 0)
    {
        field_meta_data.m_data = all_data + current_field_offset;
        current_field_offset += field_meta_data.m_bytes_per_entity * capacity;

        // initialize field data
        if (init_val != NULL)
        {
            for (size_t j = 0; j < capacity; ++j)
            {
                std::memcpy( field_meta_data.m_data + j * field_meta_data.m_bytes_per_entity, init_val, field_meta_data.m_bytes_per_entity );
            }
        }
        else
        {
            std::memset( field_meta_data.m_data, 0, capacity * field_meta_data.m_bytes_per_entity );
        }
    }
}

void setInitialValue(unsigned char* data_location, const FieldBase& field, const int numBytesPerEntity)
{
    const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
    if (init_val != NULL)
    {
        std::memcpy( data_location, init_val, numBytesPerEntity );
    }
    else
    {
        std::memset( data_location, 0, numBytesPerEntity );
    }
}

void DefaultFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
        const std::vector< FieldBase * > & fields, const PartVector& superset_parts, const size_t capacity)
{
    if ( m_num_bytes_allocated_per_field.empty() )
    {
        m_num_bytes_allocated_per_field.resize(fields.size(), 0);
    }

    size_t num_fields = fields.size();
    // Sizing loop
    size_t total_field_data_size = 0;

    for(size_t i = 0; i < num_fields; ++i)
    {
        FieldMetaData field_meta_data = {NULL, 0};

        const FieldBase & field = *fields[i];
        if(static_cast<unsigned>(field.entity_rank()) == rank)
        {
            size_t num_bytes_per_entity = getNumBytesForField(field, rank, superset_parts);
            if(num_bytes_per_entity > 0)
            {
                field_meta_data.m_bytes_per_entity = num_bytes_per_entity;
                total_field_data_size += num_bytes_per_entity * capacity;
                m_num_bytes_allocated_per_field[i] += num_bytes_per_entity * capacity;
            }
            fields[i]->get_meta_data_for_field().push_back(field_meta_data);
        }
    }

    // Allocate all field data for this bucket
    if(total_field_data_size > 0)
    {
        unsigned char* all_data = field_data_allocator().allocate(total_field_data_size);
        m_field_raw_data[rank].push_back(all_data);

        // Set data ptrs in field meta datas
        size_t current_field_offset = 0;
        for(size_t i = 0; i < fields.size(); ++i)
        {
            const FieldBase & field = *fields[i];
            if(static_cast<unsigned>(field.entity_rank()) == rank)
            {
                const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
                FieldMetaData& field_meta_data = const_cast<FieldMetaData&>(field.get_meta_data_for_field().back());
                initializeField(field_meta_data, init_val, capacity, current_field_offset, all_data);
            }
        }
    }
    else
    {
        m_field_raw_data[rank].push_back(NULL);
    }
}

void DefaultFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
        const std::vector<FieldBase*>  fields)
{
    if(fields.empty())
        return;

    if(m_field_raw_data[rank][bucket_id] != NULL)
    {
        size_t bytes_to_delete = 0;
        for(unsigned int i = 0; i < fields.size(); ++i)
        {
            if(fields[i] == NULL || static_cast<unsigned>(fields[i]->entity_rank()) != rank)
                continue;
            FieldMetaData& field_data = fields[i]->get_meta_data_for_field()[bucket_id];
            if(field_data.m_data != NULL)
            {
                bytes_to_delete += field_data.m_bytes_per_entity * capacity;
                field_data.m_bytes_per_entity = 0;
                field_data.m_data = NULL;
            }
        }
        field_data_allocator().deallocate(m_field_raw_data[rank][bucket_id], bytes_to_delete);
        m_field_raw_data[rank][bucket_id] = NULL;
    }
}

void DefaultFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds)
{
    std::vector<unsigned char*> field_raw_data(reorderedBucketIds.size());
    for(unsigned m = 0, e = reorderedBucketIds.size(); m < e; ++m)
    {
        field_raw_data[m] = m_field_raw_data[rank][reorderedBucketIds[m]];
    }
    m_field_raw_data[rank].swap(field_raw_data);

    for(size_t i = 0; i < fields.size(); ++i)
    {
        if(static_cast<unsigned>(fields[i]->entity_rank()) == rank)
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
        allocate_bucket_field_data(rank, fields, superset_parts, buckets[i]->capacity());
    }
}

void DefaultFieldDataManager::remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields)
{

}
void DefaultFieldDataManager::reinitialize_removed_entity_field_data(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields)
{
    // bucket of bucket_id shrinks by one
    for(size_t i = 0; i < fields.size(); ++i)
    {
        const FieldBase & field = *fields[i];
        if(static_cast<unsigned>(field.entity_rank()) == rank)
        {
            FieldMetaData field_meta_data = fields[i]->get_meta_data_for_field()[bucket_id];
            const int num_bytes_per_entity = field_meta_data.m_bytes_per_entity;

            if(num_bytes_per_entity > 0)
            {
                setInitialValue(field_meta_data.m_data + bucket_ord * num_bytes_per_entity, field, num_bytes_per_entity);
            }
        }
    }
}

void DefaultFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord )
{

}

void resetFieldMetaDataPointers(const size_t bucket_index_begin, const size_t bucket_index_end,
        FieldMetaDataVector& field_meta_data_vector, unsigned char* old_field_data, unsigned char* new_field_data)
{
    for (size_t j=bucket_index_begin;j<bucket_index_end;j++)
    {
        FieldMetaData &field_meta_data = field_meta_data_vector[j];
        if ( field_meta_data.m_bytes_per_entity > 0 )
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
        field_data_allocator().deallocate(m_field_raw_data[i], m_num_bytes_allocated_per_field[i]);
        m_field_raw_data[i] = NULL;
        m_num_bytes_allocated_per_field[i] = 0;
        m_num_bytes_used_per_field[i] = 0;
    }
}

void ContiguousFieldDataManager::reinitialize_removed_entity_field_data(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields)
{

}
void ContiguousFieldDataManager::allocate_bucket_field_data(const EntityRank rank,
        const std::vector< FieldBase * > & fields, const PartVector& superset_parts, const size_t capacity)
{
    if(m_field_raw_data.empty())
    {
        m_field_raw_data.resize(fields.size(), NULL);
        m_num_entities_in_field_for_bucket.resize(fields.size());
        m_num_bytes_allocated_per_field.resize(fields.size(), 0);
        m_num_bytes_used_per_field.resize(fields.size(), 0);
    }

    for(size_t i = 0; i < fields.size(); i++)
    {
        if(fields[i]->entity_rank() == rank)
        {
            unsigned field_ordinal = fields[i]->mesh_meta_data_ordinal();
            size_t num_bytes_per_entity = getNumBytesForField(*fields[i], rank, superset_parts);
            FieldMetaData field_meta_data = {NULL, 0};
            if(num_bytes_per_entity > 0)
            {
                field_meta_data.m_bytes_per_entity = num_bytes_per_entity;
                field_meta_data.m_data = m_field_raw_data[field_ordinal] + m_num_bytes_used_per_field[field_ordinal];
            }
            fields[i]->get_meta_data_for_field().push_back(field_meta_data);
            m_num_entities_in_field_for_bucket[field_ordinal].push_back(0);
        }
    }
}

void ContiguousFieldDataManager::clear_bucket_field_data(const EntityRank rm_rank, const unsigned rm_bucket_id, const std::vector<FieldBase*>  &allFields)
{
    for(size_t field_index = 0; field_index < allFields.size(); field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();

        if(field.entity_rank() == rm_rank && m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] > 0)
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[rm_bucket_id].m_bytes_per_entity;
            size_t sizeFieldForRemovedBucket = numBytesPerEntity*m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id];
            if(numBytesPerEntity > 0)
            {
                unsigned char* new_field_data = m_field_raw_data[field_ordinal];

                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());

                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[rm_bucket_id];
                size_t sizeOfPreviousBuckets = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t leftHalf = sizeOfPreviousBuckets;
                size_t rightHalf = m_num_bytes_used_per_field[field_ordinal] - leftHalf - sizeFieldForRemovedBucket;

                field_meta_data_for_modified_bucket.m_data = NULL;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(rm_bucket_id+1, numBucketsOfRank, field_meta_data_vector, m_field_raw_data[field_ordinal], new_field_data-sizeFieldForRemovedBucket);

                std::memmove(new_field_data + leftHalf, m_field_raw_data[field_ordinal] + leftHalf + sizeFieldForRemovedBucket, rightHalf);

                m_num_bytes_used_per_field[field_ordinal] -= sizeFieldForRemovedBucket;
                m_field_raw_data[field_ordinal] = new_field_data;
                m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] = 0;
            }
        }
    }
}


void ContiguousFieldDataManager::deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
        const std::vector<FieldBase*>  fields)
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
            field_data.m_bytes_per_entity = 0;
            field_data.m_data = NULL;

            ThrowRequireMsg(m_num_entities_in_field_for_bucket[field_ordinal][bucket_id] == 0, "Bucket not empty!");
        }
    }
}

void ContiguousFieldDataManager::reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds)
{
    for(size_t field_index = 0; field_index < fields.size(); ++field_index)
    {
        if(fields[field_index]->entity_rank() == rank)
        {
            FieldMetaDataVector &oldMetaData = const_cast<FieldMetaDataVector &>(fields[field_index]->get_meta_data_for_field());
            size_t newFieldSizeInBytes = 0;
            unsigned field_ordinal = fields[field_index]->mesh_meta_data_ordinal();
            for(size_t i=0; i<reorderedBucketIds.size(); i++)
            {
                newFieldSizeInBytes += m_num_entities_in_field_for_bucket[field_ordinal][reorderedBucketIds[i]]*oldMetaData[reorderedBucketIds[i]].m_bytes_per_entity;
            }

            unsigned char* new_field_data = field_data_allocator().allocate(newFieldSizeInBytes + m_extra_capacity);

            FieldMetaData fieldMeta = {NULL, 0};
            FieldMetaDataVector newMetaData(reorderedBucketIds.size(), fieldMeta);
            std::vector<size_t> size_bucket_in_field(reorderedBucketIds.size(), 0);
            unsigned new_offset = 0;
            for(unsigned bucket_index = 0, bucket_end = reorderedBucketIds.size(); bucket_index < bucket_end; ++bucket_index)
            {
                unsigned oldBucketIndex = reorderedBucketIds[bucket_index];
                const unsigned char* bucket_start_ptr = oldMetaData[oldBucketIndex].m_data;

                if (oldMetaData[oldBucketIndex].m_bytes_per_entity > 0)
                {
                    const unsigned char* end_of_field = m_field_raw_data[field_ordinal]+m_num_bytes_used_per_field[field_ordinal];
                    unsigned bucket_size = getFieldBucketSizeInBytes(oldMetaData, oldBucketIndex, end_of_field);
                    newMetaData[bucket_index].m_data = &new_field_data[new_offset];

                    int bytes_per_entity = oldMetaData[oldBucketIndex].m_bytes_per_entity;

                    ThrowRequire(m_num_entities_in_field_for_bucket[field_ordinal][oldBucketIndex] == bucket_size/bytes_per_entity);
                    size_bucket_in_field[bucket_index] = bucket_size/bytes_per_entity;
                    newMetaData[bucket_index].m_bytes_per_entity = bytes_per_entity;

                    std::memcpy(new_field_data+new_offset, bucket_start_ptr, bucket_size);
                    new_offset += bucket_size;
                }
            }

            field_data_allocator().deallocate(m_field_raw_data[field_ordinal], m_num_bytes_allocated_per_field[field_ordinal]);
            ThrowRequire(new_offset == newFieldSizeInBytes);
            m_num_bytes_used_per_field[field_ordinal] = newFieldSizeInBytes;
            m_num_bytes_allocated_per_field[field_ordinal] = newFieldSizeInBytes + m_extra_capacity;
            m_field_raw_data[field_ordinal] = new_field_data;
            m_num_entities_in_field_for_bucket[field_ordinal].swap(size_bucket_in_field);
            newMetaData.swap(oldMetaData);
        }
    }
}

void ContiguousFieldDataManager::remove_field_data_for_entity(EntityRank rm_rank, unsigned rm_bucket_id, Bucket::size_type rm_bucket_ord, const std::vector<FieldBase *> &allFields)
{
    for(size_t field_index = 0; field_index < allFields.size(); field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();
        if(field.entity_rank() == rm_rank)
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[rm_bucket_id].m_bytes_per_entity;
            if(numBytesPerEntity > 0)
            {
                unsigned char* new_field_data = m_field_raw_data[field_ordinal];

                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());
                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[rm_bucket_id];
                size_t sizeOfPreviousBuckets = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t leftHalf = sizeOfPreviousBuckets + rm_bucket_ord * numBytesPerEntity;
                size_t rightHalf = m_num_bytes_used_per_field[field_ordinal] - leftHalf - numBytesPerEntity;

                field_meta_data_for_modified_bucket.m_data = new_field_data + sizeOfPreviousBuckets;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(rm_bucket_id+1, numBucketsOfRank, field_meta_data_vector, m_field_raw_data[field_ordinal], new_field_data-numBytesPerEntity);

                std::memmove(new_field_data + leftHalf, m_field_raw_data[field_ordinal] + leftHalf + numBytesPerEntity, rightHalf);
                m_num_bytes_used_per_field[field_ordinal] -= numBytesPerEntity;
                m_field_raw_data[field_ordinal] = new_field_data;
                m_num_entities_in_field_for_bucket[field_ordinal][rm_bucket_id] -= 1;
            }
        }
    }
}

void ContiguousFieldDataManager::allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & fields)
{
    m_field_raw_data.resize(fields.size(), NULL);
    m_num_entities_in_field_for_bucket.resize(fields.size());
    for (size_t i=0;i<fields.size();i++)
    {
        if(fields[i]->entity_rank() == rank)
        {
            unsigned field_ordinal = fields[i]->mesh_meta_data_ordinal();
            m_num_entities_in_field_for_bucket[field_ordinal].resize(buckets.size(),0);
        }
    }

    m_num_bytes_allocated_per_field.resize(fields.size(), m_extra_capacity);
    m_num_bytes_used_per_field.resize(fields.size(), 0);
    for(size_t fieldIndex = 0; fieldIndex != fields.size(); fieldIndex++)
    {
        const FieldBase& field = *fields[fieldIndex];
        if(field.entity_rank() == rank)
        {
            unsigned field_ordinal = fields[fieldIndex]->mesh_meta_data_ordinal();
            for(size_t i = 0; i < buckets.size(); ++i)
            {
                const PartVector& superset_parts = buckets[i]->supersets();
                int numBytesPerEntity = getNumBytesForField(field, rank, superset_parts);
                if ( numBytesPerEntity > 0 )
                {
                    m_num_entities_in_field_for_bucket[field_ordinal][i] = buckets[i]->size();
                    m_num_bytes_used_per_field[field_ordinal] += numBytesPerEntity * buckets[i]->size();
                }
                FieldMetaData fieldMetaData = {NULL, numBytesPerEntity};
                const_cast<FieldBase&>(field).get_meta_data_for_field().push_back(fieldMetaData);
            }

            m_num_bytes_allocated_per_field[field_ordinal] += m_num_bytes_used_per_field[field_ordinal];
            m_field_raw_data[field_ordinal] = field_data_allocator().allocate(m_num_bytes_allocated_per_field[field_ordinal]);

            size_t offset = 0;
            for(size_t i = 0; i < buckets.size(); ++i)
            {
                const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
                FieldMetaData& field_meta_data = const_cast<FieldMetaData&>(field.get_meta_data_for_field()[i]);
                initializeField(field_meta_data, init_val, buckets[i]->size(), offset, m_field_raw_data[field_ordinal]);
            }
        }
    }
}

size_t getFieldBucketSizeInBytes(FieldMetaDataVector& field_meta_data_vector, const unsigned bucket_id, const unsigned char *end_of_field)
{
    size_t sizeFieldThisBucketInBytes = end_of_field - field_meta_data_vector[bucket_id].m_data;

    for (unsigned nextBucket=bucket_id+1;nextBucket<field_meta_data_vector.size();nextBucket++)
    {
        if ( field_meta_data_vector[nextBucket].m_bytes_per_entity > 0 )
        {
            sizeFieldThisBucketInBytes = field_meta_data_vector[nextBucket].m_data - field_meta_data_vector[bucket_id].m_data;
            break;
        }
    }

    return sizeFieldThisBucketInBytes;
}

void ContiguousFieldDataManager::add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord )
{
    for (size_t field_index=0;field_index<allFields.size();field_index++)
    {
        const FieldBase& field = *allFields[field_index];
        unsigned field_ordinal = field.mesh_meta_data_ordinal();
        if ( field.entity_rank() == dst_rank )
        {
            int numBytesPerEntity = field.get_meta_data_for_field()[dst_bucket_id].m_bytes_per_entity;
            if ( numBytesPerEntity > 0 )
            {
                size_t currentFieldSize = m_num_bytes_used_per_field[field_ordinal];
                size_t newFieldSize = currentFieldSize + numBytesPerEntity;
                unsigned char* new_field_data = m_field_raw_data[field_ordinal];
                bool requiresNewAllocation = newFieldSize > m_num_bytes_allocated_per_field[field_ordinal];
                if (requiresNewAllocation)
                {
                    new_field_data = field_data_allocator().allocate(newFieldSize + m_extra_capacity);
                }

                FieldMetaDataVector& field_meta_data_vector = const_cast<FieldMetaDataVector&>(field.get_meta_data_for_field());
                if (requiresNewAllocation)
                {
                    resetFieldMetaDataPointers(0, dst_bucket_id, field_meta_data_vector, m_field_raw_data[field_ordinal], new_field_data);
                }

                FieldMetaData &field_meta_data_for_modified_bucket = field_meta_data_vector[dst_bucket_id];
                size_t sizeOfPreviousBuckets = field_meta_data_for_modified_bucket.m_data - m_field_raw_data[field_ordinal];
                size_t leftHalf =  sizeOfPreviousBuckets + dst_bucket_ord*numBytesPerEntity;
                size_t rightHalf = m_num_bytes_used_per_field[field_ordinal] - leftHalf;
                field_meta_data_for_modified_bucket.m_data = new_field_data + sizeOfPreviousBuckets;

                size_t numBucketsOfRank = field_meta_data_vector.size();
                resetFieldMetaDataPointers(dst_bucket_id+1, numBucketsOfRank, field_meta_data_vector, m_field_raw_data[field_ordinal], new_field_data+numBytesPerEntity);

                if (requiresNewAllocation)
                {
                    std::memcpy(new_field_data, m_field_raw_data[field_ordinal],  leftHalf);
                    std::memcpy(new_field_data+leftHalf+numBytesPerEntity, m_field_raw_data[field_ordinal]+leftHalf, rightHalf);
                    field_data_allocator().deallocate(m_field_raw_data[field_ordinal], m_num_bytes_allocated_per_field[field_ordinal]);
                    m_num_bytes_allocated_per_field[field_ordinal] = newFieldSize + m_extra_capacity;
                }
                else
                {
                    std::memmove(new_field_data+leftHalf+numBytesPerEntity, m_field_raw_data[field_ordinal]+leftHalf, rightHalf);
                }

                setInitialValue(new_field_data+leftHalf, field, numBytesPerEntity);

                m_num_bytes_used_per_field[field_ordinal] += numBytesPerEntity;
                m_field_raw_data[field_ordinal] = new_field_data;
                ThrowRequire(dst_bucket_ord == m_num_entities_in_field_for_bucket[field_ordinal][dst_bucket_id]);
                m_num_entities_in_field_for_bucket[field_ordinal][dst_bucket_id] += 1;
            }
        }
    }
}

void ContiguousFieldDataManager::swap_fields(const int field1, const int field2)
{
    std::swap(m_field_raw_data[field1], m_field_raw_data[field2]);
    ThrowRequire(m_num_bytes_allocated_per_field[field1] == m_num_bytes_allocated_per_field[field2]);
    ThrowRequire(m_num_bytes_used_per_field[field1] == m_num_bytes_used_per_field[field2]);
    ThrowRequire(m_num_entities_in_field_for_bucket[field1].size() == m_num_entities_in_field_for_bucket[field2].size());
}

}
}
