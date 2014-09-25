// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef stk_mesh_FieldDataManager_hpp
#define stk_mesh_FieldDataManager_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Types.hpp>      // for EntityRank, PartVector
#include <stk_util/util/PageAlignedAllocator.hpp>
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::size_type
namespace stk { namespace mesh { class FieldBase; } }

namespace stk {
namespace mesh {


class FieldDataManager
{
public:
    typedef page_aligned_allocator<unsigned char, FieldDataTag> field_data_allocator;

    FieldDataManager() {}
    virtual ~FieldDataManager() {}
    virtual void allocate_bucket_field_data(const EntityRank rank,
            const std::vector< FieldBase * > & field_set, const PartVector& superset_parts, const size_t capacity) = 0;
    virtual void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
            const std::vector<FieldBase*>  fields) = 0;
    virtual void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds) = 0;
    virtual void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set) = 0;
    virtual size_t get_num_bytes_allocated_on_field(const unsigned field_index) const = 0;
    virtual void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord ) = 0;
    virtual void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &allFields) = 0;
    virtual void reinitialize_removed_entity_field_data(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields) = 0;
    virtual void swap_fields(const int field1, const int field2) = 0;
};

class DefaultFieldDataManager : public FieldDataManager
{
public:
    DefaultFieldDataManager(const size_t num_ranks) :
        m_field_raw_data(num_ranks), m_num_bytes_allocated_per_field()
    {
    }
    ~DefaultFieldDataManager() {}
    void allocate_bucket_field_data(const EntityRank rank,
            const std::vector< FieldBase * > & field_set, const PartVector& superset_parts, const size_t capacity);
    void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
            const std::vector<FieldBase*>  fields);
    void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds);
    void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set);
    size_t get_num_bytes_allocated_on_field(const unsigned field_index) const { return m_num_bytes_allocated_per_field[field_index]; }
    void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord );
    void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &allFields);
    void reinitialize_removed_entity_field_data(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields);
    void swap_fields(const int field1, const int field2) { }

private:
    std::vector<std::vector<unsigned char*> > m_field_raw_data;
    std::vector<size_t> m_num_bytes_allocated_per_field;
};

class ContiguousFieldDataManager : public FieldDataManager
{
public:
    ContiguousFieldDataManager() : m_field_raw_data(), m_num_bytes_allocated_per_field(), m_extra_capacity(8192), m_num_bytes_used_per_field(), m_num_entities_in_field_for_bucket() {}
    ~ContiguousFieldDataManager();
    void allocate_bucket_field_data(const EntityRank rank,
            const std::vector< FieldBase * > & field_set, const PartVector& superset_parts, const size_t capacity);
    void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
            const std::vector<FieldBase*>  fields);
    void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds);
    void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set);
    size_t get_num_bytes_allocated_on_field(const unsigned field_index) const { return m_num_bytes_allocated_per_field[field_index]; }
    void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord );
    void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &allFields);
    void reinitialize_removed_entity_field_data(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord, const std::vector<FieldBase *> &fields);
    void swap_fields(const int field1, const int field2);

    const std::vector<unsigned char*> &get_field_raw_data() const {return m_field_raw_data;}
    const std::vector<size_t> &get_num_bytes_allocated_per_field_array() const {return m_num_bytes_allocated_per_field;}
    const std::vector<size_t> &get_num_bytes_used_per_field_array() const {return m_num_bytes_used_per_field;}
    size_t get_extra_capacity() const { return m_extra_capacity; }

private:
    void clear_bucket_field_data(const EntityRank rm_rank, const unsigned rm_bucket_id, const std::vector<FieldBase*>  &all_fields);

    std::vector<unsigned char*> m_field_raw_data;
    std::vector<size_t> m_num_bytes_allocated_per_field;
    size_t m_extra_capacity;
    std::vector<size_t> m_num_bytes_used_per_field;
    std::vector<std::vector<size_t> >m_num_entities_in_field_for_bucket;
};

}
}

#endif
