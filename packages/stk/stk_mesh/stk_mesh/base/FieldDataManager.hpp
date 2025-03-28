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

#ifndef stk_mesh_FieldDataManager_hpp
#define stk_mesh_FieldDataManager_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Types.hpp>      // for EntityRank, PartVector
#include <stk_util/util/FieldDataAllocator.hpp>
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
namespace stk { namespace mesh { class FieldBase; } }

namespace stk {
namespace mesh {

class AllocatorAdaptorInterface {
public:
  virtual ~AllocatorAdaptorInterface() {}

  typedef unsigned char* pointer;

  virtual pointer allocate(size_t num, const void* = 0) = 0;
  virtual void deallocate(pointer p, size_t num) = 0;
};

template<class AllocatorType>
class AllocatorAdaptor : public AllocatorAdaptorInterface {
public:
  AllocatorAdaptor(){}
  virtual ~AllocatorAdaptor(){}

  pointer allocate(size_t num, const void* = 0) override {
    return AllocatorType().allocate(num);
  }

  void deallocate(pointer p, size_t num) override {
    AllocatorType().deallocate(p, num);
  }
};

struct BucketFieldSegment {
  BucketFieldSegment(int offset_, int size_)
    : offset(offset_),
      size(size_)
  {}

  int offset;
  int size;
};

class FieldDataManager
{
public:
    FieldDataManager(unsigned alignmentIncrementBytes,
                     std::unique_ptr<AllocatorAdaptorInterface> allocatorAdaptor =
                         std::unique_ptr<AllocatorAdaptorInterface>());

    virtual ~FieldDataManager() = default;
    virtual void allocate_bucket_field_data(const EntityRank rank, const std::vector< FieldBase * > & field_set,
                                            const PartVector& superset_parts, unsigned size, unsigned capacity) = 0;
    virtual void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity, const std::vector<FieldBase*>& fields) = 0;
    virtual void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds) = 0;
    virtual void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set) = 0;
    virtual void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField, const std::vector<FieldBase *> & allFields) = 0;
    virtual size_t get_num_bytes_allocated_on_field(const unsigned field_index) const = 0;
    virtual void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, unsigned dst_bucket_ord ) = 0;
    virtual void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &allFields) = 0;
    virtual void initialize_entity_field_data(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields) = 0;
    virtual void swap_fields(const int field1, const int field2) = 0;
    virtual unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const = 0;
    virtual void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                                      unsigned bucketSize, unsigned bucketCapacity) = 0;
    virtual void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                        unsigned bucketCapacity, const FieldVector & fields) = 0;

    unsigned get_alignment_bytes() const { return alignment_increment_bytes; }

protected:
    std::unique_ptr<AllocatorAdaptorInterface> m_fieldDataAllocator;
    size_t alignment_increment_bytes;
};

class DefaultFieldDataManager : public FieldDataManager
{
public:
    DefaultFieldDataManager(const unsigned num_ranks,
                            unsigned alignmentIncrementBytes = stk::impl::DEFAULT_FIELD_ALIGNMENT_BYTES)
    : FieldDataManager(alignmentIncrementBytes),
      m_field_raw_data(num_ranks),
      m_bucketCapacity(num_ranks),
      m_num_bytes_allocated_per_field()
    {
    }

    DefaultFieldDataManager(const unsigned num_ranks,
                            std::unique_ptr<AllocatorAdaptorInterface> userSpecifiedAllocatorAdaptor,
                            unsigned alignmentIncrementBytes)
    : FieldDataManager(alignmentIncrementBytes, std::move(userSpecifiedAllocatorAdaptor)),
      m_field_raw_data(num_ranks),
      m_bucketCapacity(num_ranks),
      m_num_bytes_allocated_per_field()
    {
    }

    virtual ~DefaultFieldDataManager() override = default;
    void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & field_set,
                                    const PartVector& superset_parts, unsigned size, unsigned capacity) override;
    void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
            const std::vector<FieldBase*>&  fields) override;
    void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds) override;
    void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set) override;
    void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField, const std::vector<FieldBase *> & allFields) override;
    size_t get_num_bytes_allocated_on_field(const unsigned field_index) const override { return m_num_bytes_allocated_per_field[field_index]; }
    void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, unsigned dst_bucket_ord ) override;
    void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &allFields) override;
    void initialize_entity_field_data(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields) override;
    void swap_fields(const int /*field1*/, const int /*field2*/) override { }
    unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const override { return m_bucketCapacity[rank][bucketId]; }
    void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                              unsigned bucketSize, unsigned bucketCapacity) override;
    void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                unsigned bucketCapacity, const FieldVector & fields) override;

private:
    void allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId, const std::vector<FieldBase*>& allFields);
    void reallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, FieldBase & currentField,
                                      const std::vector<FieldBase *> & allFields, const PartVector& superset_parts,
                                      unsigned bucketSize, unsigned bucketCapacity);

    std::vector<BucketFieldSegment> get_old_bucket_field_offsets(const EntityRank rank,
                                                                 const unsigned bucketId,
                                                                 const std::vector<FieldBase*>& allFields,
                                                                 const unsigned capacity) const;
    std::vector<BucketFieldSegment> get_new_bucket_field_offsets(const EntityRank rank,
                                                                 const unsigned bucketId,
                                                                 const std::vector<FieldBase*>& allFields,
                                                                 const unsigned capacity) const;
    void update_field_meta_data(const EntityRank rank, const unsigned bucketId,
                              const std::vector<FieldBase*> & allFields,
                              const PartVector & supersetParts);
    void copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                                unsigned bucketSize,
                                                unsigned bucketId,
                                                const std::vector<FieldBase*>& allFields,
                                                const std::vector<BucketFieldSegment>& oldOffsetForField,
                                                const std::vector<BucketFieldSegment>& newOffsetForField,
                                                const unsigned char* oldAllocationAllFields,
                                                unsigned char* newAllocationAllFields);
    void update_field_pointers_to_new_bucket(const EntityRank rank,
                                             const unsigned bucketId,
                                             const std::vector<FieldBase*>& allFields,
                                             const size_t capacity,
                                             unsigned char* newAllocationAllFields);
    void initialize_new_field_values(FieldBase& currentField, const EntityRank rank, const unsigned bucketId,
                                     unsigned size, unsigned capacity);

    std::vector<std::vector<unsigned char*>> m_field_raw_data;
    std::vector<std::vector<unsigned>> m_bucketCapacity;
    std::vector<size_t> m_num_bytes_allocated_per_field;
};

class ContiguousFieldDataManager : public FieldDataManager
{
public:
    ContiguousFieldDataManager(unsigned extraCapacity = 8192,
                               unsigned alignmentIncrementBytes = stk::impl::DEFAULT_FIELD_ALIGNMENT_BYTES)
    : FieldDataManager(alignmentIncrementBytes),
      m_field_raw_data(),
      m_num_bytes_allocated_per_field(),
      m_extra_capacity(extraCapacity),
      m_num_bytes_used_per_field(),
      m_num_entities_in_field_for_bucket()
    {
    }

    virtual ~ContiguousFieldDataManager() override;
    void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & field_set,
                                    const PartVector& superset_parts, unsigned size, unsigned capacity) override;
    void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucket_id, const size_t capacity,
                                      const std::vector<FieldBase*>&  fields) override;
    void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields, const std::vector<unsigned>& reorderedBucketIds) override;
    void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, const std::vector< FieldBase * > & field_set) override;
    void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField, const std::vector<FieldBase *> & allFields) override;
    size_t get_num_bytes_allocated_on_field(const unsigned field_index) const override { return m_num_bytes_allocated_per_field[field_index]; }
    void add_field_data_for_entity(const std::vector<FieldBase *> &allFields,EntityRank dst_rank,unsigned dst_bucket_id, unsigned dst_bucket_ord ) override;
    void remove_field_data_for_entity(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &allFields) override;
    void initialize_entity_field_data(EntityRank rank, unsigned bucket_id, unsigned bucket_ord, const std::vector<FieldBase *> &fields) override;
    void swap_fields(const int field1, const int field2) override;
    void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                unsigned bucketCapacity, const FieldVector & fields) override;

    const std::vector<unsigned char*> &get_field_raw_data() const {return m_field_raw_data;}
    const std::vector<size_t> &get_num_bytes_allocated_per_field_array() const {return m_num_bytes_allocated_per_field;}
    const std::vector<size_t> &get_num_bytes_used_per_field_array() const {return m_num_bytes_used_per_field;}
    size_t get_extra_capacity() const { return m_extra_capacity; }
    unsigned get_bucket_capacity(EntityRank /*rank*/, unsigned /*bucketId*/) const override { return 0; }
    void grow_bucket_capacity(const FieldVector & /*allFields*/, EntityRank /*rank*/, unsigned /*bucketId*/,
                              unsigned /*bucketSize*/, unsigned /*bucketCapacity*/) override {}

private:
    void clear_bucket_field_data(const EntityRank rm_rank, const unsigned rm_bucket_id, const std::vector<FieldBase*>  &all_fields);
    void allocate_new_field_meta_data(const EntityRank rank, const std::vector<Bucket*> & buckets, const std::vector<FieldBase*>& allFields);
    std::vector<size_t> get_field_bucket_offsets(const std::vector<Bucket*> & buckets, FieldBase & currentField) const;
    void copy_bucket_data_from_old_to_new_field(const std::vector<size_t>& oldOffsetForBucket,
                                                const std::vector<size_t>& newOffsetForBucket,
                                                const unsigned char* oldAllocationAllBuckets,
                                                unsigned char* newAllocationAllBuckets);
    void update_bucket_storage_for_field(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase& currentField);
    void update_bucket_pointers_to_new_field(const std::vector<Bucket*>& buckets, FieldBase& currentField);
    void initialize_new_bucket_values(const std::vector<Bucket*>& buckets,
                                      const std::vector<size_t> & oldOffsetForBucket,
                                      const std::vector<size_t> & newOffsetForBucket,
                                      FieldBase& currentField);

    std::vector<unsigned char*> m_field_raw_data;
    std::vector<size_t> m_num_bytes_allocated_per_field;
    size_t m_extra_capacity;
    std::vector<size_t> m_num_bytes_used_per_field;
    std::vector<std::vector<size_t>> m_num_entities_in_field_for_bucket;
};

}
}

#endif
