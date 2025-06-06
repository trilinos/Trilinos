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

#include <stk_mesh/base/Types.hpp>      // for EntityRank, PartVector
#include <stk_util/util/FieldDataAllocator.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include <cstddef>
#include <vector>                       // for vector

namespace stk { namespace mesh { class FieldBase; } }

namespace stk {
namespace mesh {

class AllocatorAdaptorInterface {
public:
  virtual ~AllocatorAdaptorInterface() {}

  using pointer = std::byte*;

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
  virtual void allocate_bucket_field_data(const EntityRank rank, const std::vector< FieldBase * > & fieldSet,
                                          const PartVector& supersetParts, unsigned size, unsigned capacity) = 0;
  virtual void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                            const std::vector<FieldBase*>& fields) = 0;
  virtual void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                         const std::vector<unsigned>& reorderedBucketIds) = 0;
  virtual void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                                   const std::vector< FieldBase * > & fieldSet) = 0;
  virtual void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField,
                                     const std::vector<FieldBase *> & allFields) = 0;
  virtual size_t get_num_bytes_allocated_on_field(const unsigned fieldIndex) const = 0;
  virtual void add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                         unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize) = 0;
  virtual void remove_field_data_for_entity(EntityRank rank, unsigned bucketId, unsigned bucketOrd,
                                            unsigned newBucketSize, const std::vector<FieldBase *> &allFields) = 0;
  virtual void initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd,
                                            unsigned newBucketSize, const std::vector<FieldBase *> &fields) = 0;
  virtual void swap_fields(const int field1, const int field2) = 0;
  virtual unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const = 0;
  virtual void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                                    unsigned bucketSize, unsigned bucketCapacity) = 0;
  virtual void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                      unsigned bucketCapacity, const FieldVector & fields) = 0;

  unsigned get_alignment_bytes() const { return m_alignmentIncrementBytes; }

protected:
  std::unique_ptr<AllocatorAdaptorInterface> m_fieldDataAllocator;
  size_t m_alignmentIncrementBytes;
};

class DefaultFieldDataManager : public FieldDataManager
{
public:
  DefaultFieldDataManager(const unsigned num_ranks,
                          unsigned alignmentIncrementBytes = stk::impl::DEFAULT_FIELD_ALIGNMENT_BYTES)
    : FieldDataManager(alignmentIncrementBytes),
      m_fieldRawData(num_ranks),
      m_bucketCapacity(num_ranks),
      m_numBytesAllocatedPerField()
  {
  }

  DefaultFieldDataManager(const unsigned num_ranks,
                          std::unique_ptr<AllocatorAdaptorInterface> userSpecifiedAllocatorAdaptor,
                          unsigned alignmentIncrementBytes)
    : FieldDataManager(alignmentIncrementBytes, std::move(userSpecifiedAllocatorAdaptor)),
      m_fieldRawData(num_ranks),
      m_bucketCapacity(num_ranks),
      m_numBytesAllocatedPerField()
  {
  }

  virtual ~DefaultFieldDataManager() override = default;
  void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & fieldSet,
                                  const PartVector& supersetParts, unsigned size, unsigned capacity) override;
  void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                    const std::vector<FieldBase*>&  fields) override;
  void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                 const std::vector<unsigned>& reorderedBucketIds) override;
  void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                           const std::vector< FieldBase * > & fieldSet) override;
  void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField,
                             const std::vector<FieldBase *> & allFields) override;
  size_t get_num_bytes_allocated_on_field(const unsigned fieldIndex) const override {
    return m_numBytesAllocatedPerField[fieldIndex];
  }
  void add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                 unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize) override;
  void remove_field_data_for_entity(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &allFields) override;
  void initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &fields) override;
  void swap_fields(const int /*field1*/, const int /*field2*/) override { }
  unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const override {
    return m_bucketCapacity[rank][bucketId];
  }
  void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                            unsigned bucketSize, unsigned bucketCapacity) override;
  void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                              unsigned bucketCapacity, const FieldVector & fields) override;

private:
  void allocate_new_field_meta_data(const EntityRank rank, const unsigned bucketId,
                                    const std::vector<FieldBase*>& allFields);
  void reallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, FieldBase & currentField,
                                    const std::vector<FieldBase *> & allFields, const PartVector& supersetParts,
                                    unsigned bucketSize, unsigned bucketCapacity);
  void allocate_bucket_ordinal_field_data(const EntityRank rank, const std::vector<FieldBase *> & fieldSet,
                                          const PartVector& supersetParts, unsigned bucketOrd, unsigned size,
                                          unsigned capacity);
  void resize_bucket_arrays(const EntityRank rank, const std::vector<FieldBase*>& allFields, int newNumBuckets);

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
                              const PartVector & supersetParts,
                              unsigned bucketSize,
                              unsigned bucketCapacity);
  void copy_field_data_from_old_to_new_bucket(EntityRank rank,
                                              unsigned bucketSize,
                                              unsigned bucketId,
                                              const std::vector<FieldBase*>& allFields,
                                              const std::vector<BucketFieldSegment>& oldOffsetForField,
                                              const std::vector<BucketFieldSegment>& newOffsetForField,
                                              const std::byte* oldAllocationAllFields,
                                              std::byte* newAllocationAllFields);
  void update_field_pointers_to_new_bucket(const EntityRank rank,
                                           const unsigned bucketId,
                                           const std::vector<FieldBase*>& allFields,
                                           const size_t capacity,
                                           std::byte* newAllocationAllFields);
  void initialize_new_field_values(FieldBase& currentField, const EntityRank rank, const unsigned bucketId,
                                   unsigned size, unsigned capacity);

  std::vector<std::vector<std::byte*>> m_fieldRawData;
  std::vector<std::vector<unsigned>> m_bucketCapacity;
  std::vector<size_t> m_numBytesAllocatedPerField;
};

class ContiguousFieldDataManager : public FieldDataManager
{
public:
  ContiguousFieldDataManager(unsigned extraCapacity = 8192,
                             unsigned alignmentIncrementBytes = stk::impl::DEFAULT_FIELD_ALIGNMENT_BYTES)
    : FieldDataManager(alignmentIncrementBytes),
      m_fieldRawData(),
      m_numBytesAllocatedPerField(),
      m_extraCapacity(extraCapacity),
      m_numBytesUsedPerField(),
      m_numEntitiesInFieldForBucket()
  {
  }

  virtual ~ContiguousFieldDataManager() override;
  void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & fieldSet,
                                  const PartVector& supersetParts, unsigned size, unsigned capacity) override;
  void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                    const std::vector<FieldBase*>&  fields) override;
  void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                 const std::vector<unsigned>& reorderedBucketIds) override;
  void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                           const std::vector< FieldBase * > & fieldSet) override;
  void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField,
                             const std::vector<FieldBase *> & allFields) override;
  size_t get_num_bytes_allocated_on_field(const unsigned fieldIndex) const override {
    return m_numBytesAllocatedPerField[fieldIndex];
  }
  void add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                 unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize) override;
  void remove_field_data_for_entity(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &allFields) override;
  void initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &fields) override;
  void swap_fields(const int field1, const int field2) override;
  void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                              unsigned bucketCapacity, const FieldVector & fields) override;

  const std::vector<std::byte*> &get_field_raw_data() const {return m_fieldRawData;}
  const std::vector<size_t> &get_num_bytes_allocated_per_field_array() const {
    return m_numBytesAllocatedPerField;
  }
  const std::vector<size_t> &get_num_bytes_used_per_field_array() const {return m_numBytesUsedPerField;}
  size_t get_extra_capacity() const { return m_extraCapacity; }
  unsigned get_bucket_capacity(EntityRank /*rank*/, unsigned /*bucketId*/) const override { return 0; }
  void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                            unsigned bucketSize, unsigned bucketCapacity) override;

private:
  void clear_bucket_field_data(const EntityRank rmRank, const unsigned rmBucketId,
                               const std::vector<FieldBase*>& allFields);
  void allocate_new_field_meta_data(const EntityRank rank, const std::vector<Bucket*> & buckets,
                                    const std::vector<FieldBase*>& allFields);
  std::vector<size_t> get_field_bucket_offsets(const std::vector<Bucket*> & buckets, FieldBase & currentField) const;
  void copy_bucket_data_from_old_to_new_field(const std::vector<size_t>& oldOffsetForBucket,
                                              const std::vector<size_t>& newOffsetForBucket,
                                              const std::byte* oldAllocationAllBuckets,
                                              std::byte* newAllocationAllBuckets);
  void update_bucket_storage_for_field(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase& currentField);
  void update_bucket_pointers_to_new_field(const std::vector<Bucket*>& buckets, FieldBase& currentField);
  void initialize_new_bucket_values(const std::vector<Bucket*>& buckets,
                                    const std::vector<size_t> & oldOffsetForBucket,
                                    const std::vector<size_t> & newOffsetForBucket,
                                    FieldBase& currentField);

  std::vector<std::byte*> m_fieldRawData;
  std::vector<size_t> m_numBytesAllocatedPerField;
  size_t m_extraCapacity;
  std::vector<size_t> m_numBytesUsedPerField;
  std::vector<std::vector<size_t>> m_numEntitiesInFieldForBucket;
};

}
}

#endif
