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
  FieldDataManager(const unsigned num_ranks,
                   unsigned alignmentIncrementBytes = stk::impl::DEFAULT_FIELD_ALIGNMENT_BYTES,
                   std::unique_ptr<AllocatorAdaptorInterface> allocatorAdaptor =
                       std::unique_ptr<AllocatorAdaptorInterface>());

  ~FieldDataManager() = default;

  void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & fieldSet,
                                  const PartVector& supersetParts, unsigned size, unsigned capacity);
  void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                    const std::vector<FieldBase*>&  fields);
  void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                 const std::vector<unsigned>& reorderedBucketIds);
  void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                           const std::vector< FieldBase * > & fieldSet);
  void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField,
                             const std::vector<FieldBase *> & allFields);
  size_t get_num_bytes_allocated_on_field(const unsigned fieldIndex) const {
    return m_numBytesAllocatedPerField[fieldIndex];
  }
  void add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                 unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize);
  void remove_field_data_for_entity(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &allFields);
  void initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &fields);
  void swap_fields(const int /*field1*/, const int /*field2*/) { }
  unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const {
    return m_bucketCapacity[rank][bucketId];
  }
  void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                            unsigned bucketSize, unsigned bucketCapacity);
  void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                              unsigned bucketCapacity, const FieldVector & fields);

  unsigned get_alignment_bytes() const { return m_alignmentIncrementBytes; }

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

  std::unique_ptr<AllocatorAdaptorInterface> m_fieldDataAllocator;
  unsigned m_alignmentIncrementBytes;
  std::vector<std::vector<std::byte*>> m_fieldRawData;
  std::vector<std::vector<unsigned>> m_bucketCapacity;
  std::vector<size_t> m_numBytesAllocatedPerField;
};

}
}

#endif
