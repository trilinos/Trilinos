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
                   unsigned alignmentPaddingSize = STK_ALIGNMENT_PADDING_SIZE);

  ~FieldDataManager() = default;

  // Called when creating new Bucket.  Initializes storage for all Fields.
  void allocate_bucket_field_data(const EntityRank rank, const std::vector<FieldBase *> & fieldSet,
                                  const PartVector& supersetParts, unsigned size, unsigned capacity);

  // Called when removing existing Bucket
  void deallocate_bucket_field_data(const EntityRank rank, const unsigned bucketId, const size_t capacity,
                                    const std::vector<FieldBase*>&  fields);

  // Called when reordering buckets (leaving contents alone)
  void reorder_bucket_field_data(EntityRank rank, const std::vector<FieldBase*> & fields,
                                 const std::vector<unsigned>& reorderedBucketIds);

  // Called when initially allocating all Buckets.  Initializes storage for all Fields.
  void allocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets,
                           const std::vector< FieldBase * > & fieldSet);

  // Called when adding a late Field.  Initializes storage for new Field.
  void reallocate_field_data(EntityRank rank, const std::vector<Bucket*>& buckets, FieldBase & currentField,
                             const std::vector<FieldBase *> & allFields);

  // Called when adding an entity to a Bucket.  Initializes storage for new Entity.
  void add_field_data_for_entity(const std::vector<FieldBase *> &allFields, EntityRank dstRank,
                                 unsigned dstBucketId, unsigned dstBucketOrd, unsigned newBucketSize);

  // Called when removing an entity from a Bucket.
  void remove_field_data_for_entity(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &allFields);

  // Called when growing a Bucket before adding an Entity
  void grow_bucket_capacity(const FieldVector & allFields, EntityRank rank, unsigned bucketId,
                            unsigned bucketSize, unsigned bucketCapacity);

  // Updates STK_FIELD_ASAN memory poisoning around current active data
  void reset_empty_field_data(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                              unsigned bucketCapacity, const FieldVector & fields);

  unsigned get_bucket_capacity(EntityRank rank, unsigned bucketId) const {
    return m_bucketCapacity[rank][bucketId];
  }

  size_t get_num_bytes_allocated_on_field(const unsigned fieldIndex) const {
    return m_numBytesAllocatedPerField[fieldIndex];
  }

  unsigned get_alignment_padding_size() const { return m_alignmentPaddingSize; }

  const FieldDataAllocator<std::byte>& get_field_data_allocator() const { return m_fieldDataAllocator; }

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
  void initialize_entity_field_data(EntityRank rank, unsigned bucketId, unsigned bucketOrd, unsigned newBucketSize,
                                    const std::vector<FieldBase *> &fields);

  std::vector<BucketFieldSegment> get_old_bucket_field_offsets(const EntityRank rank,
                                                               const unsigned bucketId,
                                                               const std::vector<FieldBase*>& allFields,
                                                               const unsigned capacity) const;
  std::vector<BucketFieldSegment> get_new_bucket_field_offsets(const EntityRank rank,
                                                               const unsigned bucketId,
                                                               const std::vector<FieldBase*>& allFields,
                                                               const unsigned capacity) const;

  using AllocationType = FieldDataAllocator<std::byte>::HostAllocationType;

  FieldDataAllocator<std::byte> m_fieldDataAllocator;
  unsigned m_alignmentPaddingSize;
  std::vector<std::vector<AllocationType>> m_fieldRawData;
  std::vector<std::vector<unsigned>> m_bucketCapacity;
  std::vector<size_t> m_numBytesAllocatedPerField;
};

}
}

#endif
