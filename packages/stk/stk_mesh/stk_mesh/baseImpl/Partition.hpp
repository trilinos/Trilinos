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

#ifndef STK_MESH_IMPL_PARTITION_HPP_
#define STK_MESH_IMPL_PARTITION_HPP_

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for lower_bound
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PartOrdinal, EntityRank
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include <stk_mesh/base/EntitySorterBase.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { namespace impl { class BucketRepository; } } }

namespace stk {
namespace mesh {
namespace impl {

class Partition
{
public:
  enum RemoveMode {FILL_HOLE_THEN_SORT, TRACK_THEN_SLIDE};

  Partition(BulkData& mesh, BucketRepository *repo, EntityRank rank,
            const PartOrdinal* keyBegin, const PartOrdinal* keyEnd);

  virtual ~Partition();

  EntityRank get_rank() const { return m_rank; }

  bool empty() const { return m_size == 0; }

  size_t size() const { return m_size; }

  std::pair<const unsigned *, const unsigned *> superset_part_ordinals() const { return m_partOrdsBeginEnd; }

  const std::vector<PartOrdinal> &get_legacy_partition_id() const { return m_extPartitionKey; }

  const unsigned * key() const { return m_extPartitionKey.data(); }

  /// Add an entity to this partition.  The entity must not be a member
  /// of another partition.
  bool add(Entity entity);

  /// Move an entity from this partition to a destination partition.
  bool move_to(Entity entity, Partition &dst_partition);

  /// Remove an entity from this partition.
  bool remove(Entity entity);

  /// Sort the entities in this partition by EntityKey without changing
  /// the number or sizes of buckets.
  void default_sort_if_needed(bool mustSortFacesByNodeIds=false);
  void sort(const EntitySorterBase& sorter);

  void set_flag_needs_to_be_sorted(bool flag) { m_updated_since_sort = flag; }

  bool needs_to_be_sorted() const { return m_updated_since_sort; }

  ////
  //// This part of the interface exposes the Buckets that are currently a part of
  //// the implementation.
  ////

  /// Does the given bucket belong to this partition.
  bool belongs(const Bucket &bkt) const { return bkt.getPartition() == this;}

  size_t num_buckets() const { return m_buckets.size();}

  inline BucketVector::iterator begin() { return m_buckets.begin(); }
  inline BucketVector::iterator end() { return m_buckets.end(); }

  inline BucketVector::const_iterator begin() const { return m_buckets.begin(); }
  inline BucketVector::const_iterator end() const { return m_buckets.end(); }

  ////
  //// The following are used internally and for unit testing.
  ////

  size_t compute_size()
  {
    size_t partition_size = 0;
    for (BucketVector::const_iterator b_i = begin(), b_e = end(); b_i != b_e; ++b_i)
    {
      partition_size += (*b_i)->size();
    }
    m_size = partition_size;

    return partition_size;
  }

  // Enables sidestepping compiler fussiness wrt overloaded operator<<(..).
  std::ostream &streamit(std::ostream &os) const;

  // Output including Entities.
  std::ostream &dumpit(std::ostream &os) const;
  std::string dumpit() const;

  void delete_bucket(Bucket* bucket);

  void remove_bucket(Bucket* bucket);

  void add_bucket(Bucket* bucket);

  void reset_partition_key(const std::vector<unsigned>& newKey);

  void set_remove_mode(RemoveMode removeMode);
  RemoveMode get_remove_mode() const { return m_removeMode; }

private:
  BulkData& m_mesh;
  BucketRepository *m_repository;

  EntityRank m_rank;

  void check_sorted(const std::string& prefixMsg);

  std::vector<PartOrdinal> m_extPartitionKey;
  std::pair<const PartOrdinal*,const PartOrdinal*> m_partOrdsBeginEnd;

  // Used if the set of buckets (not just bucket contents) are being modified.
  BucketVector m_buckets;

  // Number of entities in this partition.
  size_t m_size;

  bool m_updated_since_sort;

  RemoveMode m_removeMode;

  std::vector<FastMeshIndex> m_removedEntities;
  //
  // Internal methods
  //

  void remove_internal(Bucket& bucket, unsigned bucketOrd);
  void remove_impl();

  void clear_pending_removes_by_filling_from_end();
  void finalize_pending_removes_by_sliding_memory();

  // The partition has no buckets, not even an empty one left after removing all its
  // entities.
  bool no_buckets() const { return m_buckets.empty(); }

  // Make sure that the last bucket has room for an entity to be added to it, adding
  // an empty bucket if necessary.
  Bucket *get_bucket_for_adds();

  void internal_check_invariants() const
  {
#ifndef NDEBUG
    internal_check_no_null_buckets_invariant();
#endif
  }

  void internal_check_no_null_buckets_invariant() const
  {
#ifndef NDEBUG
    for (int i = 0, e = m_buckets.size(); i < e; ++i) {
      STK_ThrowAssert(m_buckets[i] != NULL);
    }
#endif
  }

  // Overwrite the location defined by the input arguments with the entity
  // at the very end of this entire partition
  void overwrite_from_end( Bucket& bucket, unsigned ordinal);
};

std::ostream &operator<<(std::ostream &, const stk::mesh::impl::Partition &);

struct PartitionLess {
  bool operator()( const Partition * lhs_Partition , const OrdinalVector& rhs ) const
  {
    return lhs_Partition->get_legacy_partition_id().size() != rhs.size() ?
           lhs_Partition->get_legacy_partition_id().size() < rhs.size() :
           lhs_Partition->get_legacy_partition_id() < rhs;
  }

  bool operator()( const OrdinalVector& lhs , const Partition * rhs_Partition ) const
  {
    return lhs.size() != rhs_Partition->get_legacy_partition_id().size() ?
           lhs.size() < rhs_Partition->get_legacy_partition_id().size() :
           lhs < rhs_Partition->get_legacy_partition_id();
  }
};

inline
std::vector<Partition*>::iterator
upper_bound( std::vector<Partition*> & v , const OrdinalVector& key )
{ return std::upper_bound( v.begin() , v.end() , key , PartitionLess() ); }

} // impl
} // mesh
} // stk

#endif /* STK_MESH_IMPL_PARTITION_HPP_ */

