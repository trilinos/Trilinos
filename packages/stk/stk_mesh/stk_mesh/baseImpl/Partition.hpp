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

/*
 * Partition.hpp
 */

#ifndef STK_MESH_IMPL_PARTITION_HPP_
#define STK_MESH_IMPL_PARTITION_HPP_

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for lower_bound
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PartOrdinal, EntityRank
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { namespace impl { class BucketRepository; } } }
namespace stk { namespace mesh { namespace utest { struct SyncToPartitions; } } }
namespace stk { namespace mesh { struct Entity; } }

namespace stk {
namespace mesh {

namespace utest {
}

namespace impl {


class Partition
{
  friend class BucketRepository;
  friend struct stk::mesh::utest::SyncToPartitions;

public:

  Partition(BulkData& mesh, BucketRepository *repo, EntityRank rank,
            const std::vector<PartOrdinal> &key);

  virtual ~Partition();

  ////
  //// The main part of the interface is Bucket-free.
  ////

  /// Rank of the entities in this partition.
  EntityRank get_rank() const { return m_rank; }

  bool empty() const { return m_size == 0; }

  size_t size() const { return m_size; }

  /// Returns the representation used by BucketRepository to identify a bucket,
  /// including the parts it corresponds to.
  const std::vector<PartOrdinal> &get_legacy_partition_id() const { return m_extPartitionKey; }

  const unsigned * key() const { return &m_extPartitionKey[0]; }

  /// Add an entity to this partition.  The entity must not be a member
  /// of another partition.
  bool add(Entity entity);

  /// Move an entity from this partion to a destination partition.
  bool move_to(Entity entity, Partition &dst_partition);

  /// Remove an entity from this partition.
  bool remove(Entity entity);

  /// Compress this partion into a single bucket of sorted Entities.
  void compress(bool force = false);

  /// Sort the entities in this partition by EntityKey without changing
  /// the number or sizes of buckets.
  void sort(bool force = false);

  void set_flag_needs_to_be_sorted(bool flag) { m_updated_since_sort = flag; }

  size_t field_data_footprint(const FieldBase &f) const;

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

  // Just for unit testing.  Remove after refactor.
  static BucketRepository &getRepository(stk::mesh::BulkData &mesh);

private:

  //
  // Members
  //

  BulkData& m_mesh;
  BucketRepository *m_repository;

  EntityRank m_rank;

  // Identifies the partition, borrowing the representation from BucketRepository.
  std::vector<PartOrdinal> m_extPartitionKey;

  // Used if the set of buckets (not just bucket contents) are being modified.
  BucketVector m_buckets;

  // Number of entities in this partition.
  size_t m_size;

  bool m_updated_since_compress;

  bool m_updated_since_sort;

  //
  // Internal methods
  //

  void remove_impl();

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
    internal_check_size_invariant();
#endif
  }

  void internal_check_size_invariant() const
  {
#ifndef NDEBUG
    size_t sum = 0;
    for (size_t i = 0, e = m_buckets.size(); i < e; ++i) {
      sum += m_buckets[i]->size();
      m_buckets[i]->check_size_invariant();
    }
    ThrowAssertMsg(sum == m_size, "Inconsistent sizes, bucket sum is " << sum << ", m_size is " << m_size);
#endif
  }

  void internal_check_no_null_buckets_invariant() const
  {
#ifndef NDEBUG
    for (int i = 0, e = m_buckets.size(); i < e; ++i) {
      ThrowAssert(m_buckets[i] != NULL);
    }
#endif
  }

  void internal_swap_to_end(Entity entity);

  // Overwrite the location defined by the input arguments with the entity
  // at the very end of this entire partition
  void overwrite_from_end( Bucket& bucket, unsigned ordinal);
};

std::ostream &operator<<(std::ostream &, const stk::mesh::impl::Partition &);

struct PartitionLess {
  bool operator()( const Partition * lhs_Partition , const unsigned * rhs ) const ;
  bool operator()( const unsigned * lhs , const Partition * rhs_Partition ) const ;
};

inline
bool partition_key_less( const unsigned * lhs , const unsigned * rhs )
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

// The part count and part ordinals are less
inline bool PartitionLess::operator()( const Partition * lhs_partition ,
                                       const unsigned * rhs ) const
{ return partition_key_less( lhs_partition->key() , rhs ); }

inline bool PartitionLess::operator()( const unsigned * lhs ,
                                       const Partition * rhs_partition ) const
{ return partition_key_less( lhs , rhs_partition->key() ); }

inline
std::vector<Partition*>::iterator
lower_bound( std::vector<Partition*> & v , const unsigned * key )
{ return std::lower_bound( v.begin() , v.end() , key , PartitionLess() ); }

} // impl
} // mesh
} // stk

#endif /* PartitionFAMILY_HPP_ */
