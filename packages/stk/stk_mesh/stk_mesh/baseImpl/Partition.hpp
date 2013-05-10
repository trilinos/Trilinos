/*
 * Partition.hpp
 */

#ifndef STK_MESH_IMPL_PARTITION_HPP_
#define STK_MESH_IMPL_PARTITION_HPP_

#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

namespace utest {
struct SyncToPartitions;
}

namespace impl {

class BucketRepository;

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
  const EntityRank get_rank() const { return m_rank; }

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

  size_t field_data_footprint(const FieldBase &f) const;

  ////
  //// This part of the interface exposes the Buckets that are currently a part of
  //// the implementation.
  ////

  /// Does the given bucket belong to this partition.
  bool belongs(const Bucket &bkt) const { return bkt.getPartition() == this;}

  size_t num_buckets() const { return m_buckets.size();}

  inline std::vector<Bucket *>::iterator begin() { return m_buckets.begin(); }
  inline std::vector<Bucket *>::iterator end() { return m_buckets.end(); }

  inline std::vector<Bucket *>::const_iterator begin() const { return m_buckets.begin(); }
  inline std::vector<Bucket *>::const_iterator end() const { return m_buckets.end(); }

  ////
  //// The following are used internally and for unit testing.
  ////

  size_t compute_size()
  {
    size_t partition_size = 0;
    for (std::vector<Bucket *>::const_iterator b_i = begin(), b_e = end(); b_i != b_e; ++b_i)
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
  std::vector<Bucket *> m_buckets;

  // Number of entities in this partition.
  size_t m_size;

  bool m_updated_since_compress;

  bool m_updated_since_sort;

  //
  // Internal methods
  //

  void remove_impl(Entity entity, bool due_to_move=false);

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

  void internal_check_size_invariant() const;

  void internal_check_no_null_buckets_invariant() const;

  void internal_swap_to_end(Entity entity);
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

inline
void Partition::internal_check_size_invariant() const
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

inline
void Partition::internal_check_no_null_buckets_invariant() const
{
#ifndef NDEBUG
  for (int i = 0, e = m_buckets.size(); i < e; ++i) {
    ThrowAssert(m_buckets[i] != NULL);
  }
#endif
}

} // impl
} // mesh
} // stk

#endif /* PartitionFAMILY_HPP_ */
