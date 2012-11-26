/*
 * Partition.hpp
 *
 */

#ifndef STK_MESH_IMPL_PARTITION_HPP_
#define STK_MESH_IMPL_PARTITION_HPP_

#include <stk_mesh/base/Types.hpp>


namespace stk {
namespace mesh {

namespace impl {

class BucketRepository;

class Partition
{
    friend class BucketRepository;

public:

    Partition(BucketRepository *repo, EntityRank rank,
              const std::vector<PartOrdinal> &key);

    virtual ~Partition();

    ////
    //// The main part of the interface is Bucket-free.
    ////

    /// Rank of the entities in this partition.
    const EntityRank get_rank() const { return m_rank; }

    /// Is this partition empty.
    inline bool empty() const;

    /// Returns the representation used by BucketRepository to identify a bucket,
    /// including the parts it corresponds to.
    const std::vector<PartOrdinal> &get_legacy_partition_id() const
    {
        return m_extPartitionKey;
    }

    const unsigned * key() const { return &m_extPartitionKey[0]; }

    /// Add an entity to this partition.  The entity must not be a member
    /// of another partition.
    bool add(Entity entity);

    /// Move an entity from this partion to a destination partition.
    void move_to(Entity entity, Partition &dst_partition);

    /// Remove an entity from this partition.
    bool remove(Entity entity, bool not_in_move_to = false);

    /// Compress this partion into a single bucket of sorted Entities.
    void compress(bool force = false);

    /// Sort the entities in this partition by EntityKey without changing
    /// the number or sizes of buckets.
    void sort(bool force = false);
    
    void update_state() const;

    ////
    //// This part of the interface exposes the Buckets that are currently a part of
    //// the implementation.
    ////

    /// Does the given bucket belong to this partition.
    inline bool belongs(Bucket *bkt) const;

    size_t num_buckets() const
    {
        return (m_modifyingBucketSet ? m_buckets.size() : m_endBucketIndex - m_beginBucketIndex);
    }

    inline std::vector<Bucket *>::iterator begin();
    inline std::vector<Bucket *>::iterator end();

    inline std::vector<Bucket *>::const_iterator begin() const;
    inline std::vector<Bucket *>::const_iterator end() const;

    size_t size() const { return m_size; }

    ////
    //// The following are used internally and for unit testing.
    ////

    size_t compute_size()
    {
        size_t partition_size = 0;
        std::vector<Bucket *>::iterator buckets_end = end();
        for (std::vector<Bucket *>::iterator b_i = begin(); b_i != buckets_end; ++b_i)
        {
            partition_size += (*b_i)->size();
        }
        m_size = partition_size;

        return partition_size;
    }

    // Enables sidestepping compiler fussiness wrt overloaded operator<<(..).
    std::ostream &streamit(std::ostream &os) const;

    // Just for unit testing.  Remove after refactor.
    static BucketRepository &getRepository(stk::mesh::BulkData &mesh);

    // Just for unit testing. DOES NOT SYNC DATA.  Within each bucket locally reverse
    // the order of the entities.
    void reverseEntityOrderWithinBuckets();

private:

    BucketRepository *m_repository;

    EntityRank m_rank;

    // Identifies the partition, borrowing the representation from BucketRepository.
    std::vector<PartOrdinal> m_extPartitionKey;

    // Used if the set of buckets (not just bucket contents) are being modified.
    std::vector<Bucket *> m_buckets;

    size_t m_size;

    // Used when the vector of Bucket * in the BucketRepository is the representation
    // being used.
    unsigned m_beginBucketIndex;
    unsigned m_endBucketIndex;

    // Flag that the set of buckets, and not just their contents, is being modified.
    bool m_modifyingBucketSet;

    bool m_updated_since_compress;

    bool m_updated_since_sort;

    // The partition has no buckets, not even an empty one left after removing all its
    // entities.
    bool no_buckets() const;

    bool modifying_bucket_set() const { return m_modifyingBucketSet; }

    // Take control of this partition's buckets away from the BucketRepository.
    bool take_bucket_control();

    // Make sure that the last bucket has room for an entity to be added to it, adding
    // an empty bucket if necessary.
    Bucket *get_bucket_for_adds();

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
{
    return partition_key_less( lhs_partition->key() , rhs ); }

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
