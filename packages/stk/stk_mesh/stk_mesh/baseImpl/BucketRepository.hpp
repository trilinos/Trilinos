/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_BucketRepository_hpp
#define stk_mesh_BucketRepository_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Iterators.hpp>

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository;

class BucketRepository
{
public:
  ~BucketRepository();
  BucketRepository(
      BulkData & mesh,
      unsigned bucket_capacity,
      unsigned entity_rank_count,
      EntityRepository & entity_repo
      );

  /** \brief  Query to get all buckets of a given entity rank */
  const std::vector<Bucket*> & buckets( EntityRank rank ) const
  {
    ThrowAssertMsg( rank < m_buckets.size(), "Invalid entity rank " << rank );

    if (m_need_sync_from_partitions[rank])
    {
        const_cast<BucketRepository *>(this)->sync_from_partitions(rank);
    }

    return m_buckets[ rank ];
  }

  //------------------------------------
  /** \brief  Query the upper bound on the number of mesh entities
    *         that may be associated with a single bucket.
    */
  unsigned bucket_capacity() const { return m_bucket_capacity; }

  BulkData& mesh() const { return m_mesh; }

  //------------------------------------

  /** \brief  Rotate the field data of multistate fields.
   *
   *  <PRE>
   *  Rotation of states:
   *    StateN   <- StateNP1 (StateOld <- StateNew)
   *    StateNM1 <- StateN   (StateNM1 <- StateOld)
   *    StateNM2 <- StateNM1
   *    StateNM3 <- StateNM2
   *    StateNM3 <- StateNM2
   *  </PRE>
   */
  void update_field_data_states() const ;

  void copy_fields( Bucket & k_dst , unsigned i_dst ,
                           Bucket & k_src , unsigned i_src )
  { k_dst.replace_fields(i_dst,k_src,i_src); }

  void internal_sort_bucket_entities();

  void optimize_buckets();

  AllBucketsRange get_bucket_range() const
  {
    const_cast<BucketRepository *>(this)->sync_from_partitions();
    return stk::mesh::get_bucket_range(m_buckets);
  }

  AllBucketsRange get_bucket_range(EntityRank entity_rank) const
  {
    if (m_need_sync_from_partitions[entity_rank])
    {
        const_cast<BucketRepository *>(this)->sync_from_partitions(entity_rank);
    }
    std::vector< std::vector<Bucket*> >::const_iterator itr = m_buckets.begin() + entity_rank;
    return stk::mesh::get_bucket_range(m_buckets, itr);
  }

  ////
  //// Partitions are now the primary location of buckets.
  ////

  friend class Partition;

  Partition *get_or_create_partition(const unsigned arg_entity_rank ,
                                     const PartVector &parts);

  Partition *get_or_create_partition(const unsigned arg_entity_rank ,
                                     const OrdinalVector &parts);

  // Update m_buckets from the partitions.
  void sync_from_partitions();
  void sync_from_partitions(EntityRank rank);

  // Used in unit tests.  Returns the current partitions.
  std::vector<Partition *> get_partitions(EntityRank rank);

  // Used in unit tests. Delete the Partitions in m_partitions, clear it, and then (re-)build
  // the Partitions from the m_buckets.
  void sync_to_partitions();


private:
  BucketRepository();

  BulkData                            & m_mesh ; // Associated Bulk Data Aggregate
  unsigned                              m_bucket_capacity ; // Maximum number of entities per bucket

  // Vector of bucket pointers by rank.  This is now a cache and no longer the primary
  // location of Buckets when USE_STK_MESH_IMPL_PARTITION is #defined.
  std::vector< std::vector<Bucket*> >   m_buckets ;

  EntityRepository                    & m_entity_repo ;

  std::vector<std::vector<Partition *> >         m_partitions;
  std::vector<bool>                              m_need_sync_from_partitions;
};

} // namespace impl
} // namespace mesh
} // namespace stk


#endif // stk_mesh_BucketRepository_hpp
