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

class BucketRepository {
public:
  ~BucketRepository();
  BucketRepository(
      BulkData & mesh,
      unsigned bucket_capacity,
      unsigned entity_rank_count,
      EntityRepository & entity_repo
      );

  /** \brief  Query all buckets of a given entity rank */
  const std::vector<Bucket*> & buckets( EntityRank rank ) const ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void remove_entity( Bucket * , unsigned );

  //------------------------------------
  /** \brief  Query the upper bound on the number of mesh entities
    *         that may be associated with a single bucket.
    */
  unsigned bucket_capacity() const { return m_bucket_capacity; }


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

  // Destroy the last empty bucket in a family:
  void destroy_bucket( const unsigned & entity_rank , Bucket * last );
  void destroy_bucket( Bucket * bucket );
  void declare_nil_bucket();
  Bucket * get_nil_bucket() const { return m_nil_bucket; }
  Bucket * declare_bucket(
      const unsigned entity_rank ,
      const unsigned part_count ,
      const unsigned part_ord[] ,
      const std::vector< FieldBase * > & field_set
      );
  void copy_fields( Bucket & k_dst , unsigned i_dst ,
                           Bucket & k_src , unsigned i_src );
  void zero_fields( Bucket & k_dst , unsigned i_dst );

  void internal_sort_bucket_entities();

  void add_entity_to_bucket( Entity & entity, Bucket & bucket )
  {
    bucket.m_bucketImpl.replace_entity( bucket.size() , & entity ) ;
    bucket.m_bucketImpl.increment_size();
  }

  void internal_propagate_relocation( Entity & );

  AllBucketsRange get_bucket_range() const
  {
    return stk::mesh::get_bucket_range(m_buckets);
  }

  AllBucketsRange get_bucket_range(EntityRank entity_rank) const
  {
    std::vector< std::vector<Bucket*> >::const_iterator itr = m_buckets.begin() + entity_rank;
    return stk::mesh::get_bucket_range(m_buckets, itr);
  }

private:
  BucketRepository();

  BulkData                            & m_mesh ; // Associated Bulk Data Aggregate
  unsigned                              m_bucket_capacity ; // Maximum number of entities per bucket
  std::vector< std::vector<Bucket*> >   m_buckets ; // Vector of bucket pointers by rank
  Bucket                              * m_nil_bucket ; // nil bucket

  EntityRepository                    & m_entity_repo ;
};



} // namespace impl
} // namespace mesh
} // namespace stk


#endif // stk_mesh_BucketRepository_hpp
