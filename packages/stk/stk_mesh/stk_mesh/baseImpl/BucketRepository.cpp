/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <cstdlib>
#include <stdexcept>

#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------


BucketRepository::BucketRepository(
    BulkData & mesh,
    unsigned bucket_capacity,
    unsigned entity_rank_count,
    EntityRepository & entity_repo
    )
  :m_mesh(mesh),
   m_bucket_capacity(bucket_capacity),
   m_buckets(entity_rank_count),
   m_nil_bucket(NULL),
   m_entity_repo(entity_repo)
{
}


BucketRepository::~BucketRepository()
{
  // Destroy buckets, which were *not* allocated by the set.

  try {
    for ( std::vector< std::vector<Bucket*> >::iterator
          i = m_buckets.end() ; i != m_buckets.begin() ; ) {
      try {
        std::vector<Bucket*> & kset = *--i ;

        while ( ! kset.empty() ) {
          try { destroy_bucket( kset.back() ); } catch(...) {}
          kset.pop_back();
        }
        kset.clear();
      } catch(...) {}
    }
    m_buckets.clear();
  } catch(...) {}

  try { if ( m_nil_bucket ) destroy_bucket( m_nil_bucket ); } catch(...) {}
}


//----------------------------------------------------------------------
// The current 'last' bucket in a family is to be deleted.
// The previous 'last' bucket becomes the new 'last' bucket in the family.

void BucketRepository::destroy_bucket( const unsigned & entity_rank , Bucket * bucket_to_be_deleted )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::destroy_bucket", LOG_BUCKET, bucket_to_be_deleted);

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(entity_rank),
                  "Entity rank " << entity_rank << " is invalid");

  std::vector<Bucket *> & bucket_set = m_buckets[entity_rank];

  // Get the first bucket in the same family as the bucket being deleted
  Bucket * const first = bucket_to_be_deleted->m_bucketImpl.first_bucket_in_family();

  ThrowRequireMsg( bucket_to_be_deleted->equivalent(*first), "Logic error - bucket_to_be_deleted is not in same family as first_bucket_in_family");
  ThrowRequireMsg( first->equivalent(*bucket_to_be_deleted), "Logic error - first_bucket_in_family is not in same family as bucket_to_be_deleted");

  ThrowRequireMsg( bucket_to_be_deleted->size() == 0,
      "Destroying non-empty bucket " << *(bucket_to_be_deleted->key()) );

  ThrowRequireMsg( bucket_to_be_deleted == first->m_bucketImpl.get_bucket_family_pointer(),
                   "Destroying bucket family") ;

  std::vector<Bucket*>::iterator ik = lower_bound(bucket_set, bucket_to_be_deleted->key());
  ThrowRequireMsg( ik != bucket_set.end() && bucket_to_be_deleted == *ik,
      "Bucket not found in bucket set for entity rank " << entity_rank );

  ik = bucket_set.erase( ik );

  if ( first != bucket_to_be_deleted ) {

    ThrowRequireMsg( ik != bucket_set.begin(),
                     "Where did first bucket go?" );

    first->m_bucketImpl.set_last_bucket_in_family( *--ik );

    ThrowRequireMsg ( first->m_bucketImpl.get_bucket_family_pointer()->size() != 0,
                      "TODO: Explain" );
  }

  destroy_bucket( bucket_to_be_deleted );
}

//----------------------------------------------------------------------
void BucketRepository::destroy_bucket( Bucket * bucket )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::destroy_bucket", LOG_BUCKET, bucket);

  delete bucket;
}

//
//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.
void
BucketRepository::declare_nil_bucket()
{
  TraceIf("stk::mesh::impl::BucketRepository::declare_nil_bucket", LOG_BUCKET);

  if (m_nil_bucket == NULL) {
    // Key layout:
    // { part_count + 1 , { part_ordinals } , family_count }

    std::vector<unsigned> new_key(2);
    new_key[0] = 1 ; // part_count + 1
    new_key[1] = 0 ; // family_count

    Bucket * bucket =
      new Bucket(m_mesh, InvalidEntityRank, new_key, 0);

    bucket->m_bucketImpl.set_bucket_family_pointer( bucket );

    //----------------------------------

    m_nil_bucket = bucket;
  }
}


/** 11/9/10 Discussion between Kendall, Alan, Todd:
 *  Kendall is confused about why presto would run faster simply by removing
 *  several fields that are not even used.  We considered this and posed the
 *  following possibility.  The current bucket allocation system guarantees
 *  that all the fields for a bucket are layed out contiguously in memory so
 *  that they can be accessed in a fast cache-friendly manner.  This also
 *  guarantees means that if a field is allocated but not used, it will still
 *  be chopped up and carried around in the bucket field data as part of the
 *  contiguous block of memory and that it will have to be skipped over as the
 *  computations progress over that block of data.  This would result in cache
 *  misses and reduced performance.  When they're removed, it makes sense that
 *  the performance might get better.
 *
 *  This leads to the idea that maybe we should test this in a use-case or
 *  performance test case and that we should include this in the performance
 *  comparison of the up-and-coming pluggable data module for the Bucket memory
 *  allocation.
 *
 *  It may be that a flat-array style data allocation for field data would
 *  eliminate this issue.
 **/

//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.
Bucket *
BucketRepository::declare_bucket(
                        const unsigned arg_entity_rank ,
                        const unsigned part_count ,
                        const unsigned part_ord[] ,
                        const std::vector< FieldBase * > & field_set
                              )
{
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  TraceIf("stk::mesh::impl::BucketRepository::declare_bucket", LOG_BUCKET);

  const unsigned max = static_cast<unsigned>(-1);

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
                  "Entity rank " << arg_entity_rank << " is invalid");

  ThrowRequireMsg( !m_buckets.empty(),
    "m_buckets is empty! Did you forget to initialize MetaData before creating BulkData?");
  std::vector<Bucket *> & bucket_set = m_buckets[ arg_entity_rank ];


  std::vector<unsigned> key(2+part_count) ;

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , family_count }
  // Thus family_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key.

  key[0] = part_count+1;
  key[ key[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i ) { key[i+1] = part_ord[i] ; }
  }

  //----------------------------------
  // Bucket family has all of the same parts.
  // Look for the last bucket in this family:

  const std::vector<Bucket*>::iterator ik = lower_bound( bucket_set , &key[0] );

  //----------------------------------
  // If a member of the bucket family has space, it is the last one
  // since buckets are kept packed.
  const bool bucket_family_exists =
    ik != bucket_set.begin() && bucket_part_equal( ik[-1]->key() , &key[0] );

  Bucket * const last_bucket = bucket_family_exists ? ik[-1] : NULL ;

  Bucket          * bucket    = NULL ;

  if ( last_bucket == NULL ) { // First bucket in this family
    key[ key[0] ] = 0 ; // Set the key's family count to zero
  }
  else { // Last bucket present, can it hold one more entity?

    ThrowRequireMsg( last_bucket->size() != 0,
                     "Last bucket should not be empty.");

    //field_map = last_bucket->m_bucketImpl.get_field_map();

    const unsigned last_count = last_bucket->key()[ key[0] ];

    const unsigned cap = last_bucket->capacity();

    if ( last_bucket->size() < cap ) {
      bucket = last_bucket ;
    }
    else if ( last_count < max ) {
      key[ key[0] ] = 1 + last_count ; // Increment the key's family count.
    }
    else {
      // ERROR insane number of buckets!
      ThrowRequireMsg( false, "Insanely large number of buckets" );
    }
  }


  //----------------------------------

  //Required bucket does not exist
  if ( NULL == bucket )
  {
    bucket = new Bucket( m_mesh, arg_entity_rank, key, m_bucket_capacity);

    Bucket * first_bucket = last_bucket ? last_bucket->m_bucketImpl.first_bucket_in_family() : bucket ;

    bucket->m_bucketImpl.set_first_bucket_in_family(first_bucket); // Family members point to first bucket

    first_bucket->m_bucketImpl.set_last_bucket_in_family(bucket); // First bucket points to new last bucket

    bucket_set.insert( ik , bucket );
  }

  //----------------------------------

  ThrowRequireMsg( bucket->equivalent(*bucket->m_bucketImpl.first_bucket_in_family()), "Logic error - new bucket is not in same family as first_bucket_in_family");
  ThrowRequireMsg( bucket->m_bucketImpl.first_bucket_in_family()->equivalent(*bucket), "Logic error - first_bucket_in_family is not in same family as new bucket");

  return bucket ;
}

//----------------------------------------------------------------------

void BucketRepository::initialize_fields( Bucket & k_dst , unsigned i_dst )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::initialize_fields", LOG_BUCKET, &k_dst);
  k_dst.m_bucketImpl.initialize_fields(i_dst);
}

//----------------------------------------------------------------------

void BucketRepository::update_field_data_states() const
{
  TraceIf("stk::mesh::impl::BucketRepository::update_field_data_states", LOG_BUCKET);

  for ( std::vector< std::vector<Bucket*> >::const_iterator
        i = m_buckets.begin() ; i != m_buckets.end() ; ++i ) {

    const std::vector<Bucket*> & kset = *i ;

    for ( std::vector<Bucket*>::const_iterator
          ik = kset.begin() ; ik != kset.end() ; ++ik ) {
      (*ik)->m_bucketImpl.update_state();
    }
  }
}


//----------------------------------------------------------------------


void BucketRepository::internal_sort_bucket_entities()
{
  TraceIf("stk::mesh::impl::BucketRepository::internal_sort_bucket_entities", LOG_BUCKET);

  for ( EntityRank entity_rank = 0 ;
        entity_rank < m_buckets.size() ; ++entity_rank ) {

    std::vector<Bucket*> & buckets = m_buckets[ entity_rank ];

    size_t bk = 0 ; // Offset to first bucket of the family
    size_t ek = 0 ; // Offset to end   bucket of the family

    for ( ; bk < buckets.size() ; bk = ek ) {
      Bucket * b_scratch = NULL ;
      Bucket * ik_vacant = buckets[bk]->m_bucketImpl.last_bucket_in_family();
      unsigned ie_vacant = ik_vacant->size();

      if ( ik_vacant->capacity() <= ie_vacant ) {
        // Have to create a bucket just for the scratch space...
        const unsigned * const bucket_key = buckets[bk]->key() ;
        const unsigned         part_count = bucket_key[0] - 1 ;
        const unsigned * const part_ord   = bucket_key + 1 ;

        b_scratch = declare_bucket( entity_rank ,
            part_count , part_ord ,
            MetaData::get(m_mesh).get_fields() );

        ik_vacant = b_scratch ;
        ie_vacant = 0 ;
      }

      ik_vacant->m_bucketImpl.replace_entity( ie_vacant , NULL ) ;

      // Determine offset to the end bucket in this family:
      while ( ek < buckets.size() && ik_vacant != buckets[ek] ) { ++ek ; }
      if (ek < buckets.size()) ++ek ;

      unsigned count = 0 ;
      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        count += buckets[ik]->size();
      }

      std::vector<Entity*> entities( count );

      std::vector<Entity*>::iterator j = entities.begin();

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          *j = & b[i] ;
        }
      }

      std::sort( entities.begin() , entities.end() , EntityLess() );

      j = entities.begin();

      bool change_this_family = false ;

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          Entity * const current = & b[i] ;

          if ( current != *j ) {

            if ( current ) {
              // Move current entity to the vacant spot
              copy_fields( *ik_vacant , ie_vacant , b, i );
              m_entity_repo.change_entity_bucket(*ik_vacant, *current, ie_vacant);
              ik_vacant->m_bucketImpl.replace_entity( ie_vacant , current ) ;
            }

            // Set the vacant spot to where the required entity is now.
            ik_vacant = & ((*j)->bucket()) ;
            ie_vacant = (*j)->bucket_ordinal() ;
            ik_vacant->m_bucketImpl.replace_entity( ie_vacant , NULL ) ;

            // Move required entity to the required spot
            copy_fields( b, i, *ik_vacant , ie_vacant );
            m_entity_repo.change_entity_bucket( b, **j, i);
            b.m_bucketImpl.replace_entity( i, *j );

            change_this_family = true ;
          }

          // Once a change has occured then need to propagate the
          // relocation for the remainder of the family.
          // This allows the propagation to be performed once per
          // entity as opposed to both times the entity is moved.

          if ( change_this_family ) { internal_propagate_relocation( **j ); }
        }
      }

      if ( b_scratch ) {
        // Created a last bucket, now have to destroy it.
        destroy_bucket( entity_rank , b_scratch );
        --ek ;
      }
    }
  }
}

void BucketRepository::optimize_buckets()
{
  TraceIf("stk::mesh::impl::BucketRepository::optimize_buckets", LOG_BUCKET);

  for ( EntityRank entity_rank = 0 ;
      entity_rank < m_buckets.size() ; ++entity_rank )
  {

    std::vector<Bucket*> & buckets = m_buckets[ entity_rank ];

    std::vector<Bucket*> tmp_buckets;

    size_t begin_family = 0 ; // Offset to first bucket of the family
    size_t end_family = 0 ; // Offset to end   bucket of the family

    //loop over families
    for ( ; begin_family < buckets.size() ; begin_family = end_family ) {
      Bucket * last_bucket_in_family  = buckets[begin_family]->m_bucketImpl.last_bucket_in_family();

      // Determine offset to the end bucket in this family:
      while ( end_family < buckets.size() && last_bucket_in_family != buckets[end_family] ) { ++end_family ; }
      if (end_family < buckets.size())  ++end_family ; //increment past the end

      //if compressed and sorted go to the next family
      const bool is_compressed =    (end_family-begin_family == 1)
                                 && (buckets[begin_family]->size() == buckets[begin_family]->capacity());
      if (is_compressed) {
        const Bucket & b = *buckets[begin_family];
        bool is_sorted = true;
        for (size_t i=0, end=b.size()-1; i<end && is_sorted; ++i)
        {
          if(b[i].key() >= b[i+1].key()) is_sorted = false;
        }
        if (is_sorted) {
          tmp_buckets.push_back(buckets[begin_family]);
          continue;
        }
      }

      std::vector<unsigned> new_key = buckets[begin_family]->m_bucketImpl.key_vector();
      //index of bucket in family
      new_key[ new_key[0] ] = 0;

      unsigned new_capacity = 0 ;
      for ( size_t i = begin_family ; i != end_family ; ++i ) {
        new_capacity += buckets[i]->m_bucketImpl.size();
      }

      std::vector<Entity*> entities;
      entities.reserve(new_capacity);

      for ( size_t i = begin_family ; i != end_family ; ++i ) {
        Bucket& b = *buckets[i];
        for(size_t j=0; j<b.size(); ++j) {
          entities.push_back(&b[j]);
        }
      }

      std::sort( entities.begin(), entities.end(), EntityLess() );

      Bucket * new_bucket = new Bucket( m_mesh,
          entity_rank,
          new_key,
          new_capacity
          );

      new_bucket->m_bucketImpl.set_first_bucket_in_family(new_bucket); // Family members point to first bucket
      new_bucket->m_bucketImpl.set_last_bucket_in_family(new_bucket); // First bucket points to new last bucket

      tmp_buckets.push_back(new_bucket);

      for(size_t new_ordinal=0; new_ordinal<entities.size(); ++new_ordinal) {
        //increase size of the new_bucket
        new_bucket->m_bucketImpl.increment_size();

        Entity & entity = *entities[new_ordinal];
        Bucket& old_bucket = entity.bucket();
        unsigned old_ordinal = entity.bucket_ordinal();

        //copy field data from old to new
        copy_fields( *new_bucket, new_ordinal, old_bucket, old_ordinal);
        m_entity_repo.change_entity_bucket( *new_bucket, entity, new_ordinal);
        new_bucket->m_bucketImpl.replace_entity( new_ordinal , &entity ) ;
        internal_propagate_relocation(entity);
      }

      for (size_t ik = begin_family; ik != end_family; ++ik) {
        delete buckets[ik];
        buckets[ik] = NULL;
      }
    }

    buckets.swap(tmp_buckets);
  }
}
//----------------------------------------------------------------------

void BucketRepository::remove_entity( Bucket * k , unsigned i )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::remove_entity", LOG_BUCKET, k);

  ThrowRequireMsg( k != m_nil_bucket, "Cannot remove entity from nil_bucket" );

  const EntityRank entity_rank = k->entity_rank();

  // Last bucket in the family of buckets with the same parts.
  // The last bucket is the only non-full bucket in the family.

  Bucket * const last = k->m_bucketImpl.last_bucket_in_family();

  ThrowRequireMsg( last->equivalent(*k), "Logic error - last bucket in family not equivalent to bucket");
  ThrowRequireMsg( k->equivalent(*last), "Logic error - bucket not equivalent to last bucket in family");

  // Fill in the gap if it is not the last entity being removed

  if ( last != k || k->size() != i + 1 ) {

    // Copy last entity in last bucket to bucket *k slot i

    Entity & entity = (*last)[ last->size() - 1 ];

    copy_fields( *k , i , *last , last->size() - 1 );

    k->m_bucketImpl.replace_entity(i, & entity ) ;
    m_entity_repo.change_entity_bucket( *k, entity, i);

    // Entity field data has relocated

    internal_propagate_relocation( entity );
  }

  last->m_bucketImpl.decrement_size();

  last->m_bucketImpl.replace_entity( last->size() , NULL ) ;

  if ( 0 == last->size() ) {
    destroy_bucket( entity_rank , last );
  }
}

//----------------------------------------------------------------------

void BucketRepository::internal_propagate_relocation( Entity & entity )
{
  TraceIf("stk::mesh::impl::BucketRepository::internal_propagate_relocation", LOG_BUCKET);

  const EntityRank erank = entity.entity_rank();
  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const EntityRank rel_rank = rel->entity_rank();
    if ( rel_rank < erank ) {
      Entity & e_to = * rel->entity();

      set_field_relations( entity, e_to, rel->identifier() );
    }
    else if ( erank < rel_rank ) {
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel->identifier() );
    }
  }
}


} // namespace impl
} // namespace mesh
} // namespace stk


