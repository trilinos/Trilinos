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

    for ( std::vector<std::vector<Partition *> >::iterator pv_i = m_partitions.begin();
            pv_i != m_partitions.end(); ++pv_i)
    {
        for (std::vector<Partition *>::iterator p_j = pv_i->begin();
                p_j != pv_i->end(); ++p_j)
        {
            delete *p_j;
        }
        pv_i->clear();
    }
    m_partitions.clear();

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
// The current 'last' bucket in a partition is to be deleted.
// The previous 'last' bucket becomes the new 'last' bucket in the partition.

void BucketRepository::destroy_bucket( const unsigned & entity_rank , Bucket * bucket_to_be_deleted )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::destroy_bucket", LOG_BUCKET, bucket_to_be_deleted);

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(entity_rank),
                  "Entity rank " << entity_rank << " is invalid");

  std::vector<Bucket *> & bucket_set = m_buckets[entity_rank];

  // Get the first bucket in the same partition as the bucket being deleted
  Bucket * const first = bucket_to_be_deleted->first_bucket_in_partition();

  ThrowRequireMsg( bucket_to_be_deleted->in_same_partition(*first), "Logic error - bucket_to_be_deleted is not in same partition as first bucket in partition");
  ThrowRequireMsg( first->in_same_partition(*bucket_to_be_deleted), "Logic error - first bucket in partition is not in same partition as bucket_to_be_deleted");

  ThrowRequireMsg( bucket_to_be_deleted->size() == 0,
      "Destroying non-empty bucket " << *(bucket_to_be_deleted->key()) );

  ThrowRequireMsg( bucket_to_be_deleted == first->get_partition_pointer(),
                   "Destroying partition") ;

  std::vector<Bucket*>::iterator ik = lower_bound(bucket_set, bucket_to_be_deleted->key());
  ThrowRequireMsg( ik != bucket_set.end() && bucket_to_be_deleted == *ik,
      "Bucket not found in bucket set for entity rank " << entity_rank );

  ik = bucket_set.erase( ik );

  if ( first != bucket_to_be_deleted ) {

    ThrowRequireMsg( ik != bucket_set.begin(),
                     "Where did first bucket go?" );

    first->set_last_bucket_in_partition( *--ik );

    ThrowRequireMsg ( first->get_partition_pointer()->size() != 0,
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
    // { part_count + 1 , { part_ordinals } , partition_count }

    std::vector<unsigned> new_key(2);
    new_key[0] = 1 ; // part_count + 1
    new_key[1] = 0 ; // partition_count

    Bucket * bucket =
      new Bucket(m_mesh, InvalidEntityRank, new_key, 0);

    bucket->set_partition_pointer( bucket );

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
  // { part_count + 1 , { part_ordinals } , partition_count }
  // Thus partition_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key.

  key[0] = part_count+1;
  key[ key[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i ) { key[i+1] = part_ord[i] ; }
  }

  //----------------------------------
  // partition has all of the same parts.
  // Look for the last bucket in this partition:

  const std::vector<Bucket*>::iterator ik = lower_bound( bucket_set , &key[0] );

  //----------------------------------
  // If a member of the partition has space, it is the last one
  // since buckets are kept packed.
  const bool partition_exists =
    ik != bucket_set.begin() && bucket_part_equal( ik[-1]->key() , &key[0] );

  Bucket * const last_bucket = partition_exists ? ik[-1] : NULL ;

  Bucket          * bucket    = NULL ;

  if ( last_bucket == NULL ) { // First bucket in this partition
    key[ key[0] ] = 0 ; // Set the partition key's bucket count to zero
  }
  else { // Last bucket present, can it hold one more entity?

    ThrowRequireMsg( last_bucket->size() != 0,
                     "Last bucket should not be empty.");

    //field_map = last_bucket->get_field_map();

    const unsigned last_count = last_bucket->key()[ key[0] ];

    const unsigned cap = last_bucket->capacity();

    if ( last_bucket->size() < cap ) {
      bucket = last_bucket ;
    }
    else if ( last_count < max ) {
      key[ key[0] ] = 1 + last_count ; // Increment the partition key's bucket count.
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

    Bucket * first_bucket = last_bucket ? last_bucket->first_bucket_in_partition() : bucket ;

    bucket->set_first_bucket_in_partition(first_bucket); // partition members point to first bucket

    first_bucket->set_last_bucket_in_partition(bucket); // First bucket points to new last bucket

    bucket_set.insert( ik , bucket );
  }

  //----------------------------------

  ThrowRequireMsg( bucket->in_same_partition(*bucket->first_bucket_in_partition()), "Logic error - new bucket is not in same partition as first bucket in partition");
  ThrowRequireMsg( bucket->first_bucket_in_partition()->in_same_partition(*bucket), "Logic error - first bucket in partition is not in same partition as new bucket");

  return bucket ;
}

//----------------------------------------------------------------------

void BucketRepository::initialize_fields( Bucket & k_dst , unsigned i_dst )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::initialize_fields", LOG_BUCKET, &k_dst);
  k_dst.initialize_fields(i_dst);
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
      (*ik)->update_state();
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

    size_t bk = 0 ; // Offset to first bucket of the partition
    size_t ek = 0 ; // Offset to end   bucket of the partition

    for ( ; bk < buckets.size() ; bk = ek ) {
      Bucket * b_scratch = NULL ;
      Bucket * ik_vacant = buckets[bk]->last_bucket_in_partition();
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

      ik_vacant->replace_entity( ie_vacant , Entity() ) ;

      // Determine offset to the end bucket in this partition:
      while ( ek < buckets.size() && ik_vacant != buckets[ek] ) { ++ek ; }
      if (ek < buckets.size()) ++ek ;

      unsigned count = 0 ;
      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        count += buckets[ik]->size();
      }

      std::vector<Entity> entities( count );

      std::vector<Entity>::iterator j = entities.begin();

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          *j = b[i];
        }
      }

      std::sort( entities.begin() , entities.end() , EntityLess() );

      j = entities.begin();

      bool change_this_partition = false ;

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          Entity const current = b[i];

          if ( current != *j ) {

            if ( current.is_valid() ) {
              // Move current entity to the vacant spot
              copy_fields( *ik_vacant , ie_vacant , b, i );
              m_entity_repo.change_entity_bucket(*ik_vacant, current, ie_vacant);
              ik_vacant->replace_entity( ie_vacant , current ) ;
            }

            // Set the vacant spot to where the required entity is now.
            ik_vacant = & (j->bucket()) ;
            ie_vacant = j->bucket_ordinal() ;
            ik_vacant->replace_entity( ie_vacant , Entity() ) ;

            // Move required entity to the required spot
            copy_fields( b, i, *ik_vacant , ie_vacant );
            m_entity_repo.change_entity_bucket( b, *j, i);
            b.replace_entity( i, *j );

            change_this_partition = true ;
          }

          // Once a change has occured then need to propagate the
          // relocation for the remainder of the partition.
          // This allows the propagation to be performed once per
          // entity as opposed to both times the entity is moved.

          if ( change_this_partition ) { internal_propagate_relocation( *j ); }
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

    size_t begin_partition = 0 ; // Offset to first bucket of the partition
    size_t end_partition = 0 ; // Offset to end   bucket of the partition

    //loop over families
    for ( ; begin_partition < buckets.size() ; begin_partition = end_partition ) {
      Bucket * last_bucket_in_partition  = buckets[begin_partition]->last_bucket_in_partition();

      // Determine offset to the end bucket in this partition:
      while ( end_partition < buckets.size() && last_bucket_in_partition != buckets[end_partition] ) { ++end_partition ; }
      if (end_partition < buckets.size())  ++end_partition ; //increment past the end

      //only one bucket in the partition
      //go to the next partition
      if (end_partition - begin_partition == 1) {
        tmp_buckets.push_back(buckets[begin_partition]);
        continue;
      }

      std::vector<unsigned> new_key = buckets[begin_partition]->key_vector();
      //index of bucket in partition
      new_key[ new_key[0] ] = 0;

      unsigned new_capacity = 0 ;
      for ( size_t i = begin_partition ; i != end_partition ; ++i ) {
        new_capacity += buckets[i]->capacity();
      }

      std::vector<Entity> entities;
      entities.reserve(new_capacity);

      for ( size_t i = begin_partition ; i != end_partition ; ++i ) {
        Bucket& b = *buckets[i];
        for(size_t j=0; j<b.size(); ++j) {
          entities.push_back(b[j]);
        }
      }

      std::sort( entities.begin(), entities.end(), EntityLess() );

      Bucket * new_bucket = new Bucket( m_mesh,
          entity_rank,
          new_key,
          new_capacity
          );

      new_bucket->set_first_bucket_in_partition(new_bucket); // partition members point to first bucket
      new_bucket->set_last_bucket_in_partition(new_bucket); // First bucket points to new last bucket

      tmp_buckets.push_back(new_bucket);

      for(size_t new_ordinal=0; new_ordinal<entities.size(); ++new_ordinal) {
        //increase size of the new_bucket
        new_bucket->increment_size();

        Entity entity = entities[new_ordinal];
        Bucket& old_bucket = entity.bucket();
        unsigned old_ordinal = entity.bucket_ordinal();

        //copy field data from old to new
        copy_fields( *new_bucket, new_ordinal, old_bucket, old_ordinal);
        m_entity_repo.change_entity_bucket( *new_bucket, entity, new_ordinal);
        new_bucket->replace_entity( new_ordinal , entity ) ;
        internal_propagate_relocation(entity);
      }

      for (size_t ik = begin_partition; ik != end_partition; ++ik) {
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

  // Last bucket in the partition of buckets with the same parts.
  // The last bucket is the only non-full bucket in the partition.

  Bucket * const last = k->last_bucket_in_partition();

  ThrowRequireMsg( last->in_same_partition(*k), "Logic error - last bucket's partition is not bucket's");
  ThrowRequireMsg( k->in_same_partition(*last), "Logic error - bucket's partition is not last bucket's");

  // Fill in the gap if it is not the last entity being removed

  if ( last != k || k->size() != i + 1 ) {

    // Copy last entity in last bucket to bucket *k slot i

    Entity entity = (*last)[ last->size() - 1 ];

    copy_fields( *k , i , *last , last->size() - 1 );

    k->replace_entity(i, entity ) ;
    m_entity_repo.change_entity_bucket( *k, entity, i);

    // Entity field data has relocated

    internal_propagate_relocation( entity );
  }

  last->decrement_size();

  last->replace_entity( last->size() , Entity() ) ;

  if ( 0 == last->size() ) {
    destroy_bucket( entity_rank , last );
  }
}

//----------------------------------------------------------------------

void BucketRepository::internal_propagate_relocation( Entity entity )
{
  TraceIf("stk::mesh::impl::BucketRepository::internal_propagate_relocation", LOG_BUCKET);

  const EntityRank erank = entity.entity_rank();
  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const EntityRank rel_rank = rel->entity_rank();
    if ( rel_rank < erank ) {
      Entity e_to = rel->entity();

      set_field_relations( entity, e_to, rel->identifier() );
    }
    else if ( erank < rel_rank ) {
      Entity e_from = rel->entity();

      set_field_relations( e_from, entity, rel->identifier() );
    }
  }
}

void BucketRepository::sync_to_partitions()
{
    for ( std::vector<std::vector<Partition *> >::iterator pv_i = m_partitions.begin();
            pv_i != m_partitions.end(); ++pv_i)
    {
        for (std::vector<Partition *>::iterator p_j = pv_i->begin();
                p_j != pv_i->end(); ++p_j)
        {
            delete *p_j;
        }
        pv_i->clear();
    }
    m_partitions.clear();
    m_partitions.resize(4);

    for ( EntityRank entity_rank = 0 ;
        entity_rank < m_buckets.size() ; ++entity_rank )
    {
      std::vector<stk::mesh::impl::Partition *> &partitions = m_partitions[entity_rank];
      std::vector<Bucket*> & buckets = m_buckets[ entity_rank ];

      size_t begin_partition = 0 ; // Offset to first bucket of the partition
      size_t end_partition = 0 ; // Offset to end   bucket of the partition

      //loop over families (of entity_rank) via buckets
      for ( ; begin_partition < buckets.size() ; begin_partition = end_partition )
      {
        partitions.push_back(new Partition(this, entity_rank));
        Bucket * last_bucket_in_partition = buckets[begin_partition]->last_bucket_in_partition();

        // Determine offset to the end bucket in this partition:
        while ((end_partition < buckets.size()) && (last_bucket_in_partition != buckets[end_partition]))
        {
            ++end_partition ;
        }
        if (end_partition < buckets.size())
        {
            ++end_partition ; //increment past the last bucket in the partition.
        }
        Partition *partition = partitions.back();
        partition->m_extPartitionKey = buckets[begin_partition]->key_vector();
        // partition.m_stkPartition.pop_back();
        partition->m_beginBucketIndex = static_cast<unsigned>(begin_partition);
        partition->m_endBucketIndex = static_cast<unsigned>(end_partition);
      }

      // Let all the buckets of entity_rank know about their bucket families.
      size_t num_partitions = partitions.size();
      for (size_t i = 0; i < num_partitions; ++i)
      {
          stk::mesh::impl::Partition *partition = partitions[i];
          std::vector<stk::mesh::Bucket *>::iterator bkt= partition->begin();
          for (;bkt != partition->end(); ++bkt)
          {
              (*bkt)->m_partition = partitions[i];
          }
          // std::cout << partition << std::endl;
      }
    }
}

void BucketRepository::sync_from_partitions()
{
    for (size_t rank = 0; rank < m_partitions.size(); ++rank)
    {
        bool need_sync = false;

        std::vector<Partition *> &partitions = m_partitions[rank];
        size_t num_partitions = partitions.size();
        for (size_t p_i = 0; p_i < num_partitions; ++p_i)
        {
            if (partitions[p_i]->need_sync_to_repository())
            {
                need_sync = true;
            }
        }

        if (need_sync)
        {
            size_t num_buckets = 0;
            for (size_t p_i = 0; p_i < num_partitions; ++p_i)
            {
                partitions[p_i]->take_bucket_control();
                num_buckets += partitions[p_i]->num_buckets();
            }
            if (!num_buckets)
            {
                return;
            }
            m_buckets[rank].resize(num_buckets);

            std::vector<Bucket *>::iterator bkts_i = m_buckets[rank].begin();
            for (size_t p_i = 0; p_i < num_partitions; ++p_i)
            {
                size_t num_bkts_in_partition = partitions[p_i]->num_buckets();
                std::copy(partitions[p_i]->begin(), partitions[p_i]->end(),
                          bkts_i);
                bkts_i += num_bkts_in_partition;
            }
        }
    }
}


std::vector<Partition *> BucketRepository::get_partitions(EntityRank rank)
{
    if (m_mesh.synchronized_state() != BulkData::SYNCHRONIZED)
    {
        std::vector<Partition *>();
    }
    std::vector<Partition *> retval;
    std::vector<Partition *> &bf_vec = m_partitions[rank];
    for (size_t i = 0; i < bf_vec.size(); ++i)
    {
        retval.push_back(bf_vec[i]);
    }
    return retval;
}

} // namespace impl
} // namespace mesh
} // namespace stk


