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
#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {
namespace impl {

BucketRepository::BucketRepository(
    BulkData & mesh,
    unsigned bucket_capacity,
    unsigned entity_rank_count,
    EntityRepository & entity_repo
    )
  :m_mesh(mesh),
   m_bucket_capacity(bucket_capacity),
   m_buckets(entity_rank_count),
   m_entity_repo(entity_repo),
   m_partitions(entity_rank_count),
   m_need_sync_from_partitions(entity_rank_count, false)
{
    // Nada.
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
    m_buckets.clear();
  } catch(...) {}
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

void BucketRepository::destroy_bucket( Bucket * bucket )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::destroy_bucket", LOG_BUCKET, bucket);

  delete bucket;
}


void BucketRepository::initialize_fields( Bucket & k_dst , unsigned i_dst )
{
  TraceIfWatching("stk::mesh::impl::BucketRepository::initialize_fields", LOG_BUCKET, &k_dst);
  k_dst.initialize_fields(i_dst);
}


void BucketRepository::update_field_data_states() const
{
  TraceIf("stk::mesh::impl::BucketRepository::update_field_data_states", LOG_BUCKET);

  for (std::vector<std::vector<Partition *> >::const_iterator
        i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
  {
      const std::vector<Partition *> & pset = *i ;

      for ( std::vector<Partition*>::const_iterator
                  ip = pset.begin() ; ip != pset.end() ; ++ip )
      {
          (*ip)->update_state();
      }
  }
}


void BucketRepository::internal_sort_bucket_entities()
{
    for (std::vector<std::vector<Partition *> >::const_iterator
          i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
    {
        const std::vector<Partition *> & pset = *i ;
        for ( std::vector<Partition*>::const_iterator
                    ip = pset.begin() ; ip != pset.end() ; ++ip )
        {
            (*ip)->sort();
        }
    }
}


void BucketRepository::optimize_buckets()
{
    for (std::vector<std::vector<Partition *> >::const_iterator
          i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
    {
        const std::vector<Partition *> & pset = *i ;
        for ( std::vector<Partition*>::const_iterator
                    ip = pset.begin() ; ip != pset.end() ; ++ip )
        {
            (*ip)->compress();
        }
    }
}


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

////
//// Note that in both versions of get_or_create_partition(..) we need to construct a
//// key vector that the particular format so we can use the lower_bound(..) function to
//// lookup the partition.  Because we are using partitions now instead of buckets, it
//// should be possible to do without that vector and instead do the lookup directly from
//// the PartVector or OrdinalVector.
////

Partition *BucketRepository::get_or_create_partition(
                                const unsigned arg_entity_rank ,
                                const PartVector &parts)
{
    enum { KEY_TMP_BUFFER_SIZE = 64 };

    TraceIf("stk::mesh::impl::BucketRepository::get_or_create_partition", LOG_BUCKET);

    ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
            "Entity rank " << arg_entity_rank << " is invalid");

    // Somehow, this can happen.
    ThrowRequireMsg( !m_buckets.empty(),
      "m_buckets is empty! Did you forget to initialize MetaData before creating BulkData?");

    std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

    const size_t part_count = parts.size();
    std::vector<unsigned> key(2 + part_count) ;

    //----------------------------------
    // Key layout:
    // { part_count + 1 , { part_ordinals } , partition_count }
    // Thus partition_count = key[ key[0] ]
    //
    // for upper bound search use the maximum key for a bucket in the partition.
    const unsigned max = static_cast<unsigned>(-1);
    key[0] = part_count+1;
    key[ key[0] ] = max ;

    {
        for ( unsigned i = 0 ; i < part_count ; ++i )
        {
            key[i+1] = parts[i]->mesh_meta_data_ordinal();
        }
    }

    // If the partition is found, the iterator will be right after it, thanks to the
    // trickiness above.
    const std::vector<Partition *>::iterator ik = lower_bound( partitions , &key[0] );
    const bool partition_exists =
            (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , &key[0] );

    if (partition_exists)
    {
        return ik[-1];
    }

    key[key[0]] = 0;
    Partition *partition = new Partition(this, arg_entity_rank, key);
    m_need_sync_from_partitions[arg_entity_rank] = true;
    partitions.insert( ik , partition );

    return partition ;
}



Partition *BucketRepository::get_or_create_partition(
                                const unsigned arg_entity_rank ,
                                const OrdinalVector &parts)
{
    enum { KEY_TMP_BUFFER_SIZE = 64 };

    TraceIf("stk::mesh::impl::BucketRepository::get_or_create_partition", LOG_BUCKET);

    ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
            "Entity rank " << arg_entity_rank << " is invalid");

    // Somehow, this can happen.
    ThrowRequireMsg( !m_buckets.empty(),
      "m_buckets is empty! Did you forget to initialize MetaData before creating BulkData?");

    std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

    const size_t part_count = parts.size();
    std::vector<unsigned> key(2 + part_count) ;

    //----------------------------------
    // Key layout:
    // { part_count + 1 , { part_ordinals } , partition_count }
    // Thus partition_count = key[ key[0] ]
    //
    // for upper bound search use the maximum key for a bucket in the partition.
    const unsigned max = static_cast<unsigned>(-1);
    key[0] = part_count+1;
    key[ key[0] ] = max ;

    {
        for ( unsigned i = 0 ; i < part_count ; ++i ) { key[i+1] = parts[i] ; }
    }

    // If the partition is found, the iterator will be right after it, thanks to the
    // trickiness above.
    const std::vector<Partition *>::iterator ik = lower_bound( partitions , &key[0] );
    const bool partition_exists =
            (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , &key[0] );

    if (partition_exists)
    {
        return ik[-1];
    }

    key[key[0]] = 0;
    Partition *partition = new Partition(this, arg_entity_rank, key);
    m_need_sync_from_partitions[arg_entity_rank] = true;
    partitions.insert( ik , partition );

    return partition ;
}

void BucketRepository::sync_to_partitions()
{
    for ( std::vector<std::vector<Partition *> >::iterator pv_i = m_partitions.begin();
            pv_i != m_partitions.end(); ++pv_i)
    {
        for (std::vector<Partition *>::iterator p_j = pv_i->begin();
                p_j != pv_i->end(); ++p_j)
        {
            (*p_j)->m_buckets.clear();  // Since we assume that m_buckets is up-to-date; happens in testing.
            delete *p_j;
        }
        pv_i->clear();
    }
    m_partitions.clear();
    m_partitions.resize(m_buckets.size());
    m_need_sync_from_partitions.resize(m_buckets.size());

    for ( EntityRank entity_rank = 0; entity_rank < m_buckets.size() ; ++entity_rank )
    {
        std::vector<stk::mesh::impl::Partition *> &partitions = m_partitions[entity_rank];
        std::vector<Bucket*> & buckets = m_buckets[ entity_rank ];

        size_t begin_partition = 0 ; // Offset to first bucket of the partition
        size_t end_partition = 0 ; // Offset to end   bucket of the partition

        //loop over families (of entity_rank) via buckets
        for ( ; begin_partition < buckets.size() ; begin_partition = end_partition )
        {
            Partition *partition =
                    new Partition(this, entity_rank, buckets[begin_partition]->key_vector());
            partitions.push_back(partition);
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

            size_t num_buckets_in_partition = end_partition - begin_partition;
            partition->m_buckets.resize(num_buckets_in_partition);
            std::copy(&buckets[begin_partition],&buckets[end_partition], &partition->m_buckets[0]);

            partition->compute_size();
        }

        // Let each bucket of entity_rank know about its partition.
        size_t num_partitions = partitions.size();
        for (size_t i = 0; i < num_partitions; ++i)
        {
            stk::mesh::impl::Partition *partition = partitions[i];
            std::vector<stk::mesh::Bucket *>::iterator bkt= partition->begin();
            for (;bkt != partition->end(); ++bkt)
            {
                (*bkt)->m_partition = partitions[i];
            }
        }

        m_need_sync_from_partitions[entity_rank] = false;
    }
}

void BucketRepository::sync_from_partitions()
{
    for (EntityRank rank = stk::mesh::MetaData::NODE_RANK; rank < m_partitions.size(); ++rank)
    {
        sync_from_partitions(rank);
    }
}

inline bool isNull(stk::mesh::impl::Partition *p) { return (p ? false : true);}

void BucketRepository::sync_from_partitions(EntityRank rank)
{
    if (!m_need_sync_from_partitions[rank])
    {
        return;
    }

    std::vector<Partition *> &partitions = m_partitions[rank];

    size_t num_partitions = partitions.size();
    size_t num_buckets = 0;
    for (size_t p_i = 0; p_i < num_partitions; ++p_i)
    {
        if (!partitions[p_i]->empty())
        {
            num_buckets += partitions[p_i]->num_buckets();
        }
    }

    m_buckets[rank].resize(num_buckets);

    bool has_hole = false;
    std::vector<Bucket *>::iterator bkts_i = m_buckets[rank].begin();
    for (size_t p_i = 0; p_i < num_partitions; ++p_i)
    {
        Partition &partition = *partitions[p_i];

        if (partition.empty())
        {
            delete partitions[p_i];
            partitions[p_i] = 0;
            has_hole = true;
            continue;
        }
        size_t num_bkts_in_partition = partition.num_buckets();
        std::copy(partition.begin(), partition.end(), bkts_i);
        bkts_i += num_bkts_in_partition;
    }

    if (has_hole)
    {
        std::vector<Partition *>::iterator new_end;
        new_end = std::remove_if(partitions.begin(), partitions.end(), isNull);
        size_t new_size = new_end - partitions.begin();  // OK because has_hole is true.
        partitions.resize(new_size);
    }

    m_need_sync_from_partitions[rank] = false;
}


void BucketRepository::babbleForEntity(EntityRank entity_rank, EntityId entity_id)
{
    Entity entity = m_entity_repo.get_entity( EntityKey(entity_rank, entity_id));

    if (!entity.is_valid())
        return;

    const Bucket *bkt = entity.bucket_ptr();

    if (!bkt)
        return;

    std::cout << m_mesh.synchronized_count()
              << ": P" << m_mesh.parallel_rank() << " Entity(" << entity_rank << "):" << entity_id;
    const std::vector<unsigned> &bucket_key = bkt->key_vector();
    std::cout << " [ ";
    for (size_t i = 0; i < bucket_key.size(); ++i)
    {
        std::cout << bucket_key[i] << " ";
    }
    std::cout << "]  ord = " << entity.bucket_ordinal() << std::endl;
    print(std::cout, " ", *bkt);


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
