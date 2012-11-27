/*
 * Partition.cpp
 *
 *  Created on: Oct 22, 2012
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>

#include <iostream>

#define PARTITION_DONT_COMPRESS_SMALL_PARTITIONS
// #define PARTITION_DONT_SORT_SMALL_PARTITIONS_ON_COMPRESS

namespace stk {
namespace mesh {
namespace impl {


std::ostream &operator<<(std::ostream &os, const stk::mesh::impl::Partition &bf)
{
    return bf.streamit(os);
}

} // impl
} // mesh
} // stk

using namespace stk::mesh::impl;

std::ostream &Partition::streamit(std::ostream &os) const
{
    const MetaData & mesh_meta_data = m_repository->m_mesh.mesh_meta_data();

    os << "{Partition " << this << " in m_repository = " << m_repository << " (m_rank) = " << m_rank
       << " on P" << m_repository->m_mesh.parallel_rank();

    os << "  legacy partition id : {";
    const std::vector<unsigned> &family_key = get_legacy_partition_id();
    for (size_t i = 0; i < family_key.size(); ++i)
    {
        os << " " << family_key[i];
        if ((i > 0) && (i < family_key.size() - 1))
        {
            const Part & part = mesh_meta_data.get_part( family_key[i] );
            os << " " << part.name();
        }
    }
    os << " }}";

    return os;
}


Partition::Partition(BucketRepository *repo, EntityRank rank,
                     const std::vector<PartOrdinal> &key)
    : m_repository(repo)
    , m_rank(rank)
    , m_extPartitionKey(key)
    , m_size(0)
    , m_updated_since_compress(false)
    , m_updated_since_sort(false)
{
    // Nada.
}

Partition::~Partition()
{
    size_t num_bkts = m_buckets.size();
    for (size_t i = 0; i < num_bkts; ++i)
    {
        delete m_buckets[i];
        m_buckets[i] = 0;
    }
}


BucketRepository &Partition::getRepository(stk::mesh::BulkData &mesh)
{
    return mesh.m_bucket_repository;
}


bool Partition::add(Entity entity)
{
    if (entity.bucket_ptr())
    {
        // If an entity already belongs to a partition, it cannot be added to one.
        return false;
    }

    // If the last bucket is full, automatically create a new one.
    Bucket *bucket = get_bucket_for_adds();

    unsigned dst_ordinal = bucket->size();
    bucket->initialize_fields(dst_ordinal);
    bucket->replace_entity(dst_ordinal, entity);
    entity.m_entityImpl->log_modified_and_propagate();
    entity.m_entityImpl->set_bucket_and_ordinal(bucket, dst_ordinal);
    bucket->increment_size();
    ++m_size;

    m_updated_since_compress = m_updated_since_sort = true;
    m_repository->internal_propagate_relocation(entity);
    entity.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());

    return true;
}


void Partition::move_to(Entity entity, Partition &dst_partition)
{
    Bucket *src_bucket = entity.m_entityImpl->bucket_ptr();
    if (src_bucket && (src_bucket->getPartition() == &dst_partition))
        return;
    
    ThrowRequireMsg(src_bucket && (src_bucket->getPartition() == this),
                    "Partition::move_to  cannot move an entity that does not belong to it.");

    unsigned src_ordinal = entity.m_entityImpl->bucket_ordinal();
    
    // If the last bucket is full, automatically create a new one.
    Bucket *dst_bucket = dst_partition.get_bucket_for_adds();
    
    // Move the entity's data to the new bucket before removing the entity from its old bucket.
    unsigned dst_ordinal = dst_bucket->size();
    dst_bucket->replace_fields(dst_ordinal, *src_bucket, src_ordinal);
    remove(entity, false);

    entity.m_entityImpl->log_modified_and_propagate();
    entity.m_entityImpl->set_bucket_and_ordinal(dst_bucket, dst_ordinal);
    dst_bucket->replace_entity(dst_ordinal, entity) ;
    dst_bucket->increment_size();
    dst_partition.m_updated_since_compress = dst_partition.m_updated_since_sort = true;
    dst_partition.m_size++;
    
    m_updated_since_compress = m_updated_since_sort = true;
    m_repository->internal_propagate_relocation(entity);
    entity.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());
}


bool Partition::remove(Entity e_k, bool not_in_move_to)
{
    Bucket *bucket_k = e_k.bucket_ptr();
    if (!belongs(bucket_k) || empty())
    {
        return false;
    }
    unsigned ord_k = e_k.m_entityImpl->bucket_ordinal();
    Bucket *last = *(end() - 1);

    if ( (last != bucket_k) || (bucket_k->size() != ord_k + 1))
    {
        // Copy last entity to spot being vacated.
        Entity e_swap = (*last)[ last->size() - 1 ];
        bucket_k->replace_fields(ord_k, *last , last->size() - 1 );
        bucket_k->replace_entity(ord_k, e_swap ) ;
        e_swap.m_entityImpl->set_bucket_and_ordinal(bucket_k, ord_k);

        // Entity field data has relocated.
        m_repository->internal_propagate_relocation(e_swap);
    }

    e_k.m_entityImpl->set_bucket_and_ordinal(0, 0);

    last->decrement_size();
    last->replace_entity( last->size() , Entity() ) ;

    if ( 0 == last->size() )
    {
        size_t num_buckets = m_buckets.size();

        // Don't delete the last bucket now --- might want it later in this modification cycle.
        if (num_buckets > 1)
        {
            Bucket *new_last = m_buckets[num_buckets - 2];
            m_repository->m_need_sync_from_partitions[m_rank] = true;
            m_buckets[0]->set_last_bucket_in_partition(new_last);
            delete m_buckets.back();
            m_buckets.pop_back();
        }
#ifndef PARTITION_HOLD_EMPTY_BUCKET_OPTIMIZATION
        else
        {
            m_repository->m_need_sync_from_partitions[m_rank] = true;
            delete m_buckets.back();
            m_buckets.pop_back();
        }
#endif

    }

    e_k.m_entityImpl->set_bucket_and_ordinal(0, 0);
    e_k.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());
    m_updated_since_compress = m_updated_since_sort = true;
    --m_size;

    return true;
}


void Partition::compress(bool force)
{
    if (!force && (empty() || !m_updated_since_compress))
        return;

#ifdef PARTITION_DONT_COMPRESS_SMALL_PARTITIONS
    if (!force && (m_buckets.size() <= 1))
    {
        return;
    }
#endif

    std::vector<unsigned> partition_key = get_legacy_partition_id();
    //index of bucket in partition
    partition_key[ partition_key[0] ] = 0;

    std::vector<Entity> entities(m_size);

    // Copy the entities (but not their data) into a vector, where they will be sorted.
    //
    std::vector<Bucket *>::iterator buckets_begin, buckets_end;
    buckets_begin = begin();
    buckets_end = end();
    size_t new_i = 0;
    for (std::vector<Bucket *>::iterator b_i = buckets_begin; b_i != buckets_end; ++b_i)
    {
        Bucket &b = **b_i;
        size_t b_size = b.size();
        std::copy(&b.m_entities[0], &b.m_entities[b_size], &entities[new_i]);
        new_i += b_size;
    }

#ifdef PARTITION_DONT_SORT_SMALL_PARTITIONS_ON_COMPRESS
    if (num_buckets() > 1)
    {
        std::sort( entities.begin(), entities.end(), EntityLess() );
    }
#else
    std::sort( entities.begin(), entities.end(), EntityLess() );
#endif

    // Now that the entities hav been sorted, we copy them and their field data into a
    // single bucket in sorted order.
    //
    Bucket * new_bucket = new Bucket( m_repository->m_mesh, m_rank, partition_key, m_size);
    new_bucket->set_first_bucket_in_partition(new_bucket); // partition members point to first bucket
    new_bucket->m_partition = this;

    for(size_t new_ordinal = 0; new_ordinal < entities.size(); ++new_ordinal)
    {
        Entity entity = entities[new_ordinal];
        Bucket& old_bucket = entity.bucket();
        size_t old_ordinal = entity.bucket_ordinal();

        new_bucket->replace_fields(new_ordinal, old_bucket, old_ordinal);
        entity.m_entityImpl->set_bucket_and_ordinal(new_bucket, new_ordinal);
        new_bucket->replace_entity( new_ordinal , entity ) ;
        m_repository->internal_propagate_relocation(entity);
    }
    new_bucket->m_size = m_size;

    for (std::vector<Bucket *>::iterator b_i = buckets_begin; b_i != buckets_end; ++b_i)
    {
        delete *b_i;
    }

    m_buckets.resize(1);
    m_buckets[0] = new_bucket;
    m_repository->m_need_sync_from_partitions[m_rank] = true;
    m_updated_since_compress = m_updated_since_sort = false;
}


void Partition::sort(bool force)
{
    if (!force && (empty() || !m_updated_since_sort))
        return;

    std::vector<unsigned> partition_key = get_legacy_partition_id();
    //index of bucket in partition
    partition_key[ partition_key[0] ] = 0;

    std::vector<Entity> entities(m_size);

    std::vector<Bucket *>::iterator buckets_begin, buckets_end;
    buckets_begin = begin();
    buckets_end = end();

    // Make sure that there is a vacancy somewhere.
    //
    stk::mesh::Bucket *vacancy_bucket = *(buckets_end - 1);
    size_t vacancy_ordinal = vacancy_bucket->size();
    stk::mesh::Bucket *tmp_bucket = 0;

    if (vacancy_ordinal >= vacancy_bucket->capacity())
    {
        // If we need a temporary bucket, it only needs to hold one entity and
        // the corresponding field data.
        tmp_bucket = new Bucket(m_repository->m_mesh, m_rank, partition_key, 1);
        vacancy_bucket = tmp_bucket;
        vacancy_ordinal = 0;
    }

    // Copy all the entities in the Partition into a vector for sorting.
    size_t new_i = 0;
    for (std::vector<Bucket *>::iterator b_i = begin(); b_i != buckets_end; ++b_i)
    {
        Bucket &b = **b_i;
        size_t b_size = b.size();
        std::copy(&b.m_entities[0], &b.m_entities[b_size], &entities[new_i]);
        new_i += b_size;
    }

    std::sort( entities.begin(), entities.end(), EntityLess() );

    // Now that we have the entities sorted, we need to put them and their data
    // in the right order in the buckets.

    std::vector<Entity>::iterator j = entities.begin();
    bool change_this_partition = false ;

    for (std::vector<Bucket *>::iterator b_i = begin(); b_i != buckets_end; ++b_i)
    {
        Bucket *b = *b_i;
        const unsigned n = b->size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j )
        {
            Entity const current = (*b)[i];

            if ( current != *j )
            {
                if ( current.is_valid() )
                {
                    // Move current entity to the vacant spot
                    vacancy_bucket->replace_fields(vacancy_ordinal, *b, i );
                    current.m_entityImpl->set_bucket_and_ordinal(vacancy_bucket, vacancy_ordinal);
                    vacancy_bucket->replace_entity( vacancy_ordinal , current ) ;
                }
                // Set the vacant spot to where the required entity is now.
                vacancy_bucket = & (j->bucket()) ;
                vacancy_ordinal = j->bucket_ordinal() ;
                vacancy_bucket->replace_entity( vacancy_ordinal , Entity() ) ;

                // Move required entity to the required spot
                b->replace_fields(i, *vacancy_bucket , vacancy_ordinal );
                j->m_entityImpl->set_bucket_and_ordinal(b, i);
                b->replace_entity( i, *j );
                change_this_partition = true ;
            }

            // Once a change has occurred then need to propagate the
            // relocation for the remainder of the partition.
            // This allows the propagation to be performed once per
            // entity as opposed to both times the entity is moved.
            if ( change_this_partition )
            {
                m_repository->internal_propagate_relocation( *j );
            }
        }
    }
    m_updated_since_sort = false;

    if (tmp_bucket)
        delete tmp_bucket;
}


void Partition::update_state() const
{
    std::vector<Bucket*>::const_iterator b_e = end();
    for ( std::vector<Bucket*>::const_iterator ik = begin() ; ik != b_e; ++ik )
    {
      (*ik)->update_state();
    }
}


stk::mesh::Bucket *Partition::get_bucket_for_adds()
{
    if (no_buckets())
    {
        m_repository->m_need_sync_from_partitions[m_rank] = true;

        std::vector<unsigned> partition_key = get_legacy_partition_id();
        partition_key[ partition_key[0] ] = 0;
        Bucket *bucket = new Bucket( m_repository->m_mesh, m_rank, partition_key,
                             m_repository->bucket_capacity());
        bucket->m_partition = this;
        bucket->set_first_bucket_in_partition(bucket);
        m_buckets.push_back(bucket);

        return bucket;
    }

    Bucket *bucket = *(end() - 1);  // Last bucket of the partition.

    if (bucket->size() == bucket->capacity())
    {
        m_repository->m_need_sync_from_partitions[m_rank] = true;

        std::vector<unsigned> partition_key = get_legacy_partition_id();
        partition_key[ partition_key[0] ] = m_buckets.size();
        bucket = new Bucket( m_repository->m_mesh, m_rank, partition_key,
                             m_repository->bucket_capacity());
        bucket->m_partition = this;
        bucket->set_first_bucket_in_partition(m_buckets[0]);
        m_buckets[0]->set_last_bucket_in_partition(bucket);
        m_buckets.push_back(bucket);
    }

    return bucket;
}


void Partition::reverseEntityOrderWithinBuckets()
{
    for (std::vector<stk::mesh::Bucket *>::iterator b_i = begin(); b_i != end(); ++b_i)
    {
        stk::mesh::Bucket & b = **b_i;
        std::reverse(b.m_entities.begin(), b.m_entities.end());
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i)
        {
            b.m_entities[i].m_entityImpl->set_bucket_and_ordinal(*b_i, i);
        }
    }
}
