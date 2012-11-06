/*
 * Partition.cpp
 *
 *  Created on: Oct 22, 2012
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>

#include <iostream>

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
    os << "{Partition " << std::endl
       << "  m_repository = " << m_repository << "  m_rank = " << m_rank << std::endl
       << "  (" << m_beginBucketIndex << ", " << m_endBucketIndex << ")" << std::endl;

    os << "  legacy partition id : {";
    const std::vector<unsigned> &family_key = get_legacy_partition_id();
    for (size_t i = 0; i < family_key.size(); ++i)
    {
        os << " " << family_key[i];
    }
    os << " }"  << std::endl << "}";

    return os;
}

Partition::~Partition()
{
    // TODO Auto-generated destructor stub
}

BucketRepository &Partition::getRepository(stk::mesh::BulkData &mesh)
{
    return mesh.m_bucket_repository;
}

void Partition::compress()
{
    if (m_endBucketIndex >= m_beginBucketIndex )
        return;

    std::vector<unsigned> partition_key = get_legacy_partition_id();
    //index of bucket in partition
    partition_key[ partition_key[0] ] = 0;

    size_t partition_size = compute_size();
    std::vector<Entity> entities(partition_size);

    std::vector<Bucket *>::iterator buckets_end = end();
    size_t new_i = 0;
    for (std::vector<Bucket *>::iterator b_i = begin(); b_i != buckets_end; ++b_i)
    {
        Bucket &b = **b_i;
        size_t b_size = b.size();
        std::copy(&b.m_entities[0], &b.m_entities[b_size], &entities[new_i]);
        new_i += b_size;
    }

    std::sort( entities.begin(), entities.end(), EntityLess() );

    Bucket * new_bucket = new Bucket( m_repository->m_mesh, m_rank, partition_key, partition_size);
    new_bucket->set_first_bucket_in_partition(new_bucket); // partition members point to first bucket
    new_bucket->set_last_bucket_in_partition(new_bucket);  // First bucket points to new last bucket
    new_bucket->m_partition = this;

    EntityRepository &entity_repository = m_repository->m_entity_repo;
    for(size_t new_ordinal = 0; new_ordinal < entities.size(); ++new_ordinal)
    {
        Entity entity = entities[new_ordinal];
        Bucket& old_bucket = entity.bucket();
        size_t old_ordinal = entity.bucket_ordinal();

        //copy field data from old to new
        new_bucket->replace_fields(new_ordinal, old_bucket, old_ordinal);

        entity_repository.change_entity_bucket( *new_bucket, entity, new_ordinal);
        new_bucket->replace_entity( new_ordinal , entity ) ;
        m_repository->internal_propagate_relocation(entity);
    }
    new_bucket->m_size = partition_size;

    std::vector<Bucket *> &repo_buckets = m_repository->m_buckets[m_rank];
    for (size_t ik = m_beginBucketIndex; ik < m_endBucketIndex; ++ik)
    {
        delete repo_buckets[ik];
        repo_buckets[ik] = NULL;
    }

    repo_buckets[m_beginBucketIndex] = new_bucket;
    m_endBucketIndex = m_beginBucketIndex + 1;
}


void Partition::sort()
{
    if (m_beginBucketIndex >= m_endBucketIndex )
        return;

    std::vector<unsigned> partition_key = get_legacy_partition_id();
    //index of bucket in partition
    partition_key[ partition_key[0] ] = 0;

    size_t partition_size = compute_size();
    std::vector<Entity> entities(partition_size);

    std::vector<Bucket *>::iterator buckets_end = end();

    // Make sure that there is a vacancy somewhere.
    stk::mesh::Bucket *vacancy_bucket = *(buckets_end - 1);
    size_t vacancy_offset = vacancy_bucket->size();
    stk::mesh::Bucket *tmp_bucket = 0;
    if (vacancy_offset >= vacancy_bucket->capacity())
    {
        tmp_bucket = new Bucket(m_repository->m_mesh, m_rank, partition_key, 1);
        vacancy_bucket = tmp_bucket;
        vacancy_offset = 0;
    }

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

    EntityRepository &entity_repository = m_repository->m_entity_repo;
    std::vector<Entity>::iterator j = entities.begin();
    bool change_this_partition = false ;

    for (std::vector<Bucket *>::iterator b_i = begin(); b_i != buckets_end; ++b_i)
    {
        Bucket & b = **b_i;
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j )
        {
            Entity const current = b[i];

            if ( current != *j )
            {
                if ( current.is_valid() )
                {
                    // Move current entity to the vacant spot
                    vacancy_bucket->replace_fields(vacancy_offset, b, i );
                    entity_repository.change_entity_bucket(*vacancy_bucket, current, vacancy_offset);
                    vacancy_bucket->replace_entity( vacancy_offset , current ) ;
                }
                // Set the vacant spot to where the required entity is now.
                vacancy_bucket = & (j->bucket()) ;
                vacancy_offset = j->bucket_ordinal() ;
                vacancy_bucket->replace_entity( vacancy_offset , Entity() ) ;

                // Move required entity to the required spot
                b.replace_fields(i, *vacancy_bucket , vacancy_offset );
                entity_repository.change_entity_bucket( b, *j, i);
                b.replace_entity( i, *j );
                change_this_partition = true ;
            }

            // Once a change has occured then need to propagate the
            // relocation for the remainder of the partition.
            // This allows the propagation to be performed once per
            // entity as opposed to both times the entity is moved.
            if ( change_this_partition )
            {
                m_repository->internal_propagate_relocation( *j );
            }
        }
    }
    if (tmp_bucket)
        delete tmp_bucket;
}
