/*
 * Partition.cpp
 *
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/diag/Env.hpp>
#include <stk_topology/topology.hpp>

#include <iostream>

// No more than one of these 4 should be set to 1.
#define PARTITION_DONT_COMPRESS_SMALL_PARTITIONS                     0
#define PARTITION_DONT_COMPRESS_SINGLE_BUCKET_PARTITIONS             0
#define PARTITION_DONT_SORT_SMALL_PARTITIONS_ON_COMPRESS             0
#define PARTITION_DONT_SORT_SINGLE_BUCKET_PARTITIONS_ON_COMPRESS     0

// Somewhat contradicts a reason for compress, but useful for debugging
// and experimentation.
#define PARTITION_MAINTAIN_SINGLE_BUCKET_CAPACITY_ON_COMPRESS        0

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

  os << "{Partition m_rank = " << m_rank << ", m_size = " << m_size;

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

std::ostream &Partition::dumpit(std::ostream &os) const
{
  os << "{ Partition (rank = " << m_rank << ")  \n";
  for (std::vector<Bucket *>::const_iterator b_i = begin(); b_i != end(); ++b_i)
  {
    Bucket &b = **b_i;
    print(os, "  ", b );
  }
  os << "}\n";
  return os;
}

std::string Partition::dumpit() const
{
  std::ostringstream output;
  dumpit(output);

  return output.str();
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
  TraceIf("stk::mesh::impl::Partition::add", LOG_PARTITION);
  DiagIf(LOG_PARTITION, "Adding entity: " << print_entity_key(MetaData::get(BulkData::get(*m_repository)), entity.key()));
  TraceIfWatchingDec("stk::mesh::impl::Partition::add", LOG_ENTITY, entity.key(), extra);

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
  entity.m_entityImpl->modified();
  entity.m_entityImpl->set_bucket_and_ordinal(bucket, dst_ordinal);
  bucket->increment_size();
  ++m_size;

  m_updated_since_compress = m_updated_since_sort = true;
  internal_propagate_relocation(entity);
  entity.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());

  DiagIfWatching(LOG_ENTITY, entity.key(),
                 " Bucket: " << *bucket << ", ordinal: " << dst_ordinal);
  DiagIf(LOG_PARTITION, "After add, state is: " << *this);

  return true;
}

void Partition::move_to(Entity entity, Partition &dst_partition)
{
  TraceIf("stk::mesh::impl::Partition::move_to", LOG_PARTITION);
  DiagIf(LOG_PARTITION, "Moving entity: " << print_entity_key(entity));
  TraceIfWatchingDec("stk::mesh::impl::Partition::move_to", LOG_ENTITY, entity.key(), extra);

  Bucket *src_bucket = entity.m_entityImpl->bucket_ptr();
  if (src_bucket && (src_bucket->getPartition() == &dst_partition))
  {
    DiagIfWatching(LOG_ENTITY, entity.key(),
                   " Already on destination partition at bucket " << *src_bucket << ", ordinal "
                   << entity.m_entityImpl->bucket_ordinal());
    return;
  }

  ThrowRequireMsg(src_bucket && (src_bucket->getPartition() == this),
                  "Partition::move_to  cannot move an entity that does not belong to it.");
  unsigned src_ordinal = entity.m_entityImpl->bucket_ordinal();
  DiagIfWatching(LOG_ENTITY, entity.key(),
                 " src_bucket: " << *src_bucket << ", src_ordinal: " << src_ordinal);

  // If the last bucket is full, automatically create a new one.
  Bucket *dst_bucket = dst_partition.get_bucket_for_adds();

  ThrowErrorMsgIf(src_bucket && src_bucket->topology().is_valid() && (src_bucket->topology() != dst_bucket->topology()),
                  "Error: cannot change topology of entity (rank: "
                  << static_cast<stk::topology::rank_t>(entity.entity_rank()) << ", global_id: " << entity.identifier() << ") from "
                  << src_bucket->topology() << "to " << dst_bucket->topology() << "."
                  );

  // Move the entity's data to the new bucket before removing the entity from its old bucket.
  unsigned dst_ordinal = dst_bucket->size();
  dst_bucket->replace_fields(dst_ordinal, *src_bucket, src_ordinal);
  remove(entity, false);

  DiagIfWatching(LOG_ENTITY, entity.key(),
                 " dst_bucket: " << *dst_bucket << ", dst_ordinal: " << dst_ordinal);

  entity.m_entityImpl->modified();
  entity.m_entityImpl->set_bucket_and_ordinal(dst_bucket, dst_ordinal);
  dst_bucket->replace_entity(dst_ordinal, entity) ;
  dst_bucket->increment_size();
  dst_partition.m_updated_since_compress = dst_partition.m_updated_since_sort = true;
  dst_partition.m_size++;

  m_updated_since_compress = m_updated_since_sort = true;
  internal_propagate_relocation(entity);
  entity.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());

  DiagIf(LOG_PARTITION, "After move_to, src state is: " << *this);
  DiagIf(LOG_PARTITION, "After move_to, dst state is: " << dst_partition);
}

bool Partition::remove(Entity e_k, bool not_in_move_to)
{
  TraceIf("stk::mesh::impl::Partition::remove", LOG_PARTITION);
  DiagIf(LOG_PARTITION, "Removing entity: " << print_entity_key(e_k));
  TraceIfWatchingDec("stk::mesh::impl::Partition::remove", LOG_ENTITY, e_k.key(), extra);

  Bucket *bucket_k = e_k.bucket_ptr();
  if (!belongs(bucket_k) || empty())
  {
    return false;
  }
  unsigned ord_k = e_k.m_entityImpl->bucket_ordinal();
  Bucket *last = *(end() - 1);

  DiagIfWatching(LOG_ENTITY, e_k.key(), "  Bucket: " << *bucket_k << ", ordinal: " << ord_k);

  const bool NOT_last_entity_in_last_bucket =
    (last != bucket_k) || (bucket_k->size() != ord_k + 1);
  if ( NOT_last_entity_in_last_bucket )
  {
    // Copy last entity to spot being vacated.
    Entity e_swap = (*last)[ last->size() - 1 ];
    bucket_k->replace_fields(ord_k, *last , last->size() - 1 );
    bucket_k->replace_entity(ord_k, e_swap ) ;
    e_swap.m_entityImpl->set_bucket_and_ordinal(bucket_k, ord_k);

    // Entity field data has relocated.
    internal_propagate_relocation(e_swap);
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
      // The current 'last' bucket in a partition is to be deleted.
      // The previous 'last' bucket becomes the new 'last' bucket in the partition.
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

  e_k.m_entityImpl->set_sync_count(m_repository->m_mesh.synchronized_count());
  m_updated_since_compress = m_updated_since_sort = true;
  --m_size;

  DiagIf(LOG_PARTITION, "After remove, state is: " << *this);

  return true;
}

void Partition::compress(bool force)
{
  TraceIf("stk::mesh::impl::Partition::compress", LOG_PARTITION);

  // An easy optimization.
  if (!force && (empty() || !m_updated_since_compress))
    return;

  const bool dont_compress_small_partitions =                 PARTITION_DONT_COMPRESS_SMALL_PARTITIONS;
  const bool dont_sort_small_partitions_on_compress =         PARTITION_DONT_SORT_SMALL_PARTITIONS_ON_COMPRESS;
  const bool dont_compress_single_bucket_partitions =         PARTITION_DONT_COMPRESS_SINGLE_BUCKET_PARTITIONS;
  const bool dont_sort_single_bucket_partitions_on_compress = PARTITION_DONT_SORT_SINGLE_BUCKET_PARTITIONS_ON_COMPRESS;

  const bool maintain_single_bucket_capacity_on_compress =    PARTITION_MAINTAIN_SINGLE_BUCKET_CAPACITY_ON_COMPRESS;

  const bool single_bucket_case = (m_buckets.size() <= 1);
  const bool small_partition_case = (m_size <= m_repository->bucket_capacity() );

  // Default values
  bool dont_sort_or_compress = false;    // Policy for shortcutting out.
  bool sort_partition = true;

  if (dont_compress_single_bucket_partitions)
  {
    dont_sort_or_compress = single_bucket_case;
  }
  if (dont_compress_small_partitions)
  {
    dont_sort_or_compress = small_partition_case;
  }

  if (dont_sort_single_bucket_partitions_on_compress)
  {
    sort_partition = !single_bucket_case;
  }
  if (dont_sort_small_partitions_on_compress)
  {
    sort_partition = !small_partition_case;
  }

  if (!force && dont_sort_or_compress)
  {
    DiagIf(LOG_PARTITION, "Before compress, state is: " << *this << "\n" << dumpit());
    m_updated_since_compress = false;  // correct in the sense that compress(.) has been called.
    return;
  }

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

  if (sort_partition)
  {
    DiagIf(LOG_PARTITION, "Partition::compress is sorting "
           << m_size << " entities in " << m_buckets.size() << " buckets");
    std::sort( entities.begin(), entities.end(), EntityLess() );
    m_updated_since_sort = false;
  }
  else
  {
    DiagIf(LOG_PARTITION, "Partition::compress is NOT sorting");
  }

  const size_t new_bucket_capacity = (single_bucket_case && maintain_single_bucket_capacity_on_compress
                                      ? m_buckets[0]->capacity() : m_size);

  // Now that the entities hav been sorted, we copy them and their field data into a
  // single bucket in sorted order.
  //
  Bucket * new_bucket = new Bucket( m_repository->m_mesh, m_rank, partition_key, new_bucket_capacity);
  new_bucket->set_first_bucket_in_partition(new_bucket); // partition members point to first bucket
  new_bucket->m_partition = this;

  for(size_t new_ordinal = 0; new_ordinal < entities.size(); ++new_ordinal)
  {
    Entity entity = entities[new_ordinal];
    Bucket& old_bucket = entity.bucket();
    size_t old_ordinal = entity.bucket_ordinal();

    TraceIfWatching("stk::mesh::impl::Partition::compress affects", LOG_ENTITY, entity.key());
    DiagIfWatching(LOG_ENTITY, entity.key(),
                   "  old_bucket: " << old_bucket << ", old_ordinal: " << old_ordinal << '\n'
                   << dumpit()
                   << "\n  new_bucket: " << *new_bucket << ", new_ordinal: " << new_ordinal << " of " << m_size);

    new_bucket->replace_fields(new_ordinal, old_bucket, old_ordinal);
    entity.m_entityImpl->set_bucket_and_ordinal(new_bucket, new_ordinal);
    new_bucket->replace_entity( new_ordinal , entity ) ;
    internal_propagate_relocation(entity);
  }
  new_bucket->m_size = m_size;

  for (std::vector<Bucket *>::iterator b_i = buckets_begin; b_i != buckets_end; ++b_i)
  {
    delete *b_i;
  }

  m_buckets.resize(1);
  m_buckets[0] = new_bucket;
  m_repository->m_need_sync_from_partitions[m_rank] = true;
  m_updated_since_compress = false;
  m_updated_since_sort = false;

  DiagIf(LOG_PARTITION, "After compress, state is: " << *this << "\n" << dumpit());
}

void Partition::sort(bool force)
{
  TraceIf("stk::mesh::impl::Partition::sort", LOG_PARTITION);

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
          TraceIfWatching("stk::mesh::impl::Partition::sort affects", LOG_ENTITY, current.key());
          DiagIfWatching(LOG_ENTITY, current.key(), "  old_bucket: " << j->bucket() << ", old_ordinal: " << i);

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

        TraceIfWatching("stk::mesh::impl::Partition::sort affects", LOG_ENTITY, j->key());
        DiagIfWatching(LOG_ENTITY, j->key(), "  new_bucket: " << *b << ", new_ordinal: " << i);
      }

      // Once a change has occurred then need to propagate the
      // relocation for the remainder of the partition.
      // This allows the propagation to be performed once per
      // entity as opposed to both times the entity is moved.
      if ( change_this_partition )
      {
        internal_propagate_relocation( *j );
      }
    }
  }
  m_updated_since_sort = false;

  if (tmp_bucket)
    delete tmp_bucket;

  DiagIf(LOG_PARTITION, "After sort, state is: " << *this << "\n" << dumpit());
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

void Partition::internal_propagate_relocation( Entity entity )
{
  const EntityRank erank = entity.entity_rank();
  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel )
  {
    const EntityRank rel_rank = rel->entity_rank();
    if ( rel_rank < erank )
    {
      Entity e_to = rel->entity();

      set_field_relations( entity, e_to, rel->identifier() );
    }
    else if ( erank < rel_rank )
    {
      Entity e_from = rel->entity();

      set_field_relations( e_from, entity, rel->identifier() );
    }
  }
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
