/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.               */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/util/StaticAssert.hpp>

#include <stk_util/diag/Trace.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <boost/foreach.hpp>

namespace stk {
namespace mesh {

namespace {

std::vector< parallel::DistributedIndex::KeySpan>
convert_entity_keys_to_spans( const MetaData & meta )
{
  // Make sure the distributed index can handle the EntityKey

  enum { OK = StaticAssert<
                SameType< EntityKey::raw_key_type,
                          parallel::DistributedIndex::KeyType >::value >::OK };

  // Default constructed EntityKey has all bits set.

  const EntityKey invalid_key ;
  const EntityId  min_id = 1 ;
  const EntityId  max_id = invalid_key.id();

  const size_t rank_count = meta.entity_rank_count();

  std::vector< parallel::DistributedIndex::KeySpan> spans( rank_count );

  for ( size_t rank = 0 ; rank < rank_count ; ++rank ) {
    EntityKey key_min( rank , min_id );
    EntityKey key_max( rank , max_id );
    spans[rank].first  = key_min.raw_key();
    spans[rank].second = key_max.raw_key();
  }

  return spans ;
}

void ensure_part_superset_consistency( const Entity entity )
{
  std::ostringstream errs;
  PartVector parts;
  Bucket& bucket = entity.bucket();
  bucket.supersets(parts);
  BOOST_FOREACH(Part* part, parts) {
    const PartVector& supersets = part->supersets();
    BOOST_FOREACH(Part* superset, supersets) {
      if (!bucket.member(*superset)) {
        errs << "  Due to being a member part " << part->name() << ", should have been a member of " << superset->name() << std::endl;
      }
    }
  }
  ThrowRequireMsg( errs.str() == "",
                   "Entity " << print_entity_key(entity) << " has bad part list:\n" << errs.str() );
}

}

//----------------------------------------------------------------------

BulkData::BulkData( MetaData & mesh_meta_data ,
                    ParallelMachine parallel ,
                    unsigned bucket_max_size ,
                    bool use_memory_pool )
  : m_entities_index( parallel, convert_entity_keys_to_spans(mesh_meta_data) ),
    m_entity_repo(use_memory_pool),
    m_bucket_repository(
        *this, bucket_max_size,
        mesh_meta_data.entity_rank_count(),
        m_entity_repo
        ),
    m_entity_comm_list(),
    m_ghosting(),

    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_sync_count( 0 ),
    m_sync_state( MODIFIABLE ),
    m_meta_data_verified( false ),
    m_optimize_buckets(false),
    m_mesh_finalized(false)
{
  create_ghosting( "shared" );
  create_ghosting( "shared_aura" );

  m_sync_state = SYNCHRONIZED ;
}

BulkData::~BulkData()
{
  while ( ! m_ghosting.empty() ) {
    delete m_ghosting.back();
    m_ghosting.pop_back();
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::require_ok_to_modify() const
{
  ThrowRequireMsg( m_sync_state != SYNCHRONIZED,
                   "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity entity ,
                                     unsigned owner ) const
{
  if (parallel_size() == 1) {
    //no error-check if running in serial
    return;
  }

  const bool error_not_owner = owner != entity.owner_rank() ;

  ThrowRequireMsg( !error_not_owner,
      "Entity " << print_entity_key(entity) << " owner is " <<
                   entity.owner_rank() << ", expected " << owner);
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = entity_id_valid(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == MetaData::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  ThrowRequireMsg( ok_rank,
                   "Bad key rank: " << ent_rank << " for id " << ent_id );

  ThrowRequireMsg( ok_id, "Bad key id for key: " <<
      print_entity_key(m_mesh_meta_data, EntityKey(ent_rank, ent_id) ) );
}

void BulkData::require_metadata_committed()
{
  if (!m_mesh_meta_data.is_commit()) {
    m_mesh_meta_data.commit();
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_begin()
{
  Trace_("stk::mesh::BulkData::modification_begin");

  parallel_machine_barrier( m_parallel_machine );

  ThrowRequireMsg( m_mesh_finalized == false, "Unable to modifiy, BulkData has been finalized.");

  if ( m_sync_state == MODIFIABLE && m_mesh_finalized == false ) return false ;

  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    if (parallel_size() > 1) {
      verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );
    }

    m_meta_data_verified = true ;
  }
  else {
    ++m_sync_count ;

    // Clear out the previous transaction information
    // m_transaction_log.flush();

    m_entity_repo.clean_changes();
  }

  // // It might be overkill to call this on every modification cycle.
  // m_bucket_repository.sync_to_partitions();

  m_sync_state = MODIFIABLE ;

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.

//----------------------------------------------------------------------

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );
  TraceIfWatching("stk::mesh::BulkData::declare_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "declaring entity with parts " << parts);

  std::pair< Entity , bool > result = m_entity_repo.internal_create_entity( key );

  Entity declared_entity = result.first;

  if ( result.second ) {
    // A new application-created entity
    m_entity_repo.set_entity_owner_rank(declared_entity, m_parallel_rank);
    m_entity_repo.set_entity_sync_count( declared_entity, m_sync_count);
    DiagIfWatching(LOG_ENTITY, key, "new entity: " << declared_entity);
  }
  else {
    // An existing entity, the owner must match.
    require_entity_owner( declared_entity , m_parallel_rank );
    DiagIfWatching(LOG_ENTITY, key, "existing entity: " << declared_entity);
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  PartVector rem ;
  PartVector add( parts );
  add.push_back( owns );

  change_entity_parts( declared_entity , add , rem );

  // m_transaction_log.insert_entity ( *(result.first) );

  return declared_entity ;
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    Part& universal = mesh_meta_data().universal_part();
    return declare_entity(ent_rank, ent_id, universal);
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
  ThrowAssertMsg(parallel_size() == 1,
                 "change_entity_id only supported in serial");

  require_ok_to_modify();
  require_good_rank_and_id(entity.entity_rank(),id);

  EntityKey new_key(entity.entity_rank(),id);
  EntityKey old_key = entity.key();

  m_entity_repo.update_entity_key(new_key, old_key);
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity entity )
{
  TraceIfWatching("stk::mesh::BulkData::destroy_entity", LOG_ENTITY, entity.key());
  DiagIfWatching(LOG_ENTITY, entity.key(), "entity state: " << entity);

  require_ok_to_modify();

  if (!entity.is_valid()) {
    return false;
  }

  bool has_upward_relation = false ;

  for ( PairIterRelation
        irel = entity.relations() ;
        ! irel.empty() && ! has_upward_relation ; ++irel ) {

    has_upward_relation = entity.entity_rank() <= irel->entity_rank();
  }

  if ( has_upward_relation ) { return false ; }

  //------------------------------
  // Immediately remove it from relations and buckets.
  // Postpone deletion until modification_end to be sure that
  // 1) No attempt is made to re-create it.
  // 2) Parallel index is cleaned up.
  // 3) Parallel sharing is cleaned up.
  // 4) Parallel ghosting is cleaned up.
  //
  // Must clean up the parallel lists before fully deleting the entity.

  // It is important that relations be destroyed in reverse order so that
  // the higher (back) relations are destroyed first.
  while ( ! entity.relations().empty() ) {
    destroy_relation( entity ,
                      entity.relations().back().entity(),
                      entity.relations().back().relation_ordinal());
  }

  // We need to save these items and call remove_entity AFTER the call to
  // destroy_later because remove_entity may destroy the bucket
  // which would cause problems in m_entity_repo.destroy_later because it
  // makes references to the entity's original bucket.

  // Need to invalidate Entity handles in comm-list
  std::vector<EntityCommListInfo>::iterator lb_itr =
    std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), entity.key());
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == entity.key()) {
    lb_itr->entity = Entity();
  }

  m_entities_index.register_removed_key( entity.key().raw_key() );

  entity.bucket().getPartition()->remove(entity);
  m_entity_repo.destroy_entity( entity );

  return true ;
}

//----------------------------------------------------------------------

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity>& requested_entities)
{
  Trace_("stk::mesh::BulkData::generate_new_entities");

  typedef stk::parallel::DistributedIndex::KeyType KeyType;
  std::vector< std::vector<KeyType> >
    requested_key_types;
  m_entities_index.generate_new_keys(requests, requested_key_types);

  //generating 'owned' entities
  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  std::vector<Part*> rem ;
  std::vector<Part*> add;
  add.push_back( owns );

  requested_entities.clear();
  unsigned cnt=0;
  for (std::vector< std::vector<KeyType> >::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr) {
    const std::vector<KeyType>& key_types = *itr;
    for (std::vector<KeyType>::const_iterator
        kitr = key_types.begin(); kitr != key_types.end(); ++kitr) {
      ++cnt;
    }
  }
  requested_entities.reserve(cnt);

  for (std::vector< std::vector<KeyType> >::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr) {
    const std::vector<KeyType>& key_types = *itr;
    for (std::vector<KeyType>::const_iterator
        kitr = key_types.begin(); kitr != key_types.end(); ++kitr) {
      EntityKey key(&(*kitr));
      require_good_rank_and_id(key.rank(), key.id());
      std::pair<Entity , bool> result = m_entity_repo.internal_create_entity(key);

      //if an entity is declare with the declare_entity function in
      //the same modification cycle as the generate_new_entities
      //function, and it happens to generate a key that was declare
      //previously in the same cycle it is an error
      ThrowErrorMsgIf( ! result.second,
                       "Generated " << print_entity_key(m_mesh_meta_data, key) <<
                       " which was already used in this modification cycle.");

      // A new application-created entity

      Entity new_entity = result.first;

      m_entity_repo.set_entity_owner_rank( new_entity, m_parallel_rank);
      m_entity_repo.set_entity_sync_count( new_entity, m_sync_count);

      //add entity to 'owned' part
      change_entity_parts( new_entity , add , rem );
      requested_entities.push_back(new_entity);
    }
  }
}

bool BulkData::in_shared(EntityKey key, unsigned proc) const
{
  PairIterEntityComm sharing = entity_comm_sharing(key);
  for ( ; !sharing.empty(); ++sharing ) {
    if ( proc == sharing->proc ) {
      return true ;
    }
  }
  return false ;
}

bool BulkData::in_send_ghost( EntityKey key , unsigned proc ) const
{
  const unsigned owner_rank = entity_comm_owner(key);
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id != 0 &&
         ec->proc     != owner_rank &&
         ec->proc     == proc ) {
      return true;
    }
  }
  return false;
}

bool BulkData::in_ghost( const Ghosting & ghost , EntityKey key , unsigned p ) const
{
  // Ghost communication from owner.
  EntityCommInfo tmp( ghost.ordinal() , p );

  PairIterEntityComm ec = entity_comm(key);
  std::vector<EntityCommInfo>::const_iterator i =
    std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

void BulkData::comm_procs( EntityKey key, std::vector<unsigned> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<unsigned>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}

void BulkData::comm_procs( const Ghosting & ghost ,
                           EntityKey key, std::vector<unsigned> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id == ghost.ordinal() ) {
      procs.push_back( ec->proc );
    }
  }
}

void BulkData::internal_change_owner_in_comm_data(const EntityKey& key, unsigned new_owner)
{
  const bool changed = m_entity_comm_map.change_owner_rank(key, new_owner);
  if (changed) {
    std::vector<EntityCommListInfo>::iterator lb_itr = std::lower_bound(m_entity_comm_list.begin(),
                                                                        m_entity_comm_list.end(),
                                                                        key);
    if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
      lb_itr->owner = new_owner;
    }
  }
}

void BulkData::internal_sync_comm_list_owners()
{
  for (size_t i = 0, e = m_entity_comm_list.size(); i < e; ++i) {
    m_entity_comm_list[i].owner = m_entity_comm_list[i].entity.owner_rank();
  }
}

} // namespace mesh
} // namespace stk
