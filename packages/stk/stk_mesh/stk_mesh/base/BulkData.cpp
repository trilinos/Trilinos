/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
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

}

//----------------------------------------------------------------------

BulkData::BulkData( MetaData & mesh_meta_data ,
                    ParallelMachine parallel ,
                    unsigned bucket_max_size )
  : m_entities_index( parallel, convert_entity_keys_to_spans(mesh_meta_data) ),
    m_entity_repo(),
    m_bucket_repository(
        *this, bucket_max_size,
        mesh_meta_data.entity_rank_count(),
        m_entity_repo
        ),
    m_entity_comm(),
    m_ghosting(),

    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_sync_count( 0 ),
    m_sync_state( MODIFIABLE ),
    m_meta_data_verified( false )
{
  create_ghosting( std::string("shared") );
  create_ghosting( std::string("shared_aura") );

  m_sync_state = SYNCHRONIZED ;
}

BulkData::~BulkData()
{
  try {
    while ( ! m_ghosting.empty() ) {
      delete m_ghosting.back();
      m_ghosting.pop_back();
    }
  } catch(...){}

  try { m_entity_comm.clear(); } catch(...){}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::require_ok_to_modify() const
{
  ThrowRequireMsg( m_sync_state != SYNCHRONIZED,
                   "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity & entity ,
                                     unsigned owner ) const
{
  const bool error_not_owner = owner != entity.owner_rank() ;

  ThrowRequireMsg( !error_not_owner,
      "Entity " << print_entity_key(entity) << " owner is " <<
                   entity.owner_rank() << ", expected " << owner);
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = entity_id_valid(ent_id);
  const bool ok_rank = ent_rank < rank_count ;

  ThrowRequireMsg( ok_rank,
                   "Bad key rank: " << ent_rank << " for id " << ent_id );

  ThrowRequireMsg( ok_id, "Bad key id for key: " <<
      print_entity_key(m_mesh_meta_data, EntityKey(ent_rank, ent_id) ) );
}

void BulkData::require_metadata_committed() const
{
  ThrowRequireMsg( m_mesh_meta_data.is_commit(), "MetaData not committed." );
}

//----------------------------------------------------------------------

bool BulkData::modification_begin()
{
  Trace_("stk::mesh::BulkData::modification_begin");

  parallel_machine_barrier( m_parallel_machine );

  if ( m_sync_state == MODIFIABLE ) return false ;

  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    if (parallel_size() > 1) {
      verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );
    }

    m_meta_data_verified = true ;

    m_bucket_repository.declare_nil_bucket();
  }
  else {
    ++m_sync_count ;

    // Clear out the previous transaction information
    // m_transaction_log.flush();

    m_entity_repo.clean_changes();
  }

  m_sync_state = MODIFIABLE ;

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.

//----------------------------------------------------------------------

Entity & BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                   const PartVector & parts )
{
  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );
  TraceIfWatching("stk::mesh::BulkData::declare_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "declaring entity with parts " << parts);

  std::pair< Entity * , bool > result = m_entity_repo.internal_create_entity( key );

  Entity* declared_entity = result.first;

  if ( result.second ) {
    // A new application-created entity
    m_entity_repo.set_entity_owner_rank( *declared_entity, m_parallel_rank);
    m_entity_repo.set_entity_sync_count( *declared_entity, m_sync_count);
    DiagIfWatching(LOG_ENTITY, key, "new entity: " << *declared_entity);
  }
  else {
    // An existing entity, the owner must match.
    require_entity_owner( *declared_entity , m_parallel_rank );
    DiagIfWatching(LOG_ENTITY, key, "existing entity: " << *declared_entity);
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  std::vector<Part*> rem ;
  std::vector<Part*> add( parts );
  add.push_back( owns );

  change_entity_parts( *declared_entity , add , rem );

  // m_transaction_log.insert_entity ( *(result.first) );

  return *declared_entity ;
}

//----------------------------------------------------------------------

// TODO Change the methods below to requirements (private, const invariant checkers)

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 2) PartRelation target part
// 3) Part that does not match the entity rank.

void BulkData::internal_verify_change_parts( const MetaData   & meta ,
                                             const Entity     & entity ,
                                             const PartVector & parts ) const
{
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  const EntityRank undef_rank  = InvalidEntityRank;
  const EntityRank entity_rank = entity.entity_rank();

  bool ok = true ;
  std::ostringstream msg ;

  for ( PartVector::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part * const p = *i ;
    const unsigned part_rank = p->primary_entity_rank();

    // The code below is coupled with the code in quick_verify_change_part.
    // If we change what it means for a part to be valid, code will need to be
    // changed in both places unfortunately.
    const bool error_intersection = ! p->intersection_of().empty();
    const bool error_rel_target   = ! p->relations().empty() &&
                                    p == p->relations().begin()->m_target ;
    const bool error_rank = entity_rank != part_rank &&
                            undef_rank  != part_rank ;

    if ( error_intersection || error_rel_target || error_rank ) {
      if ( ok ) {
        ok = false ;
        msg << "change parts for entity " << print_entity_key( entity );
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }

      msg << p->name() << "[" ;
      if ( part_rank < rank_names.size() ) {
        msg << rank_names[ part_rank ];
      }
      else {
        msg << part_rank ;
      }
      msg << "] " ;
      if ( error_intersection ) { msg << "is_intersection " ; }
      if ( error_rel_target )   { msg << "is_relation_target " ; }
      if ( error_rank )         { msg << "is_bad_rank " ; }
    }
  }

  ThrowErrorMsgIf( !ok, msg.str() << "}" );
}

//----------------------------------------------------------------------

namespace {

void filter_out( std::vector<unsigned> & vec ,
                 const PartVector & parts ,
                 PartVector & removed )
{
  std::vector<unsigned>::iterator i , j ;
  i = j = vec.begin();

  PartVector::const_iterator ip = parts.begin() ;

  while ( j != vec.end() && ip != parts.end() ) {
    Part * const p = *ip ;
    if      ( p->mesh_meta_data_ordinal() < *j ) { ++ip ; }
    else if ( *j < p->mesh_meta_data_ordinal() ) { *i = *j ; ++i ; ++j ; }
    else {
      removed.push_back( p );
      ++j ;
      ++ip ;
    }
  }

  if ( i != j ) { vec.erase( i , j ); }
}

void merge_in( std::vector<unsigned> & vec , const PartVector & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  PartVector::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = (*ip)->mesh_meta_data_ordinal();

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      // Now have: ord == *i
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    vec.push_back( ord );
  }
}

}

//  The 'add_parts' and 'remove_parts' are complete and disjoint.
//  Changes need to have parallel resolution during
//  modification_end.

void BulkData::internal_change_entity_parts(
  Entity & entity ,
  const PartVector & add_parts ,
  const PartVector & remove_parts )
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity.key());
  DiagIfWatching(LOG_ENTITY, entity.key(), "entity state: " << entity);
  DiagIfWatching(LOG_ENTITY, entity.key(), "add_parts: " << add_parts);
  DiagIfWatching(LOG_ENTITY, entity.key(), "remove_parts: " << remove_parts);

  Bucket * const k_old = m_entity_repo.get_entity_bucket( entity );

  const unsigned i_old = entity.bucket_ordinal() ;

  if ( k_old && k_old->member_all( add_parts ) &&
              ! k_old->member_any( remove_parts ) ) {
    // Is already a member of all add_parts,
    // is not a member of any remove_parts,
    // thus nothing to do.
    return ;
  }

  PartVector parts_removed ;

  std::vector<unsigned> parts_total ; // The final part list

  //--------------------------------

  if ( k_old ) {
    // Keep any of the existing bucket's parts
    // that are not a remove part.
    // This will include the 'intersection' parts.
    //
    // These parts are properly ordered and unique.

    const std::pair<const unsigned *, const unsigned*>
      bucket_parts = k_old->superset_part_ordinals();

    const unsigned * parts_begin = bucket_parts.first;
    const unsigned * parts_end   = bucket_parts.second;

    const unsigned num_bucket_parts = parts_end - parts_begin;
    parts_total.reserve( num_bucket_parts + add_parts.size() );
    parts_total.insert( parts_total.begin(), parts_begin , parts_end);

    if ( !remove_parts.empty() ) {
      parts_removed.reserve(remove_parts.size());
      filter_out( parts_total , remove_parts , parts_removed );
    }
  }
  else {
    parts_total.reserve(add_parts.size());
  }

  if ( !add_parts.empty() ) {
    merge_in( parts_total , add_parts );
  }

  if ( parts_total.empty() ) {
    // Always a member of the universal part.
    const unsigned univ_ord =
      m_mesh_meta_data.universal_part().mesh_meta_data_ordinal();
    parts_total.push_back( univ_ord );
  }

  //--------------------------------
  // Move the entity to the new bucket.

  Bucket * k_new =
    m_bucket_repository.declare_bucket(
        entity.entity_rank(),
        parts_total.size(),
        & parts_total[0] ,
        m_mesh_meta_data.get_fields()
        );

  // If changing buckets then copy its field values from old to new bucket

  if ( k_old ) {
    m_bucket_repository.copy_fields( *k_new , k_new->size() , *k_old , i_old );
  }
  else {
    m_bucket_repository.zero_fields( *k_new , k_new->size() );
  }

  // Set the new bucket
  m_entity_repo.change_entity_bucket( *k_new, entity, k_new->size() );
  m_bucket_repository.add_entity_to_bucket( entity, *k_new );

  // If changing buckets then remove the entity from the bucket,
  if ( k_old ) { m_bucket_repository.remove_entity( k_old , i_old ); }

  // Update the change counter to the current cycle.
  m_entity_repo.set_entity_sync_count( entity, m_sync_count );

  // Propagate part changes through the entity's relations.

  internal_propagate_part_changes( entity , parts_removed );
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity * & entity_in )
{
  Entity & entity = *entity_in ;

  TraceIfWatching("stk::mesh::BulkData::destroy_entity", LOG_ENTITY, entity.key());
  DiagIfWatching(LOG_ENTITY, entity.key(), "entity state: " << entity);

  require_ok_to_modify( );

  bool has_upward_relation = false ;

  for ( PairIterRelation
        irel = entity.relations() ;
        ! irel.empty() && ! has_upward_relation ; ++irel ) {

    has_upward_relation = entity.entity_rank() <= irel->entity_rank();
  }

  if ( has_upward_relation ) { return false ; }

  if (  EntityLogDeleted == entity.log_query() ) {
    // Cannot already be destroyed.
    return false ;
  }
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
                      * entity.relations().back().entity(),
                      entity.relations().back().identifier());
  }

  // We need to save these items and call remove_entity AFTER the call to
  // destroy_later because remove_entity may destroy the bucket
  // which would cause problems in m_entity_repo.destroy_later because it
  // makes references to the entity's original bucket.
  Bucket& orig_bucket = entity.bucket();
  unsigned orig_bucket_ordinal = entity.bucket_ordinal();

  // Set the bucket to 'bucket_nil' which:
  //   1) has no parts at all
  //   2) has no field data
  //   3) has zero capacity
  //
  // This keeps the entity-bucket methods from catastrophically failing
  // with a bad bucket pointer.

  m_entity_repo.destroy_later( entity, m_bucket_repository.get_nil_bucket() );

  m_bucket_repository.remove_entity( &orig_bucket , orig_bucket_ordinal );

  // Add destroyed entity to the transaction
  // m_transaction_log.delete_entity ( *entity_in );

  // Set the calling entity-pointer to NULL;
  // hopefully the user-code will clean up any outstanding
  // references to this entity.

  entity_in = NULL ;

  return true ;
}

//----------------------------------------------------------------------

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity *>& requested_entities)
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
      std::pair<Entity *, bool> result = m_entity_repo.internal_create_entity(key);

      //if an entity is declare with the declare_entity function in
      //the same modification cycle as the generate_new_entities
      //function, and it happens to generate a key that was declare
      //previously in the same cycle it is an error
      ThrowErrorMsgIf( ! result.second,
                       "Generated " << print_entity_key(m_mesh_meta_data, key) <<
                       " which was already used in this modification cycle.");

      // A new application-created entity

      Entity* new_entity = result.first;

      m_entity_repo.set_entity_owner_rank( *new_entity, m_parallel_rank);
      m_entity_repo.set_entity_sync_count( *new_entity, m_sync_count);

      //add entity to 'owned' part
      change_entity_parts( *new_entity , add , rem );
      requested_entities.push_back(new_entity);
    }
  }
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

