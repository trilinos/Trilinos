/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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

BulkData::BulkData( const MetaData & mesh_meta_data ,
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
  parallel_machine_barrier( m_parallel_machine );

  if ( m_sync_state == MODIFIABLE ) return false ;

  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );

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
                                   const std::vector<Part*> & parts )
{
  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );

  std::pair< Entity * , bool > result = m_entity_repo.internal_create_entity( key );

  if ( result.second ) {
    // A new application-created entity
    m_entity_repo.set_entity_owner_rank( *(result.first), m_parallel_rank);
    m_entity_repo.set_entity_sync_count( *(result.first), m_sync_count);
  }
  else {
    // An existing entity, the owner must match.
    require_entity_owner( * result.first , m_parallel_rank );
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  std::vector<Part*> rem ;
  std::vector<Part*> add( parts );
  add.push_back( owns );

  change_entity_parts( * result.first , add , rem );

  // m_transaction_log.insert_entity ( *(result.first) );

  return * result.first ;
}

//----------------------------------------------------------------------

namespace {

// TODO Change the methods below to requirements (private, const invariant checkers)

// Returns false if there is a problem. It is expected that
// verify_change_parts will be called if quick_verify_change_part detects
// a problem, therefore we leave the generation of an exception to
// verify_change_parts. We want this function to be as fast as
// possible.
inline bool quick_verify_change_part(const Entity& entity,
                                     const Part* part,
                                     const unsigned entity_rank,
                                     const unsigned undef_rank)
{
  const unsigned part_rank = part->primary_entity_rank();

  // The code below is coupled with the code in verify_change_parts. If we
  // change what it means for a part to be valid, code will need to be
  // changed in both places unfortunately.
  const bool intersection_ok = part->intersection_of().empty();
  const bool rel_target_ok   = ( part->relations().empty() ||
                                 part != part->relations().begin()->m_target );
  const bool rank_ok         = ( entity_rank == part_rank ||
                                 undef_rank  == part_rank );

  return intersection_ok && rel_target_ok && rank_ok;
}

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 2) PartRelation target part
// 3) Part that does not match the entity rank.

void verify_change_parts( const MetaData   & meta ,
                          const Entity     & entity ,
                          const PartVector & parts )
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

}

void BulkData::change_entity_parts(
  Entity & entity ,
  const PartVector & add_parts ,
  const PartVector & remove_parts )
{
  require_ok_to_modify();

  require_entity_owner( entity , m_parallel_rank );

  const EntityRank entity_rank = entity.entity_rank();
  const EntityRank undef_rank  = InvalidEntityRank;

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  PartVector a_parts( add_parts );
  bool quick_verify_check = true;

  for ( PartVector::const_iterator
        ia = add_parts.begin(); ia != add_parts.end() ; ++ia ) {
    quick_verify_check = quick_verify_check &&
      quick_verify_change_part(entity, *ia, entity_rank, undef_rank);
    a_parts.insert( a_parts.end(), (*ia)->supersets().begin(),
                                   (*ia)->supersets().end() );
  }

  order( a_parts );

  PartVector r_parts ;

  for ( PartVector::const_iterator
        ir = remove_parts.begin(); ir != remove_parts.end() ; ++ir ) {

    // The following guards should be in the public interface to
    // changing parts.  However, internal mechanisms such as changing
    // ownership calls this function to add or remove an entity from
    // the three special parts.  Without refactoring, these guards
    // cannot be put in place.
    /*
    ThrowErrorMsgIf( m_mesh_meta_data.universal_part() == **ir,
                     "Cannot remove entity from universal part" );
    ThrowErrorMsgIf( m_mesh_meta_data.locally_owned_part() == **ir,
                     "Cannot remove entity from locally owned part" );
    ThrowErrorMsgIf( m_mesh_meta_data.globally_shared_part() == **ir,
                     "Cannot remove entity from globally shared part" );
    */

    quick_verify_check = quick_verify_check &&
      quick_verify_change_part(entity, *ir, entity_rank, undef_rank);

    if ( ! contain( a_parts , **ir ) ) {
      r_parts.push_back( *ir );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( entity.bucket().member ( **cur_part ) )
          r_parts.push_back ( *cur_part );
    }
  }

  order( r_parts );

  // If it looks like we have a problem, run the full check and we should
  // expect to see an exception thrown; otherwise, only do the full check in
  // debug mode because it incurs significant overhead.
  if ( ! quick_verify_check ) {
    verify_change_parts( m_mesh_meta_data , entity , a_parts );
    verify_change_parts( m_mesh_meta_data , entity , r_parts );
    ThrowRequireMsg(false, "Expected throw from verify methods above.");
  }
  else {
#ifndef NDEBUG
    verify_change_parts( m_mesh_meta_data , entity , a_parts );
    verify_change_parts( m_mesh_meta_data , entity , r_parts );
#endif
  }

  internal_change_entity_parts( entity , a_parts , r_parts );

  return ;
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

    parts_total.assign( bucket_parts.first , bucket_parts.second );

    filter_out( parts_total , remove_parts , parts_removed );
  }

  merge_in( parts_total , add_parts );

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

  while ( ! entity.relations().empty() ) {
    destroy_relation( entity , * entity.relations().back().entity() );
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
      m_entity_repo.set_entity_owner_rank( *(result.first), m_parallel_rank);
      m_entity_repo.set_entity_sync_count( *(result.first), m_sync_count);

      //add entity to 'owned' part
      change_entity_parts( * result.first , add , rem );
      requested_entities.push_back(result.first);
    }
  }
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

