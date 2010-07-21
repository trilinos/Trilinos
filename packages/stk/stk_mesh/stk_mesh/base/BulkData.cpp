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
#include <assert.h>

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

void BulkData::assert_ok_to_modify( const char * method ) const
{
  if ( m_sync_state == SYNCHRONIZED ) {
    std::string msg( method );
    msg.append( ": FAILED, NOT in the ok-to-modify state" );
    throw std::runtime_error( msg );
  }
}

void BulkData::assert_entity_owner( const char * method ,
                                    const Entity & e ,
                                    unsigned owner ) const
{
  const bool error_not_owner = owner != e.owner_rank() ;

  if ( error_not_owner ) {
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , m_mesh_meta_data , e.key() );
    msg << " ) FAILED" ;

    msg << " : Owner( " << e.owner_rank()
        << " ) != Required( " << owner << " )" ;

    throw std::runtime_error( msg.str() );
  }
}

void BulkData::assert_good_key( const char * method ,
                                const EntityKey & key ) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = 0 < entity_id( key );
  const bool ok_type = entity_rank( key ) < rank_count ;

  if ( ! ok_type || ! ok_id ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    if ( ! ok_type ) {
      msg << entity_rank( key ) << "-"
          << entity_id( key ) << " : BAD KEY TYPE" ;
    }
    else {
      print_entity_key( msg , m_mesh_meta_data , key );
      msg << " : BAD KEY ID" ;
    }
    msg << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_begin()
{
  static const char method[] = "stk::mesh::BulkData::modification_begin" ;

  parallel_machine_barrier( m_parallel_machine );

  if ( m_sync_state == MODIFIABLE ) return false ;

  if ( ! m_meta_data_verified ) {
    m_mesh_meta_data.assert_committed( method );

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


void BulkData::verify_type_and_id(const char* calling_method,
                                  EntityRank ent_type, EntityId ent_id) const
{
  m_mesh_meta_data.assert_entity_rank( calling_method , ent_type );

  if (!entity_id_valid(ent_id)) {
    std::ostringstream msg;
    msg << calling_method << ": ent_id not valid";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.



//----------------------------------------------------------------------

Entity & BulkData::declare_entity( EntityRank ent_type , EntityId ent_id ,
				   const std::vector<Part*> & parts )
{
  static const char method[] = "stk::mesh::BulkData::declare_entity" ;

  assert_ok_to_modify( method );

  verify_type_and_id("BulkData::declare_entity", ent_type, ent_id);

  EntityKey key( ent_type , ent_id );

  assert_good_key( method , key );

  std::pair< Entity * , bool > result = m_entity_repo.internal_create_entity( key );

  if ( result.second ) {
    // A new application-created entity
    m_entity_repo.set_entity_owner_rank( *(result.first), m_parallel_rank);
    m_entity_repo.set_entity_sync_count( *(result.first), m_sync_count);
  }
  else {
    // An existing entity, the owner must match.
    assert_entity_owner( method , * result.first , m_parallel_rank );
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

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 2) PartRelation target part
// 3) Part that does not match the entity type.

void verify_change_parts( const char * method ,
                          const MetaData   & meta ,
                          const Entity     & entity ,
                          const PartVector & parts )
{
  const std::vector<std::string> & type_names = meta.entity_rank_names();
  const unsigned undef_rank  = std::numeric_limits<unsigned>::max();
  const unsigned entity_rank = entity.entity_rank();

  bool ok = true ;
  std::ostringstream msg ;

  for ( PartVector::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part * const p = *i ;
    const unsigned part_rank = p->primary_entity_rank();

    const bool error_intersection = ! p->intersection_of().empty();
    const bool error_rel_target   = ! p->relations().empty() &&
                                    p == p->relations().begin()->m_target ;
    const bool error_rank = entity_rank != part_rank &&
                            undef_rank  != part_rank ;

    if ( error_intersection || error_rel_target || error_rank ) {
      if ( ok ) {
        ok = false ;
        msg << method ;
        msg << "( " ;
        print_entity_key( msg , meta , entity.key() );
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }

      msg << p->name() << "[" ;
      if ( part_rank < type_names.size() ) {
        msg << type_names[ part_rank ];
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

  if ( ! ok ) {
    msg << " }" ;
    throw std::runtime_error( msg.str() );
  }
}

}

void BulkData::change_entity_parts(
  Entity & e ,
  const std::vector<Part*> & add_parts ,
  const std::vector<Part*> & remove_parts )
{
  static const char method[] = "stk::mesh::BulkData::change_entity_parts" ;

  assert_ok_to_modify( method );

  assert_entity_owner( method , e , m_parallel_rank );

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  PartVector a_parts( add_parts );

  for ( PartVector::const_iterator
        ia = add_parts.begin(); ia != add_parts.end() ; ++ia ) {
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
    if ( m_mesh_meta_data.universal_part() == **ir )
      throw std::runtime_error ( "Cannot remove entity from universal part" );
    if ( m_mesh_meta_data.locally_owned_part() == **ir )
      throw std::runtime_error ( "Cannot remove entity from locally owned part" );
    if ( m_mesh_meta_data.globally_shared_part() == **ir )
      throw std::runtime_error ( "Cannot remove entity from globally shared part" );
    */

    if ( ! contain( a_parts , **ir ) ) {
      r_parts.push_back( *ir );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( e.bucket().member ( **cur_part ) )
          r_parts.push_back ( *cur_part );
    }
  }

  order( r_parts );

  verify_change_parts( method , m_mesh_meta_data , e , a_parts );
  verify_change_parts( method , m_mesh_meta_data , e , r_parts );

  internal_change_entity_parts( e , a_parts , r_parts );

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
  Entity & e ,
  const PartVector & add_parts ,
  const PartVector & remove_parts )
{
  // TODO 10-07-21 figure out how to solve this without using is_bucket_valid
  Bucket * const k_old = e.is_bucket_valid()
                         ? &(e.bucket())
                         : static_cast<Bucket*>(NULL);

  const unsigned i_old = e.bucket_ordinal() ;

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
        e.entity_rank(),
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
  m_entity_repo.change_entity_bucket( *k_new, e, k_new->size() );
  m_bucket_repository.add_entity_to_bucket( e, *k_new );

  // If changing buckets then remove the entity from the bucket,
  if ( k_old ) { m_bucket_repository.remove_entity( k_old , i_old ); }

  // Update the change counter to the current cycle.
  m_entity_repo.set_entity_sync_count( e, m_sync_count );

  // Propagate part changes through the entity's relations.

  internal_propagate_part_changes( e , parts_removed );
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity * & e )
{
  static const char method[] = "stk::mesh::BulkData::destroy_entity" ;

  Entity & entity = *e ;

  assert_ok_to_modify( method );

  bool has_upward_relation = false ;

  for ( PairIterRelation
        irel = entity.relations() ;
        ! irel.empty() && ! has_upward_relation ; ++irel ) {

    has_upward_relation = entity.entity_rank() <= irel->entity_rank();
  }

  if ( has_upward_relation ) { return false ; }

  if (  entity.marked_for_destruction() ) {
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

  m_bucket_repository.remove_entity( &(entity.bucket()) , entity.bucket_ordinal() );

  // Set the bucket to 'bucket_nil' which:
  //   1) has no parts at all
  //   2) has no field data
  //   3) has zero capacity
  //
  // This keeps the entity-bucket methods from catastrophically failing
  // with a bad bucket pointer.

  m_entity_repo.destroy_later( entity, m_bucket_repository.get_nil_bucket() );

  // Add destroyed entity to the transaction
  // m_transaction_log.delete_entity ( *e );

  // Set the calling entity-pointer to NULL;
  // hopefully the user-code will clean up any outstanding
  // references to this entity.

  e = NULL ;

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

      if ( result.second ) {
        // A new application-created entity
        m_entity_repo.set_entity_owner_rank( *(result.first), m_parallel_rank);
        m_entity_repo.set_entity_sync_count( *(result.first), m_sync_count);
      }
      else {
        //if an entity is declare with the declare_entity function in the same
        //modification cycle as the generate_new_entities function, and it happens to
        //generate a key that was declare previously in the same cycle it is an error
        std::ostringstream msg;
        msg << "stk::mesh::BulkData::generate_new_entities ERROR:";
        msg << " generated ";
        print_entity_key(msg, m_mesh_meta_data, key);
        msg << " which was already used in this modification cycle.";
        throw std::runtime_error(msg.str());
      }
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

