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

  const size_t rank_count = meta.entity_type_count();

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
    m_buckets( mesh_meta_data.entity_type_count() ),
    m_entities(),
    m_entity_comm(),
    m_ghosting(),
    m_new_entities(),
    m_bucket_nil( NULL ),

    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_bucket_capacity( bucket_max_size ),
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

  // Remove entities from the buckets.
  // Destroy entities, which were allocated by the set itself.
  // Destroy buckets, which were *not* allocated by the set.

  try {
    for ( std::vector< std::vector<Bucket*> >::iterator
          i = m_buckets.end() ; i != m_buckets.begin() ; ) {
      try {
        std::vector<Bucket*> & kset = *--i ;

        while ( ! kset.empty() ) {
          try { Bucket::destroy_bucket( kset.back() ); } catch(...) {}
          kset.pop_back();
        }
        kset.clear();
      } catch(...) {}
    }
    m_buckets.clear();
  } catch(...) {}

  try { if ( m_bucket_nil ) Bucket::destroy_bucket( m_bucket_nil ); } catch(...) {}

  try {
    while ( ! m_entities.empty() ) {
      internal_expunge_entity( m_entities.begin() );
    }
  } catch(...){}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::update_field_data_states() const
{
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
  const size_t rank_count = m_mesh_meta_data.entity_type_count();
  const bool ok_id   = 0 < entity_id( key );
  const bool ok_type = entity_type( key ) < rank_count ;

  if ( ! ok_type || ! ok_id ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    if ( ! ok_type ) {
      msg << entity_type( key ) << "-"
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

  if ( ! m_meta_data_verified ) {
    m_mesh_meta_data.assert_committed( method );

    verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );

    m_meta_data_verified = true ;

    m_bucket_nil =
      Bucket::declare_nil_bucket(*this,m_mesh_meta_data.get_fields().size());
  }
  else {
    ++m_sync_count ;
  }

  parallel_machine_barrier( m_parallel_machine );

  if ( m_sync_state == MODIFIABLE ) return false ;

  m_sync_state = MODIFIABLE ;

  // Clear out the previous transaction information
  // m_transaction_log.flush();

  // Clear out the entities destroyed in the previous modification.
  // They were retained for change-logging purposes.

  for ( EntitySet::iterator i = m_entities.begin() ; i != m_entities.end() ; ) {
    const EntitySet::iterator j = i ; ++i ;
    if ( j->second->m_bucket == m_bucket_nil ) {
      internal_expunge_entity( j );
    }
  }

  return true ;
}

//----------------------------------------------------------------------

const std::vector<Bucket*> & BulkData::buckets( unsigned type ) const
{
  const char method[]= "stk::mesh::BulkData::buckets" ;

  m_mesh_meta_data.assert_entity_type( method , type );

  return m_buckets[ type ];
}

//----------------------------------------------------------------------

void BulkData::verify_type_and_id(const char* calling_method,
                                  EntityRank ent_type, EntityId ent_id) const
{
  m_mesh_meta_data.assert_entity_type( calling_method , ent_type );

  if (!entity_id_valid(ent_id)) {
    std::ostringstream msg;
    msg << calling_method << ": ent_id not valid";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }
}

Entity * BulkData::get_entity( EntityRank ent_type, EntityId ent_id ,
			       const char * /*required_by*/ ) const
{
  verify_type_and_id("BulkData::get_entity", ent_type, ent_id);

  EntityKey key(ent_type, ent_id);

  return get_entity(key);
}

Entity * BulkData::get_entity( EntityKey key) const
{
  const bool valid_key = entity_key_valid( key );

  const EntitySet::const_iterator i = m_entities.find( key );

  if ( ! valid_key ) {
    static const char method[] = "stk::mesh::BulkData::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , m_mesh_meta_data , key );
    msg << " INVALID KEY" ;
    msg << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  return i != m_entities.end() ? i->second : NULL ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.

std::pair<Entity*,bool>
BulkData::internal_create_entity( const EntityKey & key )
{
  EntitySet::value_type tmp(key,NULL);

  const std::pair< EntitySet::iterator , bool >
    insert_result = m_entities.insert( tmp );

  std::pair<Entity*,bool>
    result( insert_result.first->second , insert_result.second );

  if ( insert_result.second )  { // A new entity
    insert_result.first->second = result.first = new Entity( key );
    //result.first->m_owner_rank   = ~0u ;
    result.first->m_owner_rank = m_parallel_rank ;
    result.first->m_sync_count   = m_sync_count ;
  }

  return result ;
}


void BulkData::internal_expunge_entity( BulkData::EntitySet::iterator i )
{
  const bool ok_ptr = i->second != NULL ;
  const bool ok_key = ok_ptr ? i->first == i->second->key() : true ;

  if ( ! ok_ptr || ! ok_key ) {
    std::ostringstream msg ;
    msg << "stk::mesh::BulkData::internal_expunge_entity( " ;
    print_entity_key( msg , m_mesh_meta_data , i->first );
    if ( ! ok_ptr ) {
      msg << "NULL" ;
    }
    else {
      msg << " != " ;
      print_entity_key( msg , m_mesh_meta_data , i->second->key() );
    }
    msg << ") FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  delete i->second ;
  i->second = NULL ;
  m_entities.erase( i );
}

//----------------------------------------------------------------------

Entity & BulkData::declare_entity( EntityRank ent_type , EntityId ent_id ,
				   const std::vector<Part*> & parts )
{
  const char method[] = "stk::mesh::BulkData::declare_entity" ;

  assert_ok_to_modify( method );

  verify_type_and_id("BulkData::declare_entity", ent_type, ent_id);

  EntityKey key( ent_type , ent_id );

  assert_good_key( method , key );

  std::pair< Entity * , bool > result = internal_create_entity( key );

  if ( ! result.second ) {
    // An existing entity, the owner must match.
    assert_entity_owner( method , * result.first , m_parallel_rank );
  }

  if ( result.second ) { // A new entity
    m_new_entities.push_back( result.first );
    //result.first->m_owner_rank = m_parallel_rank ;
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
  const std::vector<std::string> & type_names = meta.entity_type_names();
  const unsigned undef_rank  = std::numeric_limits<unsigned>::max();
  const unsigned entity_rank = entity.entity_type();

  bool ok = true ;
  std::ostringstream msg ;

  for ( PartVector::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part * const p = *i ;
    const unsigned part_rank = p->primary_entity_type();

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
  const char method[] = "stk::mesh::BulkData::change_entity_parts" ;

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
    if ( m_mesh_meta_data.locally_used_part() == **ir )
      throw std::runtime_error ( "Cannot remove entity from locally used part" );
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
  Bucket * const k_old = e.m_bucket != m_bucket_nil
                       ? e.m_bucket : (Bucket *) NULL ;
  const unsigned i_old = e.m_bucket_ord ;


  if ( k_old && k_old->member_all( add_parts ) &&
              ! k_old->member_any( remove_parts ) ) {
    // Is already a member of all add_parts,
    // is not a member of any remove_parts,
    // thus nothing to do.
    return ;
  }

  // m_transaction_log.modify_entity ( e );
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
    Bucket::declare_bucket( *this ,
                            e.entity_type(),
                            parts_total.size(),
                            & parts_total[0] ,
                            m_bucket_capacity ,
                            m_mesh_meta_data.get_fields() ,
                            m_buckets[ e.entity_type() ] );

  // If changing buckets then copy its field values from old to new bucket

  if ( k_old ) {
    Bucket::copy_fields( *k_new , k_new->m_size , *k_old , i_old );
  }
  else {
    Bucket::zero_fields( *k_new , k_new->m_size );
  }

  // Set the new bucket
  e.m_bucket     = k_new ;
  e.m_bucket_ord = k_new->m_size ;
  k_new->m_entities[ k_new->m_size ] = & e ;
  ++( k_new->m_size );

  // If changing buckets then remove the entity from the bucket,
  if ( k_old ) { remove_entity( k_old , i_old ); }

  // Update the change counter to the current cycle.
  e.m_sync_count = m_sync_count ;

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

    has_upward_relation = entity.entity_type() <= irel->entity_rank();
  }

  if ( has_upward_relation ) { return false ; }

  if ( 0 == entity.m_bucket->capacity() ) {
    // Cannot already be destroyed.
    return false ;
  }
  // Add destroyed entity to the transaction

  //------------------------------
  // Immediately remove it from relations and buckets.
  // Postpone deletion until modification_end to be sure that
  // 1) No attempt is made to re-create it.
  // 2) Parallel index is cleaned up.
  // 3) Parallel sharing is cleaned up.
  // 4) Parallel ghosting is cleaned up.
  //
  // Must clean up the parallel lists before fully deleting the entity.

  while ( ! entity.m_relation.empty() ) {
    destroy_relation( entity , * entity.m_relation.back().entity() );
  }

  // m_transaction_log.delete_entity ( *e );

  remove_entity( entity.m_bucket , entity.m_bucket_ord );

  // Set the bucket to 'bucket_nil' which:
  //   1) has no parts at all
  //   2) has no field data
  //   3) has zero capacity
  //
  // This keeps the entity-bucket methods from catastrophically failing
  // with a bad bucket pointer.

  entity.m_bucket     = m_bucket_nil ;
  entity.m_bucket_ord = 0 ;

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
      std::pair<Entity *, bool> result = internal_create_entity(key);

      //if an entity is declare with the declare_entity function in the same
      //modification cycle as the generate_new_entities function, and it happens to
      //generate a key that was declare previously in the same cycle it is an error
      if (! result.second) {
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

void BulkData::remove_entity( Bucket * k , unsigned i )
{
  Bucket * const first = bucket_counter( k->m_key ) ? k->m_bucket : k ;
  Bucket * const last  = first->m_bucket ;

  // Only move if not the last entity being removed

  if ( last != k || k->m_size != i + 1 ) {

    // Not the same bucket or not the last entity

    // Copy last entity in last to ik slot i

    Entity * const entity = last->m_entities[ last->m_size - 1 ];

    Bucket::copy_fields( *k , i , *last , last->m_size - 1 );

    k->m_entities[i]     = entity ;
    entity->m_bucket     = k ;
    entity->m_bucket_ord = i ;

    // Entity field data has relocated

    internal_propagate_relocation( *entity );
  }

  --( last->m_size );

  if ( last->m_size != 0 ) {
    last->m_entities[ last->m_size ] = NULL ;
  }
  else {

    // The current 'last' bucket is to be deleted.
    // The previous 'last' bucket becomes the
    // new 'last' bucket in the family:

    std::vector<Bucket*> & bucket_set = m_buckets[ last->entity_type() ];

    std::vector<Bucket*>::iterator ik = lower_bound(bucket_set, last->m_key);

    if ( ik == bucket_set.end() || last != *ik ) {
      if ( ik == bucket_set.end() )
        std::cout << "Case 1 is met" << std::endl;
      if ( *ik != last )
        std::cout << "Case 2 is met" << std::endl;
      throw std::runtime_error(
        std::string("stk::mesh::BulkData::remove_entity INTERNAL FAILURE") );
    }

    ik = bucket_set.erase( ik );

    if ( first != last ) { first->m_bucket = *--ik ; }

    Bucket::destroy_bucket( last );
  }
}

//----------------------------------------------------------------------

void BulkData::internal_sort_bucket_entities()
{
  for ( unsigned entity_type = 0 ;
                 entity_type < m_buckets.size() ; ++entity_type ) {

    std::vector<Bucket*> & buckets = m_buckets[ entity_type ];

    size_t bk = 0 ; // Offset to first bucket of the family
    size_t ek = 0 ; // Offset to end   bucket of the family

    for ( ; bk < buckets.size() ; bk = ek ) {
      Bucket * ik_vacant = buckets[bk]->m_bucket ; // Last bucket, need space
      unsigned ie_vacant = ik_vacant->size();

      if ( ik_vacant->capacity() <= ie_vacant ) {
        // Have to create a bucket just for the scratch space...
        const unsigned * const bucket_key = buckets[bk]->m_key ;
        const unsigned         part_count = bucket_key[0] - 1 ;
        const unsigned * const part_ord   = bucket_key + 1 ;

        ik_vacant = Bucket::declare_bucket( *this , entity_type ,
                                            part_count , part_ord ,
                                            m_bucket_capacity ,
                                            m_mesh_meta_data.get_fields() ,
                                            buckets );
        ie_vacant = 0 ;
      }

      ik_vacant->m_entities[ ie_vacant ] = NULL ;

      // Determine offset to the end bucket in this family:
      while ( ek < buckets.size() && ik_vacant != buckets[ek] ) { ++ek ; }
      ++ek ;

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
          *j = b.m_entities[i] ;
        }
      }

      std::sort( entities.begin() , entities.end() , EntityLess() );

      j = entities.begin();

      bool change_this_family = false ;

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          Entity * const current = b.m_entities[i] ;

          if ( current != *j ) {

            if ( current ) {
              // Move current entity to the vacant spot
              Bucket::copy_fields( *ik_vacant , ie_vacant , b, i );
              current->m_bucket     = ik_vacant ;
              current->m_bucket_ord = ie_vacant ;
              ik_vacant->m_entities[ ie_vacant ] = current ;
            }

            // Set the vacant spot to where the required entity is now.
            ik_vacant = (*j)->m_bucket ;
            ie_vacant = (*j)->m_bucket_ord ;
            ik_vacant->m_entities[ ie_vacant ] = NULL ;

            // Move required entity to the required spot
            Bucket::copy_fields( b, i, *ik_vacant , ie_vacant );
            (*j)->m_bucket     = & b ;
            (*j)->m_bucket_ord = i ;
            b.m_entities[i]    = *j ;

            change_this_family = true ;
          }

          // Once a change has occured then need to propagate the
          // relocation for the remainder of the family.
          // This allows the propagation to be performed once per
          // entity as opposed to both times the entity is moved.

          if ( change_this_family ) { internal_propagate_relocation( **j ); }
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

