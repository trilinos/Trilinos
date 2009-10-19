/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

BulkData::BulkData( const MetaData & mesh_meta_data ,
            ParallelMachine parallel ,
             unsigned bucket_max_size )
  : m_buckets(),
    m_entities(),
    m_shares_all(),
    m_ghosting(),
    m_entities_owner_index(),
    m_new_entities(),
    m_del_entities(),
    m_bucket_nil( NULL ),

    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_bucket_capacity( bucket_max_size ),
    m_sync_count( 0 ),
    m_sync_state( false )
{
  static const char method[] = "stk::mesh::BulkData::Mesh" ;

  m_mesh_meta_data.assert_committed( method );

  verify_parallel_consistency( mesh_meta_data , parallel );

  m_buckets.resize( m_mesh_meta_data.entity_type_count() );

  m_bucket_nil =
    Bucket::declare_nil_bucket( *this , mesh_meta_data.get_fields().size() );

  create_ghosting( std::string("shared_aura") );
}

BulkData::~BulkData()
{
  try {
    while ( ! m_ghosting.empty() ) {
      delete m_ghosting.back();
      m_ghosting.pop_back();
    }
  } catch(...){}

  try { m_shares_all.clear(); } catch(...){}

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

  try { Bucket::destroy_bucket( m_bucket_nil ); } catch(...) {}

  try {
    for ( EntitySet::iterator
          it = m_entities.begin(); it != m_entities.end(); ++it) {
      it->second->m_bucket     = NULL ;
      it->second->m_bucket_ord = 0 ;
      delete it->second;
    }
    m_entities.clear();
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
  if ( m_sync_state ) {
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
  const bool error_destroyed = NULL != e.m_bucket &&
                                  0 == e.m_bucket->m_capacity ;

  if ( error_not_owner || error_destroyed ) {
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , m_mesh_meta_data , e.key() );
    msg << " ) FAILED" ;

    if ( error_not_owner ) {
      msg << " : Owner( " << e.owner_rank()
          << " ) != Required( " << owner << " )" ;
    }

    if ( error_destroyed ) {
       msg << " : Entity has been destroyed" ;
    }
    throw std::runtime_error( msg.str() );
  }
}

void BulkData::assert_good_key( const char * method ,
                                const EntityKey & key ) const
{
  const bool bad_key  = ! entity_key_valid( key );
  const bool bad_type = m_buckets.size() <= entity_type( key );

  if ( bad_key || bad_type ) { 
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    if ( bad_type ) {
      msg << entity_type( key ) << "-"
          << entity_id( key ) << " : BAD KEY TYPE" ;
    }
    else {
      print_entity_key( msg , m_mesh_meta_data , key );
      if ( bad_key ) { msg << " : BAD KEY" ; }
    }
    msg << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_begin()
{
  parallel_machine_barrier( m_parallel_machine );

  if ( ! m_sync_state ) return false ;

  m_sync_state = false ;

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
                                  EntityType ent_type, EntityId ent_id) const
{
  m_mesh_meta_data.assert_entity_type( calling_method , ent_type );

  if (!entity_id_valid(ent_id)) {
    std::ostringstream msg;
    msg << calling_method << ": ent_id not valid";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }
}

Entity * BulkData::get_entity( EntityType ent_type, EntityId ent_id ,
			       const char * required_by ) const
{
  verify_type_and_id("BulkData::get_entity", ent_type, ent_id);

  EntityKey key(ent_type, ent_id);

  const bool valid_key = entity_key_valid( key );

  const EntitySet::const_iterator i = m_entities.find( key );

  if ( ! valid_key || ( required_by && i == m_entities.end() ) ) {
    static const char method[] = "stk::mesh::BulkData::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , m_mesh_meta_data , key );
    if ( valid_key ) { msg << " , " << required_by ; }
    else { msg << " INVALID KEY" ; }
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
BulkData::internal_create_entity( const EntityKey & key ,
                                  const unsigned    owner )
{
  const char method[] = "stk::mesh::BulkData::internal_create_entity" ;

  EntitySet::value_type tmp(key,NULL);

  const std::pair< EntitySet::iterator , bool >
    insert_result = m_entities.insert( tmp );

  std::pair<Entity*,bool>
    result( insert_result.first->second , insert_result.second );

  if ( insert_result.second )  { // A new entity
    insert_result.first->second = result.first = new Entity( key );
    result.first->m_owner_rank = owner ;
  }
  else { // An existing entity, the owner must match.
    assert_entity_owner( method , * result.first , owner );
  }

  return result ;
}

//----------------------------------------------------------------------

Entity & BulkData::declare_entity( EntityType ent_type , EntityId ent_id ,
				   const std::vector<Part*> & parts )
{
  const char method[] = "stk::mesh::BulkData::declare_entity" ;

  assert_ok_to_modify( method );

  verify_type_and_id("BulkData::declare_entity", ent_type, ent_id);

  EntityKey key( ent_type , ent_id );

  assert_good_key( method , key );

  std::pair< Entity * , bool > result =
    internal_create_entity( key , m_parallel_rank );

  if ( result.second ) { // A new entity
    m_new_entities.push_back( result.first );
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  std::vector<Part*> rem ;
  std::vector<Part*> add( parts );
  add.push_back( owns );

  change_entity_parts( * result.first , add , rem );

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
    if ( ! contain( a_parts , **ir ) ) {
      r_parts.push_back( *ir );
      r_parts.insert( r_parts.end(), (*ir)->subsets().begin(),
                                     (*ir)->subsets().end() );
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
  Bucket * const k_old = e.m_bucket ;
  const unsigned i_old = e.m_bucket_ord ;

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

    has_upward_relation = entity.entity_type() <= irel->entity_type();
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

  while ( ! entity.m_relation.empty() ) {
    destroy_relation( entity , * entity.m_relation.back().entity() );
  }

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

  // Remember what entities have been destroyed for
  // modification_end_syncronize clean up
  m_del_entities.push_back( e );

  // Set the calling entity-pointer to NULL;
  // hopefully the user-code did not keep any outstanding
  // references to this entity.

  e = NULL ;

  return true ;
}

//----------------------------------------------------------------------

void BulkData::internal_destroy_entity( Entity * e )
{
  while ( ! e->m_relation.empty() ) {
    destroy_relation( * e , * e->m_relation.back().entity() );
  }

  if ( e->m_bucket->capacity() ) { // Not the 'nil' bucket.
    remove_entity( e->m_bucket , e->m_bucket_ord );
  }

  e->m_bucket     = NULL ;
  e->m_bucket_ord = 0 ;

  m_entities.erase( e->key() );

  delete e ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

