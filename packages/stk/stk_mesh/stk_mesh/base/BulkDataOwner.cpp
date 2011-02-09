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
#include <set>
#include <vector>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Bucket.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

// Given an entity, if it's a ghost, insert the closure of the entity
// into work_list.
void insert_closure_ghost( Entity * const entity ,
                           const unsigned proc_local ,
                           std::set<Entity*,EntityLess> & work_list )
{
  if ( ! in_owned_closure( *entity , proc_local ) ) {
    // This entity is a ghost, put it on the work_list
    // along with all ghosts in its closure

    std::pair< std::set<Entity*,EntityLess>::iterator , bool >
      result = work_list.insert( entity );

    if ( result.second ) {
      // This ghost entity is new to the list, traverse its closure.

      const unsigned entity_rank = entity->entity_rank();

      // Recurse over downward relations
      for ( PairIterRelation
            irel = entity->relations() ; ! irel.empty() ; ++irel ) {
        if ( irel->entity_rank() < entity_rank ) {
          insert_closure_ghost( irel->entity() , proc_local ,work_list );
        }
      }
    }
  }
}

// Given an entity, insert the closures of every entity that has this entity
// in its closure. Only ghosts will be inserted.
void insert_transitive_ghost( Entity * const entity ,
                              const unsigned proc_local ,
                              std::set<Entity*,EntityLess> & work_list )
{
  insert_closure_ghost( entity , proc_local , work_list );

  // Transitive:
  // If this entity is a member of another entity's closure
  // then that other entity is part of the traversal.

  const unsigned entity_rank = entity->entity_rank();

  // Recurse over upward relations
  for ( PairIterRelation rel = entity->relations(); ! rel.empty() ; ++rel ) {
    if ( entity_rank < rel->entity_rank() ) {
      insert_transitive_ghost( rel->entity() , proc_local , work_list );
    }
  }
}

//----------------------------------------------------------------------

// Add EntityProc pairs to send_list for every entity in the closure of the
// entity in send_entry. All these entities will be sent to the same proc as
// the original send_entry.
void insert_closure_send(
  const EntityProc                  send_entry ,
  std::set<EntityProc,EntityLess> & send_list )
{
  ThrowRequireMsg( send_entry.first->log_query() != EntityLogDeleted,
      "Cannot send destroyed entity " << print_entity_key(send_entry.first));

  std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
    result = send_list.insert( send_entry );

  if ( result.second ) {
    // First time this entity was inserted into the send_list.

    const unsigned erank  = send_entry.first->entity_rank();
    PairIterRelation irel = send_entry.first->relations();

    // Recurse over downward relations
    for ( ; ! irel.empty() ; ++irel ) {
      if ( irel->entity_rank() < erank ) {
        const EntityProc rel_send_entry( irel->entity(), send_entry.second );

        insert_closure_send( rel_send_entry , send_list );
      }
    }
  }
}

//----------------------------------------------------------------------

bool member_of_owned_closure( const Entity & e , const unsigned p_rank )
{
  bool result = p_rank == e.owner_rank();

  const unsigned entity_rank = e.entity_rank();

  // Any higher ranking entities locally owned?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = entity_rank  <  irel->entity_rank() &&
             p_rank == irel->entity()->owner_rank();
  }

  // Any higher ranking entity member of an owned closure?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = entity_rank < irel->entity_rank() &&
             member_of_owned_closure( * irel->entity() , p_rank );
  }

  return result ;
}

//----------------------------------------------------------------------

// Given a vector of local ownership changes, remove duplicates and
// sanity check.
void clean_and_verify_parallel_change(
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change )
{
  const MetaData      & meta   = MetaData::get(mesh);
  const unsigned        p_rank = mesh.parallel_rank();
  const unsigned        p_size = mesh.parallel_size();
  const ParallelMachine p_comm = mesh.parallel();

  size_t error_count = 0 ;

  std::ostringstream error_msg ;

  // Order and eliminate redundancies:
  {
    std::vector<EntityProc>::iterator i = local_change.begin() ,
                                      j = local_change.end() ;
    std::sort( i , j , EntityLess() );
    i = std::unique( i , j );
    local_change.erase( i , j );
  }

  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    std::vector<EntityProc>::iterator next = i+1 ;
    Entity * const entity    = i->first ;
    const unsigned new_owner = i->second ;

    // Verification:
    // 1) Cannot change the ownership of an entity that you've already marked as deleted
    // 2) Cannot change the ownership of an entity you do not own
    // 3) New owner must exist
    // 4) Cannot grant ownership to two different owners

    const bool bad_null = NULL == entity ;

    // Cannot change the ownership of an entity that you've already marked as deleted
    const bool bad_marked_deleted = ! bad_null && EntityLogDeleted == entity->log_query();

    // Cannot change the ownership of an entity you do not own
    const bool bad_process_not_entity_owner = ! bad_null && entity->owner_rank() != p_rank ;

    // New owner must exist
    const bool bad_new_owner_does_not_exist = p_size <= new_owner ;

    // Cannot grant ownership to two different owners
    const bool bad_inconsistent_change = ! bad_null && next != local_change.end() && entity == next->first ;

    if ( bad_null ||
         bad_process_not_entity_owner ||
         bad_new_owner_does_not_exist ||
         bad_inconsistent_change ||
         bad_marked_deleted)
    {
      ++error_count ;

      error_msg << "  P" << p_rank << ": " ;
      if ( bad_null ) { error_msg << " NULL ENTITY" ; }
      else { print_entity_key( error_msg , meta , entity->key() ); }
      if ( bad_marked_deleted ) { error_msg << " HAS_BEEN_DELETED" ; }
      if ( bad_process_not_entity_owner ) { error_msg << " NOT_CURRENT_OWNER" ; }
      if ( bad_new_owner_does_not_exist ) {
        error_msg << " BAD_NEW_OWNER( " << new_owner << " )" ;
      }
      if ( bad_inconsistent_change ) {
        error_msg << " CONFLICTING_NEW_OWNER( " << new_owner ;
        error_msg << " != " << next->second << " )" ;
      }
      error_msg << std::endl ;
    }
    else if ( new_owner == p_rank ) {
      // Eliminate non-changes
      i->first = NULL ;
      i->second = 0 ;
    }
  }

  all_reduce( p_comm , ReduceSum<1>( & error_count ) );

  if ( error_count ) {
    all_write_string( p_comm , std::cerr , error_msg.str() );

    ThrowErrorMsg("Bad change ownership directives\n");
  }

  // Filter out non-changes (entity will be NULL
  {
    std::vector<EntityProc>::iterator i = local_change.begin(),
                                      j = local_change.end();
    i = std::remove( i , j , EntityProc(NULL,0) );
    local_change.erase( i , j );
  }
}

//----------------------------------------------------------------------
// Generate a parallel consistent list of ownership changes:
// 1) Shared entities (not owned but in closure of an owned entity),
// 2) Ghosted entities (not owned and not in closure of an owned entity), and
// 3) Parallel index.

void generate_parallel_change( const BulkData & mesh ,
                               const std::vector<EntityProc> & local_change ,
                                     std::vector<EntityProc> & shared_change ,
                                     std::vector<EntityProc> & ghosted_change )
{
  const unsigned p_size = mesh.parallel_size();

  CommAll comm( mesh.parallel() );

  std::vector<unsigned> procs ;

  // pack and communicate change owner information to all
  // processes that know about the entity
  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::vector<EntityProc>::const_iterator
          ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
      Entity & entity    = * ip->first ;
      unsigned new_owner = ip->second;
      comm_procs( entity , procs );
      for ( std::vector<unsigned>::iterator
            j = procs.begin() ; j != procs.end() ; ++j )
      {
        comm.send_buffer( *j )
          .pack<EntityKey>( entity.key() )
          .pack<unsigned>(  new_owner );
      }
    }
    if (phase == 0) { // allocation phase
      comm.allocate_buffers( p_size / 4 , 0 );
    }
    else { // communication phase
      comm.communicate();
    }
  }

  // unpack communicated owner information into the
  // ghosted and shared change vectors.
  for ( unsigned ip = 0 ; ip < p_size ; ++ip ) {
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      EntityProc entry ;
      EntityKey key ;
      buf.unpack<EntityKey>( key )
         .unpack<unsigned>( entry.second );

      entry.first = mesh.get_entity( key );

      if ( in_receive_ghost( * entry.first ) ) {
        ghosted_change.push_back( entry );
      }
      else {
        shared_change.push_back( entry );
      }
    }
  }

  std::sort( shared_change.begin() , shared_change.end() , EntityLess() );
  std::sort( ghosted_change.begin() , ghosted_change.end() , EntityLess() );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::change_entity_owner( const std::vector<EntityProc> & arg_change )
{
  const MetaData  & meta = m_mesh_meta_data ;
  const unsigned  p_rank = m_parallel_rank ;
  const unsigned  p_size = m_parallel_size ;
  ParallelMachine p_comm = m_parallel_machine ;

  //------------------------------
  // Verify the input changes, generate a clean local change list, and
  // generate the remote change list so that all processes know about
  // pending changes.

  std::vector<EntityProc> local_change( arg_change );

  // Parallel synchronous clean up and verify the requested changes:
  clean_and_verify_parallel_change( *this , local_change );

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  generate_parallel_change( *this , local_change ,
                            shared_change , ghosted_change );

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  typedef std::set<EntityProc,EntityLess> EntityProcSet;
  typedef std::set<Entity*,EntityLess> EntitySet;

  // Compute the closure of all the locally changing entities
  EntityProcSet send_closure ;
  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    insert_closure_send( *i , send_closure );
  }

  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
    EntitySet modified_ghosts ;

    for ( std::vector<EntityProc>::const_iterator
          i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , modified_ghosts );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , modified_ghosts );
    }

    for ( EntityProcSet::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , modified_ghosts );
    }

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty_add ;
    std::vector<Entity*> remove_modified_ghosts( modified_ghosts.begin() ,
                                                 modified_ghosts.end() );

    // Skip 'm_ghosting[0]' which is the shared subset.
    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() + 1; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty_add , remove_modified_ghosts );
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership

  {
    PartVector owned;
    owned.push_back(& meta.locally_owned_part());

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      change_entity_parts( * i->first , PartVector() , owned );

      m_entity_repo.set_entity_owner_rank( *(i->first), i->second);
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      m_entity_repo.set_entity_owner_rank( *(i->first), i->second);
      if ( p_rank == i->second ) { // I receive ownership
        change_entity_parts( * i->first , owned , PartVector() );
      }
    }
  }

  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    CommAll comm( p_comm );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity     & entity = * i->first ;
      pack_entity_info( buffer , entity );
      pack_field_values( buffer , entity );
    }

    comm.allocate_buffers( p_size / 4 );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity     & entity = * i->first ;
      pack_entity_info( buffer , entity );
      pack_field_values( buffer , entity );
    }

    comm.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        unsigned  owner = ~0u ;

        unpack_entity_info( buf, *this, key, owner, parts, relations );

        // Received entity information will be correct,
        // modulo the owned and shared parts

        remove( parts , meta.globally_shared_part() );

        if ( owner == p_rank ) {
          // Must have the locally_owned_part
          insert( parts , meta.locally_owned_part() );
        }
        else {
          // Must not have the locally_owned_part
          remove( parts , meta.locally_owned_part() );
        }

        std::pair<Entity*,bool> result =
          m_entity_repo.internal_create_entity( key );

        m_entity_repo.log_created_parallel_copy( *(result.first) );

        // The entity was copied and not created.

        m_entity_repo.set_entity_owner_rank( *(result.first), owner);

        internal_change_entity_parts( *result.first , parts , PartVector() );

        declare_relation( *result.first , relations );

        if ( ! unpack_field_values( buf , * result.first , error_msg ) ) {
          ++error_count ;
        }
      }
    }

    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
    ThrowErrorMsgIf( error_count, error_msg.str() );

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.

    {
      Entity * entity = NULL ;

      for ( std::set<EntityProc,EntityLess>::iterator
            i = send_closure.end() ; i != send_closure.begin() ; ) {

        Entity * e = (--i)->first ;

        // The same entity may be sent to more than one process.
        // Only evaluate it once.

        if ( entity != e ) {
          entity = e ;
          if ( ! member_of_owned_closure( *e , p_rank ) ) {
            ThrowRequireMsg( destroy_entity( e ),
                "Failed to destroy entity " << print_entity_key(e) );
          }
        }
      }
    }

    send_closure.clear(); // Has been invalidated
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

