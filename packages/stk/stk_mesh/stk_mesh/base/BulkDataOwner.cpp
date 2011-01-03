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

void insert_closure_ghost( Entity * const entity ,
                           const unsigned proc_local ,
                           std::set<Entity*,EntityLess> & remove_list )
{
  if ( ! in_owned_closure( *entity , proc_local ) ) {
    // This entity is a ghost, put it on the remove_list
    // along with all ghosts in its closure

    std::pair< std::set<Entity*,EntityLess>::iterator , bool >
      result = remove_list.insert( entity );

    if ( result.second ) {
      // This ghost entity is new to the list, traverse its closure.

      const unsigned etype = entity->entity_rank();

      for ( PairIterRelation
            irel = entity->relations() ; ! irel.empty() ; ++irel ) {

        if ( irel->entity_rank() < etype ) {
          insert_closure_ghost( irel->entity() , proc_local ,remove_list );
        }
      }
    }
  }
}

void insert_transitive_ghost( Entity * const entity ,
                              const unsigned proc_local ,
                              std::set<Entity*,EntityLess> & remove_list )
{
  insert_closure_ghost( entity , proc_local , remove_list );

  // Transitive:
  // If this entity is a member of another entity's closure
  // then that other entity is part of the traversal.

  const unsigned etype = entity->entity_rank();

  for ( PairIterRelation rel = entity->relations(); ! rel.empty() ; ++rel ) {
    if ( etype < rel->entity_rank() ) {
      insert_transitive_ghost( rel->entity() , proc_local , remove_list );
    }
  }
}

//----------------------------------------------------------------------

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

  const unsigned etype = e.entity_rank();

  // Any higher ranking entities locally owned?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = etype  <  irel->entity_rank() &&
             p_rank == irel->entity()->owner_rank();
  }

  // Any higher ranking entity member of an owned closure?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = etype < irel->entity_rank() &&
             member_of_owned_closure( * irel->entity() , p_rank );
  }

  return result ;
}

//----------------------------------------------------------------------

void clean_and_verify_parallel_change(
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change )
{
  const MetaData      & meta   = mesh.mesh_meta_data() ;
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
        i = local_change.begin() ; i != local_change.end() ; ) {
    std::vector<EntityProc>::iterator j = i ; ++i ;
    Entity * const entity    = j->first ;
    const unsigned new_owner = j->second ;

    // Verification:
    // 1) If bucket has no capacity then is destined for deletion
    // 2) If not locally owned then not allowed grant ownership
    // 3) New owner must be legit
    // 4) Cannot grant to two different owners

    const bool bad_null   = NULL == entity ;
    const bool bad_delete = ! bad_null && EntityLogDeleted == entity->log_query();
    const bool bad_entity = ! bad_null && entity->owner_rank() != p_rank ;
    const bool bad_owner  = p_size <= new_owner ;
    const bool bad_dup    = ! bad_null && i != local_change.end() && entity == i->first ;

    if ( bad_null || bad_entity || bad_owner || bad_dup || bad_delete ) {
      ++error_count ;

      error_msg << "  P" << p_rank << ": " ;
      if ( bad_null ) { error_msg << " NULL ENTITY" ; }
      else { print_entity_key( error_msg , meta , entity->key() ); }
      if ( bad_delete ) { error_msg << " HAS_BEEN_DELETED" ; }
      if ( bad_entity ) { error_msg << " NOT_CURRENT_OWNER" ; }
      if ( bad_owner ) {
        error_msg << " BAD_NEW_OWNER( " << new_owner << " )" ;
      }
      if ( bad_dup ) {
        error_msg << " CONFLICTING_NEW_OWNER( " << new_owner ;
        error_msg << " != " << i->second << " )" ;
      }
      error_msg << std::endl ;
    }
    else if ( new_owner == p_rank ) {
      // Eliminate non-changes
      j->first = NULL ;
      j->second = 0 ;
    }
  }

  all_reduce( p_comm , ReduceSum<1>( & error_count ) );

  if ( error_count ) {
    all_write_string( p_comm , std::cerr , error_msg.str() );

    ThrowErrorMsg("Bad change ownership directives\n");
  }

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

  // Sizing:

  std::vector<unsigned> procs ;

  for ( std::vector<EntityProc>::const_iterator
        ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
    Entity & entity = * ip->first ;
    comm_procs( entity , procs );
    for ( std::vector<unsigned>::iterator
           j = procs.begin() ; j != procs.end() ; ++j ) {
      comm.send_buffer( *j ).skip<EntityKey>(1).skip<unsigned>(1);
    }
  }

  // Allocation:

  comm.allocate_buffers( p_size / 4 , 0 );

  // Packing new owner info:

  for ( std::vector<EntityProc>::const_iterator
        ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
    Entity & entity = * ip->first ;
    comm_procs( entity , procs );
    for ( std::vector<unsigned>::iterator
           j = procs.begin() ; j != procs.end() ; ++j ) {
      comm.send_buffer( *j )
          .pack<EntityKey>( entity.key() )
          .pack<unsigned>(  ip->second );
    }
  }

  comm.communicate();

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
  // Parallel synchronous determination of changing
  // shared and ghosted.

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

  // Closure of the owner change for impacted ghost entities.

  EntityProcSet send_closure ;

  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    insert_closure_send( *i , send_closure );
  }

  {
    EntitySet work ;

    for ( std::vector<EntityProc>::const_iterator
          i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , work );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , work );
    }

    for ( EntityProcSet::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      insert_transitive_ghost( i->first , m_parallel_rank , work );
    }

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty ;
    std::vector<Entity*> effected_ghosts( work.begin() , work.end() );

    // Skip 'm_ghosting[0]' which is the shared subset.
    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() + 1 ; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty , effected_ghosts );
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership

  {
    PartVector owned( 1 );
    owned[0] = & meta.locally_owned_part();

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

