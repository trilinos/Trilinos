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

//----------------------------------------------------------------------

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

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

bool comm_mesh_verify_parallel_consistency(
  BulkData & M , std::ostream & error_log );

//----------------------------------------------------------------------

unsigned BulkData::determine_new_owner( Entity & entity ) const
{
  const bool keeping = entity.bucket().capacity();

  unsigned new_owner = keeping ? m_parallel_rank : ~0u ;

  for ( PairIterEntityComm
        share = entity.sharing(); ! share.empty() ; ++share ) {
    if ( share->proc < m_parallel_size &&
         ( new_owner < share->proc || m_parallel_size <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------
// Send to shared and ghost-receive,

bool send_to_shared_and_ghost_recv( const BulkData                & mesh ,
                                    const std::vector<Entity*>    & send ,
                                          std::vector<EntityProc> & recv )
{
  const unsigned p_size = mesh.parallel_size();
  const unsigned p_rank = mesh.parallel_rank();

  CommAll comm( mesh.parallel() );

  std::vector<unsigned> procs ;

  bool local = false ;
  for ( std::vector<Entity*>::const_iterator
        i = send.begin(); i != send.end() ; ++i ) {
    Entity & entity = **i ;
    if ( p_rank == entity.owner_rank() || in_shared( entity ) ) {
      comm_procs( entity , procs );
      for ( std::vector<unsigned>::iterator
            ip = procs.begin() ; ip != procs.end() ; ++ip ) {
        comm.send_buffer( *ip ).pack<EntityKey>( entity.key() );
        local = true ;
      }
    }
  }

  const bool global = comm.allocate_buffers( p_size/4, false, local );

  if ( ! global ) { return false ; }

  for ( std::vector<Entity*>::const_iterator
        i = send.begin(); i != send.end() ; ++i ) {
    Entity & entity = **i ;
    if ( p_rank == entity.owner_rank() || in_shared( entity ) ) {
      comm_procs( entity , procs );
      for ( std::vector<unsigned>::iterator
            ip = procs.begin() ; ip != procs.end() ; ++ip ) {
        comm.send_buffer( *ip ).pack<EntityKey>( entity.key() );
      }
    }
  }

  comm.communicate();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm.recv_buffer( p );

    EntityProc ep ; ep.second = p ;

    while ( buf.remaining() ) {
      EntityKey key ; buf.unpack<EntityKey>( key );

      ep.first = mesh.get_entity( key );

      recv.push_back( ep );
    }
  }

  std::sort( recv.begin() , recv.end() , EntityLess() );

  return true ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_update_parallel_index(
  const std::vector<Entity*> & del_entities ,
        std::vector<Entity*> & shared_new )
{
  static const char method[] =
    "stk::mesh::BulkData::internal_update_parallel_index" ;

  const size_t new_count = m_new_entities.size();
  const size_t del_count =   del_entities.size();

  std::vector< parallel::DistributedIndex::KeyType >
    new_entities_keys(new_count),
    del_entities_keys(del_count);

  std::vector< parallel::DistributedIndex::KeyProc >
    new_entities_keyprocs ;

  for ( size_t i = 0 ; i < new_count ; ++i ) {
    new_entities_keys[i] = m_new_entities[i]->key().raw_key();
  }

  for ( size_t i = 0 ; i < del_count ; ++i ) {
    del_entities_keys[i] = del_entities[i]->key().raw_key();
  }

  //------------------------------
  // Verify created did not already exist:
  {
    m_entities_index.query( new_entities_keys , new_entities_keyprocs );

    unsigned error_count = new_entities_keyprocs.size();

    all_reduce( m_parallel_machine , ReduceSum<1>( & error_count ) );

    if ( error_count ) {
      std::ostringstream msg ;
      msg << method ;
      msg << " FAILED: Declare entities already existed" ;

      for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
            i =  new_entities_keyprocs.begin() ;
            i != new_entities_keyprocs.end() ; ++i ) {
        EntityKey key( & i->first );
        msg << std::endl << "  ( " ;
        print_entity_key( msg , m_mesh_meta_data , key );
        msg << "  Already exists on P" << i->second << " )" ;
      }
 
      throw std::runtime_error( msg.str() );
    }
  }

  //------------------------------
  // Update distributed index.

  m_entities_index.update_keys( new_entities_keys , del_entities_keys );

  //------------------------------
  // Inform creating processes:

  m_entities_index.query( new_entities_keys , new_entities_keyprocs );

  {
    Entity * entity = NULL ;

    for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
          i =  new_entities_keyprocs.begin() ;
          i != new_entities_keyprocs.end() ; ++i ) {

      EntityKey key( & i->first );

      if ( (int) m_parallel_rank != i->second ) {
        // Another process also created this entity.

        if ( entity == NULL || entity->key() != key ) {
          // Have not looked this entity up by key 
          const EntitySet::iterator j = m_entities.find( key );
          entity = j->second ;

          shared_new.push_back( entity );
        }

        // Add the other_process to the entity's sharing info.
        entity->insert( EntityCommInfo( 0 , i->second ) );
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_resolve_parallel_create_delete(
  const std::vector<Entity*> & del_entities )
{
  // 'del_entities' is guaranteed unique and sorted

  {
    std::vector<int> local_flags( m_ghosting.size() , 0 );
    std::vector<int> global_flags( m_ghosting.size() , 0 );
    std::vector<EntityProc> del_entities_remote ;
 
    bool global_delete_flag =
      send_to_shared_and_ghost_recv( *this, del_entities, del_entities_remote );

    if ( global_delete_flag ) {
      // At least one process has deleted an entity that is
      // not owned, shared, or ghosted.

      // Clear communication lists of remotely deleted entities:

      for ( ; ! del_entities_remote.empty() ; del_entities_remote.pop_back() ) {
        Entity *       entity = del_entities_remote.back().first ;
        const unsigned proc   = del_entities_remote.back().second ;
        const bool     owner  = entity->owner_rank() == proc ;
        const bool     shared = in_shared( *entity , proc );
        const bool     g_recv = in_receive_ghost( *entity );
        const bool     g_send = in_send_ghost( *entity , proc );

        if ( shared ) {

          // A shared entity is being deleted on the remote process.
          // Remove it from the sharing.

          for ( std::vector<EntityCommInfo>::iterator
                j = entity->m_comm.begin();
                j != entity->m_comm.end() && j->ghost_id == 0 ; ++j ) {
            if ( j->proc == proc ) {
              entity->m_comm.erase( j );
              break ;
            }
          }

          if ( owner ) {
            // Remote owner is deleting a shared entity.
            // Have to determine new owner from remaining sharing processes.

            const unsigned new_owner = determine_new_owner( *entity );

            entity->m_owner_rank = new_owner ;
            entity->m_sync_count = m_sync_count ;

            if ( new_owner == m_parallel_rank ) {
              // Changing remotely owned to locally owned
              PartVector add( 1 );
              add[0] = & m_mesh_meta_data.locally_owned_part();
              internal_change_entity_parts( *entity , add , PartVector() );
            }
          }
        }

        if ( g_send ) {
          // Remotely ghosted entity is being destroyed,
          // remove from ghosting list
          for ( std::vector<EntityCommInfo>::iterator
                j = entity->m_comm.end(); j != entity->m_comm.begin() ; ) {
            --j ;
            if ( j->ghost_id != 0 && j->proc == proc ) {
              local_flags[ j->ghost_id ] = 1 ;
              j = entity->m_comm.erase( j );
            }
          }
        }
        else if ( owner && g_recv ) {
          // Remotely owned entity is being destroyed.
          // This receive ghost then it must be destroyed.
          while ( ! entity->m_comm.empty() ) {
            local_flags[ entity->m_comm.back().ghost_id ] = 1 ;
            entity->m_comm.pop_back();
          }
          destroy_entity( entity );
        }
      }

      // Clear communication list of locally deleted entities:

      for ( std::vector<Entity*>::const_iterator
            i = del_entities.end() ; i != del_entities.begin() ; ) {

        Entity * entity = *--i ;
        const bool shared = in_shared( *entity );

        if ( shared ) {
          const unsigned new_owner = determine_new_owner( *entity );

          entity->m_owner_rank = new_owner ;
          entity->m_sync_count = m_sync_count ;
        }

        while ( ! entity->m_comm.empty() ) {
          local_flags[ entity->m_comm.back().ghost_id ] = 1 ;
          entity->m_comm.pop_back();
        }
      }

      all_reduce_sum( m_parallel_machine ,
                      & local_flags[0] ,
                      & global_flags[0] , local_flags.size() );

      for ( unsigned ic = 0 ; ic < local_flags.size() ; ++ic ) {
        if ( global_flags[ic] ) {
          m_ghosting[ic]->m_sync_count = m_sync_count ;
        }
      }

      // Clean up the entity communication list

      {
        std::vector<Entity*>::iterator i = m_entity_comm.begin();
        bool changed = false ;
        for ( ; i != m_entity_comm.end() ; ++i ) {
          if ( (*i)->comm().empty() ) { *i = NULL ; changed = true ; }
        }
        if ( changed ) {
          i = std::remove( m_entity_comm.begin() ,
                           m_entity_comm.end() , (Entity *) NULL );
          m_entity_comm.erase( i , m_entity_comm.end() );
        }
      }
    }
  }

  //------------------------------

  {
    std::vector<Entity*> shared_created ;

    // Update the parallel index and
    // output sharing of created entities.
    internal_update_parallel_index( del_entities , shared_created );

    // Update shared created entities.
    // - Revise ownership to selected processor
    // - Update sharing.

    Part * const owns_part = & m_mesh_meta_data.locally_owned_part();
     
    const std::vector<Part*> add_parts , remove_parts( 1 , owns_part );
     
    for ( std::vector<Entity*>::const_iterator 
          ip = shared_created.begin() ; ip != shared_created.end() ; ++ip ) {
     
      Entity * entity = *ip ;

      const unsigned new_owner = determine_new_owner( *entity );

      entity->m_owner_rank = new_owner ;

      if ( new_owner != m_parallel_rank ) {
        // Removing locally owned.
        entity->m_sync_count = m_sync_count ;
        internal_change_entity_parts( *entity , add_parts , remove_parts );
      }
    }

    const size_t n_old = m_entity_comm.size();
     
    m_entity_comm.insert( m_entity_comm.end() ,  
                          shared_created.begin() , shared_created.end() );
     
    std::inplace_merge( m_entity_comm.begin() ,  
                        m_entity_comm.begin() + n_old ,
                        m_entity_comm.end() ,
                        EntityLess() );
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_end()
{
  return internal_modification_end( true );
}

bool BulkData::internal_modification_end( bool regenerate_aura )
{
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  std::vector<Entity*> del_entities ;

  for ( EntitySet::iterator
        i = m_entities.begin() ; i != m_entities.end() ; ++i ) {
    if ( i->second->m_bucket == m_bucket_nil ) {
      del_entities.push_back( i->second );
    }
  }

  int local_change_count[2] ;
  int global_change_count[2] ;
  local_change_count[0] = m_new_entities.size();
  local_change_count[1] =   del_entities.size();

  all_reduce_sum( m_parallel_machine ,
                  local_change_count , global_change_count , 2 );

  if ( global_change_count[0] || global_change_count[1] ) {

    if ( 1 < m_parallel_size ) {
      internal_resolve_parallel_create_delete( del_entities );
    }

    m_new_entities.clear();
  }

  // Parallel distribution considerations:
  if ( 1 < m_parallel_size ) {

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura ) { internal_regenerate_shared_aura(); }

    // Resolve changed part membership:
    internal_resolve_shared_membership();

    {
      std::ostringstream msg ;
      if ( ! comm_mesh_verify_parallel_consistency( *this , msg ) ) {
        std::cerr << msg.str() << std::endl ;
        throw std::runtime_error( msg.str() );
      }
    }
  }
  // Parallel consistency is restored and verified.
  // ------------------------------
  // The very last operation performed is to sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.

  internal_sort_bucket_entities();

  // ------------------------------

  m_sync_state = SYNCHRONIZED ;

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { OWNERSHIP_PART_COUNT = 3 };

namespace {

void pack_induced_memberships( CommAll & comm ,
                               const std::vector<Entity*> & entity_comm )
{
  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( in_shared( entity , entity.owner_rank() ) ) {
      // Is shared with owner, send to owner.

      PartVector empty , induced ;

      induced_part_membership( entity , empty , induced );

      CommBuffer & buf = comm.send_buffer( entity.owner_rank() );

      unsigned tmp = induced.size();

      buf.pack<unsigned>( tmp );

      for ( PartVector::iterator
            j = induced.begin() ; j != induced.end() ; ++j ) {
        tmp = (*j)->mesh_meta_data_ordinal();
        buf.pack<unsigned>( tmp );
      }
    }
  }
}

void generate_send_list( const size_t sync_count ,
                         const unsigned p_rank ,
                         const std::vector<Entity*>    & entity_comm ,
                               std::vector<EntityProc> & send_list )
{
  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( entity.owner_rank() == p_rank &&
         entity.synchronized_count() == sync_count ) {

      for ( PairIterEntityComm ec = entity.comm() ; ! ec.empty() ; ++ec ) {
        EntityProc tmp( & entity , ec->proc );
        send_list.push_back( tmp );
      }
    }
  }

  {
    std::sort( send_list.begin() , send_list.end() , EntityLess() );
    std::vector<EntityProc>::iterator i =
      std::unique( send_list.begin() , send_list.end() );
    send_list.erase( i , send_list.end() );
  }
}

void pack_memberships( CommAll & comm ,
                       const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity & entity = * i->first ;

    std::pair<const unsigned *, const unsigned *>
      part_ord = entity.bucket().superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    const unsigned count =
      ( part_ord.second - part_ord.first ) - OWNERSHIP_PART_COUNT ;

    const unsigned * const start = part_ord.first + OWNERSHIP_PART_COUNT ;

    comm.send_buffer( i->second ).pack<EntityKey>( entity.key() )
                                 .pack<unsigned>( count )
                                 .pack<unsigned>( start , count );
  }
}

}

//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.

void BulkData::internal_resolve_shared_membership()
{
  static const char method[] =
    "stk::mesh::BulkData::internal_resolve_shared_membership" ;

  const MetaData & meta  = m_mesh_meta_data ;
  ParallelMachine p_comm = m_parallel_machine ;
  const unsigned  p_rank = m_parallel_rank ;
  const unsigned  p_size = m_parallel_size ;
  const PartVector & all_parts = meta.get_parts();

  // Quick verification of part ordinal assumptions

  if ( OWNERSHIP_PART_COUNT <= meta.universal_part().mesh_meta_data_ordinal() ||
       OWNERSHIP_PART_COUNT <= meta.locally_used_part().mesh_meta_data_ordinal() ||
       OWNERSHIP_PART_COUNT <= meta.locally_owned_part().mesh_meta_data_ordinal() ) {
    throw std::logic_error( std::string( method ) );
  }

  //  Shared entities may have been modified due to relationship changes.
  //  Send just the current induced memberships from the sharing to
  //  the owning processes.
  {
    CommAll comm( p_comm );

    pack_induced_memberships( comm , m_entity_comm );

    comm.allocate_buffers( p_size / 4 );

    pack_induced_memberships( comm , m_entity_comm );

    comm.communicate();

    for ( std::vector<Entity*>::iterator
          i = m_entity_comm.begin() ; i != m_entity_comm.end() ; ++i ) {

      Entity & entity = **i ;

      if ( entity.owner_rank() == p_rank ) {
        // Receiving from all sharing processes

        PartVector empty , induced_parts , current_parts , remove_parts ;

        induced_part_membership( entity , empty , induced_parts );

        for ( PairIterEntityComm
              ec = entity.sharing() ; ! ec.empty() ; ++ec ) {

          CommBuffer & buf = comm.recv_buffer( ec->proc );

          unsigned count = 0 ; buf.unpack<unsigned>( count );
          for ( unsigned j = 0 ; j < count ; ++j ) {
            unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
            insert( induced_parts , * all_parts[ part_ord ] );
          }
        }

        // Remove any part that is an induced part but is not
        // in the induced parts list.

        entity.bucket().supersets( current_parts );

        for ( PartVector::iterator
              p = current_parts.begin() ; p != current_parts.end() ; ++p ) {
          if ( membership_is_induced( **p , entity.entity_rank() ) &&
               ! contain( induced_parts , **p ) ) {
            remove_parts.push_back( *p );
          }
        }

        internal_change_entity_parts( entity, induced_parts, remove_parts );
      }
    }
  }

  //------------------------------
  // The owners have complete knowledge of memberships.
  // Send membership information to sync the shared and ghosted copies.
  // Only need to do this for entities that have actually changed.

  {
    std::vector<EntityProc> send_list ;

    generate_send_list( m_sync_count, p_rank, m_entity_comm, send_list);
    
    CommAll comm( p_comm );

    pack_memberships( comm , send_list );

    comm.allocate_buffers( p_size / 4 );

    pack_memberships( comm , send_list );

    comm.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer( p );
      while ( buf.remaining() ) {

        PartVector owner_parts , current_parts , remove_parts ;

        EntityKey key ; buf.unpack<EntityKey>( key );
        unsigned count = 0 ; buf.unpack<unsigned>( count );
        for ( unsigned j = 0 ; j < count ; ++j ) {
          unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
          insert( owner_parts , * all_parts[ part_ord ] );
        }

        // Any current part that is not a member of owners_parts
        // must be removed.

        Entity * const entity = get_entity( key );

        entity->bucket().supersets( current_parts );

        for ( PartVector::iterator
              ip = current_parts.begin() ; ip != current_parts.end() ; ++ip ) {
          if ( OWNERSHIP_PART_COUNT <= (*ip)->mesh_meta_data_ordinal() &&
               ! contain( owner_parts , **ip ) ) {
            remove_parts.push_back( *ip );
          }
        }

        internal_change_entity_parts( *entity , owner_parts , remove_parts );
      }
    }
  }
}

} // namespace mesh
} // namespace stk

