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
#include <stk_mesh/base/Bucket.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

bool comm_mesh_verify_parallel_consistency(
  BulkData & M , std::ostream & error_log );

//----------------------------------------------------------------------

unsigned BulkData::determine_new_owner( Entity & entity ) const
{
  unsigned new_owner = entity.marked_for_destruction() ? ~0u : m_parallel_rank;

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

void BulkData::internal_update_distributed_index(
        std::vector<Entity*> & shared_new )
{
  std::vector< parallel::DistributedIndex::KeyType >
    new_entities_keys , del_entities_keys ;

  std::vector< parallel::DistributedIndex::KeyProc >
    new_entities_keyprocs ;

  for ( impl::EntityRepository::iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i ) {

    Entity & entity = * i->second ;

    if ( entity.marked_for_destruction() ) {
      // Has been destroyed
      del_entities_keys.push_back( entity.key().raw_key() );
    }
    else if ( EntityLogNoChange != entity.log_query() &&
              in_owned_closure( entity , m_parallel_rank ) ) {
      // Has been changed and is in owned closure, may be shared
      new_entities_keys.push_back( entity.key().raw_key() );
    }
  }

  //------------------------------
  // Update distributed index.

  m_entities_index.update_keys( new_entities_keys , del_entities_keys );

  //------------------------------
  // Inform creating processes:

  m_entities_index.query_to_usage( new_entities_keys , new_entities_keyprocs );

  {
    Entity * entity = NULL ;

    for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
          i =  new_entities_keyprocs.begin() ;
          i != new_entities_keyprocs.end() ; ++i ) {

      EntityKey key( & i->first );

      if ( (int) m_parallel_rank != i->second ) {
        // Another process also created or updated this entity.

        if ( entity == NULL || entity->key() != key ) {
          // Have not looked this entity up by key
          entity = get_entity( key );

          shared_new.push_back( entity );
        }

        // Add the other_process to the entity's sharing info.
        m_entity_repo.insert_comm_info( *entity, EntityCommInfo( 0 , i->second ) );
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_resolve_parallel_create_delete(
  const std::vector<Entity*> & del_entities )
{
  static const char method[] =
    "stk::mesh::BulkData::internal_resolve_parallel_create_delete" ;

  // 'del_entities' is guaranteed unique and sorted

  {
    const size_t ghosting_count = m_ghosting.size();

    std::vector<int> local_flags(  ghosting_count , 0 );
    std::vector<int> global_flags( ghosting_count , 0 );
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
        const bool     destroyed = entity->marked_for_destruction();
        const bool     remote_owner = entity->owner_rank() == proc ;
        const bool     shared = in_shared( *entity , proc );
        const bool     g_recv = in_receive_ghost( *entity );
        const bool     g_send = in_send_ghost( *entity , proc );

        if ( shared ) {

          // A shared entity is being deleted on the remote process.
          // Remove it from the sharing.

          m_entity_repo.erase_comm_info( *entity, EntityCommInfo( 0 , proc ) );

          PartVector add_part , remove_part ;

          if ( entity->sharing().empty() ) {
            // No longer shared, remove the shared part.
            remove_part.push_back( & m_mesh_meta_data.globally_shared_part() );
          }

          if ( remote_owner ) {
            // Remote owner is deleting a shared entity.
            // Have to determine new owner from remaining sharing processes.

            const unsigned new_owner = determine_new_owner( *entity );

            m_entity_repo.set_entity_owner_rank( *entity, new_owner);
            m_entity_repo.set_entity_sync_count( *entity, m_sync_count);

            if ( new_owner == m_parallel_rank ) {
              // Changing remotely owned to locally owned
              add_part.push_back( & m_mesh_meta_data.locally_owned_part() );
            }
          }

          if ( ! destroyed && ( ! add_part.empty() || ! remove_part.empty() ) ) {
            internal_change_entity_parts( *entity , add_part , remove_part );
          }
        }

        if ( g_send ) {
          // Remotely ghosted entity is being destroyed,
          // remove from ghosting list
          for ( size_t j = ghosting_count ; j-- ; ) {
            if ( m_entity_repo.erase_comm_info( *entity, EntityCommInfo( j , proc ) ) ) {
              local_flags[ j ] = 1 ;
            }
          }
        }
        else if ( remote_owner && g_recv ) {
          // Remotely owned entity is being destroyed.
          // This receive ghost must be destroyed.
          for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
            local_flags[ ec->ghost_id ] = 1 ;
          }
          m_entity_repo.comm_clear( *entity);
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

          m_entity_repo.set_entity_owner_rank( *entity, new_owner);
          m_entity_repo.set_entity_sync_count( *entity, m_sync_count);
        }

        for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
          local_flags[ ec->ghost_id ] = 1 ;
        }
        m_entity_repo.comm_clear( *entity);
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
    std::vector<Entity*> shared_modified ;

    // Update the parallel index and
    // output shared and modified entities.
    internal_update_distributed_index( shared_modified );

    // A shared_modified entity which was not created this
    // cycle and is locally owned send claim on ownership
    // to the sharing processes.

    {
      CommAll comm_all( m_parallel_machine );

      for ( std::vector<Entity*>::iterator
            i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
        Entity & entity = **i ;
        if ( entity.owner_rank() == m_parallel_rank &&
             entity.log_query()  != EntityLogCreated ) {

          for ( PairIterEntityComm
                jc = entity.sharing() ; ! jc.empty() ; ++jc ) {
            comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity.key() );
          }
        }
      }

      comm_all.allocate_buffers( m_parallel_size / 4 );

      for ( std::vector<Entity*>::iterator
            i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
        Entity & entity = **i ;
        if ( entity.owner_rank() == m_parallel_rank &&
             entity.log_query()  != EntityLogCreated ) {

          for ( PairIterEntityComm
                jc = entity.sharing() ; ! jc.empty() ; ++jc ) {
            comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity.key() );
          }
        }
      }

      comm_all.communicate();

      for ( unsigned p = 0 ; p < m_parallel_size ; ++p ) {
        CommBuffer & buf = comm_all.recv_buffer( p );
        EntityKey key ;
        while ( buf.remaining() ) {
          buf.unpack<EntityKey>( key );

          Entity & entity = * get_entity( key );

          // Set owner, will correct part membership later
          m_entity_repo.set_entity_owner_rank( entity, p);
        }
      }
    }

    // Update shared created entities.
    // - Revise ownership to selected processor
    // - Update sharing.
    // - Work backward so the 'in_owned_closure' function
    //   can evaluate related higher ranking entities.

    std::ostringstream error_msg ;
    int error_flag = 0 ;

    PartVector shared_part , owned_part ;
    shared_part.push_back( & m_mesh_meta_data.globally_shared_part() );
    owned_part.push_back(  & m_mesh_meta_data.locally_owned_part() );

    for ( std::vector<Entity*>::const_iterator
          i = shared_modified.end() ; i != shared_modified.begin() ; ) {

      Entity * entity = *--i ;

      if ( entity->owner_rank() == m_parallel_rank &&
           entity->log_query() == EntityLogCreated ) {

        // Created and not claimed by an existing owner

        const unsigned new_owner = determine_new_owner( *entity );

        m_entity_repo.set_entity_owner_rank( *entity, new_owner);
      }

      if ( entity->owner_rank() != m_parallel_rank ) {
        // Do not own it and still have it.
        // Remove the locally owned, add the globally_shared
        m_entity_repo.set_entity_sync_count( *entity, m_sync_count);
        internal_change_entity_parts( *entity , shared_part , owned_part );
      }
      else if ( ! entity->sharing().empty() ) {
        // Own it and has sharing information.
        // Add the globally_shared
        internal_change_entity_parts( *entity , shared_part , PartVector() );
      }
      else {
        // Own it and does not have sharing information.
        // Remove the globally_shared
        internal_change_entity_parts( *entity , PartVector() , shared_part );
      }

      // Newly created shared entity had better be in the owned closure
      if ( ! in_owned_closure( *entity , m_parallel_rank ) ) {
        if ( 0 == error_flag ) {
          error_flag = 1 ;
          error_msg
            << std::endl
            << "P" << m_parallel_rank << ": " << method << " FAILED"
            << std::endl
            << "  The following entities were declared on multiple processors,"
            << std::endl
            << "  cannot be parallel-shared, and were declared with"
            << "  parallel-ghosting information. { "
            << std::endl ;
        }
        error_msg << "    " ;
        print_entity_key( error_msg , m_mesh_meta_data , entity->key() );
        error_msg << " also declared on" ;
        for ( PairIterEntityComm ec = entity->sharing(); ! ec.empty() ; ++ec ) {
          error_msg << " P" << ec->proc ;
        }
        error_msg << std::endl ;
      }
    }
    if ( error_flag ) { error_msg << "}" << std::endl ; }

    all_reduce( m_parallel_machine , ReduceMax<1>( & error_flag ) );

    if ( error_flag ) {
      throw std::runtime_error( error_msg.str() );
    }

    const size_t n_old = m_entity_comm.size();

    m_entity_comm.insert( m_entity_comm.end() ,
                          shared_modified.begin() , shared_modified.end() );

    std::inplace_merge( m_entity_comm.begin() ,
                        m_entity_comm.begin() + n_old ,
                        m_entity_comm.end() ,
                        EntityLess() );

    {
      std::vector<Entity*>::iterator i =
        std::unique( m_entity_comm.begin() , m_entity_comm.end() );

      m_entity_comm.erase( i , m_entity_comm.end() );
   }
  }
}

//----------------------------------------------------------------------

namespace {

// Enforce that shared entities must be in the owned closure:

void destroy_dependent_ghosts( BulkData & mesh , Entity * entity )
{
  static const char method[] =
    "stk::mesh::BulkData::modification_end() { destroy_dependent_ghosts() }" ;

  for ( ; ; ) {
    PairIterRelation rel = entity->relations();

    if ( rel.empty() ) { break ; }

    Entity * e = rel.back().entity();

    if ( e->entity_rank() < entity->entity_rank() ) { break ; }

    if ( in_owned_closure( *e , mesh.parallel_rank() ) ) {
      throw std::logic_error( std::string(method) );
    }

    destroy_dependent_ghosts( mesh , e );
  }

  mesh.destroy_entity( entity );
}

void resolve_modified_shared( BulkData & mesh )
{
  for ( std::vector<Entity*>::const_iterator
        i =  mesh.entity_comm().end() ;
        i != mesh.entity_comm().begin() ; ) {

    Entity * entity = *--i ;

    if ( ! entity->sharing().empty() &&
         ! in_owned_closure( *entity , mesh.parallel_rank() ) ) {
      // An entity with sharing information is not in
      // the owned closure; therefore, the entity cannot
      // be shared and must be destroyed.  If it should be
      // ghosted it will be subsequently re-ghosted.

      destroy_dependent_ghosts( mesh , entity );
    }
  }
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

  // If a shared entity is modified such that
  // it is no longer in the owned closure
  // then destroy that entity and all ghost
  // entities that depend upon it.

  resolve_modified_shared( *this );

  int local_change_count[2] = { 0 , 0 };
  int global_change_count[2] = { 0 , 0 };

  std::vector<Entity*> del_entities ;

  for ( impl::EntityRepository::iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i ) {
   Entity * entity = i->second ;
    if ( EntityLogCreated == entity->log_query() ) {
      ++ local_change_count[0] ; // Created
    }
    else if ( entity->marked_for_destruction() ) {
      del_entities.push_back( entity );
      ++ local_change_count[1] ; // Deleted
    }
  }

  all_reduce_sum( m_parallel_machine ,
                  local_change_count , global_change_count , 2 );

  if ( global_change_count[0] || global_change_count[1] ) {

    internal_resolve_parallel_create_delete( del_entities );
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

  m_bucket_repository.internal_sort_bucket_entities();

  // ------------------------------

  m_sync_state = SYNCHRONIZED ;

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { PART_ORD_UNIVERSAL = 0 };
enum { PART_ORD_OWNED     = 1 };
enum { PART_ORD_SHARED    = 2 };

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

void pack_part_memberships( CommAll & comm ,
                            const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity & entity = * i->first ;

    std::pair<const unsigned *, const unsigned *>
      part_ord = entity.bucket().superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == PART_ORD_SHARED ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( entity.key() )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
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

  const Part & part_universal = meta.universal_part();
  const Part & part_owned  = meta.locally_owned_part();
  const Part & part_shared = meta.globally_shared_part();

  // Quick verification of part ordinal assumptions

  if ( PART_ORD_UNIVERSAL != part_universal.mesh_meta_data_ordinal() ||
       PART_ORD_OWNED     != part_owned.mesh_meta_data_ordinal() ||
       PART_ORD_SHARED    != part_shared.mesh_meta_data_ordinal() ) {
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

    pack_part_memberships( comm , send_list );

    comm.allocate_buffers( p_size / 4 );

    pack_part_memberships( comm , send_list );

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
          Part * const part = *ip ;
          const unsigned part_ord = part->mesh_meta_data_ordinal();
          if ( PART_ORD_UNIVERSAL != part_ord &&
               PART_ORD_OWNED     != part_ord &&
               PART_ORD_SHARED    != part_ord &&
               ! contain( owner_parts , *part ) ) {
            remove_parts.push_back( part );
          }
        }

        internal_change_entity_parts( *entity , owner_parts , remove_parts );
      }
    }
  }
}

} // namespace mesh
} // namespace stk

