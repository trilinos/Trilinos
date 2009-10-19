#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Comm.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

enum { OWNERSHIP_PART_COUNT = 3 };

namespace {

void pack_induced_memberships( CommAll & comm ,
                               const std::vector<EntityProc> & shared )
{
  for ( std::vector<EntityProc>::const_iterator
        i = shared.begin() ; i != shared.end() ; ++i ) {

    Entity & entity = * i->first ;

    if ( entity.owner_rank() == i->second ) { // Only send to owner.

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
                         const std::vector<EntityProc> & shared ,
                         const std::vector<Ghosting*>  & ghosting ,
                               std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = shared.begin() ; i != shared.end() ; ++i ) {

    Entity & entity = * i->first ;

    if ( entity.owner_rank() == p_rank &&
         entity.synchronized_count() == sync_count ) {
      send_list.push_back( *i );
    }
  }

  for ( std::vector<Ghosting*>::const_iterator
        i = ghosting.begin() ; i != ghosting.end() ; ++i ) {

    const std::vector<EntityProc> & gh_list = (*i)->send();

    for ( std::vector<EntityProc>::const_iterator
          ig = gh_list.begin() ; ig != gh_list.end() ; ++ig ) {

      Entity & entity = * ig->first ;

      if ( entity.synchronized_count() == sync_count ) {
        send_list.push_back( *ig );
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

    pack_induced_memberships( comm , m_shares_all );

    comm.allocate_buffers( p_size / 4 );

    pack_induced_memberships( comm , m_shares_all );

    comm.communicate();

    for ( std::vector<EntityProc>::iterator
          i = m_shares_all.begin() ; i != m_shares_all.end() ; ) {

      Entity * const entity = i->first ;

      if ( entity->owner_rank() != p_rank ) {
        // Skip not-owned
        for ( ++i ; i != m_shares_all.end() && i->first == entity ; ++i );
      }
      else {
        PartVector empty , induced_parts , current_parts , remove_parts ;

        induced_part_membership( *entity , empty , induced_parts );

        for ( ; i != m_shares_all.end() && i->first == entity ; ++i ) {

          CommBuffer & buf = comm.recv_buffer( i->second );

          unsigned count = 0 ; buf.unpack<unsigned>( count );
          for ( unsigned j = 0 ; j < count ; ++j ) {
            unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
            insert( induced_parts , * all_parts[ part_ord ] );
          }
        }

        // Remove any part that is an induced part but is not
        // in the induced parts list.

        entity->bucket().supersets( current_parts );

        for ( PartVector::iterator
              p = current_parts.begin() ; p != current_parts.end() ; ++p ) {
          if ( membership_is_induced( **p , entity->entity_type() ) &&
               ! contain( induced_parts , **p ) ) {
            remove_parts.push_back( *p );
          }
        }

        internal_change_entity_parts( *entity, induced_parts, remove_parts );
      }
    }
  }

  //------------------------------
  // The owners have complete knowledge of memberships.
  // Send membership information to sync the shared and ghosted copies.
  // Only need to do this for entities that have actually changed.

  {
    std::vector<EntityProc> send_list ;

    generate_send_list( m_sync_count, p_rank, m_shares_all, m_ghosting, send_list);
    
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

        Entity * const entity = get_entity( entity_type(key), entity_id(key) );

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

