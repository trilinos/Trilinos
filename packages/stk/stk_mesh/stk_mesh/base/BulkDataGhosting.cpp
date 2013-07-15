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

#include <cstring>
#include <set>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Trace.hpp>

using std::strcmp;

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

Ghosting & BulkData::create_ghosting( const std::string & name )
{
  require_ok_to_modify();

  // Verify name is the same on all processors,
  // if not then throw an exception on all processors.
  if (parallel_size() > 1) {
    CommBroadcast bc( parallel() , 0 );

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().skip<char>( name.size() + 1 );
    }

    bc.allocate_buffer();

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack<char>( name.c_str() , name.size() + 1 );
    }

    bc.communicate();

    const char * const bc_name =
      reinterpret_cast<const char *>( bc.recv_buffer().buffer() );

    int error = 0 != strcmp( bc_name , name.c_str() );

    all_reduce( parallel() , ReduceMax<1>( & error ) );

    ThrowErrorMsgIf( error, "Parallel name inconsistency");
  }

  Ghosting * const g =
    new Ghosting( *this , name , m_ghosting.size() , m_sync_count );

  m_ghosting.push_back( g );

  return *g ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void insert_transitive_closure( BulkData& bulk_data,
                                std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry );

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & new_recv ,
        std::set< EntityProc , EntityLess > & new_send );

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< EntityKey > & new_recv );

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::destroy_ghosting( Ghosting& ghost_layer )
{
  std::vector<EntityKey> receive_list;
  ghost_layer.receive_list(receive_list);
  change_ghosting(ghost_layer, std::vector<stk::mesh::EntityProc>(), receive_list);
}

//----------------------------------------------------------------------

void BulkData::destroy_all_ghosting()
{
  Trace_("stk::mesh::BulkData::destroy_all_ghosting");

  require_ok_to_modify();

  // Clear Ghosting data

  for ( std::vector<Ghosting*>::iterator
        ig = m_ghosting.begin() ; ig != m_ghosting.end() ; ++ig ) {
    Ghosting & gh = **ig ;
    gh.m_sync_count = m_sync_count ;
  }

  // Iterate backwards so as not to invalidate a closure.

  for ( EntityCommListInfoVector::reverse_iterator
        i =  m_entity_comm_list.rbegin() ;
        i != m_entity_comm_list.rend() ; ++i) {

    if ( in_receive_ghost( i->key ) ) {
      entity_comm_clear( i->key );
      destroy_entity( i->entity );
      i->key = EntityKey();
    }
    else {
      entity_comm_clear_ghosting(i->key);
      if ( entity_comm(i->key).empty() ) {
        i->key = EntityKey();
      }
    }
  }

  EntityCommListInfoVector::iterator i =
    std::remove_if( m_entity_comm_list.begin() ,
                    m_entity_comm_list.end() , IsInvalid() );

  m_entity_comm_list.erase( i , m_entity_comm_list.end() );
}

//----------------------------------------------------------------------

void BulkData::change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
  Trace_("stk::mesh::BulkData::change_ghosting");

  //----------------------------------------
  // Verify inputs:

  require_ok_to_modify();

  const bool ok_mesh  = & BulkData::get(ghosts) == this ;
  const bool ok_ghost = 1 < ghosts.ordinal();
  bool ok_add    = true ;
  bool ok_remove = true ;

  // Verify all 'add' are locally owned.

  for ( std::vector<EntityProc>::const_iterator
        i = add_send.begin() ; ok_add && i != add_send.end() ; ++i ) {
    ok_add = parallel_owner_rank(i->first) == parallel_rank();
  }

  // Verify all 'remove' are members of the ghosting.

  for ( std::vector<EntityKey>::const_iterator
        i = remove_receive.begin() ;
        ok_remove && i != remove_receive.end() ; ++i ) {
    ok_remove = in_receive_ghost( ghosts , *i );
  }

  int ok = ok_mesh && ok_ghost && ok_add && ok_remove ;

  all_reduce( parallel() , ReduceMin<1>( & ok ) );

  if ( 0 == ok ) {
    std::ostringstream msg ;
    msg << "For ghosts " << ghosts.name() << ", " ;
    if ( ! ok_mesh )  { msg << " : Mesh does not own this ghosting" ; }
    if ( ! ok_ghost ) { msg << " : Cannot modify this ghosting" ; }
    if ( ! ok_add ) {
      msg << " : Not owned add {" ;
      for ( std::vector<EntityProc>::const_iterator
            i = add_send.begin() ; i != add_send.end() ; ++i ) {
        if ( parallel_owner_rank(i->first) != parallel_rank() ) {
          msg << " " << identifier(i->first);
        }
      }
      msg << " }" ;
    }
    if ( ! ok_remove ) {
      msg << " : Not in ghost receive {" ;
      for ( std::vector<EntityKey>::const_iterator
            i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
        if ( ! in_receive_ghost( ghosts , *i ) ) {
          msg << " " << i->id();
        }
      }
    }

    ThrowErrorMsg( msg.str() );
  }
  //----------------------------------------
  // Change the ghosting:

  internal_change_ghosting( ghosts , add_send , remove_receive );
}

//----------------------------------------------------------------------

void BulkData::internal_change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive,
  bool is_full_regen)
{
  Trace_("stk::mesh::BulkData::internal_change_ghosting");

  const MetaData & meta = m_mesh_meta_data ;
  const unsigned rank_count = meta.entity_rank_count();
  const int p_size = m_parallel_size ;

  //------------------------------------
  // Copy ghosting lists into more efficiently editted container.
  // The send and receive lists must be in entity rank-order.

  std::set<EntityProc , EntityLess> new_send(EntityLess(*this));
  std::set<EntityKey>               new_recv;

  //------------------------------------
  // Insert the current ghost receives and then remove from that list.

  // This if-check is an optimization; if doing a full regen
  // then we are removing all ghosting information and new_recv should
  // be left empty.
  if ( !is_full_regen ) {

    // Iterate over all entities with communication information, adding
    // the entity if it's a ghost on this process. new_recv will contain
    // all ghosts on this process by the end of the loop.
    for ( EntityCommListInfoVector::const_iterator
          i = comm_list().begin() ; i != comm_list().end() ; ++i ) {
      if ( in_receive_ghost( ghosts , i->key ) ) {
        new_recv.insert( i->key );
      }
    }

    // Remove any entities that are in the remove list.

    for ( std::vector<EntityKey>::const_iterator
          i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
      new_recv.erase( *i );
    }

    // Keep the closure of the remaining received ghosts.
    // Working from highest-to-lowest key (rank entity type)
    // results in insertion of the closure because
    // inserted entities will get looped over after they are inserted.

    // Insertion will not invalidate the associative container's iterator.

    for ( std::set<EntityKey>::reverse_iterator
          i = new_recv.rbegin() ; i != new_recv.rend() ; ++i) {
      const unsigned erank = i->rank();

      Entity e = get_entity(*i); // Could be performance issue? Not if you're just doing full regens

      for (EntityRank irank = stk::topology::BEGIN_RANK;
          irank < erank;
          ++irank)
      {
        Entity const *rels_i = begin(e, irank);
        Entity const *rels_e = end(e, irank);
        for (; rels_i != rels_e; ++rels_i)
        {
          if ( is_valid(*rels_i) && in_receive_ghost( ghosts , entity_key(*rels_i) ) )
          {
            new_recv.insert( entity_key(*rels_i) );
          }
        }
      }
    }
  }

  //  Initialize the new_send from the new_recv
  comm_recv_to_send( *this , new_recv , new_send );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting

  for ( std::vector< EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
    insert_transitive_closure( *this, new_send , *i );
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  comm_sync_send_recv( *this , new_send , new_recv );

  // The new_send list is now parallel complete and parallel accurate
  // The new_recv has those ghost entities that are to be kept.
  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  bool removed = false ;

  for ( EntityCommListInfoVector::reverse_iterator
        i = m_entity_comm_list.rbegin() ; i != m_entity_comm_list.rend() ; ++i) {

    const bool is_owner = i->owner == m_parallel_rank ;
    const bool remove_recv = ( ! is_owner ) &&
                             0 == new_recv.count(i->key);

    if ( is_owner ) {
      // Is owner, potentially removing ghost-sends
      // Have to make a copy

      std::vector<EntityCommInfo> comm_ghost ;
      const PairIterEntityComm ec = entity_comm(i->key, ghosts);
      comm_ghost.assign( ec.first , ec.second );

      for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
        const EntityCommInfo tmp = comm_ghost.back();

        if ( 0 == new_send.count( EntityProc( i->entity , tmp.proc ) ) ) {
          entity_comm_erase(i->key, tmp);
        }
      }
    }
    else if ( remove_recv ) {
      entity_comm_erase(i->key, ghosts);
    }

    if ( entity_comm(i->key).empty() ) {
      removed = true ;
      i->key = EntityKey(); // No longer communicated
      if ( remove_recv ) {
        ThrowRequireMsg( destroy_entity( i->entity ),
                         " FAILED attempt to destroy entity: " << identifier(i->entity) );
      }
    }
  }

  if ( removed ) {
    EntityCommListInfoVector::iterator i =
      std::remove_if( m_entity_comm_list.begin() ,
                      m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
  }

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the comm list.
  // Unpacking must proceed in entity-rank order so that higher ranking
  // entities that have relations to lower ranking entities will have
  // the lower ranking entities unpacked first.  The higher and lower
  // ranking entities may be owned by different processes,
  // as such unpacking must be performed in rank order.

  {
    const size_t entity_comm_size = m_entity_comm_list.size();

    CommAll comm( m_parallel_machine );

    for ( std::set< EntityProc , EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {

      Entity entity = j->first;
      const int proc = j->second;

      if ( ! in_ghost( ghosts , entity_key(entity) , proc ) ) {
        // Not already being sent , must send it.
        CommBuffer & buf = comm.send_buffer( proc );
        buf.pack<unsigned>( entity_rank(entity) );
        pack_entity_info(*this, buf , entity );
        pack_field_values(*this, buf , entity );
      }
    }

    comm.allocate_buffers( p_size / 4 );

    for ( std::set< EntityProc , EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {

      Entity entity = j->first;
      const int proc = j->second;

      if ( ! in_ghost( ghosts , entity_key(entity) , proc ) ) {
        // Not already being sent , must send it.
        CommBuffer & buf = comm.send_buffer( proc );
        buf.pack<unsigned>( entity_rank(entity) );
        pack_entity_info(*this,  buf , entity );
        pack_field_values(*this, buf , entity );

        entity_comm_insert(entity, EntityCommInfo(ghosts.ordinal(), proc));

        EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity)};
        m_entity_comm_list.push_back( comm_info );
      }
    }

    comm.communicate();

    std::ostringstream error_msg ;
    int error_count = 0 ;

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {
      for ( int p = 0 ; p < p_size ; ++p ) {
        CommBuffer & buf = comm.recv_buffer(p);
        while ( buf.remaining() ) {
          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            unsigned this_rank = ~0u ;
            buf.peek<unsigned>( this_rank );

            if ( this_rank != rank ) break ;

            buf.unpack<unsigned>( this_rank );
          }

          PartVector parts ;
          std::vector<Relation> relations ;
          EntityKey key ;
          int owner = ~0u ;

          unpack_entity_info( buf, *this, key, owner, parts, relations );

          // Must not have the locally_owned_part or globally_shared_part

          remove( parts , meta.locally_owned_part() );
          remove( parts , meta.globally_shared_part() );

          std::pair<Entity ,bool> result =
            m_entity_repo.internal_create_entity( key );

          Entity entity = result.first;
          const bool created   = result.second ;

          require_entity_owner( entity , owner );

          internal_change_entity_parts( entity , parts , PartVector() );

          if ( created ) {
            log_created_parallel_copy( entity );
            this->set_parallel_owner_rank( entity, owner);
          }

          declare_relation( entity , relations );

          if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }

          const EntityCommInfo tmp( ghosts.ordinal() , owner );

          if ( entity_comm_insert(entity, tmp) ) {
            EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity)};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }
    }

    if (parallel_size() > 1) {
      all_reduce( m_parallel_machine , ReduceSum<1>( & error_count ) );
    }

    ThrowErrorMsgIf( error_count, error_msg.str() );

    if ( entity_comm_size < m_entity_comm_list.size() ) {
      // Added new ghosting entities to the list,
      // must now sort and merge.

      EntityCommListInfoVector::iterator i = m_entity_comm_list.begin();
      i += entity_comm_size ;
      std::sort( i , m_entity_comm_list.end() );
      std::inplace_merge( m_entity_comm_list.begin() , i ,
                          m_entity_comm_list.end() );
      m_entity_comm_list.erase( std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() ) ,
                                m_entity_comm_list.end() );

      internal_sync_comm_list_owners();
    }
  }

  ghosts.m_sync_count = m_sync_count ;
}

//----------------------------------------------------------------------

namespace {

void insert_transitive_closure( BulkData& bulk_data,
                                std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry )
{
  // Do not insert if I can determine that this entity is already
  // owned or shared by the receiving processor.

  if ( entry.second != bulk_data.parallel_owner_rank(entry.first) &&
       ! bulk_data.in_shared( bulk_data.entity_key(entry.first) , entry.second ) ) {

    std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
      result = new_send.insert( entry );

    if ( result.second ) {
      // A new insertion, must also insert the closure

      const unsigned erank = bulk_data.entity_rank(entry.first);

      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
      {
        Entity const *irels_j = bulk_data.begin(entry.first, irank);;
        Entity const *irels_e = bulk_data.end(entry.first, irank);
        for (; irels_j != irels_e; ++irels_j)
        {
          if (bulk_data.is_valid(*irels_j)) {
            EntityProc tmp( *irels_j , entry.second );
            insert_transitive_closure( bulk_data, new_send , tmp );
          }
        }
      }
    }
  }
}

// Fill a new send list from the receive list.

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & new_recv,
        std::set< EntityProc , EntityLess > & new_send )
{
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::set<EntityKey>::const_iterator
            i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
      Entity e = mesh.get_entity(*i); // Performance issue? Not if you're just doing full regens
      const int owner = mesh.parallel_owner_rank(e);
      const EntityKey key = mesh.entity_key(e);
      all.send_buffer( owner ).pack<EntityKey>( key );
    }
    if (phase == 0) { //allocation phase
      all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );
    }
    else { //communication phase
      all.communicate();
    }
  }

  for ( int proc_rank = 0 ; proc_rank < parallel_size ; ++proc_rank ) {
    CommBuffer & buf = all.recv_buffer(proc_rank);
    while ( buf.remaining() ) {
      EntityKey key ;
      buf.unpack<EntityKey>( key );
      EntityProc tmp( mesh.get_entity( key ) , proc_rank );
      new_send.insert( tmp );
    }
  }
}

// Synchronize the send list to the receive list.

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< EntityKey > & new_recv )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< EntityProc , EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(i->first);
    all.send_buffer( i->second ).skip<EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<EntityKey>(1).skip<int>(1);
    }
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Communication packing (with message content comments):
  for ( std::set< EntityProc , EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ) {
    const int owner = mesh.parallel_owner_rank(i->first);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey entity_key = mesh.entity_key(i->first);
    const int proc = i->second;

    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< EntityProc , EntityLess >::iterator jrem = i ; ++i ;
      new_send.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  all.communicate();

  // Communication unpacking:
  for ( int p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      Entity const e = mesh.get_entity( entity_key );

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        ThrowRequireMsg(mesh.is_valid(e),
            "Unknown entity key: " <<
            MetaData::get(mesh).entity_rank_name(entity_key.rank()) <<
            "[" << entity_key.id() << "]");
        EntityProc tmp( e , proc );
        new_send.insert( tmp );
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        new_recv.insert( mesh.entity_key(e) );
      }
    }
  }
}

void insert_upward_relations(const BulkData& bulk_data, Entity rel_entity,
                             const EntityRank rank_of_orig_entity,
                             const int my_rank,
                             const int share_proc,
                             std::vector<EntityProc>& send)
{
  EntityRank rel_entity_rank = bulk_data.entity_rank(rel_entity);

  // If related entity is higher rank, I own it, and it is not
  // already shared by proc, ghost it to the sharing processor.
  if ( rank_of_orig_entity < rel_entity_rank &&
       bulk_data.parallel_owner_rank(rel_entity) == my_rank &&
       ! bulk_data.in_shared( bulk_data.entity_key(rel_entity) , share_proc ) ) {

    EntityProc entry( rel_entity , share_proc );
    send.push_back( entry );

    // There may be even higher-ranking entities that need to be ghosted, so we must recurse
    const EntityRank end_rank = bulk_data.mesh_meta_data().entity_rank_count();
    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
    {
      Entity const *irels_j = bulk_data.begin(rel_entity, irank);
      Entity const *irels_e = bulk_data.end(rel_entity, irank);
      for ( ; irels_j != irels_e; ++irels_j)
      {
        Entity const rel_of_rel_entity = *irels_j;
        if (bulk_data.is_valid(rel_of_rel_entity)) {
          insert_upward_relations(bulk_data, rel_of_rel_entity, rel_entity_rank, my_rank, share_proc, send);
        }
      }
    }
  }
}

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_regenerate_shared_aura()
{
  Trace_("stk::mesh::BulkData::internal_regenerate_shared_aura");

  require_ok_to_modify();

  std::vector<EntityProc> send ;

  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();

  // Iterate over all entities with communication info, get the sharing
  // comm info for each entity, and ensure that upwardly related
  // entities to the shared entity are ghosted on the sharing proc.
  for ( EntityCommListInfoVector::const_iterator
      i = comm_list().begin() ; i != comm_list().end() ; ++i ) {

    const unsigned erank = i->key.rank();

    const PairIterEntityComm sharing = entity_comm_sharing(i->key);

    for ( size_t j = 0 ; j < sharing.size() ; ++j ) {

      const int share_proc = sharing[j].proc ;

      for (EntityRank k_rank = stk::topology::BEGIN_RANK; k_rank < end_rank; ++k_rank)
      {
        Entity const *rels_itr = begin(i->entity, k_rank);
        Entity const *rels_end = end(i->entity, k_rank);
        for ( ; rels_itr != rels_end; ++rels_itr)
        {
          if (is_valid(*rels_itr)) {
            insert_upward_relations(*this, *rels_itr, erank, m_parallel_rank, share_proc, send);
          }
        }
      }
    }
  }

  // Add new aura, remove all of the old aura.
  // The change_ghosting figures out what to actually delete and add.
  internal_change_ghosting( shared_aura() , send , std::vector<EntityKey>(), true /*full regen*/ );
}

} // namespace mesh
} // namespace stk
