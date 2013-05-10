/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <set>
#include <stdexcept>
#include <sstream>

#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

namespace stk {
namespace mesh {

typedef std::set<Entity , EntityLess> EntitySet;
typedef std::set<EntityKeyProc> EntityProcSet;

namespace {

void construct_transitive_closure(const BulkData& mesh, std::set<Entity,EntityLess> & closure , Entity entry )
{

  std::pair< std::set<Entity,EntityLess>::const_iterator , bool >
    result = closure.insert( entry );

  // A new insertion, must also insert the closure
  if ( result.second ) {

    const EntityRank erank = mesh.entity_rank(entry);

    // insert entities with relations of lower rank into the closure
    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
    {
      Entity const *irels_j = mesh.begin_entities(entry, irank);
      Entity const *irels_e = mesh.end_entities(entry, irank);
      for (; irels_j != irels_e; ++irels_j)
      {
        if (mesh.is_valid(*irels_j)) {
          construct_transitive_closure(mesh, closure, *irels_j);
        }
      }
    }
  }
}

void find_local_closure (const BulkData& mesh, std::set<Entity,EntityLess> & closure, const EntityVector & entities)
{
  for (EntityVector::const_iterator i = entities.begin();
      i != entities.end(); ++i)
  {
    construct_transitive_closure(mesh, closure, *i);
  }
}

void construct_communication_set( const BulkData & bulk, const std::set<Entity,EntityLess> & closure, EntityProcSet & communication_set)
{
  if (bulk.parallel_size() < 2) return;

  for ( std::set<Entity,EntityLess>::const_iterator
        i = closure.begin(); i != closure.end(); ++i) {

    Entity entity = *i;

    const bool owned = bulk.parallel_rank() == bulk.parallel_owner_rank(entity);

    // Add sharing processes and ghost-send processes to communication_set

    for ( PairIterEntityComm ec = bulk.entity_comm(bulk.entity_key(entity)); ! ec.empty() ; ++ec ) {
      if ( owned || ec->ghost_id == 0 ) {
        EntityKeyProc tmp( bulk.entity_key(entity) , ec->proc );
        communication_set.insert( tmp );
      }
    }
  }
}

size_t count_non_used_entities( const BulkData & bulk, const EntityVector & entities)
{
  const unsigned proc_local = bulk.parallel_rank();
  size_t non_used_entities = 0;

  for ( EntityVector::const_iterator
        i = entities.begin(); i != entities.end(); ++i ) {
    if ( ! in_owned_closure( bulk, *i , proc_local ) ) {
      ++non_used_entities;
    }
  }

  return non_used_entities;
}

}



void find_closure( const BulkData & bulk,
    const std::vector< Entity> & entities,
    std::vector< Entity> & entities_closure)
{

  entities_closure.clear();


  EntityProcSet send_list;
  EntityLess entless(bulk);
  std::set<Entity,EntityLess>     temp_entities_closure(entless);

  const bool bulk_not_synchronized = bulk.synchronized_state() != BulkData::SYNCHRONIZED;
  const size_t non_used_entities = bulk_not_synchronized ? 0 : count_non_used_entities(bulk, entities);

  const bool local_bad_input = bulk_not_synchronized || (0 < non_used_entities);

  //Can skip if error on input
  if ( !local_bad_input) {

    find_local_closure(bulk, temp_entities_closure, entities);

    construct_communication_set(bulk, temp_entities_closure, send_list);
  }


  CommAll all( bulk.parallel() );

  //pack send_list for sizing
  for ( EntityProcSet::const_iterator
      ep = send_list.begin() ; ep != send_list.end() ; ++ep ) {
    all.send_buffer( ep->second).pack<EntityKey>(ep->first);
  }


  const bool global_bad_input = all.allocate_buffers( bulk.parallel_size() / 4 , false, local_bad_input );

  if (global_bad_input) {

    std::ostringstream msg;
    //parallel consisent throw
    if (bulk_not_synchronized) {
      msg << "stk::mesh::find_closure( const BulkData & bulk, ... ) bulk is not synchronized";
    }
    else if ( 0 < non_used_entities) {
      msg << "stk::mesh::find_closure( const BulkData & bulk, std::vector<Entity> entities, ... ) \n"
          << "entities contains " << non_used_entities << " non locally used entities \n";
    }

    throw std::runtime_error(msg.str());
  }


  //pack send_list
  for ( EntityProcSet::const_iterator
      ep = send_list.begin() ; ep != send_list.end() ; ++ep ) {
    all.send_buffer( ep->second).pack<EntityKey>(ep->first);
  }


  all.communicate();

  //unpack the send_list into the temp entities closure set
  for ( int p = 0 ; p < bulk.parallel_size() ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    EntityKey k ;
    while ( buf.remaining() ) {
      buf.unpack<EntityKey>( k );
      Entity e = bulk.get_entity(k);
      temp_entities_closure.insert(e);
    }
  }

  //copy the set into the entities_closure vector
  entities_closure.assign(temp_entities_closure.begin(), temp_entities_closure.end());
}

} // namespace mesh
} // namespace stk

