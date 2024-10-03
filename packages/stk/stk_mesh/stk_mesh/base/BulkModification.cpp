// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/BulkModification.hpp>
#include <stddef.h>                     // for size_t
#include <set>                          // for _Rb_tree_const_iterator, etc
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer
#include <stk_util/parallel/CommSparse.hpp>
#include <utility>                      // for pair
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityKeyProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter




namespace stk {
namespace mesh {

typedef std::set<Entity , EntityLess> EntitySet;
typedef std::set<EntityKeyProc> EntityKeyProcSet;

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
      Entity const *irels_j = mesh.begin(entry, irank);
      Entity const *irels_e = mesh.end(entry, irank);
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
  for (Entity entity : entities)
  {
    if (mesh.bucket(entity).owned() || mesh.bucket(entity).shared())
    {
      construct_transitive_closure(mesh, closure, entity);
    }
  }
}

void construct_communication_set( const BulkData & bulk, const std::set<Entity,EntityLess> & closure, EntityKeyProcSet & communication_set)
{
  if (bulk.parallel_size() < 2) return;

  std::vector<int> commProcs;
  for ( Entity entity : closure ) {
    const bool owned = bulk.parallel_rank() == bulk.parallel_owner_rank(entity);
    const bool shared = bulk.bucket(entity).shared();

    // Add sharing processes and ghost-send processes to communication_set

    bulk.comm_procs(entity, commProcs);
    for(size_t j=0; j<commProcs.size(); j++)
    {
      if ( owned || shared ) {
        EntityKeyProc tmp( bulk.entity_key(entity) , commProcs[j] );
        communication_set.insert( tmp );
      }
    }
  }
}

}



void find_closure( const BulkData & bulk,
    const EntityVector& entities,
    EntityVector& entities_closure)
{

  entities_closure.clear();


  EntityKeyProcSet send_list;
  EntityLess entless(bulk);
  std::set<Entity,EntityLess>     temp_entities_closure(entless);

  STK_ThrowRequireMsg(bulk.in_synchronized_state(), "ERROR find_closure requires mesh in synchronized state");

  find_local_closure(bulk, temp_entities_closure, entities);

  construct_communication_set(bulk, temp_entities_closure, send_list);


  CommSparse commSparse( bulk.parallel() );

  //pack send_list for sizing
  for ( const EntityKeyProc& ep : send_list ) {
    commSparse.send_buffer( ep.second).pack<EntityKey>(ep.first);
  }


  commSparse.allocate_buffers();

  //pack send_list
  for ( const EntityKeyProc& ep : send_list ) {
    commSparse.send_buffer( ep.second).pack<EntityKey>(ep.first);
  }

  commSparse.communicate();

  //unpack the send_list into the temp entities closure set
  for ( int p = 0 ; p < bulk.parallel_size() ; ++p ) {
    CommBuffer & buf = commSparse.recv_buffer( p );
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

