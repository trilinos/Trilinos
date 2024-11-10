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

#include <stk_mesh/baseImpl/AuraGhosting.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>

namespace stk {
namespace mesh {
namespace impl {

AuraGhosting::AuraGhosting()
: m_sendAura(),
  m_scratchSpace()
{
}

AuraGhosting::~AuraGhosting()
{
}

void AuraGhosting::generate_aura(BulkData& bulkData)
{
  m_sendAura.reset(bulkData.get_size_of_entity_index_space());
  fill_send_aura_entities(bulkData, m_sendAura);

  change_ghosting(bulkData, m_sendAura);
}

void AuraGhosting::remove_aura(BulkData& bulkData)
{
  EntityVector recvAuraEntitiesToRemove;
  bulkData.aura_ghosting().receive_list(recvAuraEntitiesToRemove);
  bulkData.internal_change_ghosting(bulkData.aura_ghosting(), {}, recvAuraEntitiesToRemove);
}

void AuraGhosting::fill_send_aura_entities(BulkData& bulkData,
                                           EntityProcMapping& sendAuraEntityProcs)
{
  const EntityRank endRank = static_cast<EntityRank>(bulkData.mesh_meta_data().entity_rank_count());
  const EntityRank maxRank = static_cast<EntityRank>(endRank-1);

  // Iterate over all shared entities, ensure that upwardly related
  // entities to each shared entity will be ghosted to the sharing proc.
  Selector shared = bulkData.mesh_meta_data().globally_shared_part();

  std::vector<int> sharingProcs;
  impl::for_each_selected_entity_run_no_threads(bulkData, stk::topology::NODE_RANK, shared,
    [&sendAuraEntityProcs, &sharingProcs, &endRank, &maxRank]
    (const BulkData& bulk, const MeshIndex& meshIndex) {
      const Bucket& bucket = *meshIndex.bucket;
      const unsigned bucketOrd = meshIndex.bucket_ordinal;

      bulk.comm_shared_procs(bucket[bucketOrd], sharingProcs);

      static constexpr EntityRank nextHigherRank = stk::topology::EDGE_RANK;
      for (EntityRank higherRank = nextHigherRank; higherRank < endRank; ++higherRank) {
        const unsigned num_rels = bucket.num_connectivity(bucketOrd, higherRank);
        if (num_rels > 0) {
          const Entity* rels     = bucket.begin(bucketOrd, higherRank);

          for (unsigned r = 0; r < num_rels; ++r) {
            if (bulk.parallel_rank() == bulk.parallel_owner_rank(rels[r])) {
              stk::mesh::impl::insert_upward_relations_for_owned(bulk, rels[r], higherRank, maxRank, sharingProcs, sendAuraEntityProcs);
            }
          }
        }
      }
    }
  ); // for_each_entity_run
}

void AuraGhosting::change_ghosting(BulkData& bulkData,
                                   EntityProcMapping& sendAuraEntityProcs)
{
  std::vector<EntityProc>& sendAuraGhosts = m_scratchSpace;
  sendAuraEntityProcs.fill_vec(sendAuraGhosts);

  //------------------------------------
  // Add the specified entities and their closure to sendAuraEntityProcs

  impl::StoreInEntityProcMapping storeEntity(bulkData, sendAuraEntityProcs);
  impl::NotAlreadyShared entityBelongsInAura(bulkData);
  for ( const EntityProc& entityProc : sendAuraGhosts ) {
    entityBelongsInAura.proc = entityProc.second;
    storeEntity.proc = entityProc.second;
    const EntityRank entityRank = bulkData.entity_rank(entityProc.first);
    if (entityRank > stk::topology::NODE_RANK) {
      if (entityRank > stk::topology::ELEM_RANK) {
        VisitClosureGeneral(bulkData, entityProc.first, entityRank, storeEntity, entityBelongsInAura);
      }
      else {
        VisitClosureBelowEntityNoRecurse(bulkData, entityProc.first, entityRank, storeEntity, entityBelongsInAura);
      }
    }
  }

  std::vector<EntityProc>& nonOwnedSendAuraGhosts = m_scratchSpace;
  nonOwnedSendAuraGhosts.clear();
  sendAuraEntityProcs.visit_entity_procs(
    [&bulkData,&nonOwnedSendAuraGhosts](Entity ent, int p)
    {
      if (!bulkData.bucket(ent).owned()) {
        nonOwnedSendAuraGhosts.emplace_back(ent,p);
      }
    });

  impl::comm_sync_nonowned_sends(bulkData, nonOwnedSendAuraGhosts, sendAuraEntityProcs);

  //------------------------------------
  // Remove send-ghost entities from the comm-list that no longer need to be sent.

  bool removed = false ;
  const unsigned auraGhostingOrdinal = bulkData.aura_ghosting().ordinal();

  const EntityCommDatabase& commDB = bulkData.internal_comm_db();
  EntityCommInfoVector comm_ghost ;
  for ( EntityCommListInfoVector::reverse_iterator
        i = bulkData.m_entity_comm_list.rbegin() ; i != bulkData.m_entity_comm_list.rend() ; ++i) {

    if (i->entity_comm == -1) {
      continue;
    }

    EntityCommListInfo& entityComm = *i;
    PairIterEntityComm commInfo = ghost_info_range(commDB.comm(entityComm.entity_comm), auraGhostingOrdinal);
    if (commInfo.empty()) {
      continue;
    }

    const bool is_owner = bulkData.parallel_owner_rank(entityComm.entity) == bulkData.parallel_rank() ;
    if ( is_owner ) {
      // Is owner, potentially removing ghost-sends
      // Have to make a copy

      comm_ghost.clear();
      for(; !commInfo.empty(); ++commInfo) {
        if (commInfo->ghost_id == auraGhostingOrdinal) {
          comm_ghost.push_back(*commInfo);
        }
      }

      if (sendAuraEntityProcs.get_num_procs(entityComm.entity) == 0) {
        for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
          const EntityCommInfo tmp = comm_ghost.back();
          bulkData.entity_comm_map_erase(entityComm.key, tmp);
        }
      }
      else {
        for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
          const EntityCommInfo tmp = comm_ghost.back();

          if (!sendAuraEntityProcs.find(entityComm.entity, tmp.proc) ) {
            bulkData.entity_comm_map_erase(entityComm.key, tmp);
          }
          else {
            sendAuraEntityProcs.eraseEntityProc(entityComm.entity, tmp.proc);
          }
        }
      }
    }

    if ( bulkData.internal_entity_comm_map(entityComm.entity).empty() ) {
      removed = true ;
      entityComm.key = EntityKey(); // No longer communicated
    }
  }

  // if an entry in the comm_list has the EntityKey() value, it is invalid,
  // and removed from the comm_list

  if ( removed ) {
    bulkData.delete_unneeded_entries_from_the_comm_list();
  }

  const std::vector<std::pair<EntityKey,EntityCommInfo>>& allRemovedGhosts = bulkData.m_removedGhosts;
  std::vector<EntityProc> removedSendGhosts;
  removedSendGhosts.reserve(allRemovedGhosts.size());
  for(const std::pair<EntityKey,EntityCommInfo>& rmGhost : allRemovedGhosts) {
    Entity rmEnt = bulkData.get_entity(rmGhost.first);
    if (bulkData.is_valid(rmEnt) &&
        rmGhost.second.ghost_id == auraGhostingOrdinal &&
        bulkData.parallel_owner_rank(rmEnt) == bulkData.parallel_rank() &&
        !sendAuraEntityProcs.find(rmEnt, rmGhost.second.proc)) {
      removedSendGhosts.push_back(EntityProc(rmEnt,rmGhost.second.proc));
    }
  }
  EntityLess entityLess(bulkData);
  {
    EntityProcVec().swap(sendAuraGhosts);
  }
  sendAuraEntityProcs.swap_vec(sendAuraGhosts);
  sendAuraEntityProcs.deallocate();
  stk::util::sort_and_unique(sendAuraGhosts, entityLess);

  const bool isFullRegen = true;
  bulkData.ghost_entities_and_fields(bulkData.aura_ghosting(), std::move(sendAuraGhosts), isFullRegen, removedSendGhosts);
}

}}} // end namepsace stk mesh impl

