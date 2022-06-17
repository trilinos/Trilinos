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
#include <stk_mesh/baseImpl/Visitors.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>

namespace stk {
namespace mesh {
namespace impl {

AuraGhosting::AuraGhosting()
{
}

AuraGhosting::~AuraGhosting()
{
}

void AuraGhosting::generate_aura(BulkData& bulkData)
{
  EntityProcMapping entitySharing(bulkData.get_size_of_entity_index_space());
  std::vector<EntityRank> ranks = {stk::topology::NODE_RANK, stk::topology::EDGE_RANK};
  const MetaData& meta = bulkData.mesh_meta_data();
  if (meta.side_rank() > stk::topology::EDGE_RANK) {
    ranks.push_back(meta.side_rank());
  }
  std::vector<int> sharingProcs;
  for(EntityRank rank : ranks) {
    impl::for_each_selected_entity_run_no_threads(bulkData, rank, meta.globally_shared_part(),
      [&entitySharing, &sharingProcs](const BulkData& bulk, const MeshIndex& meshIndex) {
        Entity entity = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
        bulk.comm_shared_procs(entity, sharingProcs);
        for(int p : sharingProcs) {
          entitySharing.addEntityProc(entity, p);
        }
      });  
  }

  EntityProcMapping sendAuraEntityProcs(bulkData.get_size_of_entity_index_space());
  fill_send_aura_entities(bulkData, sendAuraEntityProcs, entitySharing);

  change_ghosting(bulkData, sendAuraEntityProcs, entitySharing);
}

void AuraGhosting::remove_aura(BulkData& bulkData)
{
  EntityVector recvAuraEntitiesToRemove;
  bulkData.aura_ghosting().receive_list(recvAuraEntitiesToRemove);
  bulkData.internal_change_ghosting(bulkData.aura_ghosting(), {}, recvAuraEntitiesToRemove);
}

void AuraGhosting::fill_send_aura_entities(BulkData& bulkData,
                                           EntityProcMapping& sendAuraEntityProcs,
                                           const EntityProcMapping& entitySharing)
{
  const EntityRank end_rank = static_cast<EntityRank>(bulkData.mesh_meta_data().entity_rank_count());

  // Iterate over all shared entities, ensure that upwardly related
  // entities to each shared entity will be ghosted to the sharing proc.
  Selector shared = bulkData.mesh_meta_data().globally_shared_part();

  std::vector<int> sharingProcs;
  impl::for_each_selected_entity_run_no_threads(bulkData, stk::topology::NODE_RANK, shared,
    [&sendAuraEntityProcs, &entitySharing, &sharingProcs, &end_rank]
    (const BulkData& bulk, const MeshIndex& meshIndex) {
      const Bucket& bucket = *meshIndex.bucket;
      const unsigned bucketOrd = meshIndex.bucket_ordinal;
      const EntityRank nextHigherRank = stk::topology::EDGE_RANK;

      bulk.comm_shared_procs(bucket[bucketOrd], sharingProcs);
      for (const int sharingProc : sharingProcs) {

        for (EntityRank higherRank = nextHigherRank; higherRank < end_rank; ++higherRank) {
          const unsigned num_rels = bucket.num_connectivity(bucketOrd, higherRank);
          const Entity* rels     = bucket.begin(bucketOrd, higherRank);

          for (unsigned r = 0; r < num_rels; ++r) {
            stk::mesh::impl::insert_upward_relations(bulk, entitySharing, rels[r], stk::topology::NODE_RANK, sharingProc, sendAuraEntityProcs);
          }
        }
      }    
    }    
  ); // for_each_entity_run
}

void AuraGhosting::change_ghosting(BulkData& bulkData,
                                   EntityProcMapping& sendAuraEntityProcs,
                                   const EntityProcMapping& entitySharing)
{
  std::vector<EntityProc> add_send;
  sendAuraEntityProcs.fill_vec(add_send);

  //------------------------------------
  // Add the specified entities and their closure to sendAuraEntityProcs

  impl::StoreInEntityProcMapping siepm(bulkData, sendAuraEntityProcs);
  EntityProcMapping epm(bulkData.get_size_of_entity_index_space());
  impl::OnlyGhostsEPM og(bulkData, epm, entitySharing);
  for ( const EntityProc& entityProc : add_send ) {
      og.proc = entityProc.second;
      siepm.proc = entityProc.second;
      impl::VisitClosureGeneral(bulkData,entityProc.first,siepm,og);
  }

  sendAuraEntityProcs.fill_vec(add_send);

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  std::vector<bool> ghostStatus(bulkData.get_size_of_entity_index_space(), false);

  stk::mesh::impl::comm_sync_aura_send_recv(bulkData, add_send,
                                            sendAuraEntityProcs, ghostStatus );

  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  OrdinalVector addParts;
  OrdinalVector removeParts(1, bulkData.m_ghost_parts[BulkData::AURA]->mesh_meta_data_ordinal());
  OrdinalVector scratchOrdinalVec, scratchSpace;
  bool removed = false ;

  std::vector<EntityCommInfo> comm_ghost ;
  for ( EntityCommListInfoVector::reverse_iterator
        i = bulkData.m_entity_comm_list.rbegin() ; i != bulkData.m_entity_comm_list.rend() ; ++i) {

    if (!i->entity_comm) {
      continue;
    }

    EntityCommListInfo& entityComm = *i;
    if (!entityComm.entity_comm->isGhost) {
      continue;
    }

    const bool is_owner = bulkData.parallel_owner_rank(entityComm.entity) == bulkData.parallel_rank() ;
    const bool remove_recv = ( ! is_owner ) &&
                             !ghostStatus[entityComm.entity.local_offset()] && bulkData.in_receive_ghost(bulkData.aura_ghosting(), entityComm.entity);

    if(bulkData.is_valid(entityComm.entity))
    {
      if ( is_owner ) {
        // Is owner, potentially removing ghost-sends
        // Have to make a copy

          const PairIterEntityComm ec = ghost_info_range(entityComm.entity_comm->comm_map, bulkData.aura_ghosting());
          comm_ghost.assign( ec.first , ec.second );

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
      else if ( remove_recv ) {
          bulkData.entity_comm_map_erase(entityComm.key, bulkData.aura_ghosting());
          bulkData.internal_change_entity_parts(entityComm.entity, addParts, removeParts, scratchOrdinalVec, scratchSpace);
      }

      if ( bulkData.internal_entity_comm_map(entityComm.entity).empty() ) {
        removed = true ;
        entityComm.key = EntityKey(); // No longer communicated
        if ( remove_recv ) {
          ThrowRequireMsg( bulkData.internal_destroy_entity_with_notification( entityComm.entity, remove_recv ),
                           "P[" << bulkData.parallel_rank() << "]: FAILED attempt to destroy entity: "
                           << bulkData.entity_key(entityComm.entity) );
        }
      }
    }
  }

  // if an entry in the comm_list has the EntityKey() value, it is invalid,
  // and removed from the comm_list

  if ( removed ) {
    bulkData.delete_unneeded_entries_from_the_comm_list();
  }

  EntityLess entityLess(bulkData);
  std::set<EntityProc , EntityLess> finalSendGhosts(entityLess);
  sendAuraEntityProcs.fill_set(finalSendGhosts);

  const bool isFullRegen = true;
  bulkData.ghost_entities_and_fields(bulkData.aura_ghosting(), finalSendGhosts, isFullRegen);
}

}}} // end namepsace stk mesh impl

