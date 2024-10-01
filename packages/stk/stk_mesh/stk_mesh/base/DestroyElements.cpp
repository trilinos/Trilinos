#include "DestroyElements.hpp"
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>

namespace stk {
namespace mesh {

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk, stk::mesh::Entity connectedEntity, stk::mesh::EntityRank conRank)
{
  EntityVector scratchSpace;
  impl::destroy_upward_connected_aura_entities(bulk, connectedEntity, scratchSpace);
}

void pack_ghosts(BulkData &bulk, EntityVector &elemsAndRelatedEntities, stk::CommSparse& commSparse)
{
  std::vector<int> commProcs;
  for(int phase = 0; phase < 2; ++phase) {
    for(Entity ent : elemsAndRelatedEntities) {
      if (!bulk.is_valid(ent)) {
        continue;
      }

      bulk.comm_procs(ent, commProcs);
      for(int p : commProcs) {
        commSparse.send_buffer(p).pack<EntityKey>(bulk.entity_key(ent));
      }
    }

    if (phase == 0) {
      commSparse.allocate_buffers();
    }
  }
}

void unpack_recv_ghosts(stk::mesh::BulkData& bulk, stk::CommSparse& commSparse, EntityVector& recvGhostsToRemove)
{
  for(int p=0; p<commSparse.parallel_size(); ++p) {
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      EntityKey key;
      buf.unpack<EntityKey>(key);
      Entity entity = bulk.get_entity(key);
      if (bulk.is_valid(entity) && bulk.in_receive_ghost(entity)) {
        recvGhostsToRemove.push_back(entity);
      }
    }
  }
}

void remove_ghosts_from_remote_procs(stk::mesh::BulkData &bulk, EntityVector& recvGhostsToRemove)
{
  if (!recvGhostsToRemove.empty()) {
    stk::util::sort_and_unique(recvGhostsToRemove, EntityLess(bulk));

    impl::StoreEntity storeEntity(bulk);
    impl::VisitUpwardClosure(bulk, recvGhostsToRemove.begin(), recvGhostsToRemove.end(), storeEntity);
    storeEntity.store_visited_entities_in_vec(recvGhostsToRemove);

    stk::util::sort_and_unique(recvGhostsToRemove, EntityLess(bulk));
  }

  std::vector<EntityProc> emptyAdd;
  EntityVector removesForThisGhosting;
  removesForThisGhosting.reserve(recvGhostsToRemove.size());
  const bool notAddingSendGhosts = true;

  const std::vector<Ghosting*>& ghostings = bulk.ghostings();

  for(unsigned ig=0; ig<ghostings.size()-1; ++ig) {
    const unsigned reverseIdx = ghostings.size() - 1 - ig;
    Ghosting* ghosting = ghostings[reverseIdx];
    removesForThisGhosting.clear();
    for(unsigned i=0; i<recvGhostsToRemove.size(); ++i) {
      Entity ent = recvGhostsToRemove[i];
      if (bulk.is_valid(ent) && bulk.in_receive_ghost(*ghosting, ent)) {
        removesForThisGhosting.push_back(ent);
      }
    }

    bulk.internal_change_ghosting(*ghosting, emptyAdd, removesForThisGhosting, notAddingSendGhosts);
  }
}

void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy)
{
    stk::mesh::Selector orphansToDelete(bulk.mesh_meta_data().universal_part());
    destroy_elements(bulk, elementsToDestroy, orphansToDelete);
}

void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, const stk::mesh::Selector& orphansToDelete)
{
    bulk.modification_begin();
    bulk.m_bucket_repository.set_remove_mode_tracking();
    destroy_elements_no_mod_cycle(bulk, elementsToDestroy, orphansToDelete);
    bulk.modification_end();
    bulk.m_bucket_repository.set_remove_mode_fill_and_sort();
}

void get_all_related_entities(BulkData& bulk, EntityVector& elements, const Selector& orphansToDelete, EntityVector& relatedEntities)
{
  impl::StoreEntity storeEntity(bulk);

  auto ifSelected = [&](Entity ent) { return bulk.is_valid(ent) && orphansToDelete(bulk.bucket(ent)); };

  impl::VisitClosureGeneral(bulk, elements.begin(), elements.end(), storeEntity, ifSelected);
  storeEntity.store_visited_entities_in_vec(relatedEntities);

  auto ifSharedOrRecvGhost = [&](Entity ent) { return bulk.is_valid(ent) && (bulk.in_shared(ent) || bulk.in_receive_ghost(ent)); };

  const EntityRank endRank = static_cast<EntityRank>(bulk.mesh_meta_data().entity_rank_count());
  impl::VisitUpwardClosureGeneral(bulk, relatedEntities.begin(), relatedEntities.end(), endRank, storeEntity, ifSharedOrRecvGhost);
  storeEntity.store_visited_entities_in_vec(relatedEntities);
  relatedEntities.insert(relatedEntities.end(), elements.begin(), elements.end());

  stk::util::sort_and_unique(relatedEntities, stk::mesh::EntityLess(bulk));
}

void destroy_elements_no_mod_cycle(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, const stk::mesh::Selector& orphansToDelete)
{
  for(stk::mesh::Entity element : elementsToDestroy) {
    if(!bulk.is_valid(element))
        continue;

    STK_ThrowRequireMsg(!impl::has_upward_connectivity(bulk, element), "Element to be destroyed cannot have upward connectivity");
    STK_ThrowRequireMsg(bulk.entity_rank(element) == stk::topology::ELEM_RANK, "Entity to be destroyed must be an element");
  }

  stk::mesh::EntityVector elemsAndRelatedEntities;
  elemsAndRelatedEntities.reserve(2*elementsToDestroy.size());
  get_all_related_entities(bulk, elementsToDestroy, orphansToDelete, elemsAndRelatedEntities);

  stk::CommSparse commSparse(bulk.parallel());

  pack_ghosts(bulk, elemsAndRelatedEntities, commSparse);

  commSparse.communicate();

  EntityVector localEntitiesToRemove;
  EntityVector recvGhostsToRemove;
  for(Entity ent : elemsAndRelatedEntities) {
    if (bulk.is_valid(ent)) {
      if (bulk.in_receive_ghost(ent)) {
        recvGhostsToRemove.push_back(ent);
      }
      else {
        localEntitiesToRemove.push_back(ent);
      }
    }
  }

  unpack_recv_ghosts(bulk, commSparse, recvGhostsToRemove);

  remove_ghosts_from_remote_procs(bulk, recvGhostsToRemove);

  for(unsigned i=0; i<localEntitiesToRemove.size(); ++i) {
    const unsigned reverseIdx = localEntitiesToRemove.size() - 1 - i;
    Entity ent = localEntitiesToRemove[reverseIdx];
    if (bulk.is_valid(ent)) {
      bulk.destroy_entity(ent);
    }
  }
}

}}
