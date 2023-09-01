#include "DestroyElements.hpp"
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>

namespace stk {
namespace mesh {

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk, stk::mesh::Entity connectedEntity, stk::mesh::EntityRank conRank,
                                            EntityVector& scratchSpace)
{
  impl::StoreEntity storeEntity(bulk);
  impl::VisitUpwardClosure(bulk, connectedEntity, storeEntity);

  storeEntity.store_visited_entities_in_vec(scratchSpace);
  stk::util::sort_and_unique(scratchSpace, EntityLess(bulk));

  for(unsigned i=0; i<scratchSpace.size(); ++i) {
    int reverseIdx = scratchSpace.size() - 1 - i;
    Entity upwardEntity = scratchSpace[reverseIdx];

    if (bulk.is_valid(upwardEntity) && bulk.bucket(upwardEntity).in_aura()) {
      bulk.destroy_entity(upwardEntity);
    }
  }
}

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk, stk::mesh::Entity connectedEntity, stk::mesh::EntityRank conRank)
{
  EntityVector scratchSpace;
  destroy_upward_connected_aura_entities(bulk, connectedEntity, conRank, scratchSpace);
}

void remove_ghosts_on_remote_procs(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy,
                   std::vector<stk::mesh::EntityVector>& downwardConnectivity,
                   const Selector& orphansToDelete)
{
  stk::CommSparse commSparse(bulk.parallel());

  std::vector<int> commProcs;
  for(int phase = 0; phase < 2; ++phase) {
    for(Entity elem : elementsToDestroy) {
      if (!bulk.is_valid(elem) || !bulk.bucket(elem).owned()) {
        continue;
      }

      bulk.comm_procs(elem, commProcs);
      for(int p : commProcs) {
        if (!bulk.in_shared(elem, p)) {
          commSparse.send_buffer(p).pack<EntityKey>(bulk.entity_key(elem));
        }
      }
    }

    for(EntityVector& downward : downwardConnectivity) {
      for(Entity entity : downward) {
        if (bulk.is_valid(entity) && bulk.bucket(entity).owned() && orphansToDelete(bulk.bucket(entity))) {
          bulk.comm_procs(entity, commProcs);
          for(int p : commProcs) {
            if (!bulk.in_shared(entity, p)) {
              commSparse.send_buffer(p).pack<EntityKey>(bulk.entity_key(entity));
            }
          }
        }
      }
    }

    if (phase == 0) {
      commSparse.allocate_buffers();
    }
  }

  commSparse.communicate();

  EntityVector recvGhostsToRemove;
  for(int p=0; p<commSparse.parallel_size(); ++p) {
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      EntityKey key;
      buf.unpack<EntityKey>(key);
      Entity entity = bulk.get_entity(key);
      if (bulk.in_receive_ghost(entity)) {
        recvGhostsToRemove.push_back(entity);
      }
    }
  }

  if (!recvGhostsToRemove.empty()) {
    stk::util::sort_and_unique(recvGhostsToRemove, EntityLess(bulk));
    for(unsigned i=0; i<recvGhostsToRemove.size(); ++i) {
      const unsigned reverseIdx = recvGhostsToRemove.size() - 1 - i;
      bulk.destroy_entity(recvGhostsToRemove[reverseIdx]);
    }
  }
}

void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy)
{
    stk::mesh::Selector orphansToDelete(bulk.mesh_meta_data().universal_part());
    destroy_elements(bulk, elementsToDestroy, orphansToDelete);
}

void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, stk::mesh::Selector orphansToDelete)
{
    bulk.modification_begin();
    destroy_elements_no_mod_cycle(bulk, elementsToDestroy, orphansToDelete);
    bulk.modification_end();
}

void destroy_elements_no_mod_cycle(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, stk::mesh::Selector orphansToDelete)
{
    std::vector<stk::mesh::EntityVector> downwardConnectivity(bulk.mesh_meta_data().entity_rank_count());

    for(stk::mesh::Entity element : elementsToDestroy) {
        if(!bulk.is_valid(element))
            continue;

        STK_ThrowRequireMsg(!impl::has_upward_connectivity(bulk, element), "Element to be destroyed cannot have upward connectivity");
        STK_ThrowRequireMsg(bulk.entity_rank(element) == stk::topology::ELEM_RANK, "Entity to be destroyed must be an element");

        for(stk::mesh::EntityRank conRank = stk::topology::NODE_RANK; conRank < stk::topology::ELEM_RANK; ++conRank) {
            unsigned numConnected = bulk.num_connectivity(element, conRank);
            const stk::mesh::Entity* connectedEntities = bulk.begin(element, conRank);
            for(unsigned i = 0; i < numConnected; i++)
                downwardConnectivity[conRank].push_back(connectedEntities[i]);
        }

    }

    remove_ghosts_on_remote_procs(bulk, elementsToDestroy, downwardConnectivity, orphansToDelete);

    for(stk::mesh::Entity element : elementsToDestroy) {
        bulk.destroy_entity(element);
    }

    EntityVector scratchSpace;
    for(stk::mesh::EntityRank conRank : {stk::topology::FACE_RANK, stk::topology::EDGE_RANK, stk::topology::NODE_RANK})
    {
        stk::util::sort_and_unique(downwardConnectivity[conRank], stk::mesh::EntityLess(bulk));

        int numConnections = downwardConnectivity[conRank].size();
        for(int i = numConnections-1; i >= 0; --i)
        {
            stk::mesh::Entity connectedEntity =  downwardConnectivity[conRank][i];
            if(bulk.is_valid(connectedEntity) && orphansToDelete(bulk.bucket(connectedEntity)))
            {
                if (bulk.in_shared(connectedEntity)) {
                  destroy_upward_connected_aura_entities(bulk, connectedEntity, conRank, scratchSpace);
                }
                if(!impl::has_upward_connectivity(bulk, connectedEntity))
                    bulk.destroy_entity(connectedEntity);
            }
        }
    }
}

}}
