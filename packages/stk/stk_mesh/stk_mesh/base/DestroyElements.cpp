#include "DestroyElements.hpp"
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

bool has_upward_connectivity(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity)
{
    if(!bulk.is_valid(entity))
        return false;

    const stk::mesh::Bucket &bucket = bulk.bucket(entity);

    for(stk::mesh::EntityRank conRank = static_cast<stk::mesh::EntityRank>(bucket.entity_rank() + 1); conRank <= stk::topology::CONSTRAINT_RANK; ++conRank)
    {
        unsigned numConnected = bulk.num_connectivity(entity, conRank);
        if(numConnected > 0)
            return true;
    }

    return false;
}

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk, stk::mesh::Entity connectedEntity, stk::mesh::EntityRank conRank)
{
    for(stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count()-1); rank > conRank; rank--)
    {
        unsigned numUpward = bulk.num_connectivity(connectedEntity, rank);
        const stk::mesh::Entity *entities = bulk.begin(connectedEntity, rank);
        for (unsigned j=0;j<numUpward;++j)
        {
            stk::mesh::Entity upwardEntity = entities[numUpward-j-1];
            if(bulk.is_valid(upwardEntity) && bulk.bucket(upwardEntity).in_aura())
            {
                destroy_upward_connected_aura_entities(bulk, upwardEntity, rank);
                bulk.destroy_entity(upwardEntity);
            }
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
    std::vector<stk::mesh::EntityVector> downwardConnectivity(bulk.mesh_meta_data().entity_rank_count());

    for(stk::mesh::Entity element : elementsToDestroy) {
        if(!bulk.is_valid(element))
            continue;

        ThrowRequireMsg(!has_upward_connectivity(bulk, element), "Element to be destroyed cannot have upward connectivity");
        ThrowRequireMsg(bulk.entity_rank(element) == stk::topology::ELEM_RANK, "Entity to be destroyed must be an element");

        for(stk::mesh::EntityRank conRank = stk::topology::NODE_RANK; conRank < stk::topology::ELEM_RANK; ++conRank) {
            unsigned numConnected = bulk.num_connectivity(element, conRank);
            const stk::mesh::Entity* connectedEntities = bulk.begin(element, conRank);
            for(unsigned i = 0; i < numConnected; i++)
                downwardConnectivity[conRank].push_back(connectedEntities[i]);
        }

        bulk.destroy_entity(element);
    }

    for(stk::mesh::EntityRank conRank : {stk::topology::FACE_RANK, stk::topology::EDGE_RANK, stk::topology::NODE_RANK})
    {
        stk::util::sort_and_unique(downwardConnectivity[conRank], stk::mesh::EntityLess(bulk));

        int numConnections = downwardConnectivity[conRank].size();
        for(int i = numConnections-1; i >= 0; --i)
        {
            stk::mesh::Entity connectedEntity =  downwardConnectivity[conRank][i];
            if(orphansToDelete(bulk.bucket(connectedEntity)))
            {
                destroy_upward_connected_aura_entities(bulk, connectedEntity, conRank);
                if(!has_upward_connectivity(bulk, connectedEntity))
                    bulk.destroy_entity(connectedEntity);
            }
        }
    }
    bulk.modification_end();
}

}}
