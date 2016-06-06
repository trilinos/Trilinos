#include "DestroyElements.hpp"
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

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

void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy)
{
    bulk.modification_begin();
    stk::mesh::EntityVector downwardConnectivity;

    for(stk::mesh::Entity element : elementsToDestroy) {
        if(!bulk.is_valid(element))
            continue;

        ThrowRequireMsg(!has_upward_connectivity(bulk, element), "Element to be destroyed cannot have upward connectivity");
        ThrowRequireMsg(bulk.entity_rank(element) == stk::topology::ELEM_RANK, "Entity to be destroyed must be an element");

        for(stk::mesh::EntityRank conRank = stk::topology::NODE_RANK; conRank < stk::topology::ELEM_RANK; ++conRank) {
            unsigned numConnected = bulk.num_connectivity(element, conRank);
            const stk::mesh::Entity* connectedEntities = bulk.begin(element, conRank);
            for(unsigned i = 0; i < numConnected; i++)
                downwardConnectivity.push_back(connectedEntities[i]);
        }

        bulk.destroy_entity(element);
    }

    stk::util::sort_and_unique(downwardConnectivity, stk::mesh::EntityLess(bulk));

    int numConnections = downwardConnectivity.size();
    for(int i = numConnections-1; i >= 0; --i)
    {
        stk::mesh::Entity connectedEntity =  downwardConnectivity[i];
        if(!has_upward_connectivity(bulk, connectedEntity))
            bulk.destroy_entity(connectedEntity);
    }
    bulk.modification_end();
}

}}
