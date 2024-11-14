#include <stk_mesh/baseImpl/ElementTopologyDeletions.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

ElementTopologyDeletions::ElementTopologyDeletions(stk::mesh::BulkData &bulkData, stk::topology topologyToDelete) :
    bulk(bulkData),
    topoSelector(bulk.mesh_meta_data().get_topology_root_part(topologyToDelete)),
    selectsOtherElemTopologies(get_all_other_element_topologies_selector(bulk.mesh_meta_data(), topologyToDelete)),
    bucketsToDelete(bulk.get_buckets(stk::topology::ELEM_RANK, topoSelector))
{
}

void ElementTopologyDeletions::find_buckets_to_destroy_and_relations_to_sever()
{
    for(EntityRank irank = stk::topology::NODE_RANK; irank <= stk::topology::FACE_RANK; ++irank)
    {
        for(stk::mesh::Bucket *b : bulk.get_buckets(irank, topoSelector))
        {
            if(!selectsOtherElemTopologies(*b))
                bucketsToDelete.push_back(b);
            else
                append_upward_relations_to_sever(*b);
        }
    }
}

const std::vector<RelationEntityToNode> & ElementTopologyDeletions::get_relations_to_sever() const
{
    return relationsToDestroy;
}

const stk::mesh::BucketVector & ElementTopologyDeletions::get_buckets_to_delete() const
{
    return bucketsToDelete;
}

stk::mesh::Selector ElementTopologyDeletions::get_all_other_element_topologies_selector(const stk::mesh::MetaData &meta,
                                                                                        stk::topology topology)
{
    stk::mesh::Selector otherTopoSelector;
    for(stk::mesh::Part * part : meta.get_parts())
        if(is_topology_root_part(*part) && part->primary_entity_rank() == stk::topology::ELEM_RANK)
            if(meta.get_topology(*part) != topology)
                otherTopoSelector |= *part;
    return otherTopoSelector;
}

void ElementTopologyDeletions::append_upward_relations_to_sever(stk::mesh::Bucket &bucket)
{
    for(Entity entity : bucket)
    {
        for(EntityRank conRank = static_cast<EntityRank>(bucket.entity_rank() + 1); conRank <= stk::topology::ELEM_RANK;
                ++conRank)
        {
            unsigned numConnected = bulk.num_connectivity(entity, conRank);
            const Entity* connectedEntities = bulk.begin(entity, conRank);
            const ConnectivityOrdinal* ordinals = bulk.begin_ordinals(entity, conRank);
            for(unsigned i = 0; i < numConnected; i++)
                if(topoSelector(bulk.bucket(connectedEntities[i])) && !selectsOtherElemTopologies(bulk.bucket(connectedEntities[i])))
                    relationsToDestroy.push_back( {connectedEntities[i], entity, ordinals[i]});
        }
    }
}

}
}
}
