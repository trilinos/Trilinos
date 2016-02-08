
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/FaceCreator.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>

namespace stk
{
namespace mesh
{


SkinMeshUtil::SkinMeshUtil(stk::mesh::BulkData& bulkData, const ElemElemGraph& elemElemGraph, const stk::mesh::PartVector& skinParts, const stk::mesh::Selector& inputSkinSelector, const stk::mesh::Selector* inputAirSelector)
: eeGraph(elemElemGraph), skinSelector(inputSkinSelector), airSelector(inputAirSelector)
{
}

std::vector<int> SkinMeshUtil::get_exposed_sides(stk::mesh::impl::LocalId localId, int maxSidesThisElement)
{
    std::vector<int> exposedSides;
    impl::add_exposed_sides(localId, maxSidesThisElement, eeGraph.get_graph(), exposedSides);
    if(airSelector != nullptr)
        add_exposed_sides_due_to_air_selector(localId, exposedSides);
    return exposedSides;
}

void SkinMeshUtil::add_exposed_sides_due_to_air_selector(impl::LocalId local_id, std::vector<int> &exposedSides)
{
    std::vector<int> exposedSidesDueToAirSelector;
    stk::mesh::Entity element = eeGraph.get_entity_from_local_id(local_id);
    std::vector<bool> areAllConnectedElemsAir(eeGraph.get_mesh().bucket(element).topology().num_sides(), true);

    const GraphEdgesForElement& graphEdges = eeGraph.get_edges_for_element(local_id);
    for(const GraphEdge & graphEdge : graphEdges)
    {
        if(is_connected_element_air(graphEdge))
            exposedSidesDueToAirSelector.push_back(graphEdge.side1);
        else
            areAllConnectedElemsAir[graphEdge.side1] = false;
    }

    for(int exposedSide : exposedSidesDueToAirSelector)
        if(areAllConnectedElemsAir[exposedSide])
            exposedSides.push_back(exposedSide);
}

bool SkinMeshUtil::is_remote_element_air(const ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge)
{
    const stk::mesh::impl::parallel_info &parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    return parInfo.m_is_air;
}

bool SkinMeshUtil::is_connected_element_air(const stk::mesh::GraphEdge &graphEdge)
{
    if(impl::is_local_element(graphEdge.elem2))
        return is_element_selected_and_can_have_side(*airSelector, eeGraph.get_entity_from_local_id(graphEdge.elem2));
    else
        return is_remote_element_air(eeGraph.get_parallel_info_for_graph_edges(), graphEdge);
}

bool SkinMeshUtil::is_element_selected_and_can_have_side(const stk::mesh::Selector& selector, stk::mesh::Entity otherElement)
{
    return selector(eeGraph.get_mesh().bucket(otherElement)) && impl::does_element_have_side(eeGraph.get_mesh(), otherElement);
}

std::vector<int> SkinMeshUtil::get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId, int numElemSides)
{
    std::vector<bool> isConnectedToRemoteElementInBodyToSkin(numElemSides, false);
    std::vector<bool> isOnlyConnectedRemotely(numElemSides, true);
    for(const stk::mesh::GraphEdge &graphEdge : eeGraph.get_edges_for_element(localId))
    {
        if(impl::is_local_element(graphEdge.elem2))
        {
            stk::mesh::Entity other_element = eeGraph.get_entity_from_local_id(graphEdge.elem2);
            if(is_element_selected_and_can_have_side(skinSelector, other_element))
                isOnlyConnectedRemotely[graphEdge.side1] = false;
        }
        else
        {
            const impl::parallel_info &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
            bool is_other_element_selected = parallel_edge_info.m_in_body_to_be_skinned;
            if(is_other_element_selected)
                isConnectedToRemoteElementInBodyToSkin[graphEdge.side1] = true;
        }
    }

    std::vector<int> exposedSides;
    for(int side = 0; side < numElemSides; side++)
    {
        if(isConnectedToRemoteElementInBodyToSkin[side] && isOnlyConnectedRemotely[side])
        {
            exposedSides.push_back(side);
        }
    }
    return exposedSides;
}

std::vector<int> SkinMeshUtil::get_sides_for_skinning(const stk::mesh::Selector& skinSelector,
                                        const stk::mesh::Bucket& bucket,
                                        stk::mesh::Entity element,
                                        stk::mesh::impl::LocalId localId,
                                        const stk::mesh::Selector* airSelector)
{
    int maxSidesThisElement = bucket.topology().num_sides();
    std::vector<int> exposedSides;
    if(stk::mesh::impl::does_element_have_side(bucket.mesh(), element))
    {
        if(skinSelector(bucket))
            exposedSides = get_exposed_sides(localId, maxSidesThisElement);
        else if(airSelector != nullptr && (*airSelector)(bucket))
            exposedSides = get_sides_exposed_on_other_procs(localId, maxSidesThisElement);
    }
    return exposedSides;
}

std::vector<SideSetEntry> SkinMeshUtil::extract_skinned_sideset()
{
    std::vector<SideSetEntry> skinnedSideSet;

    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());

    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();++j)
        {
            stk::mesh::Entity element = bucket[j];

            stk::mesh::impl::LocalId localId = eeGraph.get_local_element_id(element);
            std::vector<int> exposedSides = get_sides_for_skinning(skinSelector, bucket, element, localId, airSelector);

            for(size_t k=0; k<exposedSides.size(); ++k)
            {
                skinnedSideSet.push_back(SideSetEntry(element, static_cast<ConnectivityOrdinal> (exposedSides[k])));
            }
        }
    }

    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));

    return skinnedSideSet;
}

std::vector<SideSetEntry> SkinMeshUtil::extract_interior_sideset()
{
    std::vector<SideSetEntry> skinnedSideSet;
    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());

    for(const stk::mesh::Bucket* bucket : buckets)
    {
        if(bucket->owned())
        {
            for(size_t i=0;i<bucket->size();++i)
            {
                stk::mesh::Entity element = (*bucket)[i];
                impl::LocalId elementId = eeGraph.get_local_element_id(element);
                ;

                for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elementId))
                {
                    stk::mesh::EntityId otherEntityId = 0;
                    stk::mesh::Entity otherElement;
                    bool isElement1InSelector = skinSelector(bucket);
                    bool isElement2InSelector = false;

                    bool isParallelEdge = graphEdge.elem2<0;
                    bool should_add_side = false;

                    if(isParallelEdge)
                    {
                        const impl::parallel_info &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
                        isElement2InSelector = parallel_edge_info.m_in_body_to_be_skinned;
                        if(!isElement1InSelector && !isElement2InSelector) continue;
                        should_add_side = !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, parallel_edge_info.m_part_ordinals);
                    }
                    else
                    {
                        otherElement = eeGraph.get_entity_from_local_id(graphEdge.elem2);
                        otherEntityId = bulkData.identifier(otherElement);
                        isElement2InSelector = skinSelector(bulkData.bucket(otherElement));
                        if(!isElement1InSelector && !isElement2InSelector) continue;
                        if(!isParallelEdge && bulkData.identifier(element) < otherEntityId)
                            should_add_side = !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, otherElement);
                    }

                    if(should_add_side)
                    {
                        skinnedSideSet.push_back(SideSetEntry(element, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side1)));
                        if(!isParallelEdge)
                            skinnedSideSet.push_back(SideSetEntry(otherElement, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side2)));
                    }
                }
            }
        }
    }

    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));

    return skinnedSideSet;
}

}
}



