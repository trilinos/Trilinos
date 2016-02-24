
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

struct InteriorSidesetFilter
{
    virtual ~InteriorSidesetFilter() {};

    virtual bool operator()(const stk::mesh::BulkData& bulkData,
                            stk::mesh::Entity element,
                            const std::vector<stk::mesh::PartOrdinal>& otherElementOrdinals) const
    {
        return true;
    }

    virtual bool operator()(const stk::mesh::BulkData& bulkData,
                            stk::mesh::Entity element,
                            stk::mesh::Entity otherElement) const
    {
        return true;
    }
};

struct InteriorBlockBoundaryFilter : public InteriorSidesetFilter
{
    virtual bool operator()(const stk::mesh::BulkData& bulkData,
                            stk::mesh::Entity element,
                            const std::vector<stk::mesh::PartOrdinal>& otherElementOrdinals) const
    {
        return !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, otherElementOrdinals);
    }

    virtual bool operator()(const stk::mesh::BulkData& bulkData,
                            stk::mesh::Entity element,
                            stk::mesh::Entity otherElement) const
    {
        return !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, element, otherElement);
    }
};

SkinMeshUtil::SkinMeshUtil(const ElemElemGraph& elemElemGraph,
                           const stk::mesh::PartVector& skinParts,
                           const stk::mesh::Selector& inputSkinSelector,
                           const stk::mesh::Selector* inputAirSelector)
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
    const stk::mesh::impl::ParallelInfo &parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    return parInfo.is_considered_air();
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

void SkinMeshUtil::mark_local_connections(const stk::mesh::GraphEdge &graphEdge,
                                           std::vector<bool> &isOnlyConnectedRemotely)
{
    if(impl::is_local_element(graphEdge.elem2))
    {
        stk::mesh::Entity other_element = eeGraph.get_entity_from_local_id(graphEdge.elem2);
        if(is_element_selected_and_can_have_side(skinSelector, other_element))
            isOnlyConnectedRemotely[graphEdge.side1] = false;
    }
}

void SkinMeshUtil::mark_remote_connections(const stk::mesh::GraphEdge &graphEdge,
                                            std::vector<bool> &isConnectedToRemoteElementInBodyToSkin)
{
    if(!impl::is_local_element(graphEdge.elem2))
    {
        const impl::ParallelInfo &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
        bool is_other_element_selected = parallel_edge_info.is_in_body_to_be_skinned();
        if(is_other_element_selected)
            isConnectedToRemoteElementInBodyToSkin[graphEdge.side1] = true;
    }
}

void SkinMeshUtil::mark_sides_exposed_on_other_procs(const stk::mesh::GraphEdge &graphEdge,
                                                      std::vector<bool> &isConnectedToRemoteElementInBodyToSkin,
                                                      std::vector<bool> &isOnlyConnectedRemotely)
{
    mark_local_connections(graphEdge, isOnlyConnectedRemotely);
    mark_remote_connections(graphEdge, isConnectedToRemoteElementInBodyToSkin);
}

std::vector<int> SkinMeshUtil::get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId, int numElemSides)
{
    std::vector<bool> isConnectedToRemoteElementInBodyToSkin(numElemSides, false);
    std::vector<bool> isOnlyConnectedRemotely(numElemSides, true);
    for(const stk::mesh::GraphEdge &graphEdge : eeGraph.get_edges_for_element(localId))
        mark_sides_exposed_on_other_procs(graphEdge, isConnectedToRemoteElementInBodyToSkin, isOnlyConnectedRemotely);

    for(const stk::mesh::GraphEdge graphEdge : eeGraph.get_coincident_edges_for_element(localId))
        mark_local_connections(graphEdge, isOnlyConnectedRemotely);

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
        else if(!eeGraph.get_coincident_edges_for_element(localId).empty())
        {
            const std::vector<GraphEdge>& edges = eeGraph.get_coincident_edges_for_element(localId);
            for(const stk::mesh::GraphEdge& edge : edges)
            {
                if(!stk::mesh::impl::is_local_element(edge.elem2))
                {
                    const impl::ParallelInfo& pinfo = eeGraph.get_parallel_info_for_graph_edge(edge);
                    if(pinfo.is_in_body_to_be_skinned())
                    {
                        exposedSides = get_exposed_sides(localId, maxSidesThisElement);
                        break;
                    }
                }
            }
        }
    }
    return exposedSides;
}

void SkinMeshUtil::extract_skinned_sideset_for_element(stk::mesh::Entity element, const stk::mesh::Bucket &bucket, std::vector<SideSetEntry> &skinnedSideSet)
{
    stk::mesh::impl::LocalId localId = eeGraph.get_local_element_id(element);
    std::vector<int> exposedSides = get_sides_for_skinning(skinSelector, bucket, element, localId, airSelector);

    for(size_t k=0; k<exposedSides.size(); ++k)
    {
        skinnedSideSet.push_back(SideSetEntry(element, static_cast<ConnectivityOrdinal> (exposedSides[k])));
    }
}

void SkinMeshUtil::extract_skinned_sideset_for_bucket(const stk::mesh::Bucket &bucket, std::vector<SideSetEntry> &skinnedSideSet)
{
    for(size_t j=0; j<bucket.size(); ++j)
    {
        stk::mesh::Entity element = bucket[j];
        extract_skinned_sideset_for_element(element, bucket, skinnedSideSet);
    }
}

std::vector<SideSetEntry> SkinMeshUtil::extract_skinned_sideset()
{
    std::vector<SideSetEntry> skinnedSideSet;

    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());

    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        extract_skinned_sideset_for_bucket(bucket, skinnedSideSet);
    }

    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));

    return skinnedSideSet;
}


bool SkinMeshUtil::is_parallel_graph_edge_candidate_for_interior_sideset(stk::mesh::Entity element, const stk::mesh::GraphEdge & graphEdge, const InteriorSidesetFilter& filter, bool &isElement1InSelector)
{
    const impl::ParallelInfo &parallel_edge_info = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
    bool isElement2InSelector = parallel_edge_info.is_in_body_to_be_skinned();
    if(!isElement1InSelector && !isElement2InSelector) return false;
    return filter(eeGraph.get_mesh(), element, parallel_edge_info.get_part_ordinals());
}

bool SkinMeshUtil::is_local_graph_edge_candidate_for_interior_sideset(stk::mesh::Entity element, const stk::mesh::GraphEdge & graphEdge, const InteriorSidesetFilter& filter, bool &isElement1InSelector)
{
    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
    stk::mesh::Entity otherElement = eeGraph.get_entity_from_local_id(graphEdge.elem2);
    bool isElement2InSelector = skinSelector(bulkData.bucket(otherElement));

    if(!isElement1InSelector && !isElement2InSelector) return false;

    if(bulkData.identifier(element) < bulkData.identifier(otherElement))
        return filter(bulkData, element, otherElement);

    return false;
}

bool SkinMeshUtil::is_graph_edge_candidate_for_interior_sideset(stk::mesh::Entity element, const stk::mesh::Bucket &bucket, const InteriorSidesetFilter& filter, const stk::mesh::GraphEdge & graphEdge)
{
    bool isElement1InSelector = skinSelector(bucket);

    if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
        return is_parallel_graph_edge_candidate_for_interior_sideset(element, graphEdge, filter, isElement1InSelector);

    return is_local_graph_edge_candidate_for_interior_sideset(element, graphEdge, filter, isElement1InSelector);
}

void SkinMeshUtil::extract_interior_sideset_for_element(stk::mesh::Entity element, const stk::mesh::Bucket &bucket, const InteriorSidesetFilter& filter, std::vector<SideSetEntry> &sideSet)
{
    impl::LocalId elementId = eeGraph.get_local_element_id(element);

    for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elementId))
    {
        bool should_add_side = is_graph_edge_candidate_for_interior_sideset(element, bucket, filter, graphEdge);

        if(should_add_side)
        {
            sideSet.push_back(SideSetEntry(element, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side1)));
            if(stk::mesh::impl::is_local_element(graphEdge.elem2))
            {
                stk::mesh::Entity otherElement = eeGraph.get_entity_from_local_id(graphEdge.elem2);
                sideSet.push_back(SideSetEntry(otherElement, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side2)));
            }
        }
    }
}

void SkinMeshUtil::extract_interior_sideset_for_bucket(const stk::mesh::Bucket &bucket, const InteriorSidesetFilter& filter, std::vector<SideSetEntry> &sideSet)
{
    for(size_t i=0;i<bucket.size();++i)
    {
        stk::mesh::Entity element = bucket[i];
        extract_interior_sideset_for_element(element, bucket, filter, sideSet);
    }
}

std::vector<SideSetEntry> SkinMeshUtil::extract_interior_sideset()
{
    std::vector<SideSetEntry> skinnedSideSet;
    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    InteriorBlockBoundaryFilter filter;

    for(const stk::mesh::Bucket* bucket : buckets) {
        extract_interior_sideset_for_bucket(*bucket, filter, skinnedSideSet);
    }

    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));
    return skinnedSideSet;
}

std::vector<SideSetEntry> SkinMeshUtil::extract_all_sides_sideset()
{
    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();

    unsigned maxNumSides = eeGraph.get_graph().get_num_edges();
    std::vector<SideSetEntry> sideSet;
    sideSet.reserve(maxNumSides);

    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    InteriorSidesetFilter filter;

    for (const stk::mesh::Bucket* bucketPtr : buckets) {
        const stk::mesh::Bucket & bucket = *bucketPtr;
        for (size_t i=0; i<bucket.size(); ++i) {
            extract_skinned_sideset_for_bucket(bucket, sideSet);
            extract_interior_sideset_for_bucket(bucket, filter, sideSet);
        }
    }

    stk::util::sort_and_unique(sideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));
    return sideSet;
}

}
}



