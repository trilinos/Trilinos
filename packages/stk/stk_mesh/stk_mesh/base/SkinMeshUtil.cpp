
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetEntryCompare.hpp>
#include <stk_mesh/base/FaceCreator.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{

SkinMeshUtil::SkinMeshUtil(ElemElemGraph& elemElemGraph,
                           const stk::mesh::Selector& inputSkinSelector,
                           const stk::mesh::Selector* inputAirSelector)
: eeGraph(elemElemGraph), skinSelector(inputSkinSelector), useAirSelector(inputAirSelector != nullptr)
{
    if (useAirSelector) airSelector = *inputAirSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(eeGraph.get_mesh(), eeGraph, skinSelector, remoteSkinSelector);
    if (useAirSelector)
        impl::populate_selected_value_for_remote_elements(eeGraph.get_mesh(), eeGraph, airSelector, remoteAirSelector);
}

std::vector<SideSetEntry> SkinMeshUtil::get_skinned_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector)
{
  bulk.initialize_face_adjacent_element_graph();
  ElemElemGraph& elemElemGraph = bulk.get_face_adjacent_element_graph();
  SkinMeshUtil skinMesh(elemElemGraph, skinSelector, nullptr);
  return skinMesh.extract_skinned_sideset();
}

std::vector<SideSetEntry> SkinMeshUtil::get_skinned_sideset_excluding_region(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector, const stk::mesh::Selector& exclusionRegionSelector)
{
  bulk.initialize_face_adjacent_element_graph();
  ElemElemGraph& elemElemGraph = bulk.get_face_adjacent_element_graph();
  SkinMeshUtil skinMesh(elemElemGraph, skinSelector, &exclusionRegionSelector);
  return skinMesh.extract_skinned_sideset();
}

std::vector<SideSetEntry> SkinMeshUtil::get_interior_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector)
{
  bulk.initialize_face_adjacent_element_graph();
  ElemElemGraph& elemElemGraph = bulk.get_face_adjacent_element_graph();
  SkinMeshUtil skinMesh(elemElemGraph, skinSelector, nullptr);
  return skinMesh.extract_interior_sideset();
}

std::vector<SideSetEntry> SkinMeshUtil::get_all_sides_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector, bool includeAuraElementSides)
{
  bulk.initialize_face_adjacent_element_graph();
  ElemElemGraph& elemElemGraph = bulk.get_face_adjacent_element_graph();
  SkinMeshUtil skinMesh(elemElemGraph, skinSelector, nullptr);
  return skinMesh.extract_all_sides_sideset(includeAuraElementSides);
}

void SkinMeshUtil::get_exposed_sides(stk::mesh::impl::LocalId localId, int maxSidesThisElement,
                                     std::vector<int>& exposedSides)
{
    exposedSides.clear();
    impl::add_exposed_sides(localId, maxSidesThisElement, eeGraph.get_graph(), exposedSides);
    if(useAirSelector)
        add_exposed_sides_due_to_air_selector(localId, exposedSides);
}

void SkinMeshUtil::add_exposed_sides_due_to_air_selector(impl::LocalId local_id, std::vector<int> &exposedSides)
{
    std::vector<int> exposedSidesDueToAirSelector;
    stk::mesh::Entity element = eeGraph.get_entity(local_id);
    std::vector<bool> areAllConnectedElemsAir(eeGraph.get_mesh().bucket(element).topology().num_sides(), true);

    GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(local_id);
    for(const GraphEdge & graphEdge : graphEdges)
    {
        if(is_connected_element_air(graphEdge))
            exposedSidesDueToAirSelector.push_back(graphEdge.side1());
        else
            areAllConnectedElemsAir[graphEdge.side1()] = false;
    }

    for(int exposedSide : exposedSidesDueToAirSelector)
        if(areAllConnectedElemsAir[exposedSide])
            exposedSides.push_back(exposedSide);
}

bool SkinMeshUtil::is_remote_element_air(const ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge)
{
    return remoteAirSelector[graphEdge.elem2()];
}

bool SkinMeshUtil::is_connected_element_air(const stk::mesh::GraphEdge &graphEdge)
{
    if(impl::is_local_element(graphEdge.elem2()))
        return is_element_selected_and_can_have_side(airSelector, eeGraph.get_entity(graphEdge.elem2()));
    else
        return is_remote_element_air(eeGraph.get_parallel_info_for_graph_edges(), graphEdge);
}

bool SkinMeshUtil::is_element_selected_and_can_have_side(const stk::mesh::Selector& selector, stk::mesh::Entity otherElement)
{
    return selector(eeGraph.get_mesh().bucket(otherElement)) && impl::does_element_have_side(eeGraph.get_mesh(), otherElement);
}

void SkinMeshUtil::mark_local_connections(const stk::mesh::GraphEdge &graphEdge,
                                          std::vector<int> &isOnlyConnectedRemotely)
{
    if(impl::is_local_element(graphEdge.elem2()))
    {
        stk::mesh::Entity other_element = eeGraph.get_entity(graphEdge.elem2());
        if(is_element_selected_and_can_have_side(skinSelector, other_element))
            isOnlyConnectedRemotely[graphEdge.side1()] = false;
    }
}

void SkinMeshUtil::mark_remote_connections(const stk::mesh::GraphEdge &graphEdge,
                                           std::vector<int> &isConnectedToRemoteElementInBodyToSkin)
{
    if(!impl::is_local_element(graphEdge.elem2()))
    {
        bool is_other_element_selected = remoteSkinSelector[graphEdge.elem2()];
        if(is_other_element_selected)
            isConnectedToRemoteElementInBodyToSkin[graphEdge.side1()] = true;
    }
}

void SkinMeshUtil::mark_sides_exposed_on_other_procs(const stk::mesh::GraphEdge &graphEdge,
                                                     std::vector<int> &isConnectedToRemoteElementInBodyToSkin,
                                                     std::vector<int> &isOnlyConnectedRemotely)
{
    mark_local_connections(graphEdge, isOnlyConnectedRemotely);
    mark_remote_connections(graphEdge, isConnectedToRemoteElementInBodyToSkin);
}

void SkinMeshUtil::get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId,
                                                    int numElemSides,
                                                    std::vector<int>& exposedSides)
{
    m_isConnectedToRemoteElem.assign(numElemSides, false);
    m_isOnlyConnectedRemotely.assign(numElemSides, true);
    for(const stk::mesh::GraphEdge &graphEdge : eeGraph.get_edges_for_element(localId))
        mark_sides_exposed_on_other_procs(graphEdge, m_isConnectedToRemoteElem, m_isOnlyConnectedRemotely);

    for(const stk::mesh::GraphEdge& graphEdge : eeGraph.get_coincident_edges_for_element(localId))
        mark_local_connections(graphEdge, m_isOnlyConnectedRemotely);

    exposedSides.clear();
    for(int side = 0; side < numElemSides; side++)
    {
        if(m_isConnectedToRemoteElem[side] && m_isOnlyConnectedRemotely[side])
        {
            exposedSides.push_back(side);
        }
    }
}


void SkinMeshUtil::get_sides_for_skinning(const stk::mesh::Bucket& bucket,
                                          stk::mesh::Entity element,
                                          stk::mesh::impl::LocalId localId,
                                          std::vector<int>& exposedSides)
{
    stk::topology elemTopo = bucket.topology();
    int maxSidesThisElement = elemTopo.num_sides();
    exposedSides.clear();
    const bool elemHasSide = elemTopo.dimension() == bucket.mesh().mesh_meta_data().spatial_dimension();
    if(elemHasSide)
    {
        if(skinSelector(bucket))
            get_exposed_sides(localId, maxSidesThisElement, exposedSides);
        else if(useAirSelector && airSelector(bucket))
            get_sides_exposed_on_other_procs(localId, maxSidesThisElement, exposedSides);
        else if(!eeGraph.get_coincident_edges_for_element(localId).empty())
        {
            const std::vector<GraphEdge>& edges = eeGraph.get_coincident_edges_for_element(localId);
            for(const stk::mesh::GraphEdge& edge : edges)
            {
                if(!stk::mesh::impl::is_local_element(edge.elem2()))
                {
                    if(remoteSkinSelector[edge.elem2()])
                    {
                        get_exposed_sides(localId, maxSidesThisElement, exposedSides);
                        break;
                    }
                }
            }
        }
    }
}

bool checkIfSideIsNotCollapsed(stk::mesh::EntityVector& sideNodes, const stk::mesh::Bucket& bucket, const stk::mesh::BulkData& bulkData, stk::mesh::Entity element, int sideOrdinal)
{
    unsigned dim = bulkData.mesh_meta_data().spatial_dimension();
    if(dim==1) return true;

    stk::topology sideTopo = bucket.topology().sub_topology(bulkData.mesh_meta_data().side_rank(), sideOrdinal);
    if (sideTopo.num_vertices() < dim) {
      return false; //side is "collapsed", possibly because this is an edge-side of a 3D shell...
    }

    sideNodes.resize(bucket.topology().sub_topology(bulkData.mesh_meta_data().side_rank(), sideOrdinal).num_nodes());
    bucket.topology().side_nodes(bulkData.begin_nodes(element), sideOrdinal, sideNodes.data());
    stk::util::sort_and_unique(sideNodes);
    return sideNodes.size() >= dim;
}

std::vector<SideSetEntry> SkinMeshUtil::extract_skinned_sideset()
{
    std::vector<SideSetEntry> skinnedSideSet;

    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();

    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());

    stk::mesh::EntityVector sideNodes;
    std::vector<int> exposedSides;

    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();++j)
        {
            stk::mesh::Entity element = bucket[j];

            stk::mesh::impl::LocalId localId = eeGraph.get_local_element_id(element);
            get_sides_for_skinning(bucket, element, localId, exposedSides);

            for(size_t k=0; k<exposedSides.size(); ++k)
            {
                if(checkIfSideIsNotCollapsed(sideNodes, bucket, bulkData, element, exposedSides[k]))
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

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    stk::mesh::impl::populate_part_ordinals_for_remote_edges(bulkData, eeGraph, parallelPartInfo);
    stk::mesh::EntityVector sideNodes;
    std::vector<stk::mesh::PartOrdinal> elemPartOrdinals, scratchOrdinals2;

    stk::mesh::Selector ownedSkinSelector = bulkData.mesh_meta_data().locally_owned_part() & skinSelector;
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, ownedSkinSelector);
    for(const stk::mesh::Bucket* bucket : buckets)
    {
        for(stk::mesh::Entity element : *bucket)
        {
            impl::LocalId elementId = eeGraph.get_local_element_id(element);
            stk::mesh::impl::get_element_block_part_ordinals(element, bulkData, elemPartOrdinals);

            for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elementId))
            {
                stk::mesh::EntityId otherEntityId = 0;
                stk::mesh::Entity otherElement;
                bool isElement2InSelector = false;

                bool isParallelEdge = !impl::is_local_element(graphEdge.elem2());
                bool should_add_side = false;

                if(isParallelEdge)
                {
                    isElement2InSelector = remoteSkinSelector[graphEdge.elem2()];
                    if(!isElement2InSelector) continue;
                    should_add_side = !(elemPartOrdinals == parallelPartInfo[graphEdge.elem2()].elementPartOrdinals);
                }
                else
                {
                    otherElement = eeGraph.get_entity(graphEdge.elem2());
                    otherEntityId = bulkData.identifier(otherElement);
                    if(bulkData.identifier(element) < otherEntityId && skinSelector(bulkData.bucket(otherElement))) {
                        should_add_side = !stk::mesh::impl::are_entity_element_blocks_equivalent(bulkData, otherElement, elemPartOrdinals, scratchOrdinals2);
                    }
                }

                if(should_add_side && checkIfSideIsNotCollapsed(sideNodes, *bucket, bulkData, element, graphEdge.side1()))
                {
                    skinnedSideSet.emplace_back(element, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side1()));
                    if(!isParallelEdge)
                        skinnedSideSet.emplace_back(otherElement, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side2()));
                }
            }
        }
    }

    stk::util::sort_and_unique(skinnedSideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));

    return skinnedSideSet;
}

void get_aura_element_sides_for_skinning(const stk::mesh::BulkData& bulk, stk::mesh::Entity element,
                                 std::vector<int>& exposedSides, std::vector<stk::mesh::Entity>& elems)
{
    exposedSides.clear();
    stk::topology topo = bulk.bucket(element).topology();
    const stk::mesh::Entity* elemNodes = bulk.begin_nodes(element);

    const unsigned maxNodesPerSide = 128;
    stk::mesh::Entity sideNodes[maxNodesPerSide];
    const unsigned numSides = topo.num_sides();
    for(unsigned side = 0; side<numSides; ++side) {
        topo.side_nodes(elemNodes, side, sideNodes);
        unsigned numNodesThisSide = topo.side_topology(side).num_nodes();
        stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEM_RANK,
                                   numNodesThisSide, sideNodes, elems);
        if (elems.size() == 1) {
            exposedSides.push_back(side);
        }
        else {
            bool connectedElementHasSmallerId = false;
            bool connectedElementIsOwned = false;

            for(unsigned elemIndex=0; elemIndex<elems.size(); ++elemIndex) {
                if (element != elems[elemIndex]) {
                    if (bulk.bucket(elems[elemIndex]).owned()) {
                        connectedElementIsOwned = true;
                        break;
                    }
                    if (bulk.identifier(elems[elemIndex]) < bulk.identifier(element)) {
                        connectedElementHasSmallerId = true;
                        break;
                    }
                }
            }

            if (!connectedElementHasSmallerId && !connectedElementIsOwned) {
                exposedSides.push_back(side);
            }
        }
    }
}

std::vector<SideSetEntry> SkinMeshUtil::extract_all_sides_sideset(bool includeAuraElementSides)
{
    const stk::mesh::BulkData& bulkData = eeGraph.get_mesh();

    unsigned maxNumSides = eeGraph.get_graph().get_num_edges();
    std::vector<SideSetEntry> sideSet;
    sideSet.reserve(maxNumSides);

    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    stk::mesh::EntityVector sideNodes;
    stk::mesh::EntityVector elems;
    std::vector<int> exposedSides;

    for (const stk::mesh::Bucket* bucketPtr : buckets) {
        const stk::mesh::Bucket & bucket = *bucketPtr;
        for (size_t i=0; i<bucket.size(); ++i) {
            stk::mesh::Entity element = bucket[i];
            impl::LocalId elementId = eeGraph.get_local_element_id(element);

            get_sides_for_skinning(bucket, element, elementId, exposedSides);

            for (size_t k=0; k<exposedSides.size(); ++k)
            {
                if(checkIfSideIsNotCollapsed(sideNodes, bucket, bulkData, element, exposedSides[k]))
                        sideSet.push_back(SideSetEntry(element, static_cast<ConnectivityOrdinal> (exposedSides[k])));
            }

            for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elementId)) {
                stk::mesh::EntityId otherEntityId = 0;
                stk::mesh::Entity otherElement;
                bool isElement1InSelector = skinSelector(bucket);
                bool isElement2InSelector = false;

                bool isParallelEdge = !impl::is_local_element(graphEdge.elem2());
                bool should_add_side = false;

                if (isParallelEdge) {
                    isElement2InSelector = remoteSkinSelector[graphEdge.elem2()];
                    if (!isElement1InSelector && !isElement2InSelector) continue;
                    should_add_side = true;
                    if (includeAuraElementSides && isElement2InSelector) {
                        stk::mesh::EntityId otherId = eeGraph.convert_negative_local_id_to_global_id(graphEdge.elem2());
                        stk::mesh::Entity auraElement = bulkData.get_entity(stk::topology::ELEM_RANK, otherId);
                        if (bulkData.is_valid(auraElement) && bulkData.bucket(auraElement).in_aura()) {
                            get_aura_element_sides_for_skinning(bulkData, auraElement, exposedSides, elems);
                            for (size_t k=0; k<exposedSides.size(); ++k)
                            {
                                 sideSet.push_back(SideSetEntry(auraElement, static_cast<ConnectivityOrdinal>(exposedSides[k])));
                            }
                        }
                    }
                }
                else {
                    otherElement = eeGraph.get_entity(graphEdge.elem2());
                    otherEntityId = bulkData.identifier(otherElement);
                    isElement2InSelector = skinSelector(bulkData.bucket(otherElement));
                    if (!isElement1InSelector && !isElement2InSelector) continue;
                    if (bulkData.identifier(element) < otherEntityId) {
                        should_add_side = true;
                    }
                }

                if(should_add_side && checkIfSideIsNotCollapsed(sideNodes, bucket, bulkData, element, graphEdge.side1()))
                {
                    sideSet.push_back(SideSetEntry(element, static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side1())));
                    if (!isParallelEdge) {
                        if (bulkData.identifier(otherElement) < bulkData.identifier(element)) {
                            sideSet.back().element = otherElement;
                            sideSet.back().side = static_cast<stk::mesh::ConnectivityOrdinal>(graphEdge.side2());
                        }
                    }
                }
            }
        }
    }

    if (includeAuraElementSides)
    {
        exposedSides.clear();
        const stk::mesh::BucketVector& auraBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().aura_part());
        for (const stk::mesh::Bucket* bucketPtr : auraBuckets)
        {
            const stk::mesh::Bucket & bucket = *bucketPtr;
            if (!skinSelector(bucket)) {
                continue;
            }

            for (size_t i=0; i<bucket.size(); ++i)
            {
                stk::mesh::Entity element = bucket[i];
                get_aura_element_sides_for_skinning(bulkData, element, exposedSides, elems);
                for (size_t k=0; k<exposedSides.size(); ++k)
                {
                     sideSet.push_back(SideSetEntry(element, static_cast<ConnectivityOrdinal>(exposedSides[k])));
                }
            }
        }
    }

    stk::util::sort_and_unique(sideSet, SideSetEntryLess(bulkData), SideSetEntryEquals(bulkData));

    return sideSet;
}

}
}
