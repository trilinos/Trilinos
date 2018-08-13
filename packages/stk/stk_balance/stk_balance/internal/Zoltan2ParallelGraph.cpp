#include "Zoltan2ParallelGraph.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/human_bytes.hpp>

void Zoltan2ParallelGraph::adjust_vertex_weights(const stk::balance::BalanceSettings& balanceSettings,
                           stk::mesh::BulkData& stkMeshBulkData,
                           const std::vector<stk::mesh::Selector>& selectors,
                           const stk::mesh::impl::LocalIdMapper& localIds)
{
    if(balanceSettings.areVertexWeightsProvidedInAVector())
    {
        size_t previousSize = mVertexWeights.size();
        mVertexWeights = balanceSettings.getVertexWeightsViaVector();
        size_t newSize = mVertexWeights.size();
        ThrowRequireWithSierraHelpMsg(newSize == previousSize);
    }
    else if(balanceSettings.areVertexWeightsProvidedViaFields())
    {
        stk::mesh::EntityVector entitiesToBalance;
        stk::mesh::get_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK), entitiesToBalance);

        fillFieldVertexWeights(balanceSettings, stkMeshBulkData, selectors, entitiesToBalance);
    }
    else if ( !stk::balance::internal::is_geometric_method(balanceSettings.getDecompMethod()) && balanceSettings.setVertexWeightsBasedOnNumberAdjacencies() )
    {
        const stk::mesh::BucketVector &buckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            const stk::mesh::Bucket &bucket = *buckets[i];
            if(bucket.owned() && bucket.topology() == stk::topology::PARTICLE)
            {
                for(size_t j = 0; j < bucket.size(); j++)
                {
                    unsigned localId = stk::balance::internal::get_local_id(localIds, bucket[j]);
                    unsigned numAdjElements = mOffsets[localId+1] - mOffsets[localId];
                    mVertexWeights[localId] = 0.5*numAdjElements;
                }
            }
        }
    }
}

void Zoltan2ParallelGraph::adjust_weights_for_small_meshes()
{
    for(size_t i=0;i<mVertexWeights.size();++i)
        mVertexWeights[i] += 0.05;
}

BalanceGlobalNumber getAdjacencyId(const stk::mesh::BulkData &stkMeshBulkData,
                                 stk::mesh::Entity adjacentVertex,
                                 bool useLocalIds,
                                 const stk::mesh::impl::LocalIdMapper& localIds)
{
    BalanceGlobalNumber adjacencyId = -1;
    if(useLocalIds)
    {
        adjacencyId = stk::balance::internal::get_local_id(localIds, adjacentVertex);
    }
    else
    {
        adjacencyId = stkMeshBulkData.identifier(adjacentVertex);
    }
    return adjacencyId;
}

bool areElementsConnected(const stk::mesh::BulkData& stkMeshBulkData,
                          const size_t numNodesRequiredForConnection,
                          const stk::mesh::Entity element1,
                          const stk::mesh::Entity element2)
{
    size_t nodesSharedBetweenElements = stk::balance::internal::getNumSharedNodesBetweenElements(stkMeshBulkData, element1, element2);
    bool elementsAreConnected = (nodesSharedBetweenElements >= numNodesRequiredForConnection);
    return elementsAreConnected;
}

std::vector<stk::mesh::Entity> getAllPossiblyConnectedElements(const stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element)
{
    unsigned numNodes = stkMeshBulkData.num_nodes(element);
    const stk::mesh::Entity* nodes = stkMeshBulkData.begin_nodes(element);
    std::vector<stk::mesh::Entity> allElementsPossiblyConnected;
    for (unsigned k = 0; k < numNodes; k++)
    {
        unsigned numElementsAttachedToNode = stkMeshBulkData.num_elements(nodes[k]);
        const stk::mesh::Entity* elements = stkMeshBulkData.begin_elements(nodes[k]);
        allElementsPossiblyConnected.insert(allElementsPossiblyConnected.end(), elements, elements+numElementsAttachedToNode);
    }
    stk::util::sort_and_unique(allElementsPossiblyConnected);
    return allElementsPossiblyConnected;
}

void createGraphEdgesForElement(stk::mesh::BulkData &stkMeshBulkData,
                                const stk::balance::BalanceSettings &balanceSettings,
                                stk::mesh::Entity elementOfConcern,
                                std::vector<stk::balance::GraphEdge> &graphEdges,
                                const stk::mesh::impl::LocalIdMapper& localIds)
{
    std::vector<stk::mesh::Entity> allElementsPossiblyConnected = getAllPossiblyConnectedElements(stkMeshBulkData, elementOfConcern);
    bool useLocalIds = balanceSettings.getGraphOption() == stk::balance::BalanceSettings::COLORING;

    for (stk::mesh::Entity possiblyConnectedElement : allElementsPossiblyConnected)
    {
        if (elementOfConcern != possiblyConnectedElement)
        {
            bool doingLoadBalancingConsiderAnyElement = (balanceSettings.getGraphOption() == stk::balance::BalanceSettings::LOADBALANCE);
            bool doingColoringConsiderOnlyOwnedElement = (balanceSettings.getGraphOption() == stk::balance::BalanceSettings::COLORING && stkMeshBulkData.bucket(possiblyConnectedElement).owned());
            if(doingLoadBalancingConsiderAnyElement || doingColoringConsiderOnlyOwnedElement)
            {
                stk::topology element1Topology = stkMeshBulkData.bucket(elementOfConcern).topology();
                stk::topology element2Topology = stkMeshBulkData.bucket(possiblyConnectedElement).topology();
                size_t numNodesRequiredForConnection = balanceSettings.getNumNodesRequiredForConnection(element1Topology, element2Topology);
                if ( areElementsConnected(stkMeshBulkData, numNodesRequiredForConnection, elementOfConcern, possiblyConnectedElement) &&
                     !stk::balance::internal::shouldOmitSpiderElement(stkMeshBulkData, balanceSettings, possiblyConnectedElement) )
                {
                    double edgeWeight = balanceSettings.getGraphEdgeWeight(element1Topology, element2Topology);
                    int vertex2ParallelOwner = 0;
                    if ( useLocalIds )
                    {
                        graphEdges.push_back(stk::balance::GraphEdge(elementOfConcern, stk::balance::internal::get_local_id(localIds, possiblyConnectedElement), vertex2ParallelOwner, edgeWeight));
                    }
                    else
                    {
                        vertex2ParallelOwner = stkMeshBulkData.parallel_owner_rank(possiblyConnectedElement);
                        graphEdges.push_back(stk::balance::GraphEdge(elementOfConcern, stkMeshBulkData.identifier(possiblyConnectedElement), vertex2ParallelOwner, edgeWeight));
                    }
                }
            }
        }
    }
}

void Zoltan2ParallelGraph::convertGraphEdgesToZoltanGraph(const stk::mesh::BulkData& stkMeshBulkData,
                                                          const std::vector<stk::balance::GraphEdge> &graphEdges,
                                                          const unsigned numElements,
                                                          std::vector<int>& adjacencyProcs,
                                                          const stk::mesh::impl::LocalIdMapper& localIds)
{
    mOffsets.clear();
    mAdjacency.clear();
    mEdgeWeights.clear();

    mOffsets.resize(numElements+1,0);

    // Need to resize
    for (size_t i=0;i<graphEdges.size();i++)
    {
        unsigned localId = stk::balance::internal::get_local_id(localIds, graphEdges[i].vertex1());
        mOffsets[localId+1]++;
        // Only need if graph is NOT symmetric here
        // offsets[graphEdges[i].vertex2().local_offset()]++;
    }

    for (size_t i=1;i<mOffsets.size();i++)
    {
        mOffsets[i] += mOffsets[i-1];
    }

    std::ostringstream os;
    size_t aboutToAllocate = mOffsets[numElements]*(sizeof(BalanceGlobalNumber) + sizeof(int) + sizeof(double));
    os << "About to allocate " << stk::human_bytes(aboutToAllocate) << " bytes";
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(), os.str());

    mAdjacency.resize(mOffsets[numElements]);
    adjacencyProcs.resize(mOffsets[numElements]);
    mEdgeWeights.resize(mOffsets[numElements]);

    std::vector<int> numElementsConnectedSoFar(numElements,0);

    for (size_t i=0;i<graphEdges.size();i++)
    {
        unsigned localId = stk::balance::internal::get_local_id(localIds, graphEdges[i].vertex1());
        int offset = numElementsConnectedSoFar[localId];
        int index = mOffsets[localId] + offset;

        mAdjacency[index] = graphEdges[i].vertex2();
        adjacencyProcs[index] = graphEdges[i].vertex2OwningProc();
        mEdgeWeights[index] = graphEdges[i].weight();
        numElementsConnectedSoFar[localId]++;
    }
}

void Zoltan2ParallelGraph::createGraphEdgesUsingNodeConnectivity(stk::mesh::BulkData &stkMeshBulkData,
                                           const stk::balance::BalanceSettings &balanceSettings,
                                           size_t numElements,
                                           std::vector<stk::balance::GraphEdge> &graphEdges,
                                           const stk::mesh::impl::LocalIdMapper& localIds)
{
    const stk::mesh::BucketVector &buckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);

    mVertexIds.clear();
    mVertexIds.resize(numElements);
    mVertexWeights.clear();
    mVertexWeights.resize(numElements, 0);
    mVertexCoordinates.clear();
    unsigned spatialDimension = stkMeshBulkData.mesh_meta_data().spatial_dimension();
    mVertexCoordinates.resize(numElements*spatialDimension, 0);

    stk::mesh::FieldBase const * coord = stk::balance::internal::get_coordinate_field(stkMeshBulkData.mesh_meta_data(),
                                                                               balanceSettings.getCoordinateFieldName());

    bool useLocalId = (balanceSettings.getGraphOption() == stk::balance::BalanceSettings::COLORING);
    for (size_t i=0; i<buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if(bucket.owned())
        {
            for (size_t j=0; j<bucket.size(); j++)
            {
                stk::mesh::Entity elementOfConcern = bucket[j];
                unsigned local_id = stk::balance::internal::get_local_id(localIds, elementOfConcern);
                mVertexIds[local_id] = getAdjacencyId(stkMeshBulkData, elementOfConcern, useLocalId, localIds);

                if (!stk::balance::internal::shouldOmitSpiderElement(stkMeshBulkData, balanceSettings, elementOfConcern))
                {
                    mVertexWeights[local_id] = balanceSettings.getGraphVertexWeight(bucket.topology());
                    stk::balance::internal::fillEntityCentroid(stkMeshBulkData, coord, elementOfConcern, &mVertexCoordinates[local_id*spatialDimension]);
                    createGraphEdgesForElement(stkMeshBulkData, balanceSettings, elementOfConcern, graphEdges, localIds);
                }
            }
        }
    }
}

void Zoltan2ParallelGraph::fillZoltan2AdapterDataFromStkMesh(stk::mesh::BulkData &stkMeshBulkData,
                                       const stk::balance::BalanceSettings &balanceSettings,
                                       std::vector<int>& adjacencyProcs,
                                       const stk::mesh::Selector& searchSelector,
                                       const stk::mesh::impl::LocalIdMapper& localIds)
{
    size_t numElements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                                    stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

    std::vector<stk::balance::GraphEdge> graphEdges;

    stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Create graph edges using node connectivity");

    // Set vertexWeights based on topology of entity
    createGraphEdgesUsingNodeConnectivity(stkMeshBulkData, balanceSettings, numElements,
                                          graphEdges, localIds);

    if ( !this->amCheckingForMechanisms() && balanceSettings.includeSearchResultsInGraph() )
    {
        stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Create graph edges using search results");

        stk::balance::internal::createGraphEdgesUsingBBSearch(stkMeshBulkData, balanceSettings, graphEdges, searchSelector);

        std::sort(graphEdges.begin(), graphEdges.end());
        std::vector<stk::balance::GraphEdge>::iterator iter = std::unique(graphEdges.begin(), graphEdges.end());
        graphEdges.erase(iter, graphEdges.end());

        std::vector<unsigned> localIdsSearch(graphEdges.size(),0);
        size_t counter = 0;
        for (size_t i=0;i<graphEdges.size();i++)
        {
            if ( graphEdges[i].isEdgeFromSearch() )
            {
                localIdsSearch[counter] = stk::balance::internal::get_local_id(localIds, graphEdges[i].vertex1());
                counter++;
            }
        }

        localIdsSearch.resize(counter);
        stk::util::sort_and_unique(localIdsSearch);

        // Adjust vertex weights if part of search
        for(size_t i=0;i<localIdsSearch.size();++i)
        {
            mVertexWeights[localIdsSearch[i]] *= balanceSettings.getVertexWeightMultiplierForVertexInSearch();
        }
    }

    stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Convert edges to a graph");

    convertGraphEdgesToZoltanGraph(stkMeshBulkData,
                                   graphEdges,
                                   numElements,
                                   adjacencyProcs,
                                   localIds);

    stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Finished creating graph");
}


