#include "Zoltan2ParallelGraph.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/human_bytes.hpp>

void Zoltan2ParallelGraph::adjust_vertex_weights(const stk::balance::BalanceSettings& balanceSettings,
                                                 stk::mesh::BulkData& stkMeshBulkData,
                                                 const stk::mesh::Selector& selector,
                                                 const stk::mesh::impl::LocalIdMapper& localIds)
{
  const stk::mesh::Part & locallyOwnedPart =  stkMeshBulkData.mesh_meta_data().locally_owned_part();

#ifndef STK_HIDE_DEPRECATED_CODE  // Delete after April 2021
  if (balanceSettings.areVertexWeightsProvidedInAVector())
  {
    size_t previousSize = mVertexWeights.size();
    mVertexWeights = balanceSettings.getVertexWeightsViaVector();
    size_t newSize = mVertexWeights.size();
    ThrowRequireWithSierraHelpMsg(newSize == previousSize);
  }
  else if (balanceSettings.areVertexWeightsProvidedViaFields())
#else
  if (balanceSettings.areVertexWeightsProvidedViaFields())
#endif
  {
    stk::mesh::EntityVector entitiesToBalance;
    const bool sortById = true;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::ELEM_RANK,
                            stkMeshBulkData.mesh_meta_data().locally_owned_part(), entitiesToBalance, sortById);

    fillFieldVertexWeights(balanceSettings, stkMeshBulkData, {selector}, entitiesToBalance);
  }
  else if (!stk::balance::internal::is_geometric_method(balanceSettings.getDecompMethod()) &&
           balanceSettings.setVertexWeightsBasedOnNumberAdjacencies())
  {
    const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, selector & locallyOwnedPart);
    for (size_t i = 0; i < buckets.size(); i++) {
      const stk::mesh::Bucket &bucket = *buckets[i];
      if (bucket.topology() == stk::topology::PARTICLE) {
        for (size_t j = 0; j < bucket.size(); j++) {
          unsigned localId = stk::balance::internal::get_local_id(localIds, bucket[j]);
          unsigned numAdjElements = mOffsets[localId+1] - mOffsets[localId];
          mVertexWeights[localId] = 0.5*numAdjElements;
        }
      }
    }
  }

  for (const auto & blockWeightMultiplier : balanceSettings.getVertexWeightBlockMultipliers()) {
    const std::string & blockName = blockWeightMultiplier.first;
    const double weightMultiplier = blockWeightMultiplier.second;

    const stk::mesh::Part * block = stkMeshBulkData.mesh_meta_data().get_part(blockName);
    ThrowRequireMsg(block != nullptr, "Mesh does not contain a block named '" + blockName + "'");

    const stk::mesh::BucketVector & buckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK,
                                                                          *block & selector & locallyOwnedPart);
    for (const stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & entity : *bucket) {
        const unsigned localId = stk::balance::internal::get_local_id(localIds, entity);
        mVertexWeights[localId] *= weightMultiplier;
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
                          const stk::mesh::Entity & element1,
                          const stk::mesh::Entity & element2)
{
    size_t nodesSharedBetweenElements = stk::balance::internal::get_num_common_nodes_between_elements(stkMeshBulkData, element1, element2);
    bool elementsAreConnected = (nodesSharedBetweenElements >= numNodesRequiredForConnection);
    return elementsAreConnected;
}

std::vector<stk::mesh::Entity> getAllPossiblyConnectedElements(const stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element)
{
  unsigned numNodes = stkMeshBulkData.num_nodes(element);
  const stk::mesh::Entity* nodes = stkMeshBulkData.begin_nodes(element);
  std::vector<stk::mesh::Entity> allElementsPossiblyConnected;
  for (unsigned k = 0; k < numNodes; ++k) {
    unsigned numElementsAttachedToNode = stkMeshBulkData.num_elements(nodes[k]);
    const stk::mesh::Entity* elements = stkMeshBulkData.begin_elements(nodes[k]);
    allElementsPossiblyConnected.insert(allElementsPossiblyConnected.end(), elements, elements+numElementsAttachedToNode);
  }
  stk::util::sort_and_unique(allElementsPossiblyConnected);
  return allElementsPossiblyConnected;
}

bool is_valid_graph_connectivity(const stk::mesh::BulkData &stkMeshBulkData,
                                 const stk::mesh::Selector& selector,
                                 const stk::balance::BalanceSettings &balanceSettings,
                                 const stk::mesh::Entity & elementOfConcern,
                                 const stk::mesh::Entity & possiblyConnectedElement)
{
  stk::topology element1Topology = stkMeshBulkData.bucket(elementOfConcern).topology();
  stk::topology element2Topology = stkMeshBulkData.bucket(possiblyConnectedElement).topology();
  size_t numNodesRequiredForConnection = balanceSettings.getNumNodesRequiredForConnection(element1Topology, element2Topology);

  return (areElementsConnected(stkMeshBulkData, numNodesRequiredForConnection, elementOfConcern, possiblyConnectedElement) &&
          selector(stkMeshBulkData.bucket(elementOfConcern)) &&
          selector(stkMeshBulkData.bucket(possiblyConnectedElement)));
}

stk::balance::GraphEdge create_graph_edge(const stk::mesh::BulkData &bulk,
                                          const stk::balance::BalanceSettings &balanceSettings,
                                          const stk::mesh::impl::LocalIdMapper& localIds,
                                          const stk::mesh::Entity & element1,
                                          const stk::mesh::Entity & element2)
{
  const stk::topology element1Topology = bulk.bucket(element1).topology();
  const stk::topology element2Topology = bulk.bucket(element2).topology();
  const stk::mesh::EntityId element2Id = stk::balance::internal::get_local_id(localIds, element2);
  double edgeWeight = balanceSettings.getGraphEdgeWeight(element1Topology, element2Topology);
  int vertex2ParallelOwner = 0;

  return stk::balance::GraphEdge(element1, element2Id, vertex2ParallelOwner, edgeWeight);
}

stk::balance::GraphEdge create_graph_edge(const stk::mesh::BulkData &bulk,
                                          const stk::balance::BalanceSettings &balanceSettings,
                                          const stk::mesh::Entity & element1,
                                          const stk::mesh::Entity & element2)
{
  const stk::topology element1Topology = bulk.bucket(element1).topology();
  const stk::topology element2Topology = bulk.bucket(element2).topology();
  const stk::mesh::EntityId element2Id = bulk.identifier(element2);
  double edgeWeight = balanceSettings.getGraphEdgeWeight(element1Topology, element2Topology);
  int vertex2ParallelOwner = bulk.parallel_owner_rank(element2);

  return stk::balance::GraphEdge(element1, element2Id, vertex2ParallelOwner, edgeWeight);
}

bool should_create_graph_edge(const stk::mesh::BulkData &bulk,
                              const stk::mesh::Selector& selector,
                              const stk::balance::BalanceSettings &balanceSettings,
                              const stk::mesh::Entity & element1,
                              const stk::mesh::Entity & element2)
{
  const bool isValidGraphEdge = is_valid_graph_connectivity(bulk, selector, balanceSettings, element1, element2);
  const bool notSpiderElement1 = !stk::balance::internal::shouldOmitSpiderElement(bulk, balanceSettings, element1);
  const bool notSpiderElement2 = !stk::balance::internal::shouldOmitSpiderElement(bulk, balanceSettings, element2);

  return isValidGraphEdge && notSpiderElement1 && notSpiderElement2;
}

void createGraphEdgesForElement(const stk::mesh::BulkData &bulk,
                                const stk::mesh::Selector& selector,
                                const stk::balance::BalanceSettings &balanceSettings,
                                const stk::mesh::Entity & elementOfConcern,
                                std::vector<stk::balance::GraphEdge> &graphEdges,
                                const stk::mesh::impl::LocalIdMapper& localIds)
{
  std::vector<stk::mesh::Entity> allElementsPossiblyConnected = getAllPossiblyConnectedElements(bulk, elementOfConcern);
  bool usingColoring = balanceSettings.usingColoring();

  if (usingColoring) {
    for (stk::mesh::Entity possiblyConnectedElement : allElementsPossiblyConnected) {
      if (elementOfConcern != possiblyConnectedElement) {
        bool considerOnlySelectedOwnedElement = (bulk.bucket(possiblyConnectedElement).owned() &&
                                                 selector(bulk.bucket(possiblyConnectedElement)));
        if (considerOnlySelectedOwnedElement) {
          if (should_create_graph_edge(bulk, selector, balanceSettings, elementOfConcern, possiblyConnectedElement)) {
            graphEdges.push_back(create_graph_edge(bulk, balanceSettings, localIds, elementOfConcern, possiblyConnectedElement));
          }
        }
      }
    }
  }
  else {
    for (stk::mesh::Entity possiblyConnectedElement : allElementsPossiblyConnected) {
      if (elementOfConcern != possiblyConnectedElement) {
        if (bulk.bucket(possiblyConnectedElement).owned()) {
          const stk::mesh::EntityId elementOfConcernId = bulk.identifier(elementOfConcern);
          const stk::mesh::EntityId possiblyConnectedElementId = bulk.identifier(possiblyConnectedElement);
          if (elementOfConcernId < possiblyConnectedElementId) {
            if (should_create_graph_edge(bulk, selector, balanceSettings, elementOfConcern, possiblyConnectedElement)) {
              graphEdges.push_back(create_graph_edge(bulk, balanceSettings, elementOfConcern, possiblyConnectedElement));
              graphEdges.push_back(create_graph_edge(bulk, balanceSettings, possiblyConnectedElement, elementOfConcern));
            }
          }
        }
        else {
          if (should_create_graph_edge(bulk, selector, balanceSettings, elementOfConcern, possiblyConnectedElement)) {
            graphEdges.push_back(create_graph_edge(bulk, balanceSettings, elementOfConcern, possiblyConnectedElement));
          }
        }
      }
    }
  }
}

void Zoltan2ParallelGraph::convertGraphEdgesToZoltanGraph(const stk::mesh::BulkData& stkMeshBulkData,
                                                          const stk::balance::BalanceSettings &balanceSettings,
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

    bool usingColoring = balanceSettings.usingColoring();
    for (size_t i=0;i<graphEdges.size();i++)
    {
        unsigned localId = stk::balance::internal::get_local_id(localIds, graphEdges[i].vertex1());
        int offset = numElementsConnectedSoFar[localId];
        int index = mOffsets[localId] + offset;

        if (usingColoring) {
            ThrowRequireMsg(graphEdges[i].vertex2_id() < numElements, "Adding invalid element to ZoltanII graph");
        }
        mAdjacency[index] = graphEdges[i].vertex2_id();
        adjacencyProcs[index] = graphEdges[i].vertex2_owning_proc();
        mEdgeWeights[index] = graphEdges[i].weight();
        numElementsConnectedSoFar[localId]++;
    }
}

void Zoltan2ParallelGraph::createGraphEdgesUsingNodeConnectivity(stk::mesh::BulkData &stkMeshBulkData,
                                                                 const stk::mesh::Selector& selector,
                                                                 const stk::balance::BalanceSettings &balanceSettings,
                                                                 size_t numElements,
                                                                 std::vector<stk::balance::GraphEdge> &graphEdges,
                                                                 const stk::mesh::impl::LocalIdMapper& localIds)
{
  bool useLocalIds = balanceSettings.usingColoring();
  const stk::mesh::BucketVector& buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, selector);

  mVertexIds.clear();
  mVertexIds.resize(numElements);
  mVertexWeights.clear();
  mVertexWeights.resize(numElements, 0);
  mVertexCoordinates.clear();
  unsigned spatialDimension = stkMeshBulkData.mesh_meta_data().spatial_dimension();
  mVertexCoordinates.resize(numElements*spatialDimension, 0);

  const stk::mesh::FieldBase * coord = stkMeshBulkData.mesh_meta_data().coordinate_field();

  for (size_t i = 0; i < buckets.size(); ++i) {
    const stk::mesh::Bucket &bucket = *buckets[i];
    if (bucket.owned()) {
      for (const stk::mesh::Entity & elementOfConcern : bucket) {
        unsigned local_id = stk::balance::internal::get_local_id(localIds, elementOfConcern);
        mVertexIds[local_id] = getAdjacencyId(stkMeshBulkData, elementOfConcern, useLocalIds, localIds);

        if (!stk::balance::internal::shouldOmitSpiderElement(stkMeshBulkData, balanceSettings, elementOfConcern)) {
          mVertexWeights[local_id] = balanceSettings.getGraphVertexWeight(bucket.topology());
          stk::balance::internal::fillEntityCentroid(stkMeshBulkData, coord, elementOfConcern, &mVertexCoordinates[local_id*spatialDimension]);
          createGraphEdgesForElement(stkMeshBulkData, selector, balanceSettings, elementOfConcern, graphEdges, localIds);
        }
      }
    }
  }
}

void Zoltan2ParallelGraph::fillZoltan2AdapterDataFromStkMesh(stk::mesh::BulkData &stkMeshBulkData,
                                                             const stk::balance::BalanceSettings &balanceSettings,
                                                             std::vector<int>& adjacencyProcs,
                                                             const stk::mesh::Selector& selector,
                                                             const stk::mesh::impl::LocalIdMapper& localIds)
{
  const stk::mesh::Selector localSelected = selector & stkMeshBulkData.mesh_meta_data().locally_owned_part();
  size_t numElements = stk::mesh::count_entities(stkMeshBulkData, stk::topology::ELEM_RANK, localSelected);

  std::vector<stk::balance::GraphEdge> graphEdges;

  stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Create graph edges using node connectivity");

  // Set vertexWeights based on topology of entity
  createGraphEdgesUsingNodeConnectivity(stkMeshBulkData, selector, balanceSettings, numElements,
                                        graphEdges, localIds);

  if (!this->amCheckingForMechanisms() && balanceSettings.includeSearchResultsInGraph()) {
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Create graph edges using search results");

    stk::balance::internal::addGraphEdgesUsingBBSearch(stkMeshBulkData, balanceSettings, graphEdges, selector);

    std::vector<unsigned> localIdsSearch(graphEdges.size(),0);
    size_t counter = 0;
    for (size_t i = 0; i < graphEdges.size(); ++i) {
      if (graphEdges[i].is_edge_from_search()) {
        localIdsSearch[counter] = stk::balance::internal::get_local_id(localIds, graphEdges[i].vertex1());
        counter++;
      }
    }

    localIdsSearch.resize(counter);
    stk::util::sort_and_unique(localIdsSearch);

    // Adjust vertex weights if part of search
    for (size_t i = 0; i < localIdsSearch.size(); ++i) {
      mVertexWeights[localIdsSearch[i]] *= balanceSettings.getVertexWeightMultiplierForVertexInSearch();
    }
  }

  stk::util::sort_and_unique(graphEdges);

  stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Convert edges to a graph");

  convertGraphEdgesToZoltanGraph(stkMeshBulkData,
                                 balanceSettings,
                                 graphEdges,
                                 numElements,
                                 adjacencyProcs,
                                 localIds);

  stk::balance::internal::logMessage(stkMeshBulkData.parallel(), "Finished creating graph");
}


