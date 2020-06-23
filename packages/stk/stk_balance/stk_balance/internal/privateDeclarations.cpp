#include <cmath>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/GeometricVertices.hpp>
#include <stk_balance/internal/StkGeometricMethodViaZoltan.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_balance/internal/SideGeometry.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/FieldBase.hpp>  // for field_data
#include <stk_mesh/base/GetEntities.hpp>  // for field_data
#include "stk_mesh/base/FEMHelpers.hpp"
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/LogWithTimeAndMemory.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <zoltan.h>
#include <Zoltan2_Version.hpp>
#include <map>

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_tools/mesh_tools/CustomAura.hpp"
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include <stk_util/environment/Env.hpp>

namespace stk {
namespace balance {

class DecompositionChangeList::Impl
{
public:
    Impl(stk::mesh::BulkData & stkMeshBulkData, const stk::mesh::EntityProcVec & decomposition)
      : m_stkMeshBulkData(stkMeshBulkData)
    {
        fill_change_list_from_raw_decomposition(decomposition);
    }

    bool has_entity(stk::mesh::Entity entity) const
    {
        return (m_dataMap.find(entity) != m_dataMap.end());
    }

    int get_entity_destination(stk::mesh::Entity entity) const
    {
        std::map<stk::mesh::Entity, int>::const_iterator entry = m_dataMap.find(entity);
        if (entry != m_dataMap.end()) {
            return entry->second;
        }
        return -1;
    }

    void set_entity_destination(stk::mesh::Entity entity, const int destination)
    {
        m_dataMap[entity] = destination;
    }

    void delete_entity(stk::mesh::Entity entity)
    {
        std::map<stk::mesh::Entity, int>::iterator entry = m_dataMap.find(entity);

        if (entry != m_dataMap.end()) {
            m_dataMap.erase(entry);
        }
    }

    stk::mesh::EntityProcVec get_all_partition_changes()
    {
        stk::mesh::EntityProcVec entityProcPairs;
        for (auto & entry : m_dataMap) {
            entityProcPairs.push_back(entry);
        }
        return entityProcPairs;
    }

    stk::mesh::BulkData & get_bulk() { return m_stkMeshBulkData; }

    void get_decomposition_with_full_closure(stk::mesh::EntityProcVec & finalDecomposition)
    {
        finalDecomposition.clear();
        for (auto & entry : m_dataMap) {
            finalDecomposition.push_back(entry);
            add_downward_closure_for_entity(entry, finalDecomposition);
        }
    }

    size_t get_num_global_entity_migrations() const
    {
        size_t num_local_entity_migrations = m_dataMap.size();
        size_t num_global_entity_migrations = 0;
        stk::all_reduce_sum(m_stkMeshBulkData.parallel(), &num_local_entity_migrations, &num_global_entity_migrations, 1);
        return num_global_entity_migrations;
    }

    size_t get_max_global_entity_migrations() const
    {
        size_t num_local_entity_migrations = m_dataMap.size();
        size_t max_global_entity_migrations = 0;
        stk::all_reduce_max(m_stkMeshBulkData.parallel(), &num_local_entity_migrations, &max_global_entity_migrations, 1);
        return max_global_entity_migrations;
    }

private:
    stk::mesh::BulkData &m_stkMeshBulkData;
    std::map<stk::mesh::Entity, int> m_dataMap;

    void fill_change_list_from_raw_decomposition(const stk::mesh::EntityProcVec& decomposition) {
        for(const stk::mesh::EntityProc& entity_proc : decomposition)
        {
            if(m_stkMeshBulkData.is_valid(entity_proc.first) && entity_proc.second != m_stkMeshBulkData.parallel_rank())
                m_dataMap[entity_proc.first] = entity_proc.second;
        }
    }

    void add_downward_closure_for_entity(const std::pair<const stk::mesh::Entity, int> & entityProc, stk::mesh::EntityProcVec & finalDecomposition)
    {
        const stk::topology::rank_t entityRank = m_stkMeshBulkData.entity_rank(entityProc.first);
        for (int rank = entityRank-1; rank >= stk::topology::NODE_RANK; --rank)
            internal::add_connected_entities_of_rank(m_stkMeshBulkData, entityProc.first, entityProc.second, static_cast<stk::topology::rank_t>(rank), finalDecomposition);
    }

};


DecompositionChangeList::DecompositionChangeList(stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::EntityProcVec& decomposition)
    : pImpl(new Impl(stkMeshBulkData, decomposition))
{ }
DecompositionChangeList::~DecompositionChangeList() { delete pImpl;  }

bool DecompositionChangeList::has_entity(stk::mesh::Entity entity)                                    { return pImpl->has_entity(entity); }
int  DecompositionChangeList::get_entity_destination(stk::mesh::Entity entity)                        { return pImpl->get_entity_destination(entity); }
void DecompositionChangeList::set_entity_destination(stk::mesh::Entity entity, const int destination) { pImpl->set_entity_destination(entity, destination); }
void DecompositionChangeList::delete_entity(stk::mesh::Entity entity)                                 { pImpl->delete_entity(entity); }
stk::mesh::BulkData & DecompositionChangeList::get_bulk()                                             { return pImpl->get_bulk(); }
stk::mesh::EntityProcVec DecompositionChangeList::get_all_partition_changes()                         { return pImpl->get_all_partition_changes(); }
size_t DecompositionChangeList::get_num_global_entity_migrations() const                              { return pImpl->get_num_global_entity_migrations(); }
size_t DecompositionChangeList::get_max_global_entity_migrations() const                              { return pImpl->get_max_global_entity_migrations(); }

namespace internal
{


unsigned get_local_id(const stk::mesh::impl::LocalIdMapper& localIds, stk::mesh::Entity entity)
{
    return localIds.entity_to_local(entity);
}

void addBoxForFace(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity face, const double eps, SearchBoxIdentProcs &faceBoxes, const stk::mesh::FieldBase* coord)
{

    unsigned numElements = stkMeshBulkData.num_elements(face);
    const stk::mesh::Entity *element = stkMeshBulkData.begin_elements(face);

    ThrowRequireWithSierraHelpMsg(numElements <= 1);
    if ( element != NULL && stkMeshBulkData.bucket(*element).owned() )
    {
        unsigned numNodes = stkMeshBulkData.num_nodes(face);
        const stk::mesh::Entity *nodes = stkMeshBulkData.begin_nodes(face);

        addBoxForNodes(stkMeshBulkData, numNodes, nodes, coord, eps, stkMeshBulkData.identifier(*element), faceBoxes);
    }
}

void addEdgeAndVertexWeightsForSearchResult(stk::mesh::BulkData& stkMeshBulkData,
                                            const BalanceSettings &balanceSettings,
                                            stk::mesh::EntityId element1Id,
                                            stk::mesh::EntityId element2Id,
                                            unsigned owningProcElement2,
                                            std::vector<GraphEdge>& graphEdges)
{
    stk::mesh::EntityKey entityKeyElement1(stk::topology::ELEMENT_RANK, element1Id);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(entityKeyElement1);
    ThrowRequireWithSierraHelpMsg(stkMeshBulkData.entity_rank(element1) == stk::topology::ELEMENT_RANK);

    if(stkMeshBulkData.is_valid(element1) && stkMeshBulkData.bucket(element1).owned() && element1Id != element2Id)
    {
        stk::mesh::EntityKey entityKeyElement2(stk::topology::ELEMENT_RANK, element2Id);
        stk::mesh::Entity element2 = stkMeshBulkData.get_entity(entityKeyElement2);

        unsigned anyIntersections = 0;
        if ( stkMeshBulkData.is_valid(element2) )
        {
            anyIntersections = stk::balance::internal::getNumSharedNodesBetweenElements(stkMeshBulkData, element1, element2);
        }

        bool elementIsNotConnectedViaNodes = anyIntersections == 0;
        if ( elementIsNotConnectedViaNodes )
        {
            double edge_weight = balanceSettings.getGraphEdgeWeightForSearch();
            graphEdges.push_back(GraphEdge(element1, element2Id, owningProcElement2, edge_weight, true));
        }
    }
}

void
addSearchResultsToGraphEdges(stk::mesh::BulkData & bulk,
                             const BalanceSettings & balanceSettings,
                             const SearchElemPairs & searchResults,
                             std::vector<GraphEdge> & graphEdges)
{
  for (const auto & searchResult : searchResults) {
    stk::mesh::EntityId element1id = searchResult.first.id();
    stk::mesh::EntityId element2id = searchResult.second.id();
    stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEM_RANK, element1id);
    int owningProcElement2 = searchResult.second.proc();

    double edge_weight = balanceSettings.getGraphEdgeWeightForSearch();
    graphEdges.push_back(GraphEdge(element1, element2id, owningProcElement2, edge_weight, true));
  }
}

void filterOutNonLocalResults(const stk::mesh::BulkData & bulk, SearchElemPairs & searchResults)
{
  const int myRank = bulk.parallel_rank();
  size_t numFiltered = 0;

  for (const auto & searchResult : searchResults) {
    if (searchResult.first.proc() == myRank) {
      searchResults[numFiltered] = searchResult;
      numFiltered++;
    }
  }

  searchResults.resize(numFiltered);
}

void filterOutConnectedElements(const stk::mesh::BulkData & bulk, SearchElemPairs & searchResults)
{
  size_t numFiltered = 0;

  for (const auto & searchResult : searchResults) {
    stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEM_RANK, searchResult.first.id());
    stk::mesh::Entity element2 = bulk.get_entity(stk::topology::ELEM_RANK, searchResult.second.id());
    int owningProcElement1 = searchResult.first.proc();
    int owningProcElement2 = searchResult.second.proc();

    ThrowRequireWithSierraHelpMsg(owningProcElement1 == bulk.parallel_rank() ||
                                  owningProcElement2 == bulk.parallel_rank());

    int numIntersections = 0;

    if (element1 == element2) {
      numIntersections = 1;
    }
    else if (bulk.is_valid(element1) && bulk.is_valid(element2)) {
      numIntersections = internal::getNumSharedNodesBetweenElements(bulk, element1, element2);
    }

    if (numIntersections == 0) {
      searchResults[numFiltered] = searchResult;
      numFiltered++;
    }
  }

  searchResults.resize(numFiltered);
}

stk::mesh::OrdinalVector getExposedSideOrdinals(const stk::mesh::BulkData & bulk, stk::mesh::EntityId elemId)
{
  stk::mesh::OrdinalVector exposedSideOrdinals;

  const stk::mesh::ElemElemGraph & graph = bulk.get_face_adjacent_element_graph();
  const stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
  const stk::topology elemTopology = bulk.bucket(element).topology();

  for (stk::mesh::Ordinal sideOrd = 0; sideOrd < elemTopology.num_sides(); ++sideOrd) {
    if (!graph.is_connected_to_other_element_via_side_ordinal(element, sideOrd)) {
      exposedSideOrdinals.push_back(sideOrd);
    }
  }

  return exposedSideOrdinals;
}

struct SideInfo {
  SideInfo()
    : sideTopology(stk::topology::INVALID_TOPOLOGY),
      sideSearchTolerance(0.0)
  {}

  SideInfo(stk::topology::topology_t _sideTopology,
           double _sideSearchTolerance,
           const std::vector<stk::math::Vector3d> & _nodeCoordinates)
    : sideTopology(_sideTopology),
      sideSearchTolerance(_sideSearchTolerance),
      nodeCoordinates(_nodeCoordinates)
  {}

  ~SideInfo() = default;

  stk::topology::topology_t sideTopology;
  double sideSearchTolerance;
  std::vector<stk::math::Vector3d> nodeCoordinates;
};

using SideInfoMap = std::map<stk::mesh::EntityId, std::vector<SideInfo>>;
using SideGeometryPtr = std::unique_ptr<SideGeometry>;

bool is_line_side(const stk::topology::topology_t & t)
{
  return (t == stk::topology::LINE_2) || (t == stk::topology::LINE_3);
}

bool is_point_side(const stk::topology::topology_t & t)
{
  return t == stk::topology::NODE;
}

SideGeometryPtr makeSideGeometry(const SideInfo & sideInfo)
{
  if (stk::is_quad_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new QuadGeometry(sideInfo.nodeCoordinates[0],
                                              sideInfo.nodeCoordinates[1],
                                              sideInfo.nodeCoordinates[2],
                                              sideInfo.nodeCoordinates[3]));
  }
  else if (stk::is_tri_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new TriGeometry(sideInfo.nodeCoordinates[0],
                                             sideInfo.nodeCoordinates[1],
                                             sideInfo.nodeCoordinates[2]));
  }
  else if (is_line_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new LineGeometry(sideInfo.nodeCoordinates[0],
                                            sideInfo.nodeCoordinates[1]));
  }
  else if (is_point_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new PointGeometry(sideInfo.nodeCoordinates[0]));
  }
  else {
    ThrowErrorMsg("Unsupported side topology: " << stk::topology(sideInfo.sideTopology).name());
    return SideGeometryPtr();
  }
}

bool isAdjacent(const SideInfoMap & sideInfoMap,
                const SearchElemPair & elemPair)
{
  stk::mesh::EntityId element1Id = elemPair.first.id();
  stk::mesh::EntityId element2Id = elemPair.second.id();
  const std::vector<SideInfo> & sides1 = sideInfoMap.at(element1Id);
  const std::vector<SideInfo> & sides2 = sideInfoMap.at(element2Id);

  for (const SideInfo & sideInfo1 : sides1) {
    const SideGeometryPtr side1 = makeSideGeometry(sideInfo1);
    const double tol1 = sideInfo1.sideSearchTolerance;

    for (const SideInfo & sideInfo2 : sides2) {
      const SideGeometryPtr side2 = makeSideGeometry(sideInfo2);
      const double tol2 = sideInfo2.sideSearchTolerance;

      if (side1->are_nodes_close_to_side(*side2, tol1) || side2->are_nodes_close_to_side(*side1, tol2)) {
        return true;
      }
    }
  }

  return false;
}

std::vector<stk::math::Vector3d> nodeCoordinates(const stk::mesh::EntityVector& sideNodes,
                                                        const stk::mesh::FieldBase& coords,
                                                        int spatialDim)
{
    std::vector<stk::math::Vector3d> nodeCoordinates;
    for (const stk::mesh::Entity & node : sideNodes) {
        if (spatialDim == 3) {
            nodeCoordinates.emplace_back(reinterpret_cast<double*>(stk::mesh::field_data(coords, node)));
        }
        else if (spatialDim == 2) {
            double* nodeFieldData = reinterpret_cast<double*>(stk::mesh::field_data(coords, node));
            nodeCoordinates.emplace_back(nodeFieldData[0], nodeFieldData[1], 0.0);
        }
        else {
            ThrowErrorMsg("Problem dimensionality " << spatialDim << " not supported");
        }
    }
    return nodeCoordinates;
}

std::vector<SideInfo> getElementExposedFaceInfo(const stk::mesh::BulkData & bulk,
                                                const BalanceSettings & balanceSettings,
                                                stk::mesh::EntityId elementId)
{
  stk::mesh::OrdinalVector sideOrdinals = getExposedSideOrdinals(bulk, elementId);
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elementId);
  const stk::mesh::FieldBase & coords = *bulk.mesh_meta_data().coordinate_field();
  int spatialDim = bulk.mesh_meta_data().spatial_dimension();
  const stk::topology elemTopology = bulk.bucket(element).topology();

  stk::mesh::EntityVector sideNodes;
  std::vector<SideInfo> sideInfoVec;

  for (unsigned ord : sideOrdinals) {
    const stk::topology sideTopology = elemTopology.side_topology(ord);
    stk::mesh::get_subcell_nodes(bulk, element, bulk.mesh_meta_data().side_rank(), ord, sideNodes);
    const double tol = balanceSettings.getToleranceForFaceSearch(bulk, coords,
                                                                 sideNodes.data(), sideNodes.size());

    sideInfoVec.emplace_back(sideTopology, tol, nodeCoordinates(sideNodes, coords, spatialDim));
  }

  return sideInfoVec;
}

std::vector<SideInfo> getParticleSideInfo(const stk::mesh::BulkData & bulk,
                                          const BalanceSettings & balanceSettings,
                                          stk::mesh::EntityId elementId)
{
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elementId);
  const stk::mesh::FieldBase & coords = *bulk.mesh_meta_data().coordinate_field();
  int spatialDim = bulk.mesh_meta_data().spatial_dimension();

  const stk::topology sideTopology = stk::topology::NODE;
  const double tol = balanceSettings.getAbsoluteToleranceForParticleSearch(element);
  const stk::mesh::Entity node = bulk.begin_nodes(element)[0];

  std::vector<SideInfo> sideInfoVec;
  sideInfoVec.emplace_back(sideTopology, tol, nodeCoordinates({node}, coords, spatialDim));

  return sideInfoVec;
}

std::vector<SideInfo> getElementExposedSideInfo(const stk::mesh::BulkData & bulk,
                                                const BalanceSettings & balanceSettings,
                                                stk::mesh::EntityId elementId)
{
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elementId);
  const stk::topology elemTopology = bulk.bucket(element).topology();

  if (elemTopology == stk::topology::PARTICLE) {
    return getParticleSideInfo(bulk, balanceSettings, elementId);
  }
  else {
    return getElementExposedFaceInfo(bulk, balanceSettings, elementId);
  }
}

void insertSideInfoIfNew(const stk::mesh::BulkData & bulk,
                         const BalanceSettings & balanceSettings,
                         SideInfoMap & sideInfoMap,
                         stk::mesh::EntityId elemId)
{
  if (sideInfoMap.count(elemId) == 0) {
    sideInfoMap.emplace(std::make_pair(elemId, getElementExposedSideInfo(bulk, balanceSettings, elemId)));
  }
}

void insertLocalSideInfo(const stk::mesh::BulkData & bulk,
                         const BalanceSettings & balanceSettings,
                         const SearchElemPairs & searchResults,
                         SideInfoMap & sideInfoMap)
{
  for (const SearchElemPair & elemPair : searchResults) {
    insertSideInfoIfNew(bulk, balanceSettings, sideInfoMap, elemPair.first.id());
    if (elemPair.second.proc() == bulk.parallel_rank()) {
      insertSideInfoIfNew(bulk, balanceSettings, sideInfoMap, elemPair.second.id());
    }
  }
}

void packNonLocalSideInfo(const stk::mesh::BulkData & bulk,
                          stk::CommSparse & comm,
                          const SearchElemPairs & searchResults,
                          const SideInfoMap & sideInfoMap)
{
  for (const SearchElemPair & elemPair : searchResults) {
    const int secondProc = elemPair.second.proc();
    if (secondProc != bulk.parallel_rank()) {
      const std::vector<SideInfo> & firstSideInfo = sideInfoMap.at(elemPair.first.id());

      comm.send_buffer(secondProc).pack<stk::mesh::EntityId>(elemPair.first.id());
      comm.send_buffer(secondProc).pack<unsigned>(firstSideInfo.size());
      for (const SideInfo & sideInfo : firstSideInfo) {
        comm.send_buffer(secondProc).pack<stk::topology::topology_t>(sideInfo.sideTopology);
        comm.send_buffer(secondProc).pack<double>(sideInfo.sideSearchTolerance);
        comm.send_buffer(secondProc).pack<unsigned>(sideInfo.nodeCoordinates.size());
        for (const stk::math::Vector3d & nodeCoords : sideInfo.nodeCoordinates) {
          comm.send_buffer(secondProc).pack<stk::math::Vector3d>(nodeCoords);
        }
      }
    }
  }
}

void unpackAndInsertNonLocalSideInfo(const stk::mesh::BulkData & bulk,
                                     stk::CommSparse & comm,
                                     SideInfoMap & sideInfoMap)
{
  for (int proc = 0; proc < bulk.parallel_size(); ++proc) {
    while (comm.recv_buffer(proc).remaining()) {
      stk::mesh::EntityId id;
      unsigned numSides;
      comm.recv_buffer(proc).unpack(id);
      comm.recv_buffer(proc).unpack(numSides);
      std::vector<SideInfo> elementSideInfo(numSides);

      for (SideInfo & si : elementSideInfo) {
        unsigned numSideNodes;
        comm.recv_buffer(proc).unpack(si.sideTopology);
        comm.recv_buffer(proc).unpack(si.sideSearchTolerance);
        comm.recv_buffer(proc).unpack(numSideNodes);
        si.nodeCoordinates.resize(numSideNodes);
        for (stk::math::Vector3d & nodeCoords : si.nodeCoordinates) {
          comm.recv_buffer(proc).unpack(nodeCoords);
        }
      }
      sideInfoMap[id] = elementSideInfo;
    }
  }
}

SideInfoMap fillSideInfo(const stk::mesh::BulkData & bulk,
                         const BalanceSettings & balanceSettings,
                         const SearchElemPairs & searchResults)
{
  SideInfoMap sideInfoMap;

  insertLocalSideInfo(bulk, balanceSettings, searchResults, sideInfoMap);

  stk::CommSparse comm(bulk.parallel());

  stk::pack_and_communicate(comm, [&](){
    packNonLocalSideInfo(bulk, comm, searchResults, sideInfoMap);
  });

  unpackAndInsertNonLocalSideInfo(bulk, comm, sideInfoMap);

  return sideInfoMap;
}

void filterOutNonAdjacentElements(const stk::mesh::BulkData & bulk,
                                  const BalanceSettings & balanceSettings,
                                  SearchElemPairs & searchResults)
{
  SideInfoMap sideInfoMap = fillSideInfo(bulk, balanceSettings, searchResults);

  size_t numFiltered = 0;

  for (const SearchElemPair & elemPair : searchResults) {
    if (isAdjacent(sideInfoMap, elemPair)) {
      searchResults[numFiltered] = elemPair;
      numFiltered++;
    }
  }

  searchResults.resize(numFiltered);
}

void printGraphEdgeCounts(const stk::mesh::BulkData& stkMeshBulkData,
                          size_t edgeCounts,
                          const std::string& message)
{
  std::ostringstream os;
  size_t max = 0, min = 0, avg = 0;
  stk::get_max_min_avg(stkMeshBulkData.parallel(), edgeCounts, max, min, avg);

  os << message << ", have following distribution of graph edges: min="
     << min << ", avg=" << avg << ", max=" << max;

  logMessage(stkMeshBulkData.parallel(), os.str());
}

void addGraphEdgesUsingBBSearch(stk::mesh::BulkData & stkMeshBulkData,
                                const BalanceSettings & balanceSettings,
                                std::vector<GraphEdge> & graphEdges,
                                const stk::mesh::Selector & searchSelector)
{
  printGraphEdgeCounts(stkMeshBulkData, graphEdges.size(), "Starting search");

  SearchElemPairs searchResults = getBBIntersectionsForFacesParticles(stkMeshBulkData, balanceSettings, searchSelector);
  printGraphEdgeCounts(stkMeshBulkData, searchResults.size(), "Finished search");

  filterOutNonLocalResults(stkMeshBulkData, searchResults);
  filterOutConnectedElements(stkMeshBulkData, searchResults);
  filterOutNonAdjacentElements(stkMeshBulkData, balanceSettings, searchResults);

  addSearchResultsToGraphEdges(stkMeshBulkData, balanceSettings, searchResults, graphEdges);
  printGraphEdgeCounts(stkMeshBulkData, graphEdges.size(), "After search");
}

std::vector<int> getLocalIdsOfEntitiesNotSelected(const stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Selector selector, const stk::mesh::impl::LocalIdMapper& localIds)
{
    selector = (!selector) & stkMeshBulkData.mesh_meta_data().locally_owned_part();
    std::vector<int> local_ids;
    size_t num_total_elements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

    const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, selector);
    for (size_t i=0; i<buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        for (size_t j=0;j<bucket.size();j++)
        {
            unsigned local_id = get_local_id(localIds, bucket[j]);
            ThrowRequireWithSierraHelpMsg(local_id < num_total_elements);
            local_ids.push_back(local_id);
        }
    }
    return local_ids;
}

void fillEntityCentroid(const stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::FieldBase* coord, stk::mesh::Entity entityOfConcern, double *entityCentroid)
{
    unsigned spatialDimension = stkMeshBulkData.mesh_meta_data().spatial_dimension();
    if(stkMeshBulkData.entity_rank(entityOfConcern)==stk::topology::NODE_RANK)
    {
        double *nodeCoord = static_cast<double *>(stk::mesh::field_data(*coord, entityOfConcern));
        for(unsigned k=0; k < spatialDimension; k++)
        {
            entityCentroid[k] = nodeCoord[k];
        }
    }
    else
    {
        stk::mesh::Entity const * nodes = stkMeshBulkData.begin_nodes(entityOfConcern);
        unsigned numNodes = stkMeshBulkData.num_nodes(entityOfConcern);
        for(unsigned nodeIndex=0; nodeIndex < numNodes; nodeIndex++)
        {
            double *nodeCoord = static_cast<double *>(stk::mesh::field_data(*coord, nodes[nodeIndex]));
            for(unsigned k=0; k < spatialDimension; k++)
            {
                entityCentroid[k] += nodeCoord[k];
            }
        }
        for(unsigned k=0; k < spatialDimension; k++)
        {
            entityCentroid[k] /= numNodes;
        }
    }
}

int num_beams_connected_to_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity node)
{
  const unsigned numElems = bulk.num_elements(node);
  const stk::mesh::Entity* elems = bulk.begin_elements(node);

  int numBeams = 0;
  for(unsigned elemIndex=0; elemIndex<numElems; ++elemIndex) {
    if (bulk.bucket(elems[elemIndex]).topology() == stk::topology::BEAM_2) {
      ++numBeams;
    }
  }

  return numBeams;
}

bool is_not_part_of_spider(stk::topology::topology_t elemTopology)
{
  return ((elemTopology != stk::topology::PARTICLE) &&
          (elemTopology != stk::topology::BEAM_2));
}

int num_volume_elements_connected_to_beam(const stk::mesh::BulkData& bulk, stk::mesh::Entity beam)
{
  const unsigned numNodes = bulk.num_nodes(beam);
  const stk::mesh::Entity* nodes = bulk.begin_nodes(beam);

  int numVolumeElems = 0;
  for (unsigned nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex) {
    const unsigned numElems = bulk.num_elements(nodes[nodeIndex]);
    const stk::mesh::Entity* elems = bulk.begin_elements(nodes[nodeIndex]);
    for (unsigned elemIndex = 0; elemIndex < numElems; ++elemIndex) {
      const stk::topology::topology_t elemTopology = bulk.bucket(elems[elemIndex]).topology();
      if (is_not_part_of_spider(elemTopology)) {
        ++numVolumeElems;
      }
    }
  }

  return numVolumeElems;
}

void fill_spider_connectivity_count_fields(stk::mesh::BulkData & bulk, const BalanceSettings & balanceSettings)
{
  if (balanceSettings.shouldFixSpiders()) {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

    const stk::mesh::Field<int> * beamConnectivityCountField = balanceSettings.getSpiderBeamConnectivityCountField(bulk);
    stk::mesh::Selector selectBeamNodes(meta.locally_owned_part() &
                                        meta.get_topology_root_part(stk::topology::BEAM_2));
    stk::mesh::EntityVector nodes;
    mesh::get_selected_entities(selectBeamNodes, bulk.buckets(stk::topology::NODE_RANK), nodes);
    for(stk::mesh::Entity node : nodes) {
      int * beamConnectivityCount = stk::mesh::field_data(*beamConnectivityCountField, node);
      *beamConnectivityCount = num_beams_connected_to_node(bulk, node);
    }

    const stk::mesh::Field<int> * volumeConnectivityCountField = balanceSettings.getSpiderVolumeConnectivityCountField(bulk);
    stk::mesh::Selector selectBeamElements(meta.locally_owned_part() &
                                           meta.get_topology_root_part(stk::topology::BEAM_2));
    stk::mesh::EntityVector beams;
    mesh::get_selected_entities(selectBeamElements, bulk.buckets(stk::topology::ELEM_RANK), beams);
    for(stk::mesh::Entity beam : beams) {
      int * volumeConnectivityCount = stk::mesh::field_data(*volumeConnectivityCountField, beam);
      *volumeConnectivityCount = num_volume_elements_connected_to_beam(bulk, beam);
    }

    stk::mesh::communicate_field_data(bulk, {beamConnectivityCountField, volumeConnectivityCountField});
  }
}

void logMessage(MPI_Comm communicator, const std::string &message)
{
    stk::log_with_time_and_memory(communicator, message);
}



void fill_zoltan2_graph(const BalanceSettings& balanceSettings,
                        stk::mesh::BulkData& stkMeshBulkData,
                        Zoltan2ParallelGraph& zoltan2Graph,
                        const stk::mesh::Selector& searchSelector,
                        const stk::mesh::impl::LocalIdMapper& localIds)
{
    std::vector<int> adjacencyProcs;
    zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                   balanceSettings,
                                                   adjacencyProcs,
                                                   searchSelector,
                                                   localIds);

    logMessage(stkMeshBulkData.parallel(), "Finished filling in graph data");
}


//#define WRITE_OUT_DEBUGGING_INFO
//#define WRITE_OUT_DECOMP_METRICS#include "StkGeometricMethodViaZoltan.hpp"


stk::mesh::EntityVector get_entities_to_balance(stk::mesh::Selector selector, stk::mesh::EntityRank primaryRank, const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector entitiesToBalance;
    selector = selector & bulkData.mesh_meta_data().locally_owned_part();
    stk::mesh::get_selected_entities(selector, bulkData.buckets(primaryRank), entitiesToBalance);
    return entitiesToBalance;
}


void get_multicriteria_decomp_using_selectors_as_segregation(const stk::mesh::BulkData& mesh, const std::vector<stk::mesh::Selector>& criterions, const BalanceSettings& balanceSettings,
                                                             const int numSubdomainsToCreate, stk::mesh::EntityProcVec& decomp, const stk::mesh::impl::LocalIdMapper& localIds)
{
    stk::mesh::Selector unionSelector = stk::mesh::selectUnion(criterions);
    stk::mesh::EntityVector entitiesToBalance = get_entities_to_balance(unionSelector, stk::topology::ELEM_RANK, mesh);
    size_t num_entities = entitiesToBalance.size();
    size_t num_entities_across_procs = 0;
    stk::all_reduce_sum(mesh.parallel(), &num_entities, &num_entities_across_procs, 1);
    if(num_entities_across_procs > 0)
    {
        stk::balance::internal::GeometricVertices vertexInfo(balanceSettings, mesh, entitiesToBalance, criterions);
        std::vector<unsigned> processorOntoWhichEntityBelongs = stk::balance::get_decomposition(vertexInfo, balanceSettings, numSubdomainsToCreate, mesh.parallel());
        for(size_t i=0;i<entitiesToBalance.size();++i)
        {
            int local_id = get_local_id(localIds, entitiesToBalance[i]);
            decomp[local_id] = std::make_pair(entitiesToBalance[i], processorOntoWhichEntityBelongs[i]);
        }
    }
}


void fill_decomp_using_geometric_method(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate, stk::mesh::EntityProcVec &decomp,
                                        stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds)
{
    logMessage(stkMeshBulkData.parallel(), "Using Zoltan2 version: " + Zoltan2::Zoltan2_Version());
    logMessage(stkMeshBulkData.parallel(), "Filling in vertex data for decomp method = " + balanceSettings.getDecompMethod());

    size_t numEntities = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

    decomp.clear();
    decomp.resize(numEntities, std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));


    if (balanceSettings.isMultiCriteriaRebalance())
        get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, selectors, balanceSettings, numSubdomainsToCreate, decomp, localIds);
    else
        for(const stk::mesh::Selector selector : selectors)
            get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, std::vector<stk::mesh::Selector>{selector}, balanceSettings, numSubdomainsToCreate, decomp, localIds);

    logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");
}

void get_multicriteria_parmetis_decomp(const stk::mesh::BulkData &mesh, const BalanceSettings& balanceSettings, Zoltan2ParallelGraph &zoltan2Graph, Teuchos::ParameterList &params,
                                       stk::mesh::Selector selector, stk::mesh::EntityProcVec &decomp, const stk::mesh::impl::LocalIdMapper& localIds)
{
    StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

    logMessage(mesh.parallel(), "Setting up partitioning problem");

    Zoltan2::PartitioningProblem<StkMeshZoltanAdapter> problem(&stkMeshAdapter, &params, mesh.parallel());

    logMessage(mesh.parallel(), "Solving");

    if(balanceSettings.shouldPrintMetrics())
    {
        internal::print_statistics(stkMeshAdapter, mesh.parallel(), mesh.parallel_rank());
    }

    std::srand(mesh.parallel_rank()); // KHP: Temporary until an API is added to Zoltan2 for random seeds.
    problem.solve();

    if(balanceSettings.shouldPrintMetrics())
        internal::print_solution_statistics(stkMeshAdapter, problem.getSolution(), mesh.parallel(), mesh.parallel_rank());

#if defined(WRITE_OUT_DEBUGGING_INFO)
    std::vector<int> local_ids_of_elements_to_balance;
    std::set_difference(all_local_ids.begin(), all_local_ids.end(), localIdsOfNotSelectedEntities.begin(), localIdsOfNotSelectedEntities.end(),
                            std::inserter(local_ids_of_elements_to_balance, local_ids_of_elements_to_balance.begin()));

    std::ostringstream filename;
    filename << "rebalance_proc_" << mesh.parallel_rank() << "_for_selector_" << i << "_for_step_" << step << ".txt";
    std::ofstream out(filename.str().c_str());
    out << "For this selector: " << selectors[i] << " the following entities are being balanced: ";
    for(size_t j=0;j<local_ids_of_elements_to_balance.size();++j)
    {
        stk::mesh::Entity element = entities[local_ids_of_elements_to_balance[j]];
        out << "Element " << j << " has: " << mesh.entity_key(element) << " with topology " << mesh.bucket(element).topology() << std::endl;
    }
    stkMeshAdapter.debuggingInfo(mesh.parallel_rank(), out);
    out.close();
#endif

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(mesh.mesh_meta_data().locally_owned_part() & selector, mesh.buckets(stk::topology::ELEM_RANK), elements);

    // local entity j is on which processor
    const StkMeshZoltanAdapter::part_t *processorOntoWhichEntityBelongs = problem.getSolution().getPartListView();
    for(size_t j = 0; j < elements.size(); ++j)
    {
        int local_id = get_local_id(localIds, elements[j]);
        int dest_proc = processorOntoWhichEntityBelongs[local_id];
        decomp[local_id] = std::make_pair(elements[j], dest_proc);
    }
}

Teuchos::ParameterList getGraphBasedParameters(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate)
{
    Teuchos::ParameterList params("test params");
    params.set("debug_level", "no_status");
//    params.set("debug_level", "basic_status");

    int nparts = numSubdomainsToCreate;
    double imbalance_allowed = balanceSettings.getImbalanceTolerance();
    params.set("imbalance_tolerance", imbalance_allowed);
    params.set("num_global_parts", nparts);
    params.set("algorithm", balanceSettings.getDecompMethod());

//    params.set("partitioning_objective", "balance_object_weight");
//    params.set("partitioning_objective", "multicriteria_minimize_total_weight");
//    params.set("partitioning_objective", "multicriteria_minimize_maximum_weight");
//    params.set("partitioning_objective", "multicriteria_balance_total_maximum");

    if (balanceSettings.isIncrementalRebalance())
    {
        params.set("partitioning_approach", "repartition");
        params.set("remap_parts", true);
    }

    // should not hurt other methods, only affects RCB.
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
    zparams.set("debug_level","0");
    zparams.set("LB_APPROACH","PARTITION");
    zparams.set("LB_METHOD","GRAPH");
    zparams.set("GRAPH_PACKAGE","ParMETIS");
    zparams.set("GRAPH_SYMMETRIZE","None");
    zparams.set("PARMETIS_METHOD","PartKway");
    //zparams.set("PARMETIS_METHOD","AdaptiveRepart");
    //zparams.set("PARMETIS_ITR",1000);
    return params;
}

bool isNodePartOfSpider(const stk::mesh::BulkData& stkMeshBulkData,
                        const stk::mesh::Field<int>& beamConnectivityCountField,
                        stk::mesh::Entity node)
{
    const int spiderConnectivityThreshold = 5;
    const int connectivityCount = *stk::mesh::field_data(beamConnectivityCountField, node);
    return (connectivityCount > spiderConnectivityThreshold);
}

bool isElementPartOfSpider(const stk::mesh::BulkData& stkMeshBulkData,
                           const stk::mesh::Field<int>& beamConnectivityCountField,
                           stk::mesh::Entity element)
{
    const stk::mesh::Entity* nodes = stkMeshBulkData.begin_nodes(element);
    const unsigned numNodes = stkMeshBulkData.num_nodes(element);
    for (unsigned i = 0; i < numNodes; ++i) {
      if (isNodePartOfSpider(stkMeshBulkData, beamConnectivityCountField, nodes[i])) {
        return true;
      }
    }
    return false;
}

bool shouldOmitSpiderElement(const stk::mesh::BulkData & stkMeshBulkData,
                             const stk::balance::BalanceSettings & balanceSettings,
                             stk::mesh::Entity elem)
{
    bool omitConnection = false;
    if (balanceSettings.shouldFixSpiders()) {
        const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderBeamConnectivityCountField(stkMeshBulkData);
        stk::topology elemTopology = stkMeshBulkData.bucket(elem).topology();

        if (elemTopology == stk::topology::BEAM_2 || elemTopology == stk::topology::PARTICLE) {
            omitConnection = isElementPartOfSpider(stkMeshBulkData, beamConnectivityCountField, elem);
        }
    }

    return omitConnection;
}

bool found_valid_new_owner(const stk::mesh::BulkData & bulk, int newOwner)
{
  return (newOwner < bulk.parallel_size());
}

bool spider_body_element_exists(const stk::mesh::BulkData & bulk, const stk::mesh::Entity & bodyElement)
{
  return bulk.is_valid(bodyElement);
}

stk::mesh::Entity get_spider_particle_body_for_leg(const stk::mesh::BulkData & bulk, const stk::mesh::Entity & bodyNode)
{
  const stk::mesh::Entity * elems = bulk.begin_elements(bodyNode);
  const unsigned numElements = bulk.num_elements(bodyNode);
  for (unsigned i = 0; i < numElements; ++i) {
    const bool isParticle = (bulk.bucket(elems[i]).topology() == stk::topology::PARTICLE);
    if (isParticle) {
      return elems[i];
    }
  }
  return stk::mesh::Entity();
}

std::pair<stk::mesh::Entity, stk::mesh::Entity> get_spider_beam_body_and_node_for_leg(const stk::mesh::BulkData & bulk,
                                                                                      const BalanceSettings & balanceSettings,
                                                                                      const stk::mesh::Entity & bodyNode)
{
  const stk::mesh::Entity * elems = bulk.begin_elements(bodyNode);
  const unsigned numElements = bulk.num_elements(bodyNode);
  const stk::mesh::Field<int> & volumeElemConnectivityCountField = *balanceSettings.getSpiderVolumeConnectivityCountField(bulk);

  for (unsigned i = 0; i < numElements; ++i) {
    if (bulk.bucket(elems[i]).topology() == stk::topology::BEAM_2) {
      const int volumeElemConnectivityCount = *stk::mesh::field_data(volumeElemConnectivityCountField, elems[i]);
      if (volumeElemConnectivityCount == 0) {
        const stk::mesh::Entity * beamNodes = bulk.begin_nodes(elems[i]);
        const stk::mesh::Entity nonBodyNode = (beamNodes[0] == bodyNode) ? beamNodes[1] : beamNodes[0];
        return std::make_pair(elems[i], nonBodyNode);
      }
    }
  }

  return std::make_pair(stk::mesh::Entity(), stk::mesh::Entity());
}

void update_new_spider_entity_owner(stk::mesh::EntityProcMap & newSpiderEntityOwners,
                                    stk::mesh::Entity entity, int candidateOwner)
{
  if (newSpiderEntityOwners.find(entity) == newSpiderEntityOwners.end()) {
    newSpiderEntityOwners[entity] = candidateOwner;
  }
  newSpiderEntityOwners[entity] = std::min(candidateOwner, newSpiderEntityOwners[entity]);
}

void pack_new_spider_entity_owners(const stk::mesh::BulkData & bulk,
                                   stk::CommSparse & comm,
                                   const stk::mesh::EntityProcMap & newSpiderEntityOwners)
{
  for (const stk::mesh::EntityProc & entityNewOwner : newSpiderEntityOwners) {
    const int currentOwner = bulk.parallel_owner_rank(entityNewOwner.first);
    const stk::mesh::EntityKey entityKey = bulk.entity_key(entityNewOwner.first);
    comm.send_buffer(currentOwner).pack<stk::mesh::EntityKey>(entityKey);
    comm.send_buffer(currentOwner).pack<int>(entityNewOwner.second);
  }
}

stk::mesh::EntityProcMap unpack_new_spider_entity_owners(const stk::mesh::BulkData & bulk,
                                                         stk::CommSparse & comm)
{
  stk::mesh::EntityProcMap newSpiderEntityOwners;
  for (int proc = 0; proc < bulk.parallel_size(); ++proc) {
    while (comm.recv_buffer(proc).remaining()) {
      stk::mesh::EntityKey entityKey;
      int newOwner;
      comm.recv_buffer(proc).unpack(entityKey);
      comm.recv_buffer(proc).unpack(newOwner);
      update_new_spider_entity_owner(newSpiderEntityOwners, bulk.get_entity(entityKey), newOwner);
    }
  }
  return newSpiderEntityOwners;
}


stk::mesh::EntityProcMap determine_global_new_owner(const stk::mesh::BulkData & bulk,
                                                    const stk::mesh::EntityProcMap & newSpiderEntityOwners)
{
  stk::CommSparse comm(bulk.parallel());

  stk::pack_and_communicate(comm, [&](){
    pack_new_spider_entity_owners(bulk, comm, newSpiderEntityOwners);
  });

  return unpack_new_spider_entity_owners(bulk, comm);
}

void fix_spider_elements(const BalanceSettings & balanceSettings, stk::mesh::BulkData & bulk)
{
  stk::mesh::Ghosting * customAura = stk::tools::create_custom_aura(bulk, bulk.mesh_meta_data().globally_shared_part(), "customAura");

  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderBeamConnectivityCountField(bulk);

  stk::mesh::EntityVector beams;
  stk::mesh::Part & beamPart = meta.get_topology_root_part(stk::topology::BEAM_2);
  stk::mesh::get_selected_entities(beamPart & meta.locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), beams);

  stk::mesh::EntityVector particles;
  stk::mesh::EntityProcMap newSpiderEntityOwners;

  stk::mesh::EntityProcVec entitiesToMove;
  for (stk::mesh::Entity spiderLeg : beams) {
    if (isElementPartOfSpider(bulk, beamConnectivityCountField, spiderLeg)) {
      const stk::mesh::Entity* nodes = bulk.begin_nodes(spiderLeg);
      const int node1ConnectivityCount = *stk::mesh::field_data(beamConnectivityCountField, nodes[0]);
      const int node2ConnectivityCount = *stk::mesh::field_data(beamConnectivityCountField, nodes[1]);

      const stk::mesh::Entity footNode = (node1ConnectivityCount < node2ConnectivityCount) ? nodes[0] : nodes[1];
      const stk::mesh::Entity bodyNode = (node1ConnectivityCount < node2ConnectivityCount) ? nodes[1] : nodes[0];
      const stk::mesh::Entity* elements = bulk.begin_elements(footNode);
      const unsigned numElements = bulk.num_elements(footNode);
      int newLegOwner = std::numeric_limits<int>::max();

      for (unsigned i = 0; i < numElements; ++i) {
        const stk::topology::topology_t elemTopology = bulk.bucket(elements[i]).topology();
        if (is_not_part_of_spider(elemTopology)) {
          newLegOwner = std::min(newLegOwner, bulk.parallel_owner_rank(elements[i]));
        }
      }

      if (found_valid_new_owner(bulk, newLegOwner)) {
        if (newLegOwner != bulk.parallel_rank()) {
          entitiesToMove.push_back(std::make_pair(spiderLeg, newLegOwner));
        }

        update_new_spider_entity_owner(newSpiderEntityOwners, bodyNode, newLegOwner);

          const stk::mesh::Entity spiderParticleBody = get_spider_particle_body_for_leg(bulk, bodyNode);
          if (spider_body_element_exists(bulk, spiderParticleBody)) {
            update_new_spider_entity_owner(newSpiderEntityOwners, spiderParticleBody, newLegOwner);
          }
          else {
            const std::pair<stk::mesh::Entity, stk::mesh::Entity> spiderBeamBodyAndNode =
                get_spider_beam_body_and_node_for_leg(bulk, balanceSettings, bodyNode);
            if (spider_body_element_exists(bulk, spiderBeamBodyAndNode.first)) {
              update_new_spider_entity_owner(newSpiderEntityOwners, spiderBeamBodyAndNode.first, newLegOwner);
              update_new_spider_entity_owner(newSpiderEntityOwners, spiderBeamBodyAndNode.second, newLegOwner);
            }

        }

        if (newLegOwner != bulk.parallel_owner_rank(footNode)) {
          update_new_spider_entity_owner(newSpiderEntityOwners, footNode, newLegOwner);
        }
      }
    }
  }

  stk::mesh::EntityProcMap globalNewSpiderEntityOwners = determine_global_new_owner(bulk, newSpiderEntityOwners);

  for (const auto & entityOwner : globalNewSpiderEntityOwners) {
    if (entityOwner.second != bulk.parallel_rank()) {
      entitiesToMove.push_back(entityOwner);
    }
  }

  stk::tools::destroy_custom_aura(bulk, customAura);
  bulk.change_entity_owner(entitiesToMove);
}

void keep_spiders_on_original_proc(stk::mesh::BulkData &bulk, const stk::balance::BalanceSettings & balanceSettings, DecompositionChangeList &changeList)
{
    // Need to keep spiders on the original proc until the remaining elements have moved,
    // so that we can properly determine the final ownership of the elements on the end.
    // Then, we can move them.
    //
    const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderBeamConnectivityCountField(bulk);

    stk::mesh::EntityProcVec entityProcs = changeList.get_all_partition_changes();
    for (const stk::mesh::EntityProc & entityProc : entityProcs) {
        stk::mesh::Entity entity = entityProc.first;
        const stk::topology entityTopology = bulk.bucket(entity).topology();
        if (entityTopology == stk::topology::BEAM_2 || entityTopology == stk::topology::PARTICLE) {
            if (isElementPartOfSpider(bulk, beamConnectivityCountField, entity)) {
                changeList.delete_entity(entity);
            }
        }
    }
}

void createZoltanParallelGraph(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData,
                               const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds,
                               Zoltan2ParallelGraph& zoltan2Graph)
{
    std::vector<size_t> counts;
    stk::mesh::Selector locallyOwnedSelector(stkMeshBulkData.mesh_meta_data().locally_owned_part());

    stk::mesh::comm_mesh_counts(stkMeshBulkData, counts, &locallyOwnedSelector);
    zoltan2Graph.set_num_global_elements(counts[stk::topology::ELEM_RANK]);
    zoltan2Graph.set_spatial_dim(stkMeshBulkData.mesh_meta_data().spatial_dimension());
    if (balanceSettings.isMultiCriteriaRebalance())
        zoltan2Graph.set_num_field_criteria(selectors.size()*balanceSettings.getNumCriteria());

    stk::mesh::Selector selectUnion = stk::mesh::selectUnion(selectors);

    // set vertex weights using entity's topology and if search is part of algorithm, use multiplier
    fill_zoltan2_graph(balanceSettings, stkMeshBulkData, zoltan2Graph, selectUnion, localIds);

    // now can reset those vertex weights based on fields or other critieria
    zoltan2Graph.adjust_vertex_weights(balanceSettings, stkMeshBulkData, selectors, localIds);

    if (balanceSettings.allowModificationOfVertexWeightsForSmallMeshes())
    {
        bool isSmallMesh = (counts[stk::topology::ELEM_RANK] / stkMeshBulkData.parallel_size()) <= 10;
        if(isSmallMesh)
        {
            logMessage(stkMeshBulkData.parallel(), "Changing weights since mesh is small");
            zoltan2Graph.adjust_weights_for_small_meshes();
        }
    }
}

void fill_decomp_using_parmetis(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate, stk::mesh::EntityProcVec &decomp, stk::mesh::BulkData& stkMeshBulkData,
                                const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds)
{
    Teuchos::ParameterList params = getGraphBasedParameters(balanceSettings, numSubdomainsToCreate);

    std::ostringstream os;
    os << "Using Zoltan2 version: " << Zoltan2::Zoltan2_Version();
    logMessage(stkMeshBulkData.parallel(), os.str());
    logMessage(stkMeshBulkData.parallel(), "Filling in graph data");


    Zoltan2ParallelGraph zoltan2Graph;
    createZoltanParallelGraph(balanceSettings, stkMeshBulkData, selectors, localIds, zoltan2Graph);

    std::vector<double> copyOrigWeights = zoltan2Graph.get_vertex_weights();
    std::vector<int> all_local_ids(copyOrigWeights.size());
    for(size_t i=0;i<all_local_ids.size();++i)
    {
        all_local_ids[i]=i;
    }

    decomp.clear();
    decomp.resize(copyOrigWeights.size(), std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));

    #if defined(WRITE_OUT_DEBUGGING_INFO) || defined(WRITE_OUT_DECOMP_METRICS)
    static int step = 0;
    stk::mesh::EntityVector entities = save_for_debugging_local_ids(stkMeshBulkData);
    #endif

    if (balanceSettings.isMultiCriteriaRebalance())
    {
        stk::mesh::Selector selectUnion = stk::mesh::selectUnion(selectors);
        get_multicriteria_parmetis_decomp(stkMeshBulkData, balanceSettings, zoltan2Graph, params, selectUnion, decomp, localIds);
    }
    else
    {
        for(size_t i=0;i<selectors.size();++i)
        {
            zoltan2Graph.set_vertex_weights(copyOrigWeights);

            std::vector<int> localIdsOfNotSelectedEntities = getLocalIdsOfEntitiesNotSelected(stkMeshBulkData, selectors[i], localIds);

            size_t numItemsToRebalance = copyOrigWeights.size() - localIdsOfNotSelectedEntities.size();
            size_t numVerticesGlobal = 0;
            stk::all_reduce_sum(stkMeshBulkData.parallel(), &numItemsToRebalance, &numVerticesGlobal, 1);

            if(numVerticesGlobal>0)
            {
                for(size_t j=0;j<localIdsOfNotSelectedEntities.size();++j)
                {
                    zoltan2Graph.set_vertex_weight(localIdsOfNotSelectedEntities[j], 0.0);
                }

                get_multicriteria_parmetis_decomp(stkMeshBulkData, balanceSettings, zoltan2Graph, params, selectors[i], decomp, localIds);
            }
        }
    }

    logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");

    #if defined(WRITE_OUT_DECOMP_METRICS)
    std::vector<double> weights_per_proc(numSubdomainsToCreate,0);
    double max = 0;
    double sum = 0;
    for(size_t i=0; i<zoltan2Graph.mVertexIds.size();++i)
    {
        sum += copyOrigWeights[i];
        max = std::max(max, copyOrigWeights[i]);
        weights_per_proc[decomp[i]] += copyOrigWeights[i];
    }

    double global_sum = 0;
    stk::all_reduce_sum(stkMeshBulkData.parallel(), &sum, &global_sum, 1);
    double global_max = 0;
    stk::all_reduce_max(stkMeshBulkData.parallel(), &max, &global_max, 1);
    double max_weight_any_proc = 0;
    stk::all_reduce_max(stkMeshBulkData.parallel(), &sum, &max_weight_any_proc, 1);

    std::vector<double> weights_per_proc_sum(weights_per_proc.size(),0);
    stk::all_reduce_sum(stkMeshBulkData.parallel(), weights_per_proc.data(), weights_per_proc_sum.data(), weights_per_proc.size());

    {
        std::ostringstream filename;
        filename << "balance_info_" << stkMeshBulkData.parallel_rank() << ".txt";
        std::ofstream out(filename.str().c_str(), std::ofstream::app);
        std::ostringstream os;
        os << "=========================== For Step " << step << " ====================================" << std::endl;
        os << "Max element weight: " << global_max << " and average weight per proc: " << global_sum/stkMeshBulkData.parallel_size() << " and max total weight any proc = " << max_weight_any_proc << std::endl;
        os << "Max element weight/Ave weight per proc: " << global_max/std::max(0.000001, global_sum/stkMeshBulkData.parallel_size()) << std::endl;
        os << "Weight before this proc: " << sum << " and weight after decomp: " << weights_per_proc_sum[stkMeshBulkData.parallel_rank()] << std::endl;
        os << "Weight/average before: " << sum/(std::max(0.000001, global_sum/stkMeshBulkData.parallel_size())) << " and after: " << weights_per_proc_sum[stkMeshBulkData.parallel_rank()]/(std::max(0.000001, global_sum/stkMeshBulkData.parallel_size())) << std::endl;
        out << os.str();
        out.close();
    }
    #endif

    #if defined(WRITE_OUT_DEBUGGING_INFO) || defined(WRITE_OUT_DECOMP_METRICS)
    step++;
    #endif
}

bool is_geometric_method(const std::string& method)
{
  return (method=="rcb" ||
          method=="rib" ||
          method=="multijagged");
}

bool is_graph_based_method(const std::string& method)
{
  return (method == "parmetis");
}

void calculateGeometricOrGraphBasedDecomp(const BalanceSettings& balanceSettings,
                                              const int numSubdomainsToCreate,
                                              stk::mesh::EntityProcVec &decomp,
                                              stk::mesh::BulkData& stkMeshBulkData,
                                              const std::vector<stk::mesh::Selector>& selectors)
{
    ThrowRequireWithSierraHelpMsg(numSubdomainsToCreate > 0);
    ThrowRequireWithSierraHelpMsg(is_geometric_method(balanceSettings.getDecompMethod()) ||
                                  is_graph_based_method(balanceSettings.getDecompMethod()));

    if (is_geometric_method(balanceSettings.getDecompMethod()))
    {
        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);
        fill_decomp_using_geometric_method(balanceSettings, numSubdomainsToCreate,
                                           decomp, stkMeshBulkData, selectors, localIds);
    }
    else if (is_graph_based_method(balanceSettings.getDecompMethod()))
    {
        stk::mesh::Ghosting * customAura = stk::tools::create_custom_aura(stkMeshBulkData,
                                                                          stkMeshBulkData.mesh_meta_data().globally_shared_part(),
                                                                          "customAura");
        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);
        internal::fill_spider_connectivity_count_fields(stkMeshBulkData, balanceSettings);
        fill_decomp_using_parmetis(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors, localIds);
        stk::tools::destroy_custom_aura(stkMeshBulkData, customAura);
    }
}

bool compareEntityEqualityOnly(const std::pair<stk::mesh::Entity,int> &a, const std::pair<stk::mesh::Entity,int> &b)
{
    return (a.first == b.first);
}

void add_if_owned(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity, int newOwningProc, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    if(stkMeshBulkData.bucket(entity).owned())
        entityProcPairs.emplace_back(entity, newOwningProc);
}

void add_connected_entities_of_rank(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element, int newOwningProc, stk::mesh::EntityRank rank, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    unsigned numEntities = stkMeshBulkData.num_connectivity(element, rank);
    const stk::mesh::Entity *entities = stkMeshBulkData.begin(element, rank);
    for(unsigned int i = 0; i < numEntities; i++)
    {
        add_if_owned(stkMeshBulkData, entities[i], newOwningProc, entityProcPairs);
    }
}

void fillEntityProcPairsForEntityAndItsClosure(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity elementToMove, int newOwningProc,
        std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    entityProcPairs.emplace_back(elementToMove, newOwningProc);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::FACE_RANK, entityProcPairs);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::EDGE_RANK, entityProcPairs);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::NODE_RANK, entityProcPairs);
}

void performModifications(stk::mesh::BulkData& stkMeshBulkData, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    std::sort(entityProcPairs.begin(), entityProcPairs.end());
    std::vector<std::pair<stk::mesh::Entity, int> >::iterator iter = std::unique(entityProcPairs.begin(), entityProcPairs.end(), compareEntityEqualityOnly);
    entityProcPairs.erase(iter, entityProcPairs.end());

    stkMeshBulkData.change_entity_owner(entityProcPairs);
}

void rebalance(DecompositionChangeList &changeList)
{
    stk::mesh::EntityProcVec entityProcPairs;
    changeList.pImpl->get_decomposition_with_full_closure(entityProcPairs);
    performModifications(changeList.get_bulk(), entityProcPairs);
}

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec entityProcPairs;
    for(const stk::mesh::EntityProc& entity_proc : decomposition)
    {
        if(entity_proc.second != stkMeshBulkData.parallel_rank())
            fillEntityProcPairsForEntityAndItsClosure(stkMeshBulkData, entity_proc.first, entity_proc.second, entityProcPairs);
    }
    performModifications(stkMeshBulkData, entityProcPairs);
}

stk::mesh::EntityProcVec get_mapped_decomp(const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec mapped_decomp(decomposition.size());
    for(size_t i=0;i<decomposition.size();++i)
    {
        mapped_decomp[i].first = decomposition[i].first;
        mapped_decomp[i].second = mappings[decomposition[i].second];
    }
    return mapped_decomp;
}

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec mapped_decomp = get_mapped_decomp(mappings, decomposition);
    rebalance(stkMeshBulkData, mapped_decomp);
}

void print_rebalance_metrics(const size_t num_global_entity_migrations, const size_t max_global_entity_migrations, stk::mesh::BulkData & stkMeshBulkData)
{
    std::vector<size_t> meshCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, meshCounts);
    double fraction_of_mesh_moved = static_cast<double>(num_global_entity_migrations)/meshCounts[stk::topology::ELEM_RANK];
    double avg_global_entity_migrations = static_cast<double>(num_global_entity_migrations)/stkMeshBulkData.parallel_size();

    std::ostringstream oss;
    oss << "Percentage of total mesh moved = " << fraction_of_mesh_moved*100.0 << "%";
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(),oss.str());

    oss.str("");
    oss << "Max/Avg global entity migrations = " << max_global_entity_migrations/avg_global_entity_migrations;
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(),oss.str());
}

} //internal

}
}
