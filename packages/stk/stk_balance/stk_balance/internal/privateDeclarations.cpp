#include <cmath>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/GeometricVertices.hpp>
#include <stk_balance/internal/StkGeometricMethodViaZoltan.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_balance/internal/Diagnostics.hpp>
#include <stk_balance/internal/DiagnosticsContainer.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>  // for field_data
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Comm.hpp"
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
#include <stk_math/SideGeometry.hpp>
#include <zoltan.h>
#include <Zoltan2_Version.hpp>
#include <map>

#include "stk_mesh/base/FieldParallel.hpp"
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <algorithm>
#include <tuple>
#include <vector>
#include <numeric>

namespace stk {
namespace balance {

DecompositionChangeList::DecompositionChangeList(stk::mesh::BulkData & bulk)
  : m_bulk(&bulk)
{
}

DecompositionChangeList::DecompositionChangeList(stk::mesh::BulkData & bulk,
                                                 const stk::mesh::EntityProcVec & decomposition)
  : m_bulk(&bulk)
{
  fill_change_list_from_raw_decomposition(decomposition);
}

bool
DecompositionChangeList::has_entity(stk::mesh::Entity entity) const
{
  return (m_dataMap.find(entity) != m_dataMap.end());
}

int
DecompositionChangeList::get_entity_destination(stk::mesh::Entity entity) const
{
  std::map<stk::mesh::Entity, int>::const_iterator entry = m_dataMap.find(entity);
  if (entry != m_dataMap.end()) {
    return entry->second;
  }
  return -1;
}

void
DecompositionChangeList::set_entity_destination(stk::mesh::Entity entity, const int destination)
{
  m_dataMap[entity] = destination;
}

void
DecompositionChangeList::delete_entity(stk::mesh::Entity entity)
{
  std::map<stk::mesh::Entity, int>::iterator entry = m_dataMap.find(entity);

  if (entry != m_dataMap.end()) {
    m_dataMap.erase(entry);
  }
}

stk::mesh::EntityProcVec
DecompositionChangeList::get_all_partition_changes()
{
  stk::mesh::EntityProcVec entityProcPairs;
  for (auto & entry : m_dataMap) {
    entityProcPairs.push_back(entry);
  }
  return entityProcPairs;
}

stk::mesh::BulkData &
DecompositionChangeList::get_bulk()
{
  return *m_bulk;
}

struct ElementWillChangeProcessor {
  ElementWillChangeProcessor(const stk::mesh::BulkData & bulk)
    : m_bulk(bulk)
  {  }

  bool operator()(const std::pair<stk::mesh::Entity, int> & entityProcPair) {
    return (m_bulk.entity_rank(entityProcPair.first) == stk::topology::ELEM_RANK) &&
        (entityProcPair.second != m_bulk.parallel_rank());
  }

  const stk::mesh::BulkData & m_bulk;
};

size_t
DecompositionChangeList::get_num_global_entity_migrations() const
{
  size_t num_local_entity_migrations = std::count_if(m_dataMap.begin(), m_dataMap.end(), ElementWillChangeProcessor(*m_bulk));

  return stk::get_global_sum(m_bulk->parallel(), num_local_entity_migrations);
}

size_t
DecompositionChangeList::get_max_global_entity_migrations() const
{
  size_t num_local_entity_migrations = std::count_if(m_dataMap.begin(), m_dataMap.end(), ElementWillChangeProcessor(*m_bulk));

  return stk::get_global_max(m_bulk->parallel(), num_local_entity_migrations);
}

void
DecompositionChangeList::fill_change_list_from_raw_decomposition(const stk::mesh::EntityProcVec& decomposition)
{
  for (const stk::mesh::EntityProc& entity_proc : decomposition) {
    if (m_bulk->is_valid(entity_proc.first) && (entity_proc.second != m_bulk->parallel_rank())) {
      m_dataMap[entity_proc.first] = entity_proc.second;
    }
  }
}

stk::mesh::EntityProcVec DecompositionChangeList::get_decomposition() const
{
  stk::mesh::EntityProcVec finalDecomposition;
  finalDecomposition.reserve(m_dataMap.size());
  for (auto & entry : m_dataMap) {
    finalDecomposition.push_back(entry);
  }
  return finalDecomposition;
}

stk::mesh::EntityProcVec DecompositionChangeList::get_decomposition_with_full_closure() const
{
  stk::mesh::EntityProcVec finalDecomposition;
  finalDecomposition.reserve(m_dataMap.size());
  for (auto & entry : m_dataMap) {
    finalDecomposition.push_back(entry);
    add_downward_closure_for_entity(entry, finalDecomposition);
  }
  return finalDecomposition;
}

void
DecompositionChangeList::add_downward_closure_for_entity(const std::pair<const stk::mesh::Entity, int> & entityProc,
                                                         stk::mesh::EntityProcVec & finalDecomposition) const
{
  const stk::topology::rank_t entityRank = m_bulk->entity_rank(entityProc.first);
  for (int rank = entityRank-1; rank >= stk::topology::NODE_RANK; --rank)
    internal::add_connected_entities_of_rank(*m_bulk, entityProc.first, entityProc.second,
                                             static_cast<stk::topology::rank_t>(rank), finalDecomposition);
}


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

  STK_ThrowRequireWithSierraHelpMsg(numElements <= 1);
  if ( element != NULL && stkMeshBulkData.bucket(*element).owned() )
  {
    unsigned numNodes = stkMeshBulkData.num_nodes(face);
    const stk::mesh::Entity *nodes = stkMeshBulkData.begin_nodes(face);

    addBoxForNodes(stkMeshBulkData, numNodes, nodes, coord, eps, stkMeshBulkData.identifier(*element), faceBoxes);
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

    double edge_weight = balanceSettings.getGraphEdgeWeightForSearch() * balanceSettings.getGraphEdgeWeightMultiplier();
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

    STK_ThrowRequireWithSierraHelpMsg(owningProcElement1 == bulk.parallel_rank() ||
                                      owningProcElement2 == bulk.parallel_rank());

    bool anyIntersections = false;

    if (element1 == element2) {
      anyIntersections = true;
    }
    else if (bulk.is_valid(element1) && bulk.is_valid(element2)) {
      anyIntersections = internal::has_common_nodes_between_elements(bulk, element1, element2);
    }

    if (!anyIntersections) {
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
using SideGeometryPtr = std::unique_ptr<stk::math::SideGeometry>;

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
    return SideGeometryPtr(new stk::math::QuadGeometry(sideInfo.nodeCoordinates[0],
                                                       sideInfo.nodeCoordinates[1],
                                                       sideInfo.nodeCoordinates[2],
                                                       sideInfo.nodeCoordinates[3]));
  }
  else if (stk::is_tri_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new stk::math::TriGeometry(sideInfo.nodeCoordinates[0],
                                                      sideInfo.nodeCoordinates[1],
                                                      sideInfo.nodeCoordinates[2]));
  }
  else if (is_line_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new stk::math::LineGeometry(sideInfo.nodeCoordinates[0],
                                                       sideInfo.nodeCoordinates[1]));
  }
  else if (is_point_side(sideInfo.sideTopology)) {
    return SideGeometryPtr(new stk::math::PointGeometry(sideInfo.nodeCoordinates[0]));
  }
  else {
    STK_ThrowErrorMsg("Unsupported side topology: " << stk::topology(sideInfo.sideTopology).name());
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

  stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(coords,
    [&](auto& coordsData) {
      for (const stk::mesh::Entity & node : sideNodes) {
        auto nodeFieldData = coordsData.entity_values(node);
        if (spatialDim == 3) {
          nodeCoordinates.emplace_back(nodeFieldData(0_comp), nodeFieldData(1_comp), nodeFieldData(2_comp));
        }
        else if (spatialDim == 2) {
          nodeCoordinates.emplace_back(nodeFieldData(0_comp), nodeFieldData(1_comp), 0.0);
        }
        else {
          STK_ThrowErrorMsg("Problem dimensionality " << spatialDim << " not supported");
        }
      }
    }
  );

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

    // FIXME SHELL_SIDE_TOPOLOGY
    if (elemTopology.is_shell() && sideTopology != elemTopology.side_topology()) { break; }

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

void fillEntityCentroid(const stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::FieldBase* coord,
                        stk::mesh::Entity entityOfConcern, double *entityCentroid)
{
  int spatialDimension = stkMeshBulkData.mesh_meta_data().spatial_dimension();

  stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*coord,
    [&](auto& coordData) {
      if (stkMeshBulkData.entity_rank(entityOfConcern)==stk::topology::NODE_RANK) {
        auto nodeCoord = coordData.entity_values(entityOfConcern);
        for (stk::mesh::ComponentIdx k(0); k < spatialDimension; ++k) {
          entityCentroid[k] += nodeCoord(k);
        }
      }
      else {
        stk::mesh::Entity const* nodes = stkMeshBulkData.begin_nodes(entityOfConcern);
        unsigned numNodes = stkMeshBulkData.num_nodes(entityOfConcern);
        for (unsigned nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex) {
          auto nodeCoord = coordData.entity_values(nodes[nodeIndex]);

          for (stk::mesh::ComponentIdx k(0); k < spatialDimension; ++k) {
            entityCentroid[k] += nodeCoord(k);
          }
        }

        for (int k = 0; k < spatialDimension; ++k) {
          entityCentroid[k] /= numNodes;
        }
      }
    }
  );
}

bool is_not_part_of_spider(const stk::mesh::BulkData & bulk, const stk::mesh::Part & spiderPart,
                           stk::mesh::Entity element)
{
  return not bulk.bucket(element).member(spiderPart);
}

bool is_not_spider_topology(stk::topology::topology_t elemTopology)
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
    if (is_not_spider_topology(elemTopology)) {
        ++numVolumeElems;
      }
    }
  }

  return numVolumeElems;
}



void register_internal_fields_and_parts(stk::mesh::BulkData & bulk, const stk::balance::BalanceSettings & balanceSettings)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  const bool hasSpiderField = meta.get_field(stk::topology::ELEM_RANK,
                                             balanceSettings.getSpiderVolumeConnectivityCountFieldName()) != nullptr;

  if (balanceSettings.shouldFixSpiders() && not hasSpiderField) {
    meta.enable_late_fields();
    stk::mesh::Field<int> & volumeField = meta.declare_field<int>(stk::topology::ELEM_RANK,
                                                                  balanceSettings.getSpiderVolumeConnectivityCountFieldName());
    stk::mesh::put_field_on_mesh(volumeField, meta.universal_part(), nullptr);

    stk::mesh::Field<int> & outputSubdomainField = meta.declare_field<int>(stk::topology::ELEM_RANK,
                                                                           balanceSettings.getOutputSubdomainFieldName());
    stk::mesh::put_field_on_mesh(outputSubdomainField, meta.universal_part(), nullptr);

    meta.declare_part(balanceSettings.getSpiderPartName(), stk::topology::ELEM_RANK);
  }

  const bool hasWeightField = meta.get_field(stk::topology::ELEM_RANK,
                                             balanceSettings.getDiagnosticElementWeightFieldName()) != nullptr;

  if (balanceSettings.shouldPrintDiagnostics() && not hasWeightField) {
    meta.enable_late_fields();
    stk::mesh::Field<double> & weightField = meta.declare_field<double>(stk::topology::ELEM_RANK,
                                                                        balanceSettings.getDiagnosticElementWeightFieldName());
    stk::mesh::put_field_on_mesh(weightField, meta.universal_part(), balanceSettings.getNumCriteria(), nullptr);
  }

  const bool hasConnectivityWeightField = meta.get_field(stk::topology::ELEM_RANK,
                                                         balanceSettings.getVertexConnectivityWeightFieldName()) != nullptr;

  if ((balanceSettings.shouldPrintDiagnostics() ||
       balanceSettings.getVertexWeightMethod() == VertexWeightMethod::CONNECTIVITY) && not hasConnectivityWeightField) {
    meta.enable_late_fields();
    stk::mesh::Field<double> & connectivityWeightField = meta.declare_field<double>(stk::topology::ELEM_RANK,
                                                                                    balanceSettings.getVertexConnectivityWeightFieldName());
    stk::mesh::put_field_on_mesh(connectivityWeightField, meta.universal_part(), nullptr);
  }
}

unsigned number_of_connected_beams(const stk::mesh::BulkData & bulk, stk::mesh::Entity node)
{
  unsigned numBeams = 0;

  const stk::mesh::Entity * elements = bulk.begin_elements(node);
  const unsigned numElems = bulk.num_elements(node);
  for (unsigned i = 0; i < numElems; ++i) {
    stk::topology elemTopology = bulk.bucket(elements[i]).topology();
    if (elemTopology == stk::topology::BEAM_2) {
      ++numBeams;
    }
  }

  return numBeams;
}

bool element_has_an_unconnected_node(const stk::mesh::BulkData & bulk, stk::mesh::Entity element)
{
  const stk::mesh::Entity * nodes = bulk.begin_nodes(element);
  const unsigned numNodes = bulk.num_nodes(element);
  for (unsigned i = 0; i < numNodes; ++i) {
    const unsigned numElems = bulk.num_elements(nodes[i]);
    if (numElems < 2) {  // Only has self
      return true;
    }
  }
  return false;
}

void fill_spider_connectivity_count_fields_and_parts(stk::mesh::BulkData & bulk,
                                                     const BalanceSettings & balanceSettings)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Selector selectBeamsAndParticles(meta.locally_owned_part() &
                                              (meta.get_topology_root_part(stk::topology::BEAM_2) |
                                               meta.get_topology_root_part(stk::topology::PARTICLE)));

  constexpr unsigned spiderConnectivityThreshold = 5;

  stk::mesh::EntityVector spiderElements;
  const auto beamsAndParticles = mesh::get_entities(bulk, stk::topology::ELEM_RANK, selectBeamsAndParticles);
  for (stk::mesh::Entity candidateElem : beamsAndParticles) {
    if (element_has_an_unconnected_node(bulk, candidateElem)) {
      continue;
    }

    const stk::mesh::Entity * nodes = bulk.begin_nodes(candidateElem);
    const unsigned numNodes = bulk.num_nodes(candidateElem);
    for (unsigned i = 0; i < numNodes; ++i) {
      if (number_of_connected_beams(bulk, nodes[i]) > spiderConnectivityThreshold) {
        spiderElements.push_back(candidateElem);
        break;
      }
    }
  }

  stk::mesh::Part * spiderPart = balanceSettings.getSpiderPart(bulk);
  bulk.batch_change_entity_parts(spiderElements, {spiderPart}, {});


  const stk::mesh::Field<int> * volumeConnectivityCountField = balanceSettings.getSpiderVolumeConnectivityCountField(bulk);
  auto volumeConnectivityCountFieldData = volumeConnectivityCountField->data<stk::mesh::ReadWrite>();
  stk::mesh::Selector selectBeamElements(meta.locally_owned_part() &
                                         meta.get_topology_root_part(stk::topology::BEAM_2));
  stk::mesh::EntityVector beams;
  mesh::get_entities(bulk, stk::topology::ELEM_RANK, selectBeamElements, beams);
  for(stk::mesh::Entity beam : beams) {
    auto volumeConnectivityCount = volumeConnectivityCountFieldData.entity_values(beam);
    volumeConnectivityCount() = num_volume_elements_connected_to_beam(bulk, beam);
  }

  stk::mesh::communicate_field_data(bulk, {volumeConnectivityCountField});
}

void fill_output_subdomain_field(const stk::mesh::BulkData & bulk, const BalanceSettings & balanceSettings,
                                 stk::mesh::EntityProcVec & decomp)
{
  if (balanceSettings.shouldFixSpiders()) {
    const stk::mesh::Field<int> * outputSubdomainField = balanceSettings.getOutputSubdomainField(bulk);
    auto outputSubdomainFieldData = outputSubdomainField->data<stk::mesh::ReadWrite>();

    for (const stk::mesh::EntityProc & entityProc : decomp) {
      const stk::mesh::Entity elem = entityProc.first;
      const int proc = entityProc.second;
      if (bulk.is_valid(elem)) {
        auto outputSubdomain = outputSubdomainFieldData.entity_values(elem);
        outputSubdomain() = proc;
      }
    }

    stk::mesh::communicate_field_data(bulk, {outputSubdomainField});
  }
}

void logMessage(MPI_Comm communicator, const std::string &message)
{
  stk::log_with_time_and_memory(communicator, message, stk::outputP0());
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

stk::mesh::EntityVector get_entities_to_balance(stk::mesh::Selector selector,
                                                stk::mesh::EntityRank primaryRank,
                                                const stk::mesh::BulkData& bulkData)
{
  const bool sortById = true;
  stk::mesh::EntityVector entitiesToBalance;
  selector = selector & bulkData.mesh_meta_data().locally_owned_part();
  stk::mesh::get_entities(bulkData, primaryRank, selector, entitiesToBalance, sortById);
  return entitiesToBalance;
}

size_t count_decomp_work_in_this_comm(const stk::mesh::BulkData & bulk,
                                      const stk::ParallelMachine & comm,
                                      const stk::mesh::Selector & selector)
{
  const stk::mesh::Selector locallyOwnedAndSelected = selector & bulk.mesh_meta_data().locally_owned_part();
  const size_t numSelectedLocal = stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, locallyOwnedAndSelected);
  size_t numSelectedGlobal = 0;
  stk::all_reduce_sum(comm, &numSelectedLocal, &numSelectedGlobal, 1);
  return numSelectedGlobal;
}

bool has_decomp_work_in_this_comm(const stk::mesh::BulkData & bulk,
                                  const stk::ParallelMachine & comm,
                                  const stk::mesh::Selector & selector)
{
  return (count_decomp_work_in_this_comm(bulk, comm, selector) > 0);
}

void store_diagnostic_element_weights(const stk::mesh::BulkData & bulk,
                                      const BalanceSettings & balanceSettings,
                                      const stk::mesh::Selector & /*selector*/,
                                      const Vertices & vertices)
{
  if (stk::balance::get_diagnostic<TotalElementWeightDiagnostic>()) {
    const stk::mesh::Field<double> & weightField = *balanceSettings.getDiagnosticElementWeightField(bulk);
    auto weightFieldData = weightField.data<stk::mesh::ReadWrite>();
    const std::vector<double> & vertexWeights = vertices.get_vertex_weights();
    const int numSelectors = 1;
    const int selectorIndex = 0;
    const int numCriteria = balanceSettings.getNumCriteria();

    stk::mesh::EntityVector entitiesToBalance;
    const bool sortById = true;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK,
                            bulk.mesh_meta_data().locally_owned_part(), entitiesToBalance, sortById);

    for (unsigned entityIndex = 0; entityIndex < entitiesToBalance.size(); ++entityIndex) {
      auto weight = weightFieldData.entity_values(entitiesToBalance[entityIndex]);
      for (stk::mesh::ComponentIdx criterion = 0_comp; criterion < numCriteria; ++criterion) {
        unsigned index = stk::balance::internal::get_index(numSelectors, numCriteria, entityIndex,
                                                           selectorIndex, criterion);
        weight(criterion) = vertexWeights[index];
      }
    }
  }
}


void get_multicriteria_decomp_using_selectors_as_segregation(const stk::mesh::BulkData& bulk,
                                                             const std::vector<stk::mesh::Selector>& criterions,
                                                             const stk::ParallelMachine & decompCommunicator,
                                                             const BalanceSettings& balanceSettings,
                                                             const int numSubdomainsToCreate,
                                                             stk::mesh::EntityProcVec& decomp,
                                                             const stk::mesh::impl::LocalIdMapper& localIds)
{
  stk::mesh::Selector unionSelector = stk::mesh::selectUnion(criterions);
  const stk::mesh::Selector locallyOwnedAndSelected = unionSelector & bulk.mesh_meta_data().locally_owned_part();

  stk::mesh::EntityVector entitiesToBalance = get_entities_to_balance(locallyOwnedAndSelected, stk::topology::ELEM_RANK, bulk);
  size_t num_entities = entitiesToBalance.size();
  size_t num_entities_across_procs = 0;
  stk::all_reduce_sum(decompCommunicator, &num_entities, &num_entities_across_procs, 1);

  if (num_entities_across_procs > 0) {
    stk::balance::internal::GeometricVertices vertexInfo(balanceSettings, bulk, entitiesToBalance, criterions);

    store_diagnostic_element_weights(bulk, balanceSettings, unionSelector, vertexInfo);

    std::vector<unsigned> processorOntoWhichEntityBelongs = stk::balance::get_decomposition(vertexInfo, balanceSettings,
                                                                                            numSubdomainsToCreate, decompCommunicator);
    for (size_t i = 0; i < entitiesToBalance.size(); ++i) {
      int local_id = get_local_id(localIds, entitiesToBalance[i]);
      decomp[local_id] = std::make_pair(entitiesToBalance[i], processorOntoWhichEntityBelongs[i]);
    }
  }
}

void fill_decomp_using_geometric_method(stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors,
                                        const stk::ParallelMachine & decompCommunicator, const int numSubdomainsToCreate,
                                        const BalanceSettings& balanceSettings, const stk::mesh::impl::LocalIdMapper& localIds,
                                        stk::mesh::EntityProcVec &decomp)
{
  logMessage(stkMeshBulkData.parallel(), "Using Zoltan2 version: " + Zoltan2::Zoltan2_Version());
  logMessage(stkMeshBulkData.parallel(), "Filling in vertex data for decomp method = " + balanceSettings.getDecompMethod());

  size_t numEntities = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

  decomp.clear();
  decomp.resize(numEntities, std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));


  if (balanceSettings.isMultiCriteriaRebalance())
    get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, selectors, decompCommunicator, balanceSettings, numSubdomainsToCreate, decomp, localIds);
  else
    for(const stk::mesh::Selector& selector : selectors)
      get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, std::vector<stk::mesh::Selector>{selector}, decompCommunicator, balanceSettings, numSubdomainsToCreate, decomp, localIds);

  logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");
}


void get_multicriteria_graph_based_decomp(const stk::mesh::BulkData &mesh,
                                          stk::mesh::Selector selector,
                                          const stk::ParallelMachine & decompCommunicator,
                                          const BalanceSettings& balanceSettings,
                                          const stk::mesh::impl::LocalIdMapper& domainLocalIds,
                                          Zoltan2ParallelGraph &zoltan2Graph,
                                          Teuchos::ParameterList &params,
                                          stk::mesh::EntityProcVec &decomp)

{
  StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

  logMessage(mesh.parallel(), "Setting up partitioning problem");

  Zoltan2::PartitioningProblem<StkMeshZoltanAdapter> problem(&stkMeshAdapter, &params, decompCommunicator);

  logMessage(mesh.parallel(), "Solving");

  if (has_decomp_work_in_this_comm(mesh, decompCommunicator, selector)) {
    std::srand(mesh.parallel_rank()); // KHP: Temporary until an API is added to Zoltan2 for random seeds.
    problem.solve();

    stk::mesh::impl::LocalIdMapper graphLocalIds(mesh, stk::topology::ELEM_RANK, selector);
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(mesh, stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part() & selector, elements);

    const StkMeshZoltanAdapter::part_t *processorOntoWhichEntityBelongs = problem.getSolution().getPartListView();
    for (size_t j = 0; j < elements.size(); ++j) {
      int domainLocalId = get_local_id(domainLocalIds, elements[j]);
      int graphLocalId = get_local_id(graphLocalIds, elements[j]);
      int dest_proc = processorOntoWhichEntityBelongs[graphLocalId];
      decomp[domainLocalId] = std::make_pair(elements[j], dest_proc);
    }
  }

  if (balanceSettings.shouldPrintMetrics()) {
    internal::print_statistics(stkMeshAdapter, mesh.parallel(), mesh.parallel_rank());
    internal::print_solution_statistics(stkMeshAdapter, problem.getSolution(), mesh.parallel(), mesh.parallel_rank());
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
  auto volumeElemConnectivityCountFieldData = volumeElemConnectivityCountField.data();

  for (unsigned i = 0; i < numElements; ++i) {
    if (bulk.bucket(elems[i]).topology() == stk::topology::BEAM_2) {
      auto volumeElemConnectivityCount = volumeElemConnectivityCountFieldData.entity_values(elems[i]);
      if (volumeElemConnectivityCount() == 0) {
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
  for (const stk::mesh::EntityProcMap::value_type & entityNewOwner : newSpiderEntityOwners) {
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

unsigned count_connected_spider_elements(const stk::mesh::BulkData & bulk,
                                         const stk::mesh::Part & spiderPart,
                                         const stk::mesh::Entity node)
{
  const stk::mesh::Entity * connectedElems = bulk.begin_elements(node);
  const unsigned numConnectedElems = bulk.num_elements(node);
  unsigned numSpiderElems = 0;
  for (unsigned i = 0; i < numConnectedElems; ++i) {
    if (bulk.bucket(connectedElems[i]).member(spiderPart)) {
      ++numSpiderElems;
    }
  }
  return numSpiderElems;
}

std::pair<stk::mesh::Entity, stk::mesh::Entity> get_foot_and_body_nodes(const stk::mesh::BulkData & bulk,
                                                                        const stk::mesh::Part & spiderPart,
                                                                        const stk::mesh::Entity spiderLeg)
{
  const stk::mesh::Entity* nodes = bulk.begin_nodes(spiderLeg);
  const int node1ConnectivityCount = count_connected_spider_elements(bulk, spiderPart, nodes[0]);
  const int node2ConnectivityCount = count_connected_spider_elements(bulk, spiderPart, nodes[1]);

  return (node1ConnectivityCount < node2ConnectivityCount) ? std::make_pair(nodes[0], nodes[1])
                                                           : std::make_pair(nodes[1], nodes[0]);
}

stk::mesh::EntityVector get_local_spider_legs(const stk::mesh::BulkData & bulk,
                                              const stk::mesh::Part & spiderPart)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Selector localSpiderLegsSelector = spiderPart &
                                                meta.get_topology_root_part(stk::topology::BEAM_2) &
                                                meta.locally_owned_part();
  return stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, localSpiderLegsSelector);
}

int get_new_spider_leg_owner(const stk::mesh::BulkData & bulk, const stk::mesh::Entity & footNode,
                             const stk::mesh::Part & spiderPart, const stk::mesh::Field<int> & outputSubdomainField)
{
  auto outputSubdomainFieldData = outputSubdomainField.data();
  const stk::mesh::Entity* footElements = bulk.begin_elements(footNode);
  const unsigned numElements = bulk.num_elements(footNode);
  int newLegOwner = std::numeric_limits<int>::max();

  for (unsigned i = 0; i < numElements; ++i) {
    if (is_not_part_of_spider(bulk, spiderPart, footElements[i])) {
      auto outputSubdomain = outputSubdomainFieldData.entity_values(footElements[i]);
      newLegOwner = std::min(newLegOwner, outputSubdomain());
    }
  }

  return newLegOwner;
}

void update_spider_body_component_owners(const stk::mesh::BulkData & bulk,
                                         const BalanceSettings & balanceSettings,
                                         const stk::mesh::Entity & bodyNode,
                                         const stk::mesh::Entity & footNode,
                                         int newLegOwner,
                                         stk::mesh::EntityProcMap & newBodyComponentOwners)
{
  update_new_spider_entity_owner(newBodyComponentOwners, bodyNode, newLegOwner);

  const stk::mesh::Entity spiderParticleBody = get_spider_particle_body_for_leg(bulk, bodyNode);
  if (spider_body_element_exists(bulk, spiderParticleBody)) {
    update_new_spider_entity_owner(newBodyComponentOwners, spiderParticleBody, newLegOwner);
  }
  else {
    const std::pair<stk::mesh::Entity, stk::mesh::Entity> spiderBeamBodyAndNode =
        get_spider_beam_body_and_node_for_leg(bulk, balanceSettings, bodyNode);
    if (spider_body_element_exists(bulk, spiderBeamBodyAndNode.first)) {
      update_new_spider_entity_owner(newBodyComponentOwners, spiderBeamBodyAndNode.first, newLegOwner);
      update_new_spider_entity_owner(newBodyComponentOwners, spiderBeamBodyAndNode.second, newLegOwner);
    }
  }

  if (newLegOwner != bulk.parallel_owner_rank(footNode)) {
    update_new_spider_entity_owner(newBodyComponentOwners, footNode, newLegOwner);
  }
}

void update_output_subdomain_field(const stk::mesh::BulkData & bulk,
                                   const stk::mesh::Field<int> & outputSubdomainField,
                                   stk::mesh::Entity spiderEntity,
                                   int newLegOwner)
{
  auto outputSubdomainFieldData = outputSubdomainField.data<stk::mesh::ReadWrite>();
  if (bulk.entity_rank(spiderEntity) == stk::topology::ELEM_RANK) {
    auto outputSubdomain = outputSubdomainFieldData.entity_values(spiderEntity);
    outputSubdomain() = newLegOwner;
  }
}

void fix_spider_elements(const BalanceSettings & balanceSettings, stk::mesh::BulkData & bulk,
                         DecompositionChangeList & changeList)
{
  const stk::mesh::Field<int> & outputSubdomainField = *balanceSettings.getOutputSubdomainField(bulk);
  const stk::mesh::Part & spiderPart = *balanceSettings.getSpiderPart(bulk);

  stk::mesh::EntityVector spiderLegs = get_local_spider_legs(bulk, spiderPart);
  stk::mesh::EntityProcMap newSpiderComponentOwners;

  for (stk::mesh::Entity spiderLeg : spiderLegs) {
    stk::mesh::Entity footNode;
    stk::mesh::Entity bodyNode;
    std::tie(footNode, bodyNode) = get_foot_and_body_nodes(bulk, spiderPart, spiderLeg);

    const int newLegOwner = get_new_spider_leg_owner(bulk, footNode, spiderPart, outputSubdomainField);

    if (found_valid_new_owner(bulk, newLegOwner)) {
      changeList.set_entity_destination(spiderLeg, newLegOwner);
      update_output_subdomain_field(bulk, outputSubdomainField, spiderLeg, newLegOwner);
      update_spider_body_component_owners(bulk, balanceSettings, bodyNode, footNode, newLegOwner, newSpiderComponentOwners);
    }
  }

  stk::mesh::EntityProcMap globalNewSpiderComponentOwners = determine_global_new_owner(bulk, newSpiderComponentOwners);

  for (const auto & entityOwner : globalNewSpiderComponentOwners) {
    const stk::mesh::Entity spiderEntity = entityOwner.first;
    const int newOwner = entityOwner.second;
    changeList.set_entity_destination(spiderEntity, newOwner);
    update_output_subdomain_field(bulk, outputSubdomainField, spiderEntity, newOwner);
  }

  stk::mesh::communicate_field_data(bulk, {&outputSubdomainField});
}

void createZoltanParallelGraph(stk::mesh::BulkData & stkMeshBulkData,
                               const stk::mesh::Selector & selector,
                               const stk::ParallelMachine & decompCommunicator,
                               const BalanceSettings & balanceSettings,
                               Zoltan2ParallelGraph & zoltan2Graph)
{
  std::ostringstream os;
  os << "Using Zoltan2 version: " << Zoltan2::Zoltan2_Version();
  logMessage(stkMeshBulkData.parallel(), os.str());
  logMessage(stkMeshBulkData.parallel(), "Filling in graph data");

  stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK, selector);

  const size_t numSelectedElements = count_decomp_work_in_this_comm(stkMeshBulkData, decompCommunicator, selector);
  int numProcsInCommunicator = 0;
  MPI_Comm_size(decompCommunicator, &numProcsInCommunicator);

  zoltan2Graph.set_num_global_elements(numSelectedElements);
  zoltan2Graph.set_spatial_dim(stkMeshBulkData.mesh_meta_data().spatial_dimension());

  if (balanceSettings.isMultiCriteriaRebalance()) {
    zoltan2Graph.set_num_field_criteria(balanceSettings.getNumCriteria());
  }

  // set vertex weights using entity's topology and if search is part of algorithm, use multiplier
  fill_zoltan2_graph(balanceSettings, stkMeshBulkData, zoltan2Graph, selector, localIds);

  // now can reset those vertex weights based on fields or other critieria
  zoltan2Graph.adjust_vertex_weights(balanceSettings, stkMeshBulkData, selector, localIds);

  if (balanceSettings.allowModificationOfVertexWeightsForSmallMeshes())
  {
    bool isSmallMesh = (numSelectedElements / numProcsInCommunicator) <= 10;
    if (isSmallMesh) {
      logMessage(stkMeshBulkData.parallel(), "Changing weights since mesh is small");
      zoltan2Graph.adjust_weights_for_small_meshes();
    }
  }
}

void fill_decomp_using_graph_based_method(stk::mesh::BulkData & stkMeshBulkData,
                                          const std::vector<stk::mesh::Selector> & selectors,
                                          const stk::ParallelMachine & decompCommunicator,
                                          const int numSubdomainsToCreate,
                                          const BalanceSettings & balanceSettings,
                                          stk::mesh::EntityProcVec & decomp)
{
  const stk::mesh::Selector locallyOwnedSelector = stkMeshBulkData.mesh_meta_data().locally_owned_part();
  stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK, locallyOwnedSelector);

  const size_t numLocallyOwned = stk::mesh::count_entities(stkMeshBulkData, stk::topology::ELEM_RANK, locallyOwnedSelector);
  decomp.clear();
  decomp.resize(numLocallyOwned, std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));

  Teuchos::ParameterList params = getGraphBasedParameters(balanceSettings, numSubdomainsToCreate);

  if (balanceSettings.isMultiCriteriaRebalance()) {
    stk::mesh::Selector selectUnion = stk::mesh::selectUnion(selectors);
    Zoltan2ParallelGraph zoltan2Graph;
    createZoltanParallelGraph(stkMeshBulkData, selectUnion, decompCommunicator, balanceSettings, zoltan2Graph);
    store_diagnostic_element_weights(stkMeshBulkData, balanceSettings, selectUnion, zoltan2Graph);
    get_multicriteria_graph_based_decomp(stkMeshBulkData, selectUnion, decompCommunicator, balanceSettings,
                                         localIds, zoltan2Graph, params, decomp);
  }
  else {
    for (const stk::mesh::Selector & selector : selectors) {
      Zoltan2ParallelGraph zoltan2Graph;
      createZoltanParallelGraph(stkMeshBulkData, selector, decompCommunicator, balanceSettings, zoltan2Graph);
      store_diagnostic_element_weights(stkMeshBulkData, balanceSettings, selector, zoltan2Graph);
      get_multicriteria_graph_based_decomp(stkMeshBulkData, selector, decompCommunicator, balanceSettings,
                                           localIds, zoltan2Graph, params, decomp);
    }
  }

  logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");
}

void collapse_to_serial_partition(stk::mesh::BulkData & stkMeshBulkData,
                                  const std::vector<stk::mesh::Selector> & decompSelectors,
                                  stk::mesh::EntityProcVec & decomp)
{
  const stk::mesh::Selector locallyOwnedSelector = stkMeshBulkData.mesh_meta_data().locally_owned_part();
  stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK, locallyOwnedSelector);

  constexpr int RootPartition = 0;
  const size_t numLocallyOwned = stk::mesh::count_entities(stkMeshBulkData, stk::topology::ELEM_RANK, locallyOwnedSelector);
  decomp.clear();
  decomp.resize(numLocallyOwned, std::make_pair(stk::mesh::Entity(), RootPartition));

  for (const stk::mesh::Selector & decompSelector : decompSelectors) {
    const stk::mesh::Selector localSelected = locallyOwnedSelector & decompSelector;
    for (const stk::mesh::Bucket * bucket : stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, localSelected)) {
      for (const stk::mesh::Entity & element : *bucket) {
        unsigned localId = get_local_id(localIds, element);
        decomp[localId].first = element;
      }
    }
  }
}

bool is_geometric_method(const std::string& method)
{
  return (method=="rcb" ||
          method=="rib" ||
          method=="multijagged");
}

bool is_graph_based_method(const std::string& method)
{
  return (method == "parmetis" ||
          method == "scotch");
}

void calculateGeometricOrGraphBasedDecomp(stk::mesh::BulkData & stkMeshBulkData,
                                          const std::vector<stk::mesh::Selector> & selectors,
                                          const stk::ParallelMachine & decompCommunicator,
                                          const int numSubdomainsToCreate,
                                          const BalanceSettings & balanceSettings,
                                          stk::mesh::EntityProcVec & decomp)
{
  STK_ThrowRequireWithSierraHelpMsg(numSubdomainsToCreate > 0);
  STK_ThrowRequireWithSierraHelpMsg(is_geometric_method(balanceSettings.getDecompMethod()) ||
                                    is_graph_based_method(balanceSettings.getDecompMethod()));

  std::vector<stk::mesh::Selector> balanceSelectors = selectors;
  if (balanceSettings.shouldFixSpiders()) {
    internal::fill_spider_connectivity_count_fields_and_parts(stkMeshBulkData, balanceSettings);

    const stk::mesh::Part & spiderPart = *balanceSettings.getSpiderPart(stkMeshBulkData);
    for (stk::mesh::Selector & selector : balanceSelectors) {
      selector = selector & !spiderPart;
    }
  }

  compute_connectivity_weight(stkMeshBulkData, balanceSettings);

  if (numSubdomainsToCreate > 1) {
    if (is_geometric_method(balanceSettings.getDecompMethod())) {
      stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);
      fill_decomp_using_geometric_method(stkMeshBulkData, balanceSelectors, decompCommunicator, numSubdomainsToCreate,
                                         balanceSettings, localIds, decomp);
    }
    else if (is_graph_based_method(balanceSettings.getDecompMethod())) {
      fill_decomp_using_graph_based_method(stkMeshBulkData, balanceSelectors, decompCommunicator, numSubdomainsToCreate,
                                           balanceSettings, decomp);
    }
  }
  else {
    collapse_to_serial_partition(stkMeshBulkData, selectors, decomp);
  }

  internal::fill_output_subdomain_field(stkMeshBulkData, balanceSettings, decomp);
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
  stk::mesh::EntityProcVec entityProcPairs = changeList.get_decomposition_with_full_closure();
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

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings,
               const stk::mesh::EntityProcVec& decomposition)
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

void compute_element_count_diagnostic(ElementCountDiagnostic & diag, const stk::mesh::BulkData & bulk, int rank)
{
  const unsigned numElements = stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK,
                                                         bulk.mesh_meta_data().locally_owned_part());
  diag.store_value(rank, numElements);
}

void compute_total_element_weight_diagnostic(TotalElementWeightDiagnostic & diag, const stk::mesh::BulkData & bulk,
                                             const stk::balance::BalanceSettings & balanceSettings,
                                             const stk::mesh::Field<double> & weightField, int rank)
{
  auto weightFieldData = weightField.data();
  const stk::mesh::BucketVector & buckets = bulk.get_buckets(stk::topology::ELEM_RANK,
                                                             bulk.mesh_meta_data().locally_owned_part());
  const int numCriteria = balanceSettings.getNumCriteria();
  std::vector<double> weightSum(numCriteria, 0.0);

  for (const stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & element : *bucket) {
      auto weight = weightFieldData.entity_values(element);
      for (stk::mesh::ComponentIdx criterion : weight.components()) {
        weightSum[criterion] += weight(criterion);
      }
    }
  }

  for (int criterion = 0; criterion < numCriteria; ++criterion) {
    diag.store_value(criterion, rank, weightSum[criterion]);
  }
}

void compute_relative_node_interface_size_diagnostic(RelativeNodeInterfaceSizeDiagnostic & diag, const stk::mesh::BulkData & bulk, int rank)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Selector sharedNodesSelector = meta.globally_shared_part();
  stk::mesh::Selector allNodesSelector = meta.locally_owned_part() | meta.globally_shared_part();

  size_t numSharedNodes = stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, sharedNodesSelector);
  size_t numNodes = stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, allNodesSelector);

  diag.store_value(rank, (numNodes > 0) ? (static_cast<double>(numSharedNodes) / numNodes)
                                        : 0.0);
}

double getTypicalElemsPerNode(stk::topology type)
{
  switch(type)
  {
  case stk::topology::PARTICLE:
    return 1;
  case stk::topology::LINE_2_1D:
    return 1;
  case stk::topology::LINE_3_1D:
    return 1.0/2.0;
  case stk::topology::BEAM_2:
    return 1;
  case stk::topology::BEAM_3:
    return 1.0/2.0;
  case stk::topology::SHELL_LINE_2:
    return 1;
  case stk::topology::SHELL_LINE_3:
    return 1.0/2.0;
  case stk::topology::SPRING_2:
    return 1;
  case stk::topology::SPRING_3:
    return 1.0/2.0;
  case stk::topology::TRI_3_2D:
    return 2;
  case stk::topology::TRI_4_2D:
    return 2.0/3.0;
  case stk::topology::TRI_6_2D:
    return 2.0/4.0;
  case stk::topology::SHELL_TRI_3:
    return 2;
  case stk::topology::SHELL_TRI_4:
    return 2.0/3.0;
  case stk::topology::SHELL_TRI_6:
    return 2.0/4.0;
  case stk::topology::QUAD_4_2D:
    return 1;
  case stk::topology::QUAD_8_2D:
    return 1.0/3.0;
  case stk::topology::QUAD_9_2D:
    return 1.0/4.0;
  case stk::topology::SHELL_QUAD_4:
    return 1;
  case stk::topology::SHELL_QUAD_8:
    return 1.0/3.0;
  case stk::topology::SHELL_QUAD_9:
    return 1.0/4.0;
  case stk::topology::TET_4:
    return 6;
  case stk::topology::TET_8:
    return 6.0/13.0;
  case stk::topology::TET_10:
    return 6.0/8.0;
  case stk::topology::TET_11:
    return 6.0/14.0;
  case stk::topology::PYRAMID_5:
    return 6.0/2.0;
  case stk::topology::PYRAMID_13:
    return 6.0/13.0;
  case stk::topology::PYRAMID_14:
    return 6.0/19.0;
  case stk::topology::WEDGE_6:
    return 2;
  case stk::topology::WEDGE_12:
    return 2.0/4.0;
  case stk::topology::WEDGE_15:
    return 2.0/5.0;
  case stk::topology::WEDGE_18:
    return 2.0/8.0;
  case stk::topology::HEXAHEDRON_8:
    return 1;
  case stk::topology::HEXAHEDRON_20:
    return 1.0/4.0;
  case stk::topology::HEXAHEDRON_27:
    return 1.0/8.0;
  default:
    if ( type.is_superelement( ))
    {
      return 1.0/100.0;
    }
    throw("Invalid Element Type In WeightsOfElement");
  }
}

double get_connected_node_weight(const stk::mesh::BulkData & bulk, std::vector<stk::mesh::Entity> & connectedNodesBuffer,
                                 const stk::mesh::Entity node)
{
  connectedNodesBuffer.clear();
  const unsigned numElements = bulk.num_elements(node);
  const stk::mesh::Entity * elements = bulk.begin_elements(node);

  for (unsigned elemIndex = 0; elemIndex < numElements; ++elemIndex) {
    const stk::mesh::Entity element = elements[elemIndex];
    const unsigned numElementNodes = bulk.num_nodes(element);
    const stk::mesh::Entity * elementNodes = bulk.begin_nodes(element);

    for (unsigned nodeIndex = 0; nodeIndex < numElementNodes; ++nodeIndex) {
      const stk::mesh::Entity elementNode = elementNodes[nodeIndex];
      connectedNodesBuffer.push_back(elementNode);
    }
  }

  stk::util::sort_and_unique(connectedNodesBuffer);

  return connectedNodesBuffer.size();
}

void spread_weight_across_connected_elements(const stk::mesh::BulkData & bulk, const stk::mesh::Entity & node,
                                             double nodeWeight, const stk::mesh::Field<double> & elementWeights)
{
  auto elementWeightsData = elementWeights.data<stk::mesh::ReadWrite>();
  const unsigned numElements = bulk.num_elements(node);
  const stk::mesh::Entity * elements = bulk.begin_elements(node);

  for (unsigned elemIndex = 0; elemIndex < numElements; ++elemIndex) {
    const stk::mesh::Entity element = elements[elemIndex];
    if (bulk.bucket(element).owned()) {
      auto elemWeight = elementWeightsData.entity_values(element);
      const unsigned numNodes = bulk.num_nodes(element);
      const double typicalElemsPerNode = getTypicalElemsPerNode(bulk.bucket(element).topology());
      elemWeight() += nodeWeight / (numNodes * typicalElemsPerNode);
    }
  }
}

void compute_connectivity_weight(const stk::mesh::BulkData & bulk,  const BalanceSettings & balanceSettings)
{
  if (not (balanceSettings.shouldPrintDiagnostics() ||
           balanceSettings.getVertexWeightMethod() == VertexWeightMethod::CONNECTIVITY)) {
    return;
  }

  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Selector allNodesSelector = meta.locally_owned_part() | meta.globally_shared_part();
  std::vector<stk::mesh::Entity> connectedNodesBuffer;
  const stk::mesh::Field<double> & elementWeights = *balanceSettings.getVertexConnectivityWeightField(bulk);
  stk::mesh::field_fill(0.0, elementWeights);

  const stk::mesh::BucketVector & buckets = bulk.get_buckets(stk::topology::NODE_RANK, allNodesSelector);
  for (const stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & node : *bucket) {
      const double nodeWeight = get_connected_node_weight(bulk, connectedNodesBuffer, node);
      spread_weight_across_connected_elements(bulk, node, nodeWeight, elementWeights);
    }
  }
}

void compute_connectivity_weight_diagnostic(ConnectivityWeightDiagnostic & diag, const stk::mesh::BulkData & bulk,
                                            const BalanceSettings & balanceSettings, int rank)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const stk::mesh::Field<double> & elementWeights = *balanceSettings.getVertexConnectivityWeightField(bulk);
  auto elementWeightsData = elementWeights.data();

  double totalWeight = 0.0;
  const stk::mesh::BucketVector & buckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  for (const stk::mesh::Bucket * bucket : buckets) {
    for (const stk::mesh::Entity & element : *bucket) {
      auto elemWeight = elementWeightsData.entity_values(element);
      totalWeight += elemWeight();
    }
  }

  diag.store_value(rank, totalWeight);
}

void compute_balance_diagnostics(const stk::mesh::BulkData & bulk, const stk::balance::BalanceSettings & balanceSettings)
{
  auto * elemCountDiag = get_diagnostic<ElementCountDiagnostic>();
  if (elemCountDiag) {
    compute_element_count_diagnostic(*elemCountDiag, bulk, bulk.parallel_rank());
  }

  auto * totalElemWeightDiag = get_diagnostic<TotalElementWeightDiagnostic>();
  if (totalElemWeightDiag) {
    const stk::mesh::Field<double> & weightField = *balanceSettings.getDiagnosticElementWeightField(bulk);
    compute_total_element_weight_diagnostic(*totalElemWeightDiag, bulk, balanceSettings, weightField, bulk.parallel_rank());
  }

  auto * relNodeInterfaceSizeDiag = get_diagnostic<RelativeNodeInterfaceSizeDiagnostic>();
  if (relNodeInterfaceSizeDiag) {
    compute_relative_node_interface_size_diagnostic(*relNodeInterfaceSizeDiag, bulk, bulk.parallel_rank());
  }

  auto * connectivityWeightDiag = get_diagnostic<ConnectivityWeightDiagnostic>();
  if (connectivityWeightDiag) {
    compute_connectivity_weight_diagnostic(*connectivityWeightDiag, bulk, balanceSettings, bulk.parallel_rank());
  }
}

EnableAura::EnableAura(stk::mesh::BulkData & bulk)
  : m_bulk(bulk),
    m_weTurnedOnAura(false)
{
  if (not bulk.is_automatic_aura_on()) {
    const bool applyImmediately = true;
    m_bulk.set_automatic_aura_option(stk::mesh::BulkData::AUTO_AURA, applyImmediately);
    m_weTurnedOnAura = true;
  }
}

EnableAura::~EnableAura()
{
  if (m_weTurnedOnAura) {
    const bool applyImmediately = true;
    m_bulk.set_automatic_aura_option(stk::mesh::BulkData::NO_AUTO_AURA, applyImmediately);
  }
}

} //internal

}
}
