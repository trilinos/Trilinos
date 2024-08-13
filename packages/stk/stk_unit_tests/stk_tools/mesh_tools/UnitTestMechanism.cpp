#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include "stk_unit_test_utils/getOption.h"
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/GraphCycleDetector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/Comm.hpp"
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_util/util/TreeTraverser.hpp"

#include <limits>
#include <string>
#include <algorithm>
#include <unordered_map>

using stk::unit_test_util::build_mesh;

//#define DEBUG

struct ElementCluster
{
  static constexpr int64_t INVALID_ID = std::numeric_limits<int64_t>::max();

  int64_t globalId = INVALID_ID;
  stk::mesh::EntityVector elements{};

  ElementCluster() : globalId(INVALID_ID) {}
  ElementCluster(stk::mesh::EntityVector& elements_) : globalId(INVALID_ID), elements(std::move(elements_)) {}
};

struct ElementIdCluster
{
  static constexpr int64_t INVALID_ID = -1;

  int64_t globalId = INVALID_ID;
  stk::mesh::EntityIdVector elements{};

  ElementIdCluster() : globalId(INVALID_ID) {}
  ElementIdCluster(int64_t id_) : globalId(id_) {}
  ElementIdCluster(int64_t id_, const stk::mesh::EntityIdVector& elements_) : globalId(id_), elements(elements_) {}
};

struct ParallelClusterRecvInfo
{
  int64_t clusterIndex = -1;
  stk::mesh::EntityId clusterMinElemId = stk::mesh::InvalidEntityId;

  ParallelClusterRecvInfo(int64_t clusterIndex_, stk::mesh::EntityId clusterMinElemId_)
  : clusterIndex(clusterIndex_), clusterMinElemId(clusterMinElemId_) {}
};

struct ParallelClusterInfo
{
  int64_t clusterIndex = -1;
  stk::mesh::EntityId minElemId = stk::mesh::InvalidEntityId;
  std::unordered_map<int, stk::mesh::EntityIdVector> sendData;
  std::unordered_map<int, std::vector<ParallelClusterRecvInfo>> recvData;

  ParallelClusterInfo() : clusterIndex(-1) {}
  ParallelClusterInfo(int64_t index) : clusterIndex(index) {}
};

struct ClusterGraphNode
{
  int proc = -1;
  int64_t clusterIndex = -1;
  stk::mesh::EntityId clusterMinElemId = stk::mesh::InvalidEntityId;

  ClusterGraphNode() : proc(-1), clusterIndex(-1), clusterMinElemId(stk::mesh::InvalidEntityId) {}
  ClusterGraphNode(int proc_, int64_t clusterIndex_, stk::mesh::EntityId clusterMinElemId_)
  : proc(proc_), clusterIndex(clusterIndex_), clusterMinElemId(clusterMinElemId_) {}

  bool operator==(const ClusterGraphNode &rhs) const
  {
    return ((proc == rhs.proc) && (clusterIndex == rhs.clusterIndex) && (clusterMinElemId == rhs.clusterMinElemId));
  }

  bool operator!=(const ClusterGraphNode &rhs) const
  {
    return ((proc != rhs.proc) || (clusterIndex != rhs.clusterIndex) || (clusterMinElemId != rhs.clusterMinElemId));
  }

  bool operator<(const ClusterGraphNode &rhs) const
  {
    if(proc < rhs.proc)
      return true;
    else if (proc == rhs.proc && clusterIndex < rhs.clusterIndex)
      return true;
    else if(proc == rhs.proc && clusterIndex == rhs.clusterIndex) {
      STK_ThrowRequireMsg(clusterMinElemId == rhs.clusterMinElemId,
                      "Two equivalent cluster nodes {proc,clusterIndex} = {"
                      << proc << "," << clusterIndex << "} do not have the same clusterMinElemId: "
                      << clusterMinElemId << " vs " << rhs.clusterMinElemId);
    }

    return false;
  }
};

struct ClusterGraphEdge
{
  ClusterGraphNode node1;
  ClusterGraphNode node2;

  ClusterGraphEdge(const ClusterGraphNode& node1_, const ClusterGraphNode& node2_) : node1(node1_), node2(node2_) {}
  ClusterGraphEdge reciprocal_edge() const {return ClusterGraphEdge(node2, node1);}
};

inline std::ostream &operator<<(std::ostream &out, const ClusterGraphNode &t)
{
  return out << "{proc:" << t.proc << " index:" << t.clusterIndex << " minId:" << t.clusterMinElemId << "}";
}

inline std::ostream &operator<<(std::ostream &out, const ClusterGraphEdge &t)
{
  return out << t.node1 << " <-> " << t.node2;
}

template <class T>
inline void hash_combine(std::size_t & s, const T & v)
{
  std::hash<T> h;
  s ^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
}


namespace std {

template<>
struct hash<ClusterGraphNode>
{
  size_t operator()(ClusterGraphNode const& s) const
  {
    std::size_t res = 0;
    hash_combine(res,s.proc);
    hash_combine(res,s.clusterIndex);
    hash_combine(res,s.clusterMinElemId);
    return res;
  }
};

}

using ElementClusterVector = std::vector<ElementIdCluster>;
using ClusterGraphData = std::map<ClusterGraphNode, std::vector<ClusterGraphNode>>;

template <typename GraphKey, typename GraphData>
class ClusterGraph : public stk::util::TreeTraverser<GraphKey, typename GraphData::const_iterator> {
 public:
  using Key = GraphKey;
  using Iterator = typename GraphData::const_iterator;

  explicit ClusterGraph(const GraphData& clusterData)
    : stk::util::TreeTraverser<Key, Iterator>()
    , m_clusterData(clusterData)
  {
  }

  Iterator begin() const override { return m_clusterData.begin(); }

  Iterator end() const override { return m_clusterData.end(); }

  const Key& get_node(Iterator iter) const override { return iter->first; }

  size_t size() const override { return m_clusterData.size(); }

  bool node_exists(const Key& node) const override { return m_clusterData.find(node) != m_clusterData.end(); }

  const std::vector<Key>& children(const Key& node) const override
  {
    const auto& iter = m_clusterData.find(node);
    if(iter != m_clusterData.end()) {
      return iter->second;
    }

    return m_emptyChildren;
  }

 private:
  ClusterGraph() = delete;
  ClusterGraph(const ClusterGraph&) = delete;

  const GraphData& m_clusterData;
  const std::vector<Key> m_emptyChildren{};
};

class NeighborClusterConnectivityHelper
{
public:
  NeighborClusterConnectivityHelper(stk::mesh::BulkData& bulk)
  : m_bulk(bulk), m_comm(bulk.parallel())
  {

  }

  void fill_neighbor_cluster_connectivity(std::vector<ElementCluster>& elementClusters, std::vector<ParallelClusterInfo>& pcInfoVec)
  {
    std::unordered_map<stk::mesh::EntityId, int64_t> entityToClusterMap;
    for(size_t i=0; i<elementClusters.size(); i++) {
      for(stk::mesh::Entity entry : elementClusters[i].elements) {
        entityToClusterMap[m_bulk.identifier(entry)] = i;
      }
    }

    pack_and_communicate_local_cluster_info(pcInfoVec);
    unpack_and_update_cluster_connectivity(pcInfoVec, entityToClusterMap);
  }

private:
  void pack_data_for_clusters(std::vector<ParallelClusterInfo>& pcInfoVec)
  {
    for(const auto& pcInfo : pcInfoVec)
    {
      for(const auto& item : pcInfo.sendData) {
        const int procRank = item.first;
        const stk::mesh::EntityIdVector &procElemIds = item.second;

        m_comm.send_buffer(procRank).pack<int64_t>(pcInfo.clusterIndex);
        m_comm.send_buffer(procRank).pack<stk::mesh::EntityId>(pcInfo.minElemId);

        m_comm.send_buffer(procRank).pack<unsigned>(procElemIds.size());
        for(stk::mesh::EntityId procElemId : procElemIds) {
          m_comm.send_buffer(procRank).pack<stk::mesh::EntityId>(procElemId);
        }
      }
    }
  }

  void pack_and_communicate_local_cluster_info(std::vector<ParallelClusterInfo>& pcInfoVec)
  {
    stk::pack_and_communicate(m_comm, [this, &pcInfoVec]()
                                    {
                                      pack_data_for_clusters(pcInfoVec);
                                    });
  }

  void unpack_data_for_clusters(int rank, std::vector<ParallelClusterInfo>& pcInfoVec,
                                std::unordered_map<stk::mesh::EntityId, int64_t>& entityToClusterMap)
  {
    int64_t remoteClusterIndex;
    stk::mesh::EntityId remoteClusterMinElemId;

    m_comm.recv_buffer(rank).unpack<int64_t>(remoteClusterIndex);
    m_comm.recv_buffer(rank).unpack<stk::mesh::EntityId>(remoteClusterMinElemId);

    unsigned numLocalElemIds;
    m_comm.recv_buffer(rank).unpack<unsigned>(numLocalElemIds);

    stk::mesh::EntityId localElemId;
    for(unsigned i=0; i<numLocalElemIds; i++) {
      m_comm.recv_buffer(rank).unpack<stk::mesh::EntityId>(localElemId);

#if defined(DEBUG)
      debugStream_ << "\tUnpacked local element: " << localElemId
                   << " from proc: " << rank
                   << " for remote cluster index: " << remoteClusterIndex
                   << " and remote min elem id: " << remoteClusterMinElemId
                   << std::endl;
#endif

      auto iter = entityToClusterMap.find(localElemId);
      STK_ThrowRequireMsg(iter != entityToClusterMap.end(),
                      "Could not find local element " << localElemId << " from proc: " << rank
                      << " for remote cluster index: " << remoteClusterIndex);

      int64_t localClusterIndex = iter->second;
      ParallelClusterRecvInfo entry(remoteClusterIndex, remoteClusterMinElemId);
      pcInfoVec[localClusterIndex].recvData[rank].push_back(entry);
    }
  }

  void unpack_and_update_cluster_connectivity(std::vector<ParallelClusterInfo>& pcInfoVec,
                                  std::unordered_map<stk::mesh::EntityId, int64_t>& entityToClusterMap)
  {
    stk::unpack_communications(m_comm, [this, &pcInfoVec, &entityToClusterMap](int rank)
                                     {
                                       unpack_data_for_clusters(rank, pcInfoVec, entityToClusterMap);
                                     });
  }

  NeighborClusterConnectivityHelper() = delete;
  NeighborClusterConnectivityHelper(const NeighborClusterConnectivityHelper&) = delete;

  const stk::mesh::BulkData& m_bulk;
  stk::CommSparse m_comm;
};

class GlobalClusterConnectivityHelper
{
public:
  GlobalClusterConnectivityHelper(stk::mesh::BulkData& bulk)
  : m_bulk(bulk), m_comm(bulk.parallel())
  {

  }

  void fill_global_cluster_graph_edges(const std::vector<ParallelClusterInfo>& pcInfoVec,
                                       std::vector<ClusterGraphEdge>& globalClusterEdges)
  {
    std::vector<ClusterGraphEdge> localClusterEdges;

    for(const ParallelClusterInfo& pcInfo : pcInfoVec) {
      ClusterGraphNode localClusterNode(m_bulk.parallel_rank(), pcInfo.clusterIndex, pcInfo.minElemId);
      std::set<ClusterGraphNode> remoteClusterNodesSet;

      for(const auto& item : pcInfo.recvData) {
        int proc = item.first;
        for(const ParallelClusterRecvInfo& recvInfo : item.second) {
          ClusterGraphNode remoteClusterNode(proc, recvInfo.clusterIndex, recvInfo.clusterMinElemId);
          remoteClusterNodesSet.insert(remoteClusterNode);
        }
      }

      for(const ClusterGraphNode& remoteClusterNode : remoteClusterNodesSet) {
        localClusterEdges.emplace_back(localClusterNode, remoteClusterNode);
#if defined(DEBUG)
      debugStream_ << "\t\t" << localClusterEdges.back() << std::endl;
#endif
      }
    }

    pack_and_communicate_cluster_edges(localClusterEdges);
    unpack_remote_cluster_edges(globalClusterEdges);
  }

private:
  void pack_data_for_cluster_edges(std::vector<ClusterGraphEdge>& clusterEdges)
  {
    for(int i=0; i<m_bulk.parallel_size(); i++) {
      if(i == m_bulk.parallel_rank()) continue;

      for(const auto& clusterEdge : clusterEdges)
      {
        m_comm.send_buffer(i).pack<int>(clusterEdge.node1.proc);
        m_comm.send_buffer(i).pack<int64_t>(clusterEdge.node1.clusterIndex);
        m_comm.send_buffer(i).pack<stk::mesh::EntityId>(clusterEdge.node1.clusterMinElemId);

        m_comm.send_buffer(i).pack<int>(clusterEdge.node2.proc);
        m_comm.send_buffer(i).pack<int64_t>(clusterEdge.node2.clusterIndex);
        m_comm.send_buffer(i).pack<stk::mesh::EntityId>(clusterEdge.node2.clusterMinElemId);
      }
    }
  }

  void pack_and_communicate_cluster_edges(std::vector<ClusterGraphEdge>& clusterEdges)
  {
    stk::pack_and_communicate(m_comm, [this, &clusterEdges]()
                                    {
                                      pack_data_for_cluster_edges(clusterEdges);
                                    });
  }

  void unpack_data_for_cluster_edges(int rank, std::vector<ClusterGraphEdge>& remoteClusterEdges)
  {
    ClusterGraphNode node1;
    ClusterGraphNode node2;

    m_comm.recv_buffer(rank).unpack<int>(node1.proc);
    m_comm.recv_buffer(rank).unpack<int64_t>(node1.clusterIndex);
    m_comm.recv_buffer(rank).unpack<stk::mesh::EntityId>(node1.clusterMinElemId);

    m_comm.recv_buffer(rank).unpack<int>(node2.proc);
    m_comm.recv_buffer(rank).unpack<int64_t>(node2.clusterIndex);
    m_comm.recv_buffer(rank).unpack<stk::mesh::EntityId>(node2.clusterMinElemId);

    remoteClusterEdges.emplace_back(node1, node2);

#if defined(DEBUG)
      debugStream_ << "\tUnpacked remote cluster edge: " << remoteClusterEdges.back()
                   << " from proc: " << rank << std::endl;
#endif
  }

  void unpack_remote_cluster_edges(std::vector<ClusterGraphEdge>& remoteClusterEdges)
  {
    stk::unpack_communications(m_comm, [this, &remoteClusterEdges](int rank)
                                     {
                                       unpack_data_for_cluster_edges(rank, remoteClusterEdges);
                                     });
  }

  GlobalClusterConnectivityHelper() = delete;
  GlobalClusterConnectivityHelper(const GlobalClusterConnectivityHelper&) = delete;

  const stk::mesh::BulkData& m_bulk;
  stk::CommSparse m_comm;
};

class MeshTraverser
{
public:
  static constexpr size_t INVALID_SEED = std::numeric_limits<size_t>::max();

  MeshTraverser(stk::mesh::BulkData& bulk) : m_bulk(bulk), m_currentSeed(0), m_currentNumVisitedElements(0)
  {
    m_bulk.initialize_face_adjacent_element_graph();
    create_local_element_clusters();
    assign_global_cluster_ids();

#if defined(DEBUG)
    std::cout << "P[" << m_bulk.parallel_rank() << "]\n" << debugStream_.str();
#endif
  }

  const std::vector<ElementCluster>& get_clusters() const { return m_elementClusters; }

private:
  void fill_parallel_cluster_send_info(ParallelClusterInfo& pcInfo)
  {
    if(pcInfo.clusterIndex < 0) return;

    stk::mesh::ElemElemGraph& graph = m_bulk.get_face_adjacent_element_graph();
    const stk::mesh::EntityVector& elements = m_elementClusters[pcInfo.clusterIndex].elements;

    for(stk::mesh::Entity entry : elements) {
      stk::mesh::impl::LocalId localId = graph.get_local_element_id(entry);
      for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(localId)) {
        if(!graphEdge.is_elem2_local()) {
          stk::mesh::EntityId remoteElemId = graph.convert_negative_local_id_to_global_id(graphEdge.elem2());
          int remoteProc =  graph.get_parallel_info_for_graph_edge(graphEdge).get_proc_rank_of_neighbor();
          pcInfo.sendData[remoteProc].push_back(remoteElemId);

  #if defined(DEBUG)
          debugStream_ << "\tFound remote boundary element: " << remoteElemId
                       << " on proc: " << remoteProc
                       << " for element: " << m_bulk.entity_key(entry)
                       << std::endl;
  #endif
        }
      }

      pcInfo.minElemId = std::min(pcInfo.minElemId, m_bulk.identifier(entry));
    }
  }

  void update_global_cluster_ids_from_graph(const std::vector<ParallelClusterInfo>& pcInfoVec,
                                            ClusterGraph<ClusterGraphNode, ClusterGraphData>& globalClusterGraph)
  {
    for(const ParallelClusterInfo& pcInfo : pcInfoVec) {
      ClusterGraphNode localClusterNode(m_bulk.parallel_rank(), pcInfo.clusterIndex, pcInfo.minElemId);
      std::vector<ClusterGraphNode> clusterList = globalClusterGraph.get_forward_traversal_list(localClusterNode);

      stk::mesh::EntityId globalId = m_elementClusters[pcInfo.clusterIndex].globalId;
      for(const auto& clusterNode : clusterList) {
        globalId = std::min(globalId, clusterNode.clusterMinElemId);
      }

      m_elementClusters[pcInfo.clusterIndex].globalId = globalId;
    }
  }

  void fill_global_cluster_graph_edges(const std::vector<ParallelClusterInfo>& pcInfoVec,
                                       std::vector<ClusterGraphEdge>& globalClusterEdges)
  {
    GlobalClusterConnectivityHelper helper(m_bulk);
    helper.fill_global_cluster_graph_edges(pcInfoVec, globalClusterEdges);
  }

  void fill_cluster_graph_data(const std::vector<ClusterGraphEdge>& clusterEdges,
                               ClusterGraphData& clusterGraphData)
  {
    for(const ClusterGraphEdge& clusterEdge : clusterEdges) {
      std::vector<ClusterGraphNode>& node1Entries = clusterGraphData[clusterEdge.node1];
      std::vector<ClusterGraphNode>& node2Entries = clusterGraphData[clusterEdge.node2];

      stk::util::insert_keep_sorted_and_unique(clusterEdge.node2, node1Entries);
      stk::util::insert_keep_sorted_and_unique(clusterEdge.node1, node2Entries);
    }
  }

  void update_global_cluster_ids(const std::vector<ParallelClusterInfo>& pcInfoVec)
  {
    std::vector<ClusterGraphEdge> globalClusterEdges;
    fill_global_cluster_graph_edges(pcInfoVec, globalClusterEdges);

    ClusterGraphData globalClusterGraphData;
    fill_cluster_graph_data(globalClusterEdges, globalClusterGraphData);

    ClusterGraph<ClusterGraphNode, ClusterGraphData> globalClusterGraph(globalClusterGraphData);
    update_global_cluster_ids_from_graph(pcInfoVec, globalClusterGraph);
  }

  void fill_neighbor_cluster_connectivity(std::vector<ParallelClusterInfo>& pcInfoVec)
  {
    NeighborClusterConnectivityHelper helper(m_bulk);
    helper.fill_neighbor_cluster_connectivity(m_elementClusters, pcInfoVec);
  }

  void assign_global_cluster_ids()
  {
    std::vector<ParallelClusterInfo> pcInfoVec;
    pcInfoVec.reserve(m_elementClusters.size());

    for(size_t i=0; i<m_elementClusters.size(); i++) {
      pcInfoVec.emplace_back(i);

      ParallelClusterInfo& currentInfo = pcInfoVec.back();
      fill_parallel_cluster_send_info(currentInfo);

      m_elementClusters[i].globalId = currentInfo.minElemId;
    }

    fill_neighbor_cluster_connectivity(pcInfoVec);
    update_global_cluster_ids(pcInfoVec);
  }

  void create_local_element_clusters()
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(m_bulk.mesh_meta_data().locally_owned_part(), m_bulk.buckets(stk::topology::ELEM_RANK), elements);

    for (stk::mesh::Entity element : elements) {
      m_visitedElements[element] = false;
    }

    size_t newSeed = get_new_seed(elements);

    while(newSeed != INVALID_SEED) {
      add_sorted_element_cluster_from_seed(newSeed, elements);
      newSeed = get_new_seed(elements);
    }
  }

  size_t get_new_seed(const stk::mesh::EntityVector& elements)
  {
    for(size_t i=m_currentSeed; i<elements.size(); i++) {
      if(not m_visitedElements[elements[i]]) {
        return i;
      }
    }

    return INVALID_SEED;
  }

  void update_cluster_queue_with_entry(stk::mesh::EntityVector& queue, stk::mesh::Entity entry)
  {
    if(m_visitedElements[entry] || !m_bulk.bucket(entry).owned()) return;

    stk::mesh::ElemElemGraph& graph = m_bulk.get_face_adjacent_element_graph();
    stk::mesh::impl::LocalId localId = graph.get_local_element_id(entry);
    for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(localId)) {
      if(graphEdge.is_elem2_local()) {
        stk::mesh::Entity neighbor = graph.get_entity(graphEdge.elem2());

        if(not m_visitedElements[neighbor]) {
          queue.push_back(neighbor);

#if defined(DEBUG)
          debugStream_ << "\t" << m_bulk.entity_key(neighbor)
                       << " [" << m_bulk.bucket(neighbor).topology() << "]" << std::endl;
#endif
        }
      }
    }

    m_visitedElements[entry] = true;
  }

  void add_sorted_element_cluster_from_seed(const size_t newSeed, const stk::mesh::EntityVector& elements)
  {
    if(m_visitedElements[elements[newSeed]]) {
      return;
    }

    m_currentSeed = newSeed;

    m_elementClusters.emplace_back();
    ElementCluster& cluster = m_elementClusters.back();

    cluster.elements.reserve(elements.size() - m_currentNumVisitedElements);
    cluster.elements.push_back(elements[m_currentSeed]);

#if defined(DEBUG)
    debugStream_ << "Cluster: " << m_elementClusters.size()-1 << std::endl;
    debugStream_ << "\t" << m_bulk.entity_key(cluster.elements.back())
                 << " [" << m_bulk.bucket(cluster.elements.back()).topology() << "]" << std::endl;
#endif

    for(size_t i=0; i<cluster.elements.size(); i++) {
      update_cluster_queue_with_entry(cluster.elements, cluster.elements[i]);
    }

    cluster.elements.resize(cluster.elements.size());
    std::sort(cluster.elements.begin(), cluster.elements.end(), stk::mesh::EntityLess(m_bulk));

    m_currentNumVisitedElements += cluster.elements.size();

#if defined(DEBUG)
    debugStream_ << std::endl;
#endif
  }

  MeshTraverser() = delete;
  MeshTraverser(const MeshTraverser&) = delete;

  stk::mesh::BulkData& m_bulk;
  size_t m_currentSeed = 0;
  size_t m_currentNumVisitedElements = 0;
  std::unordered_map<stk::mesh::Entity, bool> m_visitedElements;
  std::vector<ElementCluster> m_elementClusters;

#if defined(DEBUG)
  std::ostringstream debugStream_;
#endif
};

stk::mesh::EntityIdVector convert_entities_to_sorted_ids(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& entities)
{
  stk::mesh::EntityIdVector ids;
  for(stk::mesh::Entity element : entities) {
    ids.push_back(bulk.identifier(element));
  }
  std::sort(ids.begin(), ids.end());

  return ids;
}

void test_clusters(const stk::mesh::BulkData& bulk, const std::vector<ElementCluster>& clusters, const ElementClusterVector& goldIds)
{
  EXPECT_EQ(goldIds.size(), clusters.size());

  for(size_t i=0; i<goldIds.size(); i++) {
    EXPECT_EQ(goldIds[i].globalId, clusters[i].globalId);

    const stk::mesh::EntityIdVector& goldClusterIds = goldIds[i].elements;
    const stk::mesh::EntityVector& clusterElements = clusters[i].elements;

    stk::mesh::EntityIdVector clusterIds = convert_entities_to_sorted_ids(bulk, clusterElements);

    for(stk::mesh::EntityId goldClusterId : goldClusterIds) {
      EXPECT_TRUE(std::binary_search(clusterIds.begin(), clusterIds.end(), goldClusterId));
    }
  }
}

void run_traverser(const std::string& meshDesc, const ElementClusterVector& goldIds, const unsigned dim = 3)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(dim,MPI_COMM_WORLD);
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  MeshTraverser bfs(*bulk);
  test_clusters(*bulk, bfs.get_clusters(), goldIds);
}

TEST(Mechanism, twoHexConnected_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12";

  ElementIdCluster cluster(1, stk::mesh::EntityIdVector{1,2});

  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, twoHexConnected_Parallel)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "1,2,HEX_8,5,6,7,8,9,10,11,12";

  unsigned myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  ElementIdCluster cluster{1, stk::mesh::EntityIdVector{myProc+1}};

  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, threeHexWithTwoConnectedAcrossProcessorBoundary_Parallel)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,3,HEX_8,13,14,15,16,17,18,19,20";

  ElementClusterVector goldIds;
  if(stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    ElementIdCluster cluster1{1, stk::mesh::EntityIdVector{1}};
    ElementIdCluster cluster2{2, stk::mesh::EntityIdVector{2}};
    goldIds = ElementClusterVector{cluster1, cluster2};
  } else {
    ElementIdCluster cluster{2, stk::mesh::EntityIdVector{3}};
    goldIds = ElementClusterVector{cluster};
  }
  run_traverser(meshDesc, goldIds);
}

TEST(Mechanism, twoHexDisconnected_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";

  ElementIdCluster cluster1(1, stk::mesh::EntityIdVector{1});
  ElementIdCluster cluster2(2, stk::mesh::EntityIdVector{2});
  run_traverser(meshDesc, ElementClusterVector{cluster1, cluster2});
}

TEST(Mechanism, threeHexConnected_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
      "0,3,HEX_8,9,10,11,12,13,14,15,16";

  ElementIdCluster cluster(1, stk::mesh::EntityIdVector{1,2,3});
  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, threeHexConnected_Parallel)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "1,2,HEX_8,5,6,7,8,9,10,11,12\n"
      "2,3,HEX_8,9,10,11,12,13,14,15,16";

  unsigned myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  ElementIdCluster cluster{1, stk::mesh::EntityIdVector{myProc+1}};

  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, fourHexConnected_Parallel)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 4) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "1,2,HEX_8,5,6,7,8,9,10,11,12\n"
      "2,3,HEX_8,9,10,11,12,13,14,15,16\n"
      "3,4,HEX_8,13,14,15,16,17,18,19,20";

  unsigned myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);
  ElementIdCluster cluster{1, stk::mesh::EntityIdVector{myProc+1}};

  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, threeHexDisconnected_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
      "0,3,HEX_8,17,18,19,20,21,22,23,24";

  ElementIdCluster cluster1(1, stk::mesh::EntityIdVector{1});
  ElementIdCluster cluster2(2, stk::mesh::EntityIdVector{2});
  ElementIdCluster cluster3(3, stk::mesh::EntityIdVector{3});
  run_traverser(meshDesc, ElementClusterVector{cluster1, cluster2, cluster3});
}

TEST(Mechanism, connectedPyramidAndHexAndTet_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,PYRAMID_5,1,2,3,4,5\n"
      "0,2,HEX_8,1,4,3,2,6,9,8,7\n"
      "0,3,TET_4,2,3,5,10";

  ElementIdCluster cluster(1, stk::mesh::EntityIdVector{1,2,3});
  run_traverser(meshDesc, ElementClusterVector{cluster});
}

TEST(Mechanism, connectedPyramidAndHexAndDisconnectedTet_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,PYRAMID_5,1,2,3,4,5\n"
      "0,2,HEX_8,1,4,3,2,6,9,8,7\n"
      "0,3,TET_4,10,11,12,13";

  ElementIdCluster cluster1(1, stk::mesh::EntityIdVector{1,2});
  ElementIdCluster cluster2(3, stk::mesh::EntityIdVector{3});
  run_traverser(meshDesc, ElementClusterVector{cluster1, cluster2});
}

TEST(Mechanism, disconnectedPyramidAndHexAndTet_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc =
      "0,1,PYRAMID_5,1,2,3,4,5\n"
      "0,2,HEX_8,6,7,8,9,10,11,12,13\n"
      "0,3,TET_4,14,15,16,17";

  ElementIdCluster cluster1(1, stk::mesh::EntityIdVector{1});
  ElementIdCluster cluster2(2, stk::mesh::EntityIdVector{2});
  ElementIdCluster cluster3(3, stk::mesh::EntityIdVector{3});
  run_traverser(meshDesc, ElementClusterVector{cluster1, cluster2, cluster3});
}

TEST(Mechanism, fourQuadBowtie_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,5,6,7,4,block_1\n"
                         "0,3,QUAD_4_2D,4,8,9,10,block_1\n"
                         "0,4,QUAD_4_2D,4,11,12,13,block_1";

  ElementIdCluster cluster1(1, stk::mesh::EntityIdVector{1});
  ElementIdCluster cluster2(2, stk::mesh::EntityIdVector{2});
  ElementIdCluster cluster3(3, stk::mesh::EntityIdVector{3});
  ElementIdCluster cluster4(4, stk::mesh::EntityIdVector{4});
  run_traverser(meshDesc, ElementClusterVector{cluster1, cluster2, cluster3, cluster4}, 2u);
}

TEST(Mechanism, fourQuadPacman_Serial)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_1\n"
                         "0,4,QUAD_4_2D,3,4,9,10,block_1";

  ElementIdCluster cluster(1, stk::mesh::EntityIdVector{1,2,3,4});
  run_traverser(meshDesc, ElementClusterVector{cluster}, 2u);
}
