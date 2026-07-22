#include "NodeBalancer.hpp"
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include "stk_balance/internal/privateDeclarations.hpp"

namespace stk {
namespace balance {
namespace internal {

NodeBalancer::NodeBalancer(stk::mesh::BulkData& bulk)
  : m_bulkData(bulk),
    m_metaData(bulk.mesh_meta_data())
{}

bool
NodeBalancer::balance_node_entities(const double targetLoadBalance,
                                    const unsigned maxIterations)
{
  bool changedNodeOwnership = false;
  if (m_bulkData.parallel_size() == 1) return changedNodeOwnership;

  for (unsigned bal = 0; bal < maxIterations; ++bal) {
    int numLocallyOwnedNodes = 0;
    double loadFactor = 0;

    getGlobalLoadImbalance(loadFactor, numLocallyOwnedNodes);

    const bool converged = (loadFactor <= targetLoadBalance);
    if (converged) break;

    std::set<int> neighborProcessors;
    std::map<stk::mesh::Entity, std::vector<int>> interfaceNodesAndProcessors;

    getInterfaceDescription(neighborProcessors,
                            interfaceNodesAndProcessors);

    std::map<int, int> numLocallyOwnedByRank;
    exchangeLocalSizes(neighborProcessors, numLocallyOwnedNodes, numLocallyOwnedByRank);

    changedNodeOwnership = true;
    changeOwnersOfNodes(interfaceNodesAndProcessors, numLocallyOwnedByRank, numLocallyOwnedNodes);
  }
  return changedNodeOwnership;
}

void
NodeBalancer::getInterfaceDescription(std::set<int>& neighborProcessors,
                                      std::map<stk::mesh::Entity,
                                      std::vector<int> >& interfaceNodesAndProcessors)
{
  stk::mesh::Selector sharedSelector(m_metaData.globally_shared_part());
  const stk::mesh::BucketVector& buckets = m_bulkData.get_buckets(stk::topology::NODE_RANK, sharedSelector);

  for (auto ib = buckets.begin(); ib != buckets.end(); ++ib) {

    const stk::mesh::Bucket& b = **ib;
    const size_t nnodes = b.size();
    for (size_t n = 0; n < nnodes; ++n) {
      stk::mesh::Entity node = b[n];
      std::vector<int> shared_procs;

      m_bulkData.comm_shared_procs(m_bulkData.entity_key(node), shared_procs);
      neighborProcessors.insert(shared_procs.begin(), shared_procs.end());
      std::pair<stk::mesh::Entity, std::vector<int>> item(node, shared_procs);
      interfaceNodesAndProcessors.insert(item);
    }
  }
}

void
NodeBalancer::getGlobalLoadImbalance(double &loadFactor, int& numLocallyOwnedNodes)
{
  stk::mesh::Selector localSelector = m_metaData.locally_owned_part();
  numLocallyOwnedNodes = stk::mesh::count_entities(m_bulkData, stk::topology::NODE_RANK, localSelector);

  int maxLocallyOwned = 0;
  int minLocallyOwned = 0;
  stk::all_reduce_max(m_bulkData.parallel(), &numLocallyOwnedNodes, &maxLocallyOwned, 1);
  stk::all_reduce_min(m_bulkData.parallel(), &numLocallyOwnedNodes, &minLocallyOwned, 1);
  loadFactor = double(maxLocallyOwned) / double(minLocallyOwned);

  internal::logMessage(m_bulkData.parallel(),
                       "Balancing locally-owned nodes: min = " + std::to_string(minLocallyOwned) +
                       ", max = " + std::to_string(maxLocallyOwned) +
                       " (loadFactor = " + std::to_string(loadFactor) + ")");
}

void
NodeBalancer::exchangeLocalSizes(const std::set<int>& neighborProcessors,
                                 int& numLocallyOwnedNodes,
                                 std::map<int, int> &numLocallyOwnedByRank)
{
  size_t numCommunications = neighborProcessors.size();
  if (numCommunications == 0u) return;

  std::vector<int> recvBuffer(numCommunications);
  std::vector<MPI_Request> receiveRequests(numCommunications);
  std::vector<MPI_Request> sendRequests(numCommunications);

  int bufferCounter = 0;
  for (int p : neighborProcessors) {
    MPI_Irecv(&recvBuffer[bufferCounter], 1, MPI_INT, p,
              MPI_ANY_TAG, m_bulkData.parallel(), &receiveRequests[bufferCounter]);
    ++bufferCounter;
  }

  bufferCounter = 0;
  for (int p : neighborProcessors) {
    MPI_Isend(&numLocallyOwnedNodes, 1, MPI_INT, p, 0, m_bulkData.parallel(), &sendRequests[bufferCounter]);
    ++bufferCounter;
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), receiveStati.data());

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), sendRequests.data(), sendStati.data());

  int i = 0;
  for (int p : neighborProcessors) {
    int n = recvBuffer[i];
    numLocallyOwnedByRank.insert(std::pair<int, int>(p, n));
    ++i;
  }
}

void
NodeBalancer::changeOwnersOfNodes(const std::map<stk::mesh::Entity, std::vector<int> >& interfaceNodesAndProcessors,
                                  std::map<int, int>& numLocallyOwnedByRank,
                                  int numLocallyOwnedNodes)
{
  int myRank = m_bulkData.parallel_rank();

  stk::mesh::EntityProcVec nodesToMove;

  for (const auto & nodeProcs : interfaceNodesAndProcessors) {
    stk::mesh::Entity node = nodeProcs.first;
    const bool isOwned = m_bulkData.bucket(node).owned();

    if (!isOwned) continue;

    const auto & procs = nodeProcs.second;
    double maxLoad = 1;
    int destination = myRank;

    for (auto p : procs) {
      double numOwnedByP = numLocallyOwnedByRank[p];
      double load = double(numLocallyOwnedNodes) / numOwnedByP;
      if (load >  maxLoad) {  //could test for max rank
        maxLoad = load;
        destination = p;
      }
    }

    if (destination != myRank) {
      std::pair<stk::mesh::Entity, int> item(node, destination);
      nodesToMove.push_back(item);
      numLocallyOwnedNodes--;
      numLocallyOwnedByRank[destination] += 1;
    }

  }

  m_bulkData.change_entity_owner(nodesToMove);
}


}
}
}
