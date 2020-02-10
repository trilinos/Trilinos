#ifndef NODEBALANCER_HPP
#define NODEBALANCER_HPP

#include <set>
#include <map>
#include <vector>

#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace balance {
namespace internal {

class NodeBalancer {

public:
  NodeBalancer(stk::mesh::BulkData& bulk);

  bool balance_node_entities(const double targetLoadBalance,
                             const unsigned maxIterations);

private:
  void getInterfaceDescription(std::set<int>& neighborProcessors,
                               std::map<stk::mesh::Entity, std::vector<int>>& interfaceNodesAndProcessors);

  void getGlobalLoadImbalance(double &loadFactor,
                              int& numLocallyOwnedNodes);

  void exchangeLocalSizes(const std::set<int>& neighborProcessors,
                          int& numLocallyOwnedNodes,
                          std::map<int, int>& numLocallyOwnedByRank);

  void changeOwnersOfNodes(const std::map<stk::mesh::Entity, std::vector<int> >& interfaceNodesAndProcessors,
                           std::map<int, int>& numLocallyOwnedByRank,
                           int numLocallyOwnedNodes);

  stk::mesh::BulkData & m_bulkData;
  const stk::mesh::MetaData & m_metaData;
};

}
}
}

#endif // NODEBALANCER_HPP
