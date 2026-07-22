
#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_SECONDORDERMESHSPECS_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_SECONDORDERMESHSPECS_HPP_

#include <Akri_TopologyData.hpp>
#include <algorithm>
#include <array>
#include <map>

#include <Akri_MeshSpecs.hpp>

namespace krino {

template<typename BASEMESHSPEC, stk::topology::topology_t TOPO>
struct HigherOrderSimplexMeshSpec
{
  static constexpr stk::topology::topology_t TOPOLOGY = TOPO;
  using TopoData = TopologyData<TOPO>;
  static constexpr unsigned NPE = TopoData::num_nodes();
  static constexpr unsigned DIM = TopoData::spatial_dimension();
  using MidNodeMap = std::map<std::pair<unsigned, unsigned>, unsigned>;

  HigherOrderSimplexMeshSpec()
  {
    const auto & baseNodeLocs = baseMeshSpec.nodeLocs;

    for (auto loc : baseNodeLocs)
      nodeLocs.push_back(loc);

    std::array<unsigned,2> baseEdgeNodes;
    std::array<unsigned, NPE> elemConn;
    MidNodeMap midNodeMap;
    unsigned midNodeId = baseMeshSpec.nodeLocs.size();
    unsigned baseNPE = elemTopology.base().num_nodes();
    for (auto baseElemConn : baseMeshSpec.allElementConn)
    {
      for (unsigned i=0; i<baseNPE; ++i)
        elemConn[i] = baseElemConn[i];
      for (unsigned iEdge=0; iEdge<elemTopology.num_edges(); ++iEdge)
      {
        elemTopology.base().edge_nodes(baseElemConn.data(), iEdge, baseEdgeNodes.data());
        const auto edgeKey = std::minmax(baseEdgeNodes[0],baseEdgeNodes[1]);
        const auto iter = midNodeMap.find(edgeKey);
        unsigned & elemEdgeMidNode = elemConn[baseNPE+iEdge]; // assumes simplex
        if (iter == midNodeMap.end())
        {
          elemEdgeMidNode = midNodeId++;
          midNodeMap[edgeKey] = elemConn[baseNPE+iEdge];
          nodeLocs.push_back(0.5*(baseNodeLocs[baseEdgeNodes[0]] + baseNodeLocs[baseEdgeNodes[1]]));
        }
        else
        {
          elemEdgeMidNode = iter->second;
        }
      }
      allElementConn.push_back(elemConn);
    }
  }

  stk::topology elemTopology{TOPO};
  std::vector<stk::math::Vec<double,DIM>> nodeLocs;
  std::vector<std::array<unsigned, NPE>> allElementConn;

  BASEMESHSPEC baseMeshSpec;
};

using RegularTri6 = HigherOrderSimplexMeshSpec<RegularTri, stk::topology::TRIANGLE_6_2D>;
using RightTri6SurroundedByEdgeTris = HigherOrderSimplexMeshSpec<RightTriSurroundedByEdgeTris, stk::topology::TRIANGLE_6_2D>;
using RegularTet10 = HigherOrderSimplexMeshSpec<RegularTet, stk::topology::TETRAHEDRON_10>;
using UMRRegularTet10 = HigherOrderSimplexMeshSpec<UMRRegularTet, stk::topology::TETRAHEDRON_10>;

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_SECONDORDERMESHSPECS_HPP_ */
