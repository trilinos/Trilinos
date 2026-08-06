// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "stk_tools/mesh_tools/PairwiseSideInfo.hpp"
#include "stk_tools/mesh_tools/DetectHingesImpl.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk {
namespace tools {
namespace impl {

bool are_equal(const stk::mesh::Entity* nodesPtr, const stk::mesh::EntityVector& nodesVec)
{
  for(stk::mesh::Entity node : nodesVec) {
    if (node != *nodesPtr) {
      return false;
    }
    ++nodesPtr;
  }
  return true;
}

class SideFinder
{
public:
  SideFinder(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
    : m_bulk(bulkData), m_elem(element)
  { 
    const unsigned numNodes = m_bulk.num_nodes(m_elem);
    if (numNodes > MAX_NUM_NODES) {
      m_heapScratchSpace.resize(numNodes);
      m_scratchSpace = m_heapScratchSpace.data();
    }
    else {
      m_scratchSpace = m_stackScratchSpace;
    }
  }
  
  virtual ~SideFinder() {}
  
  bool put_nodes_in_side_order(stk::mesh::EntityVector& nodes)
  { 
    stk::topology elementTopo = m_bulk.bucket(m_elem).topology();
    stk::mesh::EntityRank subRank = m_bulk.mesh_meta_data().side_rank();
    
    const unsigned numSubTopo = elementTopo.num_sub_topology(subRank);
    const unsigned numNodes = nodes.size();
    const stk::mesh::Entity* elemNodes = m_bulk.begin_nodes(m_elem);
    
    for(unsigned i = 0; i < numSubTopo; i++) {
      stk::topology subTopo = elementTopo.sub_topology(subRank, i);
      if(numNodes != subTopo.num_nodes()) {
        continue;
      }
      
      elementTopo.sub_topology_nodes(elemNodes, subRank, i, m_scratchSpace);
      std::sort(m_scratchSpace, m_scratchSpace+subTopo.num_nodes());
      
      if(are_equal(m_scratchSpace, nodes)) {
        elementTopo.sub_topology_nodes(elemNodes, subRank, i, nodes.data());
        return true;
      }
    }
    return false;
  }

private:
  const stk::mesh::BulkData& m_bulk;
  stk::mesh::Entity m_elem; 
  static constexpr unsigned MAX_NUM_NODES = 32;
  stk::mesh::Entity m_stackScratchSpace[MAX_NUM_NODES];
  stk::mesh::EntityVector m_heapScratchSpace;
  stk::mesh::Entity* m_scratchSpace;
};

std::pair<stk::mesh::EntityVector,bool> get_pairwise_common_nodes(const stk::mesh::BulkData& bulk,
                                                                  stk::mesh::Entity elem1,
                                                                  stk::mesh::Entity elem2)
{
  std::pair<stk::mesh::EntityVector,bool> result;
  stk::mesh::EntityVector& commonNodes = result.first;
  commonNodes.assign(bulk.begin_nodes(elem1), bulk.end_nodes(elem1));

  const stk::mesh::Entity* elem2Nodes = bulk.begin_nodes(elem2);
  const unsigned numElem2Nodes = bulk.num_nodes(elem2);
  unsigned numIntersect = 0;
  const unsigned numElem1Nodes = commonNodes.size();
  for(unsigned n=0; n<numElem1Nodes; ++n) {
    for(unsigned m=0; m<numElem2Nodes; ++m) {
      if (commonNodes[n] == elem2Nodes[m]) {
        if (n > numIntersect) {
          commonNodes[numIntersect] = commonNodes[n];
        }
        ++numIntersect;
        break;
      }
    }
  }
  commonNodes.resize(numIntersect);
  std::sort(commonNodes.begin(), commonNodes.end());

  SideFinder sideFinder(bulk, elem1);
  bool commonNodesAreSideNodes = sideFinder.put_nodes_in_side_order(commonNodes);
  result.second = commonNodesAreSideNodes;

  return result;
}

}}} // namespace stk::tools::impl

