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

#include "stk_util/util/GraphCycleDetector.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include "stk_util/util/SortAndUnique.hpp"  // for insert_keep_sorted_and_unique
#include <memory>                           // for allocator_traits<>::value_type

namespace stk {
namespace tools {
namespace impl {

void GraphCycleDetector::add_edge(unsigned node1Id, unsigned node2Id)
{
  STK_ThrowRequire(node1Id < numNodes && node2Id < numNodes);
  nodeGraph[node1Id].push_back(node2Id);
  nodeGraph[node2Id].push_back(node1Id);
}

void GraphCycleDetector::get_cycles_helper(unsigned visitNodeId, unsigned parentNodeId)
{
  if(visited[visitNodeId] == 2) { return; }

  if(visited[visitNodeId] == 1) {
    unsigned currentNodeId = parentNodeId;
    stk::util::insert_keep_sorted_and_unique(currentNodeId, nodesInCycle);

    while(currentNodeId != visitNodeId) {
      currentNodeId = parents[currentNodeId];
      stk::util::insert_keep_sorted_and_unique(currentNodeId, nodesInCycle);
    }
    return;
  }

  parents[visitNodeId] = parentNodeId;
  visited[visitNodeId] = 1;

  for(unsigned node : nodeGraph[visitNodeId]) {
    if(node == parents[visitNodeId]) {
      continue;
    } else {
      get_cycles_helper(node, visitNodeId);
    }
  }

  visited[visitNodeId] = 2;
}

const std::vector<unsigned>& GraphCycleDetector::get_nodes_in_cycles()
{
  if(!initialized) {
    for(unsigned i = 0; i < numNodes; i++) {
      if(visited[i] != 0) { continue; }

      get_cycles_helper(i, i);
    }
    initialized = true;
  }

  return nodesInCycle;
}

}}}