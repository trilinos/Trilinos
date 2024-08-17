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

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/baseImpl/ElemDeathImpl.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool is_node_connected_to_active_element_locally(const stk::mesh::BulkData &mesh,
                                                 stk::mesh::Entity node,
                                                 const stk::mesh::Part &activePart)
{
  bool activeNode = false;
  const int numElements = mesh.num_elements(node);
  const stk::mesh::Entity * elements = mesh.begin_elements(node);
  for (int elementI=0 ; elementI<numElements ; ++elementI) {    
    stk::mesh::Entity connectedElement = elements[elementI];
    stk::mesh::Bucket &connectedElementBucket = mesh.bucket(connectedElement);
    if (connectedElementBucket.owned() && connectedElementBucket.member(activePart)) {
      activeNode = true;
      break;
    }
  }    
  return activeNode;
}

stk::mesh::EntityVector
get_nodes_to_deactivate(const stk::mesh::BulkData& bulk,
                        const stk::mesh::EntityVector & deactivatedElements,
                        const stk::mesh::Part & activePart)
{
    stk::mesh::EntityVector nodesToDeactivate;

    stk::mesh::EntityVector potentiallyDeactivatedNodes;
    for (stk::mesh::Entity element : deactivatedElements) {
        const int numNodes = bulk.num_nodes(element);
        const stk::mesh::Entity * nodes = bulk.begin_nodes(element);
        for (int nodeI=0 ; nodeI<numNodes ; ++nodeI) {
            potentiallyDeactivatedNodes.push_back(nodes[nodeI]);
        }
    }
    stk::util::sort_and_unique(potentiallyDeactivatedNodes);

    stk::mesh::EntityVector nodesToCommunicate;
    for (stk::mesh::Entity node : potentiallyDeactivatedNodes) {
        if (bulk.bucket(node).owned() || bulk.bucket(node).shared()) {
            bool activeNode = is_node_connected_to_active_element_locally(bulk, node, activePart);
            if (!activeNode) {
              if (bulk.bucket(node).shared()) {
                nodesToCommunicate.push_back(node);
              }
              else {
                nodesToDeactivate.push_back(node);
              }
            }
        }
    }

    std::vector<int> sharedProcs;
    stk::CommSparse inquiryComm(bulk.parallel());
    pack_and_communicate(inquiryComm,
        [&bulk,&inquiryComm,&nodesToCommunicate,&sharedProcs]()
        {
            for (stk::mesh::Entity node : nodesToCommunicate) {
                const stk::mesh::EntityKey nodeKey = bulk.entity_key(node);
                bulk.comm_shared_procs(nodeKey,sharedProcs);
                for (int otherProc : sharedProcs) {
                    inquiryComm.send_buffer(otherProc).pack<stk::mesh::EntityId>(nodeKey.id());
                }
            }
        }
    );
    stk::mesh::EntityVector incomingNodes;
    unpack_communications(inquiryComm,
        [&bulk,&inquiryComm,&incomingNodes](int procId)
        {
            stk::mesh::EntityId nodeId;
            inquiryComm.recv_buffer(procId).unpack<stk::mesh::EntityId>(nodeId);
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
            STK_ThrowAssertMsg(bulk.is_valid(node),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            incomingNodes.push_back(node);
        }
    );

    std::map<stk::mesh::Entity,bool> nodeToActiveStatusMap;
    stk::CommSparse answerComm(bulk.parallel());
    pack_and_communicate(answerComm,
        [&bulk,&answerComm,&incomingNodes,&nodeToActiveStatusMap,&activePart]()
        {
            for (stk::mesh::Entity incomingNode : incomingNodes) {
                std::vector<int> sharingProcs;
                bulk.comm_shared_procs(bulk.entity_key(incomingNode),sharingProcs);
                bool activeStatus = is_node_connected_to_active_element_locally(bulk, incomingNode, activePart);
                for (int otherProc : sharingProcs) {
                    answerComm.send_buffer(otherProc).pack<stk::mesh::EntityId>(bulk.identifier(incomingNode));
                    answerComm.send_buffer(otherProc).pack<bool>(activeStatus);
                }
                auto nodeLocationInMap = nodeToActiveStatusMap.find(incomingNode);
                if (nodeLocationInMap == nodeToActiveStatusMap.end()) {
                    nodeToActiveStatusMap.emplace(incomingNode, activeStatus);
                }
                else {
                    nodeLocationInMap->second = nodeLocationInMap->second || activeStatus;
                }
            }
        }
    );

    unpack_communications(answerComm,
        [&bulk,&answerComm,&nodeToActiveStatusMap](int procId)
        {
            stk::mesh::EntityId nodeId;
            answerComm.recv_buffer(procId).unpack<stk::mesh::EntityId>(nodeId);
            bool activeStatus = false;
            answerComm.recv_buffer(procId).unpack<bool>(activeStatus);
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK,nodeId);
            STK_ThrowAssertMsg(bulk.is_valid(node),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            auto nodeLocationInMap = nodeToActiveStatusMap.find(node);
            if (nodeLocationInMap == nodeToActiveStatusMap.end()) {
                nodeToActiveStatusMap.emplace(node, activeStatus);
            }
            else {
                nodeLocationInMap->second = nodeLocationInMap->second || activeStatus;
            }
        }
    );

    for (auto nodeActiveStatusPair : nodeToActiveStatusMap) {
        stk::mesh::Entity node = nodeActiveStatusPair.first;
        bool nodeIsActiveOnAnyOtherProcessors = nodeActiveStatusPair.second;
        if (!nodeIsActiveOnAnyOtherProcessors) {
            nodesToDeactivate.push_back(node);
        }
    }

    return nodesToDeactivate;
}

} // namespace impl
} // namespace mesh
} // namespace stk

