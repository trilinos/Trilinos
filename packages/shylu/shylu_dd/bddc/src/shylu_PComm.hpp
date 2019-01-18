
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_PCOMM_H
#define BDDC_PCOMM_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <mpi.h>

#include "Tpetra_Distributor.hpp"
#include "shylu_errorBDDC.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {

template <class LO, class GO, class SX> class PComm
{
public:
  
  PComm(const PComm<LO,GO,SX> & pComm,
	RCP<const Teuchos::Comm<int> > Comm)
  {
    m_distributor = rcp( new Tpetra::Distributor(Comm) );
    m_distributor->createFromSendsAndRecvs(pComm.m_distributor->getProcsTo(),
					   pComm.m_distributor->getProcsFrom());
    m_sourceLength = pComm.m_sourceLength;
    m_targetLength = pComm.m_targetLength;
    m_sendCount = pComm.m_sendCount;
    m_recvCount = pComm.m_recvCount;
    m_sendDofs = pComm.m_sendDofs;
    m_recvDofs = pComm.m_recvDofs;
    m_nodeBeginOwned = pComm.m_nodeBeginOwned;
    m_ownedDofs = pComm.m_ownedDofs;
    m_ownedDofsTarget = pComm.m_ownedDofsTarget;
    m_sendVals = pComm.m_sendVals;
    m_recvVals = pComm.m_recvVals;
    m_nodeGlobalIDsOwned = pComm.m_nodeGlobalIDsOwned;
    m_numNodesSource = pComm.m_numNodesSource;
  }

  PComm(RCP<Tpetra::Distributor> distributor,
	const LO numNodes,
	const LO* nodeBegin,
	const GO* nodeGlobalIDs,
	const int* nodeSend) :
    m_distributor(distributor),
      m_sourceLength(nodeBegin[numNodes]),
      m_numNodesSource(numNodes)
  {
    setNodeSend(numNodes, nodeSend);
    std::map<GO,LO> nodeMap;
    initializeOwnedData(numNodes, nodeBegin, nodeGlobalIDs, nodeSend, nodeMap);
    std::vector<size_t> sendCountNodes, recvCountNodes;
    std::vector<GO> sendNodes, recvNodes;
    prepareSendData(numNodes, nodeBegin, nodeGlobalIDs,
		    nodeSend, sendCountNodes, sendNodes);
    receiveData(sendCountNodes, sendNodes, recvCountNodes, recvNodes);
    setRecvCount(recvCountNodes, recvNodes);
    processReceivedData(recvNodes, nodeMap);
    checkSizes();
  }

  PComm(RCP<Tpetra::Distributor> distributor,
	const LO numNodesSource,
	const LO* nodeBeginSource,
	const GO* nodeGlobalIDsSource,
	const int* nodeSend,
	const LO numNodesTarget,
	const LO* nodeBeginTarget,
	const GO* nodeGlobalIDsTarget) :
    m_distributor(distributor),
      m_sourceLength(nodeBeginSource[numNodesSource]),
      m_targetLength(nodeBeginTarget[numNodesTarget]),
      m_numNodesSource(numNodesSource)
  {
    setNodeSend(numNodesSource, nodeSend);
    std::map<GO,LO> nodeMapSource;
    for (LO i=0; i<numNodesSource; i++) {
      nodeMapSource.emplace(nodeGlobalIDsSource[i], i);
    }
    initializeOwnedData2(numNodesTarget, nodeBeginTarget, nodeGlobalIDsTarget, 
			nodeBeginSource, nodeMapSource);
    std::vector<size_t> sendCountNodes, recvCountNodes;
    std::vector<GO> sendNodes, recvNodes;
    prepareSendData(numNodesSource, nodeBeginSource, nodeGlobalIDsSource,
		    nodeSend, sendCountNodes, sendNodes);
    receiveData(sendCountNodes, sendNodes, recvCountNodes, recvNodes);
    setRecvCount(recvCountNodes, recvNodes);
    processReceivedData(recvNodes, numNodesTarget, nodeBeginTarget,
			nodeGlobalIDsTarget);
    checkSizes();
  }

  ~PComm()
  {
  }

  LO getNumNodesSource()
  {
    return m_numNodesSource;
  }

  LO getSourceLength() const
  {
    return m_sourceLength;
  }

  LO getTargetLength() const
  {
    return m_targetLength;
  }

  LO getOwnedLength() const
  {
    return m_ownedDofs.size();
  }

  SX* getSourcePtr() {
    m_sourceVals.resize(m_sourceLength);
    return m_sourceVals.data();
  }

  SX* getSourcePtr2() {
    m_sourceVals2.resize(m_sourceLength);
    return m_sourceVals2.data();
  }

  SX* getSourcePtr3() {
    m_sourceVals3.resize(m_sourceLength);
    return m_sourceVals3.data();
  }

  SX* getOwnedPtr() {
    m_ownedVals.resize(getOwnedLength());
    return m_ownedVals.data();
  }

  SX* getOwnedPtr2() {
    m_ownedVals2.resize(getOwnedLength());
    return m_ownedVals2.data();
  }

  const std::vector<GO> & getOwnedNodesGlobalIDs() const
  {
    return m_nodeGlobalIDsOwned;
  }

  const std::vector<LO> & getNodeBeginOwned() const
  {
    return m_nodeBeginOwned;
  }

  const int* getNodeSend() const
  {
    return m_nodeSend.data();
  }

  RCP<Tpetra::Distributor> getDistributor()
  {
    return m_distributor;
  }

  void doExportLocal(const SX* sourceVals,
		     SX* targetVals)
  {
    // only load on-processor owned values
    for (size_t i=0; i<m_ownedDofs.size(); i++) {
      const LO source = m_ownedDofs[i];
      const LO target = m_ownedDofsTarget[i];
      targetVals[target] = sourceVals[source];
    }
  }

  void doExport(const SX* sourceVals,
		SX* targetVals)
  {
    // Here we sum contributions from sending processors in sourceVals
    // to "owned" data in the array targetVals.
    // first, load on-processor owned values
    for (size_t i=0; i<m_ownedDofs.size(); i++) {
      const LO source = m_ownedDofs[i];
      const LO target = m_ownedDofsTarget[i];
      targetVals[target] = sourceVals[source];
    }
    // load and communicate values to send
    for (size_t i=0; i<m_sendDofs.size(); i++) {
      const LO dof = m_sendDofs[i];
      m_sendVals[i] = sourceVals[dof];
    }

    m_distributor->doPostsAndWaits
      (Teuchos::ArrayView<const SX>(m_sendVals),
       Teuchos::ArrayView<const size_t>(m_sendCount),
       Teuchos::ArrayView<SX>(m_recvVals),
       Teuchos::ArrayView<const size_t>(m_recvCount));

    // sum contributions
    for (size_t i=0; i<m_recvDofs.size(); i++) {
      const LO dof = m_recvDofs[i];
      targetVals[dof] += m_recvVals[i];
    }
  }

  void doImport(const SX* sourceVals,
		SX* targetVals)
  {
    // Here we "reverse communicate" owned data back to sending processors
    // (no summing takes place).
    // first, load on-processor owned values
    for (size_t i=0; i<m_ownedDofs.size(); i++) {
      const LO target = m_ownedDofs[i];
      const LO source = m_ownedDofsTarget[i];
      targetVals[target] = sourceVals[source];
    }

    for (size_t i=0; i<m_recvDofs.size(); i++) {
      const LO dof = m_recvDofs[i];
      m_recvVals[i] = sourceVals[dof];
    }

    m_distributor->doReversePostsAndWaits
      (Teuchos::ArrayView<const SX>(m_recvVals),
       Teuchos::ArrayView<const size_t>(m_recvCount),
       Teuchos::ArrayView<SX>(m_sendVals),
       Teuchos::ArrayView<const size_t>(m_sendCount));

    for (size_t i=0; i<m_sendDofs.size(); i++) {
      const LO dof = m_sendDofs[i];
      targetVals[dof] = m_sendVals[i];
    }
  }

 private: // member variables
  RCP<Tpetra::Distributor> m_distributor;
  LO m_sourceLength{0}, m_targetLength{0};
  std::vector<size_t> m_sendCount, m_recvCount;
  std::vector<LO> m_sendDofs, m_recvDofs, m_nodeBeginOwned, m_ownedDofs,
    m_ownedDofsTarget;
  std::vector<SX> m_sendVals, m_recvVals, m_sourceVals, m_sourceVals2,
    m_sourceVals3, m_ownedVals, m_ownedVals2;
  std::vector<GO> m_nodeGlobalIDsOwned;
  std::vector<int> m_nodeSend;
  LO m_numNodesSource;

 private: // member functions

  void setNodeSend(const LO numNodes, 
		   const int* nodeSend)
  {
    m_nodeSend.resize(numNodes);
    for (LO i=0; i<numNodes; i++) {
      m_nodeSend[i] = nodeSend[i];
    }
  }

  void initializeOwnedData2(const LO numNodesTarget, 
			    const LO* nodeBeginTarget, 
			    const GO* nodeGlobalIDsTarget, 
			    const LO* nodeBeginSource,
			    const std::map<GO,LO> & nodeMapSource)
  {
    for (LO i=0; i<numNodesTarget; i++) {
      const GO globalID = nodeGlobalIDsTarget[i];
      auto iter = nodeMapSource.find(globalID);
      if (iter != nodeMapSource.end()) {
	const LO node = iter->second;
	const bool test = ((nodeBeginSource[node+1] - nodeBeginSource[node]) ==
			   (nodeBeginTarget[i+1] - nodeBeginTarget[i]));
	BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
				"error in initializeOwnedData2");
	for (LO j=nodeBeginSource[node]; j<nodeBeginSource[node+1]; j++) {
	  m_ownedDofs.push_back(j);
	}
	for (LO j=nodeBeginTarget[i]; j<nodeBeginTarget[i+1]; j++) {
	  m_ownedDofsTarget.push_back(j);
	}
      }
    }
  }

  void initializeOwnedData(const LO numNodes, 
			   const LO* nodeBegin, 
			   const GO* nodeGlobalIDs, 
			   const int* nodeSend,
			   std::map<GO,LO> & nodeMap)
  {
    m_nodeBeginOwned.resize(1, 0);
    LO numDof(0), numNodesOwned(0);
    for (LO i=0; i<numNodes; i++) {
      if (nodeSend[i] == -1) {
	numDof += nodeBegin[i+1] - nodeBegin[i];
	m_nodeBeginOwned.push_back(numDof);
	const GO globalID = nodeGlobalIDs[i];
	m_nodeGlobalIDsOwned.push_back(globalID);
	nodeMap.emplace(globalID, numNodesOwned++);
      }
    }
    m_ownedDofs.resize(numDof);
    m_ownedDofsTarget.resize(numDof); // needed for consistency with other constructor
    numDof = 0;
    for (LO i=0; i<numNodes; i++) {
      if (nodeSend[i] == -1) {
	for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  m_ownedDofsTarget[numDof] = numDof;
	  m_ownedDofs[numDof++] = j;
	}
      }
    }
  }

  void checkSizes()
  {
    const LO numRecvs = m_distributor->getNumReceives();
    std::vector<size_t> recvCount(numRecvs);
    m_distributor->doPostsAndWaits
      (Teuchos::ArrayView<const size_t>(m_sendCount), 
       1,
       Teuchos::ArrayView<size_t>(recvCount));
    for (LO i=0; i<numRecvs; i++) {
      BDDC_TEST_FOR_EXCEPTION(m_recvCount[i] != recvCount[i], 
			      std::runtime_error, "error in checkSizes");
    }
  }

  void setRecvCount(const std::vector<size_t> & recvCountNodes, 
		    const std::vector<GO> & recvNodes)
  {
    const LO numRecvs = recvCountNodes.size();
    m_recvCount.resize(numRecvs);
    LO index(0);
    for (LO i=0; i<numRecvs; i++) {
      LO numDofs(0);
      const LO numPairs = recvCountNodes[i]/2;
      for (LO j=0; j<numPairs; j++) {
	numDofs += recvNodes[index+1];
	index += 2;
      }
      m_recvCount[i] = numDofs;
    }
  }

  void receiveData(const std::vector<size_t> & sendCount, 
		   const std::vector<GO> & sendNodes,
		   std::vector<size_t> & recvCount,
		   std::vector<GO> & recvNodes)
  {
    const LO numRecvs = m_distributor->getNumReceives();
    recvCount.resize(numRecvs);
    m_distributor->doPostsAndWaits(Teuchos::ArrayView<const size_t>(sendCount), 
				   1,
				   Teuchos::ArrayView<size_t>(recvCount));
    LO numTermsRecv(0);
    for (LO i=0; i<numRecvs; i++) numTermsRecv += recvCount[i];
    recvNodes.resize(numTermsRecv);
    m_distributor->doPostsAndWaits(Teuchos::ArrayView<const GO>(sendNodes),
				   Teuchos::ArrayView<const size_t>(sendCount),
				   Teuchos::ArrayView<GO>(recvNodes),
				   Teuchos::ArrayView<const size_t>(recvCount));
  }

  void processReceivedData(const std::vector<GO> & recvNodes,
			   const LO numNodesTarget,
			   const LO* nodeBeginTarget,
			   const GO* nodeGlobalIDsTarget)
  {
    std::map<GO,LO> nodeMapTarget;
    for (LO i=0; i<numNodesTarget; i++) {
      nodeMapTarget.emplace(nodeGlobalIDsTarget[i], i);
    }
    const LO numPairs = recvNodes.size()/2;
    m_recvDofs.resize(0);
    for (LO i=0; i<numPairs; i++) {
      const GO nodeGlobalID = recvNodes[2*i];
      auto iter = nodeMapTarget.find(nodeGlobalID);
      BDDC_TEST_FOR_EXCEPTION(iter == nodeMapTarget.end(), 
			      std::runtime_error, "nodeGlobalID not found");
      const LO node = iter->second;
      const bool test = recvNodes[2*i+1] == 
	(nodeBeginTarget[node+1] - nodeBeginTarget[node]);
      BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			      "size error in processRecievedData");
      for (LO j=nodeBeginTarget[node]; j<nodeBeginTarget[node+1]; j++) {
	m_recvDofs.push_back(j);
      }
    }
    m_recvDofs.shrink_to_fit();
    m_recvVals.resize(m_recvDofs.size());
    
  }

  void processReceivedData(const std::vector<GO> & recvNodes,
			   std::map<GO,LO> & nodeMap)
  {
    const LO numPairs = recvNodes.size()/2;
    std::vector<LO> recvNodeLIDs(numPairs, -1);
    LO numNodesOwned = m_nodeGlobalIDsOwned.size();
    LO numDofsOwned = m_ownedDofs.size();
    for (LO i=0; i<numPairs; i++) {
      const GO nodeGlobalID = recvNodes[2*i];
      auto iter = nodeMap.find(nodeGlobalID);
      if (iter == nodeMap.end()) {
	m_nodeGlobalIDsOwned.push_back(nodeGlobalID);
	recvNodeLIDs[i] = numNodesOwned;
	nodeMap.emplace(nodeGlobalID, numNodesOwned++);
	numDofsOwned += recvNodes[2*i+1];
        m_nodeBeginOwned.push_back(numDofsOwned);
      }
      else {
	recvNodeLIDs[i] = iter->second;
      }
    }
    for (LO i=0; i<numPairs; i++) {
      BDDC_TEST_FOR_EXCEPTION(recvNodeLIDs[i] == -1, std::runtime_error, 
			      "invalid recvNodeLIDs[i]");
    }
    m_recvDofs.resize(0);
    for (LO i=0; i<numPairs; i++) {
      const LO node = recvNodeLIDs[i];
      for (LO j=m_nodeBeginOwned[node]; j<m_nodeBeginOwned[node+1]; j++) {
	m_recvDofs.push_back(j);
      }
    }
    m_recvDofs.shrink_to_fit();
    m_recvVals.resize(m_recvDofs.size());
  }

  void prepareSendData(const LO numNodes, 
		       const LO* nodeBegin,
		       const GO* nodeGlobalIDs,
		       const int* nodeSend,
		       std::vector<size_t> & sendCount,
		       std::vector<GO> & sendNodes)
  {
    const LO numSends = m_distributor->getNumSends();
    sendCount.assign(numSends, 0);
    m_sendCount.assign(numSends, 0);
    for (LO i=0; i<numNodes; i++) {
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      if (numDofNode == 0) continue;
      const int proc = nodeSend[i];
      if (proc != -1) {
	sendCount[proc] += 2;
	m_sendCount[proc] += numDofNode;
      }
    }
    std::vector<LO> sendStart(numSends+1, 0), sendStartDofs(numSends+1, 0);
    for (LO i=0; i<numSends; i++) {
      sendStart[i+1] = sendStart[i] + sendCount[i];
      sendCount[i] = 0;
      sendStartDofs[i+1] = sendStartDofs[i] + m_sendCount[i];
      m_sendCount[i] = 0;
    }
    sendNodes.resize(sendStart[numSends]);
    m_sendDofs.resize(sendStartDofs[numSends]);
    m_sendVals.resize(sendStartDofs[numSends]);
    for (LO i=0; i<numNodes; i++) {
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      if (numDofNode == 0) continue;
      const int proc = nodeSend[i];
      if (proc != -1) {
	LO index = sendStart[proc] + sendCount[proc];
	sendNodes[index+0] = nodeGlobalIDs[i];
	sendNodes[index+1] = numDofNode;
	sendCount[proc] += 2;
	for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  index = sendStartDofs[proc] + m_sendCount[proc];
	  m_sendDofs[index] = j;
	  m_sendCount[proc]++;
	}
      }
    }
  }

};

} // namespace bddc

#endif // BDDC_PCOMM_H
