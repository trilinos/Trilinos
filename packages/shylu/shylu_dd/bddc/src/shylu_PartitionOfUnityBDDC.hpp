
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

#ifndef BDDC_PARTITIONOFUNITY_H
#define BDDC_PARTITIONOFUNITY_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <mpi.h>
#include "shylu_enumsBDDC.hpp"
#include "shylu_UtilBDDC.hpp"
#include "shylu_PComm.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"  

#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Distributor.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {

template <class SX,
	  class SM,
	  class LO,
	  class GO> class PartitionOfUnity
{
 public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::Map<LO,GO>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
  typedef Tpetra::Export<LO,GO>                              Export;
  typedef Tpetra::Import<LO,GO>                              Import;

  PartitionOfUnity()
  {
  }

  PartitionOfUnity(LO numNodes,
		   const GO* nodeGlobalIDs,
		   const std::vector< std::vector<LO> > & subNodeBegin,
		   const std::vector< std::vector<LO> > & subNodes,
		   LO spatialDim,
		   RCP<Teuchos::ParameterList> Parameters,
		   RCP<const Teuchos::Comm<int> > Comm,
		   RCP<Tpetra::Distributor> & distributor,
		   const std::vector<int> & nodeSend,
		   const SM* xCoord=nullptr,
		   const SM* yCoord=nullptr,
		   const SM* zCoord=nullptr) :
  m_numNodes(numNodes), 
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_subNodeBegin(subNodeBegin),
    m_subNodes(subNodes),
    m_spatialDim(spatialDim),
    m_constructAdjacencyGraph
    (Parameters->get("Construct Subdomain Adjacency Graph", false)),
    m_useCorners(Parameters->get("Use Corners", false)),
    m_useEdges(Parameters->get("Use Edges", false)),
    m_useFaces(Parameters->get("Use Faces", false)),
    m_useVertexCoarseSpace(Parameters->get("Use Vertex Coarse Space", false)),
    m_Comm(Comm),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()),
    m_numNodesCoarse(0),
    m_startingSub(0),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord)
  {
    determineActiveNodes();
    determineStartingSubdomain();
    determineEquivalenceClasses(distributor, nodeSend);
    //  PrintSubdomainPU();
  }

  ~PartitionOfUnity()
  {
  }

  const std::vector< std::vector<LO> > & getSubdomainEquivClasses() const
  {
    return m_subdomainEquivClasses;
  }

  int getEquivFlag(LO equiv)
  {
    return m_equivFlag[equiv];
  }

  int getNumberOfMySubdomains() const
  {
    return m_subdomainEquivClasses.size();
  }

  const std::vector< std::vector<LO> > & getEquivClasses() const
  {
    return m_equivClasses;
  }

  LO getNumEquivClasses() const
  {
    return m_equivClasses.size();
  }

  const std::vector<LO> & getEquivCardinality() const
  {
    return m_equivCard;
  }

  const std::vector<GO> & getGlobalIDs() const 
  {
    return m_nodeGlobalIDsCoarse;
  }

  const std::vector< std::vector<LO> > & getCoarseNodeSubdomains() const 
  {
    return m_coarseNodeSubdomains;
  }

  const std::vector< std::vector<LO> > & getEquivCoarseNodes() const
  { 
    return m_equivCoarseNodes;
  }

  const std::vector<GO> & getEquivCoarseNodesGIDs() const
  { 
    return m_equivCoarseNodesGIDs;
  }

  const std::vector<LO> & getEquivCoarseNodesLIDs() const
  { 
    return m_equivCoarseNodesLIDs;
  }

  const std::vector<SM> & getEquivCoordX() const
  {
    return m_xCoordEquiv;
  }

  const std::vector<SM> & getEquivCoordY() const
  {
    return m_yCoordEquiv;
  }

  const std::vector<SM> & getEquivCoordZ() const
  {
    return m_zCoordEquiv;
  }

  RCP<const Teuchos::Comm<int> > getComm()
  {
    return m_Comm;
  }

  enum EquivType getEquivType(const LO equiv) const
  {
    BDDC_TEST_FOR_EXCEPTION(m_numAncestors.size() != m_equivCard.size(), 
			    std::runtime_error, "m_numAncestors size error");
    BDDC_TEST_FOR_EXCEPTION(m_equivCard[equiv] <= 1, std::runtime_error, 
			    "m_equivCard[equiv] must be greater than 1");
    enum EquivType equivType = FACE;
    if (m_equivCard[equiv] == 2) {
      if (m_spatialDim == 2) equivType = EDGE;
    }
    else {
      if (m_numAncestors[equiv] == 0) equivType = CORNER;
      else equivType = EDGE;
    }
    /*
    enum EquivType equivType = FACE;
    if (m_equivClasses[equiv].size() == 1) {
      equivType = CORNER;
    }
    else if (m_spatialDim == 2) {
      equivType = EDGE;
    }
    else if (m_spatialDim == 3) {
      if (m_equivCard[equiv] > 2) equivType = EDGE;
    }
    */
    return equivType;
  }

  const std::vector< std::vector<LO> > & getEquivClassSubdomains() const 
  {
    return m_equivClassSubdomains;
  }

  int getStartingSub() const
  {
    return m_startingSub;
  }

  int getNumActiveAncestors(LO equiv) const
  {
    LO numEquiv = m_equivClasses.size();
    const bool test = (equiv >= 0) && (equiv < numEquiv);
    BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			    "equiv out of range in getNumActiveAncestors");
    return m_numActiveAncestors[equiv];
  }

  void getConnectivityGraph(const GO* & subdomainGIDs,
			    const GO* & subdomainConnGIDs,
			    const GO* & subdomainConnProcs,
			    const LO* & subdomainConnBegin) {
    subdomainGIDs = m_subdomainGIDs.data();
    subdomainConnGIDs = m_subdomainConnGIDs.data();
    subdomainConnProcs = m_subdomainConnProcs.data();
    subdomainConnBegin = m_subdomainConnBegin.data();
  }

 private:
  void determineActiveNodes()
  {
    m_nodeIsActive.resize(m_numNodes, false);
    int numSub = m_subNodes.size();
    for (LO i=0; i<numSub; i++) {
      const std::vector<LO> & nodeBegin = m_subNodeBegin[i];
      const std::vector<LO> & subNodes = m_subNodes[i];
      LO numNodes = m_subNodes[i].size();
      for (LO j=0; j<numNodes; j++) {
	LO node = subNodes[j];
	if (nodeBegin[j+1] > nodeBegin[j]) m_nodeIsActive[node] = true;
      }
    }
  }

  void determineStartingSubdomain()
  {
    int numSub = m_subNodes.size();
    int numSubScan;
    Teuchos::scan<int, int> (*m_Comm, Teuchos::REDUCE_SUM, 1, &numSub, 
			     &numSubScan);
    m_startingSub = numSubScan - numSub;
  }

  void determineNodeSubdomains
    (RCP<Tpetra::Distributor> & distributor,
     const std::vector<int> & nodeSend,
     std::vector< std::vector<LO> > & nodeSubdomains,
     std::map<GO,LO> & subProcsMap) const
  {
    // nodeSubdomains[i] = list of subdomain globalIDs for all subdomains
    //  (across all processors) containing node with globalID m_nodeGlobalIDs[i]
    const LO myPID = m_Comm->getRank();
    // mySubData[i] = list of subdomain globalIDs of all on-processor
    //                subdomains for node with globalID m_nodeGlobalIDs[i]
    std::vector< std::vector<GO> > mySubData(m_numNodes);
    LO numSub = m_subNodes.size();
    for (LO i=0; i<numSub; i++) {
      const GO subGID = m_startingSub + i;
      for (size_t j=0; j<m_subNodes[i].size(); j++) {
	const LO node = m_subNodes[i][j];
	mySubData[node].push_back(subGID);
	if (m_constructAdjacencyGraph) {
	  mySubData[node].push_back(myPID);
	}
      }
    }
    std::vector< std::vector<GO> > unionSubData;
    bddc::unionData<LO,GO>(m_numNodes, m_nodeGlobalIDs, nodeSend.data(),
			   distributor, mySubData, unionSubData);
    nodeSubdomains.resize(m_numNodes);
    LO den(1);
    if (m_constructAdjacencyGraph) den = 2;
    for (LO i=0; i<m_numNodes; i++) {
      const LO length = unionSubData[i].size()/den;
      nodeSubdomains[i].resize(length);
      for (LO j=0; j<length; j++) {
	const LO index = den*j;
	const GO subGID = unionSubData[i][index];
	nodeSubdomains[i][j] = subGID;
	if (m_constructAdjacencyGraph) {
	  auto iter = subProcsMap.find(subGID);
	  if (iter == subProcsMap.end()) {
	  const GO proc = unionSubData[i][index+1];
	    subProcsMap.emplace(subGID, proc);
	  }
	}
      }
    }
  }

  void determineEquivalenceClasses
    (RCP<Tpetra::Distributor> & distributor,
     const std::vector<int> & nodeSend)
  {
    std::vector< std::vector<LO> > nodeSubdomains;
    std::map<GO,LO> subProcsMap;
    determineNodeSubdomains(distributor, nodeSend, nodeSubdomains, subProcsMap);
    equivalenceClasses(nodeSubdomains);
    determineEquivClassGlobalIDs();
    determineEquivClassSubdomains(nodeSubdomains);
    determineSubdomainEquivClasses(nodeSubdomains);
    std::vector< std::vector<LO> > subsForEquiv, equivsForSub;
    determineLocalSubdomainsForEquivClasses
      (nodeSubdomains, subsForEquiv, equivsForSub);
    determineCoarseNodeDependencies(subsForEquiv, equivsForSub);
    determineNumAncestors(subsForEquiv, equivsForSub, m_numAncestors);
    const bool onlyConsiderActiveAncestors = true;
    determineNumAncestors(subsForEquiv, equivsForSub, m_numActiveAncestors,
			  onlyConsiderActiveAncestors);
    constructSubdomainConnectivityGraph(nodeSubdomains, subProcsMap);
    determineCoarseNodeSubdomains(nodeSubdomains);
    markEquivClasses(distributor, nodeSend, nodeSubdomains);
  }

  void getAdjacentSubEquivs(const LO i,
			    const std::vector< std::vector<LO> > & nodeSubdomains,
			    std::vector< std::vector<LO> > & adjSubEquivs)
  {
    std::map<LO,LO> adjSubMap;
    const LO mySubGID = m_startingSub + i;
    LO numAdjSub(0);
    for (size_t j=0; j<m_subdomainEquivClasses[i].size(); j++) {
      const LO equiv = m_subdomainEquivClasses[i][j];
      const LO node = m_equivClasses[equiv][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      for (size_t k=0; k<subs.size(); k++) {
	const LO subGID = subs[k];
	if (subGID != mySubGID) {
	  auto iter = adjSubMap.find(subGID);
	  if (iter == adjSubMap.end()) {
	    adjSubMap.emplace(subGID, numAdjSub++);
	  }
	}
      }
    }
    adjSubEquivs.resize(numAdjSub);
    for (size_t j=0; j<m_subdomainEquivClasses[i].size(); j++) {
      const LO equiv = m_subdomainEquivClasses[i][j];
      const LO node = m_equivClasses[equiv][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      for (size_t k=0; k<subs.size(); k++) {
	const LO subGID = subs[k];
	if (subGID != mySubGID) {
	  auto iter = adjSubMap.find(subGID);
	  BDDC_TEST_FOR_EXCEPTION(iter == adjSubMap.end(), std::runtime_error, 
				  "subGID not found");
	  bool found(false);
	  const LO localSub = iter->second;
	  for (size_t m=0; m<adjSubEquivs[localSub].size(); m++) {
	    if (adjSubEquivs[localSub][m] == equiv) {
	      found = true;
	      break;
	    }
	  }
	  if (found == false) adjSubEquivs[localSub].push_back(equiv);
	}
      }
    }
  }

  void markEquivClasses(RCP<Tpetra::Distributor> & distributor, 
			const std::vector<int> & nodeSend, 
			const std::vector< std::vector<LO> > & nodeSubdomains)
  {
    LO numEquiv = m_equivClasses.size();
    m_equivFlag.resize(numEquiv, 0);
    LO numSub = m_subNodes.size();
    std::vector<LO> cornerEquivs;
    for (LO i=0; i<numSub; i++) {
      std::vector< std::vector<LO> > adjSubEquivs;
      getAdjacentSubEquivs(i, nodeSubdomains, adjSubEquivs);
      for (size_t j=0; j<adjSubEquivs.size(); j++) {
	LO numFace(0), numEdge(0);
	cornerEquivs.resize(0);
	for (size_t k=0; k<adjSubEquivs[j].size(); k++) {
	  const LO equiv = adjSubEquivs[j][k];
	  const enum EquivType equivT = getEquivType(equiv);
	  if (equivT == FACE) numFace++;
	  if (equivT == EDGE) numEdge++;
	  if (equivT == CORNER) cornerEquivs.push_back(equiv);
	}
	if ((numFace == 0) && (numEdge == 0)) {
	  if (cornerEquivs.size() > 1) {
	    for (size_t k=0; k<cornerEquivs.size(); k++) {
	      m_equivFlag[cornerEquivs[k]] = 1;
	    }
	  }
	}
      }
    }
    // communicate data
    std::vector<LO> nodeBeginEquiv(numEquiv+1, 0);
    std::vector<int> nodeSendEquiv(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      nodeBeginEquiv[i+1] = nodeBeginEquiv[i] + 1;
      const LO localNode = m_equivClasses[i][0];
      nodeSendEquiv[i] = nodeSend[localNode];
    }
    PComm<LO,GO,int> pComm(distributor, numEquiv, nodeBeginEquiv.data(),

			   m_nodeGlobalIDsCoarse.data(), nodeSendEquiv.data());
    int* ownedVals = pComm.getOwnedPtr();
    pComm.doExport(m_equivFlag.data(), ownedVals);
    pComm.doImport(ownedVals, m_equivFlag.data());

  }

  void determineCoarseNodeSubdomains
    (const std::vector< std::vector<LO> > & nodeSubdomains)
  {
    const LO numEquiv = m_equivClasses.size();
    m_coarseNodeSubdomains.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      const LO node = m_equivClasses[i][0];
      m_coarseNodeSubdomains[i] = nodeSubdomains[node];
      std::sort(m_coarseNodeSubdomains[i].begin(),
		m_coarseNodeSubdomains[i].end());
    }
  }

  void constructSubdomainConnectivityGraph
    (const std::vector< std::vector<LO> > & nodeSubdomains,
     const std::map<GO,LO> & subProcsMap)
  {
    if (m_constructAdjacencyGraph == false) return;
    LO numSub = m_subNodes.size();
    std::vector< std::vector< std::pair<LO, enum EquivType> > > 
      allConnections(numSub);
    determineAllSubdomainConnections(nodeSubdomains, allConnections);
    std::vector< std::vector<GO> > subdomainAdjacencies(numSub);
    determineSubdomainAdjacencies(allConnections, subdomainAdjacencies);
    m_subdomainGIDs.resize(numSub);
    m_subdomainConnBegin.resize(numSub+1, 0);
    LO numTerms(0);
    for (LO i=0; i<numSub; i++) {
      m_subdomainGIDs[i] = m_startingSub + i;
      numTerms += subdomainAdjacencies[i].size();
      m_subdomainConnBegin[i+1] = numTerms;
    }
    m_subdomainConnGIDs.resize(numTerms);
    m_subdomainConnProcs.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numSub; i++) {
      for (size_t j=0; j<subdomainAdjacencies[i].size(); j++) {
	const GO subGID = subdomainAdjacencies[i][j];
	m_subdomainConnGIDs[numTerms] = subGID;
	auto iter = subProcsMap.find(subGID);
	BDDC_TEST_FOR_EXCEPTION(iter == subProcsMap.end(), std::runtime_error, 
				"subGID not found");
	m_subdomainConnProcs[numTerms++] = iter->second;
      }
    }
  }
  /*
  static bool myLessThan(const std::pair<GO,LO> & left,
			 const GO right)
  {
    return left.first < right;
  }

  static bool myLessThan2(const std::pair<GO,LO> & left,
			  const std::pair<GO,LO> & right)
  {
    return left.first < right.first;
  }
  */
  void equivalenceClasses(const std::vector< std::vector<LO> > & nodeSubdomains)
  {
    std::vector< std::pair<double, LO> > subEncoding(m_numNodes);
    const double sqrt7(sqrt(7.0));
    for (LO i=0; i<m_numNodes; i++) {
      LO numSub = nodeSubdomains[i].size();
      double sum = 0;
      if ((numSub > 1) && (m_nodeIsActive[i] == true)) {
	for (LO j=0; j<numSub; j++) {
	  sum += sqrt(double(nodeSubdomains[i][j]+sqrt7));
	}
      }
      subEncoding[i] = std::make_pair(sum, i);
    }
    std::sort(subEncoding.begin(), subEncoding.end());
    double previousVal(0), currentVal, tol(1e-12);
    std::vector<LO> start, localNodes(m_numNodes);
    for (LO i=0; i<m_numNodes; i++) {
      currentVal = subEncoding[i].first;
      localNodes[i] = subEncoding[i].second;
      if (fabs(currentVal-previousVal) > tol) {
	start.push_back(i);
	previousVal = currentVal;
      }
    }
    start.push_back(m_numNodes);
    LO numEquivClasses = start.size()-1;
    m_equivClasses.resize(numEquivClasses);
    for (LO i=0; i<numEquivClasses; i++) {
      m_equivClasses[i].resize(start[i+1]-start[i]);
      for (LO j=start[i]; j<start[i+1]; j++) {
	m_equivClasses[i][j-start[i]] = localNodes[j];
      }
      sortByGlobalIDs(m_equivClasses[i]);
    }
  }

  void determineAllSubdomainConnections
    (const std::vector< std::vector<LO> > & nodeSubdomains,
     std::vector< std::vector< std::pair<LO, enum EquivType> > > & allConnections) const
  {
    LO numSub = m_subNodes.size();
    for (LO i=0; i<numSub; i++) {
      const LO mySubGID = m_startingSub + i;
      for (size_t j=0; j<m_subdomainEquivClasses[i].size(); j++) {
	const LO equiv = m_subdomainEquivClasses[i][j];
	const enum EquivType equivT = getEquivType(equiv);
	const LO node = m_equivClasses[equiv][0];
	const std::vector<LO> & subs = nodeSubdomains[node];
	for (size_t k=0; k<subs.size(); k++) {
	  if (subs[k] != mySubGID) {
	    allConnections[i].push_back(std::make_pair(subs[k], equivT));
	  }
	}
      }
    }
    /*
    LO numEquiv = m_equivClasses.size();
    for (LO i=0; i<numEquiv; i++) {
      enum EquivType equivT = getEquivType(i);
      LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      for (size_t j=0; j<subs.size(); j++) {
	LO localSub = subs[j] - m_startingSub;
	if ((localSub >= 0) && (localSub < numSub)) {
	  for (size_t k=0; k<subs.size(); k++) {
	    if (subs[k] != subs[j]) {
	      allConnections[localSub].push_back
		(std::make_pair(subs[k], equivT));
	    }
	  }
	}
      }
    }
    */
  }

  bool checkAdjacency(std::vector<int> & countEquivType) const
  {
    bool isAdjacent = false;
    if (m_spatialDim == 2) {
      if (countEquivType[EDGE] > 0) isAdjacent = true;
      if (countEquivType[CORNER] > 1) isAdjacent = true;
    }
    else if (m_spatialDim == 3) {
      if (countEquivType[FACE] > 0) isAdjacent = true;
      if (countEquivType[EDGE] > 1) isAdjacent = true;
      if (countEquivType[CORNER] > 2) isAdjacent = true;
    }
    return isAdjacent;
  }

  void determineSubdomainAdjacencies
    (std::vector< std::vector< std::pair<LO, enum EquivType> > > & allConnections,
     std::vector< std::vector<GO> > & subdomainAdjacencies) const
  {
    LO numSub = allConnections.size();
    for (LO i=0; i<numSub; i++) {
      std::sort(allConnections[i].begin(), allConnections[i].end());
      int previousSub = -1;
      std::vector<int> subBegin;
      for (size_t j=0; j<allConnections[i].size(); j++) {
	if (allConnections[i][j].first != previousSub) {
	  subBegin.push_back(j);
	  previousSub = allConnections[i][j].first;
	}
      }
      subBegin.push_back(allConnections[i].size());
      int numSub2 = subBegin.size() - 1;
      for (int j=0; j<numSub2; j++) {
	std::vector<int> countEquivType(3, 0);
	for (int k=subBegin[j]; k<subBegin[j+1]; k++) {
	  countEquivType[allConnections[i][k].second]++;
	}
	bool isAdjacent = checkAdjacency(countEquivType);
	if (isAdjacent == true) {
	  int index = subBegin[j];
	  GO adjSubdomain = allConnections[i][index].first;
	  subdomainAdjacencies[i].push_back(adjSubdomain);
	}
      }
    }
  }

  void sortByGlobalIDs(std::vector<LO> & localIDs) const
  {
    LO numNodes = localIDs.size();
    std::vector< std::pair<GO, LO> > globalAndLocal(numNodes);
    for (LO i=0; i<numNodes; i++) {
      LO localID = localIDs[i];
      GO globalID = m_nodeGlobalIDs[localID];
      globalAndLocal[i] = std::make_pair(globalID, localID);
    }
    std::sort(globalAndLocal.begin(), globalAndLocal.end());
    for (LO i=0; i<numNodes; i++) {
      localIDs[i] = globalAndLocal[i].second;
    }
  }

  void determineEquivClassGlobalIDs()
  {
    // Note: equivalence class nodes have already been sorted in ascending order
    //       by their globalIDs
    LO numEquiv = m_equivClasses.size();
    m_nodeGlobalIDsCoarse.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      LO localNode = m_equivClasses[i][0];
      m_nodeGlobalIDsCoarse[i] = m_nodeGlobalIDs[localNode];
    }
    m_numNodesCoarse = numEquiv;
  }

  void determineEquivClassSubdomains
    (std::vector< std::vector<LO> > & nodeSubdomains)
  {
    LO numEquiv = m_equivClasses.size();
    m_equivClassSubdomains.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      LO firstNode = m_equivClasses[i][0];
      m_equivClassSubdomains[i] = nodeSubdomains[firstNode];
    }
  }

  void determineSubdomainEquivClasses
    (std::vector< std::vector<LO> > & nodeSubdomains)
  {
    LO numSub = m_subNodes.size();
    m_subdomainEquivClasses.resize(numSub);
    LO numEquiv = m_equivClasses.size();
    m_equivCard.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      m_equivCard[i] = subs.size();
      for (size_t j=0; j<subs.size(); j++) {
	LO localSub = subs[j] - m_startingSub;
	if ((localSub >= 0) && (localSub < numSub)) {
	  m_subdomainEquivClasses[localSub].push_back(i);
	}
      }
    }
  }

  void determineLocalSubdomainsForEquivClasses
    (const std::vector< std::vector<LO> > & nodeSubdomains,
     std::vector< std::vector<LO> > & subsForEquiv,
     std::vector< std::vector<LO> > & equivsForSub)
  {
    // first, get list of subdomains for all equivalence classes
    const LO numEquiv = m_equivClasses.size();
    std::vector<LO> subList;
    for (LO i=0; i<numEquiv; i++) {
      const LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      for (size_t j=0; j<subs.size(); j++) subList.push_back(subs[j]);
    }
    std::vector<LO> uniqueSubs(subList);
    std::sort(uniqueSubs.begin(), uniqueSubs.end());
    auto iter = std::unique(uniqueSubs.begin(), uniqueSubs.end());
    uniqueSubs.erase(iter, uniqueSubs.end());
    // next, determine local subdomain numbers for each equivalence class
    // and vice versa
    const LO numSub = uniqueSubs.size();
    subsForEquiv.resize(numEquiv);
    equivsForSub.resize(numSub);
    for (LO i=0; i<numEquiv; i++) {
      const LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      subsForEquiv[i].resize(subs.size());
      for (size_t j=0; j<subs.size(); j++) {
	const LO sub = subs[j];
	const LO localSub = findEntry(uniqueSubs, sub);
	subsForEquiv[i][j] = localSub;
	equivsForSub[localSub].push_back(i);
      }
    }
  }

  void determineCoarseNodeDependencies
    (const std::vector< std::vector<LO> > & subsForEquiv,
     const std::vector< std::vector<LO> > & equivsForSub)
  {
    if (m_useVertexCoarseSpace == false) return;
    const LO numSub = equivsForSub.size();
    std::vector<bool> subFlag(numSub, false);
    const LO numEquiv = subsForEquiv.size();
    std::vector<LO> equivMap(numEquiv, -1), numAncestors(numEquiv, 0);
    // first determine coarse nodes
    LO numCoarseNodes(0), pass(1);
    m_equivCoarseNodes.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      const LO sub = subsForEquiv[i][0];
      for (size_t j=0; j<equivsForSub[sub].size(); j++) {
	const LO equiv = equivsForSub[sub][j];
	if (equiv != i) {
	  decendentCheck(subFlag, i, equiv, subsForEquiv, pass, numAncestors,
			 equivMap, m_equivCoarseNodes[i]);
	}
      }
      if (numAncestors[i] == 0) {
	equivMap[i] = numCoarseNodes;
	m_equivCoarseNodes[i].push_back(numCoarseNodes++);
      }
    }
    // determine coarse node dependencies
    pass = 2;
    for (LO i=0; i<numEquiv; i++) {
      const LO sub = subsForEquiv[i][0];
      for (size_t j=0; j<equivsForSub[sub].size(); j++) {
	const LO equiv = equivsForSub[sub][j];
	if (equiv != i) {
	  decendentCheck(subFlag, i, equiv, subsForEquiv, pass, numAncestors,
			 equivMap, m_equivCoarseNodes[i]);
	}
      }
    }
    m_equivCoarseNodesGIDs.resize(numCoarseNodes);
    m_equivCoarseNodesLIDs.resize(numCoarseNodes);
    m_xCoordEquiv.resize(numEquiv);
    m_yCoordEquiv.resize(numEquiv);
    m_zCoordEquiv.resize(numEquiv);
    numCoarseNodes = 0;
    for (LO i=0; i<numEquiv; i++) {
      getEquivCoords(i, m_xCoordEquiv[i], m_yCoordEquiv[i], m_zCoordEquiv[i]);
      if (equivMap[i] != -1) {
	const LO coarseNode = equivMap[i];
	m_equivCoarseNodesGIDs[coarseNode] = m_nodeGlobalIDsCoarse[i];
	m_equivCoarseNodesLIDs[coarseNode] = i;
      }
    }
  }

  void getEquivCoords(const LO equiv, 
		      SM & xCoordEquiv,
		      SM & yCoordEquiv,
		      SM & zCoordEquiv)
  {
    SM sumX(0), sumY(0), sumZ(0);
    const LO numNodeEquiv = m_equivClasses[equiv].size();
    for (LO i=0; i<numNodeEquiv; i++) {
      const LO node = m_equivClasses[equiv][i];
      sumX += m_xCoord[node];
      sumY += m_yCoord[node];
      sumZ += m_zCoord[node];
    }
    xCoordEquiv = sumX/numNodeEquiv;
    yCoordEquiv = sumY/numNodeEquiv;
    zCoordEquiv = sumZ/numNodeEquiv;
  }

  void decendentCheck(std::vector<bool> & subFlag, 
		      const LO i, 
		      const LO equiv, 
		      const std::vector< std::vector<LO> > & subsForEquiv,
		      const LO pass,
		      std::vector<LO> & numAncestors,
		      std::vector<LO> & equivMap,
		      std::vector<LO> & equivCoarseNodes)
  {
    const LO numSubEquiv = subsForEquiv[i].size();
    const LO numSubEquiv2 = subsForEquiv[equiv].size();
    for (LO k=0; k<numSubEquiv2; k++) {
      subFlag[subsForEquiv[equiv][k]] = true;
    }
    bool isaDecendent = true;
    for (LO k=0; k<numSubEquiv; k++) {
      if (subFlag[subsForEquiv[i][k]] == false) {
	isaDecendent = false;
	break;
      }
    }
    for (LO k=0; k<numSubEquiv2; k++) {
      subFlag[subsForEquiv[equiv][k]] = false;
    }
    if (isaDecendent) {
      if (pass == 1) numAncestors[i]++;
      else {
	if (numAncestors[equiv] == 0) {
	  const LO coarseNode = equivMap[equiv];
	  BDDC_TEST_FOR_EXCEPTION(coarseNode == -1, std::runtime_error, 
				  "invalid coarseNode");
	  equivCoarseNodes.push_back(coarseNode);
	}
      }
    }
  }

  void determineNumAncestors
    (const std::vector< std::vector<LO> > & subsForEquiv,
     const std::vector< std::vector<LO> > & equivsForSub,
     std::vector<LO> & numAncestors,
     const bool onlyConsiderActiveAncestors=false)
  {
    // determine number of active ancestors for each equivalence class
    LO numEquiv = m_equivClasses.size();
    LO numSub = equivsForSub.size();
    numAncestors.resize(numEquiv);
    std::vector<bool> subFlag(numSub, false);
    for (LO i=0; i<numEquiv; i++) {
      LO sub = subsForEquiv[i][0];
      LO numSubEquiv = subsForEquiv[i].size();
      for (size_t j=0; j<equivsForSub[sub].size(); j++) {
	LO equiv = equivsForSub[sub][j];
	if (onlyConsiderActiveAncestors) {
	  bool isActive = isActiveEquiv(equiv);
	  if (isActive == false) continue;
	}
	LO numSubEquiv2 = subsForEquiv[equiv].size();
	if (numSubEquiv2 <= numSubEquiv) continue;
	for (LO k=0; k<numSubEquiv2; k++) {
	  subFlag[subsForEquiv[equiv][k]] = true;
	}
	bool isAncestor = true;
	for (LO k=0; k<numSubEquiv; k++) {
	  if (subFlag[subsForEquiv[i][k]] == false) {
	    isAncestor = false;
	    break;
	  }
	}
	for (LO k=0; k<numSubEquiv2; k++) {
	  subFlag[subsForEquiv[equiv][k]] = false;
	}
	if (isAncestor == true) numAncestors[i]++;
      }
    }
  }

  bool isActiveEquiv(LO equiv) const
  {
    enum EquivType equivType = getEquivType(equiv);
    if ((equivType == CORNER) && (m_useCorners)) return true;
    if ((equivType == EDGE) && (m_useEdges)) return true;
    if ((equivType == FACE) && (m_useFaces)) return true;
    return false;
  }

  LO findEntry(std::vector<LO> & vector, 
	       LO value) const
  {
    auto iter = std::lower_bound(vector.begin(), vector.end(), value);
    BDDC_TEST_FOR_EXCEPTION(iter == vector.end(), std::runtime_error, 
			    "error in findEntry");
    LO index = iter - vector.begin();
    return index;
  }

 private: // variables
  LO m_numNodes;
  const GO* m_nodeGlobalIDs;
  const std::vector< std::vector<LO> > & m_subNodeBegin;
  const std::vector< std::vector<LO> > & m_subNodes;
  LO m_spatialDim;
  bool m_constructAdjacencyGraph, m_useCorners, m_useEdges, 
    m_useFaces, m_useVertexCoarseSpace;
  RCP<const Teuchos::Comm<int> > m_Comm;
  Tpetra::global_size_t m_IGO;
  LO m_numNodesCoarse;
  int m_startingSub;
  std::vector<LO> m_equivCard;
  std::vector<bool> m_nodeIsActive;
  std::vector< std::vector<LO> > m_equivClasses, m_subdomainEquivClasses,
    m_equivCoarseNodes;
  std::vector<GO> m_nodeGlobalIDsCoarse, m_equivCoarseNodesGIDs;
  std::vector<SM> m_xCoordEquiv, m_yCoordEquiv, m_zCoordEquiv;
  std::vector< std::vector<LO> > m_equivClassSubdomains;
  std::vector<LO> m_numAncestors, m_numActiveAncestors, m_equivCoarseNodesLIDs;
  RCP<CrsGraph> m_subdomainConnectivity;
  const SM *m_xCoord, *m_yCoord, *m_zCoord;
  std::vector<GO> m_subdomainGIDs, m_subdomainConnGIDs, m_subdomainConnProcs;
  std::vector<LO> m_subdomainConnBegin, m_equivFlag;
  std::vector< std::vector<LO> > m_coarseNodeSubdomains;
};

} // namespace bddc

#endif // BDDC_PARTITIONOFUNITY_H
  
