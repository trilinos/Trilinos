
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef PARTITIONOFUNITYBDDC_H
#define PARTITIONOFUNITYBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <mpi.h>
#include "shylu_enumsBDDC.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"  

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
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
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::Map<LO,GO,Node>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO,Node>                            CrsGraph;
  typedef Tpetra::Export<LO,GO,Node>                              Export;
  typedef Tpetra::Import<LO,GO,Node>                              Import;

  PartitionOfUnity()
  {
  }

  PartitionOfUnity(LO numNodes,
		   const GO* nodeGlobalIDs,
		   const std::vector< std::vector<LO> > & subNodeBegin,
		   const std::vector< std::vector<LO> > & subNodes,
		   LO** subRowBegin,
		   LO** subColumns,
		   SX** subValues,
		   LO spatialDim,
		   RCP<Teuchos::ParameterList> & Parameters,
		   RCP<const Teuchos::Comm<int> > Comm) :
  m_numNodes(numNodes), 
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_subNodeBegin(subNodeBegin),
    m_subNodes(subNodes),
    m_subRowBegin(subRowBegin),
    m_subColumns(subColumns),
    m_subValues(subValues),
    m_spatialDim(spatialDim),
    m_constructAdjacencyGraph
    (Parameters->get("Construct Subdomain Adjacency Graph", false)),
    m_addCorners(Parameters->get("Add Corners", false)),
    m_useCorners(Parameters->get("Use Corners", false)),
    m_useEdges(Parameters->get("Use Edges", false)),
    m_useFaces(Parameters->get("Use Faces", false)),
    m_Comm(Comm),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()),
    m_numNodesCoarse(0),
    m_startingSub(0)
  {
    determineActiveNodes();
    determineStartingSubdomain();
    determineEquivalenceClasses();
    //  PrintSubdomainPU();
  }

  ~PartitionOfUnity()
  {
  }

  const std::vector< std::vector<LO> > & getSubdomainEquivClasses() const
  {
    return m_subdomainEquivClasses;
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

  RCP<const Teuchos::Comm<int> > getComm()
  {
    return m_Comm;
  }

  enum EquivType getEquivType(LO equiv) const
  {
    LO numNodes = m_equivClasses[equiv].size();
    LO equivCard = m_equivCard[equiv];
    enum EquivType equivType = FACE;
    assert (equivCard > 1);
    if (numNodes == 1) equivType = CORNER;
    else if (m_spatialDim == 2) equivType = EDGE;
    else if (m_spatialDim == 3) {
      if (equivCard == 2) equivType = FACE;
      else equivType = EDGE;
    }
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
    assert ((equiv >= 0) && (equiv < numEquiv));
    return m_numActiveAncestors[equiv];
  }

  RCP<const CrsGraph> getConnectivityGraph() {
    return m_subdomainConnectivity;
  }

 private:
  void determineAdjacentNodes(LO & numAdjNodes, 
			      std::vector<LO> & adjNodes, 
			      std::vector<bool> & nodeFlag,
			      std::vector< std::vector<LO> > & subDofNodes,
			      LO sub, 
			      LO node) const
  {
    const LO* rowBegin = m_subRowBegin[sub];
    const LO* columns = m_subColumns[sub];
    const SX* values = m_subValues[sub];
    const std::vector<LO> & nodeBegin = m_subNodeBegin[sub];
    const std::vector<LO> & dofNodes = subDofNodes[sub];
    const std::vector<LO> & subNodes = m_subNodes[sub];
    for (LO i=nodeBegin[node]; i<nodeBegin[node+1]; i++) {
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	LO col = columns[j];
	LO node2 = subNodes[dofNodes[col]];
	if (std::abs(values[j]) > 0) {
	  if (nodeFlag[node2] == false) {
	    nodeFlag[node2] = true;
	    adjNodes[numAdjNodes] = node2;
	    numAdjNodes++;
	  }
	}
      }
    }
  }

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
    (std::vector< std::vector<LO> > & nodeSubdomains) const
  {
    LO numSub = m_subNodes.size();
    std::vector<GO> subGIDs(numSub);
    nodeSubdomains.resize(m_numNodes);
    for (int i=0; i<numSub; i++) {
      subGIDs[i] = m_startingSub + i;
      for (size_t j=0; j<m_subNodes[i].size(); j++) {
	nodeSubdomains[m_subNodes[i][j]].push_back(i);
      }
    }
    Teuchos::ArrayRCP<size_t> count(m_numNodes);
    for (LO i=0; i<m_numNodes; i++) count[i] = nodeSubdomains[i].size();
    RCP<const Map> nodeMap = 
      rcp(new Map(m_IGO, 
		  Teuchos::ArrayView<const GO>(m_nodeGlobalIDs,m_numNodes), 
		  0, m_Comm));
    RCP<const Map> colMap =
      rcp( new Map(m_IGO, Teuchos::ArrayView<GO>(subGIDs), 0, m_Comm) );
    CrsGraph nodeSubs(nodeMap, colMap, count, Tpetra::StaticProfile);
    for (LO i=0; i<m_numNodes; i++) {
      nodeSubs.insertLocalIndices(i, Teuchos::ArrayView<LO>(nodeSubdomains[i]));
    }
    RCP<const Map> nodeMap1to1 = Tpetra::createOneToOne<LO,GO,Node>(nodeMap);
    nodeSubs.fillComplete(colMap, nodeMap1to1);
    CrsGraph nodeSubs1to1(nodeMap1to1, 0);
    Export Exporter(nodeMap, nodeMap1to1);
    nodeSubs1to1.doExport(nodeSubs, Exporter, Tpetra::ADD);
    nodeSubs1to1.fillComplete(colMap, nodeMap1to1);
    Import Importer(nodeMap1to1, nodeMap);
    CrsGraph nodeSubsAll(nodeMap, 0);
    nodeSubsAll.doImport(nodeSubs1to1, Importer, Tpetra::INSERT);
    nodeSubsAll.fillComplete(colMap, nodeMap1to1);
    Teuchos::ArrayView<const LO> Indices;
    for (LO i=0; i<m_numNodes; i++) {
      nodeSubsAll.getLocalRowView(i, Indices);
      LO numSubNode = Indices.size();
      nodeSubdomains[i].resize(numSubNode);
      for (LO j=0; j<numSubNode; j++) {
	nodeSubdomains[i][j] = 
	  nodeSubsAll.getColMap()->getGlobalElement(Indices[j]);
      }
      std::sort(nodeSubdomains[i].begin(), nodeSubdomains[i].end());
    }
  }

  void determineNodeConnectivity
    (std::vector< std::vector<LO> > & nodeConnectivity) const
  {
    LO numSub = m_subNodes.size();
    std::vector< std::vector<LO> > subDofNodes(numSub);
    std::vector< std::vector< std::pair<LO, LO> > > subsForNodes(m_numNodes);
    for (LO i=0; i<numSub; i++) {
      LO numNode = m_subNodes[i].size();
      LO numDof = m_subNodeBegin[i][numNode];
      subDofNodes[i].resize(numDof);
      for (LO j=0; j<numNode; j++) {
	LO node = m_subNodes[i][j];
	subsForNodes[node].push_back(std::make_pair(i, j));
	for (LO k=m_subNodeBegin[i][j]; k<m_subNodeBegin[i][j+1]; k++) {
	  subDofNodes[i][k] = j;
	}
      }
    }
    std::vector<bool> nodeFlag(m_numNodes, false);
    std::vector<LO> adjNodes(m_numNodes);
    nodeConnectivity.resize(m_numNodes);
    for (LO i=0; i<m_numNodes; i++) {
      LO numAdjNodes(0);
      for (size_t j=0; j<subsForNodes[i].size(); j++) {
	LO sub = subsForNodes[i][j].first;
	LO node = subsForNodes[i][j].second;
	determineAdjacentNodes(numAdjNodes, adjNodes, nodeFlag, subDofNodes, 
			       sub, node);
      }
      nodeConnectivity[i].resize(numAdjNodes);
      for (LO j=0; j<numAdjNodes; j++) {
	nodeConnectivity[i][j] = adjNodes[j];
	nodeFlag[adjNodes[j]] = false;
      }
    }
  }

  void determineEquivalenceClasses()
  {
    std::vector< std::vector<LO> > nodeSubdomains;
    determineNodeSubdomains(nodeSubdomains);
    std::vector< std::vector<LO> > nodeConnectivity;
    determineNodeConnectivity(nodeConnectivity);
    std::vector<int> component(m_numNodes, 0);
    equivalenceClasses(nodeSubdomains, nodeConnectivity, component);
    determineEquivClassComponents(nodeConnectivity, component);
    equivalenceClasses(nodeSubdomains, nodeConnectivity, component);
    determineEquivClassGlobalIDs();
    determineEquivClassSubdomains(nodeSubdomains);
    determineSubdomainEquivClasses(nodeSubdomains);
    determineNumActiveAncestors(nodeSubdomains);
    constructSubdomainConnectivityGraph(nodeSubdomains);
  }

  void constructSubdomainConnectivityGraph
    (std::vector< std::vector<LO> > & nodeSubdomains)
  {
    if (m_constructAdjacencyGraph == false) return;
    LO numSub = m_subNodes.size();
    std::vector< std::vector< std::pair<LO, enum EquivType> > > 
      allConnections(numSub);
    determineAllSubdomainConnections(nodeSubdomains, allConnections);
    std::vector< std::vector<GO> > subdomainAdjacencies(numSub);
    determineSubdomainAdjacencies(allConnections, subdomainAdjacencies);
    std::vector<GO> subdomainGIDs(numSub);
    Teuchos::ArrayRCP<size_t> count(numSub);
    for (LO i=0; i<numSub; i++) {
      subdomainGIDs[i] = m_startingSub + i;
      count[i] = subdomainAdjacencies[i].size();
    }
    RCP<const Map> subMap =
      rcp( new Map(m_IGO, Teuchos::ArrayView<GO>(subdomainGIDs), 0, m_Comm) );
    m_subdomainConnectivity = 
      rcp( new CrsGraph(subMap, count, Tpetra::StaticProfile) );
    for (LO i=0; i<numSub; i++) {
      m_subdomainConnectivity->insertGlobalIndices
	(subdomainGIDs[i], Teuchos::ArrayView<GO>(subdomainAdjacencies[i]));
    }
    m_subdomainConnectivity->fillComplete();
  }

  void equivalenceClasses(std::vector< std::vector<LO> > & nodeSubdomains,
			  std::vector< std::vector<LO> > & nodeConnectivity,
			  std::vector<int> & component)
  {
    std::vector< std::pair<double, LO> > subEncoding(m_numNodes);
    for (LO i=0; i<m_numNodes; i++) {
      LO numSub = nodeSubdomains[i].size();
      double sum = 0;
      if ((numSub > 1) &&	(m_nodeIsActive[i] == true)) {
	for (LO j=0; j<numSub; j++) {
	  sum += sqrt(double(nodeSubdomains[i][j]+
			     sqrt(double(7))*component[i]+1.7));
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
    (std::vector< std::vector<LO> > & nodeSubdomains,
     std::vector< std::vector< std::pair<LO, enum EquivType> > > & allConnections) const
  {
    LO numSub = m_subNodes.size();
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

  void determineEquivClassComponents
    (std::vector< std::vector<LO> > & nodeConnectivity,
     std::vector<int> & component)
  {
    LO numEquivClasses = m_equivClasses.size();
    std::vector<LO> A1, A2, nodeMap(m_numNodes, -1);
    for (LO i=0; i<numEquivClasses; i++) {
      A1.resize(0);
      const std::vector<LO> & equivNodes = m_equivClasses[i];
      LO numEquivNodes = equivNodes.size();
      LO numTerms = 0;
      for (LO j=0; j<numEquivNodes; j++) nodeMap[equivNodes[j]] = j;
      A2.resize(numEquivNodes+1, 0);
      for (LO j=0; j<numEquivNodes; j++) {
	LO node1 = equivNodes[j];
	for (size_t k=0; k<nodeConnectivity[node1].size(); k++) {
	  LO node2 = nodeMap[nodeConnectivity[node1][k]];
	  if (node2 != -1) {
	    A1.push_back(node2);
	    numTerms++;
	  }
	}
	A2[j+1] = numTerms;
      }
      for (LO j=0; j<numEquivNodes; j++) nodeMap[equivNodes[j]] = -1;
      int *componentEQ(0);
      determineComponents(&A1[0], &A2[0], numEquivNodes, componentEQ);
      for (LO j=0; j<numEquivNodes; j++) {
	component[equivNodes[j]] = componentEQ[j];
      }
      delete [] componentEQ;
    }
  }

  void determineComponents(LO *A1, 
			   LO *A2, 
			   LO N, 
			   LO* & component) const
  {
    LO i;
    component = new LO[N];
    for (i=0; i<N; i++) component[i] = -1;
    componentsFunction(A1, A2, N, component);
  }

  void componentsFunction(LO *A1, 
			  LO *A2, 
			  LO N, 
			  LO* component) const
  {
    LO i, comp_num;
    comp_num = 0;
    for (i=0; i<N; i++) {
      if (component[i] == -1) {
	depthFirstSearch(i, comp_num, component, A1, A2);
	comp_num++;
      }
    }
  }

  void depthFirstSearch(const LO v, 
			const LO comp_num, 
			LO* component,
			LO* A1, 
			LO* A2) const
  {
    LO i, adj_vertex;
    component[v] = comp_num;
    for (i=A2[v]; i<A2[v+1]; i++) {
      adj_vertex = A1[i];
      if (component[adj_vertex] == -1) 
	depthFirstSearch(adj_vertex, comp_num, component, A1, A2);
    }
  }

  void determineNumActiveAncestors
    (std::vector< std::vector<LO> > & nodeSubdomains)
  {
    // first, get list of subdomains for all equivalence classes
    LO numEquiv = m_equivClasses.size();
    std::vector<LO> subList;
    for (LO i=0; i<numEquiv; i++) {
      LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      for (size_t j=0; j<subs.size(); j++) subList.push_back(subs[j]);
    }
    std::vector<LO> uniqueSubs(subList);
    std::sort(uniqueSubs.begin(), uniqueSubs.end());
    auto iter = std::unique(uniqueSubs.begin(), uniqueSubs.end());
    uniqueSubs.erase(iter, uniqueSubs.end());
    // next, determine local subdomain numbers for each equivalence class
    // and vice versa
    LO numSub = uniqueSubs.size();
    std::vector< std::vector<LO> > subsForEquiv(numEquiv);
    std::vector< std::vector<LO> > equivsForSub(numSub);
    for (LO i=0; i<numEquiv; i++) {
      LO node = m_equivClasses[i][0];
      const std::vector<LO> & subs = nodeSubdomains[node];
      subsForEquiv[i].resize(subs.size());
      for (size_t j=0; j<subs.size(); j++) {
	LO sub = subs[j];
	LO localSub = findEntry(uniqueSubs, sub);
	subsForEquiv[i][j] = localSub;
	equivsForSub[localSub].push_back(i);
      }
    }
    // determine number of active ancestors for each equivalence class
    m_numActiveAncestors.resize(numEquiv);
    std::vector<bool> subFlag(numSub, false);
    for (LO i=0; i<numEquiv; i++) {
      LO sub = subsForEquiv[i][0];
      LO numSubEquiv = subsForEquiv[i].size();
      for (size_t j=0; j<equivsForSub[sub].size(); j++) {
	LO equiv = equivsForSub[sub][j];
	bool isActive = isActiveEquiv(equiv);
	if (isActive == false) continue;
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
	if (isAncestor == true) m_numActiveAncestors[i]++;
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
    assert (iter != vector.end());
    LO index = iter - vector.begin();
    return index;
  }

 private: // variables
  LO m_numNodes;
  const GO* m_nodeGlobalIDs;
  const std::vector< std::vector<LO> > & m_subNodeBegin;
  const std::vector< std::vector<LO> > & m_subNodes;
  LO **m_subRowBegin,  **m_subColumns;
  SX **m_subValues;
  LO m_spatialDim;
  bool m_constructAdjacencyGraph, m_addCorners, m_useCorners, m_useEdges, 
    m_useFaces;
  RCP<const Teuchos::Comm<int> > m_Comm;
  Tpetra::global_size_t m_IGO;
  LO m_numNodesCoarse;
  int m_startingSub;
  std::vector<LO> m_equivCard;
  std::vector<bool> m_nodeIsActive;
  std::vector< std::vector<LO> > m_equivClasses, m_subdomainEquivClasses;
  std::vector<GO> m_nodeGlobalIDsCoarse;
  std::vector< std::vector<LO> > m_equivClassSubdomains;
  std::vector<LO> m_numActiveAncestors;
  RCP<CrsGraph> m_subdomainConnectivity;

};

} // namespace bddc

#endif // PARTITIONOFUNITYBDDC_H
  
