
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

#define TemplateArgs1NPU class SX,class SM,class LO,class GO
#define TemplateArgs2NPU SX,SM,LO,GO

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

  PartitionOfUnity();
  PartitionOfUnity(LO numNodes,
		   const GO* nodeGlobalIDs,
		   const std::vector< std::vector<LO> > & subNodeBegin,
		   const std::vector< std::vector<LO> > & subNodes,
		   const std::vector< std::vector<LO> > & subRowBegin,
		   const std::vector< std::vector<LO> > & subColumns,
		   const std::vector< std::vector<SX> > & subValues,
		   const LO* nodeOnBoundary,
		   LO spatialDim,
		   bool addCorners,
		   RCP<const Teuchos::Comm<int> > m_Comm);
  ~PartitionOfUnity();
  const std::vector< std::vector<LO> > & getSubdomainEquivClasses() const
  {
    return m_subdomainEquivClasses;
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
  enum EquivType getEquivType(LO equiv) 
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
  /*
  const std::vector< std::vector<SM> > & getEquivClassWeights() const {
    return m_equivClassWeights;
  }
  const std::vector<GO> & getAdjacentSubdomains() const {
    return m_adjacentSubdomains;
  }
  */
 private:
  void determineAdjacentNodes(LO & numAdjNodes, 
			      std::vector<LO> & adjNodes, 
			      std::vector<bool> & nodeFlag,
			      std::vector< std::vector<LO> > & subDofNodes,
			      LO sub, 
			      LO node) const;
  void determineActiveNodes();
  void determineStartingSubdomain();
  void determineNodeSubdomains
    (std::vector< std::vector<LO> > & nodeSubdomains) const;
  void determineNodeConnectivity
    (std::vector< std::vector<LO> > & nodeConnectivity) const;
  void determineEquivalenceClasses();
  void determineAdjacentSubdomains
    (std::vector< std::vector<LO> > & nodeSubdomains);
  void equivalenceClasses(std::vector< std::vector<LO> > & nodeSubdomains,
			  std::vector< std::vector<LO> > & nodeConnectivity,
			  std::vector<int> & component);
  void sortByGlobalIDs(std::vector<LO> & localIDs) const;
  void addEdgeEndpoints(std::vector< std::vector<LO> > & nodeSubdomains, 
			std::vector< std::vector<LO> > & nodeConnectivity);
  bool edgeCandidate(LO equiv,
		     std::vector< std::vector<LO> > & nodeSubdomains) const;
  void determineCornersToAdd(LO equiv, 
			     std::vector<bool> & nodeIsCorner, 
			     std::vector< std::vector<LO> > & nodeConnectivity, 
			     std::vector<LO> & nodeMap,
			     std::vector<LO> & cornerIndices, 
			     std::vector<LO> & edgeIndices);
  void determineEquivClassGlobalIDs();
  void determineEquivClassSubdomains
    (std::vector< std::vector<LO> > & nodeSubdomains);
  void determineSubdomainEquivClasses
    (std::vector< std::vector<LO> > & nodeSubdomains);
  void determineEquivClassComponents
    (std::vector< std::vector<LO> > & nodeConnectivity,
     std::vector<int> & component);
  void determineComponents(LO *A1, 
			   LO *A2, 
			   LO N, 
			   LO* & component);
  void componentsFunction(LO *A1, 
			  LO *A2, 
			  LO N, 
			  LO* component);
  void depthFirstSearch(const LO v, 
			const LO comp_num, 
			LO* component,
			LO* A1, 
			LO* A2);

 private: // variables
  LO m_numNodes;
  const GO* m_nodeGlobalIDs;
  const std::vector< std::vector<LO> > & m_subNodeBegin;
  const std::vector< std::vector<LO> > & m_subNodes;
  const std::vector< std::vector<LO> > & m_subRowBegin;
  const std::vector< std::vector<LO> > & m_subColumns;
  const std::vector< std::vector<SX> > & m_subValues;
  const LO* m_nodeOnBoundary;
  LO m_spatialDim;
  bool m_addCorners;
  RCP<const Teuchos::Comm<int> > m_Comm;
  Tpetra::global_size_t m_IGO;
  LO m_numNodesCoarse;
  int m_startingSub;
  std::vector<LO> m_equivCard;
  std::vector<bool> m_nodeIsActive;
  std::vector< std::vector<LO> > m_equivClasses, m_adjacentSubdomains,
    m_subdomainEquivClasses;
  std::vector<GO> m_nodeGlobalIDsCoarse;
  std::vector< std::vector<LO> > m_equivClassSubdomains;
};

template <TemplateArgs1NPU>
PartitionOfUnity<TemplateArgs2NPU>::
PartitionOfUnity()
{
}

template <TemplateArgs1NPU>
PartitionOfUnity<TemplateArgs2NPU>::
~PartitionOfUnity()
{
}

template <TemplateArgs1NPU>
PartitionOfUnity<TemplateArgs2NPU>::
  PartitionOfUnity(LO numNodes,
		   const GO* nodeGlobalIDs,
		   const std::vector< std::vector<LO> > & subNodeBegin,
		   const std::vector< std::vector<LO> > & subNodes,
		   const std::vector< std::vector<LO> > & subRowBegin,
		   const std::vector< std::vector<LO> > & subColumns,
		   const std::vector< std::vector<SX> > & subValues,
		   const LO* nodeOnBoundary,
		   LO spatialDim,
		   bool addCorners,
		   RCP<const Teuchos::Comm<int> > Comm) :
 m_numNodes(numNodes), 
   m_nodeGlobalIDs(nodeGlobalIDs),
   m_subNodeBegin(subNodeBegin),
   m_subNodes(subNodes),
   m_subRowBegin(subRowBegin),
   m_subColumns(subColumns),
   m_subValues(subValues),
   m_nodeOnBoundary(nodeOnBoundary),
   m_spatialDim(spatialDim),
   m_addCorners(addCorners),
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineActiveNodes()
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineStartingSubdomain()
{
  int numSub = m_subNodes.size();
  int numSubScan;
  Teuchos::scan<int, int> (*m_Comm, Teuchos::REDUCE_SUM, 1, &numSub, 
			   &numSubScan);
  m_startingSub = numSubScan - numSub;
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineNodeSubdomains(std::vector< std::vector<LO> > & nodeSubdomains) const
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
    rcp(new Map(m_IGO, Teuchos::ArrayView<const GO>(m_nodeGlobalIDs,m_numNodes), 
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineNodeConnectivity
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineAdjacentNodes(LO & numAdjNodes, 
		       std::vector<LO> & adjNodes, 
		       std::vector<bool> & nodeFlag,
		       std::vector< std::vector<LO> > & subDofNodes,
		       LO sub, 
		       LO node) const
{
  const std::vector<LO> & rowBegin = m_subRowBegin[sub];
  const std::vector<LO> & columns = m_subColumns[sub];
  const std::vector<SX> & values = m_subValues[sub];
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineAdjacentSubdomains(std::vector< std::vector<LO> > & nodeSubdomains)
{
  LO numSub = m_subNodes.size();
  m_adjacentSubdomains.resize(numSub);
  LO numEquiv = m_equivClasses.size();
  for (LO i=0; i<numEquiv; i++) {
    LO node = m_equivClasses[i][0];
    const std::vector<LO> & subs = nodeSubdomains[node];
    if (subs.size() == 2) {
      for (LO j=0; j<2; j++) {
	LO localSub = subs[j] - m_startingSub;
	if ((localSub >= 0) && (localSub < numSub)) {
	  if (j == 0) m_adjacentSubdomains[localSub].push_back(subs[1]);
	  else m_adjacentSubdomains[localSub].push_back(subs[0]);
	}
      }
    }
  }
}

template <TemplateArgs1NPU>
bool PartitionOfUnity<TemplateArgs2NPU>::
edgeCandidate(LO equiv,
	      std::vector< std::vector<LO> > & nodeSubdomains) const
{
  bool edgeCandidate = false;
  LO numNodes = m_equivClasses[equiv].size();
  if (numNodes == 1) return false; // reject corner nodes
  LO node = m_equivClasses[equiv][0];
  const std::vector<LO> & subs = nodeSubdomains[node];
  if (m_spatialDim == 2) {
    if (subs.size() == 2) edgeCandidate = true;
  }
  else if (m_spatialDim == 3) {
    if (subs.size() > 2) edgeCandidate = true;
  }
  return edgeCandidate;
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineCornersToAdd(LO equiv, 
		      std::vector<bool> & nodeIsCorner, 
		      std::vector< std::vector<LO> > & nodeConnectivity, 
		      std::vector<LO> & nodeMap,
		      std::vector<LO> & cornerIndices, 
		      std::vector<LO> & edgeIndices)
{
  const std::vector<LO> & nodes = m_equivClasses[equiv];
  LO numNodes = nodes.size();
  std::vector< std::vector<LO> > localConn(numNodes);
  for (LO i=0; i<numNodes; i++) nodeMap[nodes[i]] = i;
  for (LO i=0; i<numNodes; i++) {
    LO node = nodes[i];
    for (size_t j=0; j<nodeConnectivity[node].size(); j++) {
      LO node2 = nodeMap[nodeConnectivity[node][j]];
      if ((node2 != -1) && (node2 != i)) localConn[i].push_back(node2);
    }
  }
  for (LO i=0; i<numNodes; i++) {
    nodeMap[nodes[i]] = -1;
    bool cornerCandidate = false;
    if (m_nodeOnBoundary != 0) {
      if (m_nodeOnBoundary[nodes[i]] == 1) cornerCandidate = true;
    }
    else {
      if (localConn[i].size() < 2) cornerCandidate = true;
    }
    if (cornerCandidate == true) {
      bool adjacentToCorner = false;
      LO node = nodes[i];
      for (size_t j=0; j<nodeConnectivity[node].size(); j++) {
	LO node2 = nodeConnectivity[node][j];
	if (nodeIsCorner[node2] == true) adjacentToCorner = true;
      }
      if (adjacentToCorner == false) cornerIndices.push_back(i);
      else edgeIndices.push_back(i);
    }
    else edgeIndices.push_back(i);
  }
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
addEdgeEndpoints(std::vector< std::vector<LO> > & nodeSubdomains, 
		 std::vector< std::vector<LO> > & nodeConnectivity)
{
  if (m_addCorners == false) return;
  std::vector<bool> nodeIsCorner(m_numNodes, false);
  LO numEquiv = m_equivClasses.size();
  for (LO i=0; i<numEquiv; i++) {
    if (m_equivClasses[i].size() == 1) {
      LO node = m_equivClasses[i][0];
      nodeIsCorner[node] = true;
    }
  }
  std::vector<LO> newCorners, nodeMap(m_numNodes, -1);
  for (LO i=0; i<numEquiv; i++) {
    if (edgeCandidate(i, nodeSubdomains) == true) {
      std::vector<LO> cornerIndices, edgeIndices;
      determineCornersToAdd(i, nodeIsCorner, nodeConnectivity, nodeMap,
			    cornerIndices, edgeIndices);
      if (cornerIndices.size() > 0) {
	for (size_t j=0; j<cornerIndices.size(); j++) {
	  newCorners.push_back(m_equivClasses[i][cornerIndices[j]]);
	}
	for (size_t j=0; j<edgeIndices.size(); j++) {
	  m_equivClasses[i][j] = m_equivClasses[i][edgeIndices[j]];
	}
	m_equivClasses[i].resize(edgeIndices.size());
      }
    }
  }
  for (size_t i=0; i<newCorners.size(); i++) {
    std::vector<LO> node(1, newCorners[i]);
    m_equivClasses.push_back(node);
  }
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineSubdomainEquivClasses(std::vector< std::vector<LO> > & nodeSubdomains)
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineEquivClassGlobalIDs()
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineEquivClassSubdomains(std::vector< std::vector<LO> > & nodeSubdomains)
{
  LO numEquiv = m_equivClasses.size();
  m_equivClassSubdomains.resize(numEquiv);
  for (LO i=0; i<numEquiv; i++) {
    LO firstNode = m_equivClasses[i][0];
    m_equivClassSubdomains[i] = nodeSubdomains[firstNode];
  }
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineEquivalenceClasses()
{
  std::vector< std::vector<LO> > nodeSubdomains;
  determineNodeSubdomains(nodeSubdomains);
  std::vector< std::vector<LO> > nodeConnectivity;
  determineNodeConnectivity(nodeConnectivity);
  std::vector<int> component(m_numNodes, 0);
  equivalenceClasses(nodeSubdomains, nodeConnectivity, component);
  determineAdjacentSubdomains(nodeSubdomains);
  determineEquivClassComponents(nodeConnectivity, component);
  equivalenceClasses(nodeSubdomains, nodeConnectivity, component);
  addEdgeEndpoints(nodeSubdomains, nodeConnectivity);
  determineEquivClassGlobalIDs();
  determineEquivClassSubdomains(nodeSubdomains);
  determineSubdomainEquivClasses(nodeSubdomains);
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
equivalenceClasses(std::vector< std::vector<LO> > & nodeSubdomains,
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
sortByGlobalIDs(std::vector<LO> & localIDs) const
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineEquivClassComponents
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
determineComponents(LO *A1, 
		    LO *A2, 
		    LO N, 
		    LO* & component)
{
  LO i;
  component = new LO[N];
  for (i=0; i<N; i++) component[i] = -1;
  componentsFunction(A1, A2, N, component);
}

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
componentsFunction(LO *A1, 
		   LO *A2, 
		   LO N, 
		   LO* component)
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

template <TemplateArgs1NPU>
void PartitionOfUnity<TemplateArgs2NPU>::
depthFirstSearch(const LO v, 
		 const LO comp_num, 
		 LO* component,
		 LO* A1, 
		 LO* A2)
{
  LO i, adj_vertex;
  component[v] = comp_num;
  for (i=A2[v]; i<A2[v+1]; i++) {
    adj_vertex = A1[i];
    if (component[adj_vertex] == -1) 
      depthFirstSearch(adj_vertex, comp_num, component, A1, A2);
  }
}

} // namespace bddc

#undef TemplateArgs1NPU
#undef TemplateArgs2NPU

#endif // PARTITIONOFUNITYBDDC_H
  
