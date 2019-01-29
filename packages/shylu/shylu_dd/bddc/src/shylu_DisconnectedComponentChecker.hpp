
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

#ifndef BDDC_DISCONNECTEDCOMPONENTCHECKER_H
#define BDDC_DISCONNECTEDCOMPONENTCHECKER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "shylu_UtilBDDC.hpp"
#include "shylu_errorBDDC.hpp"

namespace bddc {

template <class LO, class GO, class SX> 
class DisconnectedComponentChecker
{
 public:

  DisconnectedComponentChecker()
  {
  }

  DisconnectedComponentChecker(const LO numNodes,
			       const LO* nodeBegin,
			       const std::vector< std::vector<LO> > & subNodes,
			       LO** subRowBegin,
			       LO** subColumns,
			       SX** subValues)
  {
    std::vector< std::vector<LO> > subNodeBegin;
    getSubNodeBegin(subNodes, nodeBegin, subNodeBegin);
    std::vector< std::vector<int> > subNodeComponents;
    identifyComponents(numNodes, nodeBegin, subNodes, subNodeBegin, 
		       subRowBegin, subColumns, subValues, subNodeComponents);
    const LO numDofs = nodeBegin[numNodes];
    splitSubdomainData(numNodes, numDofs, subNodes, subNodeBegin, subRowBegin, 
		       subColumns, subValues, subNodeComponents);
  }
  
  ~DisconnectedComponentChecker()
  {
    delete [] m_subRowBeginPtr;
    delete [] m_subColumnsPtr;
    delete [] m_subValuesPtr;
  }

  const LO getNumSubs() const
  {
    return m_numSubNew;
  }

  void adjustProblemData(std::vector< std::vector<LO> > & subNodes,
			 LO** & subRowBegin,
			 LO** & subColumns,
			 SX** & subValues)
  {
    if (m_numMatNew == 0) return;
    subNodes = m_subNodesNew;
    subRowBegin = m_subRowBeginPtr;
    subColumns = m_subColumnsPtr;
    subValues = m_subValuesPtr;
  }

 private: // variables
  LO m_numSubNew{0}, m_numMatNew{0};
  LO **m_subRowBeginPtr{0}, **m_subColumnsPtr{0};
  SX **m_subValuesPtr{0};
  std::vector< std::vector<LO> > m_subNodesNew, m_subRowBeginNew,
    m_subColumnsNew;
  std::vector< std::vector<SX> > m_subValuesNew;

 private: // functions

  void getColNodes(const LO numNodes, 
		   const LO* nodeBegin, 
		   std::vector<LO> & colNodes)
  {
    for (LO i=0; i<numNodes; i++) {
      for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	colNodes[j] = i;
      }
    }
  }

  void getSubNodeBegin(const std::vector< std::vector<LO> > & subNodes, 
		       const LO* nodeBegin, 
		       std::vector< std::vector<LO> > & subNodeBegin)
  {
    const LO numSub = subNodes.size();
    subNodeBegin.resize(numSub);
    for (LO j=0; j<numSub; j++) {
      const LO numNodes = subNodes[j].size();
      subNodeBegin[j].resize(numNodes+1, 0);
      for (LO i=0; i<numNodes; i++) {
	const LO node = subNodes[j][i];
	const LO numDofNode = nodeBegin[node+1] - nodeBegin[node];
	subNodeBegin[j][i+1] = subNodeBegin[j][i] + numDofNode;
      }
    }
  }

  void identifyComponents(const LO numNodes,
			  const LO* nodeBegin,
			  const std::vector< std::vector<LO> > & subNodes,
			  const std::vector< std::vector<LO> > & subNodeBegin,
			  LO** subRowBegin,
			  LO** subColumns,
			  SX** subValues,
			  std::vector< std::vector<int> > & subNodeComponents)
  {
    // Note: we place nodes without any degrees of freedom or not
    //       connected to any other nodes in component 0
    const LO numSub = subNodes.size();
    const LO numDof = nodeBegin[numNodes];
    std::vector<LO> nodeMap(numNodes, -1), activeNodes(numNodes),
      colNodes(numDof);
    std::vector<bool> nodeFlag(numNodes, false);
    subNodeComponents.resize(numSub);
    for (LO i=0; i<numSub; i++) {
      std::vector<int> nodeConn, nodeConnPtr;
      const LO numNodes = subNodes[i].size();
      determineNodalConnectivity
	(numNodes, subNodeBegin[i].data(), subRowBegin[i], subColumns[i],
	 activeNodes, colNodes, nodeFlag, nodeConn, nodeConnPtr);
      // make adjustments to nodeConn and nodeConnPtr to not include 
      // any disconnected nodes
      LO numConnected(0);
      for (LO j=0; j<numNodes; j++) {
	const LO numAdj = nodeConnPtr[j+1] - nodeConnPtr[j];
	if (numAdj > 0) nodeMap[j] = numConnected++;
      }
      std::vector<int> nodeConnPtrNew(numConnected+1, 0);
      numConnected = 0;
      LO numTerms(0);
      for (LO j=0; j<numNodes; j++) {
	const LO numAdj = nodeConnPtr[j+1] - nodeConnPtr[j];
	if (numAdj > 0) {
	  for (LO k=nodeConnPtr[j]; k<nodeConnPtr[j+1]; k++) {
	    const LO node = nodeMap[nodeConn[k]];
	    BDDC_TEST_FOR_EXCEPTION(node == -1, std::runtime_error, 
				    "invalid node number");
	    nodeConn[k] = node;
	  }
	  numTerms += numAdj;
	  nodeConnPtrNew[++numConnected] = numTerms;
	}
      }
      int *components(nullptr);
      UtilBDDC<double,double>::determineComponents
	(nodeConn.data(), nodeConnPtrNew.data(), numConnected, components);
      // Note: any inactive nodes are placed in the the first component
      subNodeComponents[i].resize(numNodes, 0);
      numConnected = 0;
      for (LO j=0; j<numNodes; j++) {
	if (nodeMap[j] != -1) {
	  subNodeComponents[i][j] = components[numConnected++];
	}
	nodeMap[j] = -1;
      }
      delete [] components;
    }
  }

  void determineNodalConnectivity(const LO numNodes,
				  const LO* nodeBegin,
				  const LO* rowBegin,
				  const LO* columns,
				  std::vector<LO> & activeNodes,
				  std::vector<LO> & colNodes,
				  std::vector<bool> & nodeFlag,
				  std::vector<int> & nodeConn,
				  std::vector<int> & nodeConnPtr)
  {
    getColNodes(numNodes, nodeBegin, colNodes);
    LO row(0), numTerms(0);
    nodeConnPtr.resize(numNodes+1, 0);
    for (LO i=0; i<numNodes; i++) {
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      LO numAdjacent(0);
      for (LO j=0; j<numDofNode; j++) {
	for (LO k=rowBegin[row]; k<rowBegin[row+1]; k++) {
	  const LO col = columns[k];
	  const LO node2 = colNodes[col];
	  if (nodeFlag[node2] == false) {
	    nodeFlag[node2] = true;
	    activeNodes[numAdjacent++] = node2;
	  }
	}
	row++;
      }
      numTerms += numAdjacent;
      nodeConnPtr[i+1] = numTerms;
      for (LO j=0; j<numAdjacent; j++) nodeFlag[activeNodes[j]] = false;
    }
    BDDC_TEST_FOR_EXCEPTION(row != nodeBegin[numNodes], std::runtime_error, 
			    "inconsistent value of counter row");
    nodeConn.resize(numTerms);
    row = numTerms = 0;
    for (LO i=0; i<numNodes; i++) {
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      LO numAdjacent(0);
      for (LO j=0; j<numDofNode; j++) {
	for (LO k=rowBegin[row]; k<rowBegin[row+1]; k++) {
	  const LO col = columns[k];
	  const LO node2 = colNodes[col];
	  if (nodeFlag[node2] == false) {
	    nodeFlag[node2] = true;
	    activeNodes[numAdjacent++] = node2;
	    nodeConn[numTerms++] = node2;
	  }
	}
	row++;
      }
      for (LO j=0; j<numAdjacent; j++) nodeFlag[activeNodes[j]] = false;
    }
  }

  int getNumComponents(const std::vector<int> & components)
  {
    int maxComponent(0);
    for (size_t i=0; i<components.size(); i++) {
      if (components[i] > maxComponent) {
	maxComponent = components[i];
      }
    }
    return maxComponent + 1;
  }

  void splitSubdomainData
    (const LO numNodes,
     const LO numDofs,
     const std::vector<std::vector<LO> > & subNodes,
     const std::vector<std::vector<LO> > & subNodeBegin,
     LO** subRowBegin,
     LO** subColumns,
     SX** subValues,
     const std::vector<std::vector<int> > & subNodeComponents)
  {
    const LO numSub = subNodes.size();
    m_numSubNew = m_numMatNew = 0;
    for (LO i=0; i<numSub; i++) {
      const LO numComps = getNumComponents(subNodeComponents[i]); 
      m_numSubNew += numComps;
      if (numComps > 1) m_numMatNew += numComps;
    }
    if (m_numMatNew == 0) return;
    m_subRowBeginPtr = new LO*[m_numSubNew];
    m_subColumnsPtr = new LO*[m_numSubNew];
    m_subValuesPtr = new SX*[m_numSubNew];
    m_subNodesNew.resize(m_numSubNew);
    m_subRowBeginNew.resize(m_numMatNew);
    m_subColumnsNew.resize(m_numMatNew);
    m_subValuesNew.resize(m_numMatNew);
    m_numSubNew = m_numMatNew = 0;
    std::vector<LO> nodeMap(numNodes, -1), colNodes(numDofs), colLocal(numDofs);
    for (LO i=0; i<numSub; i++) {
      const LO numComps = getNumComponents(subNodeComponents[i]);
      if (numComps == 1) {
	m_subNodesNew[m_numSubNew] = subNodes[i];
	m_subRowBeginPtr[m_numSubNew] = subRowBegin[i];
	m_subColumnsPtr[m_numSubNew] = subColumns[i];
	m_subValuesPtr[m_numSubNew++] = subValues[i];
      }
      else {
	splitSubdomainData(numComps, subNodes[i], subNodeBegin[i].data(), 
			   subRowBegin[i], subColumns[i], subValues[i], 
			   subNodeComponents[i], colNodes, colLocal,
			   nodeMap);
      }
    }
  }

  void splitSubdomainData(const LO numComps,
			  const std::vector<LO> & subNodes, 
			  const LO* nodeBegin, 
			  const LO* rowBegin, 
			  const LO* columns, 
			  const SX* values,
			  const std::vector<int> & components,
			  std::vector<LO> & colNodes,
			  std::vector<LO> & colLocal,
			  std::vector<LO> & nodeMap)
  {
    std::vector< std::vector<LO> > nodesComp(numComps);
    std::vector< std::vector<LO> > nodeBeginComp(numComps);
    std::vector<LO> count(numComps, 0), countDof(numComps, 0);
    const LO numNodes = subNodes.size();
    for (LO i=0; i<numNodes; i++) {
      const LO comp = components[i];
      count[comp]++;
      for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	colNodes[j] = i;
	colLocal[j] = j - nodeBegin[i];
	countDof[comp]++;
      }
    }
    std::vector< std::vector<LO> > dofsComp(numComps);
    for (LO i=0; i<numComps; i++) {
      nodesComp[i].resize(count[i]);
      nodeBeginComp[i].resize(count[i]+1, 0);
      dofsComp[i].resize(countDof[i]);
      count[i] = countDof[i] = 0;
    }
    for (LO i=0; i<numNodes; i++) {
      const LO comp = components[i];
      nodesComp[comp][count[comp]] = i;
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      const LO index = count[comp]++;
      nodeBeginComp[comp][index+1] = nodeBeginComp[comp][index] + 
	                             numDofNode;
      for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	dofsComp[comp][countDof[comp]++] = j;
      }
    }
    for (LO i=0; i<numComps; i++) {
      std::vector<LO> & nodes = nodesComp[i];
      std::vector<LO> & rowBeginComp = m_subRowBeginNew[m_numMatNew];
      const LO numNode = nodes.size();
      for (LO j=0; j<numNode; j++) nodeMap[nodes[j]] = j;
      const LO numRows = dofsComp[i].size();
      rowBeginComp.resize(numRows+1, 0);
      LO rowComp(0), numTerms(0);
      for (LO j=0; j<numNode; j++) {
	const LO node = nodes[j];
	const LO numDofNode = nodeBegin[node+1] - nodeBegin[node];
	for (LO k=0; k<numDofNode; k++) {
	  const LO rowSource = nodeBegin[node] + k;
	  for (LO m=rowBegin[rowSource]; m<rowBegin[rowSource+1]; m++) {
	    const LO node2 = nodeMap[colNodes[columns[m]]];
	    BDDC_TEST_FOR_EXCEPTION(node2 == -1, std::runtime_error, 
				    "node2 in another component");
	    numTerms++;
	  }
	  rowBeginComp[rowComp+1] = numTerms;
	  rowComp++;
	}
      }
      std::vector<LO> & columnsComp = m_subColumnsNew[m_numMatNew];
      std::vector<SX> & valuesComp = m_subValuesNew[m_numMatNew];
      columnsComp.resize(numTerms);
      valuesComp.resize(numTerms);
      rowComp = numTerms = 0;
      for (LO j=0; j<numNode; j++) {
	const LO node = nodes[j];
	const LO numDofNode = nodeBegin[node+1] - nodeBegin[node];
	for (LO k=0; k<numDofNode; k++) {
	  const LO rowSource = nodeBegin[node] + k;
	  for (LO m=rowBegin[rowSource]; m<rowBegin[rowSource+1]; m++) {
	    const LO colSource = columns[m];
	    const LO node2 = nodeMap[colNodes[colSource]];
	    const LO col = nodeBeginComp[i][node2]+colLocal[colSource];
	    columnsComp[numTerms] = col;
	    valuesComp[numTerms++] = values[m];
	  }
	}
      }
      for (LO j=0; j<numNode; j++) {
	nodeMap[nodes[j]] = -1;
	nodes[j] = subNodes[nodes[j]]; // to get back to higher level numbers
      }
      m_subNodesNew[m_numSubNew] = nodes;
      m_subRowBeginPtr[m_numSubNew] = m_subRowBeginNew[m_numMatNew].data();
      m_subColumnsPtr[m_numSubNew] = m_subColumnsNew[m_numMatNew].data();
      m_subValuesPtr[m_numSubNew] = m_subValuesNew[m_numMatNew].data();
      m_numSubNew++;
      m_numMatNew++;
    }
  }

  };

} // namespace bddc

#endif // BDDC_DISCONNECTEDCOMPONENTCHECKER_H
  
