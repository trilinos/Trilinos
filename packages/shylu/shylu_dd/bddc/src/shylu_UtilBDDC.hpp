
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

#ifndef BDDC_UTIL_H
#define BDDC_UTIL_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <vector>

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Distributor.hpp"
#include "shylu_PComm.hpp"
#include "shylu_errorBDDC.hpp"

namespace bddc {

  template <class LO, class GO>
    void unionData(const LO numRows,
		   const GO* rowGlobalIDs,
		   const int* rowSend,
		   const RCP<Tpetra::Distributor> & distributor,
		   const std::vector< std::vector<GO> > & inputData,
		   std::vector< std::vector<GO> > & unionData)
  {
    // load data to send to other processors
    const LO numSends = distributor->getNumSends();
    std::vector<size_t> sendCount(numSends, 0);
    for (LO i=0; i<numRows; i++) {
      const LO localProc = rowSend[i];
      if (localProc == -1) continue;
      sendCount[localProc] += 2*inputData[i].size();
    }
    std::vector<LO> sendStart(numSends+1, 0);
    for (LO i=0; i<numSends; i++) {
      sendStart[i+1] = sendStart[i] + sendCount[i];
      sendCount[i] = 0;
    }
    std::vector<GO> sendData(sendStart[numSends]);
    unionData.resize(numRows);
    for (LO i=0; i<numRows; i++) {
      const LO length = inputData[i].size();
      const LO localProc = rowSend[i];
      if (localProc == -1) {
	unionData[i] = inputData[i];
      }
      else {
	const GO globalID = rowGlobalIDs[i];
	for (LO j=0; j<length; j++) {
	  const LO index = sendStart[localProc] + sendCount[localProc];
	  sendData[index+0] = globalID;
	  sendData[index+1] = inputData[i][j];
	  sendCount[localProc] += 2;
	}
      }
    }
    // send data sizes to owners
    const LO numRecvs = distributor->getNumReceives();
    std::vector<size_t> recvCount(numRecvs);
    distributor->doPostsAndWaits
      (Teuchos::ArrayView<const size_t>(sendCount), 1,
       Teuchos::ArrayView<size_t>(recvCount));
    LO numTermsRecv(0);
    for (LO i=0; i<numRecvs; i++) numTermsRecv += recvCount[i];
    std::vector<GO> recvData(numTermsRecv);
    // send data to owners
    distributor->doPostsAndWaits
      (Teuchos::ArrayView<const GO>(sendData),
       Teuchos::ArrayView<const size_t>(sendCount),
       Teuchos::ArrayView<GO>(recvData),
       Teuchos::ArrayView<const size_t>(recvCount));
    // process received data
    std::map<GO,LO> rowMap;
    for (LO i=0; i<numRows; i++) {
      rowMap.emplace(rowGlobalIDs[i], i);
    }
    LO index(0);
    std::vector< std::vector<LO> > sendBackRows(numRecvs);
    for (LO i=0; i<numRecvs; i++) {
      std::set<LO> nodeSet;
      const LO numTerms = recvCount[i]/2;
      for (LO j=0; j<numTerms; j++) {
	const LO globalID = recvData[index++];
	const GO value = recvData[index++];
	auto iter = rowMap.find(globalID);
	BDDC_TEST_FOR_EXCEPTION(iter == rowMap.end(), std::runtime_error, 
				"globalID not found");
	const LO row = iter->second;
	unionData[row].push_back(value);
	auto iter2 = nodeSet.find(row);
	if (iter2 == nodeSet.end()) {
	  nodeSet.emplace(row);
	  sendBackRows[i].push_back(row);
	}
      }
    }
    // reload recvData with data to send back to sending processors
    recvCount.assign(numRecvs, 0);
    LO numTerms(0);
    for (LO i=0; i<numRecvs; i++) {
      for (size_t j=0; j<sendBackRows[i].size(); j++) {
	const LO row = sendBackRows[i][j];
	recvCount[i] += 2*unionData[row].size();
      }
      numTerms += recvCount[i];
    }
    recvData.resize(numTerms);
    index = 0;
    for (LO i=0; i<numRecvs; i++) {
      for (size_t j=0; j<sendBackRows[i].size(); j++) {
	const LO row = sendBackRows[i][j];
	for (size_t k=0; k<unionData[row].size(); k++) {
	  const GO value = unionData[row][k];
	  recvData[index++] = rowGlobalIDs[row];
	  recvData[index++] = value;
	}
      }
    }    
    // send data sizes back to senders
    distributor->doReversePostsAndWaits
      (Teuchos::ArrayView<const size_t>(recvCount), 1,
       Teuchos::ArrayView<size_t>(sendCount));
    LO numTermsSend(0);
    for (LO i=0; i<numSends; i++) numTermsSend += sendCount[i];
    sendData.resize(numTermsSend);
    // send data back to senders
    distributor->doReversePostsAndWaits
      (Teuchos::ArrayView<const GO>(recvData),
       Teuchos::ArrayView<const size_t>(recvCount),
       Teuchos::ArrayView<GO>(sendData),
       Teuchos::ArrayView<const size_t>(sendCount));
    // process sent back data
    index = 0;
    for (LO i=0; i<numSends; i++) {
      for (size_t j=0; j<sendCount[i]/2; j++) {
	const GO globalID = sendData[index++];
	auto iter = rowMap.find(globalID);
	BDDC_TEST_FOR_EXCEPTION(iter == rowMap.end(), std::runtime_error, 
				"globalID not found");
	const LO row = iter->second;
	BDDC_TEST_FOR_EXCEPTION(rowSend[row] == -1, std::runtime_error, 
				"rowSend[row] is invalid");
	const GO value = sendData[index++];
	unionData[row].push_back(value);
      }
    }    
  }

  template <class LO>
  void constructDistributor(const LO numNodes,
			    std::vector<int> & nodeSend,
			    RCP<const Teuchos::Comm<int> > & Comm,
			    RCP<Tpetra::Distributor> & distributor)
  {
    // Note: If nodeSend[i] == -1, then the associated node does not send
    //       data to any other processors.
    //       If nodeSend[i] != -1, then the associated node sends data to
    //       processor nodeSend[i]. On exit, the entries in nodeSend not
    //       equal to -1 are adjusted to refer to local indices for the
    //       sending processors of the distributor.
    //       
    LO lengthSend(0);
    BDDC_TEST_FOR_EXCEPTION(LO(nodeSend.size()) != numNodes, 
			    std::runtime_error, "nodeSend size error");
    for (LO i=0; i<numNodes; i++) {
      if (nodeSend[i] != -1) lengthSend++;
    }
    std::vector<int> sendToProcs(lengthSend);
    lengthSend = 0;
    for (LO i=0; i<numNodes; i++) {
      if (nodeSend[i] != -1) {
	sendToProcs[lengthSend++] = nodeSend[i];
      }
    }
    std::sort(sendToProcs.begin(), sendToProcs.end());
    auto iter = std::unique(sendToProcs.begin(), sendToProcs.end());
    sendToProcs.erase(iter, sendToProcs.end());
    sendToProcs.shrink_to_fit();
    distributor = rcp( new Tpetra::Distributor(Comm) );
    distributor->createFromSends
      (Teuchos::ArrayView<const int>(sendToProcs));
    // adjust m_nodeSend to refer locally to sendToProcs in m_distributor
    std::map<int,LO> procMap;
    for (size_t i=0; i<sendToProcs.size(); i++) {
      procMap.emplace(sendToProcs[i], i);
    }
    // convert nodeSend to local proc numbers for distributor
    for (LO i=0; i<numNodes; i++) {
      if (nodeSend[i] != -1) {
	auto iter = procMap.find(nodeSend[i]);
	BDDC_TEST_FOR_EXCEPTION(iter == procMap.end(), std::runtime_error, 
				"nodeSend[i] not found");
	nodeSend[i] = iter->second;
      }
    }
  }

  template <class LO,class GO>
   void getNodeSend(const LO numNode, 
		    const GO* nodeGlobalIDs,
		    MPI_Comm mpiComm,
		    std::vector<int> & nodeSend)
  {
    typedef Tpetra::Map<LO,GO>                                 Map;
    typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
    typedef Tpetra::Export<LO,GO>                              Export;
    typedef Tpetra::Import<LO,GO>                              Import;
    RCP<const Teuchos::Comm<int> > Comm = 
      rcp( new Teuchos::MpiComm<int>(mpiComm) );
    const int myPID = Comm->getRank();
    Tpetra::global_size_t IGO = 
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    std::vector<GO> proc(1, myPID);
    RCP<const Map> nodeMap = 
      rcp(new Map(IGO, Teuchos::ArrayView<const GO>(nodeGlobalIDs, numNode), 
		  0, Comm));
    RCP<const Map> colMap =
      rcp(new Map(IGO, Teuchos::ArrayView<const GO>(proc), 0, Comm) );
    CrsGraph nodeProcs(nodeMap, colMap, 1, Tpetra::StaticProfile);
    std::vector<LO> zero(1, 0);
    for (LO i=0; i<numNode; i++) {
      nodeProcs.insertLocalIndices(i, Teuchos::ArrayView<LO>(zero));
    }
    RCP<const Map> nodeMap1to1 = Tpetra::createOneToOne<LO,GO>(nodeMap);
    nodeProcs.fillComplete(colMap, nodeMap1to1);
    CrsGraph nodeProcs1to1(nodeMap1to1, 0);
    Export Exporter(nodeMap, nodeMap1to1);
    nodeProcs1to1.doExport(nodeProcs, Exporter, Tpetra::ADD);
    nodeProcs1to1.fillComplete(colMap, nodeMap1to1);
    Import Importer(nodeMap1to1, nodeMap);
    CrsGraph nodeProcsAll(nodeMap, 0);
    nodeProcsAll.doImport(nodeProcs1to1, Importer, Tpetra::INSERT);
    nodeProcsAll.fillComplete(colMap, nodeMap1to1);
    RCP<const Map> colMap2 = nodeProcsAll.getColMap();
    nodeSend.resize(numNode);
    Teuchos::ArrayView<const LO> indices;
    for (LO i=0; i<numNode; i++) {
      nodeProcsAll.getLocalRowView(i, indices);
      int minProc = myPID;
      for (int j=0; j<indices.size(); j++) {
	const int proc = colMap2->getGlobalElement(indices[j]);
	if (proc < minProc) minProc = proc;
      }
      nodeSend[i] = -1;
      if (minProc < myPID) nodeSend[i] = minProc;
    }
  }

  template <class LO,class GO>
   void determineNodeProcs(const LO numNode, 
			   const GO* nodeGlobalIDs,
			   MPI_Comm mpiComm,
			   std::vector<int> & adjProcs,
			   std::vector< std::vector<int> > & nodeProcsVec)
  {
    typedef Tpetra::Map<LO,GO>                                 Map;
    typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
    typedef Tpetra::Export<LO,GO>                              Export;
    typedef Tpetra::Import<LO,GO>                              Import;
    int myPID;
    MPI_Comm_rank(mpiComm, &myPID);
    RCP<const Teuchos::Comm<int> > Comm = 
      rcp( new Teuchos::MpiComm<int>(mpiComm) );
    Tpetra::global_size_t IGO = 
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    std::vector<GO> proc(1, myPID);
    RCP<const Map> nodeMap = 
      rcp(new Map(IGO, Teuchos::ArrayView<const GO>(nodeGlobalIDs, numNode), 
		  0, Comm));
    RCP<const Map> colMap =
      rcp(new Map(IGO, Teuchos::ArrayView<const GO>(proc), 0, Comm) );
    CrsGraph nodeProcs(nodeMap, colMap, 1, Tpetra::StaticProfile);
    std::vector<LO> zero(1, 0);
    for (LO i=0; i<numNode; i++) {
      nodeProcs.insertLocalIndices(i, Teuchos::ArrayView<LO>(zero));
    }
    RCP<const Map> nodeMap1to1 = Tpetra::createOneToOne<LO,GO>(nodeMap);
    nodeProcs.fillComplete(colMap, nodeMap1to1);
    CrsGraph nodeProcs1to1(nodeMap1to1, 0);
    Export Exporter(nodeMap, nodeMap1to1);
    nodeProcs1to1.doExport(nodeProcs, Exporter, Tpetra::ADD);
    nodeProcs1to1.fillComplete(colMap, nodeMap1to1);
    Import Importer(nodeMap1to1, nodeMap);
    CrsGraph nodeProcsAll(nodeMap, 0);
    nodeProcsAll.doImport(nodeProcs1to1, Importer, Tpetra::INSERT);
    nodeProcsAll.fillComplete(colMap, nodeMap1to1);
    Teuchos::ArrayView<const LO> Indices;
    RCP<const Map> colMap2 = nodeProcsAll.getColMap();
    const LO numProcAll = colMap2->getNodeNumElements();
    adjProcs.resize(std::max(0, numProcAll-1)); // don't include self
    LO numOtherProc(0);
    for (LO i=0; i<numProcAll; i++) {
      const GO procGID = colMap2->getGlobalElement(i);
      if (procGID != myPID) {
	adjProcs[numOtherProc++] = procGID;
      }
    }
    BDDC_TEST_FOR_EXCEPTION(numOtherProc != std::max(0, numProcAll-1), 
			    std::runtime_error, "numOtherProc error");
    std::sort(adjProcs.begin(), adjProcs.end()); // sorted important later
    std::map<int,int> procMap;
    for (LO i=0; i<numOtherProc; i++) {
      procMap.insert(std::make_pair(adjProcs[i], i)); 
    }
    nodeProcsVec.resize(numNode);
    for (LO i=0; i<numNode; i++) {
      nodeProcsAll.getLocalRowView(i, Indices);
      LO numProc = Indices.size();
      BDDC_TEST_FOR_EXCEPTION(numProc < 1, std::runtime_error, 
			      "numProc less than 1");
      for (LO j=0; j<numProc; j++) {
	const GO procGID = colMap2->getGlobalElement(Indices[j]);
	auto iter = procMap.find(procGID);
	if (iter != procMap.end()) {
	  nodeProcsVec[i].push_back(iter->second); 
	}
      }
    }
  }

  template <class LO>
  void getNodeData(const std::vector<LO> & subNodes, 
	      const std::vector<LO> & subNodeBegin, 
	      std::vector<LO> & nodeMap, 
	      LO & numActiveNodes,
	      std::vector<LO> & colNodes)
  {
    numActiveNodes = 0;
    for (size_t i=0; i<subNodes.size(); i++) {
      const LO numDofNode = subNodeBegin[i+1] - subNodeBegin[i];
      if (numDofNode > 0) {
	nodeMap[i] = numActiveNodes++;
      }
      for (LO j=subNodeBegin[i]; j<subNodeBegin[i+1]; j++) colNodes[j] = i;
    }
  }

  template <class LO>
  void determineNodalConnectivity(const std::vector<LO> & subNodes,
				  const std::vector<LO> & subNodeBegin,
				  const LO* rowBegin,
				  const LO* columns,
				  std::vector<LO> & colNodes,
				  std::vector<LO> & nodeMap,
				  std::vector<bool> & nodeFlag,
				  std::vector<LO> & activeNodes,
				  std::vector<int> & nodeConn,
				  std::vector<int> & nodeConnPtr)
  {
    LO numActiveNodes;
    getNodeData(subNodes, subNodeBegin, nodeMap, numActiveNodes, colNodes);
    LO row(0), numTerms(0);
    nodeConnPtr.resize(numActiveNodes+1, 0);
    numActiveNodes = 0;
    for (size_t i=0; i<subNodes.size(); i++) {
      if (nodeMap[i] == -1) continue;
      const LO numDofNode = subNodeBegin[i+1] - subNodeBegin[i];
      BDDC_TEST_FOR_EXCEPTION(numDofNode <= 0, std::runtime_error, 
			      "numDofNode must be greater than zero");
      LO numAdjacent(0);
      for (LO j=0; j<numDofNode; j++) {
	for (LO k=rowBegin[row]; k<rowBegin[row+1]; k++) {
	  const LO col = columns[k];
	  const LO node2 = nodeMap[colNodes[col]];
	  if (node2 != -1) {
	    if (nodeFlag[node2] == false) {
	      nodeFlag[node2] = true;
	      activeNodes[numAdjacent++] = node2;
	    }
	  }
	}
	row++;
      }
      numTerms += numAdjacent;
      numActiveNodes++;
      nodeConnPtr[numActiveNodes] = numTerms;
      for (LO j=0; j<numAdjacent; j++) nodeFlag[activeNodes[j]] = false;
    }
    nodeConn.resize(numTerms);
    row = numTerms = 0;
    for (size_t i=0; i<subNodes.size(); i++) {
      if (nodeMap[i] == -1) continue;
      const LO numDofNode = subNodeBegin[i+1] - subNodeBegin[i];
      LO numAdjacent(0);
      for (LO j=0; j<numDofNode; j++) {
	for (LO k=rowBegin[row]; k<rowBegin[row+1]; k++) {
	  const LO col = columns[k];
	  const LO node2 = nodeMap[colNodes[col]];
	  if (node2 != -1) {
	    if (nodeFlag[node2] == false) {
	      nodeFlag[node2] = true;
	      activeNodes[numAdjacent++] = node2;
	      nodeConn[numTerms++] = node2;
	    }
	  }
	}
	row++;
      }
      for (LO j=0; j<numAdjacent; j++) nodeFlag[activeNodes[j]] = false;
    }
    numActiveNodes = 0;
    for (size_t i=0; i<subNodes.size(); i++) {
      nodeMap[i] = -1; 
      const LO numDofNode = subNodeBegin[i+1] - subNodeBegin[i];
      if (numDofNode > 0) activeNodes[numActiveNodes++] = i;
    }
  }

  template <class SX, class SM> 
  class UtilBDDC
{
public:
  UtilBDDC()
  {
  }

  static void determineComponents(int *A1, 
				  int *A2, 
				  int N, 
				  int* & component)
  {
    component = new int[N];
    for (int i=0; i<N; i++) component[i] = -1;
    componentsFunction2(A1, A2, N, component);
  }

  static void componentsFunction(int *A1, 
				 int *A2, 
				 int N, 
				 int* component)
  {
    int comp_num(0);
    for (int i=0; i<N; i++) {
      if (component[i] == -1) {
	depthFirstSearch(i, comp_num, component, A1, A2);
	comp_num++;
      }
    }
  }

  static int getStartingRow(const int* component, 
			    int & startLooking)
  {
    bool foundRow(false);
    int row(-1);
    while (foundRow == false) {
      if (component[startLooking] == -1) {
	foundRow = true;
	row = startLooking;
      }
      startLooking++;
    }
    BDDC_TEST_FOR_EXCEPTION(row == -1, std::runtime_error, "row not found");
    return row;
  }

  static void getRowsToAddToQueue(const int parent,
				  const int* rowBegin,
				  const int* columns,
				  const int* component,
				  int* rowsToAddToQueue,
				  int & numAdded)
  {
    numAdded = 0;
    for (int i=rowBegin[parent]; i<rowBegin[parent+1]; i++) {
      const int col = columns[i];
      if ((component[col] == -1) && (col != parent)) {
	rowsToAddToQueue[numAdded++] = col;
      }
    }
  }

  static void componentsFunction2(int *columns, 
				  int *rowBegin, 
				  int numRows, 
				  int* component)
  {
    std::vector<int> queue(numRows), rowsToAddToQueue(numRows);
    int startLooking(0), numProcessed(0), numComponents(0);
    for (int i=0; i<numRows; i++) component[i] = -1;
    while (numProcessed < numRows) {
      const int parent = getStartingRow(component, startLooking);
      component[parent] = numComponents;
      numProcessed++;
      int numAdded;
      getRowsToAddToQueue(parent, rowBegin, columns, component, 
			  rowsToAddToQueue.data(), numAdded);
      for (int i=0; i<numAdded; i++) {
	const int row = rowsToAddToQueue[i];
	queue[i] = row;
	component[row] = numComponents;
      }
      int qStart(0), qEnd(numAdded);
      while (qEnd > qStart) {
	const int child = queue[qStart];
	component[child] = numComponents;
	numProcessed++;
	getRowsToAddToQueue(child, rowBegin, columns, component, 
			    rowsToAddToQueue.data(), numAdded);
	for (int i=0; i<numAdded; i++) {
	  const int row = rowsToAddToQueue[i];
	  queue[qEnd+i] = row;
	  component[row] = numComponents;
	}
	qStart++;
	qEnd += numAdded;
      }
      numComponents++;
    }
  }

  static void depthFirstSearch(const int v, 
			       const int comp_num, 
			       int* component,
			       int* A1, 
			       int* A2)
  {
    component[v] = comp_num;
    for (int i=A2[v]; i<A2[v+1]; i++) {
      const int adj_vertex = A1[i];
      if (component[adj_vertex] == -1) 
	depthFirstSearch(adj_vertex, comp_num, component, A1, A2);
    }
  }

  static void calculateSchurComplement(std::vector<SX> & A, 
				       std::vector<int> & i1, 
				       std::vector<int> & i2, 
				       std::vector<SX> & Sc)
  {
    std::vector<SX> A21;
    calculateSchurComplement(A, i1, i2, Sc, A21);
  }

  static void calculateSchurComplement(std::vector<SX> & A, 
				       std::vector<int> & i1, 
				       std::vector<int> & i2, 
				       std::vector<SX> & Sc,
				       std::vector<SX> & A21)
  {
    // calculates Schur complement Sc = A11 - A12*inv(A22)*A21
    std::vector<SX> A12, A22;
    loadDenseMatrices(A, i1, i2, Sc, A12, A21, A22);
    Teuchos::BLAS<int, SX>  BLAS;
    Teuchos::LAPACK<int, SX> LAPACK;
    int numRows1 = i1.size();
    int numRows2 = i2.size();
    int INFO(0);
    if (numRows2 > 0) {
      // static condensation of "2" unknowns
      std::vector<int> IPIV(numRows2);
      LAPACK.GETRF(numRows2, numRows2, &A22[0], numRows2, &IPIV[0], &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRF error");
      int NRHS = numRows1;
      LAPACK.GETRS('N', numRows2, NRHS, &A22[0], numRows2, &IPIV[0], 
		   &A21[0], numRows2, &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
      SX ALPHA(-1), BETA(1);
      if (numRows1 > 0) {
	BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, numRows1,
		  numRows2, ALPHA, &A12[0], numRows1, &A21[0], numRows2, 
		  BETA, &Sc[0], numRows1);
      }
    }
  }
  
  static void solveSymmetricGeneralizedEigenProblem(std::vector<SX> & A,
						    std::vector<SX> & B,
						    int N,
						    std::vector<SM> & W)
  {
    int ITYPE = 1; // A*x = lambda*B*x;
    char JOBZ('V'), UPLO('U');
    std::vector<SM> RWORK(std::max(1, 3*N-2));
    Teuchos::LAPACK<int, SX> LAPACK;
    int LWORK(-1), INFO(0);
    SX dumWORK[1];
    LAPACK.HEGV(ITYPE, JOBZ, UPLO, N, &A[0], N, &B[0], N, &W[0], 
		dumWORK, LWORK, &RWORK[0], &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "HEGV error");
    LWORK = int(real(dumWORK[0])+0.1);
    std::vector<SX> WORK(LWORK);
    LAPACK.HEGV(ITYPE, JOBZ, UPLO, N, &A[0], N, &B[0], N, &W[0], 
		&WORK[0], LWORK, &RWORK[0], &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "HEGV error");
  }

  static void solveRealEigenValueProblem(std::vector<SX> & A,
					 int N,
					 std::vector<SM> & W,
					 const char JOBZ)
  {
    if (N == 0) return;
    BDDC_TEST_FOR_EXCEPTION(int(A.size()) != N*N, std::runtime_error, 
			    "A size error");
    char UPLO('U');
    W.resize(N);
    std::vector<SM> RWORK(std::max(1, 3*N));
    int LWORK = std::max(1, 3*N);
    std::vector<SX> WORK(LWORK);
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    LAPACK.HEEV(JOBZ, UPLO, N, &A[0], N, &W[0], &WORK[0], LWORK, 
		&RWORK[0], &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "HEEV error");
  }

  static void calculateRealEigenValues(std::vector<SX> & A, 
				       int N,
				       std::vector<SM> & W)
  {
    char JOBZ('N');
    solveRealEigenValueProblem(A, N, W, JOBZ);
  }

  static void calculateRealEigenValuesAndVectors(std::vector<SX> & A, 
						 int N,
						 std::vector<SM> & W)
  {
    char JOBZ('V');
    solveRealEigenValueProblem(A, N, W, JOBZ);
  }

  static void printIndices(const int numIndices, 
			   const int* indices,
			   const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int j=0; j<numIndices; j++) {
      fout << indices[j]+1 << std::endl;;
    }
    fout.close();
  }

  static void printDenseMatrix(const int numRows,
			       const int numCols,
			       const SX* A,
			       const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int j=0; j<numCols; j++) {
      for (int i=0; i<numRows; i++) {
	fout << i+1 << " ";
	fout << j+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	SX value = A[i+numRows*j];
	fout << real(value);
	if (isComplex(value) == true) {
	  fout << " " << imag(value);
	}
	fout << std::endl;
      }
    }
    fout.close();
  }

  static void printSparseMatrix(const int numRows,
				const int* rowBegin,
				const int* columns,
				const SX* values,
				const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int i=0; i<numRows; i++) {
      for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	fout << i+1 << "  " << columns[j]+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	SX value = values[j];
	fout << real(value);
	if (isComplex(value) == true) {
	  fout << " " << imag(value);
	}
	fout << std::endl;
      }
    }
    fout.close();
  }

  static void printCoords(std::vector<SM> & xCoords,
			  std::vector<SM> & yCoords,
			  std::vector<SM> & zCoords,
			  const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    int numRows = xCoords.size();
    for (int i=0; i<numRows; i++) {
      fout << std::setw(22) << std::setprecision(15);
      fout << xCoords[i] << " "
	   << yCoords[i] << " "
	   << zCoords[i] << std::endl;
    }
    fout.close();
  }

  static void printLocalDofs(std::vector<int> & localDofs,
			     const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    int numRows = localDofs.size();
    for (int i=0; i<numRows; i++) {
      fout << localDofs[i] << std::endl;
    }
    fout.close();
  }

  static void calculateSchurComplement(int numRows1, 
				       int numRows2,
				       int ScOption,
				       std::vector<SX> & A,
				       std::vector<SX> & Sc)
  {
    std::vector<SX> A21;
    calculateSchurComplementHere(numRows1, numRows2, ScOption, A, Sc, A21);
  }

  static void calculateSchurComplement(int numRows1, 
				       int numRows2,
				       int ScOption,
				       std::vector<SX> & A,
				       std::vector<SX> & Sc,
				       std::vector<SX> & ExtensionMatrix)
  {
    calculateSchurComplementHere(numRows1, numRows2, ScOption, A, Sc, 
				 ExtensionMatrix);
  }

  static void calculateSchurComplementHere(int numRows1, 
					   int numRows2,
					   int ScOption,
					   std::vector<SX> & A,
					   std::vector<SX> & Sc,
					   std::vector<SX> & A22invA21)
  {
    bool test = (ScOption == 1) || (ScOption == 2);
    BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			    "ScOption error");
    int N = numRows1 + numRows2;
    if (N == 0) return;
    BDDC_TEST_FOR_EXCEPTION(int(A.size()) != N*N, std::runtime_error, 
			    "invalid dimension of A");
    std::vector<int> i1(numRows1), i2(numRows2);
    for (int i=0; i<numRows1; i++) i1[i] = i;
    for (int i=0; i<numRows2; i++) i2[i] = i + numRows1;
    if (ScOption == 1) {
      calculateSchurComplement(A, i1, i2, Sc, A22invA21);
    }
    else {
      calculateSchurComplement(A, i2, i1, Sc, A22invA21);
    }
  }

  static void loadDenseMatrices(std::vector<SX> & A, 
				std::vector<int> & i1,
				std::vector<int> & i2, 
				std::vector<SX> & A11, 
				std::vector<SX> & A12, 
				std::vector<SX> & A21, 
				std::vector<SX> & A22)
  {
    int numRows1 = i1.size();
    int numRows2 = i2.size();
    A11.resize(numRows1*numRows1);
    A12.resize(numRows1*numRows2);
    A21.resize(numRows2*numRows1);
    A22.resize(numRows2*numRows2);
    int LDA = numRows1 + numRows2;
    for (int i=0; i<numRows1; i++) {
      int row = i1[i];
      for (int j=0; j<numRows1; j++) {
	int col = i1[j];
	A11[i+j*numRows1] = A[row+col*LDA];
      }
      for (int j=0; j<numRows2; j++) {
	int col = i2[j];
	A12[i+j*numRows1] = A[row+col*LDA];
      }
    }
    for (int i=0; i<numRows2; i++) {
      int row = i2[i];
      for (int j=0; j<numRows1; j++) {
	int col = i1[j];
	A21[i+j*numRows2] = A[row+col*LDA];
      }
      for (int j=0; j<numRows2; j++) {
	int col = i2[j];
	A22[i+j*numRows2] = A[row+col*LDA];
      }
    }
  }

  static float imag(float a) {return 0;};
  static double imag(double a) {return 0;};
  static float imag(std::complex<float> a) {return std::imag(a);};
  static double imag(std::complex<double> a) {return std::imag(a);};

  static float real(float a) {return a;};
  static double real(double a) {return a;};
  static float real(std::complex<float> a) {return std::real(a);};
  static double real(std::complex<double> a) {return std::real(a);};

  static bool isComplex(float a) {return false;};
  static bool isComplex(double a) {return false;};
  static bool isComplex(std::complex<float> a) {return true;};
  static bool isComplex(std::complex<double> a) {return true;};

  static float conj(float a) {return a;};
  static double conj(double a) {return a;};
  static std::complex<float> conj(std::complex<float> a) 
  {
    return std::conj(a);
  };
  static std::complex<double> conj(std::complex<double> a) 
  {
    return std::conj(a);
  };

private:

};

} // namespace bddc

#endif // BDDC_UTIL_H
  
