
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

#ifndef BDDC_DOFMANAGER_H
#define BDDC_DOFMANAGER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "shylu_errorBDDC.hpp"

namespace bddc {
  
template <class LO, class GO> 
  class DofManager
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::Map<LO,GO>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
  typedef Tpetra::CrsMatrix<GO,LO,GO>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO>                              Export;
  typedef Tpetra::Import<LO,GO>                              Import;

  DofManager()
  {
  }
   
  ~DofManager()
  {
  }

  static void determineGlobalIDs
    (LO numNodes,
     const GO* nodeGlobalIDs,
     const LO* nodeBegin,
     const LO* localDofs,
     RCP<const Teuchos::Comm<int> > & Comm,
     RCP<const Map> & dofMap,
     RCP<const Map> & dofMap1to1,
     const std::vector<GO>* nodeGlobalIDs1to1=nullptr)
  {
    Tpetra::global_size_t IGO =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const Map> nodeMap = 
      rcp( new Map(IGO, 
		   Teuchos::ArrayView<const GO>(nodeGlobalIDs, numNodes), 
		   0, Comm));
    RCP<const Map> nodeMap1to1;
    if (nodeGlobalIDs1to1 == nullptr) {
      nodeMap1to1 = Tpetra::createOneToOne<LO,GO>(nodeMap);
    }
    else {
      nodeMap1to1 = 
	rcp( new Map(IGO, Teuchos::ArrayView<const GO>(*nodeGlobalIDs1to1), 
		     0, Comm) );
    }
    // The rows and columns of the graph nodeGraphSubdomain are node global IDs 
    // and local degrees of freedom (dofs), respectively, prior to assembly.
    RCP<CrsGraph> nodeGraphSubdomain;
    determineNodalDofGraph(numNodes, nodeBegin, localDofs, Comm, nodeMap, 
			   nodeMap1to1, nodeGraphSubdomain);
    // The graph nodeGraph1to1 is the counterpart of graph nodeGraphSubdomain, 
    // but with uniquely owned rows. Note: it is assumed that the local dofs
    // for each node are identical across all processors
    RCP<CrsGraph> nodeGraph1to1 = rcp( new CrsGraph(nodeMap1to1, 0) );
    Export exporter(nodeMap, nodeMap1to1);
    nodeGraph1to1->doExport(*nodeGraphSubdomain, exporter, Tpetra::ADD);
    RCP<const Map> domainMap = nodeGraphSubdomain->getDomainMap();
    nodeGraph1to1->fillComplete(domainMap, nodeMap1to1);
    GO numDof = nodeGraph1to1->getNodeNumEntries();
    GO numDofSS;
    Teuchos::scan<int, GO> (*Comm, Teuchos::REDUCE_SUM, 1, &numDof, &numDofSS);
    GO baseGID = numDofSS - numDof;
    // The matrix nodeMatrix1to1 is the counterpart of the graph nodeGraph1to1
    // with globalIDs as entries
    Teuchos::ArrayView<const LO> Indices;
    std::vector<GO> Values;
    CrsMatrixGO nodeMatrix1to1(nodeGraph1to1);
    for (size_t i=0; i<nodeMap1to1->getNodeNumElements(); i++) {
      nodeGraph1to1->getLocalRowView(i, Indices);
      Values.resize(Indices.size());
      for (LO j=0; j<Indices.size(); j++) Values[j] = baseGID++;
      nodeMatrix1to1.replaceLocalValues
	(i, Indices, Teuchos::ArrayView<GO>(Values));
    }
    nodeMatrix1to1.fillComplete(domainMap, nodeMap1to1);
    CrsMatrixGO nodeMatrix(nodeMap, 0);
    Import importer(nodeMap1to1, nodeMap);
    nodeMatrix.doImport(nodeMatrix1to1, importer, Tpetra::INSERT);
    nodeMatrix.fillComplete(domainMap, nodeMap1to1);
    extractDofMap(nodeMatrix1to1, dofMap1to1);
    extractDofMap(nodeMatrix, dofMap);
  }

  static void determineNodeMatrices(LO numNodes,
				    const GO* nodeGlobalIDs,
				    const LO* nodeBegin,
				    const LO* localDofs,
				    RCP<const Teuchos::Comm<int> > & Comm,
				    RCP<CrsMatrixGO> & nodeMatrix,
				    RCP<CrsMatrixGO> & nodeMatrix1to1)
  {
    Tpetra::global_size_t IGO =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const Map> nodeMap = 
      rcp( new Map(IGO, 
		   Teuchos::ArrayView<const GO>(nodeGlobalIDs, numNodes), 
		   0, Comm));
    RCP<const Map> nodeMap1to1 = Tpetra::createOneToOne<LO,GO>(nodeMap);
    // The rows and columns of the graph nodeGraphSubdomain are node global IDs 
    // and local degrees of freedom (dofs), respectively, prior to assembly.
    RCP<CrsGraph> nodeGraphSubdomain;
    determineNodalDofGraph(numNodes, nodeBegin, localDofs, Comm, nodeMap, 
			   nodeMap1to1, nodeGraphSubdomain);
    // The graph nodeGraph1to1 is the counterpart of graph nodeGraphSubdomain, 
    // but with uniquely owned rows. Note: it is assumed that the local dofs
    // for each node are identical across all processors
    RCP<CrsGraph> nodeGraph1to1 = rcp( new CrsGraph(nodeMap1to1, 0) );
    Export exporter(nodeMap, nodeMap1to1);
    nodeGraph1to1->doExport(*nodeGraphSubdomain, exporter, Tpetra::ADD);
    RCP<const Map> domainMap = nodeGraphSubdomain->getDomainMap();
    nodeGraph1to1->fillComplete(domainMap, nodeMap1to1);
    GO numDof = nodeGraph1to1->getNodeNumEntries();
    GO numDofSS;
    Teuchos::scan<int, GO> (*Comm, Teuchos::REDUCE_SUM, 1, &numDof, &numDofSS);
    GO baseGID = numDofSS - numDof;
    // The matrix nodeMatrix1to1 is the counterpart of the graph nodeGraph1to1
    // with globalIDs as entries
    Teuchos::ArrayView<const LO> Indices;
    std::vector<GO> Values;
    nodeMatrix1to1 = rcp( new CrsMatrixGO(nodeGraph1to1) );
    for (size_t i=0; i<nodeMap1to1->getNodeNumElements(); i++) {
      nodeGraph1to1->getLocalRowView(i, Indices);
      Values.resize(Indices.size());
      for (LO j=0; j<Indices.size(); j++) Values[j] = baseGID++;
      nodeMatrix1to1->replaceLocalValues
	(i, Indices, Teuchos::ArrayView<GO>(Values));
    }
    nodeMatrix1to1->fillComplete(domainMap, nodeMap1to1);
    nodeMatrix = rcp( new CrsMatrixGO(nodeMap, 0) );
    Import importer(nodeMap1to1, nodeMap);
    nodeMatrix->doImport(*nodeMatrix1to1, importer, Tpetra::INSERT);
    nodeMatrix->fillComplete(domainMap, nodeMap1to1);
  }

  static void determineNodalDofGraph
    (const LO numNode, 
     const LO* nodeBegin,
     const LO* localDofs,
     RCP<const Teuchos::Comm<int> > & Comm,
     RCP<const Map> & NodeMap,
     RCP<const Map> & NodeMap1to1,
     RCP<CrsGraph> & A)
  {
    Tpetra::global_size_t IGO =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::ArrayRCP<size_t> count(numNode);
    for (LO i=0; i<numNode; i++) count[i] = nodeBegin[i+1]-nodeBegin[i];
    LO numDof = nodeBegin[numNode];
    std::vector<GO> uniqueIndicesGO(numDof);
    for (LO i=0; i<numDof; i++) uniqueIndicesGO[i] = localDofs[i];
    determineUniqueIndices(uniqueIndicesGO);
    RCP<const Map> ColMap = 
      rcp( new Map(IGO, Teuchos::ArrayView<GO>(uniqueIndicesGO), 0, Comm) );
    A = rcp( new CrsGraph(NodeMap, ColMap, count) );
    std::vector<LO> localDofsIndex(numDof);
    for (LO i=0; i<numDof; i++) {
      localDofsIndex[i] = ColMap->getLocalElement(localDofs[i]);
    }
    for (LO i=0; i<numNode; i++) {
      if (count[i] > 0) {
	A->insertLocalIndices
	  (i, Teuchos::ArrayView<LO>(&localDofsIndex[nodeBegin[i]], count[i]));
      }
    }
    RCP<const Map> ColMap1to1 = Tpetra::createOneToOne<LO,GO>(ColMap);
    A->fillComplete(ColMap1to1, NodeMap1to1);
  }

  static void extractDofMap(CrsMatrixGO & A,
			    RCP<const Map> dofMap)
  {
    Tpetra::global_size_t IGO =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    RCP<const Teuchos::Comm<int> > Comm = A.getRowMap()->getComm();
    LO numDof = A.getNodeNumEntries();
    std::vector<GO> globalIDs(numDof);
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const GO> Values;
    numDof = 0;
    for (size_t i=0; i<A.getNodeNumRows(); i++) {
      A.getLocalRowView(i, Indices, Values);
      for (LO j=0; j<Indices.size(); j++) {
	globalIDs[numDof++] = Values[j];
      }
    }
    dofMap = 
      rcp( new Map(IGO, Teuchos::ArrayView<GO>(globalIDs), 0, Comm) );
  }

  static void determineUniqueIndices(const std::vector<GO> & indices, 
				     std::vector<GO> & uniqueIndices)
  {
    uniqueIndices = indices;
    std::sort(uniqueIndices.begin(), uniqueIndices.end());
    auto iter = std::unique(uniqueIndices.begin(), uniqueIndices.end());
    uniqueIndices.erase(iter, uniqueIndices.end());
  }

  static void determineUniqueIndices(std::vector<GO> & indices)
  {
    std::sort(indices.begin(), indices.end());
    auto iter = std::unique(indices.begin(), indices.end());
    indices.erase(iter, indices.end());
  }

};

} // namespace bddc

#endif // BDDC_DOFMANAGER_H
  
