
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

#ifndef COARSESPACEBDDC_H
#define COARSESPACEBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include "shylu_enumsBDDC.h"

#include "shylu_SolverFactoryBDDC.h"
#include "shylu_PartitionOfUnityBDDC.h"
#include "shylu_SubdomainBDDC.h"
#include "shylu_DofManager.h"
#include "shylu_UtilBDDC.h"

// SRSR : Avoid depending on stk. See below for explanation.
//#include <stk_util/environment/memory_util.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {

template <class SX, class SM, class LO, class GO>
  class CoarseSpaceBDDC
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::Map<>::node_type  Node;
  typedef Tpetra::Map<LO,GO,Node>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO,Node>                            CrsGraph;
  typedef Tpetra::CrsMatrix<SX,LO,GO,Node>                        CrsMatrix;
  typedef Tpetra::CrsMatrix<GO,LO,GO,Node>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO,Node>                              Export;
  typedef Tpetra::Import<LO,GO,Node>                              Import;
  typedef Tpetra::Vector<SX,LO,GO,Node>                           Vector;
  typedef Tpetra::MultiVector<SX,LO,GO,Node>                      MV;

  CoarseSpaceBDDC()
  {
  }
  CoarseSpaceBDDC
    (std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & Subdomain,
     RCP< PartitionOfUnity<SX,SM,LO,GO> > & Partition,
     RCP<Export> & exporterB,
     RCP<Import> & importerB,
     std::vector< std::vector<LO> > & equivBoundaryDofs,
     std::vector< std::vector<LO> > & subBoundaryDofs,
     RCP<Teuchos::ParameterList> & Parameters) :
  m_Subdomain(Subdomain),
    m_Partition(Partition),
    m_exporterB(exporterB),
    m_importerB(importerB),
    m_dofMapB(exporterB->getSourceMap()),
    m_dofMapB1to1(importerB->getSourceMap()),
    m_Comm(exporterB->getSourceMap()->getComm()),
    m_equivBoundaryDofs(equivBoundaryDofs),
    m_subBoundaryDofs(subBoundaryDofs),
    m_Parameters(Parameters),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()),
    m_spatialDim(Parameters->get("Spatial Dimension", 3)),
    m_problemType(Parameters->get("Problem Type BDDC", SCALARPDE)),
    m_numSub(Subdomain.size()),
    m_rootProc(-1),
    m_numZeroEigenValues(-1),
    m_numCoarseNodes(0),
    m_coarseSolver(0)
  {
    determineGlobalIDs();
    assembleCoarseMatrices();
    coarsen();
    factorCoarseMatrix();
  }

  ~CoarseSpaceBDDC()
  {
    delete m_coarseSolver;
  }

  LO getCoarseSpaceDimension()
  {
    return m_AcRoot->getGlobalNumRows();
  }

  void applyCoarseCorrection(RCP<Vector> & rB,
                             RCP<Vector> & solB)
  {
    m_Phi->apply(*rB, *m_vecCoarse1to1, Teuchos::CONJ_TRANS);
    m_rootVec->doImport(*m_vecCoarse1to1, *m_importerInterlevel, Tpetra::INSERT);
    if (m_Comm->getRank() == m_rootProc) {
      Teuchos::ArrayRCP<SX> values = m_rootVec->getDataNonConst();
      for (size_t i=0; i<m_rhsCoarse.size(); i++) m_rhsCoarse[i] = values[i];
      if (m_coarseSolver != 0) {
        m_coarseSolver->Solve(1, &m_rhsCoarse[0], &m_solCoarse[0]);
      }
      for (size_t i=0; i<m_rhsCoarse.size(); i++) values[i] = m_solCoarse[i];
    }
    m_vecCoarse1to1->doExport(*m_rootVec, *m_exporterInterlevel, Tpetra::INSERT);
    m_Phi->apply(*m_vecCoarse1to1, *solB, Teuchos::NO_TRANS);
  }

  LO getNumZeroEigenValues() const
  {
    return m_numZeroEigenValues;
  }

  SM checkInterpolationMatrix()
  {
    // checkInterpolationMatrix is a diagnostic tool intended for examples with
    // no essential boundary conditions. Further, the coarse space must be
    // of a simple form where the number of coarse dofs for each equivalence
    // class is either 0 or the spatial dimension for elasticity problems
    // or 1 for Laplacian problems with a single dof per node.
    LO numRowsCoarse = m_dofMapCoarse1to1->getNodeNumElements();
    LO numRowsFine = m_dofMapB1to1->getNodeNumElements();
    SM maxError(0);
    if (m_problemType == SCALARPDE) {
      Teuchos::ArrayRCP<SX> values = m_vecCoarse1to1->getDataNonConst();
      for (LO i=0; i<numRowsCoarse; i++) values[i] = 1;
      m_Phi->apply(*m_vecCoarse1to1, *m_vecBoundary1to1, Teuchos::NO_TRANS);
      SM tol(1e-10);
      values = m_vecBoundary1to1->getDataNonConst();
      for (LO i=0; i<numRowsFine; i++) {
        SM error = std::abs(values[i] -1);
        if (error > maxError) maxError = error;
      }
    }
    else if (m_problemType == ELASTICITY) {
      LO numNodeCoarse = numRowsCoarse/m_spatialDim;
      LO numNodeFine = numRowsFine/m_spatialDim;
      assert (numNodeCoarse*m_spatialDim == numRowsCoarse);
      assert (numNodeFine*m_spatialDim == numRowsFine);
      Teuchos::ArrayRCP<SX> coarseValues = m_vecCoarse1to1->getDataNonConst();
      for (LO j=0; j<m_spatialDim; j++) {
        for (LO i=0; i<numRowsCoarse; i++) coarseValues[i] = 0;
        for (LO i=0; i<numNodeCoarse; i++) coarseValues[j+m_spatialDim*i] = 1;
        m_Phi->apply(*m_vecCoarse1to1, *m_vecBoundary1to1, Teuchos::NO_TRANS);
        SM tol(1e-10);
        Teuchos::ArrayRCP<SX> fineValues = m_vecBoundary1to1->getDataNonConst();
        for (LO i=0; i<numNodeFine; i++) {
          LO row = j+m_spatialDim*i;
          SM error = std::abs(fineValues[row] - 1);
          if (error > maxError) maxError = error;
        }
      }
    }
    SM maxErrorAll;
    Teuchos::reduceAll<int, SM>(*m_Comm, Teuchos::REDUCE_MAX, 1,
                                &maxError, &maxErrorAll);
    return maxErrorAll;
  }

  static void determineNodalConnectivity
    (const std::vector< std::vector<LO> > & elemConn,
     const LO numNode,
     std::vector< std::vector<LO> > & nodalConn)
  {
    std::vector< std::vector<LO> > nodeElements(numNode);
    LO numElem = elemConn.size();
    for (LO i=0; i<numElem; i++) {
      for (size_t j=0; j<elemConn[i].size(); j++) {
        LO node = elemConn[i][j];
        nodeElements[node].push_back(i);
      }
    }
    nodalConn.resize(numNode);
    std::vector<LO> anode(numNode);
    std::vector<bool> nodeFlag(numNode, false);
    for (LO i=0; i<numNode; i++) {
      LO nanode = 0;
      for (size_t j=0; j<nodeElements[i].size(); j++) {
        LO elem = nodeElements[i][j];
        for (size_t k=0; k<elemConn[elem].size(); k++) {
          LO node = elemConn[elem][k];
          if (nodeFlag[node] == false) {
            nodalConn[i].push_back(node);
            anode[nanode] = node;
            nodeFlag[node] = true;
            nanode++;
          }
        }
      }
      for (LO j=0; j<nanode; j++) {
        nodeFlag[anode[j]] = false;
      }
    }
  }

private:
  std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & m_Subdomain;
  RCP< PartitionOfUnity<SX,SM,LO,GO> > & m_Partition;
  RCP<Export> m_exporterB;
  RCP<Import> m_importerB;
  RCP<const Map> m_dofMapB, m_dofMapB1to1;
  RCP<const Teuchos::Comm<int> > m_Comm;
  std::vector< std::vector<LO> > m_equivBoundaryDofs, m_subBoundaryDofs;
  RCP<Teuchos::ParameterList> & m_Parameters;
  Tpetra::global_size_t m_IGO;
  LO m_spatialDim;
  enum ProblemType m_problemType;
  RCP<const Map> m_dofMapCoarse, m_dofMapCoarse1to1, m_dofMapCoarseRoot;
  LO m_numSub, m_rootProc, m_numZeroEigenValues, m_numCoarseNodes;
  std::vector<LO> m_nodeBeginCoarse, m_localDofsCoarse;
  std::vector<SX> m_rhsCoarse, m_solCoarse;
  RCP<Import> m_importerInterlevel;
  RCP<Export> m_exporterCoarse, m_exporterInterlevel;
  RCP<CrsMatrix> m_Phi, m_Ac, m_AcRoot;
  RCP<Vector> m_vecBoundary1to1, m_vecCoarse1to1, m_rootVec;
  SolverBase<SX>* m_coarseSolver;

  void determineReorder(LO sub,
                        std::vector<LO> & reorder)
  {
    const std::vector<LO> & coarseDofEquiv =
      m_Subdomain[sub]->getCoarseDofEquiv();
    LO numEquiv = m_Subdomain[sub]->getNumEquiv();
    LO numDofs = coarseDofEquiv.size();
    for (LO i=0; i<numEquiv; i++) {
      for (LO j=0; j<numDofs; j++) {
        if (coarseDofEquiv[j] == i) reorder.push_back(j);
      }
    }
    assert (LO(reorder.size()) == numDofs);
  }

  void determineCoarseStiffnessMatrix
    (const std::vector< std::vector<LO> > & subEquivs,
     const std::vector< std::vector<LO> > & equivConn,
     std::vector< std::vector<SX> > & subAc)
  {
    // determine sparsity pattern
    LO numNode = m_numCoarseNodes;
    const std::vector<LO> & nodeBegin = m_nodeBeginCoarse;
    LO numDof = nodeBegin[numNode];
    std::vector< std::vector<LO> > columnsA(numDof);
    for (LO i=0; i<numNode; i++) {
      for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
        for (size_t k=0; k<equivConn[i].size(); k++) {
          LO node2 = equivConn[i][k];
          for (LO m=nodeBegin[node2]; m<nodeBegin[node2+1]; m++) {
            columnsA[j].push_back(m);
          }
        }
      }
    }
    // coarse element assembly
    std::vector< std::vector<SX> > valuesA(numDof);
    for (LO i=0; i<numDof; i++) {
      valuesA[i].resize(columnsA[i].size());
    }
    std::vector<SX> rowValues;
    std::vector<LO> colMap(numDof, -1), subDofs;
    LO numSub = subEquivs.size();
    for (LO i=0; i<numSub; i++) {
      std::vector<SX> & Ac = subAc[i];
      std::vector<LO> reorder;
      determineReorder(i, reorder);
      determineSubDofs(subEquivs[i], nodeBegin, subDofs);
      for (size_t j=0; j<subDofs.size(); j++) {
        LO dof = subDofs[j];
        for (size_t k=0; k<columnsA[dof].size(); k++) {
          colMap[columnsA[dof][k]] = k;
        }
        getSubdomainMatrixRow(Ac, j, reorder, rowValues);
        for (size_t k=0; k<subDofs.size(); k++) {
          LO col = colMap[subDofs[k]];
          assert (col != -1);
          valuesA[dof][col] += rowValues[k];
        }
        for (size_t k=0; k<columnsA[dof].size(); k++) {
          colMap[columnsA[dof][k]] = -1;
        }
      }
    }
    Teuchos::ArrayRCP<size_t> count(numDof);
    for (LO i=0; i<numDof; i++) count[i] = columnsA[i].size();
    CrsMatrix ALocal(m_dofMapCoarse, m_dofMapCoarse, count,
                     Tpetra::StaticProfile);
    for (LO i=0; i<numDof; i++) {
      ALocal.insertLocalValues(i, Teuchos::ArrayView<LO>(columnsA[i]),
                               Teuchos::ArrayView<SX>(valuesA[i]));
    }
    ALocal.fillComplete(m_dofMapCoarse1to1, m_dofMapCoarse1to1);
    m_Ac = rcp( new CrsMatrix(m_dofMapCoarse1to1, 0) );
    m_Ac->doExport(ALocal, *m_exporterCoarse, Tpetra::ADD);
    m_Ac->fillComplete(m_dofMapCoarse1to1, m_dofMapCoarse1to1);
  }

  void getPhiRow(std::vector<SX> & Phi,
                 LO row,
                 LO numRows,
                 std::vector<LO> & reorder,
                 std::vector<SX> & rowValues) const
  {
    LO numCols = reorder.size();
    rowValues.resize(numCols);
    for (LO i=0; i<numCols; i++) {
      LO colUse = reorder[i];
      rowValues[i] = Phi[row+colUse*numRows];
    }
  }

  void determineSubDofs(const std::vector<LO> & subEquivs,
                        const std::vector<LO> & nodeBegin,
                        std::vector<LO> & subDofs) const
  {
    LO numSubDofs = 0;
    for (size_t i=0; i<subEquivs.size(); i++) {
      LO equiv = subEquivs[i];
      numSubDofs += nodeBegin[equiv+1] - nodeBegin[equiv];
    }
    subDofs.resize(numSubDofs);
    numSubDofs = 0;
    for (size_t i=0; i<subEquivs.size(); i++) {
      LO equiv = subEquivs[i];
      for (LO j=nodeBegin[equiv]; j<nodeBegin[equiv+1]; j++) {
        subDofs[numSubDofs++] = j;
      }
    }
  }

  void scalePhi(LO sub,
                std::vector<SX> & Phi) const
  {
    LO numRows = m_Subdomain[sub]->getNumBoundaryDofs();
    LO numCols = m_Subdomain[sub]->getCoarseSpaceDimension();
    std::vector<SX> work(numRows);
    bool applyTranspose = false;
    for (LO i=0; i<numCols; i++) {
      SX* colPhi = &Phi[i*numRows];
      for (LO j=0; j<numRows; j++) work[j] = colPhi[j];
      m_Subdomain[sub]->applyWeights(&work[0], applyTranspose, colPhi);
    }
  }

  void determineInterpolationMatrix
    (const std::vector< std::vector<LO> > & subEquivs,
     const std::vector< std::vector<LO> > & equivConn,
     std::vector< std::vector<SX> > & subPhi)
  {
    LO numEquiv = m_equivBoundaryDofs.size();
    LO numDofB = m_dofMapB->getNodeNumElements();
    std::vector< std::vector<LO> > columnsPhi(numDofB);
    const std::vector<LO> & nodeBegin = m_nodeBeginCoarse;
    for (LO i=0; i<numEquiv; i++) {
      LO numAdjEquiv = equivConn[i].size();
      LO numEquivDof = m_equivBoundaryDofs[i].size();
      for (LO j=0; j<numAdjEquiv; j++) {
        LO equiv = equivConn[i][j];
        for (LO k=nodeBegin[equiv]; k<nodeBegin[equiv+1]; k++) {
          LO dofCoarse = k;
          for (LO m=0; m<numEquivDof; m++) {
            LO dofB = m_equivBoundaryDofs[i][m];
            columnsPhi[dofB].push_back(dofCoarse);
          }
        }
      }
    }
    std::vector< std::vector<SX> > valuesPhi(numDofB);
    for (LO i=0; i<numDofB; i++) {
      valuesPhi[i].resize(columnsPhi[i].size());
    }
    LO numSub = m_Subdomain.size();
    LO numCols = nodeBegin[numEquiv];
    std::vector<LO> colMap(numCols, -1), subDofs;
    std::vector<SX> rowValues;
    for (LO i=0; i<numSub; i++) {
      std::vector<SX> & Phi = subPhi[i];
      scalePhi(i, Phi);
      std::vector<LO> reorder;
      determineReorder(i, reorder);
      determineSubDofs(subEquivs[i], nodeBegin, subDofs);
      LO numBoundaryDofs = m_Subdomain[i]->getNumBoundaryDofs();
      for (size_t j=0; j<m_subBoundaryDofs[i].size(); j++) {
        LO dofB = m_subBoundaryDofs[i][j];
        for (size_t k=0; k<columnsPhi[dofB].size(); k++) {
          colMap[columnsPhi[dofB][k]] = k;
        }
        getPhiRow(Phi, j, numBoundaryDofs, reorder, rowValues);
        for (size_t k=0; k<subDofs.size(); k++) {
          LO col = colMap[subDofs[k]];
          assert (col != -1);
          valuesPhi[dofB][col] += rowValues[k];
        }
        for (size_t k=0; k<columnsPhi[dofB].size(); k++) {
          colMap[columnsPhi[dofB][k]] = -1;
        }
      }
    }
    Teuchos::ArrayRCP<size_t> count(numDofB);
    for (LO i=0; i<numDofB; i++) count[i] = columnsPhi[i].size();
    CrsMatrix PhiLocal(m_dofMapB, m_dofMapCoarse, count, Tpetra::StaticProfile);
    for (LO i=0; i<numDofB; i++) {
      PhiLocal.insertLocalValues(i, Teuchos::ArrayView<LO>(columnsPhi[i]),
                                 Teuchos::ArrayView<SX>(valuesPhi[i]));
    }
    PhiLocal.fillComplete(m_dofMapCoarse1to1, m_dofMapB1to1);
    m_Phi = rcp( new CrsMatrix(m_dofMapB1to1, 0) );
    m_Phi->doExport(PhiLocal, *m_exporterB, Tpetra::ADD);
    m_Phi->fillComplete(m_dofMapCoarse1to1, m_dofMapB1to1);
  }

  void assembleCoarseMatrices()
  {
    LO numSub = m_Subdomain.size();
    std::vector< std::vector<SX> > subPhi(numSub), subAc(numSub);
    bool restrictPhiToBoundary = true;
    for (LO i=0; i<numSub; i++) {
      m_Subdomain[i]->calculateCoarseMatrices(subPhi[i], subAc[i],
                                              restrictPhiToBoundary);
    }
    const std::vector< std::vector<LO> > & subEquivs =
      m_Partition->getSubdomainEquivClasses();
    std::vector< std::vector<LO> > equivConn;
    determineNodalConnectivity(subEquivs, m_numCoarseNodes, equivConn);
    determineCoarseStiffnessMatrix(subEquivs, equivConn, subAc);
    determineInterpolationMatrix(subEquivs, equivConn, subPhi);
  }

  void getMinMemoryProcs(int numCoarseProcs,
                         std::vector<int> & minMemoryProcs)
  {
    long procLoad = getProcessorLoad();
    int numProc = m_Comm->getSize();
    std::vector<long> procLoadAll(numProc);
    Teuchos::gatherAll<int, long> (*m_Comm, 1, &procLoad, numProc,
                                   &procLoadAll[0]);
    std::vector< std::pair<long, int> > procSize(numProc);
    for (int i=0; i<numProc; i++) {
      procSize[i] = std::make_pair(procLoadAll[i], i);
    }
    std::sort(procSize.begin(), procSize.end());
    minMemoryProcs.resize(numCoarseProcs);
    typename std::vector< std::pair<long, int> >::const_iterator iter;
    for (int i=0; i<numCoarseProcs; i++) {
      iter = procSize.begin() + i;
      minMemoryProcs[i] = iter->second;
    }
  }

  void coarsen()
  {
    int numCoarseProcs = 1; // two level method for now
    std::vector<int> minMemoryProcs;
    getMinMemoryProcs(numCoarseProcs, minMemoryProcs);
    assert (numCoarseProcs == 1);
    m_rootProc = minMemoryProcs[0];
    DofManager<LO,GO>::generateRootMap(m_dofMapCoarse1to1, m_rootProc,
                                       m_dofMapCoarseRoot);
    m_importerInterlevel = rcp(new Import(m_dofMapCoarse1to1,
                                          m_dofMapCoarseRoot) );
    m_exporterInterlevel = rcp(new Export(m_dofMapCoarseRoot,
                                          m_dofMapCoarse1to1) );
    m_rootVec = rcp( new Vector(m_dofMapCoarseRoot) );
    m_AcRoot = rcp( new CrsMatrix(m_dofMapCoarseRoot, 0) );
    m_AcRoot->doImport(*m_Ac, *m_importerInterlevel, Tpetra::INSERT);
    m_AcRoot->fillComplete(m_dofMapCoarseRoot, m_dofMapCoarseRoot);
    getNumZeroEigenValuesCoarseStiffnessMatrix(*m_AcRoot);
  }

  void factorCoarseMatrix()
  {
    if (m_Comm->getRank() == m_rootProc) {
      LO numRows = m_AcRoot->getNodeNumRows();
      if (numRows == 0) return;
      m_rhsCoarse.resize(numRows);
      m_solCoarse.resize(numRows);
      LO numTerms = m_AcRoot->getNodeNumEntries();
      std::vector<LO> rowBegin(numRows+1), columns(numTerms);
      std::vector<SX> values(numTerms);
      Teuchos::ArrayView<const LO> Indices;
      Teuchos::ArrayView<const SX> Values;
      numTerms = 0;
      for (LO i=0; i<numRows; i++) {
        m_AcRoot->getLocalRowView(i, Indices, Values);
        for (LO j=0; j<Indices.size(); j++) {
          columns[numTerms] = Indices[j];
          values[numTerms] = Values[j];
          numTerms++;
        }
        rowBegin[i+1] = numTerms;
      }
      bool printCoarseMatrix = m_Parameters->get("Print Coarse Matrix", false);
      if (printCoarseMatrix) {
        UtilBDDC<SX,SM>::printSparseMatrix
          (numRows, &rowBegin[0], &columns[0], &values[0], "Ac.dat");
      }
      SolverFactory<SX> Factory;
      m_coarseSolver = Factory.Generate(numRows,
                                        &rowBegin[0],
                                        &columns[0],
                                        &values[0],
                                        *m_Parameters);
      m_coarseSolver->Initialize();
    }
  }

  long getProcessorLoad()
  {
    //    size_t s_now, s_hwm;
    // SRSR : We just retrun zero as the available memory. All the ranks
    // will have equal weight in the sorting that follows this call and
    // rank 0 will get picked based on tie-breaking. Consider finding the
    // current memory without depending on stk.
    //
    //stk::get_memory_usage(s_now, s_hwm);
    //long currentMemory = static_cast<long>(s_now);
    return 0;
  }

  void extractDenseMatrix(const CrsMatrix & A,
                          std::vector<SX> & AFull) const
  {
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const SX> Values;
    LO numRows = A.getNodeNumRows();
    LO numCols = A.getNodeNumCols();
    AFull.resize(numRows*numCols);
    AFull.assign(numRows*numCols, 0);
    for (LO i=0; i<numRows; i++) {
      A.getLocalRowView(i, Indices, Values);
      for (LO j=0; j<Indices.size(); j++) {
        LO col = Indices[j];
        AFull[i+col*numRows] = Values[j];
      }
    }
  }

  void getNumZeroEigenValuesCoarseStiffnessMatrix(CrsMatrix & Ac)
  {
    // getNumZeroEigenValuesCoarseStiffnessMatrix is a diagnostic tool
    // intended for smaller examples with no essential boundary conditions.
    bool checkMatrices = m_Parameters->get("Check Coarse Matrices", false);
    if (checkMatrices == false) return;
    m_numZeroEigenValues = -1;
    if (m_Comm->getRank() == m_rootProc) {
      std::vector<SX> AcFull;
      extractDenseMatrix(Ac, AcFull);
      Teuchos::LAPACK<int, SX> LAPACK;
      char JOBZ('N'), UPLO('U');
      int N = Ac.getNodeNumRows();
      std::vector<SM> W(N), RWORK(std::max(1, 3*N));
      int LWORK = std::max(1, 3*N);
      int INFO(0);
      std::vector<SX> WORK(LWORK);
      LAPACK.HEEV(JOBZ, UPLO, N, &AcFull[0], N, &W[0], &WORK[0], LWORK,
                  &RWORK[0], &INFO);
      assert (INFO == 0);
      SM tol(1e-10);
      SM compareValue = tol*W[N-1];
      if (m_problemType == SCALARPDE) {
        if (N == 1) compareValue = tol;
      }
      else if (m_problemType == ELASTICITY) {
        if (m_spatialDim == 2) {
          if (N < 4) compareValue = tol;
        }
        else if (m_spatialDim == 3) {
          if (N < 7) compareValue = tol;
        }
      }
      m_numZeroEigenValues = 0;
      for (LO i=0; i<N; i++) {
        if (std::abs(W[i]) < compareValue) m_numZeroEigenValues++;
      }
    }
    Teuchos::broadcast<int, LO> (*m_Comm, m_rootProc, 1, &m_numZeroEigenValues);
  }

  void getSubdomainMatrixRow(std::vector<SX> & Ac,
                             LO row,
                             std::vector<LO> & reorder,
                             std::vector<SX> & rowValues)
  {
    LO numElemDofs = reorder.size();
    rowValues.resize(numElemDofs);
    LO rowUse = reorder[row];
    for (LO i=0; i<numElemDofs; i++) {
      LO colUse = reorder[i];
      rowValues[i] = Ac[rowUse+colUse*numElemDofs];
    }
  }

  void determineGlobalIDs()
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses =
      m_Partition->getSubdomainEquivClasses();
    LO numCoarseNodes = m_Partition->getNumEquivClasses();
    std::vector<LO> numDofsCoarseNode(numCoarseNodes);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      for (LO j=0; j<numEquiv; j++) {
        LO equiv = subdomainEquivClasses[i][j];
        numDofsCoarseNode[equiv] = 0;
      }
      const std::vector<LO> & coarseDofEquiv =
        m_Subdomain[i]->getCoarseDofEquiv();
      for (size_t j=0; j<coarseDofEquiv.size(); j++) {
        LO equiv = subdomainEquivClasses[i][coarseDofEquiv[j]];
        numDofsCoarseNode[equiv]++;
      }
    }
    m_nodeBeginCoarse.resize(numCoarseNodes+1, 0);
    LO numCoarseDofs(0);
    for (LO i=0; i<numCoarseNodes; i++) numCoarseDofs += numDofsCoarseNode[i];
    m_localDofsCoarse.resize(numCoarseDofs);
    numCoarseDofs = 0;
    for (LO i=0; i<numCoarseNodes; i++) {
      for (LO j=0; j<numDofsCoarseNode[i]; j++) {
        m_localDofsCoarse[numCoarseDofs] = j;
        numCoarseDofs++;
      }
      m_nodeBeginCoarse[i+1] = numCoarseDofs;
    }
    const std::vector<GO> & coarseNodeGIDs = m_Partition->getGlobalIDs();
    DofManager<LO,GO>::determineGlobalIDs
      (numCoarseNodes, &coarseNodeGIDs[0], &m_nodeBeginCoarse[0],
       &m_localDofsCoarse[0], m_Comm, m_dofMapCoarse, m_dofMapCoarse1to1);
    m_numCoarseNodes = numCoarseNodes;
    m_exporterCoarse = rcp( new Export(m_dofMapCoarse, m_dofMapCoarse1to1) );
    m_vecBoundary1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_vecCoarse1to1 = rcp( new Vector(m_dofMapCoarse1to1) );
  }

};

} // namespace bddc

#endif // COARSESPACEBDDC_H

