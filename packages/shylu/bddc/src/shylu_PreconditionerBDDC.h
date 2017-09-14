
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

#ifndef PRECONDITIONERBDDC_H
#define PRECONDITIONERBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

#include "shylu_PartitionOfUnityBDDC.h"
#include "shylu_SubdomainBDDC.h"
#include "shylu_WeightsBDDC.h"
#include "shylu_ConstraintsBDDC.h"
#include "shylu_CoarseSpaceBDDC.h"
#include "shylu_DofManager.h"
#include "shylu_UtilBDDC.h"
#include "shylu_ZoltanPartition.h"
#include "shylu_SolverFactoryBDDC.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class PreconditionerBDDC
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::Map<LO,GO,Node>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO,Node>                            CrsGraph;
  typedef Tpetra::CrsMatrix<SX,LO,GO,Node>                        CrsMatrix;
  typedef Tpetra::CrsMatrix<LO,LO,GO,Node>                        CrsMatrixLO;
  typedef Tpetra::CrsMatrix<GO,LO,GO,Node>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO,Node>                              Export;
  typedef Tpetra::Import<LO,GO,Node>                              Import;
  typedef Tpetra::Vector<SM,LO,GO,Node>                           VectorSM;
  typedef Tpetra::Vector<SX,LO,GO,Node>                           Vector;
  typedef Tpetra::Vector<LO,LO,GO,Node>                           VectorLO;
  typedef Tpetra::Vector<GO,LO,GO,Node>                           VectorGO;
  typedef Tpetra::MultiVector<SX,LO,GO,Node>                      MV;
  typedef Tpetra::MultiVector<double,LO,GO,Node>                  MVD;

  PreconditionerBDDC()
  {
  }

  PreconditionerBDDC
    (LO numNodes,
     LO* nodeBegin,
     LO* localDofs,
     const GO* nodeGlobalIDs,
     const SM* xCoord,
     const SM* yCoord,
     const SM* zCoord,
     const std::vector< std::vector<LO> > subNodes,
     LO** subRowBegin,
     LO** subColumns,
     SX** subValues,
     RCP<Teuchos::ParameterList> & Parameters,
     MPI_Comm mpiComm,
     int level=0,
     std::vector<GO>* nodeGlobalIDs1to1=nullptr) :
    m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_subNodes(subNodes),
    m_subRowBegin(subRowBegin),
    m_subColumns(subColumns),
    m_subValues(subValues),
    m_Parameters(Parameters),
    m_mpiComm(mpiComm),
    m_usingSimpleInterface(false),
    m_level(level)
  {
    initializeVariables();
    initializePreconditioner(nodeGlobalIDs1to1);
  }
   
  PreconditionerBDDC
    (LO numNodes,
     const GO* nodeGlobalIDs,
     const SM* xCoord,
     const SM* yCoord,
     const SM* zCoord,
     LO* subRowBegin,
     LO* subColumns,
     SX* subValues,
     RCP<Teuchos::ParameterList> & Parameters,
     MPI_Comm mpiComm,
     int level=0,
     std::vector<GO>* nodeGlobalIDs1to1=nullptr) :
    m_numNodes(numNodes),
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_Parameters(Parameters),
    m_mpiComm(mpiComm),
    m_usingSimpleInterface(true),
    m_level(level)
  {
    convertSingleSudomain(subRowBegin, subColumns, subValues);
    initializeVariables();
    initializePreconditioner(nodeGlobalIDs1to1);
  }
   
  ~PreconditionerBDDC()
  {
    if (m_usingSimpleInterface) {
      delete [] m_subRowBegin;
      delete [] m_subColumns;
      delete [] m_subValues;
      delete [] m_nodeBegin;
      delete [] m_localDofs;
    }
    delete m_directSolver;
    for (LO i=0; i<m_numSub; i++) {
      delete m_Subdomain[i];
    }
#ifdef _OPENMP
    for (int i=0; i<m_numThreads; i++) {
      delete [] m_work1T[i];
      delete [] m_work2T[i];
    }
#endif
    delete [] m_work1T;
    delete [] m_work2T;
  }

  const std::vector<double> & getTimings() const
  {
    return m_timings;
  }

  bool useDirectSolver()
  {
    return m_useDirectSolver;
  }

  RCP<const Map> getDofMapB1to1()
  {
    return m_dofMapB1to1;
  }

  RCP<const Map> getDofMap1to1()
  {
    return m_dofMap1to1;
  }

  RCP<const Teuchos::Comm<int> > getComm() const
  {
    return m_Comm;
  }

  LO getNumSub() const
  {
    return m_numSub;
  }

  int InitialStaticCondensation(SX* rightHandSide1to1,
				SX* initialSolution1to1)
  {
    if (m_interfacePreconditioner == false) return 0;
    double startTime = GetTime();
    Teuchos::ArrayRCP<SX> rhsB = m_rhsVecB->getDataNonConst();
    for (LO i=0; i<m_numDofsB; i++) rhsB[i] = 0;
#pragma omp parallel num_threads(m_numThreads)
    {
#pragma omp for
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & interiorDofs = m_subInteriorDofs[i];
	LO numInteriorDofs = interiorDofs.size();
	SX *subRhs(0), *subSol(0);
	int threadID = getThreadID();
	getArrays(threadID, subRhs, subSol);
	for (LO j=0; j<numInteriorDofs; j++) {
	  LO row = interiorDofs[j];
	  subRhs[j] = rightHandSide1to1[row];
	}
	m_Subdomain[i]->removeInteriorResiduals(&subRhs[0], &subSol[0]);
	// no potential write conflicts because interior dofs are unique
	// to each thread
	for (LO j=0; j<numInteriorDofs; j++) {
	  LO row = interiorDofs[j];
	  initialSolution1to1[row] += subSol[j];
	}
	const std::vector<LO> & boundaryDofs = m_subBoundaryDofs[i];
	LO numBoundaryDofs = boundaryDofs.size();
	// potential write conflicts because boundary dofs are common
	// to different threads
	//#pragma omp critical(initial_static_condensation)
	{
	  for (LO j=0; j<numBoundaryDofs; j++) {
	    LO row = boundaryDofs[j];
#pragma omp atomic
	    rhsB[row] += subRhs[j];
	  }
	}
      }
    }
    m_rhsVecB1to1->doExport(*m_rhsVecB, *m_exporterB, Tpetra::ADD);
    Teuchos::ArrayRCP<SX> rhsB1to1 = m_rhsVecB1to1->getDataNonConst();
    for (LO i=0; i<m_numDofsB1to1; i++) {
      LO row = m_boundaryToAll1to1[i];
      rhsB1to1[i] += rightHandSide1to1[row];
    }
    for (LO i=0; i<m_numDofsB1to1; i++) { 
      rightHandSide1to1[i] = rhsB1to1[i];
    }
    m_timings[TIME_INIT_STATIC_COND] += GetTime() - startTime;
    return 1;
  }

  void StaticExpansion(SX* deltaSolution,
		       SX* initialSolution)
  {
    if (m_interfacePreconditioner) {
      double startTime = GetTime();
      Teuchos::ArrayRCP<SX> solB1to1 = m_solVecB1to1->getDataNonConst();
      for (LO i=0; i<m_numDofsB1to1; i++) solB1to1[i] = deltaSolution[i];
      m_solVecB->doImport(*m_solVecB1to1, *m_exporterB, Tpetra::INSERT);
      Teuchos::ArrayRCP<SX> solB = m_solVecB->getDataNonConst();
#pragma omp parallel num_threads(m_numThreads)
      {
#pragma omp for
	for (LO i=0; i<m_numSub; i++) {
	  const std::vector<LO> & subBoundaryDofs = m_subBoundaryDofs[i];
	  const std::vector<LO> & interiorDofs = m_subInteriorDofs[i];
	  LO numDofsB = subBoundaryDofs.size();
	  LO numDofsI = interiorDofs.size();
	  SX *xB(0), *xI(0);
	  int threadID = getThreadID();
	  getArrays(threadID, xB, xI);
	  for (LO j=0; j<numDofsB; j++) {
	    LO row = subBoundaryDofs[j];
	    xB[j] = solB[row];
	  }
	  m_Subdomain[i]->staticExpansion(&xB[0], &xI[0]);
	  // no potential write conflicts because interior dofs are unique
	  // to each thread
	  for (LO j=0; j<numDofsI; j++) {
	    LO row = interiorDofs[j];
	    initialSolution[row] += xI[j];
	  }
	}
      }
      for (LO i=0; i<m_numDofsB1to1; i++) {
	LO row = m_boundaryToAll1to1[i];
	initialSolution[row] += deltaSolution[i];
      }
      m_timings[TIME_STATIC_EXPANSION] += GetTime() - startTime;
    }
    else {
      LO numDofs1to1 = m_dofMap1to1->getNodeNumElements();
      for (LO i=0; i<numDofs1to1; i++) initialSolution[i] += deltaSolution[i];
    }
  }

  LO NumMyRows()
  {
    return m_dofMap1to1->getNodeNumElements();
  }

  LO NumMyRowsKrylov()
  {
    if (m_interfacePreconditioner == true) {
      return m_numDofsB1to1;
    }
    else {
      return m_dofMap1to1->getNodeNumElements();
    }
  }

  void GatherLocalDof(std::vector<LO> & localDof)
  {
    VectorLO localDofsVec(m_dofMap);
    Teuchos::ArrayRCP<LO> localDofs = localDofsVec.getDataNonConst();
    for (LO i=0; i<m_numDofs; i++) {
      localDofs[i] = m_localDofs[i];
    }
    VectorLO localDofsVec1to1(m_dofMap1to1);
    localDofsVec1to1.doExport(localDofsVec, *m_exporterAll, Tpetra::INSERT);
    Teuchos::ArrayRCP<LO> localDofs1to1 = localDofsVec1to1.getDataNonConst();
    LO numDofs1to1 = m_dofMap1to1->getNodeNumElements();
    localDof.resize(numDofs1to1);
    for (LO i=0; i<numDofs1to1; i++) {
      localDof[i] = localDofs1to1[i];
    }
  }

  void GatherScaleFactors(std::vector<SM> & scaleFactors)
  {
    scaleFactors.resize(m_dofMap1to1->getNodeNumElements(), 1);
  }

  void Apply(SX* r, 
	     SX* Pr,
	     SX* APr)
  {
    Apply(r, Pr);
    double startTime = GetTime();
    ApplyOperator(Pr, APr);
    m_timings[TIME_APPLY_OPER] += GetTime() - startTime;
    
  }

  void Apply(SX* r, 
	     SX* Pr)
  {
    if (m_useDirectSolver == true) {
      m_directSolver->Solve(1, r, Pr);
      return;
    }
    SX* PrInterface = Pr;
    SX* rUse = r;
    // apply initial static condensation correction if not at finest level
    if (m_level > 0) {
      LO numRows = NumMyRows();
      m_residualVector.resize(numRows);
      m_deltaSol.assign(numRows, 0);
      for (LO i=0; i<numRows; i++) {
	m_residualVector[i] = r[i];
	Pr[i] = 0;
      }
      InitialStaticCondensation(m_residualVector.data(), Pr);
      PrInterface = m_deltaSol.data();
      rUse = m_residualVector.data();
    }
    // assign Vectors, importer, exporter based on preconditioner type
    // and whether or not fine/coarse communicators are disjoint
    RCP<Vector> rhsVec1to1, solVec1to1, rhsVec, solVec;
    RCP<Export> exporter;
    std::vector< std::vector<LO> > & subDofs = m_subBoundaryDofs;
    setTpetraObjects(rhsVec1to1, solVec1to1, rhsVec, solVec, 
		     exporter, subDofs);
    // restrict residual to next coarser level
    Teuchos::ArrayRCP<SX> rhsData = rhsVec1to1->getDataNonConst();
    LO numDof = rhsVec1to1->getMap()->getNodeNumElements();
    for (LO i=0; i<numDof; i++) rhsData[i] = rUse[i];
    applyPhiTranspose(rhsVec1to1, m_coarseVecNextLevel, 
		      m_interfacePreconditioner);
    LO numCoarseDof = m_coarseVecNextLevel->getMap()->getNodeNumElements();
    Teuchos::ArrayRCP<const SX> rhsCoarse = m_coarseVecNextLevel->getData();
    for (LO i=0; i<numCoarseDof; i++) m_coarseRhsWork[i] = rhsCoarse[i];
    // calculate coarse and local (subdomain) corrections
    if (m_disjointCommunicators) {
      // asynchronous work across different levels
      if (m_coarseProc == true) {
	applyCoarseCorrection(m_coarseRhsWork.data(), m_coarseSolWork.data());
      }
      else {	
	changeMapsForTpetraObjects(rhsVec1to1, solVec1to1, rhsVec, solVec,
				   exporter);
	applyLocalCorrections(rhsVec1to1, rhsVec, solVec1to1, solVec,
			      exporter, subDofs);
	resetMapsForTpetraObjects(rhsVec1to1, solVec1to1, rhsVec, solVec,
				  exporter);
      }
    }
    else {
      if (m_coarseProc == true) {
	applyCoarseCorrection(m_coarseRhsWork.data(), m_coarseSolWork.data());
      }
      applyLocalCorrections(rhsVec1to1, rhsVec, solVec1to1, solVec, 
			    exporter, subDofs); 
    }
    // sum of coarse and local (subdomain) corrections
    Teuchos::ArrayRCP<SX> solCoarse = m_coarseVecNextLevel->getDataNonConst();
    for (LO i=0; i<numCoarseDof; i++) solCoarse[i] = m_coarseSolWork[i];
    SX beta(1);
    Teuchos::ArrayRCP<SX> solData1to1 = solVec1to1->getDataNonConst();
    applyPhi(m_coarseVecNextLevel, solVec1to1, beta, m_interfacePreconditioner);
    for (LO i=0; i<m_numDofsB1to1; i++) PrInterface[i] = solData1to1[i];
    // apply static expansion if not at finest level
    if (m_level > 0) {
      StaticExpansion(PrInterface, Pr);
    }
  }

  void applyCoarseCorrection(SX* coarseRhs, 
			     SX* coarseSol)
  {
    if (m_coarsePreconditioner->useDirectSolver() == true) {
      if (m_reorderCoarseDofs.size() == 0) {
	determineReorderCoarseDofs();
      }
      LO numRows = m_reorderCoarseDofs.size();
      m_coarseRhsWork2.resize(numRows);
      for (LO i=0; i<numRows; i++) {
	LO row = m_reorderCoarseDofs[i];
	m_coarseRhsWork2[i] = coarseRhs[row];
      }
      m_coarsePreconditioner->Apply(m_coarseRhsWork2.data(), coarseRhs);
      for (LO i=0; i<numRows; i++) {
	LO row = m_reorderCoarseDofs[i];
	coarseSol[row] = coarseRhs[i];
      }
    }
    else {
      m_coarsePreconditioner->Apply(coarseRhs, coarseSol);
    }
  }

  void applyCoarseCorrection(std::vector<SX> & m_coarseRhsWork, 
			     std::vector<SX> & m_coarseSolWork)
  {
    double startTime = GetTime();
    m_coarsePreconditioner->Apply(m_coarseRhsWork.data(),
				  m_coarseSolWork.data());
    m_timings[TIME_COARSE_CORR] += GetTime() - startTime;
  }

  void applyLocalCorrections(RCP<Vector> & rhsVec1to1,
			     RCP<Vector> & rhsVec,
			     RCP<Vector> & solVec1to1, 
			     RCP<Vector> & solVec,
			     RCP<Export> & exporter,
			     std::vector< std::vector<LO> > & subDofs)
  {
    double startTime = GetTime();
    rhsVec->doImport(*rhsVec1to1, *exporter, Tpetra::INSERT);
    Teuchos::ArrayRCP<SX> rhs = rhsVec->getDataNonConst();
    Teuchos::ArrayRCP<SX> sol = solVec->getDataNonConst();
    solVec->putScalar(0);
#pragma omp parallel num_threads(m_numThreads)
    {
#pragma omp for
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & subdofs = subDofs[i];
	LO numDofSub = subdofs.size();
	SX *subRhs(0), *subSol(0);
	int threadID = getThreadID();
	getArrays(threadID, subRhs, subSol);
	for (LO j=0; j<numDofSub; j++) {
	  subRhs[j] = rhs[subdofs[j]];
	}
	m_Subdomain[i]->applyNeumannCorrection(subRhs, subSol);
	//#pragma omp critical(apply)
	{
	  for (LO j=0; j<numDofSub; j++) {
	    int row = subdofs[j];
#pragma omp atomic
	    sol[row] += subSol[j];
	  }
	}
      }
    }
    solVec1to1->putScalar(0);
    solVec1to1->doExport(*solVec, *exporter, Tpetra::ADD);
    m_timings[TIME_SUB_CORR] += GetTime() - startTime;
  }

  void ApplyOperator(SX* x, 
		     SX* Ax)
  {
    if (m_interfacePreconditioner) {
      RCP<Vector> & xVec1to1 = m_rhsVecB1to1;
      Teuchos::ArrayRCP<SX> xValues1to1 = xVec1to1->getDataNonConst();
      for (LO i=0; i<m_numDofsB1to1; i++) xValues1to1[i] = x[i];
      RCP<Vector> & xVec = m_rhsVecB;
      xVec->doImport(*xVec1to1, *m_exporterB, Tpetra::INSERT);
      Teuchos::ArrayRCP<SX> xValues = xVec->getDataNonConst();
      RCP<Vector> & AxVec = m_solVecB;
      Teuchos::ArrayRCP<SX> AxValues = AxVec->getDataNonConst();
      for (LO i=0; i<m_numDofsB; i++) AxValues[i] = 0;
#pragma omp parallel num_threads(m_numThreads)
      {
#if 0
        if (omp_get_thread_num() == 0)
          std::cout << "PreconditionerBDDC::ApplyOperator nthreads "
                    << omp_get_num_threads() << "\n";
#endif
#pragma omp for
	for (LO i=0; i<m_numSub; i++) {
	  const std::vector<LO> & subBoundaryDofs = m_subBoundaryDofs[i];
	  LO numDofB = subBoundaryDofs.size();
	  SX *xSub(0), *AxSub(0);
	  int threadID = getThreadID();
	  getArrays(threadID, xSub, AxSub);
	  for (LO j=0; j<numDofB; j++) {
	    xSub[j] = xValues[subBoundaryDofs[j]];
	  }
	  m_Subdomain[i]->applyBoundaryOperator(&xSub[0], &AxSub[0]);
	  //#pragma omp critical(apply_operator)
	  {
	    for (LO j=0; j<numDofB; j++) {
	      int row = subBoundaryDofs[j]; 
#pragma omp atomic
	      AxValues[row] += AxSub[j];
	    }
	  }
	}
      }
      RCP<Vector> & AxVec1to1 = m_solVecB1to1;
      AxVec1to1->doExport(*AxVec, *m_exporterB, Tpetra::ADD);
      Teuchos::ArrayRCP<SX> AxValues1to1 = AxVec1to1->getDataNonConst();
      for (LO i=0; i<m_numDofsB1to1; i++) Ax[i] = AxValues1to1[i];
    }
    else {
      ApplyFullOperator(x, Ax);
    }
  }

  void ApplyFullOperator(SX* x, 
			 SX* Ax)
  {
    Teuchos::ArrayRCP<SX> xValues1to1 = m_xVecAll1to1->getDataNonConst();
    LO numDofs1to1 = m_dofMap1to1->getNodeNumElements();
    for (LO i=0; i<numDofs1to1; i++) xValues1to1[i] = x[i];
    m_xVecAll->doImport(*m_xVecAll1to1, *m_exporterAll, Tpetra::INSERT);
    Teuchos::ArrayRCP<SX> xValues = m_xVecAll->getDataNonConst();
    Teuchos::ArrayRCP<SX> AxValues = m_AxVecAll->getDataNonConst();
    for (LO i=0; i<m_numDofs; i++) AxValues[i] = 0;
#pragma omp parallel num_threads(m_numThreads)
    {
#if 0
      if (omp_get_thread_num() == 0)
        std::cout << "PreconditionerBDDC::ApplyFullOperator nthreads "
                  << omp_get_num_threads()
                  << " limit " << omp_get_thread_limit()
                  << " nested lvl " << omp_get_level()
                  << " ancestor tid " << omp_get_ancestor_thread_num(omp_get_level() - 1)
                  << " active lvl " << omp_get_active_level()
                  << " team sz " << omp_get_team_size(omp_get_level() - 1) << "\n";
#endif
#pragma omp for
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & subDofs = m_subDofs[i];
	LO numDofs = subDofs.size();
	SX *xSub(0), *AxSub(0);
	int threadID = getThreadID();
	getArrays(threadID, xSub, AxSub);
	for (LO j=0; j<numDofs; j++) {
	  xSub[j] = xValues[subDofs[j]];
	}
	m_Subdomain[i]->applyFullOperator(&xSub[0], &AxSub[0]);
	//#pragma omp critical(apply_full_operator)
	{
	  for (LO j=0; j<numDofs; j++) {
	    int row = subDofs[j];
#pragma omp atomic
	    AxValues[row] += AxSub[j];
	  }
	}
      }
    }
    m_AxVecAll1to1->doExport(*m_AxVecAll, *m_exporterAll, Tpetra::ADD);
    Teuchos::ArrayRCP<SX> AxValuesAll1to1 = m_AxVecAll1to1->getDataNonConst();
    for (LO i=0; i<numDofs1to1; i++) Ax[i] = AxValuesAll1to1[i];
  }

  void ReduceSum(SX* inputVec,
		 LO numTerms,
		 SX* summedVec)
  {
    Teuchos::reduceAll<int, SX> (*m_Comm, Teuchos::REDUCE_SUM, numTerms, 
				 inputVec, summedVec);
  }

  SM Norm2(SX* x, 
	   LO numTerm) 
  {
    SX dotprod = DotProd(x, x, numTerm);
    return sqrt(UtilBDDC<SX,SM>::real(dotprod));
  }

  SX DotProd(SX* x, 
	     SX* y, 
	     LO numTerm)
  {
    SX localDotProd(0), dotProd;
    for (LO i=0; i<numTerm; i++) {
      localDotProd += UtilBDDC<SX,SM>::conj(x[i])*y[i];
    }
    Teuchos::reduceAll<int, SX> (*m_Comm, Teuchos::REDUCE_SUM, 1, 
				 &localDotProd, &dotProd);
    return dotProd;
  }

  int MyPID() 
  {
    return m_myPID;
  }

  double GetTime()
  {
    return MPI_Wtime();
  }

private:
  LO m_numNodes;
  LO *m_nodeBegin, *m_localDofs;
  const GO* m_nodeGlobalIDs;
  const SM *m_xCoord, *m_yCoord, *m_zCoord;
  std::vector< std::vector<LO> > m_subNodes;
  LO **m_subRowBegin, **m_subColumns;
  SX  **m_subValues;
  RCP<Teuchos::ParameterList> & m_Parameters;
  MPI_Comm m_mpiComm, m_mpiCommSplit;
  LO m_spatialDim;
  enum bddc::ProblemType m_problemType; 
  RCP<const Teuchos::Comm<int> > m_Comm, m_CommFine;
  int m_myPID, m_numThreads;
  Tpetra::global_size_t m_IGO;
  bool m_interfacePreconditioner, m_usingSimpleInterface;
  LO m_numDofs, m_numSub, m_numDofsB, m_numDofsB1to1;
  std::vector< std::vector<LO> > m_subNodeBegin, m_subLocalDofs;
  RCP< bddc::PartitionOfUnity<SX,SM,LO,GO> > m_Partition;
  std::vector< bddc::SubdomainBDDC<SX,SM,LO,GO>* > m_Subdomain;
  RCP< bddc::WeightsBDDC<SX,SM,LO,GO> > m_Weights;
  RCP< bddc::ConstraintsBDDC<SX,SM,LO,GO> > m_Constraints;
  std::vector< std::vector<LO> > m_boundaryDofsLocal, m_subBoundaryDofs,
    m_subDofs, m_subInteriorDofs;
  std::vector<LO> m_boundaryDofs, m_globalToBoundaryMap, m_boundaryToAll1to1;
  std::vector<SM> m_diagBoundary;
  std::vector< std::vector<LO> > m_equivBoundaryDofs;
  RCP<const Map> m_dofMap, m_dofMap1to1, m_dofMapB, m_dofMapB1to1,
    m_coarseMap, m_coarseMapNextLevel;
  RCP<Export> m_exporterB, m_exporterAll, m_exporterCoarse;
  RCP<Vector> m_rhsVecB, m_solVecB, m_rhsVecB1to1, m_solVecB1to1,
    m_xVecAll, m_AxVecAll,
    m_xVecAll1to1, m_AxVecAll1to1, m_rhsVecAll, m_solVecAll,
    m_coarseVecNextLevel, m_rhsVecAll1to1, m_solVecAll1to1, m_coarseVec,
    m_fineVec, m_fineVec1to1, m_fineVecB, m_fineVecB1to1;
  std::vector<SX> m_work1, m_work2;
  SX **m_work1T, **m_work2T;
  std::vector<double> m_timings;
  // coarsening data
  LO m_numNodesCoarse;
  std::vector<LO> m_nodeBeginCoarse, m_localDofsCoarse, m_reorderCoarseDofs;
  std::vector<GO> m_nodeGlobalIDsCoarse, m_nodeGlobalIDsCoarse1to1;
  std::vector<SM> m_xCoordCoarse, m_yCoordCoarse, m_zCoordCoarse;
  std::vector<SX> m_coarseRhsWork, m_coarseSolWork, m_fine1to1Work,
    m_coarseRhsWork2, m_residualVector, m_deltaSol;
  std::vector< std::vector<LO> > m_subNodesCoarse, m_subRowBeginCoarse, 
    m_subColumnsCoarse, m_subDofsCoarse;
  std::vector< std::vector<SX> > m_subValuesCoarse;
  RCP< PreconditionerBDDC<SX,SM,LO,GO> > m_coarsePreconditioner;
  SolverBase<SX>* m_directSolver;
  int m_level, m_myCoarseMpiRank;
  bool m_useDirectSolver, m_fineProc, m_coarseProc, m_disjointCommunicators;
  // data for disjoint communicators
  RCP<const Map> m_dofMapB_dis, m_dofMapB1to1_dis, m_dofMap_dis,
    m_dofMap1to1_dis;
  RCP<Export> m_exporterB_dis, m_exporterAll_dis;
  RCP<Vector> m_rhsVecB_dis, m_solVecB_dis, m_rhsVecB1to1_dis,
    m_solVecB1to1_dis, m_rhsVecAll_dis, m_solVecAll_dis,
    m_rhsVecAll1to1_dis, m_solVecAll1to1_dis;

  void convertSingleSudomain(LO* subRowBegin, 
			     LO* subColumns, 
			     SX* subValues)
  {
    m_subNodes.resize(1);
    m_subNodes[0].resize(m_numNodes);
    for (LO i=0; i<m_numNodes; i++) m_subNodes[0][i] = i;
    m_subRowBegin = new LO*[1]; m_subRowBegin[0] = subRowBegin;
    m_subColumns = new LO*[1]; m_subColumns[0] = subColumns;
    m_subValues = new SX*[1]; m_subValues[0] = subValues;
    LO spatialDim = m_Parameters->get("Spatial Dimension", 3);
    enum bddc::ProblemType problemType = 
      m_Parameters->get("Problem Type BDDC", bddc::SCALARPDE);
    int numDofPerNode = 0;
    if (problemType == bddc::SCALARPDE) {
      numDofPerNode = 1;
    }
    else if (problemType == bddc::ELASTICITY) {
      numDofPerNode = spatialDim;
    }
    assert (numDofPerNode != 0);
    LO numDof = m_numNodes*numDofPerNode;
    m_nodeBegin = new LO[m_numNodes+1]; m_nodeBegin[0] = 0;
    m_localDofs = new LO[numDof];
    numDof = 0;
    for (LO i=0; i<m_numNodes; i++) {
      for (int j=0; j<numDofPerNode; j++) {
	m_localDofs[numDof++] = j;
      }
      m_nodeBegin[i+1] = numDof;
    }
  }

  void initializeVariables()
  {
    m_spatialDim = m_Parameters->get("Spatial Dimension", 3);
    m_problemType = m_Parameters->get("Problem Type BDDC", bddc::SCALARPDE);
    m_Comm = rcp( new Teuchos::MpiComm<int>(m_mpiComm) );
    m_myPID = m_Comm->getRank();
    m_numThreads = m_Parameters->get<int>("numThreadsOuter", 1);
    m_IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    m_interfacePreconditioner =
      m_Parameters->get("Interface Preconditioner", true);
    m_numDofs = m_nodeBegin[m_numNodes];
    m_numSub = 0;
    m_numDofsB = 0;
    m_numDofsB1to1 = 0;
    m_work1T = nullptr;
    m_work2T = nullptr;
    m_directSolver = 0;
    m_myCoarseMpiRank = -1;
    m_useDirectSolver = false;
    m_fineProc = false;
    m_coarseProc = false;
    m_disjointCommunicators = false;
  }

  void initializePreconditioner(std::vector<GO>* nodeGlobalIDs1to1)
  {
    determineSubNodeBeginAndSubLocalDofs();
    m_useDirectSolver = checkUseDirectSolver();
    if (m_useDirectSolver == true) {
      initializeDirectSolverAndVectors();
      return;
    }
    double startTime = GetTime();
    m_Partition =
      rcp( new bddc::PartitionOfUnity<SX,SM,LO,GO>
	   (m_numNodes, m_nodeGlobalIDs, m_subNodeBegin, m_subNodes, 
	    m_subRowBegin, m_subColumns, m_subValues,
	    m_spatialDim, m_Parameters, m_Comm) );
    bddc::DofManager<LO,GO>::
      determineGlobalIDs(m_numNodes, m_nodeGlobalIDs, m_nodeBegin, m_localDofs,
			 m_Comm, m_dofMap, m_dofMap1to1, nodeGlobalIDs1to1);
    determineBoundaryDofs();
    determineBoundaryMaps();
    initializeSubdomains();
    determineInteriorMaps();
    determineDiagBoundary();
    determineBaseConstraints();
    determineWeights();
    determineAuxiliaryConstraints();
    determineCoarseSpace();
    initializeVectors();
    initializeDisjointData();
    reserveMemory();
    m_timings[TIME_INITIALIZATION] += GetTime() - startTime;
  }

  void setTpetraObjects(RCP<Vector> & rhsVec1to1, 
			RCP<Vector> & solVec1to1, 
			RCP<Vector> & rhsVec, 
			RCP<Vector> & solVec, 
			RCP<Export> & exporter,
			std::vector< std::vector<LO> > & subDofs)
  {
    subDofs = m_subBoundaryDofs;
    if (m_interfacePreconditioner == true) {
      rhsVec1to1 = m_rhsVecB1to1;
      solVec1to1 = m_solVecB1to1;
      rhsVec = m_rhsVecB;
      solVec = m_solVecB;
      exporter = m_exporterB;
    }
    else {
      subDofs = m_subDofs;
      rhsVec1to1 = m_rhsVecAll1to1;
      solVec1to1 = m_solVecAll1to1;
      rhsVec = m_rhsVecAll;
      solVec = m_solVecAll;
      exporter = m_exporterAll;
    }
  }

  void resetMapsForTpetraObjects(RCP<Vector> & rhsVec1to1, 
				 RCP<Vector> & solVec1to1, 
				 RCP<Vector> & rhsVec, 
				 RCP<Vector> & solVec,
				 RCP<Export> & exporter)
  {
    if (m_interfacePreconditioner == true) {
      rhsVec1to1->replaceMap(m_dofMapB1to1);
      solVec1to1->replaceMap(m_dofMapB1to1);
      rhsVec->replaceMap(m_dofMapB);
      solVec->replaceMap(m_dofMapB);
      exporter = m_exporterB;
    }
    else {
      rhsVec1to1->replaceMap(m_dofMap1to1);
      solVec1to1->replaceMap(m_dofMap1to1);
      rhsVec->replaceMap(m_dofMap);
      solVec->replaceMap(m_dofMap);
      exporter = m_exporterAll;
    }
  }

  void changeMapsForTpetraObjects(RCP<Vector> & rhsVec1to1, 
				  RCP<Vector> & solVec1to1, 
				  RCP<Vector> & rhsVec, 
				  RCP<Vector> & solVec,
				  RCP<Export> & exporter)
  {
    if (m_interfacePreconditioner == true) {
      rhsVec1to1->replaceMap(m_dofMapB1to1_dis);
      solVec1to1->replaceMap(m_dofMapB1to1_dis);
      rhsVec->replaceMap(m_dofMapB_dis);
      solVec->replaceMap(m_dofMapB_dis);
      exporter = m_exporterB_dis;
    }
    else {
      rhsVec1to1->replaceMap(m_dofMap1to1_dis);
      solVec1to1->replaceMap(m_dofMap1to1_dis);
      rhsVec->replaceMap(m_dofMap_dis);
      solVec->replaceMap(m_dofMap_dis);
      exporter = m_exporterAll_dis;
    }
  }

  void determineReorderCoarseDofs()
  {
    LO numRows = m_localDofsCoarse.size();
    m_reorderCoarseDofs.resize(numRows);
    assert (m_subNodesCoarse.size() == 1);
    const std::vector<LO> & subNodesCoarse = m_subNodesCoarse[0];
    assert (LO(subNodesCoarse.size()) == m_numNodesCoarse);
    numRows = 0;
    for (LO i=0; i<m_numNodesCoarse; i++) {
      LO node = subNodesCoarse[i];
      for (LO j=m_nodeBeginCoarse[node]; j<m_nodeBeginCoarse[node+1]; j++) {
	m_reorderCoarseDofs[numRows++] = j;
      }
    }
  }

  int getThreadID()
  {
    int threadID(0);
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    return threadID;
  }

  void getArrays(int threadID, 
		 SX* & array1, 
		 SX* & array2)
  {
#ifdef _OPENMP
    array1 = m_work1T[threadID];
    array2 = m_work2T[threadID];
#else
    array1 = &m_work1[0];
    array2 = &m_work2[0];
#endif
  }

  void reserveMemory()
  {
    m_timings.resize(TIME_LENGTH);
    LO maxNumDof(0);
    for (LO i=0; i<m_numSub; i++) {
      LO numDof = m_Subdomain[i]->getNumDofs();
      if (numDof > maxNumDof) maxNumDof = numDof;
    }
#ifdef _OPENMP
    int nthreads_obtained(0);
#pragma omp parallel num_threads(m_numThreads)
    {
      nthreads_obtained = m_numThreads = omp_get_num_threads();
    }
    assert(nthreads_obtained == m_numThreads);
    m_work1T = new SX*[m_numThreads];
    m_work2T = new SX*[m_numThreads];
    for (int i=0; i<m_numThreads; i++) {
      m_work1T[i] = new SX[maxNumDof];
      m_work2T[i] = new SX[maxNumDof];
    }
#else
    m_work1.resize(maxNumDof);
    m_work2.resize(maxNumDof);
#endif
  }

  void determineNewMap(RCP<const Map> dofMap, 
		       RCP<const Map> & dofMap_dis)
  {
    LO numRows = dofMap->getNodeNumElements();
    std::vector<GO> globalIDs(numRows);
    for (LO i=0; i<numRows; i++) {
      globalIDs[i] = dofMap->getGlobalElement(i);
    }
    dofMap_dis = 
      rcp( new Map(m_IGO, Teuchos::ArrayView<GO>(globalIDs), 0, m_CommFine) );
  }

  void determineSubNodeBeginAndSubLocalDofs()
  {
    LO numSub = m_subNodes.size();
    m_subNodeBegin.resize(numSub);
    m_subLocalDofs.resize(numSub);
    for (LO i=0; i<numSub; i++) {
      LO numDof(0);
      LO numNode = m_subNodes[i].size();
      m_subNodeBegin[i].resize(numNode+1, 0);
      for (LO j=0; j<numNode; j++) {
	LO node = m_subNodes[i][j];
	numDof += m_nodeBegin[node+1] - m_nodeBegin[node];
	m_subNodeBegin[i][j+1] = numDof;
      }
      m_subLocalDofs[i].resize(numDof);
      numDof = 0;
      for (LO j=0; j<numNode; j++) {
	LO node = m_subNodes[i][j];
	for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	  m_subLocalDofs[i][numDof++] = m_localDofs[j];
	}
      }
    }
  }

  void initializeDisjointData()
  {
    if (m_disjointCommunicators == false) return;
    if (m_fineProc == true) {
      if (m_interfacePreconditioner == true) {
	determineNewMap(m_dofMapB, m_dofMapB_dis);
	determineNewMap(m_dofMapB1to1, m_dofMapB1to1_dis);
	m_rhsVecB_dis = rcp( new Vector(m_dofMapB_dis) );
	m_solVecB_dis = rcp( new Vector(m_dofMapB_dis) );
	m_rhsVecB1to1_dis = rcp( new Vector(m_dofMapB1to1_dis) );
	m_solVecB1to1_dis = rcp( new Vector(m_dofMapB1to1_dis) );
	m_exporterB_dis = rcp( new Export(m_dofMapB_dis, m_dofMapB1to1_dis) );
      }
      else {
	determineNewMap(m_dofMap, m_dofMap_dis);
	determineNewMap(m_dofMap1to1, m_dofMap1to1_dis);
	m_exporterAll_dis = rcp( new Export(m_dofMap1to1_dis, m_dofMap_dis) );
	m_rhsVecAll_dis = rcp( new Vector(m_dofMap_dis) );
	m_solVecAll_dis = rcp( new Vector(m_dofMap_dis) );
	m_rhsVecAll1to1_dis = rcp( new Vector(m_dofMap1to1_dis) );
	m_solVecAll1to1_dis = rcp( new Vector(m_dofMap1to1_dis) );
      }
    }
  }

  void initializeVectors()
  {
    m_rhsVecB = rcp( new Vector(m_dofMapB) );
    m_solVecB = rcp( new Vector(m_dofMapB) );
    m_rhsVecB1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_solVecB1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_xVecAll = rcp( new Vector(m_dofMap) );
    m_AxVecAll = rcp( new Vector(m_dofMap) );
    m_xVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_AxVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_rhsVecAll = rcp( new Vector(m_dofMap) );
    m_solVecAll = rcp( new Vector(m_dofMap) );
    m_rhsVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_solVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
  }

  bool checkUseDirectSolver()
  {
    bool useDirectSolver = false;
    GO numMySubs = m_subNodeBegin.size();
    GO numSubs(0);
    Teuchos::reduceAll<int, GO>(*m_Comm, Teuchos::REDUCE_SUM, 1,
				&numMySubs, &numSubs);
    assert (numSubs > 0);
    if (numSubs == 1) {
      useDirectSolver = true;
    }
    return useDirectSolver;
  }

  void initializeDirectSolverAndVectors()
  {
    assert (m_Comm->getSize() == 1);
    assert (m_subNodeBegin.size() == 1);
    LO numRows = m_numDofs;
    std::cout << "coarse space dimension = " << numRows << std::endl;
    SolverFactory<SX> Factory;
    LO* rowBegin = m_subRowBegin[0];
    LO* columns = m_subColumns[0];
    SX* values = m_subValues[0];
    m_directSolver = Factory.Generate(numRows, rowBegin, columns, values,
				      *m_Parameters);
    m_directSolver->Initialize();
  }

  void determineWeights()
  {
    m_Weights = 
      rcp( new bddc::WeightsBDDC<SX,SM,LO,GO>
	   (m_Subdomain, m_Partition, m_exporterB,
	    m_subBoundaryDofs, m_diagBoundary, m_Parameters) );
    m_Weights->determineWeights();
  }

  RCP<const Map> constructOneToOne(RCP<const Map> map, 
				   RCP<const Map> map1to1)
  {
    std::vector<GO> globalIDs;
    LO numRows = map->getNodeNumElements();
    for (LO i=0; i<numRows; i++) {
      GO globalID = map->getGlobalElement(i);
      LO localID = map1to1->getLocalElement(globalID);
      if (localID != Teuchos::OrdinalTraits<LO>::invalid()) {
	globalIDs.push_back(globalID);
      }
    }
    RCP<const Map> newMap1to1 = generateMap(globalIDs);
    return newMap1to1;
  }

  void determineBoundaryMaps()
  {
    std::vector<GO> globalIDs(m_numDofsB);
    for (LO i=0; i<m_numDofsB; i++) {
      LO dof = m_boundaryDofs[i];
      globalIDs[i] = m_dofMap->getGlobalElement(dof);
    }
    m_dofMapB = generateMap(globalIDs);
    m_dofMapB1to1 = constructOneToOne(m_dofMapB, m_dofMap1to1);
    //    Tpetra::createOneToOne<LO,GO,Node>(m_dofMapB);
    m_exporterB = rcp( new Export(m_dofMapB, m_dofMapB1to1) );
    m_numDofsB1to1 = m_dofMapB1to1->getNodeNumElements();
    m_exporterAll = rcp( new Export(m_dofMap, m_dofMap1to1) );
    m_numDofsB1to1 = m_dofMapB1to1->getNodeNumElements();
    m_boundaryToAll1to1.resize(m_numDofsB1to1);
    for (LO i=0; i<m_numDofsB1to1; i++) {
      GO globalID = m_dofMapB1to1->getGlobalElement(i);
      LO localID = m_dofMap1to1->getLocalElement(globalID);
      assert (localID != Teuchos::OrdinalTraits<LO>::invalid());
      m_boundaryToAll1to1[i] = localID;
    }
  }

  void determineInteriorMaps()
  {
    m_subInteriorDofs.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      const std::vector<LO> & interiorDofs = m_Subdomain[i]->getInteriorDofs();
      LO numInteriorDofs = interiorDofs.size();
      m_subInteriorDofs[i].resize(numInteriorDofs);
      for (LO j=0; j<numInteriorDofs; j++) {
	LO localID = m_subDofs[i][interiorDofs[j]];
	GO globalID = m_dofMap->getGlobalElement(localID);
	LO localID1to1 = m_dofMap1to1->getLocalElement(globalID);
	assert (localID1to1 != Teuchos::OrdinalTraits<LO>::invalid());
	m_subInteriorDofs[i][j] = localID1to1;
      }
    }
  }

  void determineBoundaryDofs()
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    const std::vector< std::vector<LO> > & equivClasses = 
      m_Partition->getEquivClasses();
    m_numSub = subdomainEquivClasses.size();
    m_boundaryDofsLocal.resize(m_numSub);
    m_subBoundaryDofs.resize(m_numSub);
    m_subDofs.resize(m_numSub);
    std::vector<LO> globalToLocalMap(m_numNodes, -1);
    std::vector<bool> boundaryFlag(m_numDofs, false);
    for (LO i=0; i<m_numSub; i++) {
      getBoundaryDofs(i, subdomainEquivClasses[i], equivClasses, 
		      globalToLocalMap, m_boundaryDofsLocal[i],
		      m_subBoundaryDofs[i]);
      getSubDofs(i, m_subDofs[i]);
      for (size_t j=0; j<m_subBoundaryDofs[i].size(); j++) {
	boundaryFlag[m_subBoundaryDofs[i][j]] = true;
      }
    }
    for (LO i=0; i<m_numDofs; i++) {
      if (boundaryFlag[i] == true) m_boundaryDofs.push_back(i);
    }
    m_numDofsB = m_boundaryDofs.size();
    m_globalToBoundaryMap.resize(m_numDofs, -1);
    for (LO i=0; i<m_numDofsB; i++) {
      m_globalToBoundaryMap[m_boundaryDofs[i]] = i;
    }
    for (LO i=0; i<m_numSub; i++) {
      for (size_t j=0; j<m_subBoundaryDofs[i].size(); j++) {
	LO dof = m_subBoundaryDofs[i][j];
	LO dofB = m_globalToBoundaryMap[dof];
	assert (dofB != -1);
	m_subBoundaryDofs[i][j] = dofB;
      }
    }
    determineEquivBoundaryDofs(equivClasses);
  }

  void determineEquivBoundaryDofs
    (const std::vector< std::vector<LO> > & equivClasses)
  {
    LO numEquiv = equivClasses.size();
    m_equivBoundaryDofs.resize(numEquiv);
    for (LO i=0; i<numEquiv; i++) {
      for (size_t j=0; j<equivClasses[i].size(); j++) {
	LO node = equivClasses[i][j];
	for (LO k=m_nodeBegin[node]; k<m_nodeBegin[node+1]; k++) {
	  LO dofB = m_globalToBoundaryMap[k];
	  assert (dofB != -1);
	  m_equivBoundaryDofs[i].push_back(dofB);
	}
      }
    }
  }

  void determineCoarseMpiRanks(int & myCoarseMpiRank,
			       std::vector<int> & mpiRanks)
  {
    int numSubdomainsPerCoarseSubdomain = 
      m_Parameters->get("numSubdomainsPerCoarseSubdomain", 512);
    int numCoarseSubdomainsPerMpiRank = 
      m_Parameters->get("numCoarseSubdomainsPerMpiRank", 1);
    int numSubdomainsPerMpiRank = 
      numCoarseSubdomainsPerMpiRank*numSubdomainsPerCoarseSubdomain;
    int numMySubdomains = m_Partition->getNumberOfMySubdomains();
    int totalNumSubdomains;
    Teuchos::reduceAll<int, int> 
      (*m_Comm, Teuchos::REDUCE_SUM, 1, &numMySubdomains, &totalNumSubdomains);
    int averageNumSubdomainsPerMpiRank = totalNumSubdomains/m_Comm->getSize();
    numSubdomainsPerMpiRank = std::max(numSubdomainsPerMpiRank,
				       averageNumSubdomainsPerMpiRank);
    int numMpiRanks = totalNumSubdomains/numSubdomainsPerMpiRank;
    if (totalNumSubdomains % numSubdomainsPerMpiRank != 0) numMpiRanks++;
    numMpiRanks = std::min(numMpiRanks, m_Comm->getSize());
    mpiRanks.resize(numMpiRanks);
    determineCoarseProcs(mpiRanks);
    int myRankActive(0);
    for (int i=0; i<numMpiRanks; i++) {
      if (mpiRanks[i] == m_myPID) myRankActive = 1;
    }
    int myRankActiveSS;
    Teuchos::scan<int, int> (*m_Comm, Teuchos::REDUCE_SUM, 1, &myRankActive, 
			     &myRankActiveSS);
    myCoarseMpiRank = -1;
    if (myRankActive == 1) {
      myCoarseMpiRank = myRankActiveSS - 1;
    }
  }

  SM getProcessorLoad()
  {
    SM procLoad = 0;
    for (LO i=0; i<m_numSub; i++) {
      LO numNodes = m_subNodes[i].size();
      LO numRows = m_subNodeBegin[i][numNodes];
      SM numTerms = m_subRowBegin[i][numRows];
      procLoad += numTerms;
    }
    return procLoad;
  }

  void determineCoarseProcs(std::vector<int> & mpiRanks)
  {
    SM procLoad = getProcessorLoad();
    int numProc = m_Comm->getSize();
    int numCoarseProc = mpiRanks.size();
    std::vector<SM> procLoadAll(numProc);
    Teuchos::gatherAll<int, SM> (*m_Comm, 1, &procLoad, numProc, 
				 procLoadAll.data());
    std::vector< std::pair<SM, int> > procSize(numProc);
    for (int i=0; i<numProc; i++) {
      procSize[i] = std::make_pair(procLoadAll[i], i);
    }
    std::sort(procSize.begin(), procSize.end());
    for (int i=0; i<numCoarseProc; i++) {
      auto iter = procSize.begin() + i;
      int proc = iter->second;
      mpiRanks[i] = proc;
    }
  }

  void determineSubdomainCoords(std::vector<double> & coords)
  {
    coords.resize(3*m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      double sumX(0), sumY(0), sumZ(0);
      LO numNode = m_subNodes[i].size();
      for (LO j=0; j<numNode; j++) {
	LO node = m_subNodes[i][j];
	sumX += m_xCoord[node];
	sumY += m_yCoord[node];
	sumZ += m_zCoord[node];
      }
      coords[i+0*m_numSub] = sumX/numNode;
      coords[i+1*m_numSub] = sumY/numNode;
      coords[i+2*m_numSub] = sumZ/numNode;
    }
  }

  void determineSubdomainMap(RCP<const Map> & subdomainMap)
  {
    GO startingSub = m_Partition->getStartingSub();
    std::vector<GO> subGIDs(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      subGIDs[i] = startingSub + i;
    }
    subdomainMap = generateMap(subGIDs); 
  }

  void determineSubdomainCoordsForMpiRank
    (RCP<const CrsGraph> subsForCoarseProc,
     std::vector<double> & subCoordsForCoarse)
  {
    RCP<const Map> sourceMap;
    determineSubdomainMap(sourceMap);
    MVD coordsVec(sourceMap, 3);
    Teuchos::ArrayRCP<double> xCoords = coordsVec.getDataNonConst(0);
    Teuchos::ArrayRCP<double> yCoords = coordsVec.getDataNonConst(1);
    Teuchos::ArrayRCP<double> zCoords = coordsVec.getDataNonConst(2);
    std::vector<double> subdomainCoords;
    determineSubdomainCoords(subdomainCoords);
    for (LO i=0; i<m_numSub; i++) {
      xCoords[i] = subdomainCoords[i+0*m_numSub];
      yCoords[i] = subdomainCoords[i+1*m_numSub];
      zCoords[i] = subdomainCoords[i+2*m_numSub];
    }
    RCP<const Map> targetMap = subsForCoarseProc->getColMap();
    MVD coordsVecTarget(targetMap, 3); 
    Export exporter(sourceMap, targetMap);
    coordsVecTarget.doExport(coordsVec, exporter, Tpetra::INSERT);
    Teuchos::ArrayRCP<double> xCoordsTarget = 
      coordsVecTarget.getDataNonConst(0);
    Teuchos::ArrayRCP<double> yCoordsTarget = 
      coordsVecTarget.getDataNonConst(1);
    Teuchos::ArrayRCP<double> zCoordsTarget = 
      coordsVecTarget.getDataNonConst(2);
    LO numSubTarget = targetMap->getNodeNumElements();
    subCoordsForCoarse.resize(3*numSubTarget);    
    for (LO i=0; i<numSubTarget; i++) {
      subCoordsForCoarse[i+0*numSubTarget] = xCoordsTarget[i];
      subCoordsForCoarse[i+1*numSubTarget] = yCoordsTarget[i];
      subCoordsForCoarse[i+2*numSubTarget] = zCoordsTarget[i];
    }
  }
  
  void getConnectivityGraphForCoarse
    (RCP<const CrsGraph> subsForCoarseProc,
     RCP<CrsGraph> & Graph)
  {
    RCP<const Map> targetMap = subsForCoarseProc->getColMap();
    CrsGraph targetGraph(targetMap, 0);
    RCP<const CrsGraph> sourceGraph = m_Partition->getConnectivityGraph();
    Import importer(sourceGraph->getRowMap(), targetMap);
    targetGraph.doImport(*sourceGraph, importer, Tpetra::INSERT);
    targetGraph.fillComplete(sourceGraph->getDomainMap(),
			     sourceGraph->getRangeMap());
    LO numRows = targetMap->getNodeNumElements();
    if (numRows > 0) {
      MPI_Comm mpiComm = MPI_COMM_SELF;
      RCP<const Teuchos::Comm<int> > Comm = 
	rcp( new Teuchos::MpiComm<int>(mpiComm) );
      std::vector<GO> rowIDs(numRows);
      Teuchos::ArrayRCP<size_t> count(numRows);
      Teuchos::ArrayView<const LO> Indices;
      for (LO i=0; i<numRows; i++) {
	rowIDs[i] = i;
	targetGraph.getLocalRowView(i, Indices);
	LO numLocal(0);
	for (int j=0; j<Indices.size(); j++) {
	  GO globalID = targetGraph.getColMap()->getGlobalElement(Indices[j]);
	  LO localID = targetMap->getLocalElement(globalID);
	  if (localID !=  Teuchos::OrdinalTraits<LO>::invalid()) numLocal++;
	}
	count[i] = numLocal;
      }
      RCP<const Map> rowMap = generateMap(rowIDs);
      Graph = rcp( new CrsGraph(rowMap, rowMap, count, Tpetra::StaticProfile) );
      for (LO i=0; i<numRows; i++) {
	targetGraph.getLocalRowView(i, Indices);
	std::vector<LO> indices;
	for (int j=0; j<Indices.size(); j++) {
	  GO globalID = targetGraph.getColMap()->getGlobalElement(Indices[j]);
	  LO localID = targetMap->getLocalElement(globalID);
	  if (localID !=  Teuchos::OrdinalTraits<LO>::invalid()) {
	    indices.push_back(localID);
	  }
	}
	Graph->insertLocalIndices(i, Teuchos::ArrayView<LO>(indices));
      }
      Graph->fillComplete();
    }
  }
  
  void partitionSubdomainsWithinMpiRanks
    (RCP<const CrsGraph> subsForCoarseProc,
     std::vector< std::vector<GO> > & coarseSubdomainSubs)
  {
    LO numSubdomains = subsForCoarseProc->getNodeNumCols();
    int numSubdomainsPerCoarseSubdomain = 
      m_Parameters->get("numSubdomainsPerCoarseSubdomain", 512);
    int numParts = numSubdomains/numSubdomainsPerCoarseSubdomain;
    if (numSubdomains % numSubdomainsPerCoarseSubdomain != 0) numParts++;
    std::vector<GO> subdomainGIDs(numSubdomains);
    for (LO i=0; i<numSubdomains; i++) {
      subdomainGIDs[i] = subsForCoarseProc->getColMap()->getGlobalElement(i);
    }
    int maxNumParts(0);
    Teuchos::reduceAll<int, int>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				 &numParts, &maxNumParts);
    // quick exit if no further partitioning required
    if (maxNumParts < 2) {
      if (numParts == 1) {
	coarseSubdomainSubs.resize(1);
	coarseSubdomainSubs[0] = subdomainGIDs;
      }
      return;
    }
    std::string coarsenOption = m_Parameters->get("Coarsening Option", "Graph");
    std::vector<double> elemCoords;
    RCP<CrsGraph> Graph;
    if (coarsenOption == "Graph") {
      getConnectivityGraphForCoarse(subsForCoarseProc, Graph);
    }
    else {
      determineSubdomainCoordsForMpiRank(subsForCoarseProc, elemCoords);
    }
    if (numParts == 1) {
      coarseSubdomainSubs.resize(1);
      coarseSubdomainSubs[0] = subdomainGIDs;
      return;
    }
    std::vector<LO> parts(numSubdomains, 0);
    RCP<Teuchos::ParameterList> params = rcp( new Teuchos::ParameterList() );
    params->set("Number of Parts", numParts);
    ZoltanPartition<LO, GO>* partition(nullptr);
    params->set("Coordinates", elemCoords.data());
    if (coarsenOption == "Graph") {
      params->set("LB Method", "Graph");
      partition = new ZoltanPartition<LO, GO>(Graph, params);
    }
    else {
      if (coarsenOption == "Recursive Coordinate Bisection") {
	params->set("LB Method", "RCB");
      }
      else if (coarsenOption == "Recursive Inertial Bisection") {
	params->set("LB Method", "RIB");
      }
      LO *rowBegin(0), *columns(0);
      std::vector<GO> globalIDs(numSubdomains);
      for (LO i=0; i<numSubdomains; i++) globalIDs[i] = i;
      MPI_Comm mpiComm = MPI_COMM_SELF;
      partition = new ZoltanPartition<LO, GO>
	(numSubdomains, rowBegin, columns, elemCoords.data(), 
	 globalIDs.data(), mpiComm, params);
    }
    partition->setDefaultWeights();
    partition->doPartition();
    const int* partNumbers = partition->getParts();
    coarseSubdomainSubs.resize(numParts);
    for (int i=0; i<m_numSub; i++) {
      const int part = partNumbers[i];
      coarseSubdomainSubs[part].push_back(subdomainGIDs[i]);
    }
    delete partition;
  }

  void partitionSubdomainsAcrossMpiRanks(const std::vector<int> & mpiRanks,
					 std::vector<GO> & parts)
  {
    parts.resize(m_numSub, 0);
    std::vector<double> subdomainCoords;
    determineSubdomainCoords(subdomainCoords);
    if (mpiRanks.size() > 1) {
      std::vector<GO> globalIDs(m_numSub);
      int startingSub = m_Partition->getStartingSub();
      for (int i=0; i<m_numSub; i++) globalIDs[i] = startingSub + i;
      RCP<Teuchos::ParameterList> params = rcp( new Teuchos::ParameterList() );
      int numCoarseRanks = mpiRanks.size();
      params->set("Number of Parts", numCoarseRanks);
      params->set("Coordinates", subdomainCoords.data());
      ZoltanPartition<LO, GO>* partition(nullptr);
      std::string coarsenOption = 
	m_Parameters->get("Coarsening Option", "Graph");
      if (coarsenOption == "Graph") {
	RCP<const CrsGraph> Graph = m_Partition->getConnectivityGraph();
	params->set("LB Method", "Graph");
	partition = new ZoltanPartition<LO, GO>(Graph, params);
      }
      else {
	if (coarsenOption == "Recursive Coordinate Bisection") {
	  params->set("LB Method", "RCB");
	}
	else if (coarsenOption == "Recursive Inertial Bisection") {
	  params->set("LB Method", "RIB");
	}
	LO *rowBegin(0), *columns(0);
	partition = new ZoltanPartition<LO, GO>
	  (m_numSub, rowBegin, columns, subdomainCoords.data(), 
	   globalIDs.data(), m_mpiComm, params);
      }
      partition->setDefaultWeights();
      partition->doPartition();
      const int* partNumbers = partition->getParts();
      for (int i=0; i<m_numSub; i++) parts[i] = partNumbers[i];
      delete partition;
    }
  }
  
  void determineSubdomainsInCoarsePartition
    (int myCoarseMpiRank, 
     const std::vector<GO> & parts, 
     RCP<CrsGraph> & subsForCoarseProc)
  {
    std::vector<GO> uniqueParts;
    DofManager<LO,GO>::determineUniqueIndices(parts, uniqueParts);
    RCP<const Map> partMap = generateMap(uniqueParts);
    Teuchos::ArrayRCP<size_t> count(uniqueParts.size());
    std::vector< std::vector<LO> > localIndices(uniqueParts.size());
    for (size_t i=0; i<parts.size(); i++) {
      LO index = partMap->getLocalElement(parts[i]);
      localIndices[index].push_back(i);
      count[index]++;
    }
    RCP<const Map> subMap;
    determineSubdomainMap(subMap);
    CrsGraph subsForParts(partMap, subMap, count, Tpetra::StaticProfile);
    for (size_t i=0; i<uniqueParts.size(); i++) {
      subsForParts.insertLocalIndices
	(i, Teuchos::ArrayView<LO>(localIndices[i]));
    }
    std::vector<GO> coarseParts;
    if (myCoarseMpiRank != -1) coarseParts.resize(1, myCoarseMpiRank);
    RCP<const Map> coarsePartsMap = generateMap(coarseParts);
    subsForParts.fillComplete(subMap, coarsePartsMap);
    subsForCoarseProc = rcp( new CrsGraph(coarsePartsMap, 0) );
    Export exporter(partMap, coarsePartsMap);
    subsForCoarseProc->doExport(subsForParts, exporter, Tpetra::ADD);
    subsForCoarseProc->fillComplete(subMap, coarsePartsMap);
  }
  
  void getCoarseNodeLocalDofs
    (std::vector< std::vector<LO> > & coarseNodeLocalDofs)
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    LO numCoarseNodes = m_Partition->getNumEquivClasses();
    coarseNodeLocalDofs.resize(numCoarseNodes);
    std::vector<bool> equivProcessed(numCoarseNodes, false);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      for (LO j=0; j<numEquiv; j++) {
	LO equiv = subdomainEquivClasses[i][j];
	if (equivProcessed[equiv] == false) {
	  coarseNodeLocalDofs[equiv] = m_Subdomain[i]->getEquivLocalDofs(j); 
	  equivProcessed[equiv] = true;
	}
      }
    }
  }

  void determineActiveCoarseNodeCoordinatesAndLocalDofs
    (std::vector<LO> & activeCoarseNodes,
     std::vector<LO> & activeCoarseNodesNumDofs,
     RCP<MVD> & coarseCoordsVec1to1,
     RCP<CrsGraph> & localDofsGraph1to1,
     RCP<Export> & nodeExporter)
  {
    // first, determine active coarse nodes
    const std::vector<GO> coarseNodeGIDs = m_Partition->getGlobalIDs();
    std::vector< std::vector<LO> > coarseNodeLocalDofs;
    getCoarseNodeLocalDofs(coarseNodeLocalDofs);
    assert (coarseNodeLocalDofs.size() == coarseNodeGIDs.size());
    LO numCoarseNodes = coarseNodeGIDs.size();
    std::vector<GO> activeCoarseNodeGIDs;
    for (LO i=0; i<numCoarseNodes; i++) {
      if (coarseNodeLocalDofs[i].size() > 0) {
	activeCoarseNodes.push_back(i);
	activeCoarseNodeGIDs.push_back(coarseNodeGIDs[i]);
      }
    }
    LO numActiveCoarseNodes = activeCoarseNodes.size();
    // next, determine coarse node coordinates
    std::vector<SM> coarseCoords;
    getCoordsCoarseNode(activeCoarseNodes, coarseCoords);
    RCP<const Map> coarseNodeMap = generateMap(activeCoarseNodeGIDs);
    MVD coarseCoordsVec(coarseNodeMap, 3);
    Teuchos::ArrayRCP<double> xCoords = coarseCoordsVec.getDataNonConst(0);
    Teuchos::ArrayRCP<double> yCoords = coarseCoordsVec.getDataNonConst(1);
    Teuchos::ArrayRCP<double> zCoords = coarseCoordsVec.getDataNonConst(2);
    std::vector<double> subdomainCoords;
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      xCoords[i] = coarseCoords[i+0*numActiveCoarseNodes];
      yCoords[i] = coarseCoords[i+1*numActiveCoarseNodes];
      zCoords[i] = coarseCoords[i+2*numActiveCoarseNodes];
    }
    RCP<const Map> coarseNodeMap1to1 = 
      Tpetra::createOneToOne<LO,GO,Node>(coarseNodeMap);
    coarseCoordsVec1to1 = rcp( new MVD(coarseNodeMap1to1, 3) );
    nodeExporter = rcp( new Export(coarseNodeMap, coarseNodeMap1to1) );
    coarseCoordsVec1to1->doExport(coarseCoordsVec, *nodeExporter, 
				  Tpetra::INSERT);
    // next, determine localDofs graph
    std::vector<GO> allLocalDofs, allLocalDofsUnique;
    Teuchos::ArrayRCP<size_t> count(numActiveCoarseNodes);
    LO numCoarseDof(0);
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      LO index = activeCoarseNodes[i];
      const std::vector<LO> & localDofs = coarseNodeLocalDofs[index];
      count[i] = localDofs.size();
      numCoarseDof += count[i];
      for (size_t j=0; j<localDofs.size(); j++) {
	allLocalDofs.push_back(localDofs[j]);
      }
    }
    DofManager<LO,GO>::determineUniqueIndices(allLocalDofs, allLocalDofsUnique);
    std::map<GO, LO> mapLocalDofs;
    for (size_t i=0; i<allLocalDofsUnique.size(); i++) {
      mapLocalDofs.insert(std::make_pair(allLocalDofsUnique[i], i));
    }
    RCP<const Map> colMap = generateMap(allLocalDofsUnique);
    CrsGraph localDofsGraph(coarseNodeMap, colMap, count, 
			    Tpetra::StaticProfile);
    std::vector<LO> indices;
    activeCoarseNodesNumDofs.resize(numActiveCoarseNodes);
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      LO index = activeCoarseNodes[i];
      const std::vector<LO> & localDofs = coarseNodeLocalDofs[index];
      activeCoarseNodesNumDofs[i] = localDofs.size();
      indices.resize(localDofs.size());
      for (size_t j=0; j<localDofs.size(); j++) {
	auto iter = mapLocalDofs.find(localDofs[j]);
	assert (iter != mapLocalDofs.end());
	indices[j] = iter->second;
      }
      localDofsGraph.insertLocalIndices(i, Teuchos::ArrayView<LO>(indices));
    }
    RCP<const Map> colMap1to1 =
      Tpetra::createOneToOne<LO,GO,Node>(colMap);
    localDofsGraph.fillComplete(colMap1to1, coarseNodeMap1to1);
    localDofsGraph1to1 = rcp( new CrsGraph(coarseNodeMap1to1, 0) );
    localDofsGraph1to1->doExport(localDofsGraph, *nodeExporter, Tpetra::INSERT);
    localDofsGraph1to1->fillComplete(colMap1to1, coarseNodeMap1to1);
  }
  
  RCP<const Map> generateMap(const std::vector<GO> & globalIDs)
  {
    return
      rcp( new Map(m_IGO, Teuchos::ArrayView<const GO>(globalIDs), 0, m_Comm) );
  }

  void determineCoarseMapNextLevel
    (RCP<const VectorGO> startNodeVector1to1,
     RCP<const VectorGO> numDofNodeVector1to1,
     RCP<const Map> & coarseMapNextLevel)
  {
    RCP<const Map> nodeMap = generateMap(m_nodeGlobalIDsCoarse);
    RCP<const Map> targetMap = Tpetra::createOneToOne<LO,GO,Node>(nodeMap);
    RCP<const Map> sourceMap = startNodeVector1to1->getMap();
    VectorGO startNodeVector(targetMap);
    VectorGO numDofNodeVector(targetMap);
    Import importer(sourceMap, targetMap);
    startNodeVector.doImport(*startNodeVector1to1, importer, Tpetra::INSERT);
    numDofNodeVector.doImport(*numDofNodeVector1to1, importer, Tpetra::INSERT);
    Teuchos::ArrayRCP<const GO> startValues = startNodeVector.getData();
    Teuchos::ArrayRCP<const GO> numDofValues = numDofNodeVector.getData();
    std::vector<GO> globalIDs;
    LO numTarget = targetMap->getNodeNumElements();
    m_nodeGlobalIDsCoarse1to1.resize(numTarget);
    for (LO i=0; i<numTarget; i++) {
      m_nodeGlobalIDsCoarse1to1[i] = targetMap->getGlobalElement(i);
      GO start = startValues[i];
      for (GO j=0; j<numDofValues[i]; j++) {
	globalIDs.push_back(start+j);
      }
    }
    coarseMapNextLevel = generateMap(globalIDs);
  }

  void determineCoarseMap
    (RCP<const CrsGraph> localDofsGraph1to1,
     RCP<const Export> nodeExporter,
     const std::vector<LO> & activeCoarseNodesNumDofs,
     RCP<VectorGO> & startNodeVector1to1,
     RCP<VectorGO> & numDofNodeVector1to1,
     RCP<const Map> & coarseMap)
  {
    determineCoarseNodeData(localDofsGraph1to1, startNodeVector1to1,
			    numDofNodeVector1to1);
    std::vector<GO> coarseDofGIDs;
    determineCoarseDofGIDs(startNodeVector1to1, nodeExporter, 
			   activeCoarseNodesNumDofs, coarseDofGIDs);
    coarseMap = generateMap(coarseDofGIDs);
  }

  void determineCoarseNodeData
    (RCP<const CrsGraph> localDofsGraph1to1, 
     RCP<VectorGO> & startNodeVector1to1,
     RCP<VectorGO> & numDofNodeVector1to1)
  {
    GO numDof1to1 = localDofsGraph1to1->getNodeNumEntries();
    GO numDof1to1SS;
    Teuchos::scan<int, GO> (*m_Comm, Teuchos::REDUCE_SUM, 1, &numDof1to1, 
			     &numDof1to1SS);
    GO start = numDof1to1SS - numDof1to1;
    RCP<const Map> coarseNodeMap1to1 = localDofsGraph1to1->getRowMap();
    startNodeVector1to1 = rcp( new VectorGO(coarseNodeMap1to1) );
    numDofNodeVector1to1 = rcp( new VectorGO(coarseNodeMap1to1) );
    Teuchos::ArrayRCP<GO> startVals1to1 = 
      startNodeVector1to1->getDataNonConst();
    Teuchos::ArrayRCP<GO> numDofVals1to1 = 
      numDofNodeVector1to1->getDataNonConst();
    Teuchos::ArrayView<const LO> Indices;
    for (size_t i=0; i<coarseNodeMap1to1->getNodeNumElements(); i++) {
      localDofsGraph1to1->getLocalRowView(i, Indices);
      startVals1to1[i] = start;
      numDofVals1to1[i] = Indices.size();
      start += Indices.size();
    }
  }

  void determineCoarseDofGIDs
    (RCP<const VectorGO> startNodeVector1to1, 
     RCP<const Export> nodeExporter,
     const std::vector<LO> & activeCoarseNodesNumDofs,
     std::vector<GO> & coarseDofGIDs)
  {
    RCP<const Map> coarseNodeMap = nodeExporter->getSourceMap();
    VectorGO startNodeVector(coarseNodeMap);
    startNodeVector.doImport(*startNodeVector1to1, *nodeExporter, 
			     Tpetra::INSERT);
    Teuchos::ArrayRCP<const GO> startVals = startNodeVector.getData();
    LO numCoarseDof(0);
    for (size_t i=0; i<activeCoarseNodesNumDofs.size(); i++) {
      numCoarseDof += activeCoarseNodesNumDofs[i];
    }
    coarseDofGIDs.resize(numCoarseDof);
    numCoarseDof = 0;
    for (size_t i=0; i<activeCoarseNodesNumDofs.size(); i++) {
      for (LO j=0; j<activeCoarseNodesNumDofs[i]; j++) {
	coarseDofGIDs[numCoarseDof++] = startVals[i] + j;
      }
    }
  }

  void determineCoarseNodesForSubdomainsAndSubdomainCoarseDofs
    (const std::vector<LO> & activeCoarseNodes,
     const std::vector<LO> & activeCoarseNodesNumDofs,
     RCP<CrsGraph> & coarseNodesForSubdomains)
  {
    // first, determine Map containing globalIDs of the original subdomains
    RCP<const Map> subMap;
    determineSubdomainMap(subMap);
    // next, determine CrsGraph containing active coarse nodes of the
    // original subdomains
    const std::vector< std::vector<LO> > subdomainCoarseNodes = 
      m_Partition->getSubdomainEquivClasses();
    const std::vector<GO> coarseNodeGIDs = m_Partition->getGlobalIDs();
    LO numActiveCoarseNodes = activeCoarseNodes.size();
    std::vector<GO> activeCoarseNodeGIDs(numActiveCoarseNodes);
    std::map<GO, LO> mapCoarseNodeGIDs;
    std::vector<LO> nodeBeginCoarse(numActiveCoarseNodes+1, 0);
    LO numDofCoarse(0);
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      numDofCoarse += activeCoarseNodesNumDofs[i];
      nodeBeginCoarse[i+1] = numDofCoarse;
      LO index = activeCoarseNodes[i];
      GO coarseNode = coarseNodeGIDs[index];
      activeCoarseNodeGIDs[i] = coarseNode;
      mapCoarseNodeGIDs.insert(std::make_pair(coarseNode, i));
    }
    std::vector<LO> indexMap(coarseNodeGIDs.size(), -1);
    for (size_t i=0; i<coarseNodeGIDs.size(); i++) {
      GO coarseNode = coarseNodeGIDs[i];
      auto iter = mapCoarseNodeGIDs.find(coarseNode);
      if (iter != mapCoarseNodeGIDs.end()) {
	indexMap[i] = iter->second;
      }
    }
    Teuchos::ArrayRCP<size_t> count(m_numSub);
    std::vector< std::vector<LO> > indices(m_numSub);
    m_subDofsCoarse.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      count[i] = 0;
      for (size_t j=0; j<subdomainCoarseNodes[i].size(); j++) {
	LO index = indexMap[subdomainCoarseNodes[i][j]];
	if (index != -1) {
	  indices[i].push_back(index);
	  count[i]++;
	  for (LO k=nodeBeginCoarse[index]; k<nodeBeginCoarse[index+1]; k++) {
	    m_subDofsCoarse[i].push_back(k);
	  }
	}
      }
      m_subDofsCoarse[i].shrink_to_fit();
    }
    RCP<const Map> activeCoarseNodeMap = generateMap(activeCoarseNodeGIDs);
    coarseNodesForSubdomains = rcp
      ( new CrsGraph(subMap, activeCoarseNodeMap, count, 
		     Tpetra::StaticProfile) );
    for (LO i=0; i<m_numSub; i++) {
      coarseNodesForSubdomains->insertLocalIndices
	(i, Teuchos::ArrayView<const LO>(indices[i]));
    }
    RCP<const Map> activeCoarseNodeMap1to1 = 
      Tpetra::createOneToOne<LO,GO,Node>(activeCoarseNodeMap);
    coarseNodesForSubdomains->fillComplete(activeCoarseNodeMap1to1, subMap);
  }
  
  void getCoordsCoarseNode(const std::vector<LO> & activeCoarseNodes,
			   std::vector<SM> & coordsCoarseNode)
  {
    const std::vector< std::vector<LO> > & equivClasses = 
      m_Partition->getEquivClasses();
    LO numActiveCoarseNodes = activeCoarseNodes.size();
    coordsCoarseNode.resize(3*numActiveCoarseNodes);
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      LO index = activeCoarseNodes[i];
      LO numNode = equivClasses[index].size();
      SM sumX(0), sumY(0), sumZ(0);
      for (LO j=0; j<numNode; j++) {
	LO node = equivClasses[index][j];
	sumX += m_xCoord[node];
	sumY += m_yCoord[node];
	sumZ += m_zCoord[node];
      }
      coordsCoarseNode[i+0*numActiveCoarseNodes] = sumX/numNode;
      coordsCoarseNode[i+1*numActiveCoarseNodes] = sumY/numNode;
      coordsCoarseNode[i+2*numActiveCoarseNodes] = sumZ/numNode;
    }
  }
  
  void getSubdomainParts(const std::vector<GO> & partsWithinMpiRanks, 
			 std::vector< std::vector<LO> > & partsForSubs)
  {
    LO numSubs(0);
    for (size_t i=0; i<partsWithinMpiRanks.size(); i++) {
      GO part = partsWithinMpiRanks[i];
      if (part+1 > numSubs) numSubs = part+1;
    }
    partsForSubs.resize(numSubs);
    for (size_t i=0; i<partsWithinMpiRanks.size(); i++) {
      GO part = partsWithinMpiRanks[i];
      partsForSubs[part].push_back(i);
    }
  }

  void allocateCoarseEntities()
  {
    // importer, exporter, vectors and work memory for coarse level
    m_exporterCoarse = rcp( new Export(m_coarseMap, m_coarseMapNextLevel) );
    m_coarseVec = rcp( new Vector(m_coarseMap) );
    m_coarseVecNextLevel = rcp( new Vector(m_coarseMapNextLevel) );
    if (m_interfacePreconditioner == true) {
      m_fineVecB = rcp( new Vector(m_dofMapB) );
      m_fineVecB1to1 = rcp( new Vector(m_dofMapB1to1) );
    }
    else {
      m_fineVec = rcp( new Vector(m_dofMap) );
      m_fineVec1to1 = rcp( new Vector(m_dofMap1to1) );
    }
    LO numCoarseDof = m_coarseMapNextLevel->getNodeNumElements();
    m_coarseRhsWork.resize(numCoarseDof);
    m_coarseSolWork.resize(numCoarseDof);
  }

  void applyPhi(RCP<Vector> & coarseVectorNextLevel,
		RCP<Vector> & fineVector1to1,
		SX beta,
		bool restrictToBoundary)
  {
    RCP<Vector> coarseVector = m_coarseVec;
    coarseVector->doImport(*coarseVectorNextLevel, *m_exporterCoarse,
			   Tpetra::INSERT);

    Teuchos::ArrayRCP<const SX> xVals = coarseVector->getData();
    RCP<Vector> fineVector;
    RCP<Export> exporter;
    std::vector< std::vector<LO> > & subDofs = m_subBoundaryDofs;
    if (restrictToBoundary == true) {
      fineVector = m_fineVecB;
      exporter = m_exporterB;
    }
    else {
      fineVector = m_fineVec;
      exporter = m_exporterAll;
      subDofs = m_subDofs;
    }
    fineVector->putScalar(0);
    Teuchos::ArrayRCP<SX> bVals = fineVector->getDataNonConst();
    for (LO i=0; i<m_numSub; i++) {
      SX *x(0), *b(0);
      m_Subdomain[i]->getSubdomainVectors(x, b);
      const std::vector<LO> & subDofsCoarse = m_subDofsCoarse[i];
      LO numRows = subDofsCoarse.size();
      LO numDof(0);
      for (LO j=0; j<numRows; j++) {
	LO row = subDofsCoarse[j];
	x[numDof++] = xVals[row];
      }
      bool transpose(false);
      m_Subdomain[i]->multiplyByPhi(x, b, restrictToBoundary, transpose);
      const std::vector<LO> & dofs = subDofs[i];
      LO numDofs = dofs.size();
      for (LO j=0; j<numDofs; j++) {
	bVals[dofs[j]] += b[j];
      }
    }
    // The following two lines of code doesn't do what I expect. It seems
    // like the values in fineVector1to1 get zeroed out even for beta = 1.
    // That is, the doExport seems to ignore initial values in fineVector1to1.
    // Work around is ugly, but gets the job done.
    /*
    fineVector1to1->scale(beta);
    fineVector1to1->doExport(*fineVector, *exporter, Tpetra::ADD);
    */
    Teuchos::ArrayRCP<SX> fineValues1to1 = fineVector1to1->getDataNonConst();
    LO numRows = fineVector1to1->getMap()->getNodeNumElements();
    m_fine1to1Work.resize(numRows);
    for (LO i=0; i<numRows; i++) {
      m_fine1to1Work[i] = beta*fineValues1to1[i];
    }
    fineVector1to1->putScalar(0);
    fineVector1to1->doExport(*fineVector, *exporter, Tpetra::ADD);
    for (LO i=0; i<numRows; i++) {
      fineValues1to1[i] += m_fine1to1Work[i];
    }
  }

  void applyPhiTranspose(RCP<Vector> & fineVector1to1,
			 RCP<Vector> & coarseVectorNextLevel,
			 bool restrictToBoundary)
  {
    RCP<Vector> fineVector;
    RCP<Export> exporter;
    std::vector< std::vector<LO> > & subDofs = m_subBoundaryDofs;
    if (restrictToBoundary == true) {
      fineVector = m_fineVecB;
      exporter = m_exporterB;
    }
    else {
      fineVector = m_fineVec;
      exporter = m_exporterAll;
      subDofs = m_subDofs;
    }

    Teuchos::ArrayRCP<const SX> values = fineVector1to1->getData();
    fineVector->doImport(*fineVector1to1, *exporter, Tpetra::INSERT);
    RCP<Vector> coarseVector = m_coarseVec;
    coarseVector->putScalar(0);
    Teuchos::ArrayRCP<SX> bVals = coarseVector->getDataNonConst();
    Teuchos::ArrayRCP<const SX> xVals = fineVector->getData();
    for (LO i=0; i<m_numSub; i++) {
      SX *x(0), *b(0);
      m_Subdomain[i]->getSubdomainVectors(x, b);
      const std::vector<LO> & dofs = subDofs[i];
      LO numDofs = dofs.size();
      for (LO j=0; j<numDofs; j++) {
	x[j] = xVals[dofs[j]];
      }
      bool transpose(true);
      m_Subdomain[i]->multiplyByPhi(x, b, restrictToBoundary, transpose);
      const std::vector<LO> & subDofsCoarse = m_subDofsCoarse[i];
      LO numRows = subDofsCoarse.size();
      LO numDof(0);
      for (LO j=0; j<numRows; j++) {
	LO row = subDofsCoarse[j];
	bVals[row] += b[numDof++];
      }
    }
    coarseVectorNextLevel->putScalar(0);
    coarseVectorNextLevel->doExport(*coarseVector, *m_exporterCoarse,
				    Tpetra::ADD);
  }

  void printCoarseMatrices()
  {
    LO numCoarseSub = m_subRowBeginCoarse.size();
    for (LO i=0; i<numCoarseSub; i++) {
      char fname[101];
      sprintf(fname, "Ac_%d_%d.dat", m_myPID, i);
      int numDof = m_subRowBeginCoarse[i].size() - 1;
      int* rowBegin = m_subRowBeginCoarse[i].data();
      int* columns = m_subColumnsCoarse[i].data();
      SX* values = m_subValuesCoarse[i].data();
      UtilBDDC<SX,SM>::printSparseMatrix
	(numDof, rowBegin, columns, values, fname);
    }
  }

  LO getNumDofCoarseSub(const LO sub)
  {
    LO numNode = m_subNodesCoarse[sub].size();
    LO numDof(0);
    for (LO i=0; i<numNode; i++) {
      LO node = m_subNodesCoarse[sub][i];
      numDof += m_nodeBeginCoarse[node+1] - m_nodeBeginCoarse[node];
    }
    return numDof;
  }

  void extractSubdomainMatrices(RCP<const CrsMatrix> coarseMatrix)
  {
    LO numSub = m_subNodesCoarse.size();
    m_subRowBeginCoarse.resize(numSub);
    m_subColumnsCoarse.resize(numSub);
    m_subValuesCoarse.resize(numSub);
    LO startRow(0);
    for (LO i=0; i<numSub; i++) {
      LO numDof = getNumDofCoarseSub(i);
      m_subRowBeginCoarse[i].resize(numDof+1);
      Teuchos::ArrayView<const LO> Indices;
      Teuchos::ArrayView<const SX> Values;
      LO numTerms(0);
      for (LO j=0; j<numDof; j++) {
	coarseMatrix->getLocalRowView(startRow+j, Indices, Values);
	numTerms += Indices.size();
	m_subRowBeginCoarse[i][j+1] = numTerms;
      }
      m_subValuesCoarse[i].resize(numTerms);
      m_subColumnsCoarse[i].resize(numTerms);
      startRow += numDof;
    }
    startRow = 0;
    for (LO i=0; i<numSub; i++) {
      LO numDof = getNumDofCoarseSub(i);
      Teuchos::ArrayView<const LO> Indices;
      Teuchos::ArrayView<const SX> Values;
      LO numTerms(0);
      for (LO j=0; j<numDof; j++) {
	coarseMatrix->getLocalRowView(startRow+j, Indices, Values);
	for (int k=0; k<Indices.size(); k++) {
	  LO subCol = Indices[k] - startRow;
	  assert ((subCol >= 0) && (subCol < numDof));
	  m_subValuesCoarse[i][numTerms] = Values[k];
	  m_subColumnsCoarse[i][numTerms] = subCol;
	  numTerms++;
	}
      }
      startRow += numDof;
    }
    printCoarseMatrices();
  }

  void getSubdomainCoarseDofs(LO sub, 
			      const std::vector<LO> & equivBegin,
			      std::vector<LO> & myCoarseDofs)
  {
    LO numCoarseDof = m_Subdomain[sub]->getCoarseSpaceDimension();
    myCoarseDofs.resize(numCoarseDof);
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    LO numEquiv = m_Subdomain[sub]->getNumEquiv();
    LO currentCol(0);
    for (LO j=0; j<numEquiv; j++) {
      LO equiv = subdomainEquivClasses[sub][j];
      for (LO k=equivBegin[equiv]; k<equivBegin[equiv+1]; k++) {
	myCoarseDofs[currentCol++] = k;
      }
    }
    assert (currentCol == numCoarseDof);
  }

  void getMapGIDs(RCP<const Map> map, 
		  std::vector<GO> & mapGIDs)
  {
    LO numRows = map->getNodeNumElements();
    mapGIDs.resize(numRows);
    for (LO i=0; i<numRows; i++) {
      mapGIDs[i] = map->getGlobalElement(i);
    }
  }

  void determineCoarseNodesForCoarseProc
    (const std::vector< std::vector<GO> > & coarseSubdomainSubs, 
     RCP<const CrsGraph> coarseNodesForSubdomains, 
     RCP<CrsGraph> & coarseNodesForCoarseProc)
  {
    LO numSubAll(0);
    for (size_t i=0; i<coarseSubdomainSubs.size(); i++) {
      numSubAll += coarseSubdomainSubs[i].size();
    }
    std::vector<GO> subdomainGIDs(numSubAll);
    numSubAll = 0;
    for (size_t i=0; i<coarseSubdomainSubs.size(); i++) {
      for (size_t j=0; j<coarseSubdomainSubs[i].size(); j++) {
	subdomainGIDs[numSubAll++] = coarseSubdomainSubs[i][j];
      }
    }
    RCP<const Map> subMap = generateMap(subdomainGIDs);
    coarseNodesForCoarseProc = rcp( new CrsGraph(subMap, 0) );
    Import importer(coarseNodesForSubdomains->getRowMap(), subMap);
    coarseNodesForCoarseProc->doImport(*coarseNodesForSubdomains, importer,
				       Tpetra::INSERT);
    coarseNodesForCoarseProc->fillComplete
      (coarseNodesForSubdomains->getDomainMap(), subMap);
  }

  void determineNodalDataForCoarseProblem
    (RCP<const CrsGraph> coarseNodesForCoarseProc, 
     RCP<const CrsGraph> coarseLocalDofsGraph1to1,
     RCP<MVD> coarseCoordsVec1to1)
  {
    RCP<const Map> nodeMap = coarseNodesForCoarseProc->getColMap();
    CrsGraph nodeLocalDofs(nodeMap, 0);
    Import importer(coarseLocalDofsGraph1to1->getRowMap(), nodeMap);
    nodeLocalDofs.doImport(*coarseLocalDofsGraph1to1, importer, Tpetra::INSERT);
    nodeLocalDofs.fillComplete(coarseLocalDofsGraph1to1->getDomainMap(),
			       coarseLocalDofsGraph1to1->getRangeMap());
    LO numNode = nodeMap->getNodeNumElements();
    m_numNodesCoarse = numNode;
    Teuchos::ArrayView<const LO> Indices;
    // coarse node global IDs and local dofs
    m_nodeGlobalIDsCoarse.resize(numNode);
    m_nodeBeginCoarse.resize(numNode+1, 0);
    LO numDof(0);
    for (LO i=0; i<numNode; i++) {
      m_nodeGlobalIDsCoarse[i] = nodeMap->getGlobalElement(i);
      nodeLocalDofs.getLocalRowView(i, Indices);
      numDof += Indices.size();
      m_nodeBeginCoarse[i+1] = numDof;
    }
    m_localDofsCoarse.resize(numDof);
    numDof = 0;
    RCP<const Map> colMap = nodeLocalDofs.getColMap();
    for (LO i=0; i<numNode; i++) {
      nodeLocalDofs.getLocalRowView(i, Indices);
      for (int j=0; j<Indices.size(); j++) {
	int localDof = colMap->getGlobalElement(Indices[j]);
	m_localDofsCoarse[numDof++] = localDof;
      }
    }
    // coordinates of coarse nodes
    MVD coarseCoords(nodeMap, 3);
    Import importer2(coarseCoordsVec1to1->getMap(), nodeMap);
    coarseCoords.doImport(*coarseCoordsVec1to1, importer2, Tpetra::INSERT);
    Teuchos::ArrayRCP<const SM> xvals = coarseCoords.getData(0);
    Teuchos::ArrayRCP<const SM> yvals = coarseCoords.getData(1);
    Teuchos::ArrayRCP<const SM> zvals = coarseCoords.getData(2);
    m_xCoordCoarse.resize(numNode);
    m_yCoordCoarse.resize(numNode);
    m_zCoordCoarse.resize(numNode);
    for (LO i=0; i<numNode; i++) {
      m_xCoordCoarse[i] = xvals[i];
      m_yCoordCoarse[i] = yvals[i];
      m_zCoordCoarse[i] = zvals[i];
    }
  }

  void determineNodesForEachCoarseSubdomain
    (const std::vector< std::vector<GO> > & coarseSubdomainSubs, 
     RCP<const CrsGraph> coarseNodesForCoarseProc)
  {
    LO numSubCoarse = coarseSubdomainSubs.size();
    m_subNodesCoarse.resize(numSubCoarse);
    RCP<const Map> nodeMap = coarseNodesForCoarseProc->getColMap();
    LO numNode = nodeMap->getNodeNumElements();
    std::vector<bool> nodeFlag(numNode, false);
    LO row(0);
    Teuchos::ArrayView<const LO> Indices;
    for (LO i=0; i<numSubCoarse; i++) {
      for (size_t j=0; j<coarseSubdomainSubs[i].size(); j++) {
	coarseNodesForCoarseProc->getLocalRowView(row++, Indices);
	for (int k=0; k<Indices.size(); k++) {
	  if (nodeFlag[Indices[k]] == false) {
	    m_subNodesCoarse[i].push_back(Indices[k]);
	  }
	  nodeFlag[Indices[k]] = true;
	}
      }
      for (size_t j=0; j<m_subNodesCoarse[i].size(); j++) {
	LO node = m_subNodesCoarse[i][j];
	nodeFlag[node] = false;
      }
    }
  }

  void getSubData(int row,
		  const std::vector<GO> & coarseSubdomainSubs, 
		  RCP<const CrsGraph> coarseNodesForCoarseProc,
		  std::vector<LO> & activeNodes,
		  std::vector<LO> & nodeFlag,
		  std::vector< std::vector<LO> > & localNodeNumbers,
		  std::vector<LO> & nodeBeginSub)
  {
    LO numSub = coarseSubdomainSubs.size();
    localNodeNumbers.resize(numSub);
    Teuchos::ArrayView<const LO> Indices;
    LO numActive(0);
    for (LO i=0; i<numSub; i++) {
      coarseNodesForCoarseProc->getLocalRowView(row, Indices);
      localNodeNumbers[i].resize(Indices.size());
      for (int j=0; j<Indices.size(); j++) {
	LO node = Indices[j];
	if (nodeFlag[node] == -1) {
	  activeNodes[numActive] = node;
	  nodeFlag[node] = numActive++;
	}
	localNodeNumbers[i][j] = nodeFlag[node];
      }
      row++;
    }
    nodeBeginSub.resize(numActive+1, 0);
    for (LO i=0; i<numActive; i++) {
      LO node = activeNodes[i];
      LO numDofNode = m_nodeBeginCoarse[node+1] - m_nodeBeginCoarse[node];
      nodeFlag[node] = -1;
      nodeBeginSub[i+1] = nodeBeginSub[i] + numDofNode;
    }
  }

  void determineGlobalIDsForTargetMatrix
    (std::vector<GO> & globalIDsTargetMatrix,
     GO & startRowMpiRank)
  {
    GO numRowsMpiRank(0);
    LO numSubCoarse = m_subNodesCoarse.size();
    for (LO i=0; i<numSubCoarse; i++) {
      for (size_t j=0; j<m_subNodesCoarse[i].size(); j++) {
	LO node = m_subNodesCoarse[i][j];
	LO numDofNode = m_nodeBeginCoarse[node+1] - m_nodeBeginCoarse[node];
	numRowsMpiRank += numDofNode;
      }
    }
    GO numRowsMpiRankSS;
    Teuchos::scan<int, GO> 
      (*m_Comm, Teuchos::REDUCE_SUM, 1, &numRowsMpiRank, &numRowsMpiRankSS);
    startRowMpiRank = numRowsMpiRankSS - numRowsMpiRank;
    globalIDsTargetMatrix.resize(numRowsMpiRank);
    numRowsMpiRank = 0;
    for (LO i=0; i<numSubCoarse; i++) {
      for (size_t j=0; j<m_subNodesCoarse[i].size(); j++) {
	LO node = m_subNodesCoarse[i][j];
	for (LO m=m_nodeBeginCoarse[node]; m<m_nodeBeginCoarse[node+1]; m++) {
	  globalIDsTargetMatrix[numRowsMpiRank++] = startRowMpiRank++;
	}
      }
    }
    startRowMpiRank = numRowsMpiRankSS - numRowsMpiRank;
  }
  
  void determineStartingLocationsForCoarseMatrices
    (const std::vector< std::vector<GO> > & coarseSubdomainSubs,
     RCP<const CrsGraph> coarseNodesForSubdomains, 
     RCP<const CrsGraph> coarseNodesForCoarseProc, 
     const GO startRowMpiRank,
     RCP<CrsMatrixGO> & startIndicesSub)
  {
    CrsMatrixGO coarseNodesForCoarseProcStart(coarseNodesForCoarseProc);
    LO row = 0;
    GO start = startRowMpiRank;
    std::vector<GO> startIndices;
    coarseNodesForCoarseProcStart.resumeFill();
    LO numNodeAll = coarseNodesForCoarseProc->getColMap()->getNodeNumElements();
    std::vector<LO> activeNodes(numNodeAll), nodeBeginSub(numNodeAll+1), 
      nodeFlag(numNodeAll, -1);
    Teuchos::ArrayView<const LO> Indices;
    LO numSubCoarse = m_subNodesCoarse.size();
    for (LO i=0; i<numSubCoarse; i++) {
      std::vector< std::vector<LO> > localNodeNumbers;
      getSubData(row, coarseSubdomainSubs[i], coarseNodesForCoarseProc,
		 activeNodes, nodeFlag, localNodeNumbers, nodeBeginSub);
      for (size_t j=0; j<nodeBeginSub.size()-1; j++) {
	assert (m_subNodesCoarse[i][j] == activeNodes[j]);
      }
      for (size_t j=0; j<coarseSubdomainSubs[i].size(); j++) {
	const std::vector<LO> & localNodes = localNodeNumbers[j];
	startIndices.resize(localNodes.size());
	for (size_t k=0; k<localNodes.size(); k++) {
	  startIndices[k] = start + nodeBeginSub[localNodes[k]];
	}
	coarseNodesForCoarseProc->getLocalRowView(row, Indices);
	coarseNodesForCoarseProcStart.replaceLocalValues
	  (row++, Indices, Teuchos::ArrayView<const GO>(startIndices)); 
      }
      start += nodeBeginSub.back();
    }
    coarseNodesForCoarseProcStart.fillComplete();
    startIndicesSub = rcp( new CrsMatrixGO(coarseNodesForSubdomains) );
    startIndicesSub->resumeFill();
    Import importer(coarseNodesForCoarseProcStart.getRowMap(),
		    coarseNodesForSubdomains->getRowMap());
    startIndicesSub->doImport(coarseNodesForCoarseProcStart, importer,
			      Tpetra::REPLACE);
    startIndicesSub->fillComplete();
  }

  void determineStartIndicesNodes
    (RCP<const CrsMatrixGO> startIndicesSubNodes,
     std::vector< std::vector<GO> > & startIndicesNodes)
  {
    LO numSub  = startIndicesSubNodes->getRowMap()->getNodeNumElements();
    LO numNode = startIndicesSubNodes->getColMap()->getNodeNumElements();
    assert (numSub == m_numSub);
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const GO> Values;
    std::vector<LO> count(numNode, 0);
    for (LO i=0; i<numSub; i++) {
      startIndicesSubNodes->getLocalRowView(i, Indices, Values);
      for (int j=0; j<Indices.size(); j++) {
	count[Indices[j]]++;
      }
    }
    startIndicesNodes.resize(numNode);
    for (LO i=0; i<numNode; i++) {
      startIndicesNodes[i].resize(count[i]);
      count[i] = 0;
    }
    for (LO i=0; i<numSub; i++) {
      startIndicesSubNodes->getLocalRowView(i, Indices, Values);
      for (int j=0; j<Indices.size(); j++) {
	LO index = Indices[j];
	startIndicesNodes[index][count[index]++] = Values[j];
      }
    }
  }

  void determineSourceMatrixBookkeepingData
    (const std::vector< std::vector<GO> > & startIndicesNodes, 
     const std::vector<LO> & coarseNodesNumDofs,
     std::vector<LO> & nodeStart,
     std::vector<LO> & nodeBeginAll,
     std::vector<GO> & globalIDsSourceMatrix,
     LO & maxNumDofsNode)
  {
    LO numNodeAll(0), numDofAll(0);
    maxNumDofsNode = 0;
    LO numNode = startIndicesNodes.size();
    std::vector<LO> count(numNode);
    for (LO i=0; i<numNode; i++) {
      std::vector<GO> uniqueIndices;
      DofManager<LO,GO>::determineUniqueIndices(startIndicesNodes[i], 
						uniqueIndices);
      count[i] = uniqueIndices.size();
      numNodeAll += count[i];
      numDofAll += count[i]*coarseNodesNumDofs[i];
      if (coarseNodesNumDofs[i] > maxNumDofsNode) {
	maxNumDofsNode = coarseNodesNumDofs[i];
      }
    }
    nodeBeginAll.resize(numNodeAll+1, 0);
    nodeStart.resize(numNode);
    globalIDsSourceMatrix.resize(numDofAll);
    numNodeAll = numDofAll = 0;
    for (LO i=0; i<numNode; i++) {
      nodeStart[i] = numNodeAll;
      for (LO j=0; j<count[i]; j++) {
	for (LO k=0; k<coarseNodesNumDofs[i]; k++) {
	  globalIDsSourceMatrix[numDofAll++] = startIndicesNodes[i][j] + k;
	}
	nodeBeginAll[numNodeAll+1] = numDofAll;
	numNodeAll++;
      }
    }
  }

  void determineSubdomainNodesAccountingForAnyNodeReplications
    (RCP<const CrsMatrixGO> startIndicesSubNodes, 
     const std::vector< std::vector<GO> > & startIndicesNodes,
     const std::vector<LO> & nodeStart,
     std::vector< std::vector<LO> > & subNodes)
  {
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const GO> Values;
    subNodes.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      startIndicesSubNodes->getLocalRowView(i, Indices, Values);
      subNodes[i].resize(Indices.size());
      for (int j=0; j<Indices.size(); j++) {
	LO node = Indices[j];
	std::vector<GO> uniqueIndices;
	DofManager<LO,GO>::determineUniqueIndices(startIndicesNodes[node], 
						  uniqueIndices);
	auto lower = std::find(uniqueIndices.begin(), uniqueIndices.end(), 
			       Values[j]);
	assert (lower != uniqueIndices.end());
	int offset = std::distance(uniqueIndices.begin(), lower);
	LO nodeAll = nodeStart[node] + offset;
	subNodes[i][j] = nodeAll;
      }
    }
  }

  void determineNodeSubdomains
    (const std::vector< std::vector<LO> > & subNodes, 
     const std::vector<LO> & nodeBeginAll,
     std::vector< std::vector<LO> > & nodeSubs,
     std::vector< std::vector<LO> > & nodeSubsStart)
  {
    LO numNodeAll = nodeBeginAll.size() - 1;
    std::vector<LO> count(numNodeAll, 0);
    for (LO i=0; i<m_numSub; i++) {
      for (size_t j=0; j<subNodes[i].size(); j++) {
	count[subNodes[i][j]]++;
      }
    }
    nodeSubs.resize(numNodeAll);
    nodeSubsStart.resize(numNodeAll);
    for (LO i=0; i<numNodeAll; i++) {
      nodeSubs[i].resize(count[i]);
      nodeSubsStart[i].resize(count[i]);
      count[i] = 0;
    }
    for (LO i=0; i<m_numSub; i++) {
      LO start = 0;
      for (size_t j=0; j<subNodes[i].size(); j++) {
	LO node = subNodes[i][j];
	nodeSubs[node][count[node]] = i;
	nodeSubsStart[node][count[node]++] = start;
	start += nodeBeginAll[node+1] - nodeBeginAll[node];
      }
    }
  }

  void determineNodalConnectivity
    (const std::vector< std::vector<LO> > & subNodes, 
     const std::vector< std::vector<LO> > & nodeSubs, 
     std::vector< std::vector<LO> > & nodalConn)
  {
    // determine nodal connectivity 
    LO numNodeAll = nodeSubs.size();
    nodalConn.resize(numNodeAll);
    std::vector<LO> nodeFlag(numNodeAll, -1), connNodes(numNodeAll);
    for (LO i=0; i<numNodeAll; i++) {
      LO numConn = 0;
      for (size_t j=0; j<nodeSubs[i].size(); j++) {
	LO sub = nodeSubs[i][j];
	for (size_t k=0; k<subNodes[sub].size(); k++) {
	  LO node = subNodes[sub][k];
	  if (nodeFlag[node] == -1) {
	    nodeFlag[node] = numConn;
	    connNodes[numConn++] = node;
	  }
	}
      }
      nodalConn[i].resize(numConn);
      for (LO j=0; j<numConn; j++) {
	nodalConn[i][j]= connNodes[j];
	nodeFlag[connNodes[j]] = -1;
      }
    }
  }

  void initializeSourceMatrix
    (const std::vector<GO> & globalIDsSourceMatrix, 
     const std::vector< std::vector<LO> > & nodalConn,
     const std::vector<LO> & coarseNodesNumDofs,
     const std::vector<LO> & nodeBeginAll,
     std::vector<LO> & rowCount,
     RCP<CrsMatrix> & sourceMatrix)
  {
    LO numNodeAll = nodalConn.size();
    LO numDofAll = globalIDsSourceMatrix.size();
    Teuchos::ArrayRCP<size_t> count(numDofAll);
    rowCount.resize(numNodeAll);
    for (LO i=0; i<numNodeAll; i++) {
      LO numTermsRow(0);
      for (size_t j=0; j<nodalConn[i].size(); j++) {
	LO node = nodalConn[i][j];
	numTermsRow += coarseNodesNumDofs[node];
      }
      rowCount[i] = numTermsRow;
      for (LO j=nodeBeginAll[i]; j<nodeBeginAll[i+1]; j++) {
	count[j] = numTermsRow;
      }
    }
    RCP<const Map> rowMap = generateMap(globalIDsSourceMatrix);
    sourceMatrix = 
      rcp( new CrsMatrix(rowMap, rowMap, count, Tpetra::StaticProfile) );
  }

  void getMyIndices(const std::vector<LO> & subNodes, 
		    const std::vector<LO> & nodeBeginAll, 
		    LO & startIndex,
		    std::vector<LO> & nodeFlag, 
		    std::vector<LO> & myIndices)
  {
    myIndices.resize(0);
    for (size_t k=0; k<subNodes.size(); k++) {
      LO node = subNodes[k];
      if (nodeFlag[node] == -1) {
	nodeFlag[node] = startIndex;
	startIndex += nodeBeginAll[node+1] - nodeBeginAll[node];
      }
      for (LO m=nodeBeginAll[node]; m<nodeBeginAll[node+1]; m++) {
	LO mm = m - nodeBeginAll[node];
	myIndices.push_back(nodeFlag[node]+mm);
      }
    }
  }
  
  void assembleSourceMatrix
    (const std::vector<LO> & nodeBeginAll, 
     const std::vector<LO> & rowCount,
     const LO maxNumDofsNode,
     const std::vector< std::vector<LO> > & nodalConn, 
     const std::vector< std::vector<LO> > & subNodes, 
     const std::vector< std::vector<LO> > & nodeSubs, 
     const std::vector< std::vector<LO> > & nodeSubsStart,
     RCP<CrsMatrix> sourceMatrix)
  {
    // calculate coarse subdomain matrices (this could be threaded later)
    std::vector< std::vector<SX> > subAc(m_numSub);
    bool restrictPhiToBoundary(true), scaleRows(true);  
    for (LO i=0; i<m_numSub; i++) {
      m_Subdomain[i]->calculateCoarseMatrices(subAc[i], restrictPhiToBoundary,
					      scaleRows);
    }
    // assemble source coarse matrix
    LO numNodeAll = nodeBeginAll.size() - 1;
    std::vector<LO> colIndices, myIndices, nodeFlag(numNodeAll, -1);
    std::vector< std::vector<SX> > rowValues(maxNumDofsNode);
    for (LO i=0; i<numNodeAll; i++) {
      colIndices.resize(rowCount[i]);
      for (LO j=nodeBeginAll[i]; j<nodeBeginAll[i+1]; j++) {
	LO jj = j - nodeBeginAll[i];
	rowValues[jj].assign(rowCount[i], 0);
      }
      LO numEntries(0);
      for (size_t j=0; j<nodalConn[i].size(); j++) {
	LO node = nodalConn[i][j];
	for (LO k=nodeBeginAll[node]; k<nodeBeginAll[node+1]; k++) {
	  colIndices[numEntries++] = k;
	}
      }
      LO startIndex = 0;
      for (size_t j=0; j<nodeSubs[i].size(); j++) {
	LO sub = nodeSubs[i][j];
	getMyIndices(subNodes[sub], nodeBeginAll, startIndex, nodeFlag, 
		     myIndices);
	LO numDofSub = myIndices.size();
	assert (LO(subAc[sub].size()) == numDofSub*numDofSub);
	LO subStart = nodeSubsStart[i][j];
	for (LO j=nodeBeginAll[i]; j<nodeBeginAll[i+1]; j++) {
	  LO jj = j - nodeBeginAll[i];
	  for (LO k=0; k<numDofSub; k++) {
	    rowValues[jj][myIndices[k]] += subAc[sub][subStart+k*numDofSub];
	  }
	  subStart++;
	}
      }
      for (size_t j=0; j<nodalConn[i].size(); j++) {
	LO node = nodalConn[i][j];
	nodeFlag[node] = -1;
      }
      for (LO j=nodeBeginAll[i]; j<nodeBeginAll[i+1]; j++) {
	LO jj = j - nodeBeginAll[i];
	sourceMatrix->insertLocalValues(j, Teuchos::ArrayView<LO>(colIndices),
					Teuchos::ArrayView<SX>(rowValues[jj]));
      }
    }
    sourceMatrix->fillComplete();
    printSourceMatrix(sourceMatrix);
  }

  void printSourceMatrix(RCP<CrsMatrix> sourceMatrix)
  {
    char fileName[101];
    sprintf(fileName, "sourceAc_%d.dat", m_myPID);
    std::ofstream fout;
    fout.open(fileName);
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const SX> Values;
    for (size_t i=0; i<sourceMatrix->getNodeNumRows(); i++) {
      sourceMatrix->getLocalRowView(i, Indices, Values);
      for (int j=0; j<Indices.size(); j++) {
	fout << i+1 << "  " << Indices[j]+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	fout << Values[j] << std::endl;
      }
    }
    fout.close();
  }

  void assembleTargetMatrix(RCP<const CrsMatrix> sourceMatrix, 
			    const std::vector<GO> & globalIDsTargetMatrix,
			    RCP<CrsMatrix> & targetMatrix)
  {
    RCP<const Map> sourceMap = sourceMatrix->getRowMap();
    RCP<const Map> targetMap = generateMap(globalIDsTargetMatrix);
    Export exporter(sourceMap, targetMap);
    // could do some local work here to determine sparsity pattern of the
    // targetMatrix, but let's just be lazy for now
    targetMatrix = rcp( new CrsMatrix(targetMap, targetMap, 0) );
    targetMatrix->doExport(*sourceMatrix, exporter, Tpetra::ADD);
    targetMatrix->fillComplete();
  }

  void determineSubdomainDataForCoarseProblem
    (const std::vector< std::vector<GO> > & coarseSubdomainSubs, 
     RCP<const CrsGraph> coarseNodesForSubdomains, 
     RCP<const CrsGraph> coarseNodesForCoarseProc,
     const std::vector<LO> & coarseNodesNumDofs)
  {
    determineNodesForEachCoarseSubdomain
      (coarseSubdomainSubs, coarseNodesForCoarseProc);

    std::vector<GO> globalIDsTargetMatrix;
    GO startRowMpiRank;
    determineGlobalIDsForTargetMatrix(globalIDsTargetMatrix, startRowMpiRank);

    RCP<CrsMatrixGO> startIndicesSubNodes;
    determineStartingLocationsForCoarseMatrices
      (coarseSubdomainSubs, coarseNodesForSubdomains, coarseNodesForCoarseProc,
       startRowMpiRank, startIndicesSubNodes);

    std::vector< std::vector<GO> > startIndicesNodes;
    determineStartIndicesNodes(startIndicesSubNodes, startIndicesNodes);

    // since a coarse node may appear in more than a single coarse subdomain,
    // we introduce additional "nodes" to account for this (reflected by the
    // introduction of nodeBeginAll and nodeStart vectors immediately below)
    std::vector<LO> nodeBeginAll, nodeStart;
    std::vector<GO> globalIDsSourceMatrix;
    LO maxNumDofsNode;
    determineSourceMatrixBookkeepingData
      (startIndicesNodes, coarseNodesNumDofs, nodeStart, nodeBeginAll,
       globalIDsSourceMatrix, maxNumDofsNode);

    std::vector< std::vector<LO> > subNodes;
    determineSubdomainNodesAccountingForAnyNodeReplications
      (startIndicesSubNodes, startIndicesNodes, nodeStart, subNodes);

    std::vector< std::vector<LO> > nodeSubs, nodeSubsStart;
    determineNodeSubdomains(subNodes, nodeBeginAll, nodeSubs, nodeSubsStart);

    std::vector< std::vector<LO> > nodalConn;
    determineNodalConnectivity(subNodes, nodeSubs, nodalConn);

    std::vector<LO> rowCount;
    RCP<CrsMatrix> sourceMatrix;
    initializeSourceMatrix
      (globalIDsSourceMatrix, nodalConn, coarseNodesNumDofs, nodeBeginAll, 
       rowCount, sourceMatrix);

    assembleSourceMatrix(nodeBeginAll, rowCount, maxNumDofsNode, nodalConn, 
			 subNodes, nodeSubs, nodeSubsStart, sourceMatrix);

    RCP<CrsMatrix> targetMatrix;
    assembleTargetMatrix(sourceMatrix, globalIDsTargetMatrix, targetMatrix);
    extractSubdomainMatrices(targetMatrix);
  }
  
  void determineProcessorTypes(bool & fineProc, 
			       bool & coarseProc)
  {
    fineProc = coarseProc = false;
    if (m_numSub > 0) {
      fineProc = true;
    }
    else {
      coarseProc = true;
    }
    if (m_myCoarseMpiRank != -1) coarseProc = true;
  }

  void splitCommunicator()
  {
    int color(1);
    if (m_coarseProc == true) color = 0;
    int key = 0;
    MPI_Comm_split(m_mpiComm, color, key, &m_mpiCommSplit);
    if (m_disjointCommunicators == true) {
      if (m_coarseProc == false) {
	m_CommFine = rcp( new Teuchos::MpiComm<int>(m_mpiCommSplit) );
      }
    }
  }

  bool determineIfCommunicatorsAreDisjoint()
  {
    GO fineProc(0), coarseProc(0), fineProcSum(0), coarseProcSum(0);
    if (m_fineProc == true) fineProc = 1;
    if (m_coarseProc == true) coarseProc = 1;
    Teuchos::reduceAll<int, GO> 
      (*m_Comm, Teuchos::REDUCE_SUM, 1, &fineProc, &fineProcSum);
    Teuchos::reduceAll<int, GO> 
      (*m_Comm, Teuchos::REDUCE_SUM, 1, &coarseProc, &coarseProcSum);
    GO numProc = m_Comm->getSize();
    bool disjointCommunicator = false;
    if ((fineProcSum+coarseProcSum) == numProc) {
      disjointCommunicator = true;
    }
    return disjointCommunicator;
  }

  void determineCoarseSpace()
  {
    // mpiRanks is a vector of active MPI ranks at the coarse level 
    //   (this is the same vector across all current level MPI ranks)
    // m_myCoarseMpiRank is the MPI rank at the coarse level (-1 if inactive)
    std::vector<int> mpiRanks;
    determineCoarseMpiRanks(m_myCoarseMpiRank, mpiRanks);
    // m_fineProc   = true if used to store subdomain data at current level
    // m_coarseProc = true if used to store subdomain data at coarse level
    determineProcessorTypes(m_fineProc, m_coarseProc);
    m_disjointCommunicators = determineIfCommunicatorsAreDisjoint();
    // partsAcrossMpiRanks[i] is the destination coarse MPI rank for local
    // subdomain i
    std::vector<GO> partsAcrossMpiRanks;
    partitionSubdomainsAcrossMpiRanks(mpiRanks, partsAcrossMpiRanks);
    // subsForCoarseProc is a graph with at most a single local row. The
    //  columns of subsForCoarseProc are globalIDs of subdomains for this
    //  coarse processor.
    RCP<CrsGraph> subsForCoarseProc;
    determineSubdomainsInCoarsePartition
      (m_myCoarseMpiRank, partsAcrossMpiRanks, subsForCoarseProc);
    // coarseSubdomainSubs[i] is a vector of subdomain globalIDs for part
    //   i of a decomposition of the subdomains in subsForCoarseProc.
    std::vector< std::vector<GO> > coarseSubdomainSubs;
    partitionSubdomainsWithinMpiRanks(subsForCoarseProc, coarseSubdomainSubs);
    // activeCoarseNodes is a list of active coarse node indices
    // coarseCoordsVec1to1 is a multivector containing active coarse node 
    //   coordinates
    // coarseLocalDofsGraph1to1 is a CrsGraph containing active coarse node
    //  local dof numbers
    // nodeExporter is an Export object with on-processor active coarse node
    // globalIDs as the sourceMap and the corresponding 1to1 targetMap
    std::vector<LO> activeCoarseNodes, activeCoarseNodesNumDofs;
    RCP<MVD> coarseCoordsVec1to1;
    RCP<CrsGraph> coarseLocalDofsGraph1to1;
    RCP<Export> nodeExporter;
    determineActiveCoarseNodeCoordinatesAndLocalDofs
      (activeCoarseNodes, activeCoarseNodesNumDofs,
       coarseCoordsVec1to1, coarseLocalDofsGraph1to1, nodeExporter);
    // startNodeVector1to1 contains the starting locations of all the 
    //  active coarse nodes
    // numDofNodeVector1to1 contains the number of dofs for all of the
    //  active coarse nodes
    // m_coarseMap is a Map object for the coarse dofs
    RCP<VectorGO> startNodeVector1to1, numDofNodeVector1to1;
    determineCoarseMap
      (coarseLocalDofsGraph1to1, nodeExporter, activeCoarseNodesNumDofs,
       startNodeVector1to1, numDofNodeVector1to1, m_coarseMap);
    // coarseNodesForSubdomains contains the active coarse nodes for
    //   each coarse subdomain
    RCP<CrsGraph> coarseNodesForSubdomains;
    determineCoarseNodesForSubdomainsAndSubdomainCoarseDofs
      (activeCoarseNodes, activeCoarseNodesNumDofs, coarseNodesForSubdomains);
    // coarseNodesForCoarseProc contains the coarse nodes for the subdomains
    // assigned to a coarse processor
    RCP<CrsGraph> coarseNodesForCoarseProc;
    determineCoarseNodesForCoarseProc
      (coarseSubdomainSubs, coarseNodesForSubdomains, coarseNodesForCoarseProc);

    determineNodalDataForCoarseProblem
      (coarseNodesForCoarseProc, coarseLocalDofsGraph1to1, 
       coarseCoordsVec1to1); 
    // m_coarseMapNextLevel is the counterpart of m_coarseMap, but at the
    //  next coarser level
    determineCoarseMapNextLevel(startNodeVector1to1, numDofNodeVector1to1,
				m_coarseMapNextLevel);

    determineSubdomainDataForCoarseProblem
      (coarseSubdomainSubs, coarseNodesForSubdomains, coarseNodesForCoarseProc,
       activeCoarseNodesNumDofs);

    allocateCoarseEntities();
    splitCommunicator();
    if (m_coarseProc == false) return;
    std::vector< LO* > subRowBeginPtr, subColumnsPtr;
    std::vector< SX* > subValuesPtr;
    convertToPointers
      (m_subRowBeginCoarse, m_subColumnsCoarse, m_subValuesCoarse,
       subRowBeginPtr, subColumnsPtr, subValuesPtr);
    m_coarsePreconditioner = 
      rcp( new PreconditionerBDDC<SX,SM,LO,GO>
	(m_numNodesCoarse, m_nodeBeginCoarse.data(), m_localDofsCoarse.data(), 
	 m_nodeGlobalIDsCoarse.data(), m_xCoordCoarse.data(), 
	 m_yCoordCoarse.data(), m_zCoordCoarse.data(), m_subNodesCoarse, 
	 subRowBeginPtr.data(), subColumnsPtr.data(), subValuesPtr.data(), 
	 m_Parameters, m_mpiCommSplit, m_level+1, &m_nodeGlobalIDsCoarse1to1) );
  }

  void convertToPointers(std::vector< std::vector<LO> > & subRowBeginCoarse, 
			 std::vector< std::vector<LO> > & subColumnsCoarse, 
			 std::vector< std::vector<SX> > & subValuesCoarse,
			 std::vector< LO* > & subRowBeginPtr, 
			 std::vector< LO* > & subColumnsPtr,
			 std::vector< SX* > & subValuesPtr)
  {    
    const LO numSub = subRowBeginCoarse.size();
    subRowBeginPtr.resize(numSub);
    subColumnsPtr.resize(numSub);
    subValuesPtr.resize(numSub);
    for (LO i=0; i<numSub; i++) {
      subRowBeginPtr[i] = subRowBeginCoarse[i].data();
      subColumnsPtr[i] = subColumnsCoarse[i].data();
      subValuesPtr[i] = subValuesCoarse[i].data();
    }
  }

  void initializeSubdomains()
  {
    m_Subdomain.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      m_Subdomain[i] = new bddc::SubdomainBDDC<SX,SM,LO,GO>
	(m_subNodes[i].size(), &m_subNodeBegin[i][0], &m_subLocalDofs[i][0],
	 &m_subRowBegin[i][0], &m_subColumns[i][0], &m_subValues[i][0], 
	 m_boundaryDofsLocal[i].size(), &m_boundaryDofsLocal[i][0], 
	 *m_Parameters);
    }
  }

  void determineBaseConstraints()
  {
    m_Constraints = 
      rcp( new bddc::ConstraintsBDDC<SX,SM,LO,GO>
	   (m_numNodes, m_nodeBegin, m_localDofs, m_xCoord, m_yCoord, m_zCoord,
	    m_subNodes, m_Subdomain, m_Partition, 
	    m_exporterB, m_diagBoundary, m_Parameters) );
    m_Constraints->determineBaseConstraints();
  }

  void determineAuxiliaryConstraints()
  {
    // nothing for now
  }

  void determineDiagBoundary()
  {
    m_diagBoundary.resize(m_numDofsB);
    for (LO i=0; i<m_numSub; i++) {
      LO numB = m_Subdomain[i]->getNumBoundaryDofs();
      for (LO j=0; j<numB; j++) {
	LO dofB = m_subBoundaryDofs[i][j];
	m_diagBoundary[dofB] += m_Subdomain[i]->getBoundaryDiagValue(j);
      }
    }
    VectorSM diagBoundary(m_dofMapB), diagBoundary1to1(m_dofMapB1to1);
    for (LO i=0; i<m_numDofsB; i++) {
      diagBoundary.replaceLocalValue(i, m_diagBoundary[i]);
    }
    diagBoundary1to1.doExport(diagBoundary, *m_exporterB, Tpetra::ADD);
    diagBoundary.doImport(diagBoundary1to1, *m_exporterB, Tpetra::INSERT);
    Teuchos::ArrayRCP<const SM> values = diagBoundary.getData();
    for (LO i=0; i<m_numDofsB; i++) m_diagBoundary[i] = values[i];
  }

  void getSubDofs(LO sub, 
		  std::vector<LO> & subDofs)
  {
    const std::vector<LO> & subNodes = m_subNodes[sub];
    LO numNode = subNodes.size();
    for (LO i=0; i<numNode; i++) {
      LO nodeGlobal = subNodes[i];
      for (LO k=m_nodeBegin[nodeGlobal]; k<m_nodeBegin[nodeGlobal+1]; k++) {
	subDofs.push_back(k);
      }
    }
  }

  void getBoundaryDofs(LO sub,
		       const std::vector<LO> & subEquivClasses,
		       const std::vector< std::vector<LO> > & equivClasses,
		       std::vector<LO> & globalToLocalMap,
		       std::vector<LO> & boundaryDofsLocal,
		       std::vector<LO> & boundaryDofsGlobal)
  {
    const std::vector<LO> & subNodes = m_subNodes[sub];
    const std::vector<LO> & nodeBeginSub = m_subNodeBegin[sub];
    LO numEquiv = subEquivClasses.size();
    LO numNode = subNodes.size();
    for (LO i=0; i<numNode; i++) globalToLocalMap[subNodes[i]] = i;
    for (LO i=0; i<numEquiv; i++) {
      LO equiv = subEquivClasses[i];
      for (size_t j=0; j<equivClasses[equiv].size(); j++) {
	LO nodeGlobal = equivClasses[equiv][j];
	LO nodeLocal = globalToLocalMap[nodeGlobal];
	assert (nodeLocal != -1);
	for (LO k=nodeBeginSub[nodeLocal]; k<nodeBeginSub[nodeLocal+1]; k++) {
	  boundaryDofsLocal.push_back(k);
	}
	for (LO k=m_nodeBegin[nodeGlobal]; k<m_nodeBegin[nodeGlobal+1]; k++) {
	  boundaryDofsGlobal.push_back(k);
	}
	LO numDofsLocal = nodeBeginSub[nodeLocal+1] - nodeBeginSub[nodeLocal];
	LO numDofsGlobal = m_nodeBegin[nodeGlobal+1] - m_nodeBegin[nodeGlobal];
	assert (numDofsLocal == numDofsGlobal);
      }
    }
  }

};

} // namespace bddc

#endif // PRECONDITIONERBDDC_H
  
