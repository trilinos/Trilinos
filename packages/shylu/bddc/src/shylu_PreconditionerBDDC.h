
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

#include "PartitionOfUnityBDDC.h"
#include "SubdomainBDDC.h"
#include "WeightsBDDC.h"
#include "ConstraintsBDDC.h"
#include "CoarseSpaceBDDC.h"
#include "DofManager.h"
#include "UtilBDDC.h"

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {
  
enum Timings{
  TIME_INIT_STATIC_COND, TIME_STATIC_EXPANSION, TIME_SUB_CORR, 
  TIME_COARSE_CORR, TIME_APPLY_OPER, TIME_INITIALIZATION, TIME_LENGTH
};

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
  typedef Tpetra::CrsMatrix<GO,LO,GO,Node>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO,Node>                              Export;
  typedef Tpetra::Import<LO,GO,Node>                              Import;
  typedef Tpetra::Vector<SM,LO,GO,Node>                           VectorSM;
  typedef Tpetra::Vector<SX,LO,GO,Node>                           Vector;
  typedef Tpetra::Vector<LO,LO,GO,Node>                           VectorLO;
  typedef Tpetra::MultiVector<SX,LO,GO,Node>                      MV;

  PreconditionerBDDC()
  {
  }
  PreconditionerBDDC
    (LO numNodes,
     const LO* nodeBegin,
     const LO* localDofs,
     const GO* nodeGlobalIDs,
     const SM* xCoord,
     const SM* yCoord,
     const SM* zCoord,
     const std::vector< std::vector<LO> > & subNodeBegin,
     const std::vector< std::vector<LO> > & subLocalDofs,
     const std::vector< std::vector<LO> > & subNodes,
     const std::vector< std::vector<LO> > & subRowBegin,
     const std::vector< std::vector<LO> > & subColumns,
     const std::vector< std::vector<SX> > & subValues,
     RCP<Teuchos::ParameterList> & Parameters,
     MPI_Comm mpiComm) :
    m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_subNodeBegin(subNodeBegin),
    m_subLocalDofs(subLocalDofs),
    m_subNodes(subNodes),
    m_subRowBegin(subRowBegin),
    m_subColumns(subColumns),
    m_subValues(subValues),
    m_Parameters(Parameters),
    m_mpiComm(mpiComm),
    m_spatialDim(Parameters->get("Spatial Dimension", 3)),
    m_problemType(Parameters->get("Problem Type BDDC", bddc::SCALARPDE)),
    m_Comm(rcp( new Teuchos::MpiComm<int>(mpiComm) )),
    m_myPID(m_Comm->getRank()),
    m_numThreads(Parameters->get<int>("numThreadsOuter", 1)),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()),
    m_interfacePreconditioner
    (Parameters->get("Interface Preconditioner", false)),
    m_addCorners(Parameters->get("Add Corners", true)),
    m_numDofs(nodeBegin[numNodes]),
    m_numSub(0),
    m_numDofsB(0),
    m_numDofsB1to1(0),
    m_work1T(0),
    m_work2T(0)
  {
    const LO* nodeOnBoundary = 
      Parameters->get("Node On Boundary Array", (LO*)(0));
    double startTime = GetTime();
    m_Partition =
      rcp( new bddc::PartitionOfUnity<SX,SM,LO,GO>
	   (numNodes, nodeGlobalIDs, subNodeBegin, subNodes, 
	    subRowBegin, subColumns, subValues, nodeOnBoundary,
	    m_spatialDim, m_addCorners, m_Comm) );
    bddc::DofManager<LO,GO>::
      determineGlobalIDs(numNodes, nodeGlobalIDs, nodeBegin, localDofs,
			 m_Comm, m_dofMap, m_dofMap1to1);
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
    reserveMemory();
    m_timings[TIME_INITIALIZATION] += GetTime() - startTime;
  }
   
  ~PreconditionerBDDC()
  {
    for (LO i=0; i<m_numSub; i++) {
      delete m_Subdomain[i];
    }
    for (int i=0; i<m_numThreads; i++) {
      delete [] m_work1T[i];
      delete [] m_work2T[i];
    }
    delete [] m_work1T;
    delete [] m_work2T;
  }

  const std::vector<double> & getTimings() const
  {
    return m_timings;
  }

  LO getCoarseSpaceDimension() 
  {
    return m_CoarseSpace->getCoarseSpaceDimension();
  }

  LO getNumDofsB() const
  {
    return m_numDofsB;
  }

  LO getNumDofsB1to1() const
  {
    return m_numDofsB1to1;
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
	LO numDofs = m_Subdomain[i]->getNumDofs();
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
#pragma omp critical(initial_static_condensation)
	{
	  for (LO j=0; j<numBoundaryDofs; j++) {
	    LO row = boundaryDofs[j];
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
      m_solVecB->doImport(*m_solVecB1to1, *m_importerB, Tpetra::INSERT);
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

  void GetInterfaceValues(SX* allValues,
			  SX* interfaceValues)
  {
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

  GO NumGlobalRowsKrylov()
  {
    if (m_interfacePreconditioner == true) {
      return m_dofMapB1to1->getGlobalNumElements();
    }
    else {
      return m_dofMap1to1->getGlobalNumElements();
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
    // subdomain correction
    double startTime = GetTime();
    Teuchos::ArrayRCP<SX> subRhs1to1 = m_rhsVecB1to1->getDataNonConst();
    for (LO i=0; i<m_numDofsB1to1; i++) subRhs1to1[i] = r[i];
    m_rhsVecB->doImport(*m_rhsVecB1to1, *m_importerB, Tpetra::INSERT);
    Teuchos::ArrayRCP<SX> rhsB = m_rhsVecB->getDataNonConst();
    Teuchos::ArrayRCP<SX> solB = m_solVecB->getDataNonConst();
    for (LO i=0; i<m_numDofsB; i++) solB[i] = 0;
#pragma omp parallel num_threads(m_numThreads)
    {
#pragma omp for
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & subBoundaryDofs = m_subBoundaryDofs[i];
	LO numDofB = subBoundaryDofs.size();
	SX *subRhs(0), *subSol(0);
	int threadID = getThreadID();
	getArrays(threadID, subRhs, subSol);
	for (LO j=0; j<numDofB; j++) {
	  subRhs[j] = rhsB[subBoundaryDofs[j]];
	}
	m_Subdomain[i]->applyNeumannCorrection(&subRhs[0], &subSol[0]);
#pragma omp critical(apply)
	{
	  for (LO j=0; j<numDofB; j++) {
	    solB[subBoundaryDofs[j]] += subSol[j];
	  }
	}
      }
    }
    m_solVecB1to1->doExport(*m_solVecB, *m_exporterB, Tpetra::ADD);
    m_timings[TIME_SUB_CORR] += GetTime() - startTime;
    // coarse correction
    startTime = GetTime();
    Teuchos::ArrayRCP<SX> coarseRhs1to1 = m_coarseRhsVec1to1->getDataNonConst();
    for (LO i=0; i<m_numDofsB1to1; i++) coarseRhs1to1[i] = r[i];
    m_CoarseSpace->applyCoarseCorrection(m_coarseRhsVec1to1, 
					 m_coarseSolVec1to1);
    m_timings[TIME_COARSE_CORR] += GetTime() - startTime;
    // sum of subdomain and coarse corrections
    Teuchos::ArrayRCP<SX> subSol1to1 = m_solVecB1to1->getDataNonConst();
    Teuchos::ArrayRCP<SX> coarseSol1to1 = m_coarseSolVec1to1->getDataNonConst();
    for (LO i=0; i<m_numDofsB1to1; i++) {
      Pr[i] = subSol1to1[i] + coarseSol1to1[i];
    }
    startTime = GetTime();
    ApplyOperator(Pr, APr);
    m_timings[TIME_APPLY_OPER] += GetTime() - startTime;
  }

  LO numBoundaryDofs() {
  }

  void ApplyOperator(SX* x, 
		     SX* Ax)
  {
    if (m_interfacePreconditioner) {
      RCP<Vector> & xVec1to1 = m_rhsVecB1to1;
      Teuchos::ArrayRCP<SX> xValues1to1 = xVec1to1->getDataNonConst();
      for (LO i=0; i<m_numDofsB1to1; i++) xValues1to1[i] = x[i];
      RCP<Vector> & xVec = m_rhsVecB;
      xVec->doImport(*xVec1to1, *m_importerB, Tpetra::INSERT);
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
#pragma omp critical(apply_operator)
	  {
	    for (LO j=0; j<numDofB; j++) {
	      AxValues[subBoundaryDofs[j]] += AxSub[j];
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
    m_xVecAll->doImport(*m_xVecAll1to1, *m_importerAll, Tpetra::INSERT);
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
#pragma omp critical(apply_full_operator)
	{
	  for (LO j=0; j<numDofs; j++) {
	    AxValues[subDofs[j]] += AxSub[j];
	  }
	}
      }
    }
    m_AxVecAll1to1->doExport(*m_AxVecAll, *m_exporterAll, Tpetra::ADD);
    Teuchos::ArrayRCP<SX> AxValuesAll1to1 = m_AxVecAll1to1->getDataNonConst();
    for (LO i=0; i<numDofs1to1; i++) Ax[i] = AxValuesAll1to1[i];
  }

  void MassMatrixMultiply(SX* x, 
			  SX* Mx)
  {
  }

  RCP< Teuchos::ParameterList > GetParameters()
  {
      return m_Parameters;
  }

  bool InterfacePreconditioner() {
    return m_interfacePreconditioner;
  }

  void MaxAll(SM* localMax,
	      SM* globalMax,
	      int numEntries)
  {
    Teuchos::reduceAll<int, SX> (*m_Comm, Teuchos::REDUCE_MAX, numEntries, 
				 localMax, globalMax);
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

  virtual MPI_Comm GetMPI_Comm() {
    return m_mpiComm;
  }

  double GetTime()
  {
    return MPI_Wtime();
  }

  void OutputPreconditioner()
  {
  }

  void MatrixHasChanged() 
  {
  }

  RCP< bddc::WeightsBDDC<SX,SM,LO,GO> > getWeights()
  {
    return m_Weights;
  }

  RCP< bddc::CoarseSpaceBDDC<SX,SM,LO,GO> > getCoarseSpace()
  {
    return m_CoarseSpace;
  }

private:
  LO m_numNodes;
  const LO *m_nodeBegin, *m_localDofs;
  const GO* m_nodeGlobalIDs;
  const SM *m_xCoord, *m_yCoord, *m_zCoord;
  const std::vector< std::vector<LO> > & m_subNodeBegin;
  const std::vector< std::vector<LO> > & m_subLocalDofs;
  const std::vector< std::vector<LO> > & m_subNodes;
  const std::vector< std::vector<LO> > & m_subRowBegin;
  const std::vector< std::vector<LO> > & m_subColumns;
  const std::vector< std::vector<SX> > & m_subValues;
  RCP<Teuchos::ParameterList> & m_Parameters;
  MPI_Comm m_mpiComm;
  LO m_spatialDim;
  enum bddc::ProblemType m_problemType; 
  RCP<const Teuchos::Comm<int> > m_Comm;
  int m_myPID, m_numThreads;
  Tpetra::global_size_t m_IGO;
  bool m_interfacePreconditioner, m_addCorners;
  LO m_numDofs, m_numSub, m_numDofsB, m_numDofsB1to1;
  RCP< bddc::PartitionOfUnity<SX,SM,LO,GO> > m_Partition;
  std::vector< bddc::SubdomainBDDC<SX,SM,LO,GO>* > m_Subdomain;
  RCP< bddc::CoarseSpaceBDDC<SX,SM,LO,GO> > m_CoarseSpace;
  RCP< bddc::WeightsBDDC<SX,SM,LO,GO> > m_Weights;
  RCP< bddc::ConstraintsBDDC<SX,SM,LO,GO> > m_Constraints;
  std::vector< std::vector<LO> > m_boundaryDofsLocal, m_subBoundaryDofs,
    m_subDofs, m_subInteriorDofs;
  std::vector<LO> m_boundaryDofs, m_globalToBoundaryMap, m_boundaryToAll1to1;
  std::vector<SM> m_diagBoundary;
  std::vector< std::vector<LO> > m_equivBoundaryDofs;
  RCP<const Map> m_dofMap, m_dofMap1to1, m_dofMapB, m_dofMapB1to1;
  RCP<Export> m_exporterB, m_exporterAll;
  RCP<Import> m_importerB, m_importerAll;
  RCP<Vector> m_rhsVecB, m_solVecB, m_rhsVecB1to1, m_solVecB1to1,
    m_coarseRhsVec1to1, m_coarseSolVec1to1, m_xVecAll, m_AxVecAll,
    m_xVecAll1to1, m_AxVecAll1to1, m_rhsVecAll, m_solVecAll,
    m_rhsVecAll1to1, m_solVecAll1to1;
  std::vector<SX> m_work1, m_work2;
  SX **m_work1T, **m_work2T;
  std::vector<double> m_timings;

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

  void initializeVectors()
  {
    m_rhsVecB = rcp( new Vector(m_dofMapB) );
    m_solVecB = rcp( new Vector(m_dofMapB) );
    m_rhsVecB1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_solVecB1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_coarseRhsVec1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_coarseSolVec1to1 = rcp( new Vector(m_dofMapB1to1) );
    m_xVecAll = rcp( new Vector(m_dofMap) );
    m_AxVecAll = rcp( new Vector(m_dofMap) );
    m_xVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_AxVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_rhsVecAll = rcp( new Vector(m_dofMap) );
    m_solVecAll = rcp( new Vector(m_dofMap) );
    m_rhsVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
    m_solVecAll1to1 = rcp( new Vector(m_dofMap1to1) );
  }

  void determineWeights()
  {
    m_Weights = 
      rcp( new bddc::WeightsBDDC<SX,SM,LO,GO>
	   (m_Subdomain, m_Partition, m_exporterB, m_importerB,
	    m_subBoundaryDofs, m_diagBoundary, m_Parameters) );
    m_Weights->determineWeights();
  }

  void determineBoundaryMaps()
  {
    std::vector<GO> globalIDs(m_numDofsB);
    for (LO i=0; i<m_numDofsB; i++) {
      LO dof = m_boundaryDofs[i];
      globalIDs[i] = m_dofMap->getGlobalElement(dof);
    }
    m_dofMapB = 
      rcp( new Map(m_IGO, Teuchos::ArrayView<GO>(globalIDs), 0, m_Comm) );
    m_dofMapB1to1 = Tpetra::createOneToOne<LO,GO,Node>(m_dofMapB);
    m_importerB = rcp( new Import(m_dofMapB1to1, m_dofMapB) );
    m_exporterB = rcp( new Export(m_dofMapB, m_dofMapB1to1) );
    m_numDofsB1to1 = m_dofMapB1to1->getNodeNumElements();
    m_importerAll = rcp( new Import(m_dofMap1to1, m_dofMap) );
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

  void determineCoarseSpace()
  {
    m_CoarseSpace = 
      rcp( new bddc::CoarseSpaceBDDC<SX,SM,LO,GO>
	   (m_Subdomain, m_Partition, m_exporterB, m_importerB,
	    m_equivBoundaryDofs, m_subBoundaryDofs, m_Parameters) );
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
	   (m_numNodes, m_nodeBegin, m_localDofs, m_subNodes, 
	    m_Subdomain, m_Partition, 
	    m_exporterB, m_importerB, m_diagBoundary, m_Parameters) );
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
    diagBoundary.doImport(diagBoundary1to1, *m_importerB, Tpetra::INSERT);
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
  
