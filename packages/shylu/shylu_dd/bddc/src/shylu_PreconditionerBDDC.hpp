
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

#ifndef BDDC_PRECONDITIONER_H
#define BDDC_PRECONDITIONER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include <set>

#include "Tpetra_Distributor.hpp"
#include "shylu_PartitionOfUnityBDDC.hpp"
#include "shylu_SubdomainBDDC.hpp"
#include "shylu_WeightsBDDC.hpp"
#include "shylu_ConstraintsBDDC.hpp"
#include "shylu_DofManager.hpp"
#include "shylu_UtilBDDC.hpp"
#include "shylu_ZoltanPartition.hpp"
#include "shylu_SolverFactoryBDDC.hpp"
#include "shylu_PComm.hpp"
#include "shylu_DisconnectedComponentChecker.hpp"
#include "shylu_OperatorBase.hpp"
#include "shylu_PreconditionerBase.hpp"
#include "shylu_errorBDDC.hpp"

//#define VERTEXCOARSESPACEBDDC
#ifdef VERTEXCOARSESPACEBDDC
#include "shylu_VertexCoarseSpace.hpp"
#endif


#ifdef _OPENMP
#include "omp.h"
#endif

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class PreconditionerBDDC : 
  public PreconditionerBase<SX,SM,LO,GO>
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::Map<LO,GO>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
  typedef Tpetra::CrsMatrix<SX,LO,GO>                        CrsMatrix;
  typedef Tpetra::CrsMatrix<LO,LO,GO>                        CrsMatrixLO;
  typedef Tpetra::CrsMatrix<GO,LO,GO>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO>                              Export;
  typedef Tpetra::Import<LO,GO>                              Import;
  typedef Tpetra::Vector<SM,LO,GO>                           VectorSM;
  typedef Tpetra::Vector<SX,LO,GO>                           Vector;
  typedef Tpetra::Vector<LO,LO,GO>                           VectorLO;
  typedef Tpetra::Vector<GO,LO,GO>                           VectorGO;
  typedef Tpetra::MultiVector<SX,LO,GO>                      MV;
  typedef Tpetra::MultiVector<double,LO,GO>                  MVD;

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
     OperatorBase<SX>* Operator=nullptr,     
     std::vector<int>* nodeSend=nullptr,
     std::vector<GO>* nodeGlobalIDs1to1=nullptr,
     RCP<const Map> dofMap=Teuchos::null,
     RCP<const Map> dofMap1to1=Teuchos::null
#ifdef VERTEXCOARSESPACEBDDC
     , RCP< VertexCoarseSpace<SX,SM,LO,GO> > vertexCoarseSpaceIn=Teuchos::null
#endif
     ) :
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
    m_interfacePreconditioner(Parameters->get("Interface Preconditioner", true)),
    m_usingSimpleInterface(false),
    m_useVertexCoarseSpace(Parameters->get("Use Vertex Coarse Space", false)),
    m_dofMap(dofMap), 
    m_dofMap1to1(dofMap1to1),
    m_level(level),
    m_Operator(Operator),
    m_usePComm(Parameters->get("Use PComm", true))
#ifdef VERTEXCOARSESPACEBDDC
    , m_vertexCoarseSpaceIn(vertexCoarseSpaceIn)
#endif
  {
    m_timings.resize(TIME_LENGTH);
    double startTime = GetTime();
    // check for disconnected subdomain components and make adjustments
    m_componentChecker = 
      rcp( new DisconnectedComponentChecker<LO,GO,SX>
	   (numNodes, nodeBegin, subNodes, subRowBegin, subColumns, subValues) );
    m_componentChecker->adjustProblemData
      (m_subNodes, m_subRowBegin, m_subColumns, m_subValues);
    m_Comm = rcp( new Teuchos::MpiComm<int>(mpiComm) );
    initialize(nodeGlobalIDs1to1, nodeSend);
    m_timings[TIME_INITIALIZATION] += GetTime() - startTime;
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
     std::vector<int>* nodeSend=nullptr,
     std::vector<GO>* nodeGlobalIDs1to1=nullptr) :
    m_numNodes(numNodes),
    m_nodeGlobalIDs(nodeGlobalIDs),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_Parameters(Parameters),
    m_mpiComm(mpiComm),
    m_interfacePreconditioner(Parameters->get("Interface Preconditioner", true)),
    m_usingSimpleInterface(true),
    m_useVertexCoarseSpace(Parameters->get("Use Vertex Coarse Space", false)),
    m_level(level),
    m_usePComm(Parameters->get("Use PComm", true))
  {
    m_timings.resize(TIME_LENGTH);
    double startTime = GetTime();
    m_Comm = rcp( new Teuchos::MpiComm<int>(mpiComm) );
    convertSingleSudomain(subRowBegin, subColumns, subValues);
    initialize(nodeGlobalIDs1to1, nodeSend);
    m_timings[TIME_INITIALIZATION] += GetTime() - startTime;
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
    for (size_t i=0; i<m_Subdomain.size(); i++) {
      delete m_Subdomain[i];
    }
    if (m_work1T != nullptr) {
#ifdef _OPENMP
      for (int i=0; i<m_numThreads; i++) {
	delete [] m_work1T[i];
	delete [] m_work2T[i];
      }
#endif
      delete [] m_work1T;
      delete [] m_work2T;
    }
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

  RCP<const Map> getDofMap()
  {
    return m_dofMap;
  }

  GO getCoarseSpaceDim() const
  {
    return m_numCoarseDofs;
  }

  RCP<const Teuchos::Comm<int> > getComm() const
  {
    return m_Comm;
  }

  LO getNumSub() const
  {
    return m_numSub;
  }

  GO getNumDof() const
  {
    if (m_pComm != Teuchos::null) {
      const LO numMyDof = m_pComm->getOwnedLength();
      LO numMyDofSum;
      Teuchos::reduceAll<int, LO> (*m_Comm, Teuchos::REDUCE_SUM, 1, 
				 &numMyDof, &numMyDofSum);
      return numMyDofSum;
    }
    else {
      return m_numDofs;
    }
  }

  const std::vector<LO> & getCoarseSizes()
  {
    retrieveCoarseSizesAcrossLevels();
    return m_coarseSizes;
  }

  void printTimings(const char* fileName)
  {
    retrieveTimingsAcrossLevels();
    if (m_myPID == 0) {
      std::ofstream timingsFile;
      timingsFile.open(fileName, std::ios::out);
      const LO numLevels = m_timingsMax.size()/TIME_LENGTH;
      timingsFile << "Timings (maximum) in seconds for BDDC preconditioner\n";
      for (LO i=0; i<numLevels; i++) {
	const double* timings = &m_timingsMax[i*TIME_LENGTH];
	timingsFile << "----------------------------------------------------\n";
	timingsFile << "Timings for level " << i << std::endl;
	timingsFile << "overall initialization time  = "
		    << timings[bddc::TIME_INITIALIZATION] << std::endl;
	timingsFile << "  partition of unity         = "
		    << timings[bddc::TIME_PARTITION_OF_UNITY] << std::endl;
	timingsFile << "  dof manager                = "
		    << timings[bddc::TIME_DOF_MANAGER] << std::endl;
	timingsFile << "  step 1 init group          = "
		    << timings[bddc::TIME_INIT_STEP1] << std::endl;
	timingsFile << "  coarse space preparations  = "
		    << timings[bddc::TIME_COARSE_SPACE_PREP] << std::endl;
	timingsFile << "  Note: next two times involve matrix factorizations\n";
	timingsFile << "  determine base constraints = "
		    << timings[bddc::TIME_BASE_CONSTRAINTS] << std::endl;
	timingsFile << "  factor interior matrices   = "
		    << timings[bddc::TIME_INTERIOR_FACTORIZATIONS] << std::endl;
	timingsFile << "solve phase:\n";
	timingsFile << "  apply preconditioner       = "
		    << timings[bddc::TIME_APPLY_PRECONDITIONER] << std::endl;
	timingsFile << "  initial static condensatio = "
		    << timings[bddc::TIME_INIT_STATIC_COND] << std::endl;
	timingsFile << "  static expansion           = "
		    << timings[bddc::TIME_STATIC_EXPANSION] << std::endl;
	timingsFile << "  subdomain corrections      = "
		    << timings[bddc::TIME_SUB_CORR] << std::endl;
	timingsFile << "  coarse solutions           = "
		    << timings[bddc::TIME_COARSE_CORR] << std::endl;
	timingsFile << "  apply operator             = "
		    << timings[bddc::TIME_APPLY_OPER] << std::endl;
	if (m_useVertexCoarseSpace == false) continue;
	timingsFile << "  VCS reduction              = "
		    << timings[bddc::TIME_VCS_REDUCTION] << std::endl;
	timingsFile << "  VCS export                 = "
		    << timings[bddc::TIME_VCS_EXPORT] << std::endl;
	timingsFile << "  VCS solve                  = "
		    << timings[bddc::TIME_VCS_SOLVE] << std::endl;
	timingsFile << "  VCS import                 = "
		    << timings[bddc::TIME_VCS_IMPORT] << std::endl;
	timingsFile << "  VCS expansion              = "
		    << timings[bddc::TIME_VCS_EXPANSION] << std::endl;
      }
      timingsFile.close();
    }
  }

  void printBanner(const char* fileName,
		   const bool resetFile=false)
  {
    const LO numProc = m_Comm->getSize();
    const LO numSub = getNumSub();
    LO numSubAll;
    Teuchos::reduceAll<int, LO> (*m_Comm, Teuchos::REDUCE_SUM, 1, 
				 &numSub, &numSubAll);
    LO numDofsArray[2], numDofsArraySum[2];
    numDofsArray[0] = m_numDofs;
    numDofsArray[1] = m_numDofsB;
    if (m_pComm != Teuchos::null) {
      numDofsArray[0] = m_pComm->getOwnedLength();
      numDofsArray[1] = m_pCommB->getOwnedLength();
    }
    Teuchos::reduceAll<int, LO> (*m_Comm, Teuchos::REDUCE_SUM, 2, 
				 numDofsArray, numDofsArraySum);
    const std::vector<LO> & coarseSizes = getCoarseSizes();
    if (m_myPID == 0) {
      std::ofstream outputFileDD;
      if (resetFile) {
	outputFileDD.open(fileName, std::ios::out);
      }
      else {
	outputFileDD.open(fileName, std::ios::out | std::ios::app);
      }
      outputFileDD << "number of processors     = " << numProc << std::endl;
      outputFileDD << "number of subdomains     = " << numSubAll << std::endl;
      outputFileDD << "number of unknowns       = " << numDofsArraySum[0] << std::endl;
      std::string text = "no";
      if (m_interfacePreconditioner) text = "yes";
      outputFileDD << "interface preconditioner = " << text << std::endl;
      if (m_interfacePreconditioner) {
	outputFileDD << "size of interface        = " << numDofsArraySum[1] << std::endl;
      }
      if (coarseSizes.size() > 0) {
	outputFileDD << "coarse problem sizes     = ";
	for (size_t i=0; i<coarseSizes.size()/2; i++) {
	  outputFileDD << coarseSizes[2*i];
	  const LO numVertexDof = coarseSizes[2*i+1];
	  if (numVertexDof > 0) {
	    outputFileDD << "(" << numVertexDof << ")";
	  }
	  outputFileDD << " ";
	}
	outputFileDD << " " << std::endl;
      }
      std::string solverName = 
	m_Parameters->get("Dirichlet Solver", m_defaultSolver);
      outputFileDD << "Dirichlet solver         = " <<  solverName << std::endl;
      solverName = m_Parameters->get("Neumann Solver", m_defaultSolver);
      outputFileDD << "Neumann solver           = " <<  solverName << std::endl;
      solverName = m_Parameters->get("Coarse Solver", m_defaultSolver);
      outputFileDD << "coarse solver            = " <<  solverName << std::endl;
      text = "no";
      if (m_useVertexCoarseSpace) text = "yes";
      outputFileDD << "vertex coarse space      = " << text << std::endl;
      int krylovMethod = m_Parameters->get("Krylov Method", 0);
      text = "GCR";
      if (krylovMethod == 1) text = "PCG";
      outputFileDD << "Krylov Method            = " << text << std::endl;
      SM solverTol = m_Parameters->get("Convergence Tolerance", 1e-6);
      LO maxIter = m_Parameters->get("Maximum Iterations", 200);
      outputFileDD << "solver tolerance         = " << solverTol << std::endl;
      outputFileDD << "maximum iterations       = " << maxIter << std::endl;
      outputFileDD.close();
    }
  }

  void exportAll(const SX* sourceVals,
		 SX* targetVals)
  {
    if (m_usePComm) {
      if (m_pComm == Teuchos::null) {
	memcpy(targetVals, sourceVals, m_numDofs*sizeof(SX));
	return;
      }
      BDDC_TEST_FOR_EXCEPTION(m_pComm->getSourceLength() != m_numDofs, 
			      std::runtime_error, "inconsistent length in exportAll");
      m_pComm->doExport(sourceVals, targetVals);
    }
    else {
      if (m_rhsVecAll == Teuchos::null) {
	memcpy(targetVals, sourceVals, m_numDofs*sizeof(SX));
	return;
      }
      Teuchos::ArrayRCP<SX> sourceValues = m_rhsVecAll->getDataNonConst();
      for (size_t i=0; i<m_dofMap->getNodeNumElements(); i++) {
	sourceValues[i] = sourceVals[i];
      }
      m_rhsVecAll1to1->putScalar(0);
      m_rhsVecAll1to1->doExport(*m_rhsVecAll, *m_exporterAll, Tpetra::ADD);
      Teuchos::ArrayRCP<const SX> targetValues = m_rhsVecAll1to1->getData();
      for (size_t i=0; i<m_dofMap1to1->getNodeNumElements(); i++) {
	targetVals[i] = targetValues[i];
      }
    }
  }

  void importAll(const SX* sourceVals,
		 SX* targetVals)
  {
    if (m_usePComm) {
      if (m_pComm == Teuchos::null) {
	memcpy(targetVals, sourceVals, m_numDofs*sizeof(SX));
	return;
      }
      BDDC_TEST_FOR_EXCEPTION(m_pComm->getSourceLength() != m_numDofs, 
			      std::runtime_error, "inconsistent length in importAll");
      m_pComm->doImport(sourceVals, targetVals);
    }
    else {
      if (m_rhsVecAll1to1 == Teuchos::null) {
	memcpy(targetVals, sourceVals, m_numDofs*sizeof(SX));
	return;
      }
      Teuchos::ArrayRCP<SX> sourceValues = m_rhsVecAll1to1->getDataNonConst();
      for (size_t i=0; i<m_dofMap1to1->getNodeNumElements(); i++) {
	sourceValues[i] = sourceVals[i];
      }
      m_rhsVecAll->doImport(*m_rhsVecAll1to1, *m_exporterAll, Tpetra::INSERT);
      Teuchos::ArrayRCP<const SX> targetValues = m_rhsVecAll->getData();
      for (size_t i=0; i<m_dofMap->getNodeNumElements(); i++) {
	targetVals[i] = targetValues[i];
      }
    }
  }

  int InitialStaticCondensation(SX* rightHandSide1to1,
				SX* initialSolution1to1)
  {
    if (m_coarseSub != Teuchos::null) return 0;
    double startTime = GetTime();
    SX* rhsB(0), *rhsB1to1(0);
    if (m_usePComm) {
      rhsB = m_pCommB->getSourcePtr();
      rhsB1to1 = m_pCommB->getOwnedPtr();
    }
    else {
      rhsB = m_rhsVecB->getDataNonConst().get();
      rhsB1to1 = m_rhsVecB1to1->getDataNonConst().get();
    }
    for (LO i=0; i<m_numDofsB; i++) rhsB[i] = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(m_numThreads)
#endif
    {
#ifdef _OPENMP
#pragma omp for
#endif
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & interiorDofs = m_subInteriorDofs[i];
	LO numInteriorDofs = interiorDofs.size();
	const std::vector<LO> & boundaryDofs = m_subBoundaryDofs[i];
	LO numBoundaryDofs = boundaryDofs.size();
	SX *subRhs(0), *subSol(0);
	int threadID = getThreadID();
	getArrays(threadID, subRhs, subSol);
	for (LO j=0; j<numInteriorDofs; j++) {
	  LO row = interiorDofs[j];
	  subRhs[j] = rightHandSide1to1[row];
	}
	m_Subdomain[i]->makeInteriorAdjustments(&subRhs[0], &subSol[0]);
	// no potential write conflicts below because interior dofs are unique
	// to each thread
	for (LO j=0; j<numInteriorDofs; j++) {
	  LO row = interiorDofs[j];
	  initialSolution1to1[row] += subSol[j];
	  if (m_interfacePreconditioner == false) {
	    rightHandSide1to1[row] -= subRhs[j+numBoundaryDofs];
	  }
	}
	// potential write conflicts because boundary dofs are common
	// to different threads
	for (LO j=0; j<numBoundaryDofs; j++) {
	  LO row = boundaryDofs[j];
	  BDDC_TEST_FOR_EXCEPTION(row >= m_numDofsB, std::runtime_error, 
				  "invalid row index");
#ifdef _OPENMP
#pragma omp atomic
#endif
	  rhsB[row] -= subRhs[j];
	}
      }
    }
    if (m_usePComm) {
      m_pCommB->doExport(rhsB, rhsB1to1);
    }
    else {
      m_rhsVecB1to1->putScalar(0);
      m_rhsVecB1to1->doExport(*m_rhsVecB, *m_exporterB, Tpetra::ADD);
    }
    for (LO i=0; i<m_numDofsB1to1; i++) {
      LO row = m_boundaryToAll1to1[i];
      rhsB1to1[i] += rightHandSide1to1[row];
    }
    for (LO i=0; i<m_numDofsB1to1; i++) {
      LO row = i;
      if (m_interfacePreconditioner == false) row = m_boundaryToAll1to1[i];
      rightHandSide1to1[row] = rhsB1to1[i];
    }
    m_timings[TIME_INIT_STATIC_COND] += GetTime() - startTime;
    return 1;
  }

  void addBoundaryValues(const SX* solInput, 
			 const bool interfaceValuesOnlyInput,
			 SX* solutionAll)
  {
    if (interfaceValuesOnlyInput) {
      for (LO i=0; i<m_numDofsB1to1; i++) {
	const LO row = m_boundaryToAll1to1[i];
	solutionAll[row] += solInput[i];
      }
    }
    else {
      for (LO i=0; i<m_numDofsB1to1; i++) {
	const LO row = m_boundaryToAll1to1[i];
	solutionAll[row] += solInput[row];
      }
    }
  }

  void getSolutionOnBoundary(const SX* solInput,
			     const bool interfaceValuesOnlyInput,
			     SX* & solB)
  {
    SX* solB1to1(nullptr);
    if (m_usePComm) {
      solB = m_pCommB->getSourcePtr();
      solB1to1 = m_pCommB->getOwnedPtr();
    }
    else {
      solB = m_solVecB->getDataNonConst().get();
      solB1to1 = m_solVecB1to1->getDataNonConst().get();
    }
    if (interfaceValuesOnlyInput) {
      for (LO i=0; i<m_numDofsB1to1; i++) {
	solB1to1[i] = solInput[i];
      }
    }
    else {
      for (LO i=0; i<m_numDofsB1to1; i++) {
	solB1to1[i] = solInput[m_boundaryToAll1to1[i]];
      }
    }
    if (m_usePComm) {
      m_pCommB->doImport(solB1to1, solB);
    }
    else {
      m_solVecB->doImport(*m_solVecB1to1, *m_exporterB, Tpetra::INSERT);
    }
  }

  void StaticExpansion(SX* deltaSolution,
		       const bool interfaceValuesOnlyInput,
		       const bool useCurrentInteriorValues,
		       SX* initialSolution)
  {
    // Note: The dimension of deltaSolution is as follows:
    //  interfaceValuesOnlyInput == true:  the number of owned boundary dofs
    //  interfaceValuesOnlyInput == false: the number of all owned dofs
    //    The dimension of initialSolution is always the number of owned dofs
    double startTime = GetTime();
    SX *solB(0);
    getSolutionOnBoundary(deltaSolution, interfaceValuesOnlyInput, solB);
    for (LO i=0; i<m_numSub; i++) {
      const std::vector<LO> & subBoundaryDofs = m_subBoundaryDofs[i];
      const std::vector<LO> & interiorDofs = m_subInteriorDofs[i];
      LO numDofsB = subBoundaryDofs.size();
      LO numDofsI = interiorDofs.size();
      SX *xB(0), *xI(0);
      int threadID = getThreadID();
      getArrays(threadID, xB, xI);
      for (LO j=0; j<numDofsB; j++) {
	const LO row = subBoundaryDofs[j];
	xB[j] = solB[row];
      }
      if (useCurrentInteriorValues) {
	for (LO j=0; j<numDofsI; j++) {
	  const LO row = interiorDofs[j];
	  initialSolution[row] += deltaSolution[row];
	  xI[j] = deltaSolution[row];
	}
      }
      m_Subdomain[i]->staticExpansion(&xB[0], &xI[0], useCurrentInteriorValues);
      // no potential write conflicts because interior dofs are unique
      // to each thread
      for (LO j=0; j<numDofsI; j++) {
	LO row = interiorDofs[j];
	initialSolution[row] += xI[j];
      }
    }
    addBoundaryValues(deltaSolution, interfaceValuesOnlyInput,
		      initialSolution);
    m_timings[TIME_STATIC_EXPANSION] += GetTime() - startTime;
  }

  LO NumMyRows()
  {
    if (m_pComm != Teuchos::null) {
      return m_pComm->getOwnedLength();
    }
    else {
      return m_numDofs;
    }
  }

  LO NumMyRowsKrylov()
  {
    if (m_interfacePreconditioner == true) {
      return m_numDofsB1to1;
    }
    else {
      if (m_pComm != Teuchos::null) {
	return m_pComm->getOwnedLength();
      }
      else {
	return m_numDofs;
      }
    }
  }

  void Apply(const SX* r, 
	     SX* Pr,
	     SX* APr,
	     const bool interfaceValuesOnly)
  {
    Apply2(r, Pr, interfaceValuesOnly);
    double startTime = GetTime();
    if (interfaceValuesOnly) {
      ApplyOperator(Pr, APr);
    }
    else {
      ApplyFullOperator(Pr, APr);
    }
    m_timings[TIME_APPLY_OPER] += GetTime() - startTime;
  }

  void Apply(const SX* r, 
	     SX* Pr,
	     const bool interfaceValuesOnly)
  {
    double startTime = GetTime();
    // Quick return if direct solver used for current level.
    if (m_useDirectSolver == true) {
      applyDirectSolver(r, Pr);
      m_timings[TIME_APPLY_PRECONDITIONER] += GetTime() - startTime;
      return;
    }
    // Import all right-hand-side values (rImport always has both
    // interior and interface values whether or not an interface
    // preconditioner is being used).
    SX* rImport(nullptr);
    importRhs(r, interfaceValuesOnly, rImport);
    // restrict residual to next coarser level
    applyPhiTranspose(rImport, m_coarseRhsWork);
    // The following if block may be executed by MPI processes belonging
    // to a disjoint communicator, thus permitting asynchronous work.
    if (m_coarseProc == true) {
      applyCoarseCorrection(m_coarseRhsWork, m_coarseSolWork);
    }
    // The next step performs an initial static condensation correction
    // as needed. With reference to the paper "An approximate BDDC 
    // preconditioner", Numer. Linear Algebra Appl. 2007: 14:149-168, 
    // initialCorrection applies (if needed) P_1 to the residual and updates 
    // the residual according to
    //  u_0 = P_1r
    //  r <- (I - AP_1)r (see (5) or (39) of cited paper).
    SX *deltaCorrection(nullptr), *rUse(nullptr);
    initialCorrection(rImport, rUse, Pr, deltaCorrection,
		      interfaceValuesOnly);
    applyLocalCorrections(rUse, deltaCorrection); 
    // Sum coarse and local (subdomain) corrections.
    // Note: we need a barrier here to ensure applyCoarseCorrection has
    //       completed (could have been done asynchronously)
    m_Comm->barrier();
    SX* coarseSol(nullptr);
    importCoarseSolution(m_coarseSolWork, coarseSol);
    addCoarseCorrection(coarseSol, deltaCorrection);
    if (interfaceValuesOnly) {
      restrictToInterface(deltaCorrection);
      exportBoundaryValues(deltaCorrection, Pr);
    }
    else {
      exportAllValues(deltaCorrection);
      // this part of the code applies either the operation
      // Pr = P_1r + (I - P_1A)P_2(I - AP_1)r or
      // Pr = P_1r + (I - P_1A)P_2r, depending on whether or not
      //                             an initial correction was made
      //                             see (5) or (39) of cited paper
      const bool useCurrentInteriorValues = true;
      StaticExpansion(deltaCorrection, interfaceValuesOnly, 
		      useCurrentInteriorValues, Pr);
    }
    m_timings[TIME_APPLY_PRECONDITIONER] += GetTime() - startTime;
  }
  
  void Apply2(const SX* r, 
	      SX* Pr,
	      const bool interfaceValuesOnly)
  {
    double startTime = GetTime();
    // Quick return if direct solver used for current level.
    if (m_useDirectSolver == true) {
      applyDirectSolver(r, Pr);
      m_timings[TIME_APPLY_PRECONDITIONER] += GetTime() - startTime;
      return;
    }
    // Import all right-hand-side values (rImport always has both
    // interior and interface values whether or not an interface
    // preconditioner is being used).
    SX* rImport(nullptr);
    importRhs(r, interfaceValuesOnly, rImport);
    // The next step performs an initial static condensation correction
    // as needed. With reference to the paper "An approximate BDDC 
    // preconditioner", Numer. Linear Algebra Appl. 2007: 14:149-168, 
    // initialCorrection applies (if needed) P_1 to the residual and updates 
    // the residual according to
    //  u_0 = P_1r
    //  r <- (I - AP_1)r (see (5) or (39) of cited paper).
    SX *deltaCorrection(nullptr), *rUse(nullptr);
    initialCorrection(rImport, rUse, Pr, deltaCorrection,
		      interfaceValuesOnly);
    // restrict residual to next coarser level
    applyPhiTranspose(rImport, m_coarseRhsWork);
    // The following if block may be executed by MPI processes belonging
    // to a disjoint communicator, thus permitting asynchronous work.
    if (m_coarseProc == true) {
      applyCoarseCorrection(m_coarseRhsWork, m_coarseSolWork);
    }
    applyLocalCorrections(rUse, deltaCorrection); 
    // Sum coarse and local (subdomain) corrections.
    // Note: we need a barrier here to ensure applyCoarseCorrection has
    //       completed (could have been done asynchronously)
    m_Comm->barrier();
    SX* coarseSol(nullptr);
    importCoarseSolution(m_coarseSolWork, coarseSol);
    addCoarseCorrection(coarseSol, deltaCorrection);
    if (interfaceValuesOnly) {
      restrictToInterface(deltaCorrection);
      exportBoundaryValues(deltaCorrection, Pr);
    }
    else {
      exportAllValues(deltaCorrection);
      // this part of the code applies either the operation
      // Pr = P_1r + (I - P_1A)P_2(I - AP_1)r or
      // Pr = P_1r + (I - P_1A)P_2r, depending on whether or not
      //                             an initial correction was made
      //                             see (5) or (39) of cited paper
      const bool useCurrentInteriorValues = true;
      StaticExpansion(deltaCorrection, interfaceValuesOnly, 
		      useCurrentInteriorValues, Pr);
    }
    m_timings[TIME_APPLY_PRECONDITIONER] += GetTime() - startTime;
  }
  
  void exportBoundaryValues(const SX* deltaB, 
			    SX* Pr)
  {
    if (m_usePComm) {
      m_pCommB->doExport(deltaB, Pr);
    }
    else {
      Teuchos::ArrayRCP<SX> valsB = m_rhsVecB->getDataNonConst();
      for (LO i=0; i<m_numDofsB; i++) {
	valsB[i] = deltaB[i];
      }
      m_rhsVecB1to1->putScalar(0);
      m_rhsVecB1to1->doExport(*m_rhsVecB, *m_exporterB, Tpetra::ADD);
      Teuchos::ArrayRCP<const SX> valsB1to1 = m_rhsVecB1to1->getData();
      for (LO i=0; i<m_numDofsB1to1; i++) {
	Pr[i] = valsB1to1[i];
      }
    }
  }

  void exportAllValues(SX* allVals)
  {
    if (m_usePComm) {
      const LO ownedLength = m_pComm->getOwnedLength();
      SX* ownedVals = m_pComm->getOwnedPtr();
      m_pComm->doExport(allVals, ownedVals);
      memcpy(allVals, ownedVals, ownedLength*sizeof(SX));
    }
    else {
      Teuchos::ArrayRCP<SX> rhsData = m_rhsVecAll->getDataNonConst();
      const LO numRows = m_rhsVecAll->getLocalLength();
      for (LO i=0; i<numRows; i++) rhsData[i] = allVals[i];
      m_rhsVecAll1to1->putScalar(0);
      m_rhsVecAll1to1->doExport(*m_rhsVecAll, *m_exporterAll, Tpetra::ADD);
      Teuchos::ArrayRCP<const SX> rhsData1to1 = m_rhsVecAll1to1->getData();
      const LO numRows1to1 = m_rhsVecAll1to1->getLocalLength();
      for (LO i=0; i<numRows1to1; i++) allVals[i] = rhsData1to1[i];
    }
  }

  void restrictToInterface(SX* sol)
  {
    m_rhsWork.resize(m_numDofsB);
    for (LO i=0; i<m_numDofsB; i++) {
      const LO row = m_boundaryDofs[i];
      m_rhsWork[i] = sol[row];
    }
    memcpy(sol, m_rhsWork.data(), m_numDofsB*sizeof(SX));
  }

  void applyLocalCorrections(const SX* rhs, 
			     SX* sol)
  {
    // Note: on entry sol has all zero entries
    double startTime = GetTime();
    if (m_disjointCommunicators) {
      if (m_coarseProc == false) {
	localCorrections(rhs, sol);
      }
    }
    else {
      localCorrections(rhs, sol);
    }
    m_timings[TIME_SUB_CORR] += GetTime() - startTime;
  }

  void initialCorrection(SX* rImport,
			 SX* & rUse,
			 SX* initCorrection,
			 SX* & deltaCorrection, 
			 const bool interfaceValuesOnly)
  {
    // This function applies an initial static condensation correction
    // (possibly approximate) as needed. The pointer rUse contains the
    // subdomain residual (either all dofs or only those on the interface)
    rUse = rImport;
    const LO numRows = m_pComm->getSourceLength();
    m_deltaSol.assign(numRows, 0);
    deltaCorrection = m_deltaSol.data();
    const bool needCorrection = 
      needInitialStaticCondensationCorrection(interfaceValuesOnly);
    if (needCorrection) {
      SX* initCorrectionAll = m_pComm->getSourcePtr();
      memset(initCorrectionAll, 0, numRows*sizeof(SX));
      // Note: Both rImport and initCorrectionAll can be modified below,
      //       but no parallel communications takes place.
      initialCorrection(rImport, initCorrectionAll);
      m_pComm->doExportLocal(initCorrectionAll, initCorrection);
    }
    else {
      LO numRowsOwned = m_pComm->getOwnedLength();
      if (interfaceValuesOnly) numRowsOwned = m_pCommB->getOwnedLength();
      memset(initCorrection, 0, numRowsOwned*sizeof(SX));
    }
  }

  void initialCorrection(SX* rImport,
			 SX* initCorrection)
  {
    // Note: rImport on entry has all the processor dofs (not just owned 
    //       or interface ones). Likewise, initCorrection is for all
    //       processor dofs.
    double startTime = GetTime();
    const LO numRows = m_pComm->getSourceLength();
    m_residVec2.assign(numRows, 0);
    SX* rImportDelta = m_residVec2.data();
    for (LO i=0; i<m_numSub; i++) {
      const std::vector<LO> & subDofs = m_subDofs[i];
      const std::vector<LO> & interiorDofs = m_Subdomain[i]->getInteriorDofs();
      const LO* boundaryDofs = m_Subdomain[i]->getBoundaryDofs();
      LO numInteriorDofs = interiorDofs.size();
      LO numBoundaryDofs = m_Subdomain[i]->getNumBoundaryDofs();
      SX *subRhs(0), *subSol(0);
      int threadID = getThreadID();
      getArrays(threadID, subRhs, subSol);
      for (LO j=0; j<numInteriorDofs; j++) {
	LO row = subDofs[interiorDofs[j]];
	subRhs[j] = rImport[row];
      }
      m_Subdomain[i]->makeInteriorAdjustments(&subRhs[0], &subSol[0],
					      !m_interfacePreconditioner);
      // no potential write conflicts below because interior dofs are unique
      // to each thread
      // Note: first numBoundaryDofs in subRhs are for the boundary dofs
      for (LO j=0; j<numInteriorDofs; j++) {
	LO row = subDofs[interiorDofs[j]];
	initCorrection[row] = subSol[j];
	if (m_interfacePreconditioner == false) {
	  rImportDelta[row] -= subRhs[j+numBoundaryDofs];
	}
	else {
	  rImportDelta[row] = -rImport[row]; // final internal residual is 0
	}
      }
      // potential write conflicts because boundary dofs are common
      // to different threads
      for (LO j=0; j<numBoundaryDofs; j++) {
	LO row = subDofs[boundaryDofs[j]];
#ifdef _OPENMP
#pragma omp atomic
#endif
	rImportDelta[row] -= subRhs[j];
      }
    }
    // TODO: consider removing these communications and update rImport
    //       directly in the loop above
    const LO numRowsOwned = m_pComm->getOwnedLength();
    m_residVec3.assign(numRowsOwned, 0);
    m_pComm->doExport(rImportDelta, m_residVec3.data());
    m_pComm->doImport(m_residVec3.data(), rImportDelta);
    for (LO i=0; i<numRows; i++) {
      rImport[i] += rImportDelta[i];
    }
    
    m_timings[TIME_INIT_CORRECTION] += GetTime() - startTime;
  }

  void importRhs(const SX* rhs,
		 const bool interfaceValuesOnly,
		 SX* & rhsImport)
  {
    const LO numRows = m_pComm->getSourceLength();
    m_residualVector.assign(numRows, 0);
    rhsImport = m_residualVector.data();
    if (m_usePComm) {
      if (interfaceValuesOnly) {
	m_rhsWork.resize(m_numDofsB);
	m_pCommB->doImport(rhs, m_rhsWork.data());
	for (LO i=0; i<m_numDofsB; i++) {
	  const LO row = m_boundaryDofs[i];
	  rhsImport[row] = m_rhsWork[i];
	}
      }
      else {
	m_pComm->doImport(rhs, rhsImport);
      }
    }
    else {
      RCP<Vector> rhsVec1to1, rhsVec;
      RCP<Export> exporter;
      if (interfaceValuesOnly) {
	rhsVec1to1 = m_rhsVecB1to1;
	rhsVec = m_rhsVecB;
	exporter = m_exporterB;
      }
      else {
	rhsVec1to1 = m_rhsVecAll1to1;
	rhsVec = m_rhsVecAll;
	exporter = m_exporterAll;
      }
      const LO numRows1to1 = rhsVec1to1->getLocalLength();
      Teuchos::ArrayRCP<SX> rhsData1to1 = rhsVec1to1->getDataNonConst();
      for (LO i=0; i<numRows1to1; i++) rhsData1to1[i] = rhs[i];
      rhsVec->doImport(*rhsVec1to1, *exporter, Tpetra::INSERT);
      const LO numRows = rhsVec->getLocalLength();
      Teuchos::ArrayRCP<const SX> rhsData = rhsVec->getData();
      if (interfaceValuesOnly) {
	BDDC_TEST_FOR_EXCEPTION(numRows != m_numDofsB, std::runtime_error, 
				"incorrect numRows");
	for (LO i=0; i<numRows; i++) {
	  const LO row = m_boundaryDofs[i];
	  rhsImport[row] = rhsData[i];
	}
      }
      else {
	for (LO i=0; i<numRows; i++) rhsImport[i] = rhsData[i];
      }
    }
  }

  void applyDirectSolver(const SX* rhs,
			 SX* sol)
  {
    m_residualVector.resize(m_numDofs);
    SX* rhs2 = m_residualVector.data();
    memcpy(rhs2, rhs, m_numDofs*sizeof(SX));
    if (m_useVertexCoarseSpace) {
#ifdef VERTEXCOARSESPACEBDDC
      if (m_vertexCoarseSpaceIn != Teuchos::null) {
	m_vertexCoarseSpaceIn->apply(rhs2, sol);
      }
#endif
    }
    else {
      BDDC_TEST_FOR_EXCEPTION(m_coarseSub == Teuchos::null, std::runtime_error, 
			      "m_coarseSub is nullptr");
      m_coarseSub->Solve(1, rhs2, sol);
    }
  }

  void applyCoarseCorrection(std::vector<SX> & coarseRhs, 
			     std::vector<SX> & coarseSol)
  {
    double startTime = GetTime();
    const bool interfaceValuesOnly = false;
    m_coarsePreconditioner->Apply(coarseRhs.data(), coarseSol.data(), 
				  interfaceValuesOnly);
    m_timings[TIME_COARSE_CORR] += GetTime() - startTime;
  }

  void localCorrections(const SX* rhs, 
			SX* sol)
  {
    // rhs and sol are for all proc dofs
    for (LO i=0; i<m_numSub; i++) {
      const std::vector<LO> & subDofs = m_subDofs[i];
      LO numDofSub = subDofs.size();
      SX *subRhs(0), *subSol(0);
      int threadID = getThreadID();
      getArrays(threadID, subRhs, subSol);
      for (LO j=0; j<numDofSub; j++) {
	subRhs[j] = rhs[subDofs[j]];
      }
      const bool restrictToBoundary = false;
      m_Subdomain[i]->applyNeumannCorrection(subRhs, subSol,
					     restrictToBoundary);
      for (LO j=0; j<numDofSub; j++) {
	int row = subDofs[j];
#ifdef _OPENMP
#pragma omp atomic
#endif
	sol[row] += subSol[j];
      }
    }
  }

  void ApplyOperator(SX* x, 
		     SX* Ax)
  {
    if (m_interfacePreconditioner) {
      RCP<Vector> & AxVec = m_solVecB;
      SX *xValues1to1(0), *xValues(0), *AxValues(0);
      if (m_usePComm) {
	xValues = m_pCommB->getSourcePtr();
	m_pCommB->doImport(x, xValues);
	AxValues = m_pCommB->getSourcePtr2();
      }
      else {
	RCP<Vector> & xVec1to1 = m_rhsVecB1to1;
	xValues1to1 = xVec1to1->getDataNonConst().get();
	for (LO i=0; i<m_numDofsB1to1; i++) xValues1to1[i] = x[i];
	RCP<Vector> & xVec = m_rhsVecB;
	xVec->doImport(*xVec1to1, *m_exporterB, Tpetra::INSERT);
	xValues = xVec->getDataNonConst().get();
	AxValues = AxVec->getDataNonConst().get();
      }
      for (LO i=0; i<m_numDofsB; i++) AxValues[i] = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(m_numThreads)
#endif
      {
#ifdef _OPENMP
#pragma omp for
#endif
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
#ifdef _OPENMP
#pragma omp atomic
#endif
	      AxValues[row] += AxSub[j];
	    }
	  }
	}
      }
      SX *AxValues1to1(0);
      if (m_usePComm) {
	AxValues1to1 = m_pCommB->getOwnedPtr();
	m_pCommB->doExport(AxValues, AxValues1to1);
      }
      else {
	RCP<Vector> & AxVec1to1 = m_solVecB1to1;
	AxVec1to1->putScalar(0);
	AxVec1to1->doExport(*AxVec, *m_exporterB, Tpetra::ADD);
	AxValues1to1 = AxVec1to1->getDataNonConst().get();
      }
      for (LO i=0; i<m_numDofsB1to1; i++) Ax[i] = AxValues1to1[i];
    }
    else {
      ApplyFullOperator(x, Ax);
    }
  }

  void getSubDofs(const LO sub,
		  const LO* & subInteriorDofs, 
		  LO & numSubInteriorDofs, 
		  const LO* & subBoundaryDofs,
		  LO & numSubBoundaryDofs)
  {
    const std::vector<LO> & interiorDofs = m_Subdomain[sub]->getInteriorDofs();
    subInteriorDofs = interiorDofs.data();
    numSubInteriorDofs = interiorDofs.size();
    subBoundaryDofs = m_Subdomain[sub]->getBoundaryDofs();
    numSubBoundaryDofs = m_Subdomain[sub]->getNumBoundaryDofs();
  }

  void ApplyFullOperatorSimple(SX* x, 
			       SX* Ax)
  {
    BDDC_TEST_FOR_EXCEPTION(m_numSub != 1, std::runtime_error, 
			    "number of subdomains is not 1");
    if (m_Operator != nullptr) {
      m_Operator->Apply(x, Ax);
    }
    else {
      if (m_subDofs.size() != 1) {
	m_subDofs.resize(1);
	getSubDofs(0, m_subDofs[0]); 
      }
      const LO* rowBegin = m_subRowBegin[0];
      const LO* columns = m_subColumns[0];
      const SX* values = m_subValues[0];
      const std::vector<LO> & subDofs = m_subDofs[0];
      const LO numRows = subDofs.size();
      for (LO i=0; i<numRows; i++) {
	SX sum(0);
	for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	  const LO col = subDofs[columns[j]];
	  sum += values[j]*x[col];
	}
	const LO row = subDofs[i];
	Ax[row] = sum;
      }
    }
  }

  void ApplyFullOperator(SX* x, 
			 SX* Ax)
  {
    if (m_coarseSub != Teuchos::null) {
      ApplyFullOperatorSimple(x, Ax);
      return;
    }
    SX *xValues(nullptr), *AxValues(nullptr);
    if (m_usePComm) {
      BDDC_TEST_FOR_EXCEPTION(m_pComm->getSourceLength() != m_numDofs, 
			      std::runtime_error, "incorrect length");
      xValues = m_pComm->getSourcePtr();
      AxValues = m_pComm->getSourcePtr2();
      m_pComm->doImport(x, xValues);
    }
    else {
      Teuchos::ArrayRCP<SX> xValues1to1 = m_xVecAll1to1->getDataNonConst();
      LO numDofs1to1 = m_dofMap1to1->getNodeNumElements();
      for (LO i=0; i<numDofs1to1; i++) xValues1to1[i] = x[i];
      m_xVecAll->doImport(*m_xVecAll1to1, *m_exporterAll, Tpetra::INSERT);
      xValues = m_xVecAll->getDataNonConst().get();
      AxValues = m_AxVecAll->getDataNonConst().get();
    }
    if (m_Operator != nullptr) {
      m_Operator->Apply(xValues, AxValues);
    }
    else {
      for (LO i=0; i<m_numDofs; i++) AxValues[i] = 0;
#ifdef _OPENMP
#pragma omp for
#endif
      for (LO i=0; i<m_numSub; i++) {
	const std::vector<LO> & subDofs = m_subDofs[i];
	LO numDofs = subDofs.size();
	SX *xSub(0), *AxSub(0);
	int threadID = getThreadID();
	getArrays(threadID, xSub, AxSub);
	for (LO j=0; j<numDofs; j++) {
	  xSub[j] = xValues[subDofs[j]];
	}
	m_Subdomain[i]->applyFullOperator(xSub, AxSub);
	//#pragma omp critical(apply_full_operator)
	{
	  for (LO j=0; j<numDofs; j++) {
	    int row = subDofs[j];
#ifdef _OPENMP
#pragma omp atomic
#endif
	    AxValues[row] += AxSub[j];
	  }
	}
      }
    }
    if (m_usePComm) {
      m_pComm->doExport(AxValues, Ax);
    }
    else {
      m_AxVecAll1to1->putScalar(0);
      m_AxVecAll1to1->doExport(*m_AxVecAll, *m_exporterAll, Tpetra::ADD);
      Teuchos::ArrayRCP<SX> AxValuesAll1to1 = m_AxVecAll1to1->getDataNonConst();
      LO numDofs1to1 = m_dofMap1to1->getNodeNumElements();
      for (LO i=0; i<numDofs1to1; i++) Ax[i] = AxValuesAll1to1[i];
    }
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
    return sqrt(std::abs(dotprod));
  }

  SX DotProd(SX* x, 
	     SX* y, 
	     LO numTerm)
  {
    SX localDotProd(0), dotProd;
    for (LO i=0; i<numTerm; i++) {
      //      localDotProd += std::conj(x[i])*y[i];
      localDotProd += x[i]*y[i];
    }
    Teuchos::reduceAll<int, SX> (*m_Comm, Teuchos::REDUCE_SUM, 1, 
				 &localDotProd, &dotProd);
    return dotProd;
  }

  int MyPID() 
  {
    return m_myPID;
  }

  int numProc()
  {
    return m_Comm->getSize();
  }

  double GetTime()
  {
    return MPI_Wtime();
  }

  struct SubdomainData {
    LO numNodes{0}, numRows{0}, numFaces{0};
    LO *nodeBegin{0}, *localDofs{0};
    GO *nodeGIDs{0}, *faceGIDs{0};
    SX *values{0};
    SM *coords{0};
    std::vector<LO> nodeLIDs, faceLIDs;
    int* minPart{0};
  };

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
  bool m_interfacePreconditioner, m_usingSimpleInterface,
    m_useVertexCoarseSpace;
  LO m_numDofs, m_numSub, m_numDofsB, m_numDofsB1to1;
  std::vector<int> m_coarseMpiRanks;
  std::vector< std::vector<LO> > m_subNodeBegin, m_subLocalDofs;
  RCP< bddc::PartitionOfUnity<SX,SM,LO,GO> > m_Partition;
  std::vector< bddc::SubdomainBDDC<SX,SM,LO,GO>* > m_Subdomain;
  RCP< SubdomainBDDC<SX,SM,LO,GO> > m_coarseSub;
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
  std::vector<SX> m_work1, m_work2, m_rhsWork, m_rhsWork2;
  SX **m_work1T, **m_work2T;
  std::vector<double> m_timings, m_timingsMax;
  // coarsening data
  LO m_numNodesCoarse;
  std::vector<LO> m_nodeBeginCoarse, m_localDofsCoarse, m_coarseSizes;
  std::vector<GO> m_nodeGlobalIDsCoarse, m_nodeGlobalIDsCoarse1to1;
  std::vector<SM> m_xCoordCoarse, m_yCoordCoarse, m_zCoordCoarse;
  std::vector<SX> m_coarseRhsWork, m_coarseSolWork, m_fine1to1Work,
    m_coarseRhsWork2, m_residualVector, m_residVec2, m_residVec3, m_deltaSol;
  std::vector< std::vector<LO> > m_subNodesCoarse, m_subRowBeginCoarse, 
    m_subColumnsCoarse, m_subDofsCoarse;
  std::vector< std::vector<SX> > m_subValuesCoarse;
  RCP< PreconditionerBDDC<SX,SM,LO,GO> > m_coarsePreconditioner;
  GO m_numCoarseDofs;
  const std::string m_defaultSolver{"SuperLU"};
  int m_level, m_myCoarseMpiRank;
  OperatorBase<SX>* m_Operator{nullptr};     
  bool m_useDirectSolver{false}, m_fineProc{false}, m_coarseProc{false}, 
       m_disjointCommunicators{false};
  // data for disjoint communicators
  RCP<const Map> m_dofMapB_dis, m_dofMapB1to1_dis, m_dofMap_dis,
    m_dofMap1to1_dis;
  RCP<Export> m_exporterB_dis, m_exporterAll_dis;
  RCP<Vector> m_rhsVecB_dis, m_solVecB_dis, m_rhsVecB1to1_dis,
    m_solVecB1to1_dis, m_rhsVecAll_dis, m_solVecAll_dis,
    m_rhsVecAll1to1_dis, m_solVecAll1to1_dis;

  RCP< DisconnectedComponentChecker<LO,GO,SX> > m_componentChecker;
  // extra data for better communications
  std::vector<int> m_nodeSend, m_nodeSendCoarse;
  RCP<Tpetra::Distributor> m_distributor, m_distributorSubData,
    m_distributorCoarse;
  RCP< PComm<LO,GO,SX> > m_pComm, m_pCommB, m_pCommCoarse, m_pComm_dis, m_pCommB_dis;
  bool m_usePComm;
  // vertex coarse space data
#ifdef VERTEXCOARSESPACEBDDC
  RCP< VertexCoarseSpace<SX,SM,LO,GO> > m_vertexCoarseSpace, m_vertexCoarseSpaceIn;
#endif

  void checkVertexCoarseSpaceAvailability()
  {
    if (m_useVertexCoarseSpace) {
#ifndef VERTEXCOARSESPACEBDDC
      std::cout << "Error: Vertex coarse space is not available, turning off\n";
      m_useVertexCoarseSpace = false;
      m_Parameters->set("Use Vertex Coarse Space", false);
#endif
    }
  }

  bool needInitialStaticCondensationCorrection
    (const bool interfaceValuesOnly)
  {
    if (interfaceValuesOnly) return false;
    if (m_level > 0) return true;
    const std::string DirichletSolver = 
      m_Parameters->get("Dirichlet Solver", m_defaultSolver);
    if (DirichletSolver == "MueLu") return true;
    if (DirichletSolver == "NodalAMG") return true;
    if (m_Operator != nullptr) return true;
    return false;
  }

  void checkDirichletSolver()
  {
    std::string solverDirichlet = 
      m_Parameters->get("Dirichlet Solver", m_defaultSolver);
    if ((solverDirichlet == "MueLu") && m_interfacePreconditioner) {
      if (m_myPID == 0) {
	std::cout << "Error: cannot use an interface preconditioner without a direct solver\n";
	std::cout << " Dirichlet solver = " << solverDirichlet << std::endl;
      }
      BDDC_TEST_FOR_EXCEPTION(solverDirichlet == "MueLu", 
			      std::runtime_error, "MuelLu not allowed for Dirichlet solver");
    }
  }

  void initialize(std::vector<GO>* nodeGlobalIDs1to1,
		  std::vector<int>* nodeSend)
  {
    checkVertexCoarseSpaceAvailability();
    m_myPID = m_Comm->getRank();
    checkDirichletSolver();
    determineNodeSend(nodeSend);
    bddc::constructDistributor<LO>
      (m_numNodes, m_nodeSend, m_Comm, m_distributor);
    determineSubNodeBeginAndSubLocalDofs();
    initializeVariables();
    initializePreconditioner(nodeGlobalIDs1to1);
  }

  void determineNodeSend(std::vector<int>* nodeSend)
  {
    if (nodeSend == nullptr) {
      bddc::getNodeSend<LO,GO>
	(m_numNodes, m_nodeGlobalIDs, m_mpiComm, m_nodeSend);
    }
    else {
      m_nodeSend = *nodeSend;
    }
  }

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
    BDDC_TEST_FOR_EXCEPTION(numDofPerNode == 0, std::runtime_error, 
			    "numDofPerNode must be positive");
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

  LO getNumSubRowsMax()
  {
    LO numRowsMax = 0;
    const LO numSub = m_subNodesCoarse.size();
    for (LO i=0; i<numSub; i++) {
      const LO numRows = getNumSubRows(i);
      if (numRows > numRowsMax) numRowsMax = numRows;
    }
    return numRowsMax;
  }

  LO getNumSubRows(const LO sub)
  {
    LO numRows = 0;
    for (size_t i=0; i<m_subNodesCoarse[sub].size(); i++) {
      const LO node = m_subNodesCoarse[sub][i];
      numRows += m_nodeBeginCoarse[node+1] - m_nodeBeginCoarse[node];
    }
    return numRows;
  }

  void initializeVariables()
  {
    m_spatialDim = m_Parameters->get("Spatial Dimension", 3);
    m_problemType = m_Parameters->get("Problem Type BDDC", bddc::SCALARPDE);
    m_numThreads = m_Parameters->get<int>("numThreadsOuter", 1);
    m_IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    m_numDofs = m_nodeBegin[m_numNodes];
    m_numSub = m_subNodeBegin.size();
    m_numDofsB = 0;
    m_numDofsB1to1 = 0;
    m_work1T = nullptr;
    m_work2T = nullptr;
    m_myCoarseMpiRank = -1;
    m_useDirectSolver = false;
    m_fineProc = false;
    m_coarseProc = false;
    m_disjointCommunicators = false;
    m_numCoarseDofs = -1;
    // no deluxe scaling currently enabled when using PComm
    if (m_usePComm) {
      if (m_Parameters->get("Weight Type Corner", STIFFNESS) == DELUXE) {
	m_Parameters->set("Weight Type Corner", STIFFNESS);
      }
      if (m_Parameters->get("Weight Type Edge", STIFFNESS) == DELUXE) {
	m_Parameters->set("Weight Type Edge", STIFFNESS);
      }
      if (m_Parameters->get("Weight Type Face", STIFFNESS) == DELUXE) {
	m_Parameters->set("Weight Type Face", STIFFNESS);
      }
    }
  }

  void makeOwnedMapConsistent(RCP<const Map> dofMap,
			      const LO numNodes,
			      const LO* nodeBegin,
			      const std::vector<int> & nodeSend,
			      RCP<const Map> & dofMap1to1) const
  {
    const size_t numDof = nodeBegin[numNodes];
    BDDC_TEST_FOR_EXCEPTION(numDof != dofMap->getNodeNumElements(), 
			    std::runtime_error, "inconsistent numDof");
    LO row(0);
    std::vector<GO> globalIDs1to1;
    for (LO i=0; i<numNodes; i++) {
      for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	if (nodeSend[i] == -1) {
	  GO globalID = dofMap->getGlobalElement(row);
	  globalIDs1to1.push_back(globalID);
	}
	row++;
      }
    }
    dofMap1to1 = Teuchos::null;
    dofMap1to1 = 
      rcp( new Map(m_IGO, Teuchos::ArrayView<const GO>(globalIDs1to1), 
		   0, dofMap->getComm()) );
  }

  void initializePreconditioner(std::vector<GO>* nodeGlobalIDs1to1)
  {
    m_useDirectSolver = checkUseDirectSolver();
    if (m_useDirectSolver == true) {
      initializeDirectSolverAndVectors();
      return;
    }

    double startTime = GetTime();
    m_Partition =
      rcp( new bddc::PartitionOfUnity<SX,SM,LO,GO>
	   (m_numNodes, m_nodeGlobalIDs, m_subNodeBegin, m_subNodes, 
	    m_spatialDim, m_Parameters, m_Comm, m_distributor, 
	    m_nodeSend, m_xCoord, m_yCoord, m_zCoord) );
    m_timings[TIME_PARTITION_OF_UNITY] += GetTime() - startTime;
 
    if (m_usePComm == false) {
      startTime = GetTime();
      if ((m_dofMap == Teuchos::null) || (m_dofMap1to1 == Teuchos::null)) {
	bddc::DofManager<LO,GO>::
	  determineGlobalIDs(m_numNodes, m_nodeGlobalIDs, m_nodeBegin, 
			     m_localDofs, m_Comm, m_dofMap, m_dofMap1to1, 
			     nodeGlobalIDs1to1);
      }
      m_timings[TIME_DOF_MANAGER] += GetTime() - startTime;
      
      // temporily make m_dofMap1to1 consistent with m_nodeSend for testing
      // of new PComm class
      
      makeOwnedMapConsistent(m_dofMap, m_numNodes, m_nodeBegin, m_nodeSend, 
			     m_dofMap1to1);
    }

    m_pComm = rcp( new PComm<LO,GO,SX>(m_distributor, m_numNodes, m_nodeBegin,
				       m_nodeGlobalIDs, m_nodeSend.data()) );

    startTime = GetTime();
    determineBoundaryDofs();
    determineBoundaryMaps();
    initializeSubdomains();
    determineInteriorMaps();
    determineDiagBoundary();
    m_timings[TIME_INIT_STEP1] += GetTime() - startTime;

    determineBaseConstraints();
    determineWeights();
    determineAuxiliaryConstraints();
    determineCoarseSpace(m_coarseMpiRanks);
    factorInteriorMatrices();
    initializeVectors();
    initializeDisjointData();
    reserveMemory();
  }

  void retrieveCoarseSizesAcrossLevels()
  {
    if (m_coarseMpiRanks.size() > 0) {
      Teuchos::broadcast<int, GO>
	(*m_Comm, m_coarseMpiRanks[0], 1, &m_numCoarseDofs);
    }
    if (m_level == 0) {
      m_Comm->barrier();
      if (m_coarseProc) {
	m_coarseSizes.resize(0);
	LO numDofsVertex(0);
	if (m_useVertexCoarseSpace) {
#ifdef VERTEXCOARSESPACEBDDC
	  numDofsVertex = m_vertexCoarseSpace->getNumVertexDof();
#endif
	}
	m_coarsePreconditioner->appendNumCoarseDofs(m_coarseSizes, numDofsVertex);
      }
      LO numCoarseLevels = m_coarseSizes.size();
      LO numCoarseLevelsMax;
      Teuchos::reduceAll<int, LO>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				  &numCoarseLevels, &numCoarseLevelsMax);
      if (numCoarseLevelsMax == 0) return;
      LO procFlag(-1), procFlagMax;
      if (numCoarseLevels == numCoarseLevelsMax) procFlag = m_myPID;
      Teuchos::reduceAll<int, LO>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				  &procFlag, &procFlagMax);
      m_coarseSizes.resize(numCoarseLevelsMax);
      Teuchos::broadcast<int, LO> (*m_Comm, procFlagMax, numCoarseLevelsMax, 
				   m_coarseSizes.data());
    }
  }
  
  void retrieveTimingsAcrossLevels()
  {
    if (m_level == 0) {
      m_Comm->barrier();
      m_timingsMax.resize(0);
      appendTimings(m_timingsMax);
      LO numTimings = m_timingsMax.size();
      LO numTimingsMax;
      Teuchos::reduceAll<int, LO>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				  &numTimings, &numTimingsMax);
      LO procFlag(-1), procFlagMax;
      if (numTimings == numTimingsMax) procFlag = m_myPID;
      m_timingsMax.resize(numTimingsMax);
      Teuchos::reduceAll<int, LO>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				  &procFlag, &procFlagMax);
      Teuchos::broadcast<int, double> (*m_Comm, procFlagMax, numTimingsMax, 
					 m_timingsMax.data());
    }
  }
  
  void appendNumCoarseDofs(std::vector<LO> & coarseSizes,
			   const LO numDofsVertex)
  {
    LO numDofs = getNumDof();
    coarseSizes.push_back(numDofs);
    coarseSizes.push_back(numDofsVertex);
    LO numVertexDofs(0);
#ifdef VERTEXCOARSESPACEBDDC
    if (m_vertexCoarseSpace != Teuchos::null) {
      numVertexDofs = m_vertexCoarseSpace->getNumVertexDof();
    }
#endif
    if (m_coarsePreconditioner != Teuchos::null) {
      m_coarsePreconditioner->appendNumCoarseDofs(coarseSizes, numVertexDofs);
    }
  }

  void appendTimings(std::vector<double> & timings)
  {
    const LO numTimers = TIME_LENGTH;
    std::vector<double> timingsMax(numTimers);
    Teuchos::reduceAll<int, double>
      (*m_Comm, Teuchos::REDUCE_MAX, TIME_LENGTH, m_timings.data(), 
       timingsMax.data());
    for (LO i=0; i<numTimers; i++) {
      timings.push_back(timingsMax[i]);
    }
    if (m_coarseProc) {
      m_coarsePreconditioner->appendTimings(timings);
    }
  }

  void setTpetraObjects(RCP<Vector> & rhsVec1to1, 
			RCP<Vector> & solVec1to1, 
			RCP<Vector> & rhsVec, 
			RCP<Vector> & solVec, 
			RCP<Export> & exporter,
			std::vector< std::vector<LO> > & subDofs)
  {
    if (m_interfacePreconditioner) {
      subDofs = m_subBoundaryDofs;
      if (m_usePComm == false) {
	rhsVec1to1 = m_rhsVecB1to1;
	solVec1to1 = m_solVecB1to1;
	rhsVec = m_rhsVecB;
	solVec = m_solVecB;
	exporter = m_exporterB;
      }
    }
    else {
      subDofs = m_subDofs;
      if (m_usePComm == false) {
	rhsVec1to1 = m_rhsVecAll1to1;
	solVec1to1 = m_solVecAll1to1;
	rhsVec = m_rhsVecAll;
	solVec = m_solVecAll;
	exporter = m_exporterAll;
      }
    }
  }

  void setPComm(RCP< PComm<LO,GO,SX> > & pComm)
  {
    if (m_interfacePreconditioner == true) {
      pComm = m_pCommB;
    }
    else {
      pComm = m_pComm;
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
    BDDC_TEST_FOR_EXCEPTION(nthreads_obtained != m_numThreads, 
			    std::runtime_error, "invalid nthreads_obtained");
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
    if (m_fineProc == false) return;
    if (m_usePComm) {
      if (m_interfacePreconditioner) {
	m_pCommB_dis = rcp( new PComm<LO,GO,SX>(*m_pCommB, m_CommFine) );
      }
      else {
	m_pComm_dis = rcp( new PComm<LO,GO,SX>(*m_pComm, m_CommFine) );
      }
    }
    else {
      if (m_interfacePreconditioner) {
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
	m_rhsVecAll_dis = rcp( new Vector(m_dofMap_dis) );
	m_solVecAll_dis = rcp( new Vector(m_dofMap_dis) );
	m_rhsVecAll1to1_dis = rcp( new Vector(m_dofMap1to1_dis) );
	m_solVecAll1to1_dis = rcp( new Vector(m_dofMap1to1_dis) );
	m_exporterAll_dis = rcp( new Export(m_dofMap_dis, m_dofMap1to1_dis) );
      }
    }
  }

  void initializeVectors()
  {
    if (m_usePComm) return;
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
    GO numMySubs = m_subNodeBegin.size();
    GO numSubsAll;
    Teuchos::reduceAll<int, GO>(*m_Comm, Teuchos::REDUCE_SUM, 1,
				&numMySubs, &numSubsAll);
    BDDC_TEST_FOR_EXCEPTION(numSubsAll <= 0, std::runtime_error, 
			    "numSubsAll must be positive");
    bool useDirectSolver = false;
    if (numSubsAll == 1) {
      useDirectSolver = true;
    }
    m_numSub = numMySubs;
    return useDirectSolver;
  }

  void initializeDirectSolverAndVectors()
  {
    BDDC_TEST_FOR_EXCEPTION(m_Comm->getSize() != 1, std::runtime_error, 
		 "only single rank in communicator allowed for direct solver");
    BDDC_TEST_FOR_EXCEPTION(m_subNodeBegin.size() != 1, std::runtime_error, 
		 "only single subdomain allowed for direct solver");
    LO numRows = m_numDofs;
    if (numRows == 0) return;
    std::cout << "coarse space dimension = " << numRows << std::endl;
    LO* rowBegin = m_subRowBegin[0];
    LO* columns = m_subColumns[0];
    SX* values = m_subValues[0];

    if (m_Parameters->get("Print Coarsest Matrix", false)) {
      UtilBDDC<SX,SM>::printSparseMatrix
	(numRows, rowBegin, columns, values, "Ac.dat");
    }

    if (m_useVertexCoarseSpace) {
#ifdef VERTEXCOARSESPACEBDDC
      m_vertexCoarseSpaceIn->calculateReducedMatrix(numRows, rowBegin, columns, values);
      m_vertexCoarseSpaceIn->factorReducedMatrix(m_Parameters);
      m_vertexCoarseSpaceIn->initializeGaussSeidel(m_numNodes, m_nodeBegin);
#endif
    }
    else {
      SolverFactory<SX> Factory;
      m_coarseSub = rcp
	( new SubdomainBDDC<SX,SM,LO,GO>
	  (m_numNodes, m_subNodes[0].data(), m_nodeBegin, m_localDofs,
	   rowBegin, columns, values, m_xCoord, m_yCoord, m_zCoord,
	   *m_Parameters) );
    }
  }

  void determineWeights()
  {
    if (m_usePComm) {
      m_Weights = 
	rcp( new bddc::WeightsBDDC<SX,SM,LO,GO>
	     (m_Subdomain, m_Partition, m_Comm, m_pCommB,
	      m_subBoundaryDofs, m_diagBoundary, m_Parameters) );
    }
    else {
      m_Weights = 
	rcp( new bddc::WeightsBDDC<SX,SM,LO,GO>
	     (m_Subdomain, m_Partition, m_exporterB,
	      m_subBoundaryDofs, m_diagBoundary, m_Parameters) );
    }
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

  void determineBoundaryNodeData(std::vector<LO> & nodeBeginB, 
				 std::vector<GO> & nodeGlobalIDsB,
				 std::vector<int> & nodeSendB)
  {
    std::vector<LO> dofNodes(m_numDofs, -1);
    for (LO i=0; i<m_numNodes; i++) {
      for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	dofNodes[j] = i;
      }
    }
    LO prevNode(-1), numDofB(0);
    nodeBeginB.push_back(0);
    for (size_t i=0; i<m_boundaryDofs.size(); i++) {
      const LO node = dofNodes[m_boundaryDofs[i]];
      if (node != prevNode) {
	prevNode = node;
	nodeGlobalIDsB.push_back(m_nodeGlobalIDs[node]);
	nodeSendB.push_back(m_nodeSend[node]);
	numDofB += m_nodeBegin[node+1] - m_nodeBegin[node];
	nodeBeginB.push_back(numDofB);
      }
    }
    BDDC_TEST_FOR_EXCEPTION(numDofB != m_numDofsB, std::runtime_error, 
			    "inconsistent value of numDofB");
  }

  void determineBoundaryMaps()
  {
    std::vector<LO> nodeBeginB;
    std::vector<GO> nodeGlobalIDsB;
    std::vector<int> nodeSendB;
    determineBoundaryNodeData(nodeBeginB, nodeGlobalIDsB, nodeSendB);
    const LO numNodeB(nodeSendB.size());
    m_pCommB = rcp
      ( new PComm<LO,GO,SX>(m_distributor, numNodeB, nodeBeginB.data(),
			    nodeGlobalIDsB.data(), nodeSendB.data()) );
    m_numDofsB1to1 = m_pCommB->getOwnedLength();
    m_boundaryToAll1to1.resize(m_numDofsB1to1);
    std::map<GO,LO> nodeMap;
    const std::vector<GO> nodeGlobalIDsOwned = m_pComm->getOwnedNodesGlobalIDs();
    const LO numNodesOwned = nodeGlobalIDsOwned.size();
    for (LO i=0; i<numNodesOwned; i++) {
      nodeMap.emplace(nodeGlobalIDsOwned[i], i);
    }
    const std::vector<GO> nodeGlobalIDsOwnedB = m_pCommB->getOwnedNodesGlobalIDs();
    const LO numNodesOwnedB = nodeGlobalIDsOwnedB.size();
    const std::vector<LO> & nodeBeginOwned = m_pComm->getNodeBeginOwned();
    const std::vector<LO> & nodeBeginOwnedB = m_pCommB->getNodeBeginOwned();
    for (LO i=0; i<numNodesOwnedB; i++) {
      auto iter = nodeMap.find(nodeGlobalIDsOwnedB[i]);
      BDDC_TEST_FOR_EXCEPTION(iter == nodeMap.end(), std::runtime_error, 
			      "nodeGlobalIDsOwndedB[i] not found");
      const LO node = iter->second;
      for (LO j=nodeBeginOwnedB[i]; j<nodeBeginOwnedB[i+1]; j++) {
	m_boundaryToAll1to1[j] = nodeBeginOwned[node] + j - nodeBeginOwnedB[i];
      }
    }
    if (m_usePComm == false) {
      std::vector<GO> globalIDs(m_numDofsB);
      for (LO i=0; i<m_numDofsB; i++) {
	LO dof = m_boundaryDofs[i];
	globalIDs[i] = m_dofMap->getGlobalElement(dof);
      }
      m_dofMapB = generateMap(globalIDs);
      m_dofMapB1to1 = constructOneToOne(m_dofMapB, m_dofMap1to1);
      m_exporterB = rcp( new Export(m_dofMapB, m_dofMapB1to1) );
      m_numDofsB1to1 = m_dofMapB1to1->getNodeNumElements();
      m_exporterAll = rcp( new Export(m_dofMap, m_dofMap1to1) );
      m_boundaryToAll1to1.resize(m_numDofsB1to1);
      for (LO i=0; i<m_numDofsB1to1; i++) {
	GO globalID = m_dofMapB1to1->getGlobalElement(i);
	LO localID = m_dofMap1to1->getLocalElement(globalID);
	const bool test = (localID != Teuchos::OrdinalTraits<LO>::invalid());
	BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			      "invalid localID");
	m_boundaryToAll1to1[i] = localID;
      }
    }
  }

  void determineInteriorMaps()
  {
    std::vector<LO> dofMap;
    if (m_usePComm) {
      dofMap.resize(m_numDofs, -1);
      LO numDofOwned(0);
      for (LO i=0; i<m_numNodes; i++) {
	if (m_nodeSend[i] == -1) {
	  for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	    dofMap[j] = numDofOwned++;
	  }
	} 
      }
      BDDC_TEST_FOR_EXCEPTION(numDofOwned != m_pComm->getOwnedLength(), 
			      std::runtime_error, "invalid numDofOwned");
    }
    m_subInteriorDofs.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      const std::vector<LO> & interiorDofs = m_Subdomain[i]->getInteriorDofs();
      LO numInteriorDofs = interiorDofs.size();
      m_subInteriorDofs[i].resize(numInteriorDofs);
      for (LO j=0; j<numInteriorDofs; j++) {
	LO localID = m_subDofs[i][interiorDofs[j]];
	LO localID1to1;
	if (m_usePComm) {
	  localID1to1 = dofMap[localID];
	  BDDC_TEST_FOR_EXCEPTION(localID1to1 == -1, std::runtime_error, 
				  "invalid localID1to1");
	}
	else {
	  GO globalID = m_dofMap->getGlobalElement(localID);
	  localID1to1 = m_dofMap1to1->getLocalElement(globalID);
	  const bool test = 
	    localID1to1 != Teuchos::OrdinalTraits<LO>::invalid();
	  BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
				  "invalid localID1to1");
	}
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
    BDDC_TEST_FOR_EXCEPTION(m_numSub != LO(subdomainEquivClasses.size()), 
			    std::runtime_error, "inconsistency with m_numSub");
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
	BDDC_TEST_FOR_EXCEPTION(dofB == -1, std::runtime_error, 
				"invalid dofB");
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
	  BDDC_TEST_FOR_EXCEPTION(dofB == -1, std::runtime_error, 
				  "invalid dofB");
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
    int activeProc(0), activeProcSum;
    if (m_numSub > 0) activeProc = 1;
    Teuchos::reduceAll<int, int> 
      (*m_Comm, Teuchos::REDUCE_SUM, 1, &activeProc, &activeProcSum);
    int averageNumSubdomainsPerMpiRank = totalNumSubdomains/activeProcSum;
    numSubdomainsPerMpiRank = std::max(numSubdomainsPerMpiRank,
				       averageNumSubdomainsPerMpiRank);
    int numMpiRanks = totalNumSubdomains/numSubdomainsPerMpiRank;
    if (totalNumSubdomains % numSubdomainsPerMpiRank != 0) numMpiRanks++;
    numMpiRanks = std::min(numMpiRanks, m_Comm->getSize());
    determineCoarseProcs(numMpiRanks, myCoarseMpiRank, mpiRanks);
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

  void determineCoarseProcs(const int numMpiRanks,
			    int & myCoarseMpiRank,
			    std::vector<int> & mpiRanks)
  {
    mpiRanks.resize(numMpiRanks);
    SM procLoad = getProcessorLoad();
    int numProc = m_Comm->getSize();
    int numCoarseProc = mpiRanks.size();
    const int root(0);
    std::vector<SM> procLoadAll;
    if (m_myPID == root) procLoadAll.resize(numProc);
    // gather processor loads
    Teuchos::gather<int, SM> (&procLoad, 1, procLoadAll.data(), numProc,
			      root, *m_Comm);
    if (m_myPID == root) {
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
      // the following sort is needed for consistency with m_myCoarseMpiRank
      std::sort(mpiRanks.begin(), mpiRanks.end());
    }
    Teuchos::broadcast<int, int> (*m_Comm, root, numCoarseProc, mpiRanks.data());
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

  void determineElemCoords(const std::vector<SubdomainData> & coarseSubData, 
			   std::vector<double> & elemCoords)
  {
    const LO numElems = coarseSubData.size();
    elemCoords.resize(3*numElems);
    for (LO i=0; i<numElems; i++) {
      const SubdomainData & subData = coarseSubData[i];
      double sumX(0), sumY(0), sumZ(0);
      const LO numNodes = subData.numNodes;
      BDDC_TEST_FOR_EXCEPTION(numNodes == 0, std::runtime_error, 
			      "numNodes must be positive");
      for (LO j=0; j<subData.numNodes; j++) {
	sumX += subData.coords[3*j+0];
	sumY += subData.coords[3*j+1];
	sumZ += subData.coords[3*j+2];
      }
      elemCoords[i+0*numElems] = sumX/numNodes;
      elemCoords[i+1*numElems] = sumY/numNodes;
      elemCoords[i+2*numElems] = sumZ/numNodes;
    }
  }

  void determineElemConnectivity(const std::vector<SubdomainData> & coarseSubData,
				 const LO numFaces, 
				 std::vector<LO> & rowBegin, 
				 std::vector<GO> & columns)
  {
    const LO numElems = coarseSubData.size();
    std::vector< std::vector<LO> > faceElems(numFaces);
    for (LO i=0; i<numElems; i++) {
      const SubdomainData & subData = coarseSubData[i];
      for (LO j=0; j<subData.numFaces; j++) {
	const LO face = subData.faceLIDs[j];
	faceElems[face].push_back(i);
      }
    }
    std::vector< std::vector<LO> > elemConn(numElems);
    LO numTerms = 0;
    for (LO i=0; i<numFaces; i++) {
      BDDC_TEST_FOR_EXCEPTION(faceElems[i].size() >= 3, std::runtime_error, 
			      "faceElems[i] size must be less than 3");
      if (faceElems[i].size() == 2) {
	const LO elem1 = faceElems[i][0];
	const LO elem2 = faceElems[i][1];
	elemConn[elem1].push_back(elem2);
	elemConn[elem2].push_back(elem1);
	numTerms += 2;
      }
    }
    rowBegin.resize(numElems+1, 0);
    columns.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numElems; i++) {
      for (size_t j=0; j<elemConn[i].size(); j++) {
	columns[numTerms++] = elemConn[i][j];
      }
    }
  }

  void partitionSubdomainsWithinRank
    (const std::vector<SubdomainData> & coarseSubData,
     const LO numFaces,
     std::vector< std::vector<LO> > & coarseSubdomainSubs)
  {
    LO numSubdomains = coarseSubData.size();
    int numSubdomainsPerCoarseSubdomain = 
      m_Parameters->get("numSubdomainsPerCoarseSubdomain", 512);
    int numParts = numSubdomains/numSubdomainsPerCoarseSubdomain;
    if (numSubdomains % numSubdomainsPerCoarseSubdomain != 0) numParts++;
    // quick exit if no further partitioning required
    if (numParts <= 1) {
      if (numParts == 1) {
	coarseSubdomainSubs.resize(1);
	coarseSubdomainSubs[0].resize(numSubdomains);
	for (LO i=0; i<numSubdomains; i++) coarseSubdomainSubs[0][i] = i;
      }
      return;
    }
    std::string coarsenOption = m_Parameters->get("Coarsening Option", "Graph");
    std::vector<double> elemCoords;
    std::vector<LO> rowBegin;
    std::vector<GO> columnGIDs, columnPIDs;
    determineElemCoords(coarseSubData, elemCoords);
    if (coarsenOption == "Graph") {
      columnPIDs.resize(numSubdomains, 0);
      determineElemConnectivity(coarseSubData, numFaces, rowBegin, columnGIDs);
    }
    std::vector<LO> parts(numSubdomains, 0);
    RCP<Teuchos::ParameterList> params = rcp( new Teuchos::ParameterList() );
    params->set("Number of Parts", numParts);
    ZoltanPartition<LO, GO>* partition(nullptr);
    params->set("Coordinates", elemCoords.data());
    std::vector<GO> globalIDs(numSubdomains);
    for (LO i=0; i<numSubdomains; i++) globalIDs[i] = i;
    MPI_Comm mpiComm = MPI_COMM_SELF;
    if (coarsenOption == "Graph") {
      params->set("LB Method", "Graph");
    }
    else {
      if (coarsenOption == "RCB") {
	params->set("LB Method", "RCB");
      }
      else if (coarsenOption == "RIB") {
	params->set("LB Method", "RIB");
      }
    }
    partition = new ZoltanPartition<LO, GO>
      (numSubdomains, rowBegin.data(), columnGIDs.data(), columnPIDs.data(),
       elemCoords.data(), globalIDs.data(), mpiComm, params);
    partition->setDefaultWeights();
    partition->doPartition();
    const int* partNumbers = partition->getParts();
    coarseSubdomainSubs.resize(numParts);
    for (int i=0; i<numSubdomains; i++) {
      const int part = partNumbers[i];
      coarseSubdomainSubs[part].push_back(i);
    }
    delete partition;
  }

  void partitionSubdomainsAcrossMpiRanks(const std::vector<int> & mpiRanks,
					 std::vector<GO> & parts)
  {
    // Note: parts[i] is the coarse subdomain number for local subdomain i
    //       after coarsening
    parts.resize(m_numSub, 0);
    std::vector<double> subdomainCoords;
    determineSubdomainCoords(subdomainCoords);
    if (mpiRanks.size() > 1) {
      RCP<Teuchos::ParameterList> params = rcp( new Teuchos::ParameterList() );
      int numCoarseRanks = mpiRanks.size();
      params->set("Number of Parts", numCoarseRanks);
      params->set("Coordinates", subdomainCoords.data());
      ZoltanPartition<LO, GO>* partition(nullptr);
      const GO *subGIDs(0), *subConnGIDs(0), *subConnProcs(0);
      const LO *subConnBegin(0);
      const std::string coarsenOption = 
	m_Parameters->get("Coarsening Option", "Graph");
      if (coarsenOption == "Graph") {
	m_Partition->getConnectivityGraph
	  (subGIDs, subConnGIDs, subConnProcs, subConnBegin);
	//	RCP<const CrsGraph> Graph = m_Partition->getConnectivityGraph();
	params->set("LB Method", "Graph");
	const std::string graphPackage = 
	  m_Parameters->get("Graph Package", "PHG");
	params->set("Graph Package", graphPackage);
	partition = new ZoltanPartition<LO, GO>
	  (m_numSub, subConnBegin, subConnGIDs, subConnProcs, 
	   subdomainCoords.data(), subGIDs, m_mpiComm, params);
      }
      else {
	std::vector<GO> globalIDs(m_numSub);
	int startingSub = m_Partition->getStartingSub();
	for (int i=0; i<m_numSub; i++) globalIDs[i] = startingSub + i;
	if (coarsenOption == "RCB") {
	  params->set("LB Method", "RCB");
	}
	else if (coarsenOption == "RIB") {
	  params->set("LB Method", "RIB");
	}
	partition = new ZoltanPartition<LO, GO>
	  (m_numSub, subConnBegin, subConnGIDs, subConnProcs,
	   subdomainCoords.data(), globalIDs.data(), m_mpiComm, params);
      }
      partition->setDefaultWeights();
      partition->doPartition();
      const int* partNumbers = partition->getParts();
      for (int i=0; i<m_numSub; i++) parts[i] = partNumbers[i];
      delete partition;
    }
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

  void determineCoarseSubdomainDofs(const std::vector<LO> & activeCoarseNodes,
				    const std::vector<LO> & nodeBegin,
				    std::vector< std::vector<LO> > & subDofsCoarse) const
  {
    LO numCoarseNodes = m_Partition->getNumEquivClasses();
    std::vector<LO> nodeMap(numCoarseNodes, -1);
    const LO numActiveCoarseNodes = activeCoarseNodes.size();
    for (LO i=0; i<numActiveCoarseNodes; i++) {
      const LO node = activeCoarseNodes[i];
      nodeMap[node] = i;
    }
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    subDofsCoarse.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      const LO numEquiv = m_Subdomain[i]->getNumEquiv();
      LO numDofSub(0);
      for (LO j=0; j<numEquiv; j++) {
	const LO equiv = subdomainEquivClasses[i][j];
	const LO node = nodeMap[equiv];
	if (node != -1) numDofSub += nodeBegin[node+1] - nodeBegin[node];
      }
      subDofsCoarse[i].resize(numDofSub);
      numDofSub = 0;
      for (LO j=0; j<numEquiv; j++) {
	const LO equiv = subdomainEquivClasses[i][j];
	const LO node = nodeMap[equiv];
	if (node != -1) {
	  for (LO k=nodeBegin[node]; k<nodeBegin[node+1]; k++) {
	    subDofsCoarse[i][numDofSub++] = k;
	  }
	}
      }
    }
  }

  void determineActiveCoarseNodeData(std::vector<LO> & activeCoarseNodes,
				     std::vector<GO> & activeCoarseNodeGIDs,
				     std::vector<SM> & coarseCoords, 
				     std::vector<LO> & nodeBegin, 
				     std::vector<LO> & localDofs)
  {
    const std::vector<GO> coarseNodeGIDs = m_Partition->getGlobalIDs();
    std::vector< std::vector<LO> > coarseNodeLocalDofs;
    getCoarseNodeLocalDofs(coarseNodeLocalDofs);
    BDDC_TEST_FOR_EXCEPTION(coarseNodeLocalDofs.size() != coarseNodeGIDs.size(),
		  std::runtime_error, "invalid coarseNodeLocalDofs size");
    LO numCoarseNodes = coarseNodeGIDs.size();
    LO numActiveCoarseNodes(0), numActiveCoarseDofs(0);
    for (LO i=0; i<numCoarseNodes; i++) {
      /*
      if (coarseNodeGIDs[i] == 982) {
	std::cout << "number of dofs for coarseNodeGIDs[" << i << "] = 982 is " 
		  << coarseNodeLocalDofs[i].size() << std::endl;
      }
      */
      if (coarseNodeLocalDofs[i].size() > 0) {
	numActiveCoarseNodes++;
	numActiveCoarseDofs += coarseNodeLocalDofs[i].size();
      }
    }
    activeCoarseNodes.resize(numActiveCoarseNodes);
    activeCoarseNodeGIDs.resize(numActiveCoarseNodes);
    nodeBegin.resize(numActiveCoarseNodes+1, 0);
    localDofs.resize(numActiveCoarseDofs);
    numActiveCoarseNodes = numActiveCoarseDofs = 0;
    for (LO i=0; i<numCoarseNodes; i++) {
      if (coarseNodeLocalDofs[i].size() > 0) {
	activeCoarseNodes[numActiveCoarseNodes] = i;
	activeCoarseNodeGIDs[numActiveCoarseNodes++] = coarseNodeGIDs[i];
	for (size_t j=0; j<coarseNodeLocalDofs[i].size(); j++) {
	  localDofs[numActiveCoarseDofs++] = coarseNodeLocalDofs[i][j];
	}
	nodeBegin[numActiveCoarseNodes] = numActiveCoarseDofs;
      }
    }
    getCoordsCoarseNode(activeCoarseNodes, coarseCoords);
  }

  RCP<const Map> generateMap(const std::vector<GO> & globalIDs)
  {
    return
      rcp( new Map(m_IGO, Teuchos::ArrayView<const GO>(globalIDs), 0, m_Comm) );
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
    if (m_usePComm == false) {
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
    }
    LO numCoarseDof = m_pCommCoarse->getTargetLength();
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
    std::vector< std::vector<LO> > dummyVec;
    std::vector< std::vector<LO> > & subDofs = dummyVec;
    if (restrictToBoundary == true) {
      fineVector = m_fineVecB;
      exporter = m_exporterB;
      subDofs = m_subBoundaryDofs;
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

  void importCoarseSolution(std::vector<SX> & coarseSolOwned,
			    SX* & coarseSol)
  {
    coarseSol = m_pCommCoarse->getSourcePtr();
    if (m_usePComm) {
      m_pCommCoarse->doImport(coarseSolOwned.data(), coarseSol);
    }
    else {
      const LO numCoarseOwned = m_coarseVecNextLevel->getLocalLength();
      bool test = (numCoarseOwned == m_pCommCoarse->getTargetLength());
      BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			      "invalid numCoarseOwned");
      Teuchos::ArrayRCP<SX> valsOwned = m_coarseVecNextLevel->getDataNonConst();
      for (LO i=0; i<numCoarseOwned; i++) valsOwned[i] = coarseSolOwned[i];
      m_coarseVec->doImport(*m_coarseVecNextLevel, *m_exporterCoarse,
			    Tpetra::INSERT);
      Teuchos::ArrayRCP<const SX> vals = m_coarseVec->getData();
      const LO numCoarse = m_coarseVec->getLocalLength();
      test = (numCoarse == m_pCommCoarse->getSourceLength());
      BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			      "invalid numCoarse");
      for (LO i=0; i<numCoarse; i++) coarseSol[i] = vals[i];
    }
  }

  void addCoarseCorrection(const SX* coarseSol,
			   SX* currentSol)
  {
    for (LO i=0; i<m_numSub; i++) {
      SX *x(0), *b(0);
      m_Subdomain[i]->getSubdomainVectors(x, b);
      const std::vector<LO> & subDofsCoarse = m_subDofsCoarse[i];
      LO numRows = subDofsCoarse.size();
      for (LO j=0; j<numRows; j++) {
	const LO row = subDofsCoarse[j];
	x[j] = coarseSol[row];
      }
      const bool transpose(false);
      const bool interfacePreconditioner(false);
      m_Subdomain[i]->multiplyByPhi(x, b, interfacePreconditioner, transpose);
      const std::vector<LO> & dofs = m_subDofs[i];
      LO numDofs = dofs.size();
      for (LO j=0; j<numDofs; j++) {
	const LO row = dofs[j];
	currentSol[row] += b[j];
      }
    }
  }

  void applyPhiTranspose(const SX* rhsAll, 
			 std::vector<SX> & coarseRhsTarget)
  {
    SX* bVals = m_pCommCoarse->getSourcePtr();
    const LO lengthSource = m_pCommCoarse->getSourceLength();
    memset(bVals, 0, lengthSource*sizeof(SX));
    for (LO i=0; i<m_numSub; i++) {
      SX *x(0), *b(0);
      m_Subdomain[i]->getSubdomainVectors(x, b);
      const std::vector<LO> & dofs = m_subDofs[i];
      LO numDofs = dofs.size();
      for (LO j=0; j<numDofs; j++) {
	x[j] = rhsAll[dofs[j]];
      }
      const bool transpose(true);
      const bool interfacePreconditioner(false);
      m_Subdomain[i]->multiplyByPhi(x, b, interfacePreconditioner, transpose);
      const std::vector<LO> & subDofsCoarse = m_subDofsCoarse[i];
      LO numRows = subDofsCoarse.size();
      LO numDof(0);
      for (LO j=0; j<numRows; j++) {
	LO row = subDofsCoarse[j];
	bVals[row] += b[numDof++];
      }
    }
    const LO lengthTarget = m_pCommCoarse->getTargetLength();
    memset(coarseRhsTarget.data(), 0, lengthTarget*sizeof(SX));
    if (m_usePComm) {
      m_pCommCoarse->doExport(bVals, coarseRhsTarget.data());
    }
    else {
      SX* coarseVecVals = m_coarseVec->getDataNonConst().get();
      for (LO i=0; i<lengthSource; i++) {
	coarseVecVals[i] = bVals[i];
      }
      m_coarseVecNextLevel->putScalar(0);
      m_coarseVecNextLevel->doExport(*m_coarseVec, *m_exporterCoarse,
				    Tpetra::ADD);
      const LO numCoarseDof = m_coarseVecNextLevel->getLocalLength();
      BDDC_TEST_FOR_EXCEPTION(numCoarseDof != lengthTarget, std::runtime_error, 
			      "numCoarseDof not equal to lengthTarget");
      Teuchos::ArrayRCP<const SX> rhsNextLev = m_coarseVecNextLevel->getData();
      for (LO i=0; i<numCoarseDof; i++) coarseRhsTarget[i] = rhsNextLev[i];
    }
  }

  void printMatrices
    (const std::vector< std::vector<LO> > & rowBeginIn,
     const std::vector< std::vector<LO> > & columnsIn,
     const std::vector< std::vector<SX> > & valuesIn,
     const char* baseName)
  {
    LO numMatrices = rowBeginIn.size();
    for (LO i=0; i<numMatrices; i++) {
      char fname[101];
      sprintf(fname, "%s_%d_%d.dat", baseName, m_myPID, i);
      int numDof = rowBeginIn[i].size() - 1;
      const int* rowBegin = rowBeginIn[i].data();
      const int* columns = columnsIn[i].data();
      const SX* values = valuesIn[i].data();
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
    BDDC_TEST_FOR_EXCEPTION(currentCol != numCoarseDof, std::runtime_error, 
			      "counter not numCoarseDof");
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

  void determineSubsForRank(const std::vector<GO> & partsAcrossMpiRanks,
			    const std::vector<LO> & mpiRanks,
			    std::vector<int> & targetRanks, 
			    std::vector< std::vector<LO> > & subsForRank)
  {
    std::map<GO, LO> rankMap;
    const LO numSub = partsAcrossMpiRanks.size();
    LO numRanks(0), numSubsForMyself(0);
    for (LO i=0; i<numSub; i++) {
      const GO rank = mpiRanks[partsAcrossMpiRanks[i]];
      // do not include m_myPID in targetRanks
      if (rank == m_myPID) {
	numSubsForMyself++;
	continue;
      }
      auto iter = rankMap.find(rank);
      if (iter == rankMap.end()) {
	rankMap.insert(std::make_pair(rank, numRanks));
	targetRanks.push_back(rank);
	numRanks++;
      }
    }
    // include subdomains that include me as a target
    if (numSubsForMyself > 0) {
      rankMap.insert(std::make_pair(m_myPID, numRanks++));
    }
    subsForRank.resize(numRanks);
    for (LO i=0; i<numSub; i++) {
      const GO rank = mpiRanks[partsAcrossMpiRanks[i]];
      auto iter = rankMap.find(rank);
      BDDC_TEST_FOR_EXCEPTION(iter == rankMap.end(), std::runtime_error, 
			      "rank not found");
      const LO localSub = iter->second;
      subsForRank[localSub].push_back(i);
    }
  }

  void getNumNodesRowsAndFaces(const std::map<LO, LO> & nodeMap, 
			       const std::vector<LO> & equivClasses, 
			       const std::vector<LO> & nodeBegin,
			       LO & numNodes, 
			       LO & numRows,
			       LO & numFaces)
  {
    numNodes = numRows = numFaces = 0;
    const std::vector<LO> & cardinality = m_Partition()->getEquivCardinality();
    for (size_t i=0; i<equivClasses.size(); i++) {
      const LO equiv = equivClasses[i];
      if (cardinality[equiv] == 2) numFaces++;
      auto iter = nodeMap.find(equiv);
      if (iter != nodeMap.end()) {
	const LO node = iter->second;
	numNodes++;
	numRows += nodeBegin[node+1] - nodeBegin[node];
      }
    }
  }

  void determinePointersForDataToSend
    (const std::vector<LO> & activeCoarseNodes, 
     const std::vector<LO> & nodeBegin, 
     const std::vector< std::vector<LO> > & subsForRank,
     std::vector<LO> & startDataLO,
     std::vector<LO> & startDataGO,
     std::vector<LO> & startDataSM,
     std::vector<LO> & startDataSX)
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    std::map<LO, LO> nodeMap;
    const LO numActive = activeCoarseNodes.size();
    for (LO i=0; i<numActive; i++) {
      nodeMap.insert(std::make_pair(activeCoarseNodes[i], i));
    }
    const LO numRanks = subsForRank.size();
    startDataLO.resize(numRanks+1, 0);
    startDataGO.resize(numRanks+1, 0);
    startDataSM.resize(numRanks+1, 0);
    startDataSX.resize(numRanks+1, 0);
    for (LO i=0; i<numRanks; i++) {
      LO lengthLO(0), lengthGO(0), lengthSM(0), lengthSX(0);
      for (size_t j=0; j<subsForRank[i].size(); j++) {
	const LO sub = subsForRank[i][j];
	const std::vector<LO> & equivClasses = subdomainEquivClasses[sub];
	LO numNodes(0), numRows(0), numFaces(0);
	getNumNodesRowsAndFaces(nodeMap, equivClasses, nodeBegin, numNodes, 
				numRows, numFaces);
	// LO data has numNodes, numRows, numFaces, nodeBegin, localDofs,
	// and nodeSend
	lengthLO += 1 + 1 + 1 + (numNodes+1) + numRows + numNodes;
	// GO data has globalIDs of nodes and faces
	lengthGO += numNodes + numFaces;
	// SM data has nodal coordinates
	lengthSM += 3*numNodes;
	// SX data has coarse element matrices
	lengthSX += numRows*numRows;
      }
      startDataLO[i+1] = startDataLO[i] + lengthLO;
      startDataGO[i+1] = startDataGO[i] + lengthGO;
      startDataSM[i+1] = startDataSM[i] + lengthSM;
      startDataSX[i+1] = startDataSX[i] + lengthSX;
    }
  }

  void packCoarseData
    (const std::vector<LO> & activeCoarseNodes, 
     const std::vector<GO> & activeCoarseNodeGIDs,
     const std::vector<SM> & coarseCoords, 
     const std::vector<LO> & nodeBegin, 
     const std::vector<LO> & localDofs,
     const std::vector< std::vector<LO> > & subsForRank,
     const std::vector<int> & targetRanks,
     const std::vector<LO> & startDataLO,
     const std::vector<LO> & startDataGO,
     const std::vector<LO> & startDataSM,
     const std::vector<LO> & startDataSX,
     const std::vector<int> & minPartForActiveCoarseNodes,
     std::vector<LO> & dataLO,
     std::vector<GO> & dataGO,
     std::vector<SM> & dataSM,
     std::vector<SX> & dataSX)
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    std::map<LO, LO> nodeMap;
    const LO numActive = activeCoarseNodes.size();
    for (LO i=0; i<numActive; i++) {
      nodeMap.insert(std::make_pair(activeCoarseNodes[i], i));
    }
    const LO numRanks = subsForRank.size();
    dataLO.resize(startDataLO[numRanks]);
    dataGO.resize(startDataGO[numRanks]);
    dataSM.resize(startDataSM[numRanks]);
    dataSX.resize(startDataSX[numRanks]);
    LO countLO(0), countGO(0), countSM(0), countSX(0);
    const std::vector<LO> & cardinality = m_Partition()->getEquivCardinality();
    const std::vector<GO> & equivGIDs = m_Partition()->getGlobalIDs();
    for (LO i=0; i<numRanks; i++) {
      for (size_t j=0; j<subsForRank[i].size(); j++) {
	const LO sub = subsForRank[i][j];
	const std::vector<LO> & equivClasses = subdomainEquivClasses[sub];
	LO numNodes(0), numRows(0), numFaces(0);
	getNumNodesRowsAndFaces(nodeMap, equivClasses, nodeBegin, numNodes, 
				numRows, numFaces);
	dataLO[countLO++] = numNodes;
	dataLO[countLO++] = numRows;
	dataLO[countLO++] = numFaces;
	LO numDof = 0;
	// nodeBegin followed by localDofs
	dataLO[countLO++] = numDof;
	LO startLocalDof = countLO + numNodes;
	LO startNodeSend = startLocalDof + numRows;
	for (size_t k=0; k<equivClasses.size(); k++) {
	  const LO equiv = equivClasses[k];
	  auto iter = nodeMap.find(equiv);
	  if (iter != nodeMap.end()) {
	    const LO node = iter->second;
	    numDof += nodeBegin[node+1] - nodeBegin[node];
	    dataLO[countLO++] = numDof; // nodeBegin
	    dataLO[startNodeSend++] = minPartForActiveCoarseNodes[node];
	    for (LO m=nodeBegin[node]; m<nodeBegin[node+1]; m++) {
	      dataLO[startLocalDof++] = localDofs[m];
	    }
	    dataGO[countGO++] = activeCoarseNodeGIDs[node];
	    dataSM[countSM++] = coarseCoords[node + 0*numActive];
	    dataSM[countSM++] = coarseCoords[node + 1*numActive];
	    dataSM[countSM++] = coarseCoords[node + 2*numActive];
	  }
	}
	countLO += numRows + numNodes;
	for (size_t k=0; k<equivClasses.size(); k++) {
	  const LO equiv = equivClasses[k];
	  if (cardinality[equiv] == 2) {
	    dataGO[countGO++] = equivGIDs[equiv];
	  }
	}
	bool restrictPhiToBoundary(true), scaleRows(true);
	if (m_interfacePreconditioner == false) restrictPhiToBoundary = false;
	std::vector<SX> Ac;
	m_Subdomain[sub]->calculateCoarseMatrices
	  (Ac, restrictPhiToBoundary, scaleRows);
	const LO numTerms = Ac.size();
	BDDC_TEST_FOR_EXCEPTION(numTerms != numRows*numRows, 
				std::runtime_error, "Ac size error");
	for (LO k=0; k<numTerms; k++) dataSX[countSX++] = Ac[k];
      }
    }
  }

  void sendCoarseData(const std::vector< std::vector<LO> > & subsForRank,
		      const std::vector<LO> & startDataLO,
		      const std::vector<LO> & startDataGO,
		      const std::vector<LO> & startDataSM,
		      const std::vector<LO> & startDataSX,
		      const std::vector<LO> & dataLO,
		      const std::vector<GO> & dataGO,
		      const std::vector<SM> & dataSM,
		      const std::vector<SX> & dataSX,
		      std::vector<LO> & dataReceiveLO, 
		      std::vector<GO> & dataReceiveGO, 
		      std::vector<SM> & dataReceiveSM, 
		      std::vector<SX> & dataReceiveSX,
		      LO & numElems,
		      LO & numReceives,
		      std::vector<LO> & subReceivePtr)
  {
    // communicate coarse data sizes
    const LO numSends = m_distributorSubData->getNumSends();
    std::vector<LO> arraySendLO(5*numSends);
    std::vector<size_t> sendLengthLO(numSends), sendLengthGO(numSends),
      sendLengthSM(numSends), sendLengthSX(numSends);
    for (LO i=0; i<numSends; i++) {
      arraySendLO[5*i+0] = subsForRank[i].size();
      arraySendLO[5*i+1] = startDataLO[i+1] - startDataLO[i];
      arraySendLO[5*i+2] = startDataGO[i+1] - startDataGO[i];
      arraySendLO[5*i+3] = startDataSM[i+1] - startDataSM[i];
      arraySendLO[5*i+4] = startDataSX[i+1] - startDataSX[i];
      sendLengthLO[i] = startDataLO[i+1] - startDataLO[i];
      sendLengthGO[i] = startDataGO[i+1] - startDataGO[i];
      sendLengthSM[i] = startDataSM[i+1] - startDataSM[i];
      sendLengthSX[i] = startDataSX[i+1] - startDataSX[i];
    }
    numReceives = m_distributorSubData->getNumReceives();
    const bool test = (m_distributorSubData->hasSelfMessage() == false);
    BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			    "distributor hasSelfMessage error");
    std::vector<LO> arrayReceiveLO(5*numReceives);
    m_distributorSubData->doPostsAndWaits
      (Teuchos::ArrayView<const LO>(arraySendLO), 5,
       Teuchos::ArrayView<LO>(arrayReceiveLO));
    // communicate coarse data of types LO, GO, SM, and SX
    std::vector<size_t> receiveLengthLO(numReceives), 
      receiveLengthGO(numReceives), receiveLengthSM(numReceives), 
      receiveLengthSX(numReceives);
    LO lengthLO(0), lengthGO(0), lengthSM(0), lengthSX(0);
    numElems = 0;
    subReceivePtr.resize(numReceives+1, 0);
    for (LO i=0; i<numReceives; i++) {
      numElems          += arrayReceiveLO[5*i+0];
      receiveLengthLO[i] = arrayReceiveLO[5*i+1]; lengthLO += receiveLengthLO[i];
      receiveLengthGO[i] = arrayReceiveLO[5*i+2]; lengthGO += receiveLengthGO[i];
      receiveLengthSM[i] = arrayReceiveLO[5*i+3]; lengthSM += receiveLengthSM[i];
      receiveLengthSX[i] = arrayReceiveLO[5*i+4]; lengthSX += receiveLengthSX[i];
      subReceivePtr[i+1] = numElems;
    }
    // append on-processor data as needed
    const LO numSubsForRank = subsForRank.size();
    LO startLOSelf(lengthLO), startGOSelf(lengthGO), startSMSelf(lengthSM),
      startSXSelf(lengthSX);
    if (numSubsForRank > numSends) {
      BDDC_TEST_FOR_EXCEPTION(numSubsForRank != numSends+1, std::runtime_error, 
			    "invalid numSubsForRank");
      numElems += subsForRank[numSends].size();
      lengthLO += startDataLO[numSends+1] - startDataLO[numSends];
      lengthGO += startDataGO[numSends+1] - startDataGO[numSends];
      lengthSM += startDataSM[numSends+1] - startDataSM[numSends];
      lengthSX += startDataSX[numSends+1] - startDataSX[numSends];
      subReceivePtr.push_back(numElems);
    }
    dataReceiveLO.resize(lengthLO);
    dataReceiveGO.resize(lengthGO);
    dataReceiveSM.resize(lengthSM);
    dataReceiveSX.resize(lengthSX);
    if (numSubsForRank > numSends) {
      const LO* sourceDataLO = &dataLO[startDataLO[numSends]];
      const GO* sourceDataGO = &dataGO[startDataGO[numSends]];
      const SM* sourceDataSM = &dataSM[startDataSM[numSends]];
      const SX* sourceDataSX = &dataSX[startDataSX[numSends]];
      LO* targetDataLO = &dataReceiveLO[startLOSelf];
      GO* targetDataGO = &dataReceiveGO[startGOSelf];
      SM* targetDataSM = &dataReceiveSM[startSMSelf];
      SX* targetDataSX = &dataReceiveSX[startSXSelf];
      memcpy(targetDataLO, sourceDataLO, (lengthLO-startLOSelf)*sizeof(LO));
      memcpy(targetDataGO, sourceDataGO, (lengthGO-startGOSelf)*sizeof(GO));
      memcpy(targetDataSM, sourceDataSM, (lengthSM-startSMSelf)*sizeof(SM));
      memcpy(targetDataSX, sourceDataSX, (lengthSX-startSXSelf)*sizeof(SX));
    }
    const LO numTargetRanks = numSends;
    // LO data
    LO lengthToSend = startDataLO[numTargetRanks];
    const LO* sendLO = ((lengthToSend == 0) ? nullptr : dataLO.data());
    LO* recvLO = ((startLOSelf == 0) ? nullptr : dataReceiveLO.data());
    m_distributorSubData->doPostsAndWaits
      (Teuchos::ArrayView<const LO>(sendLO, lengthToSend),
       Teuchos::ArrayView<const size_t>(sendLengthLO),
       Teuchos::ArrayView<LO>(recvLO, startLOSelf),
       Teuchos::ArrayView<const size_t>(receiveLengthLO));
    // GO data
    lengthToSend = startDataGO[numTargetRanks];
    const GO* sendGO = ((lengthToSend == 0) ? nullptr : dataGO.data());
    GO* recvGO = ((startGOSelf == 0) ? nullptr : dataReceiveGO.data());
    m_distributorSubData->doPostsAndWaits
      (Teuchos::ArrayView<const GO>(sendGO, lengthToSend),
       Teuchos::ArrayView<const size_t>(sendLengthGO),
       Teuchos::ArrayView<GO>(recvGO, startGOSelf),
       Teuchos::ArrayView<const size_t>(receiveLengthGO));
    // SM data
    lengthToSend = startDataSM[numTargetRanks];
    const SM* sendSM = ((lengthToSend == 0) ? nullptr : dataSM.data());
    SM* recvSM = ((startSMSelf == 0) ? nullptr : dataReceiveSM.data());
    m_distributorSubData->doPostsAndWaits
      (Teuchos::ArrayView<const SM>(sendSM, lengthToSend),
       Teuchos::ArrayView<const size_t>(sendLengthSM),
       Teuchos::ArrayView<SM>(recvSM, startSMSelf),
       Teuchos::ArrayView<const size_t>(receiveLengthSM));
    // SX data
    lengthToSend = startDataSX[numTargetRanks];
    const SX* sendSX = ((lengthToSend == 0) ? nullptr : dataSX.data());
    SX* recvSX = ((startSXSelf == 0) ? nullptr : dataReceiveSX.data());
    m_distributorSubData->doPostsAndWaits
      (Teuchos::ArrayView<const SX>(sendSX, lengthToSend),
       Teuchos::ArrayView<const size_t>(sendLengthSX),
       Teuchos::ArrayView<SX>(recvSX, startSXSelf),
       Teuchos::ArrayView<const size_t>(receiveLengthSX));
  }

  void convertSubData(std::vector<SubdomainData> & coarseSubData,
		      LO & numNodes,
		      LO & numFaces,
		      std::vector<LO> & nodeBegin,
		      std::vector<LO> & localDofs,
		      std::vector<GO> & nodeGlobalIDsCoarse,
		      std::vector<SM> & xCoords,
		      std::vector<SM> & yCoords,
		      std::vector<SM> & zCoords,
		      std::vector<int> & minPart)
  {
    const LO numElems = coarseSubData.size();
    std::map<GO, LO> nodeMap;
    numNodes = 0;
    nodeBegin.push_back(0);
    for (LO i=0; i<numElems; i++) {
      SubdomainData & subData = coarseSubData[i];
      subData.nodeLIDs.resize(subData.numNodes);
      for (LO j=0; j<subData.numNodes; j++) {
	const GO globalID = subData.nodeGIDs[j];
	auto found = nodeMap.find(globalID);
	if (found == nodeMap.end()) {
	  nodeMap[globalID] = numNodes;
	  minPart.push_back(subData.minPart[j]);
	  subData.nodeLIDs[j] = numNodes++;
	  xCoords.push_back(subData.coords[3*j+0]);
	  yCoords.push_back(subData.coords[3*j+1]);
	  zCoords.push_back(subData.coords[3*j+2]);
	  for (LO k=subData.nodeBegin[j]; k<subData.nodeBegin[j+1]; k++) {
	    localDofs.push_back(subData.localDofs[k]);
	  }
	  nodeBegin.push_back(localDofs.size());
	  nodeGlobalIDsCoarse.push_back(globalID);
	}
	else {
	  const LO node = found->second;
	  subData.nodeLIDs[j] = node;
	  // checks
	  BDDC_TEST_FOR_EXCEPTION(subData.minPart[j] != minPart[node], 
		      std::runtime_error, "subData.minPart[j] is invalid");
	  const LO numDofNode = subData.nodeBegin[j+1] - subData.nodeBegin[j];
	  const LO numDofNode2 = nodeBegin[node+1] - nodeBegin[node];
	  BDDC_TEST_FOR_EXCEPTION(numDofNode2 != numDofNode, 
		      std::runtime_error, "invalid numDofNode2");
	  for (LO k=nodeBegin[node]; k<nodeBegin[node+1]; k++) {
	    const LO index = subData.nodeBegin[j] + k - nodeBegin[node];
	    BDDC_TEST_FOR_EXCEPTION(localDofs[k] != subData.localDofs[index], 
		       std::runtime_error, "localDofs[k] is invalid");
	  }
	}
      }
    }
    std::map<GO, LO> faceMap;
    numFaces = 0;
    for (LO i=0; i<numElems; i++) {
      SubdomainData & subData = coarseSubData[i];
      subData.faceLIDs.resize(subData.numFaces);
      for (LO j=0; j<subData.numFaces; j++) {
	const GO globalID = subData.faceGIDs[j];
	auto found = faceMap.find(globalID);
	if (found == faceMap.end()) {
	  faceMap[globalID] = numFaces;
	  subData.faceLIDs[j] = numFaces++;
	}
	else {
	  subData.faceLIDs[j] = found->second;
	}
      }
    }
  }

  void packSendAndReceiveCoarseData
    (const std::vector<LO> & activeCoarseNodes, 
     const std::vector<GO> & activeCoarseNodeGIDs,
     const std::vector<SM> & coarseCoords, 
     const std::vector<LO> & nodeBegin, 
     const std::vector<LO> & localDofs,
     const std::vector< std::vector<LO> > & subsForRank,
     const std::vector<int> & targetRanks,
     const std::vector<int> & sourceRanks,
     const std::vector<LO> & startDataLO,
     const std::vector<LO> & startDataGO,
     const std::vector<LO> & startDataSM,
     const std::vector<LO> & startDataSX,
     const std::vector<int> & minPartForActiveCoarseNodes,
     LO & numElems,
     LO & numReceives,
     std::vector<LO> & dataReceiveLO,
     std::vector<GO> & dataReceiveGO,
     std::vector<SM> & dataReceiveSM,
     std::vector<SX> & dataReceiveSX,
     std::vector<LO> & subReceivePtr)
  {
    std::vector<LO> dataLO;
    std::vector<GO> dataGO;
    std::vector<SM> dataSM;
    std::vector<SX> dataSX;
    packCoarseData(activeCoarseNodes, activeCoarseNodeGIDs, coarseCoords,
		   nodeBegin, localDofs, subsForRank, targetRanks, startDataLO,
		   startDataGO, startDataSM, startDataSX, 
		   minPartForActiveCoarseNodes, dataLO, dataGO, dataSM, dataSX);

    m_distributorSubData = rcp( new Tpetra::Distributor(m_Comm) );
    m_distributorSubData->createFromSends
      (Teuchos::ArrayView<const int>(targetRanks));
    
    sendCoarseData(subsForRank, startDataLO, startDataGO, startDataSM, 
		   startDataSX, dataLO, dataGO, dataSM, dataSX,
		   dataReceiveLO, dataReceiveGO, dataReceiveSM, 
		   dataReceiveSX, numElems, numReceives, subReceivePtr);
  }

  void getActiveCoarseNodeLIDs
    (const std::vector<GO> & activeCoarseNodeGIDs, 
     std::vector<LO> & activeCoarseNodeLIDs)
  {
    std::map<GO,LO> nodeMap;
    for (LO i=0; i<m_numNodes; i++) {
      const GO globalID = m_nodeGlobalIDs[i];
      nodeMap.emplace(globalID, i);
    }
    const LO numActive = activeCoarseNodeGIDs.size();
    activeCoarseNodeLIDs.resize(numActive);
    for (LO i=0; i<numActive; i++) {
      const GO globalID = activeCoarseNodeGIDs[i];
      auto iter = nodeMap.find(globalID);
      BDDC_TEST_FOR_EXCEPTION(iter == nodeMap.end(), std::runtime_error, 
			      "globalID not found");
      activeCoarseNodeLIDs[i] = iter->second;
    }
  }

  void getActiveCoarseNodeParts
    (const std::vector<LO> & activeCoarseNodes,
     const std::vector<GO> & activeCoarseNodeGIDs,
     const std::vector<GO> & partsAcrossMpiRanks,
     std::vector< std::vector<GO> > & activeCoarseNodeParts)
  {
    // activeCoarseNodeParts[i] has all the parts for the active coarse
    // node with globalID activeCoarseNodeGIDs[i]
    std::vector<LO> activeCoarseNodeLIDs;
    getActiveCoarseNodeLIDs(activeCoarseNodeGIDs, activeCoarseNodeLIDs);
    const LO numRows = activeCoarseNodeGIDs.size();
    std::vector<int> rowSend(numRows);
    std::vector< std::vector<GO> > inputData(numRows);
    const std::vector< std::vector<LO> > & coarseNodeSubdomains = 
      m_Partition->getCoarseNodeSubdomains();
    const int myStartingSub = m_Partition->getStartingSub();
    for (LO i=0; i<numRows; i++) {
      rowSend[i] = m_nodeSend[activeCoarseNodeLIDs[i]];
      const LO equiv = activeCoarseNodes[i];
      const LO length = coarseNodeSubdomains[equiv].size();
      for (LO j=0; j<length; j++) {
	const GO subGID = coarseNodeSubdomains[equiv][j];
	const LO localSub = subGID - myStartingSub;
	if ((localSub >= 0) && (localSub < m_numSub)) {
	  inputData[i].push_back(partsAcrossMpiRanks[localSub]);
	}
      }
    }
    bddc::unionData<LO,GO>
      (numRows, activeCoarseNodeGIDs.data(), rowSend.data(),
       m_distributor, inputData, activeCoarseNodeParts);
    // check sizes
    for (LO i=0; i<numRows; i++) {
      const LO equiv = activeCoarseNodes[i];
      const bool test = (activeCoarseNodeParts[i].size() == 
	      coarseNodeSubdomains[equiv].size());
      BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			      "invalid acticeCoarseNodeParts[i] size");
    }
  }

  void extractCoarseData(LO numElems, 
			 const LO numReceives,
			 std::vector<LO> & dataReceiveLO,
			 std::vector<GO> & dataReceiveGO,
			 std::vector<SM> & dataReceiveSM,
			 std::vector<SX> & dataReceiveSX,
			 std::vector<SubdomainData> & coarseSubData)
  {
    coarseSubData.resize(numElems);
    const LO numElemsAll = numElems;
    numElems = 0;
    LO startLO(0), startGO(0), startSM(0), startSX(0);
    for (LO i=0; i<numElemsAll; i++) {
      extractData(startLO, startGO, startSM, startSX, dataReceiveLO,
		  dataReceiveGO, dataReceiveSM, dataReceiveSX,
		  coarseSubData[numElems++]);
    }
  }

  void extractData(LO & startLO, 
		   LO & startGO, 
		   LO & startSM, 
		   LO & startSX, 
		   std::vector<LO> & dataReceiveLO,
		   std::vector<GO> & dataReceiveGO, 
		   std::vector<SM> & dataReceiveSM, 
		   std::vector<SX> & dataReceiveSX,
		   SubdomainData & coarseSubData)
  {
    const LO numNodes = dataReceiveLO[startLO++];
    const LO numRows  = dataReceiveLO[startLO++];
    const LO numFaces = dataReceiveLO[startLO++];
    coarseSubData.numNodes = numNodes;
    coarseSubData.numRows  = numRows;
    coarseSubData.numFaces = numFaces;
    coarseSubData.nodeBegin = &dataReceiveLO[startLO];
    startLO += numNodes + 1;
    coarseSubData.localDofs = &dataReceiveLO[startLO];
    startLO += numRows;
    coarseSubData.minPart  = &dataReceiveLO[startLO];
    startLO += numNodes;
    coarseSubData.nodeGIDs = &dataReceiveGO[startGO];
    startGO += numNodes;
    coarseSubData.faceGIDs = &dataReceiveGO[startGO];
    startGO += numFaces;
    coarseSubData.coords = &dataReceiveSM[startSM];
    startSM += 3*numNodes;
    coarseSubData.values = &dataReceiveSX[startSX];
    startSX += numRows*numRows;
  }

  void setNodeMap(const std::vector<LO> & coarseSubs,
		  const std::vector<SubdomainData> & coarseSubData,
		  std::vector<LO> & nodeMap, 
		  std::vector<LO> & activeNodes, 
		  std::vector<LO> & subNodes,
		  const LO numSubs)
  {
    LO numNodeSub(0);
    for (size_t i=0; i<coarseSubs.size(); i++) {
      const LO sub = coarseSubs[i];
      const SubdomainData & subData = coarseSubData[sub];
      for (LO j=0; j<subData.numNodes; j++) {
	const LO node = subData.nodeLIDs[j];
	if (nodeMap[node] == -1) {
	  nodeMap[node] = numNodeSub;
	  activeNodes[numNodeSub++] = node;
	}
      }
    }
    subNodes.resize(numNodeSub);
    for (LO i=0; i<numNodeSub; i++) {
      subNodes[i] = activeNodes[i];
    }
    // avoid potential ordering problems downstream
    if (numSubs == 1) {
      std::sort(subNodes.begin(), subNodes.end());
      for (LO i=0; i<numNodeSub; i++) {
	const LO node = subNodes[i];
	nodeMap[node] = i;
      }
    }
  }

  void resetNodeMap(std::vector<LO> & nodeMap, 
		    const std::vector<LO> & subNodes)
  {
    for (size_t i=0; i<subNodes.size(); i++) {
      nodeMap[subNodes[i]] = -1;
    }
  }

  void determineNodeGraph(const std::vector<SubdomainData> & coarseSubData, 
			  const std::vector<LO> & nodeMap, 
			  const LO numNodes,
			  const std::vector<LO> & nodeElems,
			  const std::vector<LO> & nodeElemsPtr,
			  std::vector<LO> & activeNodes,
			  std::vector<bool> & nodeFlag,
			  std::vector<LO> & nodalConn,
			  std::vector<LO> & nodalConnPtr)
  {
    nodalConnPtr.resize(numNodes+1, 0);
    for (LO i=0; i<numNodes; i++) {
      LO numActive(0);
      for (LO j=nodeElemsPtr[i]; j<nodeElemsPtr[i+1]; j++) {
	const LO sub = nodeElems[j];
	const SubdomainData & subData = coarseSubData[sub];
	for (LO k=0; k<subData.numNodes; k++) {
	  const LO node = nodeMap[subData.nodeLIDs[k]];
	  if (nodeFlag[node] == false) {
	    nodeFlag[node] = true;
	    activeNodes[numActive++] = node;
	  }
	}
      }
      nodalConnPtr[i+1] = nodalConnPtr[i] + numActive;
      for (LO j=0; j<numActive; j++) nodeFlag[activeNodes[j]] = false;
    }
    nodalConn.resize(nodalConnPtr[numNodes]);
    for (LO i=0; i<numNodes; i++) {
      LO numActive(0);
      for (LO j=nodeElemsPtr[i]; j<nodeElemsPtr[i+1]; j++) {
	const LO sub = nodeElems[j];
	const SubdomainData & subData = coarseSubData[sub];
	for (LO k=0; k<subData.numNodes; k++) {
	  const LO node = nodeMap[subData.nodeLIDs[k]];
	  if (nodeFlag[node] == false) {
	    nodeFlag[node] = true;
	    activeNodes[numActive++] = node;
	  }
	}
      }
      for (LO j=0; j<numActive; j++) {
	const LO node = activeNodes[j];
	nodeFlag[node] = false;
	const LO index = nodalConnPtr[i] + j;
	nodalConn[index] = node;
      }
    }
  }

  void determineElemsForNodes(const std::vector<LO> & coarseSubs, 
			      const std::vector<SubdomainData> & coarseSubData, 
			      const std::vector<LO> & nodeMap, 
			      const LO numNodes,
			      std::vector<LO> & nodeElems,
			      std::vector<LO> & nodeElemsLocalNode,
			      std::vector<LO> & nodeElemsPtr)
  {
    const LO numSubs = coarseSubs.size();
    std::vector<LO> count(numNodes, 0);
    for (LO i=0; i<numSubs; i++) {
      const LO sub = coarseSubs[i];
      const SubdomainData & subData = coarseSubData[sub];
      for (LO j=0; j<subData.numNodes; j++) {
	const LO node = nodeMap[subData.nodeLIDs[j]];
	BDDC_TEST_FOR_EXCEPTION(node == -1, std::runtime_error, 
				"node is invalid");
	count[node]++;
      }
    }
    nodeElemsPtr.resize(numNodes+1, 0);
    for (LO i=0; i<numNodes; i++) {
      nodeElemsPtr[i+1] = nodeElemsPtr[i] + count[i];
      count[i] = 0;
    }
    const LO numTerms = nodeElemsPtr[numNodes];
    nodeElems.resize(numTerms);
    nodeElemsLocalNode.resize(numTerms);
    for (LO i=0; i<numSubs; i++) {
      const LO sub = coarseSubs[i];
      const SubdomainData & subData = coarseSubData[sub];
      for (LO j=0; j<subData.numNodes; j++) {
	const LO node = nodeMap[subData.nodeLIDs[j]];
	const LO index = nodeElemsPtr[node] + count[node];
	nodeElems[index] = sub;
	nodeElemsLocalNode[index] = j;
	count[node]++;
      }
    }
  }

  void determineSubdomainGraph(const std::vector<LO> & nodeBegin,
			       const std::vector<LO> & subNodes, 
			       const std::vector<LO> & nodalConn,
			       const std::vector<LO> & nodalConnPtr,
			       std::vector<LO> & nodeBeginSub,
			       std::vector<LO> & rowBegin, 
			       std::vector<LO> & columns)
  {
    const LO numNodes = subNodes.size();
    LO numRows(0);
    nodeBeginSub.resize(numNodes+1, 0);
    for (LO i=0; i<numNodes; i++) {
      const LO node = subNodes[i];
      numRows += nodeBegin[node+1] - nodeBegin[node];
      nodeBeginSub[i+1] = numRows;
    }
    rowBegin.resize(numRows+1, 0);
    numRows = 0;
    // first pass to determine number of nonzeros
    for (LO i=0; i<numNodes; i++) {
      LO numCols(0);
      for (LO j=nodalConnPtr[i]; j<nodalConnPtr[i+1]; j++) {
	const LO node = nodalConn[j];
	numCols += nodeBeginSub[node+1] - nodeBeginSub[node];
      }
      for (LO j=nodeBeginSub[i]; j<nodeBeginSub[i+1]; j++) {
	rowBegin[numRows+1] = rowBegin[numRows] + numCols;
	numRows++;
      }
    }
    LO numTerms = rowBegin[numRows];
    columns.resize(numTerms);
    numTerms = 0;
    // second pass to determine nonzero columns for each row
    for (LO i=0; i<numNodes; i++) {
      for (LO j=nodeBeginSub[i]; j<nodeBeginSub[i+1]; j++) {
	for (LO k=nodalConnPtr[i]; k<nodalConnPtr[i+1]; k++) {
	  const LO node = nodalConn[k];
	  for (LO m=nodeBeginSub[node]; m<nodeBeginSub[node+1]; m++) {
	    columns[numTerms++] = m;
	  }
	}
      }
    }
  }

  void assembleSubdomain(const std::vector<SubdomainData> & coarseSubData, 
			 const LO numNodes,
			 const std::vector<LO> & nodeMap,
			 std::vector<LO> & dofMap, 
			 const std::vector<LO> & nodeBeginSub, 
			 const std::vector<LO> & nodeElems,
			 const std::vector<LO> & nodeElemsLocalNode,
			 const std::vector<LO> & nodeElemsPtr,
			 const std::vector<LO> & rowBegin, 
			 const std::vector<LO> & columns, 
			 std::vector<SX> & values)
  {
    const LO numRows = nodeBeginSub[numNodes];
    LO numTerms = rowBegin[numRows];
    values.resize(numTerms, 0);
    for (LO i=0; i<numNodes; i++) {
      const LO firstRow = nodeBeginSub[i];
      for (LO j=rowBegin[firstRow]; j<rowBegin[firstRow+1]; j++) {
	const LO col = columns[j];
	dofMap[col] = j - rowBegin[firstRow];
      }
      for (LO row=nodeBeginSub[i]; row<nodeBeginSub[i+1]; row++) {
	const LO start = rowBegin[row];
	const LO localIndex = row - nodeBeginSub[i];
	for (LO j=nodeElemsPtr[i]; j<nodeElemsPtr[i+1]; j++) {
	  const LO sub = nodeElems[j];
	  const SubdomainData & subData = coarseSubData[sub];
	  const LO subNode = nodeElemsLocalNode[j];
	  const LO subRow = subData.nodeBegin[subNode] + localIndex;
	  for (LO k=0; k<subData.numNodes; k++) {
	    const LO subNode2 = nodeMap[subData.nodeLIDs[k]];
	    const LO numDofSubNode2 = 
	      subData.nodeBegin[k+1] - subData.nodeBegin[k];
	    for (LO m=0; m<numDofSubNode2; m++) {
	      const LO subCol = subData.nodeBegin[k] + m;
	      const LO subIndex = subRow + subData.numRows*subCol;
	      const LO col = nodeBeginSub[subNode2] + m;
	      const LO rowCol = dofMap[col];
	      BDDC_TEST_FOR_EXCEPTION(rowCol == -1, std::runtime_error, 
				      "rowCol is invalid");
	      const LO index = start + rowCol;
	      values[index] += subData.values[subIndex];
	    }
	  }
	}
      }
      for (LO j=rowBegin[firstRow]; j<rowBegin[firstRow+1]; j++) {
	const LO col = columns[j];
	dofMap[col] = -1; // needed only for BDDC_TEST above
      }
    }
  }

  void assembleSubMatrices
    (const std::vector<SubdomainData> & coarseSubData, 
     const std::vector< std::vector<LO> > & coarseSubdomainSubs, 
     const std::vector<LO> & nodeBegin,
     const std::vector<LO> & localDofs, 
     std::vector< std::vector<LO> > & rowBeginSub, 
     std::vector< std::vector<LO> > & columnsSub, 
     std::vector< std::vector<SX> > & valuesSub, 
     std::vector< std::vector<LO> > & subNodesCoarse)
  {
    const LO numSubs = coarseSubdomainSubs.size();
    rowBeginSub.resize(numSubs);
    columnsSub.resize(numSubs);
    valuesSub.resize(numSubs);
    subNodesCoarse.resize(numSubs);
    const LO numNodes = nodeBegin.size() - 1;
    std::vector<LO> nodeMap(numNodes, -1);
    std::vector<LO> activeNodes(numNodes);
    std::vector<bool> nodeFlag(numNodes, false);
    const LO numRows = localDofs.size();
    std::vector<LO> dofMap(numRows, -1), nodeBeginSub(numNodes);
    for (LO i=0; i<numSubs; i++) {
      const std::vector<LO> & coarseSubs = coarseSubdomainSubs[i];
      std::vector<LO> & subNodes = subNodesCoarse[i];
      setNodeMap(coarseSubs, coarseSubData, nodeMap, activeNodes, subNodes,
		 numSubs);
      std::vector<LO> nodeElems, nodeElemsLocalNode, nodeElemsPtr;
      determineElemsForNodes(coarseSubs, coarseSubData, nodeMap, subNodes.size(),
			     nodeElems, nodeElemsLocalNode, nodeElemsPtr);

      std::vector<LO> nodalConn, nodalConnPtr;
      determineNodeGraph(coarseSubData, nodeMap, subNodes.size(),
			 nodeElems, nodeElemsPtr, activeNodes, nodeFlag,
			 nodalConn, nodalConnPtr);
      
      determineSubdomainGraph(nodeBegin, subNodes, nodalConn, nodalConnPtr, 
			      nodeBeginSub, rowBeginSub[i], columnsSub[i]);
      
      assembleSubdomain(coarseSubData, subNodes.size(), nodeMap, dofMap, 
			nodeBeginSub, nodeElems, nodeElemsLocalNode, 
			nodeElemsPtr, rowBeginSub[i], columnsSub[i], 
			valuesSub[i]);

      resetNodeMap(nodeMap, subNodes);
    }
  }

  void constructCoarseMaps(const std::vector<GO> & nodeGlobalIDsCoarse, 
			   const std::vector<LO> & nodeBeginCoarse, 
			   const std::vector<LO> & localDofsCoarse,
			   RCP<const Map> & dofMapCoarse, 
			   RCP<const Map> & dofMapCoarse1to1)
  {
    std::vector<GO> globalIDs1to1;
    if (m_coarseProc == true) {
      RCP<const Teuchos::Comm<int> > Comm = 
	rcp( new Teuchos::MpiComm<int>(m_mpiCommSplit) );
      LO numNodesCoarse = nodeGlobalIDsCoarse.size();
      bddc::DofManager<LO,GO>::
	determineGlobalIDs(numNodesCoarse, nodeGlobalIDsCoarse.data(), 
			   nodeBeginCoarse.data(), localDofsCoarse.data(), 
			   Comm, dofMapCoarse, dofMapCoarse1to1);

      makeOwnedMapConsistent
	(dofMapCoarse, nodeGlobalIDsCoarse.size(), nodeBeginCoarse.data(),
	 m_nodeSendCoarse, dofMapCoarse1to1);

      const LO numRows = dofMapCoarse1to1->getNodeNumElements();
      globalIDs1to1.resize(numRows);
      for (LO i=0; i<numRows; i++) {
	globalIDs1to1[i] = dofMapCoarse1to1->getGlobalElement(i);
      }
    }
    m_coarseMapNextLevel = generateMap(globalIDs1to1);
  }

  void addSendData(std::vector<GO> & sendData, 
		   LO & sendLength, 
		   const std::vector< std::pair<GO,LO> > & globalAndLocalIDs, 
		   const std::vector<int> & nodeSendCoarse,
		   const std::vector<int> & mpiRanks)
  {
    const LO numNodes = globalAndLocalIDs.size();
    for (LO j=0; j<numNodes; j++) {
      sendData[sendLength++] = globalAndLocalIDs[j].first;
      const LO localNode = globalAndLocalIDs[j].second;
      GO owner = nodeSendCoarse[localNode];
      if (owner == -1) owner = m_myPID;
      else owner = mpiRanks[owner];
      sendData[sendLength++] = owner;
    }
  }

  void addSendData(std::vector<GO> & sendData, 
		   LO & sendLength, 
		   const std::vector< std::pair<GO,LO> > & globalAndLocalIDs, 
		   const std::vector<LO> & nodeBeginCoarse,
		   RCP<const Map> & dofMapCoarse)
  {
    const LO numNodes = globalAndLocalIDs.size();
    for (LO j=0; j<numNodes; j++) {
      sendData[sendLength++] = globalAndLocalIDs[j].first;
      const LO localNode     = globalAndLocalIDs[j].second;
      const LO row = nodeBeginCoarse[localNode];
      sendData[sendLength++] = dofMapCoarse->getGlobalElement(row);
    }
  }

  void getUniqueIDs(const LO i,
		    const std::vector<SubdomainData> & coarseSubData,
		    const std::vector<LO> & subReceivePtr,
		    std::vector< std::pair<GO, LO> > & globalAndLocalIDs)
  {
    LO numTerms(0);
    for (LO j=subReceivePtr[i]; j<subReceivePtr[i+1]; j++) {
      numTerms += coarseSubData[j].numNodes;
    }
    globalAndLocalIDs.resize(numTerms);
    numTerms = 0;
    for (LO j=subReceivePtr[i]; j<subReceivePtr[i+1]; j++) {
      const SubdomainData & subData = coarseSubData[j];
      for (LO k=0; k<subData.numNodes; k++) {
	globalAndLocalIDs[numTerms++] = std::make_pair(subData.nodeGIDs[k],
						       subData.nodeLIDs[k]);
      }
    }
    std::sort(globalAndLocalIDs.begin(), globalAndLocalIDs.end());
    auto iter = std::unique(globalAndLocalIDs.begin(), globalAndLocalIDs.end());
    globalAndLocalIDs.erase(iter, globalAndLocalIDs.end());
  }

  void reverseCommunicateNodeOwner
    (const std::vector<int> & nodeSendCoarse, 
     const std::vector<LO> & mpiRanks, 
     const std::vector<GO> & nodeGlobalIDsCoarse, 
     const std::vector<LO> & subReceivePtr,
     const std::vector<GO> & activeCoarseNodeGIDs,
     const std::vector<SubdomainData> & coarseSubData,
     std::vector<int> & nodeOwner)
  {
    // Note: next two lines swapped since doing reverse communications
    size_t numSends = m_distributorSubData->getNumReceives();
    size_t numReceives = m_distributorSubData->getNumSends();
    LO sendLength(0);
    std::vector< std::pair<GO, LO> > globalAndLocalIDs;
    for (size_t i=0; i<numSends; i++) {
      getUniqueIDs(i, coarseSubData, subReceivePtr, globalAndLocalIDs);
      sendLength += 2*globalAndLocalIDs.size();
    }
    std::vector<GO> sendData(sendLength);
    sendLength = 0;
    std::vector<size_t> lengthsSend(numSends, 0);
    for (size_t i=0; i<numSends; i++) {
      getUniqueIDs(i, coarseSubData, subReceivePtr, globalAndLocalIDs);
      lengthsSend[i] = 2*globalAndLocalIDs.size();
      addSendData(sendData, sendLength, globalAndLocalIDs, nodeSendCoarse,
		  mpiRanks);
    }
    std::vector<size_t> lengthsReceive(numReceives);
    m_distributorSubData->doReversePostsAndWaits
      (Teuchos::ArrayView<const size_t>(lengthsSend), 1,
       Teuchos::ArrayView<size_t>(lengthsReceive));
    LO receiveLength(0);
    for (size_t i=0; i<numReceives; i++) receiveLength += lengthsReceive[i];
    LO extraLength(0);
    if (subReceivePtr.size() > numSends+1) { // have on-proc data
      BDDC_TEST_FOR_EXCEPTION(subReceivePtr.size() != numSends+2, 
			std::runtime_error, "subReceivePtr size is invalid");
      getUniqueIDs(numSends, coarseSubData, subReceivePtr, globalAndLocalIDs);
      extraLength = 2*globalAndLocalIDs.size();
    }
    std::vector<GO> receiveData(receiveLength + extraLength);
    const GO* sendGO = ((sendLength == 0) ? nullptr : sendData.data());
    GO* recvGO = ((receiveLength == 0) ? nullptr : receiveData.data());
    m_distributorSubData->doReversePostsAndWaits
      (Teuchos::ArrayView<const GO>(sendGO, sendLength),
       Teuchos::ArrayView<const size_t>(lengthsSend),
       Teuchos::ArrayView<GO>(recvGO, receiveLength),
       Teuchos::ArrayView<const size_t>(lengthsReceive));
    if (subReceivePtr.size() > numSends+1) { // append on-proc data
      getUniqueIDs(numSends, coarseSubData, subReceivePtr, globalAndLocalIDs);
      addSendData(receiveData, receiveLength, globalAndLocalIDs,
		  nodeSendCoarse, mpiRanks);
    }
    // extract nodeOwner
    std::map<GO,LO> nodeMap;
    const LO numActive = activeCoarseNodeGIDs.size();
    for (LO i=0; i<numActive; i++) {
      nodeMap.emplace(activeCoarseNodeGIDs[i], i);
    }
    nodeOwner.resize(numActive, -1);
    for (size_t i=0; i<receiveData.size()/2; i++) {
      const GO globalID = receiveData[2*i];
      const int owner = receiveData[2*i+1];
      auto iter = nodeMap.find(globalID);
      BDDC_TEST_FOR_EXCEPTION(iter == nodeMap.end(), std::runtime_error, 
			      "globlID not found");
      const LO node = iter->second;
      if (nodeOwner[node] == -1) nodeOwner[node] = owner;
      else {
	BDDC_TEST_FOR_EXCEPTION(nodeOwner[node] != owner, std::runtime_error, 
				"nodeOwnder[node] is invalid");
      }
    }
    for (LO i=0; i<numActive; i++) {
      BDDC_TEST_FOR_EXCEPTION(nodeOwner[i] == -1, std::runtime_error, 
			      "nodeOwner[i] is invalid");
    }
  }

  void reverseCommunicateNodeGlobalIDsAndStartingLocations
    (const std::vector<GO> & nodeGlobalIDsCoarse, 
     const std::vector<LO> & nodeBeginCoarse, 
     RCP<const Map> & dofMapCoarse, 
     const std::vector<SubdomainData> & coarseSubData,
     const std::vector<LO> & subReceivePtr,
     std::vector<GO> & receiveData)
  {
    // Note: next two lines swapped since doing reverse communications
    size_t numSends = m_distributorSubData->getNumReceives();
    size_t numReceives = m_distributorSubData->getNumSends();
    LO sendLength(0);
    std::vector< std::pair<GO, LO> > globalAndLocalIDs;
    for (size_t i=0; i<numSends; i++) {
      getUniqueIDs(i, coarseSubData, subReceivePtr, globalAndLocalIDs);
      sendLength += 2*globalAndLocalIDs.size();
    }
    std::vector<GO> sendData(sendLength);
    sendLength = 0;
    std::vector<size_t> lengthsSend(numSends, 0);
    for (size_t i=0; i<numSends; i++) {
      getUniqueIDs(i, coarseSubData, subReceivePtr, globalAndLocalIDs);
      lengthsSend[i] = 2*globalAndLocalIDs.size();
      addSendData(sendData, sendLength, globalAndLocalIDs, nodeBeginCoarse,
		  dofMapCoarse);
    }
    std::vector<size_t> lengthsReceive(numReceives);
    m_distributorSubData->doReversePostsAndWaits
      (Teuchos::ArrayView<const size_t>(lengthsSend), 1,
       Teuchos::ArrayView<size_t>(lengthsReceive));
    LO receiveLength(0);
    for (size_t i=0; i<numReceives; i++) {
      receiveLength += lengthsReceive[i];
    }
    LO extraLength(0);
    if (subReceivePtr.size() > numSends+1) { // have on-proc data
      BDDC_TEST_FOR_EXCEPTION(subReceivePtr.size() != numSends+2, 
	      std::runtime_error, "invalid subReceivePtr size");
      getUniqueIDs(numSends, coarseSubData, subReceivePtr, globalAndLocalIDs);
      extraLength = 2*globalAndLocalIDs.size();
    }
    receiveData.resize(receiveLength + extraLength);
    const GO* sendGO = ((sendLength == 0) ? nullptr : sendData.data());
    GO* recvGO = ((receiveLength == 0) ? nullptr : receiveData.data());
    m_distributorSubData->doReversePostsAndWaits
      (Teuchos::ArrayView<const GO>(sendGO, sendLength),
       Teuchos::ArrayView<const size_t>(lengthsSend),
       Teuchos::ArrayView<GO>(recvGO, receiveLength),
       Teuchos::ArrayView<const size_t>(lengthsReceive));
    if (subReceivePtr.size() > numSends+1) { // append on-proc data
      getUniqueIDs(numSends, coarseSubData, subReceivePtr, globalAndLocalIDs);
      addSendData(receiveData, receiveLength, globalAndLocalIDs,
		  nodeBeginCoarse, dofMapCoarse);
    }
  }

  void determineStartingLocations(const std::vector<GO> & nodeData,
				  const std::vector<LO> & activeCoarseNodes,
				  const std::vector<GO> & activeCoarseNodeGIDs,
				  std::vector<GO> & nodeStart)
  {
    const LO numNodes = activeCoarseNodes.size();
    std::map<GO, LO> nodeMap;
    for (LO i=0; i<numNodes; i++) {
      nodeMap[activeCoarseNodeGIDs[i]] = i;
    }
    nodeStart.resize(numNodes, -1);
    for (size_t i=0; i<nodeData.size()/2; i++) {
      const GO globalID = nodeData[2*i+0];
      const GO start    = nodeData[2*i+1];
      auto found = nodeMap.find(globalID);
      BDDC_TEST_FOR_EXCEPTION(found == nodeMap.end(), std::runtime_error, 
			      "globalID not found");
      const LO node = found->second;
      if (nodeStart[node] == -1) {
	nodeStart[node] = start;
      }
      else {
	BDDC_TEST_FOR_EXCEPTION(nodeStart[node] != start, std::runtime_error, 
			      "nodeStart[node] is invalid");
      }
    }
    for (LO i=0; i<numNodes; i++) {
      BDDC_TEST_FOR_EXCEPTION(nodeStart[i] == -1, std::runtime_error, 
			      "nodeStart[i] is invalid");
    }
  }

  void determineCoarseMaps(const std::vector<GO> & nodeStart,
			   const std::vector<LO> & nodeBegin)
  {
    const LO numNodes = nodeStart.size();
    LO numRows = nodeBegin[numNodes];
    std::vector<GO> rowGIDs(numRows);
    numRows = 0;
    for (LO i=0; i<numNodes; i++) {
      const LO start = nodeStart[i];
      const LO numDofNode = nodeBegin[i+1] - nodeBegin[i];
      for (LO j=0; j<numDofNode; j++) {
	rowGIDs[numRows++] = start + j;
      }
    }
    m_coarseMap = generateMap(rowGIDs);
  }

  void constructCoarseMapsReverse
    (const std::vector<GO> & nodeGlobalIDsCoarse, 
     const std::vector<LO> & nodeBeginCoarse,
     RCP<const Map> & dofMapCoarse, 
     const std::vector<SubdomainData> & coarseSubData,
     const std::vector<LO> & activeCoarseNodes,
     const std::vector<GO> & activeCoarseNodeGIDs,
     const std::vector<LO> & nodeBegin,
     const std::vector<LO> & subReceivePtr)
  {
    const bool test = (activeCoarseNodes.size() == nodeBegin.size()-1);
    BDDC_TEST_FOR_EXCEPTION(test == false, std::runtime_error, 
			    "activeCoarseNodes size is invalid");
    std::vector<GO> nodeData;
    reverseCommunicateNodeGlobalIDsAndStartingLocations
      (nodeGlobalIDsCoarse, nodeBeginCoarse, dofMapCoarse, coarseSubData,
       subReceivePtr, nodeData);
 
    std::vector<GO> nodeStart;
    determineStartingLocations(nodeData, activeCoarseNodes,
			       activeCoarseNodeGIDs, nodeStart);

    determineCoarseMaps(nodeStart, nodeBegin);
  }

  void getMinPartForActiveCoarseNodes
    (const std::vector< std::vector<GO> > & activeCoarseNodeParts, 
     std::vector<int> & minPartForActiveCoarseNodes)
  {
    const LO numActive = activeCoarseNodeParts.size();
    minPartForActiveCoarseNodes.resize(numActive);
    for (LO i=0; i<numActive; i++) {
      const LO numPart = activeCoarseNodeParts[i].size();
      BDDC_TEST_FOR_EXCEPTION(numPart <= 0, std::runtime_error, 
			      "numPart must be positive");
      int minPart = activeCoarseNodeParts[i][0];
      for (LO j=1; j<numPart; j++) {
	const int part = activeCoarseNodeParts[i][j];
	if (part < minPart) minPart = part;
      }
      minPartForActiveCoarseNodes[i] = minPart;
    }
  }

  void determineNodeSend(const std::vector<int> & minPart, 
			 std::vector<int> & nodeSend)
  {
    const LO numNode = minPart.size();
    nodeSend.resize(numNode);
    for (LO i=0; i<numNode; i++) {
      if (m_myCoarseMpiRank == minPart[i]) {
	nodeSend[i] = -1; // own this node and only receive data
      }
      else {
	nodeSend[i] = minPart[i]; // send data to minPart[i]
      }
    }
  }

  void constructCoarseDistributorAndPComm
    (const std::vector<GO> & activeCoarseNodeGIDs, 
     const std::vector<LO> & nodeBeginActive, 
     const LO numNodesCoarse, 
     const std::vector<LO> & nodeBeginCoarse,
     const std::vector<GO> & nodeGlobalIDsCoarse, 
     const std::vector<int> & nodeOwner, 
     const std::vector<int> & nodeSendCoarse)
  {
    const LO numActive = activeCoarseNodeGIDs.size();
    std::vector<int> nodeSend(numActive);
    for (LO i=0; i<numActive; i++) {
      if (nodeOwner[i] == m_myPID) nodeSend[i] = -1;
      else nodeSend[i] = nodeOwner[i];
    }
    bddc::constructDistributor<LO>
      (numActive, nodeSend, m_Comm, m_distributorCoarse);
    std::vector<LO> nodeBeginOwned(1, 0);
    std::vector<GO> nodeGlobalIDsOwned;
    LO numDofOwned(0);
    for (LO i=0; i<numNodesCoarse; i++) {
      if (nodeSendCoarse[i] == -1) {
	nodeGlobalIDsOwned.push_back(nodeGlobalIDsCoarse[i]);
	numDofOwned += nodeBeginCoarse[i+1] - nodeBeginCoarse[i];
	nodeBeginOwned.push_back(numDofOwned);
      }
    }
    const LO numNodeOwned = nodeGlobalIDsOwned.size();
    m_pCommCoarse = 
      rcp( new PComm<LO,GO,SX>
	   (m_distributorCoarse, numActive, nodeBeginActive.data(), 
	    activeCoarseNodeGIDs.data(), nodeSend.data(), numNodeOwned,
	    nodeBeginOwned.data(), nodeGlobalIDsOwned.data()) );
  }

  void constructCoarseSubdomainsNew
    (const std::vector<GO> & partAcrossMpiRanks,
     const std::vector<LO> & mpiRanks,
     const std::vector<int> & sourceRanks,
     std::vector<LO> & activeCoarseNodes,
     RCP<const Map> & dofMapCoarse,
     RCP<const Map> & dofMapCoarse1to1,
     std::vector<GO> & activeCoarseNodeGIDs)
  {
    // activeCoarseNodes: list of active equivalence classes (coarse nodes)
    // activeCoarseNodeGIDS: globalIDs of activeCoarseNodes
    // coarseCoords: coordinates of activeCoarseNodes (all x, all y, all z)
    // nodeBegin: starting locations for activeCoarseNodes
    // localDofs: local degree of freedom numbers for activeCoarseNodes
    std::vector<SM> coarseCoords;
    std::vector<LO> nodeBegin, localDofs;
    determineActiveCoarseNodeData
      (activeCoarseNodes, activeCoarseNodeGIDs, coarseCoords,
       nodeBegin, localDofs);
    // activeCoarseNodeParts[i] = all part (coarse proc) numbers for active 
    //          coarse node which has globalID activeCoarseNodeGIDs[i]
    std::vector< std::vector<GO> > activeCoarseNodeParts;
    getActiveCoarseNodeParts
      (activeCoarseNodes, activeCoarseNodeGIDs, partAcrossMpiRanks, 
       activeCoarseNodeParts);
    // minPartForActiveCoarseNodes[i] = smallest part number which contains
    //   coarse node with globalID activeCoarseNodeGIDs[i]
    std::vector<int> minPartForActiveCoarseNodes;
    getMinPartForActiveCoarseNodes
      (activeCoarseNodeParts, minPartForActiveCoarseNodes);
    // m_subDofsCoarse[i]: coarse degree of freedom numbers (local) for 
    //                     subdomain i
    determineCoarseSubdomainDofs(activeCoarseNodes, nodeBegin, m_subDofsCoarse);

    // Note: We do not include ourselves in targetRanks (the list of 
    //       MPI process IDs to send data to) because there appears to be 
    //       a problem with the 4-parameter version of
    //       Tpretra_Distributor.doPostsAndWaits if one's own processor
    //       ID is included. Doing so may also be a little more efficient 
    //       (e.g. avoid unnecessary copies). However, subsForRank
    //       does include on-processor subdomains to be sent to self at
    //       its end.
    std::vector<int> targetRanks;
    std::vector< std::vector<LO> > subsForRank;
    determineSubsForRank
      (partAcrossMpiRanks, mpiRanks, targetRanks, subsForRank);
    // We later pack all coarse data of types LO, GO, SM and SX to be sent
    // to other processors. The pointers startDataLO, startDataGO,
    // startDataSM, and startDataSX provide offsets for this data
    std::vector<LO> startDataLO, startDataGO, startDataSM, startDataSX;
    determinePointersForDataToSend
      (activeCoarseNodes, nodeBegin, subsForRank, startDataLO, startDataGO,
       startDataSM, startDataSX);

    LO numElems, numReceives;
    std::vector<LO> dataReceiveLO;
    std::vector<GO> dataReceiveGO;
    std::vector<SM> dataReceiveSM;
    std::vector<SX> dataReceiveSX;
    std::vector<LO> subReceivePtr;
    
    packSendAndReceiveCoarseData
      (activeCoarseNodes, activeCoarseNodeGIDs, coarseCoords, nodeBegin,
       localDofs, subsForRank, targetRanks, sourceRanks,
       startDataLO, startDataGO, startDataSM, startDataSX, 
       minPartForActiveCoarseNodes, numElems, numReceives, dataReceiveLO, 
       dataReceiveGO, dataReceiveSM, dataReceiveSX, subReceivePtr);
    
    std::vector<SubdomainData> coarseSubData;
    extractCoarseData(numElems, numReceives, dataReceiveLO, dataReceiveGO,
		      dataReceiveSM, dataReceiveSX, coarseSubData);

    LO numFacesCoarse;
    std::vector<int> minPart;
    convertSubData(coarseSubData, m_numNodesCoarse, numFacesCoarse, 
		   m_nodeBeginCoarse, m_localDofsCoarse, m_nodeGlobalIDsCoarse,
		   m_xCoordCoarse, m_yCoordCoarse, m_zCoordCoarse,
		   minPart);
    determineNodeSend(minPart, m_nodeSendCoarse);

    std::vector<int> nodeOwner;
    reverseCommunicateNodeOwner
      (m_nodeSendCoarse, mpiRanks, m_nodeGlobalIDsCoarse, subReceivePtr, 
       activeCoarseNodeGIDs, coarseSubData, nodeOwner);
    constructCoarseDistributorAndPComm
      (activeCoarseNodeGIDs, nodeBegin, m_numNodesCoarse, m_nodeBeginCoarse,
       m_nodeGlobalIDsCoarse, nodeOwner, m_nodeSendCoarse);

    if (m_usePComm == false) {
      constructCoarseMaps(m_nodeGlobalIDsCoarse, m_nodeBeginCoarse, 
			  m_localDofsCoarse, dofMapCoarse, dofMapCoarse1to1);
      constructCoarseMapsReverse(m_nodeGlobalIDsCoarse, m_nodeBeginCoarse,
				 dofMapCoarse, coarseSubData, activeCoarseNodes,
				 activeCoarseNodeGIDs, nodeBegin, subReceivePtr);
    }
    m_distributorSubData = Teuchos::null;

    std::vector< std::vector<LO> > coarseSubdomainSubs;
    partitionSubdomainsWithinRank(coarseSubData, numFacesCoarse, 
				  coarseSubdomainSubs);

    std::vector< std::vector<LO> > rowBeginSub, columnsSub, subNodesCoarse;
    std::vector< std::vector<SX> > valuesSub;
    assembleSubMatrices(coarseSubData, coarseSubdomainSubs, m_nodeBeginCoarse, 
			m_localDofsCoarse, m_subRowBeginCoarse, 
			m_subColumnsCoarse, m_subValuesCoarse, m_subNodesCoarse);
  }

  void determineCoarseSpace(std::vector<LO> & mpiRanks)
  {
    double startTime = GetTime();
    // mpiRanks is a vector of active MPI ranks at the coarse level 
    //   (this is the same vector across all current level MPI ranks)
    // m_myCoarseMpiRank is the MPI rank at the coarse level (-1 if inactive)
    determineCoarseMpiRanks(m_myCoarseMpiRank, mpiRanks);
    // m_fineProc   = true if used to store subdomain data at current level
    // m_coarseProc = true if used to store subdomain data at coarse level
    determineProcessorTypes(m_fineProc, m_coarseProc);
    m_disjointCommunicators = determineIfCommunicatorsAreDisjoint();
    // partsAcrossMpiRanks[i] is the destination coarse MPI rank for local
    // subdomain i
    std::vector<GO> partsAcrossMpiRanks;
    std::vector<int> importProcsLB, importPartsLB;
    partitionSubdomainsAcrossMpiRanks(mpiRanks, partsAcrossMpiRanks);

    splitCommunicator();

    std::vector<LO> activeCoarseNodes;
    std::vector<GO> activeCoarseNodeGIDs;
    RCP<Export> nodeExporter;
    RCP<Import> nodeImporter;
    RCP<const Map> dofMapCoarse, dofMapCoarse1to1;
    std::vector<int> sourceRanks;
    constructCoarseSubdomainsNew
      (partsAcrossMpiRanks, mpiRanks, sourceRanks, activeCoarseNodes,
       dofMapCoarse, dofMapCoarse1to1, activeCoarseNodeGIDs);
    allocateCoarseEntities();
    // vertex coarse space option
    if (m_useVertexCoarseSpace && (mpiRanks.size() == 1)) {
      const int numSubCoarse = m_subNodesCoarse.size();
      int numSubCoarseMax;
      Teuchos::reduceAll<int, int>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				   &numSubCoarse, &numSubCoarseMax);
      if (numSubCoarseMax == 1) {
	bool iAmTargetProc = false;
	if (m_myPID == mpiRanks[0]) iAmTargetProc = true;
#ifdef VERTEXCOARSESPACEBDDC
	const bool useEconomicVersion = 
	  m_Parameters->get("Economic Vertex Coarse Space", false);
	m_vertexCoarseSpace = 
	  rcp( new VertexCoarseSpace<SX,SM,LO,GO>
	       (m_Partition, m_pCommCoarse, m_numNodes, m_nodeGlobalIDs, m_nodeSend,
		activeCoarseNodes, activeCoarseNodeGIDs, m_nodeGlobalIDsCoarse,
		m_nodeBeginCoarse, m_localDofsCoarse, m_xCoordCoarse, m_yCoordCoarse, 
		m_zCoordCoarse, m_problemType, iAmTargetProc, useEconomicVersion) );
#endif
	(void)(iAmTargetProc);
      }
    }

    m_timings[TIME_COARSE_SPACE_PREP] += GetTime() - startTime;
    if (m_coarseProc == true) {
      LO numNodesCoarse = m_numNodesCoarse;
      LO* nodeBeginCoarse = m_nodeBeginCoarse.data();
      LO* localDofsCoarse = m_localDofsCoarse.data();
      GO* nodeGlobalIDsCoarse = m_nodeGlobalIDsCoarse.data();
      std::vector<GO>* nodeGlobalIDsCoarse1to1 = &m_nodeGlobalIDsCoarse1to1;
      SM* xCoordCoarse = m_xCoordCoarse.data();
      SM* yCoordCoarse = m_yCoordCoarse.data();
      SM* zCoordCoarse = m_zCoordCoarse.data();
      std::vector< std::vector<LO> >* subNodesCoarse = &m_subNodesCoarse;
      std::vector< std::vector<LO> >* subRowBeginCoarse = &m_subRowBeginCoarse;
      std::vector< std::vector<LO> >* subColumnsCoarse = &m_subColumnsCoarse;
      std::vector< std::vector<SX> >* subValuesCoarse = &m_subValuesCoarse;
      std::vector<int>* nodeSend = &m_nodeSendCoarse;
      std::vector< LO* > subRowBeginPtr, subColumnsPtr;
      std::vector< SX* > subValuesPtr;
      OperatorBase<SX>* coarseOperator(nullptr);     
      convertToPointers
	(subRowBeginCoarse, subColumnsCoarse, subValuesCoarse,
	 subRowBeginPtr, subColumnsPtr, subValuesPtr);
      m_coarsePreconditioner = 
	rcp( new PreconditionerBDDC<SX,SM,LO,GO>
	 (numNodesCoarse, nodeBeginCoarse, localDofsCoarse, 
	  nodeGlobalIDsCoarse, xCoordCoarse, yCoordCoarse, zCoordCoarse, 
	  *subNodesCoarse, subRowBeginPtr.data(), subColumnsPtr.data(), 
	  subValuesPtr.data(), m_Parameters, m_mpiCommSplit, m_level+1, 
	  coarseOperator, nodeSend, nodeGlobalIDsCoarse1to1, 
	  dofMapCoarse, dofMapCoarse1to1
#ifdef VERTEXCOARSESPACEBDDC
	  , m_vertexCoarseSpace
#endif
	  ) );
      m_numCoarseDofs = m_coarsePreconditioner->getNumDof();
    }
  }
  
  void convertToPointers(std::vector< std::vector<LO> >* subRowBeginCoarse, 
			 std::vector< std::vector<LO> >* subColumnsCoarse, 
			 std::vector< std::vector<SX> >* subValuesCoarse,
			 std::vector< LO* > & subRowBeginPtr, 
			 std::vector< LO* > & subColumnsPtr,
			 std::vector< SX* > & subValuesPtr)
  {    
    const LO numSub = subRowBeginCoarse->size();
    subRowBeginPtr.resize(numSub);
    subColumnsPtr.resize(numSub);
    subValuesPtr.resize(numSub);
    for (LO i=0; i<numSub; i++) {
      subRowBeginPtr[i] = (*subRowBeginCoarse)[i].data();
      subColumnsPtr[i] = (*subColumnsCoarse)[i].data();
      subValuesPtr[i] = (*subValuesCoarse)[i].data();
    }
  }

  void initializeSubdomains()
  {
    m_Subdomain.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      m_Parameters->set("subdomain number", i);
      m_Subdomain[i] = new bddc::SubdomainBDDC<SX,SM,LO,GO>
	(m_subNodes[i].size(), m_subNodes[i].data(), &m_subNodeBegin[i][0], 
	 &m_subLocalDofs[i][0],
	 &m_subRowBegin[i][0], &m_subColumns[i][0], &m_subValues[i][0],
	 m_xCoord, m_yCoord, m_zCoord, m_boundaryDofsLocal[i].size(), 
	 &m_boundaryDofsLocal[i][0], *m_Parameters);
    }
  }

  void factorInteriorMatrices()
  {
    double startTime = GetTime();
    for (LO i=0; i<m_numSub; i++) {
      m_Subdomain[i]->factorInteriorMatrix();
    }
    m_timings[TIME_INTERIOR_FACTORIZATIONS] += GetTime() - startTime;
  }

  void determineBaseConstraints()
  {
    double startTime = GetTime();
    m_Constraints = 
      rcp( new bddc::ConstraintsBDDC<SX,SM,LO,GO>
	   (m_numNodes, m_nodeBegin, m_localDofs, m_xCoord, m_yCoord, m_zCoord,
	    m_subNodes, m_Subdomain, m_Partition, m_diagBoundary, m_Parameters) );
    m_Constraints->determineBaseConstraints();
    m_timings[TIME_BASE_CONSTRAINTS] += GetTime() - startTime;
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
    if (m_usePComm) {
      std::vector<SX> diagBoundary(m_numDofsB);
      for (LO i=0; i<m_numDofsB; i++) diagBoundary[i] = m_diagBoundary[i];
      SX* diagBoundary1to1 = m_pCommB->getOwnedPtr();
      m_pCommB->doExport(diagBoundary.data(), diagBoundary1to1);
      m_pCommB->doImport(diagBoundary1to1, diagBoundary.data());
      for (LO i=0; i<m_numDofsB; i++) m_diagBoundary[i] = diagBoundary[i];
    }
    else {
      VectorSM diagBoundary(m_dofMapB), diagBoundary1to1(m_dofMapB1to1);
      for (LO i=0; i<m_numDofsB; i++) {
	diagBoundary.replaceLocalValue(i, m_diagBoundary[i]);
      }
      diagBoundary1to1.putScalar(0);
      diagBoundary1to1.doExport(diagBoundary, *m_exporterB, Tpetra::ADD);
      diagBoundary.doImport(diagBoundary1to1, *m_exporterB, Tpetra::INSERT);
      Teuchos::ArrayRCP<const SM> values = diagBoundary.getData();
      for (LO i=0; i<m_numDofsB; i++) m_diagBoundary[i] = values[i];
    }
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
	BDDC_TEST_FOR_EXCEPTION(nodeLocal == -1, std::runtime_error, 
			      "nodeLocal is invalid");
	for (LO k=nodeBeginSub[nodeLocal]; k<nodeBeginSub[nodeLocal+1]; k++) {
	  boundaryDofsLocal.push_back(k);
	}
	for (LO k=m_nodeBegin[nodeGlobal]; k<m_nodeBegin[nodeGlobal+1]; k++) {
	  boundaryDofsGlobal.push_back(k);
	}
	LO numDofsLocal = nodeBeginSub[nodeLocal+1] - nodeBeginSub[nodeLocal];
	LO numDofsGlobal = m_nodeBegin[nodeGlobal+1] - m_nodeBegin[nodeGlobal];
	BDDC_TEST_FOR_EXCEPTION(numDofsLocal != numDofsGlobal, 
				std::runtime_error, "numDofsLocal is invalid");
      }
    }
  }

};

} // namespace bddc

#endif // BDDC_PRECONDITIONER_H
  
