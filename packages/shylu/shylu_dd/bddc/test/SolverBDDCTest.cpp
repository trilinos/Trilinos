
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

#include <gtest/gtest.h>
#include <mpi.h>
#include "ProblemMakerBDDC.hpp"
#include "OperatorStandard.hpp"
#include "ShyLU_DDBDDC_config.h"
#include "shylu_PreconditionerBDDC.hpp"
#include "shylu_KrylovSolver.hpp"
#include "shylu_UtilBDDC.hpp"
#include "shylu_enumsBDDC.hpp"
#include "setupTest.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>

#if defined(_OPENMP)
#include <omp.h>
#endif

using Teuchos::RCP;

namespace {

typedef int LO; // Local Ordinal
typedef Tpetra::Map<>::global_ordinal_type GO; // Global Ordinal
typedef double SX; // floating point data type
typedef double SM; // real (magnitude) for SX

void readInputFiles(MPI_Comm Comm,
		    const std::string fileNamePM, 
		    const std::string fileNameBDDC, 
		    const std::string fileNameMueLu, 
		    const std::string fileNameNodalAMG, 
		    RCP<Teuchos::ParameterList> & parametersPM,
		    RCP<Teuchos::ParameterList> & parametersBDDC,
		    RCP<Teuchos::ParameterList> & parametersMueLu,
		    RCP<Teuchos::ParameterList> & parametersNodalAMG)
{
  RCP<const Teuchos::Comm<int> > TComm = 
    Teuchos::rcp( new Teuchos::MpiComm<int>(Comm) );
  Teuchos::ParameterList paramsPM, paramsBDDC, paramsMueLu, paramsNodalAMG;

  Teuchos::updateParametersFromXmlFileAndBroadcast
    (fileNamePM, Teuchos::Ptr<Teuchos::ParameterList>(&paramsPM), *TComm);
  parametersPM = Teuchos::rcp( new Teuchos::ParameterList(paramsPM) );
  
  Teuchos::updateParametersFromXmlFileAndBroadcast
    (fileNameBDDC, Teuchos::Ptr<Teuchos::ParameterList>(&paramsBDDC), *TComm);
  parametersBDDC = Teuchos::rcp( new Teuchos::ParameterList(paramsBDDC) );

  if (fileNameMueLu != "") {
    Teuchos::updateParametersFromXmlFileAndBroadcast
      (fileNameMueLu, Teuchos::Ptr<Teuchos::ParameterList>(&paramsMueLu), 
       *TComm);
    parametersMueLu = Teuchos::rcp( new Teuchos::ParameterList(paramsMueLu) );
  }

  if (fileNameNodalAMG != "") {
    Teuchos::updateParametersFromXmlFileAndBroadcast
      (fileNameNodalAMG, Teuchos::Ptr<Teuchos::ParameterList>(&paramsNodalAMG),
       *TComm);
    parametersNodalAMG = 
      Teuchos::rcp( new Teuchos::ParameterList(paramsNodalAMG) );
  }
}

void runTest(RCP<Teuchos::ParameterList> & parametersPM,
	     RCP<Teuchos::ParameterList> & parametersBDDC,
	     RCP<Teuchos::ParameterList> & parametersMueLu,
	     RCP<Teuchos::ParameterList> & parametersNodalAMG,
	     const std::string meshDataFile,
	     const int myPID,
	     MPI_Comm Comm,
	     const bool resetFile,
	     int & numIterations)
{
  // setup test data
  bddc::setupTest<LO,GO,SX,SM> 
    test(parametersPM, parametersBDDC, parametersMueLu,
	 parametersNodalAMG, meshDataFile, Comm);
  LO numNode, *nodeBegin(nullptr), *localDofs(nullptr);
  const GO *nodeGlobalIDs(nullptr);
  const SM *xCoord(nullptr), *yCoord(nullptr), *zCoord(nullptr);
  std::vector<LO*> subRowBeginPtr, subColumnsPtr;
  std::vector<SX*> subValuesPtr;
  bddc::OperatorBase<SX>* Operator(nullptr);
  // extract test data
  test.getProblemData(numNode, nodeBegin, localDofs, nodeGlobalIDs,
		      xCoord, yCoord, zCoord, subRowBeginPtr,
		      subColumnsPtr, subValuesPtr, Operator);
  const std::vector< std::vector<LO> > & subNodes = test.getSubNodes();
  std::vector<int> nodeSend;
  bddc::getNodeSend(numNode, nodeGlobalIDs, Comm, nodeSend);
  // initialize preconditioner
  int level(0);
  MPI_Barrier(Comm);
  double startTimeBDDCPre = test.getTime();
  RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > Preconditioner =
    rcp( new bddc::PreconditionerBDDC<SX,SM,LO,GO>
	 (numNode, nodeBegin, localDofs, nodeGlobalIDs, xCoord, yCoord, zCoord, 
	  subNodes, subRowBeginPtr.data(), subColumnsPtr.data(), 
	  subValuesPtr.data(), parametersBDDC, Comm, level, Operator, 
	  &nodeSend) );
  MPI_Barrier(Comm);
  double stopTimeBDDCPre = test.getTime();
  const LO numMyRows = Preconditioner->NumMyRows();
  std::vector<SX> rhs;
  test.getRhs(numMyRows, rhs);
  double startTimeKrylovInit(0), startTimeKrylovSolve(0);
  double stopTimeKrylovInit(0), stopTimeKrylovSolve(0);
  test.openFiles(Preconditioner, resetFile);
  MPI_Barrier(Comm);
  // initialize Krylov solver
  startTimeKrylovInit = test.getTime();
  RCP< bddc::KrylovSolver<SX,SM,LO,GO> > Solver =
    rcp ( new bddc::KrylovSolver<SX,SM,LO,GO>(Preconditioner, parametersBDDC) );
  std::vector<SX> sol(numMyRows), Ax(numMyRows);
  MPI_Barrier(Comm);
  // solve equations
  stopTimeKrylovInit = test.getTime();
  startTimeKrylovSolve = test.getTime();
  LO numberOfSolves = parametersPM->get("Number of Solves", 1);
  for (int i=0; i<numberOfSolves; i++) {
    if (i > 0) test.updateRhs(rhs.data(), numMyRows);
    Solver->Solve(rhs.data(), sol.data());
  }
  MPI_Barrier(Comm);
  stopTimeKrylovSolve = test.getTime();
  Preconditioner->ApplyFullOperator(sol.data(), Ax.data());
  for (LO i=0; i<numMyRows; i++) {
    Ax[i] = rhs[i] - Ax[i];
  }
  // check solution
  SM normError = Preconditioner->Norm2(Ax.data(), numMyRows);
  SM normRhs = Preconditioner->Norm2(rhs.data(), numMyRows);
  double solverTolerance = parametersBDDC->get("Convergence Tolerance", 1e-6);
  EXPECT_LT(normError, 1.01*solverTolerance*normRhs);
  if (myPID == 0) {
    std::cout << "BDDC Preconditioner initialization time = "
	      << stopTimeBDDCPre - startTimeBDDCPre << std::endl;
    std::cout << "KrylovSolver initialization time        = "
	      << stopTimeKrylovInit - startTimeKrylovInit << std::endl;
    std::cout << "KrylovSolver solve solve time           = "
	      << stopTimeKrylovSolve - startTimeKrylovSolve << std::endl;
  }
  const bool outputTiming = parametersBDDC->get("Output Timing", false);
  if (outputTiming) {
    int krylovMethod = parametersBDDC->get("Krylov Method", 0);
    test.printTimings(Preconditioner, Solver, krylovMethod);
    Preconditioner->printTimings("BDDC_timers.dat");
  }
}

TEST(SolverBDDC, Test1)
{
  const std::string fileNamePM = "problemMaker.xml";
  const std::string fileNameBDDC = "bddc.xml";
  const std::string fileNameMueLu = "mueLu_SGS.xml";
  const std::string fileNameNodalAMG = "";
  const std::string meshDataFile = "";
  RCP<Teuchos::ParameterList> parametersPM, parametersBDDC, parametersMueLu,
    parametersNodalAMG;
  MPI_Comm Comm = MPI_COMM_WORLD;
  int myPID;
  MPI_Comm_rank(Comm, &myPID);  
  readInputFiles(Comm, fileNamePM, fileNameBDDC, fileNameMueLu, 
		 fileNameNodalAMG, parametersPM, parametersBDDC,
		 parametersMueLu, parametersNodalAMG);
  int numIterations;
  const bool runSingleTest = false;
  bool resetFile = true;
  if (runSingleTest) {
    runTest(parametersPM, parametersBDDC, parametersMueLu, parametersNodalAMG, 
	    meshDataFile, myPID, Comm, resetFile, numIterations);
  }
  else {
    const std::string solverDirichlet = 
      parametersBDDC->get("Dirichlet Solver", "SuperLU");
    const std::string solverNeumann = 
      parametersBDDC->get("Neumann Solver", "SuperLU");
    const std::string solverCoarse = 
      parametersBDDC->get("Coarse Solver", "SuperLU");
    for (int i=0; i<2; i++) {
      parametersBDDC->set("Dirichlet Solver", solverDirichlet);
      parametersBDDC->set("Neumann Solver", solverNeumann);
      parametersBDDC->set("Coarse Solver", solverCoarse);
      if (i == 0) {
	parametersPM->set("Problem Type String", "Poisson-3D");
	parametersBDDC->set("Krylov Solver", "PCG");
      }
      if (i == 1) {
	parametersPM->set("Problem Type String", "Elasticity-3D");
	parametersBDDC->set("Krylov Solver", "GCR");
      }
      int jMax = 2;
#if defined(HAVE_SHYLU_DDBDDC_MUELU)
      jMax++;
#endif
      for (int j=0; j<jMax; j++) {
	if (j == 0) parametersBDDC->set("Interface Preconditioner", true);
	if (j == 1) parametersBDDC->set("Interface Preconditioner", false);
	if (j == 2) {
	  parametersBDDC->set("Dirichlet Solver", "MueLu");
	  parametersBDDC->set("Neumann Solver", "MueLu");
	  parametersBDDC->set("Coarse Solver", "MueLu");
	}
	if ((i == 1) && (j == 2)) continue;
	runTest(parametersPM, parametersBDDC, parametersMueLu,
		parametersNodalAMG, meshDataFile, myPID, Comm, 
		resetFile, numIterations);
	resetFile = false;
      }
    }
  }
}

} // end namespace
