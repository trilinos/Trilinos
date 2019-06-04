
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
#else
#define omp_get_nested() 0
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

TEST(SolverBDDCSimple, Test1)
{
  const std::string fileNamePM = "problemMakerSimple.xml";
  const std::string fileNameBDDC = "bddc.xml";
  const std::string fileNameMueLu = "elasticity3D.xml";
  const std::string fileNameNodalAMG = "";
  const std::string meshDataFile = "";
  MPI_Comm Comm = MPI_COMM_WORLD;
  int myPID;
  MPI_Comm_rank(Comm, &myPID);  
  RCP<Teuchos::ParameterList> parametersPM, parametersBDDC, parametersMueLu,
    parametersNodalAMG;
  readInputFiles(Comm, fileNamePM, fileNameBDDC, fileNameMueLu, 
		 fileNameNodalAMG, parametersPM, parametersBDDC,
		 parametersMueLu, parametersNodalAMG);
  // setup test data
  bddc::setupTest<LO,GO,SX,SM> 
    test(parametersPM, parametersBDDC, parametersMueLu, parametersNodalAMG,
	 meshDataFile, Comm);
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
  BDDC_TEST_FOR_EXCEPTION(subNodes.size() != 1, std::runtime_error, 
			  "number of subdomains per MPI rank must be 1");  
  // initialize preconditioner
  int level(0);
  LO* rowBegin = subRowBeginPtr[0];
  LO* columns = subColumnsPtr[0];
  SX* values = subValuesPtr[0];
  // make adjustments for essential boundary conditions
  std::vector<SM> xUse(numNode), yUse(numNode), zUse(numNode);
  std::vector<GO> nodeGlobalIDsUse(numNode);
  LO numNodeUse(0);
  for (LO i=0; i<numNode; i++) {
    if (nodeBegin[i+1] > nodeBegin[i]) {
      xUse[numNodeUse] = xCoord[i];
      yUse[numNodeUse] = yCoord[i];
      zUse[numNodeUse] = zCoord[i];
      nodeGlobalIDsUse[numNodeUse++] = nodeGlobalIDs[i];
    }
  }
  std::vector<int> nodeSend;
  bddc::getNodeSend(numNodeUse, nodeGlobalIDsUse.data(), Comm, nodeSend);
  MPI_Barrier(Comm);
  double startTimeBDDCPre = test.getTime();
  RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > Preconditioner =
    rcp( new bddc::PreconditionerBDDC<SX,SM,LO,GO>
	 (numNodeUse, nodeGlobalIDsUse.data(), xUse.data(), yUse.data(), 
	  zUse.data(), rowBegin, columns, values, parametersBDDC, Comm, level, 
	  &nodeSend) );
  MPI_Barrier(Comm);
  double stopTimeBDDCPre = test.getTime();
  const LO numMyRows = Preconditioner->NumMyRows();
  std::vector<SX> rhs;
  test.getRhs(numMyRows, rhs);
  double startTimeKrylovInit(0), startTimeKrylovSolve(0);
  double stopTimeKrylovInit(0), stopTimeKrylovSolve(0);
  const bool resetFile = true;
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
  int krylovMethod = parametersBDDC->get("Krylov Method", 0);
  test.printTimings(Preconditioner, Solver, krylovMethod);
  if (myPID == 0) {
    std::cout << "BDDC Preconditioner initialization time = "
	      << stopTimeBDDCPre - startTimeBDDCPre << std::endl;
    std::cout << "KrylovSolver initialization time        = "
	      << stopTimeKrylovInit - startTimeKrylovInit << std::endl;
    std::cout << "KrylovSolver solve solve time           = "
	      << stopTimeKrylovSolve - startTimeKrylovSolve << std::endl;
  }
  Preconditioner->printTimings("BDDC_timers.dat");
}

} // end namespace
