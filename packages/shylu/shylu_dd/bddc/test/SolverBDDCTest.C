
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
#include "ProblemMakerBDDC.h"
#include "ShyLU_DDBDDC_config.h"
#include "shylu_PreconditionerBDDC.h"
#include "shylu_KrylovSolver.h"
#include "shylu_UtilBDDC.h"
#include "shylu_enumsBDDC.h"
#if defined(_OPENMP)

#include <omp.h>
#else
#define omp_get_nested() 0

#endif

using Teuchos::RCP;

namespace {

typedef int LO; // Local Ordinal
  //typedef int GO; // Global Ordinal
typedef long long GO; // Global Ordinal
typedef double SX; // floating point data type
typedef double SM; // real (magnitude) for SX
//
// Tpetra typedefs
//
typedef Tpetra::Map<>::node_type  Node;
typedef Tpetra::Map<LO,GO,Node>                                 Map;
typedef Tpetra::CrsMatrix<SX,LO,GO,Node>                        CrsMatrix;
typedef Tpetra::Import<LO,GO,Node>                              Import;
typedef Tpetra::Export<LO,GO,Node>                              Export;
typedef Tpetra::Vector<SX,LO,GO,Node>                           Vector;

void openFiles
(RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner,
 LO krylovMethod);

void printTimings
(RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner,
 RCP< bddc::KrylovSolver<SX,SM,LO,GO> > & Solver,
 const LO krylovMethod);

void readInputFile(double & lengthDir1,
                   double & lengthDir2,
                   double & lengthDir3,
                   LO & numSubDir1,
                   LO & numSubDir2,
                   LO & numSubDir3,
                   LO & numSubDir1PerProc,
                   LO & numSubDir2PerProc,
                   LO & numSubDir3PerProc,
                   LO & Hh,
                   LO & spatialDim,
                   LO & problemTypeInt,
                   LO & matrixType,
                   double & diagScaleFactor,
                   LO & numThreadsOuter,
                   LO & numThreadsInner,
                   SM & solverTolerance,
                   LO & maxIterations,
                   LO & krylovMethod,
                   LO & estimateConditionNumber,
                   LO & directSolver,
                   LO & cornerOption,
                   LO & edgeOption,
                   LO & faceOption,
                   LO & numCoarseSubdomainsPerMpiRank,
                   LO & numSubdomainsPerCoarseSubdomain,
                   LO & coarseningOption,
                   LO & applyOnlyPreconditioner,
                   LO & numberOfSolves);

int setParameters(double lengthDir1,
                  double lengthDir2,
                  double lengthDir3,
                  LO numSubDir1,
                  LO numSubDir2,
                  LO numSubDir3,
                  LO numSubDir1PerProc,
                  LO numSubDir2PerProc,
                  LO numSubDir3PerProc,
                  LO Hh,
                  LO spatialDim,
                  LO problemTypeInt,
                  LO matrixType,
                  double diagScaleFactor,
                  LO numThreadsOuter,
                  SM solverTolerance,
                  LO maxIterations,
                  LO krylovMethod,
                  LO estimateConditionNumber,
                  LO directSolver,
                  LO cornerOption,
                  LO edgeOption,
                  LO faceOption,
                  LO numCoarseSubdomainsPerMpiRank,
                  LO numSubdomainsPerCoarseSubdomain,
                  LO coarseningOption,
                  int numProc,
                  int myPID,
                  RCP<Teuchos::ParameterList> & Parameters);

void updateRhs(SX* rhs,
               LO numRows);

TEST(SolverBDDC, Test1)
{
  int numProc, myPID;
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Comm_size(Comm, &numProc);
  MPI_Comm_rank(Comm, &myPID);
  if (myPID == 0) {
    LO indexLO(0);
    GO indexGO(0);
    std::cout << "sizeof(LO) = " << sizeof(indexLO) << std::endl;
    std::cout << "sizeof(GO) = " << sizeof(indexGO) << std::endl;
  }
  double lengthDir1, lengthDir2, lengthDir3, diagScaleFactor, solverTolerance;
  LO numSubDir1, numSubDir2, numSubDir3, numSubDir1PerProc,
    numSubDir2PerProc, numSubDir3PerProc, Hh, spatialDim, problemTypeInt,
    numThreadsOuter, numThreadsInner, maxIterations, directSolver,
    applyOnlyPreconditioner, numberOfSolves, matrixType, cornerOption,
    edgeOption, faceOption, krylovMethod, estimateConditionNumber,
    numCoarseSubdomainsPerMpiRank, numSubdomainsPerCoarseSubdomain,
    coarseningOption;
  //
  // We currently only work with structured meshes generated on-the-fly.
  // Later, we plan to read in matrices originating from unstructured meshes.
  // Matrices originating from unstructured meshes can be read in directly
  // for testing of threaded sparse direct solvers (see SparseSolverTest.C).
  //
  readInputFile(lengthDir1, lengthDir2, lengthDir3, numSubDir1, numSubDir2,
                numSubDir3, numSubDir1PerProc, numSubDir2PerProc,
                numSubDir3PerProc, Hh, spatialDim, problemTypeInt,
                matrixType, diagScaleFactor, numThreadsOuter, numThreadsInner,
                solverTolerance, maxIterations, krylovMethod,
                estimateConditionNumber, directSolver, cornerOption,
                edgeOption, faceOption, numCoarseSubdomainsPerMpiRank,
                numSubdomainsPerCoarseSubdomain, coarseningOption,
                applyOnlyPreconditioner, numberOfSolves);
  RCP<Teuchos::ParameterList> Parameters;
  int returnValue =
  setParameters(lengthDir1, lengthDir2, lengthDir3, numSubDir1, numSubDir2,
                numSubDir3, numSubDir1PerProc, numSubDir2PerProc,
                numSubDir3PerProc, Hh, spatialDim, problemTypeInt,
                matrixType, diagScaleFactor, numThreadsOuter, solverTolerance,
                maxIterations, krylovMethod, estimateConditionNumber,
                directSolver, cornerOption, edgeOption, faceOption,
                numCoarseSubdomainsPerMpiRank, numSubdomainsPerCoarseSubdomain,
                coarseningOption,numProc, myPID, Parameters);
  if (returnValue != 0) return;
  // generate problem
  RCP< bddc::ProblemMaker<LO,GO,SX,SM> > Problem =
    rcp( new bddc::ProblemMaker<LO,GO,SX,SM>(Parameters, Comm) );
  std::vector< std::vector<LO> > subNodes, subNodeBegin, subRowBegin,
    subLocalDofs, subColumns, subElems;
  std::vector< std::vector<SX> > subValues;
  Problem->getSubDomainElements(numSubDir1PerProc, numSubDir2PerProc,
                                numSubDir3PerProc, subElems);
  Problem->getSubDomainNodeData(subElems, subNodes, subNodeBegin,
                                subLocalDofs);
  Problem->getSubdomainMatrices(subElems, subNodes, subNodeBegin, subRowBegin,
                                subColumns, subValues);
  Problem->addDiagonalStiffness(subRowBegin, subColumns, subValues,
                                diagScaleFactor);
  if (matrixType == 1) {
    Problem->addAsymmetry(subRowBegin, subColumns, subValues);
  }
  LO numNode = Problem->getNumNode();
  LO numSub = subRowBegin.size();
  std::vector< LO* > subRowBeginPtr(numSub), subColumnsPtr(numSub);
  std::vector< SX* > subValuesPtr(numSub);
  for (LO i=0; i<numSub; i++) {
    subRowBeginPtr[i] = subRowBegin[i].data();
    subColumnsPtr[i] = subColumns[i].data();
    subValuesPtr[i] = subValues[i].data();
  }
  LO* nodeBegin = const_cast<LO*>(Problem->getNodeBegin());
  LO* localDofs = const_cast<LO*>(Problem->getLocalDofs());
  const GO* nodeGlobalIDs = Problem->getNodeGlobalIDs();
  const SM* xCoord = Problem->getXcoord();
  const SM* yCoord = Problem->getYcoord();
  const SM* zCoord = Problem->getZcoord();

  if (myPID == 0) {
    std::cout << "omp_get_nested() return " << omp_get_nested() << std::endl;
  }
  // Set number of threads (will be used for inner loops, number of threads
  // for outer loops set manually using num_threads clause). Thanks for idea
  // suggested by Andrew Bradley who observed environment variable
  // specifications like OMP_NUM_THREADS=2,4 are not always supported.
#if defined(_OPENMP)
  omp_set_num_threads(numThreadsInner);
#endif
  // initialize preconditioner
  int level(0);
  RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > Preconditioner =
    rcp( new bddc::PreconditionerBDDC<SX,SM,LO,GO>
         (numNode, nodeBegin, localDofs, nodeGlobalIDs, xCoord, yCoord, zCoord,
          subNodes, subRowBeginPtr.data(), subColumnsPtr.data(),
          subValuesPtr.data(), Parameters, Comm, level) );
  RCP<const Map> dofMap = Preconditioner->getDofMap1to1();
  if (applyOnlyPreconditioner) {
    dofMap = Preconditioner->getDofMapB1to1();
  }
  LO numDofs = dofMap->getNodeNumElements();
  // initialize right hand side
  std::vector<SX> rhs(numDofs, 0);
  for (LO i=0; i<numDofs; i++) {
    //    if ((i % spatialDim) == 0) rhs[i] = 1;
    rhs[i] = 0.7*rand()/RAND_MAX;
    //    rhs[i] = 1;
  }
  if (applyOnlyPreconditioner) {
    std::vector<SX> Pr(numDofs), APr(numDofs);
    for (int i=0; i<numberOfSolves; i++) {
      updateRhs(rhs.data(), numDofs);
      Preconditioner->Apply(rhs.data(), Pr.data(), APr.data());
    }
  }
  else {
    openFiles(Preconditioner, krylovMethod);
    RCP< bddc::KrylovSolver<SX,SM,LO,GO> > Solver =
      rcp ( new bddc::KrylovSolver<SX,SM,LO,GO>(Preconditioner, Parameters) );
    std::vector<SX> sol(numDofs), Ax(numDofs);
    for (int i=0; i<numberOfSolves; i++) {
      if (i > 0) updateRhs(rhs.data(), numDofs);
      Solver->Solve(rhs.data(), sol.data());
    }
    Preconditioner->ApplyFullOperator(sol.data(), Ax.data());
    for (LO i=0; i<numDofs; i++) {
      Ax[i] = rhs[i] - Ax[i];
    }
    SM normError = Preconditioner->Norm2(Ax.data(), numDofs);
    SM normRhs = Preconditioner->Norm2(rhs.data(), numDofs);
    EXPECT_LT(normError, 1.01*solverTolerance*normRhs);
    printTimings(Preconditioner, Solver, krylovMethod);
  }
}

void updateRhs(SX* rhs,
               LO numRows)
{
  for (LO i=0; i<numRows; i++) {
    rhs[i] += 0.07*rand()/RAND_MAX;
  }
}

void readInputFile(double & lengthDir1,
                   double & lengthDir2,
                   double & lengthDir3,
                   LO & numSubDir1,
                   LO & numSubDir2,
                   LO & numSubDir3,
                   LO & numSubDir1PerProc,
                   LO & numSubDir2PerProc,
                   LO & numSubDir3PerProc,
                   LO & Hh,
                   LO & spatialDim,
                   LO & problemTypeInt,
                   LO & matrixType,
                   double & diagScaleFactor,
                   LO & numThreadsOuter,
                   LO & numThreadsInner,
                   SM & solverTolerance,
                   LO & maxIterations,
                   LO & krylovMethod,
                   LO & estimateConditionNumber,
                   LO & directSolver,
                   LO & cornerOption,
                   LO & edgeOption,
                   LO & faceOption,
                   LO & numCoarseSubdomainsPerMpiRank,
                   LO & numSubdomainsPerCoarseSubdomain,
                   LO & coarseningOption,
                   LO & applyOnlyPreconditioner,
                   LO & numberOfSolves)
{
  char buff[101];
  std::ifstream fin;
  fin.open("SolverBDDCTest.inp");
  fin >> lengthDir1; fin.getline(buff,100);
  fin >> lengthDir2; fin.getline(buff,100);
  fin >> lengthDir3; fin.getline(buff,100);
  fin >> numSubDir1; fin.getline(buff,100);
  fin >> numSubDir2; fin.getline(buff,100);
  fin >> numSubDir3; fin.getline(buff,100);
  fin >> numSubDir1PerProc; fin.getline(buff,100);
  fin >> numSubDir2PerProc; fin.getline(buff,100);
  fin >> numSubDir3PerProc; fin.getline(buff,100);
  fin >> Hh; fin.getline(buff,100);
  fin >> spatialDim; fin.getline(buff,100);
  fin >> problemTypeInt; fin.getline(buff,100);
  fin >> matrixType; fin.getline(buff,100);
  fin >> diagScaleFactor; fin.getline(buff,100);
  fin >> numThreadsOuter; fin.getline(buff,100);
  fin >> numThreadsInner; fin.getline(buff,100);
  fin >> solverTolerance; fin.getline(buff,100);
  fin >> maxIterations; fin.getline(buff,100);
  fin >> krylovMethod; fin.getline(buff,100);
  fin >> estimateConditionNumber; fin.getline(buff,100);
  fin >> directSolver; fin.getline(buff,100);
  fin >> cornerOption; fin.getline(buff,100);
  fin >> edgeOption; fin.getline(buff,100);
  fin >> faceOption; fin.getline(buff,100);
  fin >> numCoarseSubdomainsPerMpiRank; fin.getline(buff,100);
  fin >> numSubdomainsPerCoarseSubdomain; fin.getline(buff,100);
  fin >> coarseningOption; fin.getline(buff,100);
  fin >> applyOnlyPreconditioner; fin.getline(buff,100);
  fin >> numberOfSolves; fin.getline(buff,100);
  fin.close();
}

int setParameters(double lengthDir1,
                  double lengthDir2,
                  double lengthDir3,
                  LO numSubDir1,
                  LO numSubDir2,
                  LO numSubDir3,
                  LO numSubDir1PerProc,
                  LO numSubDir2PerProc,
                  LO numSubDir3PerProc,
                  LO Hh,
                  LO spatialDim,
                  LO problemTypeInt,
                  LO matrixType,
                  double diagScaleFactor,
                  LO numThreadsOuter,
                  SM solverTolerance,
                  LO maxIterations,
                  LO krylovMethod,
                  LO estimateConditionNumber,
                  LO directSolver,
                  LO cornerOption,
                  LO edgeOption,
                  LO faceOption,
                  LO numCoarseSubdomainsPerMpiRank,
                  LO numSubdomainsPerCoarseSubdomain,
                  LO coarseningOption,
                  int numProc,
                  int myPID,
                  RCP<Teuchos::ParameterList> & Parameters)
{
  LO numElemPerSubDir1 = numSubDir1PerProc*Hh;
  LO numElemPerSubDir2 = numSubDir2PerProc*Hh;
  LO numElemPerSubDir3 = numSubDir3PerProc*Hh;
  enum bddc::ProblemType problemType = bddc::SCALARPDE;
  if (problemTypeInt == 2) problemType = bddc::ELASTICITY;
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  int loadDirection(0);
  if (spatialDim == 2) {
    numSubDir3 = 1;
    numSubDir3PerProc = 1;
    numElemPerSubDir3 = 1;
  }
  int numDofPerNode = 1;
  if (problemType == bddc::ELASTICITY) {
    numDofPerNode = spatialDim;
  }
  int numExpectedProc = numSubDir1*numSubDir2*numSubDir3;
  if (numProc < numExpectedProc) {
    std::cout << "number of processors not large enough\n";
    return 1;
  }
  //
  Parameters = Teuchos::rcp( new Teuchos::ParameterList() );
  Parameters->set("numDofPerNode", numDofPerNode);
  Parameters->set("numThreadsOuter", numThreadsOuter);
  Parameters->set("Spatial Dimension", spatialDim);
  Parameters->set("Problem Type", problemType);
  Parameters->set("Problem Type BDDC", problemType);
  Parameters->set("Analysis Type", analysisType);
  if (matrixType == 0) {
    Parameters->set("Matrix Type", "Symmetric");
  }
  else {
    Parameters->set("Matrix Type", "NonSymmetric");
  }
  Parameters->set("Length Direction 1", lengthDir1);
  Parameters->set("Length Direction 2", lengthDir2);
  Parameters->set("Length Direction 3", lengthDir3);
  Parameters->set("Number of Subdomains Direction 1", numSubDir1);
  Parameters->set("Number of Subdomains Direction 2", numSubDir2);
  Parameters->set("Number of Subdomains Direction 3", numSubDir3);
  Parameters->set("Number of Elements Per Subdomain Direction 1",
                  numElemPerSubDir1);
  Parameters->set("Number of Elements Per Subdomain Direction 2",
                  numElemPerSubDir2);
  Parameters->set("Number of Elements Per Subdomain Direction 3",
                  numElemPerSubDir3);
  Parameters->set("Apply Left Side Essential BCs", true);
  Parameters->set("Apply Right Side Essential BCs", false);
  Parameters->set("Load Direction", loadDirection);
  Parameters->set("Artificial Foundation Stiffness", 0.0);
  Parameters->set("omega", 0.0);
  Parameters->set("Generate Constraint Equations", false);
  Parameters->set("Interface Preconditioner", true);
  Parameters->set("Print Interior Matrices", false);
  enum bddc::WeightType weightTypeCorner = bddc::STIFFNESS;
  enum bddc::WeightType weightTypeEdge = bddc::STIFFNESS;
  enum bddc::WeightType weightTypeFace = bddc::STIFFNESS;
  Parameters->set("Weight Type Corner", weightTypeCorner);
  Parameters->set("Weight Type Edge", weightTypeEdge);
  Parameters->set("Weight Type Face", weightTypeFace);
  Parameters->set("Convergence Tolerance", solverTolerance);
  Parameters->set("Maximum Iterations", maxIterations);
  Parameters->set("Maximum Stored Directions", maxIterations);
  Parameters->set("Krylov Method", krylovMethod);
  Parameters->set("Estimate Condition Number", estimateConditionNumber);
  Parameters->set("Print Summary", 3);
  Parameters->set("Output Timing", false);
  bool checkCoarseMatrices = true;
  bool printCoarseMatrix = true;
  Parameters->set("Check Coarse Matrices", checkCoarseMatrices);
  Parameters->set("Print Coarse Matrix", printCoarseMatrix);
  if (cornerOption) Parameters->set("Use Corners", true);
  else Parameters->set("Use Corners", false);
  if (edgeOption) Parameters->set("Use Edges", true);
  else Parameters->set("Use Edges", false);
  if (faceOption) Parameters->set("Use Faces", true);
  else Parameters->set("Use Faces", false);
  // coarse space parameters
  Parameters->set("numCoarseSubdomainsPerMpiRank",
                  numCoarseSubdomainsPerMpiRank);
  Parameters->set("numSubdomainsPerCoarseSubdomain",
                  numSubdomainsPerCoarseSubdomain);
  if (coarseningOption == 0) {
    Parameters->set("Construct Subdomain Adjacency Graph", true);
    Parameters->set("Coarsening Option", "Graph");
  }
  else {
    Parameters->set("Construct Subdomain Adjacency Graph", false);
    if (coarseningOption == 1) {
      Parameters->set("Coarsening Option", "Recursive Coordinate Bisection");
    }
    else if (coarseningOption == 2) {
      Parameters->set("Coarsening Option", "Recursive Inertial Bisection");
    }
  }
  // solver options
  Parameters->set("Solver", "Tacho");
  switch (directSolver) {
  case 0: // SuperLU
#if defined(HAVE_SHYLU_DDBDDC_SUPERLU)
    Parameters->set("Solver", "SuperLU");
#else
    if (myPID == 0) {
      std::cout << "SuperLU solver not available in this build\n";
    }
#endif
    break;
  case 1: // MKL Pardiso solver
#if defined(HAVE_SHYLU_DDBDDC_PARDISO_MKL)
    Parameters->set("Solver", "Pardiso");
#else
    if (myPID == 0) {
      std::cout << "Pardiso solver not available\n";
    }
#endif
    break;
  case 2: // Tacho solver
#if defined(HAVE_SHYLU_DDBDDC_SHYLU_NODETACHO)
    Parameters->set("Solver", "Tacho");
#else
    std::cout << "Tacho solver not available in this build\n";
#endif
    break;
  default:
    std::cout << "Error: unavailable solver option\n";
    return 1;
    break;
  }
  return 0;
}

void printTimings
(RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner,
 RCP< bddc::KrylovSolver<SX,SM,LO,GO> > & Solver,
 const LO krylovMethod)
{
  RCP<const Teuchos::Comm<int> > Comm = Preconditioner->getComm();
  LO myPID = Comm->getRank();
  if (myPID != 0) return;
  std::ofstream timingsFile;
  timingsFile.open("timingsBDDC.dat", std::ios::out);
  const std::vector<double> & timings = Preconditioner->getTimings();
  timingsFile << "Timings in seconds for different solver components\n";
  timingsFile << "initial static condensations = "
              << timings[bddc::TIME_INIT_STATIC_COND] << std::endl;
  timingsFile << "static expansions            = "
              << timings[bddc::TIME_STATIC_EXPANSION] << std::endl;
  timingsFile << "subdomain corrections        = "
              << timings[bddc::TIME_SUB_CORR] << std::endl;
  timingsFile << "coarse corrections           = "
              << timings[bddc::TIME_COARSE_CORR] << std::endl;
  timingsFile << "subdomain corrections        = "
              << timings[bddc::TIME_SUB_CORR] << std::endl;
  timingsFile << "apply operator               = "
              << timings[bddc::TIME_APPLY_OPER] << std::endl;
  if (krylovMethod == 0) {
    timingsFile << "OthogGCR:projections         = "
                << Solver->getProjectionTime() << std::endl;
    timingsFile << "OthogGCR:orthogonalizations  = "
                << Solver->getOrthogonalizationTime() << std::endl;
  }
  timingsFile << "initialization time          = "
              << timings[bddc::TIME_INITIALIZATION] << std::endl;
  timingsFile.close();
}

void openFiles
(RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner,
 LO krylovMethod)
{
  RCP<const Teuchos::Comm<int> > Comm = Preconditioner->getComm();
  LO numProc = Comm->getSize();
  LO myPID = Comm->getRank();
  LO numSub = Preconditioner->getNumSub();
  LO numSubAll;
  Teuchos::reduceAll<int, LO> (*Comm, Teuchos::REDUCE_SUM, 1,
                               &numSub, &numSubAll);
  RCP<const Map> dofMap1to1 = Preconditioner->getDofMap1to1();
  RCP<const Map> dofMapB1to1 = Preconditioner->getDofMapB1to1();
  LO numDofs = dofMap1to1->getGlobalNumElements();
  LO numDofsB = dofMapB1to1->getGlobalNumElements();
  if (myPID == 0) {
    std::ofstream outputFileDD, outputFileKrylov;
    outputFileDD.open("dd_solver.dat", std::ios::out);
    outputFileDD << "number of processors   = " << numProc << std::endl;
    outputFileDD << "number of subdomains   = " << numSubAll << std::endl;
    outputFileDD << "number of unknowns     = " << numDofs << std::endl;
    outputFileDD << "size of interface      = " << numDofsB << std::endl;
    std::string text = "GCR";
    if (krylovMethod == 1) text = "PCG";
    outputFileDD << "Krylov Method          = " << text << std::endl;
    outputFileDD.close();
    outputFileKrylov.open("krylov_solver.dat", std::ios::out);
  }
}

} // end namespace
