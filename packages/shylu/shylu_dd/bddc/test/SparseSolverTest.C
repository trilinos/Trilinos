
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
#include "ShyLUBDDC_config.h"
#include "shylu_SolverFactoryBDDC.h"

using Teuchos::RCP;

namespace {

void readInputFile(int & matrixSource,
		   int & numElemDir1, 
		   int & numElemDir2, 
		   int & numElemDir3, 
		   int & problemTypeInt,
		   int & matrixType, 
		   double & diagScaleFactor,
		   int & numberOfSolves,
		   int & precisionOption,
		   int & directSolver);

int setParameters(int numElemDir1,
		  int numElemDir2,
		  int numElemDir3,
		  int problemTypeInt,
		  int matrixType,
		  int precisionOption,
		  int directSolver,
		  int numProc,
		  RCP<Teuchos::ParameterList> & Parameters);
  
double clockIt()
{
  struct timeval start;
  gettimeofday(&start, NULL);
  double duration = 
    (double)(start.tv_sec + start.tv_usec/1000000.0);
  return duration;
}

double calculateError(int numRows,
		      const int* rowBegin,
		      const int* columns,
		      const double* values,
		      std::vector<double> & rhs,
		      std::vector<double> & sol,
		      int numRhs);

float calculateError(int numRows,
		     const int* rowBegin,
		     const int* columns,
		     const float* values,
		     std::vector<float> & rhs,
		     std::vector<float> & sol,
		     int numRhs);

void writeTime(const char* title,
	       double time);

void readMatrix(std::vector<int> & rowBegin, 
		std::vector<int> & columns, 
		std::vector<double> & values);

TEST(SparseSolverBDDC, Test1)
{
  int numProc;
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Comm_size(Comm, &numProc);
  double diagScaleFactor;
  int numElemDir1, numElemDir2, numElemDir3, problemTypeInt, matrixType,
    numberOfSolves, precisionOption, directSolver, matrixSource;
  readInputFile(matrixSource, numElemDir1, numElemDir2, numElemDir3, 
		problemTypeInt,	matrixType, diagScaleFactor, numberOfSolves, 
		precisionOption, directSolver);
  RCP<Teuchos::ParameterList> Parameters;
  int returnVal = 
    setParameters(numElemDir1, numElemDir2, numElemDir3, problemTypeInt,
		  matrixType, precisionOption, directSolver, numProc,
		  Parameters);
  if (returnVal) {
    std::cout << "Error in setParameters\n";
    return;
  }
  RCP< bddc::ProblemMaker<int,int,double,double> > Problem = 
    rcp( new bddc::ProblemMaker<int,int,double,double>(Parameters, Comm) );
  std::vector< std::vector<int> > subNodes, subNodeBegin, subRowBegin, 
    subLocalDofs, subColumns, subElems;
  std::vector< std::vector<double> > subValues;
  int numNode(0), numRows(0), *rowBegin(0), *columns(0);
  double *values(0);
  std::vector<int> rowBeginFile, columnsFile;
  std::vector<double> valuesFile;
  if (matrixSource == 0) {
    Problem->getSubDomainElements(1, 1, 1, subElems);
    Problem->getSubDomainNodeData(subElems, subNodes, subNodeBegin,
				  subLocalDofs);
    Problem->getSubdomainMatrices(subElems, subNodes, subNodeBegin, 
				  subRowBegin, subColumns, subValues);
    Problem->addDiagonalStiffness(subRowBegin, subColumns, subValues,
				  diagScaleFactor);
    if (matrixType == 1) {
      Problem->addAsymmetry(subRowBegin, subColumns, subValues);
    }
    numNode = Problem->getNumNode();
    numRows = numNode*Parameters->get("numDofPerNode", 0);
    rowBegin = subRowBegin[0].data();
    columns = subColumns[0].data();
    values = subValues[0].data();
  }
  else {
    readMatrix(rowBeginFile, columnsFile, valuesFile);
    numRows = rowBeginFile.size() - 1;
    rowBegin = rowBeginFile.data();
    columns = columnsFile.data();
    values = valuesFile.data();
  }
  int numRhs(1);
  std::vector<double> rhs(numRows*numRhs), sol(numRows*numRhs);
  for (int i=0; i<numRows*numRhs; i++) {
    rhs[i] = 0.7*rand()/RAND_MAX;
    rhs[i] = 1;
  }
  std::vector<float> rhsF, solF, valuesF;
  if (precisionOption == 1) {
    int numTerms = rowBegin[numRows];
    valuesF.resize(numTerms);
    for (int i=0; i<numTerms; i++) valuesF[i] = values[i];
    solF.resize(numRows);
    rhsF.resize(numRows);
    for (int i=0; i<numRows; i++) rhsF[i] = rhs[i];
  }
  bddc::SolverFactory<double> Factory;
  bddc::SolverBase<double>* Solver(0);
  bddc::SolverFactory<float> FactoryF;
  bddc::SolverBase<float>* SolverF(0);
  if (precisionOption == 0) {
    double startTime = clockIt();
    Solver = 
      Factory.Generate(numRows, rowBegin, columns, values, *Parameters);
    Solver->Initialize();
    double dt = clockIt() - startTime;
    writeTime("Solver(double) initialization time = ", dt);
    startTime = clockIt();
    for (int i=0; i<numberOfSolves; i++) {
      Solver->Solve(numRhs, rhs.data(), sol.data());
    }
    dt = clockIt() - startTime;
    writeTime("SuperLU(double) solve time = ", dt);
    double solError = calculateError(numRows, rowBegin, columns, values, 
				     rhs, sol, numRhs);
    std::cout << "solError = " << solError << std::endl;
    EXPECT_LT(solError, 1e-10);
    delete Solver;
  }
  else {
    double startTime = clockIt();
    SolverF = 
      FactoryF.Generate(numRows, rowBegin, columns, valuesF.data(), 
			*Parameters);
    SolverF->Initialize();
    double dt = clockIt() - startTime;
    writeTime("Solver(single) initialization time = ", dt);
    startTime = clockIt();
    for (int i=0; i<numberOfSolves; i++) {
      SolverF->Solve(numRhs, rhsF.data(), solF.data());
    }
    dt = clockIt() - startTime;
    writeTime("Solver(single) solve time = ", dt);
    float solError = calculateError(numRows, rowBegin, columns, 
				    valuesF.data(), rhsF, solF, numRhs);
    std::cout << "solError = " << solError << std::endl;
    EXPECT_LT(solError, 1e-5);
    delete SolverF;
  }
}

void writeTime(const char* title,
	       double time)
{
  std::cout << title << time << std::endl;
}

void readInputFile(int & matrixSource,
		   int & numElemDir1, 
		   int & numElemDir2, 
		   int & numElemDir3, 
		   int & problemTypeInt,
		   int & matrixType, 
		   double & diagScaleFactor,
		   int & numberOfSolves,
		   int & precisionOption,
		   int & directSolver)
{
  char buff[101];
  std::ifstream fin;
  fin.open("SparseSolverTest.inp");
  fin >> matrixSource; fin.getline(buff,100);
  fin >> numElemDir1; fin.getline(buff,100);
  fin >> numElemDir2; fin.getline(buff,100);
  fin >> numElemDir3; fin.getline(buff,100);
  fin >> problemTypeInt; fin.getline(buff,100);
  fin >> matrixType; fin.getline(buff,100);
  fin >> diagScaleFactor; fin.getline(buff,100);
  fin >> numberOfSolves; fin.getline(buff,100);
  fin >> precisionOption; fin.getline(buff,100);
  fin >> directSolver; fin.getline(buff,100);
  fin.close();
}

int setParameters(int numElemDir1,
		  int numElemDir2,
		  int numElemDir3,
		  int problemTypeInt,
		  int matrixType,
		  int precisionOption,
		  int directSolver,
		  int numProc,
		  RCP<Teuchos::ParameterList> & Parameters)
{
  int spatialDim(3);
  enum bddc::ProblemType problemType = bddc::SCALARPDE;
  if (problemTypeInt == 2) problemType = bddc::ELASTICITY;
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  int numDofPerNode(0), loadDirection(0);
  numDofPerNode = 1;
  if (problemType == bddc::ELASTICITY) {
    numDofPerNode = spatialDim;
  }
  //  
  Parameters = Teuchos::rcp( new Teuchos::ParameterList() );
  Parameters->set("numDofPerNode", numDofPerNode);
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
  double lengthDir1(numProc), lengthDir2(1), lengthDir3(1);
  int numSubDir1(numProc), numSubDir2(1), numSubDir3(1);
  Parameters->set("Length Direction 1", lengthDir1);
  Parameters->set("Length Direction 2", lengthDir2);
  Parameters->set("Length Direction 3", lengthDir3);
  Parameters->set("Number of Subdomains Direction 1", numSubDir1);
  Parameters->set("Number of Subdomains Direction 2", numSubDir2);
  Parameters->set("Number of Subdomains Direction 3", numSubDir3);
  Parameters->set("Number of Elements Per Subdomain Direction 1",
		  numElemDir1);
  Parameters->set("Number of Elements Per Subdomain Direction 2",
		  numElemDir2);
  Parameters->set("Number of Elements Per Subdomain Direction 3",
		  numElemDir3);
  Parameters->set("Apply Left Side Essential BCs", false);
  Parameters->set("Apply Right Side Essential BCs", false);
  Parameters->set("Load Direction", loadDirection);
  Parameters->set("Artificial Foundation Stiffness", 0.0);
  Parameters->set("omega", 0.0);
  Parameters->set("Generate Constraint Equations", false);
  if (precisionOption == 1) {
    Parameters->set("Precision", "single");
  }
  else {
    Parameters->set("Precision", "double");
  }
  switch (directSolver) {
  case 0: // SuperLU
#if defined(HAVE_SHYLUBDDC_SUPERLU)
    Parameters->set("Solver", "SuperLU");
#else
    std::cout << "SuperLU solver not available in this build\n";
#endif
    break;
  case 1: // MKL Pardiso solver
#if defined(HAVE_SHYLUBDDC_PARDISO_MKL)
    Parameters->set("Solver", "Pardiso");
#else
    std::cout << "Pardiso solver not available in this build\n";
#endif
    break;
  case 2: // Tacho solver
#if defined(HAVE_SHYLUBDDC_SHYLUTACHO)
    Parameters->set("Solver", "Tacho");
#else
    std::cout << "Tacho solver not available in this build\n";
#endif
    break;
  default:
    std::cout << "Error: unavailable direct solver option\n";
    return 1;
    break;
  }
  return 0;
}

double calculateError(int numRows,
		      const int* rowBegin,
		      const int* columns,
		      const double* values,
		      std::vector<double> & rhs,
		      std::vector<double> & sol,
		      int numRhs)
{
  double maxError(0);
  for (int k=0; k<numRhs; k++) {
    double* RHS = &rhs[k*numRows];
    double* SOL = &sol[k*numRows];
    for (int i=0; i<numRows; i++) {
      double resid = RHS[i];
      for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	resid -= values[j]*SOL[columns[j]];
      }
      resid = std::abs(resid);
      if (resid > maxError) maxError = resid;
    }
  }
  return maxError;
}

float calculateError(int numRows,
		     const int* rowBegin,
		     const int* columns,
		     const float* values,
		     std::vector<float> & rhs,
		     std::vector<float> & sol,
		     int numRhs)
{
  float maxError(0);
  for (int k=0; k<numRhs; k++) {
    float* RHS = &rhs[k*numRows];
    float* SOL = &sol[k*numRows];
    for (int i=0; i<numRows; i++) {
      float resid = RHS[i];
      for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	resid -= values[j]*SOL[columns[j]];
      }
      resid = std::abs(resid);
      if (resid > maxError) maxError = resid;
    }
  }
  return maxError;
}

void readMatrix(std::vector<int> & rowBegin, 
		std::vector<int> & columns, 
		std::vector<double> & values)
{  
  std::ifstream fin;
  std::string matrixFileName = "matrix.dat";
  fin.open(matrixFileName);
  int row, col, rowMax(0), colMax(0), numTerms(0);
  double value;
  fin >> row >> col >> value;
  while (!fin.eof()) {
    if (row > rowMax) rowMax = row;
    if (col > colMax) colMax = col;
    fin >> row >> col >> value;
    numTerms++;
  }
  assert (rowMax == colMax);
  int numDof = rowMax;
  fin.close();
  fin.open(matrixFileName);
  std::vector<int> count(numDof, 0);
  for (int i=0; i<numTerms; i++) {
    fin >> row >> col >> value;
    count[row-1]++;
  }
  rowBegin.resize(numDof+1, 0);
  for (int i=0; i<numDof; i++) {
    rowBegin[i+1] = rowBegin[i] + count[i];
  }
  fin.close();
  assert (numTerms == rowBegin[numDof]);
  columns.resize(numTerms);
  values.resize(numTerms);
  count.assign(numDof, 0);
  fin.open(matrixFileName);
  for (int i=0; i<numTerms; i++) {
    fin >> row >> col >> value;
    row--; col--;
    int index = rowBegin[row] + count[row];
    columns[index] = col;
    values[index] = value;
    count[row]++;
  }
  fin.close();  
}

} // end namespace
