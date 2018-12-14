
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
#include "ShyLU_DDBDDC_config.h"
#include "shylu_SolverFactoryBDDC.hpp"
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

using Teuchos::RCP;

namespace {

int setParameters(RCP<Teuchos::ParameterList> & params);

void printMatrix(const int numRows, 
		 const int* rowBegin, 
		 const int* columns,
		 const double* values, 
		 const char* fileName)
{
  std::ofstream fout;
  fout.open(fileName);
  for (int i=0; i<numRows; i++) {
    for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
      const int col = columns[j];
      fout << i+1 << " " << col+1 << " ";
      fout << std::setw(23) << std::setprecision(16);
      fout << values[j] << '\n';
    }
  }
  fout.close();
}
  
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

TEST(SparseSolverBDDC, Test1)
{
  const std::string fileName = "sparseSolverTest.xml";
  MPI_Comm Comm = MPI_COMM_WORLD;
  RCP<const Teuchos::Comm<int> > TComm = 
    Teuchos::rcp( new Teuchos::MpiComm<int>(Comm) );
  Teuchos::ParameterList parameters;
  Teuchos::updateParametersFromXmlFileAndBroadcast
    (fileName, Teuchos::Ptr<Teuchos::ParameterList>(&parameters), *TComm);
  RCP<Teuchos::ParameterList> params = 
    Teuchos::rcp( new Teuchos::ParameterList(parameters) );
  int myPID;
  MPI_Comm_rank(Comm, &myPID);
  if (myPID > 0) return;
  int returnVal = setParameters(params);
  if (returnVal) {
    std::cout << "Error in setParameters\n";
    return;
  }
  MPI_Comm Comm2 = MPI_COMM_SELF;
  RCP< bddc::ProblemMaker<int,int,double,double> > Problem = 
    rcp( new bddc::ProblemMaker<int,int,double,double>(params, Comm2) );
  std::vector< std::vector<int> > subNodes, subNodeBegin, subRowBegin, 
    subLocalDofs, subColumns;
  std::vector< std::vector<double> > subValues;
  int *rowBegin(nullptr), *columns(nullptr);
  double *values(nullptr);
  std::vector<int> rowBeginFile, columnsFile;
  std::vector<double> valuesFile;
  Problem->getSubDomainNodeData(subNodes, subNodeBegin, subLocalDofs);
  Problem->getSubdomainMatrices(subNodes, subNodeBegin, 
				subRowBegin, subColumns, subValues);
  if (params->get("Matrix Type", "Symmetric") == "NonSymmetric") {
    Problem->addAsymmetry(subRowBegin, subColumns, subValues);
  }
  const int numRows = Problem->getNumDof();
  rowBegin = subRowBegin[0].data();
  columns = subColumns[0].data();
  values = subValues[0].data();
  if (params->get("Print Matrix", false)) {
    printMatrix(numRows, rowBegin, columns, values, "A.dat");
  }
  int numRhs(1);
  std::vector<double> rhs(numRows*numRhs), sol(numRows*numRhs);
  for (int i=0; i<numRows*numRhs; i++) {
    rhs[i] = 0.7*rand()/RAND_MAX;
    //    rhs[i] = 1;
  }
  std::vector<float> rhsF, solF, valuesF;
  if (params->get("Precision", "double") == "single") {
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
  const int numberOfSolves = params->get("Number of Solves", 100);
  std::cout << "Solver = " << params->get("Solver", "SuperLU") << std::endl;
  std::cout << "number of unknowns = " << numRows << std::endl;
  if (params->get("Precision", "double") == "double") {
    double startTime = clockIt();
    Solver = 
      Factory.Generate(numRows, rowBegin, columns, values, *params);
    Solver->Initialize();
    double dt = clockIt() - startTime;
    writeTime("Solver(double) initialization time = ", dt);
    startTime = clockIt();
    for (int i=0; i<numberOfSolves; i++) {
      Solver->Solve(numRhs, rhs.data(), sol.data());
    }
    dt = clockIt() - startTime;
    writeTime("Solver(double) solve time = ", dt);
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
			*params);
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
    EXPECT_LT(solError, 1e-4);
    delete SolverF;
  }
}

void writeTime(const char* title,
	       double time)
{
  std::cout << title << time << std::endl;
}

int setParameters(RCP<Teuchos::ParameterList> & params)
{
  int numElemDir1 = params->get("numElemDir1", 7);
  int numElemDir2 = params->get("numElemDir2", 7);
  int numElemDir3 = params->get("numElemDir3", 7);
  const std::string precisionOption = params->get("Precision Option", "double");
  const std::string problemTypeString = 
    params->get("Problem Type String", "Poisson-3D");
  enum bddc::ProblemType problemType = bddc::SCALARPDE;
  if ((problemTypeString == "Elasticity-2D") ||
      (problemTypeString == "Elasticity-3D")) problemType = bddc::ELASTICITY;
  int spatialDim = 3;
  if ((problemTypeString == "Poisson-2D") || 
      (problemTypeString == "Elasticity-2D")) spatialDim = 2;  
  params->set("Problem Type", problemType);
  params->set("Spatial Dimension", spatialDim);
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  int numDofPerNode(0), loadDirection(0);
  numDofPerNode = 1;
  if (problemType == bddc::ELASTICITY) {
    numDofPerNode = spatialDim;
  }
  //  
  params->set("numDofPerNode", numDofPerNode);
  params->set("Analysis Type", analysisType);
  double lengthDir1(1), lengthDir2(1), lengthDir3(1);
  int numSubDir1(1), numSubDir2(1), numSubDir3(1);
  params->set("Length Direction 1", lengthDir1);
  params->set("Length Direction 2", lengthDir2);
  params->set("Length Direction 3", lengthDir3);
  params->set("Number of Subdomains Direction 1", numSubDir1);
  params->set("Number of Subdomains Direction 2", numSubDir2);
  params->set("Number of Subdomains Direction 3", numSubDir3);
  params->set("Number of Elements Per Subregion Direction 1", numElemDir1);
  params->set("Number of Elements Per Subregion Direction 2", numElemDir2);
  params->set("Number of Elements Per Subregion Direction 3", numElemDir3);
  params->set("Apply Left Side Essential BCs", true);
  params->set("Apply Right Side Essential BCs", false);
  params->set("Load Direction", loadDirection);
  params->set("Artificial Foundation Stiffness", 0.0);
  params->set("omega", 0.0);
  params->set("Generate Constraint Equations", false);
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

} // end namespace
