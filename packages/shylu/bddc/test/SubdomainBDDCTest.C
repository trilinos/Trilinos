
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
#include "shylu_SubdomainBDDC.h"
#include "shylu_UtilBDDC.h"
#include "shylu_enumsBDDC.h"

using Teuchos::RCP;

namespace {

typedef int LO; // Local Ordinal
typedef int GO; // Global Ordinal
typedef double SX; // floating point data type
typedef double SM; // real (magnitude) for SX
    
void getBoundaryDofs(const LO spatialDim,
		     const LO numNode,
		     const LO numDofPerNode,
		     const SM* xCoord,
		     const SM* yCoord,
		     const SM* zCoord,
		     const SM lengthDir1,
		     const SM lengthDir2,
		     const SM lengthDir3,
		     std::vector<LO> & boundaryDofs);
void getEquivClasses(const LO spatialDim,
		     const LO numNode,
		     const LO numDofPerNode,
		     const SM* xCoord,
		     const SM* yCoord,
		     const SM* zCoord,
		     const SM lengthDir1,
		     const SM lengthDir2,
		     const SM lengthDir3,
		     std::vector<LO> & boundaryDofs,
		     std::vector< std::vector<LO> > & equivClassDofs);
void getEquivClassConstraints
  (const LO spatialDim,
   const LO numDofPerNode,
   const bool useCorners,
   const bool useEdges,
   const bool useFaces,
   std::vector< std::vector<LO> > & equivClassDofs,
   std::vector< std::vector<LO> > & equivConstraintsLocalDofs,
   std::vector< std::vector<SX> > & equivConstraints);
void addConstraints(const LO numNodes,
		    const LO numDofPerNode,
		    std::vector<LO> & constraintsLocalDofs,
		    std::vector<SX> & constraints);
void getAuxiliaryConstraints
  (std::vector< std::vector<LO> > & equivClassDofs,
   const LO numDofPerNode,
   const SM* xCoord,
   const SM* yCoord,
   const SM* zCoord,
   std::vector< std::vector<LO> > & auxConstraintsLocalDofs, 
   std::vector< std::vector<SX> > & auxConstraints);
void getEquivalenceClassWeights
            (std::vector< std::vector<LO> > & equivClassDofs, 
	     std::vector< std::vector<LO> > & rowBeginWeight, 
	     std::vector< std::vector<LO> > & columnsWeight,
	     std::vector< std::vector<SX> > & valuesWeight);
void printMatrix(const LO numRows,
		 const LO numCols,
		 const SX* A,
		 const char* fileName);
void checkSymmetry(std::vector<SX> & A,
		   LO numRows,
		   const char* label);
void checkCoarseMatrices(std::vector<SX> & Phi, 
			 std::vector<SX> & Ac, 
			 LO numRows,
			 LO numCols,
			 const char* label);
void determineBoundingECs(LO spatialDim,
			  LO & equivClass, 
			  std::vector<LO> & boundingECs);
void checkNullSpace(std::vector<SX> & A,
		    LO numRows,
		    const char* label);

TEST(SubdomainBDDC, Test1)
{
  int numProc, myPID;
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Comm_size(Comm, &numProc);
  MPI_Comm_rank(Comm, &myPID);
  if (numProc != 1) return;
  double lengthDir1 = 1;
  double lengthDir2 = 1;
  double lengthDir3 = 1;
  LO numSubDir1 = 1;
  LO numSubDir2 = 1;
  LO numSubDir3 = 1;
  LO numElemPerSubDir1 = 4;
  LO numElemPerSubDir2 = 4;
  LO numElemPerSubDir3 = 4;
  
  enum bddc::ProblemType problemType = bddc::SCALARPDE; 
  // ELASTICITY of SCALARPDE
  int spatialDim(2), numDofPerNode(0), loadDirection(0);
  if (problemType == bddc::SCALARPDE) numDofPerNode = 1;
  if (problemType == bddc::ELASTICITY) numDofPerNode = spatialDim;
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  RCP<Teuchos::ParameterList> Parameters;
  Parameters = Teuchos::rcp( new Teuchos::ParameterList() );
  Parameters->set("Spatial Dimension", spatialDim);
  Parameters->set("Problem Type", problemType);
  Parameters->set("Analysis Type", analysisType);
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
  
  Parameters->set("Apply Left Side Essential BCs", false);
  Parameters->set("Apply Right Side Essential BCs", false);

  Parameters->set("Load Direction", loadDirection);
  Parameters->set("Artificial Foundation Stiffness", 0);
  Parameters->set("omega", 0.0);
  Parameters->set("Generate Constraint Equations", false);
  Parameters->set("Interface Preconditioner", true);
  Parameters->set("Print Interior Matrices", false);

  assert (numProc == numSubDir1*numSubDir2*numSubDir3);
  //
  // generate problem
  //
  RCP< bddc::ProblemMaker<LO,GO,SX,SM> > Problem = 
    rcp( new bddc::ProblemMaker<LO,GO,SX,SM>(Parameters, Comm) );
  LO numNode = Problem->getNumNode();
  const LO* nodeBegin = Problem->getNodeBegin();
  const LO* localDofs = Problem->getLocalDofs();
  //  const GO* nodeGIDs = Problem->getNodeGlobalIDs();
  const LO* rowBegin = Problem->getRowBegin();
  const LO* columns = Problem->getColumns();
  const SX* values = Problem->getValues();
  //
  // subdomain object
  //
  std::vector<LO> boundaryDofs;
  getBoundaryDofs(spatialDim, numNode, numDofPerNode, Problem->getXcoord(),
		  Problem->getYcoord(), Problem->getZcoord(),
		  lengthDir1, lengthDir2, lengthDir3, boundaryDofs);
  RCP< bddc::SubdomainBDDC<SX,SM,LO,GO> > subdomainBDDC = 
    rcp( new bddc::SubdomainBDDC<SX,SM,LO,GO>
	 (numNode, nodeBegin, localDofs, rowBegin, columns, values,
	  boundaryDofs.size(), &boundaryDofs[0], *Parameters) );
  std::vector< std::vector<LO> > equivClassDofs;
  getEquivClasses(spatialDim, numNode, numDofPerNode, Problem->getXcoord(),
		  Problem->getYcoord(), Problem->getZcoord(),
		  lengthDir1, lengthDir2, lengthDir3, boundaryDofs,
		  equivClassDofs);
  subdomainBDDC->setEquivalenceClasses(equivClassDofs);
  // original constrains
  bool useCorners(true), useEdges(true), useFaces(false);
  std::vector< std::vector<LO> > equivConstraintsLocalDofs;
  std::vector< std::vector<SX> > equivConstraints;
  getEquivClassConstraints(spatialDim, numDofPerNode, useCorners, useEdges,
			   useFaces, equivClassDofs, 
			   equivConstraintsLocalDofs, equivConstraints);
  subdomainBDDC->setInitialConstraints(equivConstraintsLocalDofs, 
				       equivConstraints);
  std::vector<SX> Phi, Ac;
  subdomainBDDC->calculateCoarseMatrices(Ac);
  subdomainBDDC->getPhi(Phi);
  LO numRows = subdomainBDDC->getNumDofs();
  LO numCols = subdomainBDDC->getCoarseSpaceDimension();
  //  printMatrix(numCols, numCols, &Ac[0], "AcOrig.dat");
  checkCoarseMatrices(Phi, Ac, numRows, numCols, "Coarse Matrices Orig");
  // auxiliary constraints
  std::vector< std::vector<LO> > auxConstraintsLocalDofs;
  std::vector< std::vector<SX> > auxConstraints;
  getAuxiliaryConstraints
    (equivClassDofs, numDofPerNode, 
     Problem->getXcoord(), Problem->getYcoord(), Problem->getZcoord(),
     auxConstraintsLocalDofs, auxConstraints);
  subdomainBDDC->setAuxiliaryConstraints(auxConstraintsLocalDofs,
					 auxConstraints);
  subdomainBDDC->calculateCoarseMatrices(Ac);
  subdomainBDDC->getPhi(Phi);
  numCols = subdomainBDDC->getCoarseSpaceDimension();
  //  printMatrix(numCols, numCols, &Ac[0], "AcAux.dat");
  checkCoarseMatrices(Phi, Ac, numRows, numCols, "Coarse Matrices Aux");
  // weight matrices
  std::vector< std::vector<LO> > rowBeginWeight, columnsWeight;
  std::vector< std::vector<SX> > valuesWeight;
  getEquivalenceClassWeights(equivClassDofs, rowBeginWeight, columnsWeight,
			     valuesWeight);
  subdomainBDDC->setEquivalenceClassWeightMatrices
    (rowBeginWeight, columnsWeight, valuesWeight);
  LO numB = boundaryDofs.size();
  std::vector<SX> MinvSub(numB*numB), g(numB, 0);
  for (LO i=0; i<numB; i++) {
    g[i] = 1;
    subdomainBDDC->applyNeumannCorrection(&g[0], &MinvSub[i*numB]);
    g[i] = 0;
  }
  checkSymmetry(MinvSub, numB, "Neumann Correction");
  // Schur complement for equivalence class and "bounding" equivalence classes
  LO equivClass;
  std::vector<LO> boundingECs;
  determineBoundingECs(spatialDim, equivClass, boundingECs);
  std::vector<SX> ScEquivBound;
  int numRowsBound;
  subdomainBDDC->calculateBoundingSchurComplement
    (equivClass, boundingECs, ScEquivBound, numRows, numRowsBound);
  //  printMatrix(numRows, numRows, &ScEquivBound[0], "ScEquivBound.dat");
  checkSymmetry(ScEquivBound, numRows, "ScEquivBound");
  checkNullSpace(ScEquivBound, numRows, "ScEquivBound");
  // basic Schur complement for equivalence class
  std::vector<SX> ScEquiv;
  subdomainBDDC->calculateSchurComplement(equivClass, ScEquiv);
  numRows = LO(std::sqrt(ScEquiv.size()));
  //  printMatrix(numRows, numRows, &ScEquiv[0], "ScEquiv.dat");
  checkSymmetry(ScEquiv, numRows, "ScEquiv");
  // full Schur complement
  std::vector<SX> ScAll(numB*numB), xB(numB);
  SX ALPHA(1), BETA(0);
  for (LO i=0; i<numB; i++) {
    xB[i] = 1;
    subdomainBDDC->applyBoundaryOperator(&xB[0], &ScAll[i*numB], ALPHA, BETA);
    xB[i] = 0;
  }
  checkSymmetry(ScAll, numB, "ScAll");
}

void checkNullSpace(std::vector<SX> & A,
		    LO numRows,
		    const char* label)
{
  SM tol(1e-10);
  for (LO i=0; i<numRows; i++) {
    SX sum(0);
    for (LO j=0; j<numRows; j++) {
      sum += A[i+j*numRows];
    }
    EXPECT_LT(std::abs(sum), tol) << label << " Null Space";
  }

}

void determineBoundingECs(LO spatialDim,
			  LO & equivClass, 
			  std::vector<LO> & boundingECs)
{
  if (spatialDim == 2) {
    // equivalence class is edge on the left side
    equivClass = 7;
    boundingECs.resize(0);
    boundingECs.push_back(0); // lower left corner
    boundingECs.push_back(3); // upper left corner
  }
}

void checkCoarseMatrices(std::vector<SX> & Phi, 
			 std::vector<SX> & Ac, 
			 LO numRows,
			 LO numCols,
			 const char* label)
{
  SM tol(1e-10);
  for (LO i=0; i<numRows; i++) {
    SX sum(0), ONE(1);
    for (LO j=0; j<numCols; j++) {
      sum += Phi[i+j*numRows];
    }
    SM diff = std::abs(sum - ONE);
    EXPECT_LT(diff, tol) << label << " Partition of Unity";
  }
  checkSymmetry(Ac, numCols, label);
  checkNullSpace(Ac, numCols, label);
}

void checkSymmetry(std::vector<SX> & A,
		   LO numRows,
		   const char* label)
{
  SM tol(1e-10), maxAbsDiag(0);
  for (LO i=0; i<numRows; i++) {
    SM absDiag = std::abs(A[i+i*numRows]);
    if (absDiag > maxAbsDiag) maxAbsDiag = absDiag;
  }
  assert (maxAbsDiag > 0);
  for (LO i=0; i<numRows; i++) {
    for (LO j=0; j<i; j++) {
      SM diff = std::abs(A[i+j*numRows] - A[j+i*numRows]);
      SM error = diff/maxAbsDiag;
      EXPECT_LT(error, tol) << label;
    }
  }
}

void getEquivalenceClassWeights
            (std::vector< std::vector<LO> > & equivClassDofs, 
	     std::vector< std::vector<LO> > & rowBeginWeight, 
	     std::vector< std::vector<LO> > & columnsWeight,
	     std::vector< std::vector<SX> > & valuesWeight)
{
  LO numEquiv = equivClassDofs.size();
  rowBeginWeight.resize(numEquiv);
  columnsWeight.resize(numEquiv);
  valuesWeight.resize(numEquiv);
  for (LO i=0; i<numEquiv; i++) {
    LO numEquivDofs = equivClassDofs[i].size();
    LO numTerms = numEquivDofs*numEquivDofs;
    rowBeginWeight[i].resize(numEquivDofs+1, 0);
    columnsWeight[i].resize(numTerms);
    valuesWeight[i].resize(numTerms);
    numTerms = 0;
    for (LO j=0; j<numEquivDofs; j++) {
      for (LO k=0; k<numEquivDofs; k++) {
	columnsWeight[i][numTerms] = k;
	valuesWeight[i][numTerms] = 0.7*rand()/RAND_MAX;
	numTerms++;
      }
      rowBeginWeight[i][j+1] = numTerms;
    }
  }
}

void getAuxiliaryConstraints
  (std::vector< std::vector<LO> > & equivClassDofs,
   const LO numDofPerNode,
   const SM* xCoord,
   const SM* yCoord,
   const SM* zCoord,
   std::vector< std::vector<LO> > & auxConstraintsLocalDofs, 
   std::vector< std::vector<SX> > & auxConstraints)
{
  LO numEquiv = equivClassDofs.size();
  auxConstraintsLocalDofs.resize(numEquiv);
  auxConstraints.resize(numEquiv);
  for (LO i=0; i<numEquiv; i++) {
    LO numDofs = equivClassDofs[i].size();
    LO numNodes = numDofs/numDofPerNode;
    if (numNodes >= 3) {
      SM xCent(0), yCent(0), zCent(0);
      for (LO j=0; j<numDofs; j++) {
	LO node = equivClassDofs[i][j]/numDofPerNode;
	xCent += xCoord[node];
	yCent += yCoord[node];
	zCent += zCoord[node];
      }
      xCent /= numDofs;
      yCent /= numDofs;
      zCent /= numDofs;
      bool centerNodePresent(false);
      LO centerNode = -1;
      for (LO j=0; j<numNodes; j++) {
	LO node = equivClassDofs[i][numDofPerNode*j]/numDofPerNode;
	if ((std::abs(xCoord[node]-xCent) < 1e-6) &&
	    (std::abs(yCoord[node]-yCent) < 1e-6) &&
	    (std::abs(zCoord[node]-zCent) < 1e-6)) {
	  centerNodePresent = true;
	  centerNode = j;
	  break;
	}
      }
      if (centerNodePresent == true) {
	auxConstraintsLocalDofs[i].resize(numDofPerNode);
	auxConstraints[i].resize(numDofPerNode*numDofs);
	for (LO j=0; j<numDofPerNode; j++) {
	  auxConstraintsLocalDofs[i][j] = j;
	  LO index = numDofs*j + numDofPerNode*centerNode + j;
	  auxConstraints[i][index] = 1;
	}
      }
    }
  }
}

void printMatrix(const LO numRows,
		 const LO numCols,
		 const SX* A,
		 const char* fileName)
{
  std::ofstream fout;
  fout.open(fileName);
  SX value(0);
  bool isComplex = bddc::UtilBDDC<SX,SM>::isComplex(value);
  for (LO j=0; j<numCols; j++) {
    for (LO i=0; i<numRows; i++) {
      fout << i+1 << " ";
      fout << j+1 << " ";
      fout << std::setw(22) << std::setprecision(15);
      value = A[i+numRows*j];
      fout << bddc::UtilBDDC<SX,SM>::real(value);
      if (isComplex == true) {
	fout << " " << bddc::UtilBDDC<SX,SM>::imag(value);
      }
      fout << std::endl;
    }
  }
  fout.close();
}

void getEquivClasses(const LO spatialDim,
		     const LO numNode,
		     const LO numDofPerNode,
		     const SM* xCoord,
		     const SM* yCoord,
		     const SM* zCoord,
		     const SM lengthDir1,
		     const SM lengthDir2,
		     const SM lengthDir3,
		     std::vector<LO> & boundaryDofs,
		     std::vector< std::vector<LO> > & equivClassDofs)
{
  LO numBoundaryNodes = boundaryDofs.size()/numDofPerNode;
  std::vector<LO> boundaryNodes(numBoundaryNodes);
  for (LO i=0; i<numBoundaryNodes; i++) {
    boundaryNodes[i] = boundaryDofs[numDofPerNode*i]/numDofPerNode;
  }
  SM tol = 1e-6;
  if (spatialDim == 2) equivClassDofs.resize(8);
  if (spatialDim == 3) equivClassDofs.resize(26);
  for (LO i=0; i<numBoundaryNodes; i++) {
    LO node = boundaryNodes[i];
    SM dx1 = std::abs(xCoord[node] - 0);
    SM dx2 = std::abs(xCoord[node] - lengthDir1);
    int ix = 1;
    if (dx1 < tol*lengthDir1) ix = 0;
    if (dx2 < tol*lengthDir1) ix = 2;
    SM dy1 = std::abs(yCoord[node] - 0);
    SM dy2 = std::abs(yCoord[node] - lengthDir2);
    int iy = 1;
    if (dy1 < tol*lengthDir2) iy = 0;
    if (dy2 < tol*lengthDir2) iy = 2;
    int ec = -1;
    if (spatialDim == 2) {
      if ((ix == 0) && (iy == 0)) ec = 0;
      if ((ix == 2) && (iy == 0)) ec = 1;
      if ((ix == 2) && (iy == 2)) ec = 2;
      if ((ix == 0) && (iy == 2)) ec = 3;
      if ((ix == 1) && (iy == 0)) ec = 4;
      if ((ix == 2) && (iy == 1)) ec = 5;
      if ((ix == 1) && (iy == 2)) ec = 6;
      if ((ix == 0) && (iy == 1)) ec = 7;
    }
    if (ec != -1) {
      for (LO j=0; j<numDofPerNode; j++) {
	equivClassDofs[ec].push_back(numDofPerNode*node+j);
      }
    }
  }
}

void getEquivClassConstraints
  (const LO spatialDim,
   const LO numDofPerNode,
   const bool useCorners,
   const bool useEdges,
   const bool useFaces,
   std::vector< std::vector<LO> > & equivClassDofs,
   std::vector< std::vector<LO> > & equivConstraintsLocalDofs,
   std::vector< std::vector<SX> > & equivConstraints)
{
  size_t numEquiv = equivClassDofs.size();
  equivConstraintsLocalDofs.resize(numEquiv);
  equivConstraints.resize(numEquiv);
  for (size_t i=0; i<numEquiv; i++) {
    if (spatialDim == 2) {
      LO numNodes = equivClassDofs[i].size()/numDofPerNode;
      bool activateConstraints = false;
      if ((numNodes == 1) && (useCorners == true)) activateConstraints = true;
      if ((numNodes > 1) && (useEdges == true)) activateConstraints = true;
      if (activateConstraints == true) {
	addConstraints(numNodes, numDofPerNode, equivConstraintsLocalDofs[i], 
		       equivConstraints[i]);
      }
    }
  }
}

void addConstraints(const LO numNodes,
		    const LO numDofPerNode,
		    std::vector<LO> & constraintsLocalDofs,
		    std::vector<SX> & constraints)
{
  LO numDof = numNodes*numDofPerNode;
  constraints.resize(numDof*numDofPerNode);
  constraintsLocalDofs.resize(numDofPerNode);
  for (LO i=0; i<numDofPerNode; i++) {
    constraintsLocalDofs[i] = i;
    LO begin = numDof*i;
    for (LO j=0; j<numNodes; j++) {
      constraints[begin+numDofPerNode*j+i] = SX(1)/numNodes;
    }
  }
}

void getBoundaryDofs(const LO spatialDim,
		     const LO numNode,
		     const LO numDofPerNode,
		     const SM* xCoord,
		     const SM* yCoord,
		     const SM* zCoord,
		     const SM lengthDir1,
		     const SM lengthDir2,
		     const SM lengthDir3,
		     std::vector<LO> & boundaryDofs)
{
  for (LO i=0; i<numNode; i++) {
    SM dx1 = std::abs(xCoord[i] - 0);
    SM dx2 = std::abs(xCoord[i] - lengthDir1);
    SM dx = std::min(dx1, dx2);
    SM dy1 = std::abs(yCoord[i] - 0);
    SM dy2 = std::abs(yCoord[i] - lengthDir2);
    SM dy = std::min(dy1, dy2);
    bool onBoundary = false;
    SM tol(1e-6);
    if ((dx < tol*lengthDir1) || (dy < tol*lengthDir2)) {
      onBoundary = true;
    }
    if (spatialDim == 3) {
      SM dz1 = std::abs(zCoord[i] - 0);
      SM dz2 = std::abs(zCoord[i] - lengthDir3);
      SM dz = std::min(dz1, dz2);
      if (dz < tol*lengthDir3) onBoundary = true;
    }
    if (onBoundary == true) {
      for (int j=0; j<numDofPerNode; j++) {
	boundaryDofs.push_back(numDofPerNode*i+j);
      }
    }
  }
}

} // end namespace
