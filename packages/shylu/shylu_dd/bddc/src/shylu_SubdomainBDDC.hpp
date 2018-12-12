
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

#ifndef BDDC_SUBDOMAIN_H
#define BDDC_SUBDOMAIN_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>

#include "shylu_SolverFactoryBDDC.hpp"
#include "Teuchos_ParameterList.hpp"  
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "shylu_UtilBDDC.hpp"
#include "shylu_enumsBDDC.hpp"
#include "shylu_errorBDDC.hpp"

#if defined(HAVE_SHYLU_DDBDDC_MUELU)
#include "shylu_SolverMueLu.hpp"
#endif

namespace bddc {
  
template <class SX, class SM, class LO, class GO> class SubdomainBDDC
{
public:
  SubdomainBDDC();
  SubdomainBDDC(LO numNodes,
		const LO* nodes,
		const LO* nodeBegin,
		const LO* localDofs,
		const LO* rowBegin,
		const LO* columns,
		const SX* values,
		const SM* xCoord,
		const SM* yCoord,
		const SM* zCoord,
		LO numBoundaryDofs,
		const LO* boundaryDofs,
		Teuchos::ParameterList & Parameters) :
  m_numDofs(nodeBegin[numNodes]),
    m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_rowBegin(rowBegin),
    m_columns(columns),
    m_values(values),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_numBoundaryDofs(numBoundaryDofs),
    m_numAuxConstraints(0),
    m_boundaryDofs(boundaryDofs),
    m_interfacePreconditioner
    (Parameters.get("Interface Preconditioner", false)),
    m_useBlockPreconditioner(Parameters.get("Use Block Preconditioner", false)),
    m_problemType(Parameters.get("Problem Type BDDC", SCALARPDE)),
    m_subNumber(Parameters.get("subdomain number", 0)),
    m_interiorSolver(0),
    m_interiorSolverOp(0),
    m_Parameters(Parameters)
  {
    m_work1.resize(m_numDofs);
    m_work2.resize(m_numDofs);
    m_workNC1.resize(m_numDofs);
    m_workNC2.resize(m_numDofs);
    extractMatricesAndFactor();
    m_subVec1.resize(m_numDofs);
    m_subVec2.resize(m_numDofs);
    setDofCoordsAndNodes(nodes);
  }

  SubdomainBDDC(LO numNodes,
		const LO* nodes,
		const LO* nodeBegin,
		const LO* localDofs,
		const LO* rowBegin,
		const LO* columns,
		const SX* values,
		const SM* xCoord,
		const SM* yCoord,
		const SM* zCoord,
		Teuchos::ParameterList & Parameters) :
  m_numDofs(nodeBegin[numNodes]),
    m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_rowBegin(rowBegin),
    m_columns(columns),
    m_values(values),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_Parameters(Parameters)
  {
    // Just factor the original matrix
    setDofCoordsAndNodes(nodes);
    determineReorderRows(nodes);
    std::vector<LO> rows(m_numDofs);
    for (LO i=0; i<m_numDofs; i++) rows[i] = i;
    const std::string solverName = Parameters.get("Coarse Solver", "SuperLU");
    std::string matrixName("Original Matrix");
    factorMatrix(rows, rowBegin, columns, values, solverName, m_subSolver,
		 matrixName);
  }

  ~SubdomainBDDC()
  {
    delete m_interiorSolver;
    delete m_NeumannSolver;
    delete m_interiorSolverOp;
    delete m_subSolver;
  }

  void Solve(int numRhs, 
	     SX* rhs, 
	     SX* sol)
  {
    for (int j=0; j<numRhs; j++) {
      SX* rhsCol = &rhs[j*m_numDofs];
      SX* solCol = &sol[j*m_numDofs];
      if (m_reorderIsIdentity) {
	m_subSolver->Solve(1, rhsCol, solCol);
      }
      else {
	for (LO i=0; i<m_numDofs; i++) {
	  LO row = m_reorderRows[i];
	  m_work1[i] = rhsCol[row];
	}
	m_subSolver->Solve(1, m_work1.data(), m_work2.data());
	for (LO i=0; i<m_numDofs; i++) {
	  LO row = m_reorderRows[i];
	  solCol[row] = m_work2[i];
	}
      }
    }
  }

  LO getNumEquiv() const
  {
    return m_equivClassDofs.size();
  }

  LO getNumDofs() const
  {
    return m_numDofs;
  }

  const std::vector<LO> & getEquivConstraintsBegin() const
  {
    return m_equivConstraintsBegin;
  }

  LO getNumBoundaryDofs() const
  {
    return m_numBoundaryDofs;
  }

  const std::vector<LO> & getInteriorDofs() const
  {
    return m_interiorDofs;
  }

  const LO* getBoundaryDofs() const
  {
    return m_boundaryDofs;
  }

  LO getCoarseSpaceDimension() const
  {
    LO numCoarse = m_pivotRows.size() + m_numAuxConstraints;
    return numCoarse;
  }

  const std::vector<LO> & getEquivLocalDofs(const LO equiv) const
  {
    return m_equivConstraintsLocalDofs[equiv];
  }

  SM getBoundaryDiagValue(LO rowB)
  {
    LO row = m_boundaryDofs[rowB];
    return diagValue(row);
  }

  void getSubdomainVectors(SX* & vec1,
			   SX* & vec2)
  {
    vec1 = m_subVec1.data();
    vec2 = m_subVec2.data();
  }

  void multiplyByPhi(const SX* x, 
		     SX* b, 
		     const bool interfaceValuesOnly,
		     const bool transpose)
  {
    LO numBaseConstraints = m_pivotRows.size();
    LO numCols = numBaseConstraints + m_numAuxConstraints;
    LO numRows(0);
    SX *Phi(0);
    if (interfaceValuesOnly) {
      numRows = m_numBoundaryDofs;
      Phi = m_PhiB.data();
    }
    else {
      numRows = m_numDofs;
      Phi = m_Phi.data();
    }
    Teuchos::BLAS<int, SX> BLAS;
    SX ALPHA(1), BETA(0);
    int INCX(1), INCY(1);
    if ((numRows > 0) && (numCols > 0)) {
      if (transpose == false) {
	BLAS.GEMV(Teuchos::NO_TRANS, numRows, numCols, ALPHA, Phi, numRows, 
		  x, INCX, BETA, b, INCY);
      }
      else {
	BLAS.GEMV(Teuchos::CONJ_TRANS, numRows, numCols, ALPHA, Phi, numRows, 
		  x, INCX, BETA, b, INCY);
      }
    }
  }

  LO getNumInitialEquivConstraints(LO equiv) const
  {
    BDDC_TEST_FOR_EXCEPTION(m_equivConstraintsBegin.size() <= 0, 
		std::runtime_error, "m_equivConstraintsBegin size error");
    LO numConstraints = m_equivConstraintsBegin[equiv+1] -
      m_equivConstraintsBegin[equiv];
    return numConstraints;
  }

  LO getNumInitialFreeDofs(LO equiv) const
  {
    LO numDofsEquiv = getNumEquivDofs(equiv);
    LO numInitialConstraintsEquiv = getNumInitialEquivConstraints(equiv);
    BDDC_TEST_FOR_EXCEPTION(numDofsEquiv < numInitialConstraintsEquiv, 
			    std::runtime_error, "invalid numDofsEquiv");
    return numDofsEquiv - numInitialConstraintsEquiv;
  }

  const std::vector<SX> & getInitialEquivConstraints(LO equiv) const
  {
    return m_equivConstraints[equiv];
  }

  const std::vector<LO> getEquivPivotRows(LO equiv) const
  {
    return m_pivotRowsLocal[equiv];

  }

  void printSparseMatrix(const char* fileName) {
    UtilBDDC<SX,SM>::printSparseMatrix(m_numDofs, m_rowBegin, m_columns,
				       m_values, fileName);
  }

  void outputData(LO equiv,
		  std::vector<int> & boundingECs)
  {
    UtilBDDC<SX,SM>::printSparseMatrix(m_numDofs, m_rowBegin, m_columns, 
				       m_values, "A.dat");
    UtilBDDC<SX,SM>::printIndices(m_equivClassDofs[equiv].size(), 
				  &m_equivClassDofs[equiv][0], "equivDofs.dat");
    UtilBDDC<SX,SM>::printIndices(m_pivotRows.size(), &m_pivotRows[0],
				  "pivotRows.dat");
    UtilBDDC<SX,SM>::printIndices(m_remainRows.size(), &m_remainRows[0],
				  "remainRows.dat");
    // bounding constraints
    LO numConstraints = 0;
    for (size_t i=0; i<boundingECs.size(); i++) {
      LO equiv = boundingECs[i];
      numConstraints += getNumInitialEquivConstraints(equiv);
    }
    std::vector<SX> Cb(numConstraints*m_numDofs);
    LO rowCount = 0;
    for (size_t i=0; i<boundingECs.size(); i++) {
      LO equiv = boundingECs[i];
      const std::vector<LO> & equivDofs = m_equivClassDofs[equiv];
      LO count = 0;
      for (LO j=m_equivConstraintsBegin[equiv]; 
	   j<m_equivConstraintsBegin[equiv+1]; j++) {
	for (size_t k=0; k<equivDofs.size(); k++) {
	  LO dof = equivDofs[k];
	  Cb[rowCount+dof*numConstraints] = m_equivConstraints[equiv][count++];
	}
	rowCount++;
      }
    }
    UtilBDDC<SX,SM>::printDenseMatrix(numConstraints, m_numDofs, &Cb[0], 
				      "Cbound.dat");
    // all base constraints
    LO numEquiv = m_equivClassDofs.size();
    numConstraints = m_equivConstraintsBegin[numEquiv];
    std::vector<SX> C(numConstraints*m_numDofs);
    for (LO i=0; i<numEquiv; i++) {
      const std::vector<LO> & equivDofs = m_equivClassDofs[i];
      LO count = 0;
      for (LO j=m_equivConstraintsBegin[i]; j<m_equivConstraintsBegin[i+1]; j++) {
	for (size_t k=0; k<equivDofs.size(); k++) {
	  LO dof = equivDofs[k];
	  C[j+dof*numConstraints] = m_equivConstraints[i][count++];
	}
      }
    }
    UtilBDDC<SX,SM>::printDenseMatrix(numConstraints, m_numDofs, &C[0], "Cb.dat");
  }

  void getBoundaryDiagValues(LO equiv,
			     std::vector<SM> & diagValues)
  {
    LO numEquivDofs = getNumEquivDofs(equiv);
    diagValues.resize(numEquivDofs);
    for (LO i=0; i<numEquivDofs; i++) {
      LO row = m_equivClassDofs[equiv][i];
      diagValues[i] = diagValue(row);
    }
  }

  size_t getNumEquivDofs(LO equiv) const
  {
    return m_equivClassDofs[equiv].size();
  }

  void setEquivalenceClasses(std::vector< std::vector<LO> > & equivClassDofs)
  {
    m_equivClassDofs = equivClassDofs;
    // check that equivalence classes are disjoint and cover boundary
    size_t numEquiv = equivClassDofs.size();
    LO numBoundaryDofs = 0;
    std::vector<bool> flagDof(m_numDofs, false);
    for (size_t i=0; i<numEquiv; i++) {
      for (size_t j=0; j<equivClassDofs[i].size(); j++) {
	LO dof = equivClassDofs[i][j];
	BDDC_TEST_FOR_EXCEPTION(flagDof[dof] == true, std::runtime_error, 
				"flagDof[dof] error");
	flagDof[dof] = true;
	numBoundaryDofs++;
      }
    }
    BDDC_TEST_FOR_EXCEPTION(numBoundaryDofs != m_numBoundaryDofs, 
			    std::runtime_error, "numBoundaryDofs error");
  }

  void setEquivalenceClassWeightMatrices
    (std::vector< std::vector<LO> > & rowBeginWeight,
     std::vector< std::vector<LO> > & columnsWeight,
     std::vector< std::vector<SX> > & valuesWeight)
  {
    m_rowBeginWeight = rowBeginWeight;
    m_columnsWeight = columnsWeight;
    m_valuesWeight = valuesWeight;
  }

  void applyNeumannCorrection(const SX* g,
			      SX* gSol,
			      const bool restrictToBoundary)
  {
    SX* gScaled = m_workNC1.data();
    bool applyTranspose = true;
    applyWeights(g, applyTranspose, gScaled, restrictToBoundary);
    m_workNC2.assign(m_numDofs, 0);
    SX* gFullScaled = m_workNC2.data();
    if (restrictToBoundary) {
      for (LO i=0; i<m_numBoundaryDofs; i++) {
	gFullScaled[m_boundaryDofs[i]] = gScaled[i];
      }
    }
    else {
      memcpy(gFullScaled, gScaled, m_numDofs*sizeof(SX));
    }
    LO numRowsa = m_numAuxConstraints;
    LO numRowsb = m_pivotRows.size();
    std::vector<SX> ea(numRowsa), solLambdaa(numRowsa), eb(numRowsb),
      solLambdab(numRowsb);
    SX* solFull = m_workNC1.data();
    solveNeumann(gFullScaled, eb.data(), ea.data(), solFull, 
		 solLambdab.data(), solLambdaa.data());
    applyTranspose = false;
    SX* sol = solFull;
    if (restrictToBoundary) {
      sol = m_workNC2.data();
      for (LO i=0; i<m_numBoundaryDofs; i++) {
	sol[i] = solFull[m_boundaryDofs[i]];
      }
    }
    applyWeights(sol, applyTranspose, gSol, restrictToBoundary);
  }

  void setInitialConstraints
    (std::vector< std::vector<LO> > & equivConstraintsLocalDofs,
     std::vector< std::vector<SX> > & equivConstraints)
  {
    m_equivConstraintsLocalDofs = equivConstraintsLocalDofs,
    m_equivConstraints = equivConstraints;
    size_t numEquiv = equivConstraintsLocalDofs.size();
    std::vector<LO> numEquivConstraints(numEquiv);
    for (size_t i=0; i<numEquiv; i++) {
      numEquivConstraints[i] = equivConstraintsLocalDofs[i].size();
    }
    m_pivotRowsLocal.resize(numEquiv);
    m_equivConstraintsBegin.resize(numEquiv+1, 0);
    for (size_t i=0; i<numEquiv; i++) {
      m_equivConstraintsBegin[i+1] = m_equivConstraintsBegin[i] +
	equivConstraintsLocalDofs[i].size();
    }
    BDDC_TEST_FOR_EXCEPTION(numEquiv != m_equivClassDofs.size(), 
			    std::runtime_error, "numEquiv is invalid");
    m_pivotRows.resize(0);
    for (size_t i=0; i<numEquiv; i++) {
      size_t numTerms = numEquivConstraints[i]*m_equivClassDofs[i].size();
      BDDC_TEST_FOR_EXCEPTION(equivConstraints[i].size() != numTerms, 
	       std::runtime_error, "equivConstraints size error");
      findPivotRows(numEquivConstraints[i], equivConstraints[i], 
		    m_equivClassDofs[i], m_pivotRowsLocal[i], m_pivotRows);
    }
    std::vector<bool> flagRows(m_numDofs, false);
    for (size_t i=0; i<m_pivotRows.size(); i++) {
      flagRows[m_pivotRows[i]] = true;
    }
    m_remainRows.resize(0);
    for (LO i=0; i<m_numDofs; i++) {
      if (flagRows[i] == false) m_remainRows.push_back(i);
    }
    std::vector<LO> mapRemain(m_numDofs, -1);
    for (size_t i=0; i<m_remainRows.size(); i++) {
      mapRemain[m_remainRows[i]] = i;
    }
    std::vector<LO> rowBeginA, columnsA;
    std::vector<SX> valuesA;
    const LO *rowBegin, *columns;
    const SX *values;
    getSparseMatrix(rowBegin, columns, values);
    extractMatrix(rowBegin, columns, values, m_remainRows.size(),
		  &m_remainRows[0], &mapRemain[0], rowBeginA, columnsA, valuesA);
    std::string solverName = m_Parameters.get("Neumann Solver", "SuperLU");
    std::string matrixName("Neumann");
    factorMatrix(m_remainRows, &rowBeginA[0], &columnsA[0], 
		 &valuesA[0], solverName, m_NeumannSolver, matrixName);
    setupSaddlePointProblem(numEquivConstraints, equivConstraints);
  }

  void getSparseMatrix(const LO* & rowBegin, 
		       const LO* & columns, 
		       const SX* & values)
  {
    if (m_useBlockPreconditioner) {
      rowBegin = &m_rowBeginBlock[0];
      columns = &m_columnsBlock[0];
      values = &m_valuesBlock[0];
    }
    else {
      rowBegin = m_rowBegin;
      columns = m_columns;
      values = m_values;
    }
  }

  void setAuxiliaryConstraints
    (std::vector< std::vector<LO> > & auxConstraintsLocalDofs,
     std::vector< std::vector<SX> > & auxConstraints)
  {
    size_t numEquiv = auxConstraintsLocalDofs.size();
    if (numEquiv == 0) return;
    LO numRowsa = 0;
    std::vector<LO> numAuxConstraints(numEquiv);
    for (size_t i=0; i<numEquiv; i++) {
      m_equivConstraintsLocalDofs.push_back(auxConstraintsLocalDofs[i]);
      numAuxConstraints[i] = auxConstraintsLocalDofs[i].size();
      numRowsa += numAuxConstraints[i];
    }
    m_numAuxConstraints = numRowsa;
    LO numRows = m_numDofs;
    std::vector<LO> imap(numRows);
    for (LO i=0; i<numRows; i++) imap[i] = i;
    std::vector<SX> CaT;
    extractConstraintMatrix(numAuxConstraints, auxConstraints, numRows,
			    numRowsa, imap, CaT);
    m_AhatInvCaT_X.resize(numRows*numRowsa);
    LO numRowsb = m_pivotRows.size();
    m_AhatInvCaT_Lambda.resize(numRowsb*numRowsa);
    std::vector<SX> eb(numRowsb, 0);
    for (LO i=0; i<numRowsa; i++) {
      SX* solX = &m_AhatInvCaT_X[numRows*i];
      SX* solLambdab = &m_AhatInvCaT_Lambda[numRowsb*i];
      solveOriginalNeumann(&CaT[i*numRows], &eb[0], solX, solLambdab);
    }
    if (numRowsa > 0) {
      Teuchos::BLAS<int, SX>  BLAS;
      SX ALPHA(-1), BETA(0);
      m_Sp2.resize(numRowsa*numRowsa);
      if (numRows > 0) {
	BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRowsa, numRowsa, 
		  numRows, ALPHA, &CaT[0], numRows, &m_AhatInvCaT_X[0], numRows, 
		  BETA, &m_Sp2[0], numRowsa);
      }
      m_Sp2IPIV.resize(numRowsa);
      Teuchos::LAPACK<int, SX> LAPACK;
      int INFO(0);
      int N = numRowsa;
      LAPACK.GETRF(N, N, &m_Sp2[0], N, &m_Sp2IPIV[0], &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRF error");
    }
  }

  void solveOriginalNeumann(const SX* g,
			    const SX* eb,
			    SX* solX,
			    SX* solLambdab)
  {
    SX* g1 = &m_work1[0];
    SX* g2 = &m_work3[0];
    for (size_t i=0; i<m_remainRows.size(); i++) {
      g1[i] = g[m_remainRows[i]];
    }
    LO numRows1 = m_remainRows.size();
    LO numRows2 = m_pivotRows.size();
    for (LO i=0; i<numRows2; i++) {
      g2[i]          = g[m_pivotRows[i]];
      g2[i+numRows2] = eb[i];
    }
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(-1), BETA(1);
    int N = 2*numRows2;
    SM maxAbsG1 = getMaxAbsValue(g1, m_remainRows.size());
    if ((maxAbsG1 > 0) && (N > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, 1, 
		numRows1, ALPHA, &m_A11invA12[0], numRows1, g1, numRows1, 
		BETA, g2, N);
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, 1, 
		numRows1, ALPHA, &m_A11invCb1T[0], numRows1, g1, numRows1, 
		BETA, &g2[numRows2], N);
    }
    Teuchos::LAPACK<int, SX> LAPACK;
    if (N > 0) {
      int INFO(0);
      LAPACK.GETRS('N', N, 1, &m_Sp1[0], N, &m_Sp1IPIV[0], 
		   g2, N, &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
    }
    SX* A11invg1 = &m_work2[0];
    m_NeumannSolver->Solve(1, g1, A11invg1);
    SM maxAbsG2_first = getMaxAbsValue(g2, numRows2);
    if ((maxAbsG2_first > 0) && (N > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, 1,
		numRows2, ALPHA, &m_A11invA12[0], numRows1, g2, N,
		BETA, A11invg1, numRows1);
    }
    SM maxAbsG2_second = getMaxAbsValue(&g2[numRows2], numRows2);
    if ((maxAbsG2_second > 0) && (N > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, 1,
		numRows2, ALPHA, &m_A11invCb1T[0], numRows1, &g2[numRows2], N,
		BETA, A11invg1, numRows1);
    }
    for (size_t i=0; i<m_remainRows.size(); i++) {
      solX[m_remainRows[i]] = A11invg1[i];
    }
    for (size_t i=0; i<m_pivotRows.size(); i++) {
      solX[m_pivotRows[i]] = g2[i];
      solLambdab[i] = g2[numRows2+i];
    }
  }

  void solveNeumann(const SX* g,
		    const SX* eb,
		    const SX* ea,
		    SX* solX,
		    SX* solLambdab,
		    SX* solLambdaa)
  {
    solveOriginalNeumann(g, eb, solX, solLambdab);
    if (m_numAuxConstraints == 0) return;
    LO numRows = m_numDofs;
    LO numRowsb = m_pivotRows.size();
    LO numRowsa = m_numAuxConstraints;    
    for (LO i=0; i<numRowsa; i++) solLambdaa[i] = ea[i];
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(-1), BETA(1);
    SM maxAbsG = getMaxAbsValue(g, numRows);
    if ((maxAbsG > 0) && (numRows > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRowsa, 1, 
		numRows, ALPHA, &m_AhatInvCaT_X[0], numRows, g, numRows, 
		BETA, solLambdaa, numRowsa);
    }
    SM maxAbsEb = getMaxAbsValue(eb, numRowsb);
    if ((maxAbsEb > 0) && (numRowsb > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRowsa, 1, 
		numRowsb, ALPHA, &m_AhatInvCaT_Lambda[0], numRowsb, eb, 
		numRowsb, BETA, solLambdaa, numRowsa);
    }
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO;
    LAPACK.GETRS('N', numRowsa, 1, &m_Sp2[0], numRowsa, 
		 &m_Sp2IPIV[0], solLambdaa, numRowsa, &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
    if (numRows > 0) {
      BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows, 1,
		numRowsa, ALPHA, &m_AhatInvCaT_X[0], numRows,
		solLambdaa, numRowsa, BETA, solX, numRows);
    }
    if (numRowsb > 0) {
      BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRowsb, 1,
		numRowsa, ALPHA, &m_AhatInvCaT_Lambda[0], numRowsb, 
		solLambdaa, numRowsa, BETA, solLambdab, numRowsb);
    }
  }

  void getBoundingMatrices(std::vector<LO> & constraintRows, 
			   std::vector<SX> & Sp1, 
			   std::vector<SX> & A11invCb1T)
  {
    LO numRows2 = m_pivotRows.size();
    LO numRowsb = constraintRows.size();
    LO dimSp1 = numRows2 + numRowsb;
    Sp1.resize(dimSp1*dimSp1);
    LO dimSp1full = 2*numRows2;
    BDDC_TEST_FOR_EXCEPTION(dimSp1full < dimSp1, std::runtime_error, 
			    "invalid dimSp1full");
    std::vector<LO> indexMap(dimSp1full, -1);
    for (LO i=0; i<numRows2; i++) indexMap[i] = i;
    for (LO i=0; i<numRowsb; i++) {
      LO row = constraintRows[i];
      indexMap[numRows2+row] = numRows2 + i;
    }
    extractSquareMatrix(&m_Sp1NotFactored[0], dimSp1full, &Sp1[0], dimSp1, 
			&indexMap[0]);
    LO numRows1 = m_remainRows.size();
    A11invCb1T.resize(numRows1*numRowsb);
    indexMap.assign(numRows2, -1);
    for (LO i=0; i<numRowsb; i++) indexMap[constraintRows[i]] = i;
    for (LO j=0; j<numRows2; j++) {
      LO col = indexMap[j];
      if (col != -1) {
	for (LO i=0; i<numRows1; i++) {
	  A11invCb1T[i+col*numRows1] = m_A11invCb1T[i+j*numRows1];
	}
      }
    }
  }

  void calculateSchurComplement(LO equivClass,
				std::vector<SX> & Sc)
  {
    if (m_interfacePreconditioner == true) {
      if (m_equivToBoundaryMap.size() == 0) initializeEquivBoundaryMaps();
      std::vector<SX> xB(m_numBoundaryDofs), SxB(m_numBoundaryDofs);
      LO numEquivDofs = m_equivClassDofs[equivClass].size();
      Sc.resize(numEquivDofs*numEquivDofs);
      for (LO i=0; i<numEquivDofs; i++) {
	LO row = m_equivToBoundaryMap[equivClass][i];
	xB[row] = 1;
	SX ALPHA(1), BETA(0);
	applyBoundaryOperator(&xB[0], &SxB[0], ALPHA, BETA, PRECONDITIONER);
	xB[row] = 0;
	for (LO j=0; j<numEquivDofs; j++) {
	  row = m_equivToBoundaryMap[equivClass][j];
	  Sc[j+i*numEquivDofs] = SxB[row];
	}
      }
    }
  }

  void getDenseWeightMatrix(LO equivClass,
			    std::vector<SX> & Weight)
  {
    const std::vector<LO> & rowBegin = m_rowBeginWeight[equivClass];
    const std::vector<LO> & columns = m_columnsWeight[equivClass];
    const std::vector<SX> & values = m_valuesWeight[equivClass];
    LO numRows = getNumEquivDofs(equivClass);
    BDDC_TEST_FOR_EXCEPTION(numRows != LO(rowBegin.size()-1), 
			    std::runtime_error, "invalid numRows");
    Weight.resize(numRows*numRows, 0);
    for (LO i=0; i<numRows; i++) {
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	LO col = columns[j];
	Weight[i+col*numRows] = values[j];
      }
    }
  }

  void determineBoundingECs(LO equivClass,
			    const std::vector< std::vector<LO> > & equivSubs,
			    std::vector<LO> & boundingECs)
  {
    boundingECs.resize(0);
    LO numEquiv = m_equivClassDofs.size();
    const std::vector<LO> subs = equivSubs[equivClass];
    LO numEquivDofs = getNumEquivDofs(equivClass);
    LO numEquivConstraints = getNumInitialEquivConstraints(equivClass);
    BDDC_TEST_FOR_EXCEPTION(numEquivConstraints > numEquivDofs, 
			    std::runtime_error, "invalid numEquivConstraints");
    if (numEquivConstraints == numEquivDofs) return;
    for (LO i=0; i<numEquiv; i++) {
      if (m_equivConstraintsBegin[i+1] > m_equivConstraintsBegin[i]) {
	const std::vector<LO> subs2 = equivSubs[i];
	std::vector<int> common;
	std::set_intersection(subs.begin(), subs.end(),
			      subs2.begin(), subs2.end(),
			      std::back_inserter(common));
	if (common.size() == subs.size()) {
	  boundingECs.push_back(i);
	}
      }
    }
  }

  void calculateBoundingSchurComplement(LO equivClass,
					std::vector<LO> & boundingECs,
					std::vector<SX> & Sc,
					LO & dimSc,
					LO & numRowsb)
  {
    std::vector<LO> constraintRows;
    for (size_t i=0; i<boundingECs.size(); i++) {
      LO ec = boundingECs[i];
      for (LO j=m_equivConstraintsBegin[ec]; j<m_equivConstraintsBegin[ec+1]; 
	   j++) {
	constraintRows.push_back(j);
      }
    }
    std::vector<SX> Sp1, A11invCb1T;
    getBoundingMatrices(constraintRows, Sp1, A11invCb1T);
    LO numRowsEquiv = m_equivClassDofs[equivClass].size();
    LO numRows1 = m_remainRows.size();
    LO numRows2 = m_pivotRows.size();
    numRowsb = constraintRows.size();
    LO dimSp1 = numRows2 + numRowsb;
    dimSc = numRowsEquiv + numRowsb;
    std::vector<SX> Rhs(dimSp1*dimSc), g(m_numDofs, 0);
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(-1), BETA(1);
    SX* g1 = &m_work1[0];
    for (LO k=0; k<dimSc; k++) {
      SX* g2 = &Rhs[k*dimSp1];
      loadG1G2(k, g, numRowsEquiv, m_equivClassDofs[equivClass], numRows1,
	       numRows2, g1, g2);
      SM maxAbsG1 = getMaxAbsValue(g1, numRows1);
      if ((maxAbsG1 > 0) && (dimSp1 > 0) && (numRows1 > 0)) {
	BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, 1, 
		  numRows1, ALPHA, &m_A11invA12[0], numRows1, g1, numRows1, 
		  BETA, g2, dimSp1);
	BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRowsb, 1, 
		  numRows1, ALPHA, &A11invCb1T[0], numRows1, g1, numRows1, 
		  BETA, &g2[numRows2], dimSp1);
      }
    }
    Teuchos::LAPACK<int, SX> LAPACK;
    std::vector<SM> S(dimSp1), RWORK(std::max(1, 5*dimSp1));
    std::vector<SX> WORK(1);
    SM RCOND(1e-10);
    int RANK, LWORK(-1), INFO(0);
    if (dimSp1 > 0) {
      LAPACK.GELSS(dimSp1, dimSp1, dimSc, &Sp1[0], dimSp1, &Rhs[0], dimSp1, 
		   &S[0], RCOND, &RANK, &WORK[0], LWORK, &RWORK[0], &INFO);
      LWORK = int(UtilBDDC<SX,SM>::real(WORK[0])+0.1);
      WORK.resize(LWORK);
      LAPACK.GELSS(dimSp1, dimSp1, dimSc, &Sp1[0], dimSp1, &Rhs[0], dimSp1, 
		   &S[0], RCOND, &RANK, &WORK[0], LWORK, &RWORK[0], &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GELSS error");
    }
    m_work3.assign(dimSp1, 0);
    std::vector<SX> Xe(numRowsEquiv*numRowsEquiv), Xb(numRowsEquiv*numRowsb),
      Lambdae(numRowsb*numRowsEquiv), Lambdab(numRowsb*numRowsb);
    for (int k=0; k<dimSc; k++) {
      SX* g2 = &m_work3[0];
      loadG1G2(k, g, numRowsEquiv, m_equivClassDofs[equivClass], numRows1,
	       numRows2, g1, g2);
      SX* A11invg1 = &m_work2[0];
      m_NeumannSolver->Solve(1, g1, A11invg1);
      SX* x2 = &Rhs[k*dimSp1];
      if ((dimSp1 > 0) && (numRows1 > 0)) {
	BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, 1,
		  numRows2, ALPHA, &m_A11invA12[0], numRows1, x2, dimSp1,
		  BETA, A11invg1, numRows1);      
	BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, 1,
		  numRowsb, ALPHA, &A11invCb1T[0], numRows1, &x2[numRows2], dimSp1,
		  BETA, A11invg1, numRows1);
      }
      SX *Xtarget(0), *LambdaTarget(0);
      if (k < numRowsEquiv) {
	Xtarget = &Xe[k*numRowsEquiv];
	LambdaTarget = &Lambdae[k*numRowsb];
      }
      else {
	Xtarget = &Xb[(k-numRowsEquiv)*numRowsEquiv];
	LambdaTarget = &Lambdab[(k-numRowsEquiv)*numRowsb];
      }
      SX* xAll = &m_work1[0];
      for (LO i=0; i<numRows1; i++) xAll[m_remainRows[i]] = A11invg1[i];
      for (LO i=0; i<numRows2; i++) xAll[m_pivotRows[i]] = x2[i];
      for (LO i=0; i<numRowsEquiv; i++) {
	LO row = m_equivClassDofs[equivClass][i];
	Xtarget[i] = xAll[row];
      }
      for (LO i=0; i<numRowsb; i++) {
	LambdaTarget[i] = x2[numRows2+i];
      }
    }
    std::vector<SX> equivRhs(numRowsEquiv*dimSc, 0);
    for (LO i=0; i<numRowsEquiv; i++) {
      equivRhs[i+numRowsEquiv*i] = 1;
    }
    for (LO i=0; i<numRowsEquiv; i++) {
      for (LO j=0; j<numRowsb; j++) {
	LO col = numRowsEquiv + j;
	equivRhs[i+numRowsEquiv*col] = Xb[i+j*numRowsEquiv];
      }
    }
    std::vector<int> IPIV(numRowsEquiv);
    LAPACK.GESV(numRowsEquiv, dimSc, &Xe[0], numRowsEquiv, &IPIV[0], 
		&equivRhs[0], numRowsEquiv, &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GESV error");
    if (numRowsb > 0) {
      BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRowsb, numRowsb,
		numRowsEquiv, ALPHA, &Lambdae[0], numRowsb, 
		&equivRhs[numRowsEquiv*numRowsEquiv], numRowsEquiv,
		BETA, &Lambdab[0], numRowsb);
    }
    Sc.resize(dimSc*dimSc);
    // bounding dofs first then equivalence class dofs next
    for (LO i=0; i<numRowsEquiv; i++) {
      LO row = numRowsb + i;
      for (LO j=0; j<numRowsEquiv; j++) {
	LO col = numRowsb + j;
	Sc[row+col*dimSc] = equivRhs[i+j*numRowsEquiv];
      }
      for (LO j=0; j<numRowsb; j++) {
	LO col = j + numRowsEquiv;
	SX value = -equivRhs[i+col*numRowsEquiv];
	Sc[row+j*dimSc] = value;
	Sc[j+row*dimSc] = UtilBDDC<SX,SM>::conj(value);
      }
    }
    for (LO i=0; i<numRowsb; i++) {
      for (LO j=0; j<numRowsb; j++) {
	Sc[i+j*dimSc] = -Lambdab[i+j*numRowsb];
      }
    }
  }

  void calculateCoarseMatrices(std::vector<SX> & Ac,
			       bool restrictPhiToBoundary = false,
			       bool scaleRows = false)
  {
    LO numRows = m_numDofs;
    LO numBaseConstraints = m_pivotRows.size();
    LO numCols = numBaseConstraints + m_numAuxConstraints;
    m_Phi.resize(numRows*numCols);
    std::vector<SX> eb(numBaseConstraints, 0), solLambdab(numBaseConstraints),
      ea(m_numAuxConstraints), solLambdaa(m_numAuxConstraints);
    m_work1.assign(numRows, 0);
    SX* g = &m_work1[0];
    for (LO i=0; i<numBaseConstraints; i++) {
      eb[i] = 1;
      solveNeumann(g, &eb[0], &ea[0], &m_Phi[i*numRows], 
		   &solLambdab[0], &solLambdaa[0]);
      eb[i] = 0;
    }
    for (LO i=0; i<m_numAuxConstraints; i++) {
      ea[i] = 1;
      solveNeumann(g, &eb[0], &ea[0], &m_Phi[(i+numBaseConstraints)*numRows],
		   &solLambdab[0], &solLambdaa[0]);
      ea[i] = 0;
    }
    std::vector<SX> APhi(numRows*numCols);
    for (LO i=0; i<numCols; i++) {
      applyFullOperator(&m_Phi[i*numRows], &APhi[i*numRows]);
    }
    Ac.resize(numCols*numCols);
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(1), BETA(0);
    if ((numRows > 0) && (numCols > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numCols, numCols, 
		numRows, ALPHA, &m_Phi[0], numRows, &APhi[0], numRows, 
		BETA, &Ac[0], numCols);
    }
    if (scaleRows == true) {
      bool applyTranspose = false;
      std::vector<SX> Phi(m_Phi);
      for (LO i=0; i<numCols; i++) {
	applyWeights(&Phi[i*numRows], applyTranspose, &m_Phi[i*numRows], false);
      }
    }
    if (restrictPhiToBoundary == true) {
      LO numRowsB = m_numBoundaryDofs;
      m_PhiB.resize(numRowsB*numCols);
      for (LO j=0; j<numCols; j++) {
	for (LO i=0; i<numRowsB; i++) {
	  m_PhiB[i+j*numRowsB] = m_Phi[m_boundaryDofs[i]+j*numRows];
	}
      }
      /*
      m_Phi.resize(0);
      m_Phi.shrink_to_fit();
      */
    }
  }

  void getPhi(std::vector<SX> & Phi,
	      bool restrictedToBoundary = false) const
  {
    if (restrictedToBoundary == true) {
      Phi = m_PhiB;
    }
    else {
      Phi = m_Phi;
    }
  }

  void applyFullOperator(SX* x,
			 SX* Ax)
  {
    for (LO i=0; i<m_numDofs; i++) {
      SX sum = 0;
      for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	LO col = m_columns[j];
	sum += m_values[j]*x[col];
      }
      Ax[i] = sum;
    }
  }

  void staticExpansion(SX* xB,
		       SX* solI,
		       const bool useCurrentInteriorValues)
  {
    const LO* rowBeginIB = &m_rowBeginIB[0];
    const LO* columnsIB = &m_columnsIB[0];
    const SX* valuesIB = &m_valuesIB[0];
    SolverBase<SX>* interiorSolver = m_interiorSolver;
    if (m_useBlockPreconditioner == true) {
      rowBeginIB = &m_rowBeginIBOp[0];
      columnsIB = &m_columnsIBOp[0];
      valuesIB = &m_valuesIBOp[0];
      interiorSolver = m_interiorSolverOp;
    }
    LO numInteriorDofs = m_interiorDofs.size();
    SX* rhsSolver = &m_work1[0];
    for (LO i=0; i<numInteriorDofs; i++) {
      rhsSolver[i] = 0;
      for (LO j=rowBeginIB[i]; j<rowBeginIB[i+1]; j++) {
	LO col = columnsIB[j];
	rhsSolver[i] -= valuesIB[j]*xB[col];
      }
      if (useCurrentInteriorValues) {
	for (LO j=m_rowBeginII[i]; j<m_rowBeginII[i+1]; j++) {
	  const LO col = m_columnsII[j];
	  rhsSolver[i] -= m_valuesII[j]*solI[col];
	}
      }
    }
    interiorSolver->Solve(1, rhsSolver, solI);
  }

  void applyBoundaryOperator(SX* x,
			     SX* Sx,
			     SX ALPHA=1,
			     SX BETA=0,
			     enum MatrixType matrixType=OPERATOR)
  {
    const LO* rowBeginBB = &m_rowBeginBB[0];
    const LO* rowBeginIB = &m_rowBeginIB[0];
    const LO* rowBeginBI = &m_rowBeginBI[0];
    const LO* columnsBB = &m_columnsBB[0];
    const LO* columnsIB = &m_columnsIB[0];
    const LO* columnsBI = &m_columnsBI[0];
    const SX* valuesBB = &m_valuesBB[0];
    const SX* valuesIB = &m_valuesIB[0];
    const SX* valuesBI = &m_valuesBI[0];
    SolverBase<SX>* interiorSolver = m_interiorSolver;
    if ((matrixType == OPERATOR) && (m_useBlockPreconditioner == true)) {
      rowBeginBB = &m_rowBeginBBOp[0];
      rowBeginIB = &m_rowBeginIBOp[0];
      rowBeginBI = &m_rowBeginBIOp[0];
      columnsBB = &m_columnsBBOp[0];
      columnsIB = &m_columnsIBOp[0];
      columnsBI = &m_columnsBIOp[0];
      valuesBB = &m_valuesBBOp[0];
      valuesIB = &m_valuesIBOp[0];
      valuesBI = &m_valuesBIOp[0];
      interiorSolver = m_interiorSolverOp;
    }
    applyBoundaryOperator(x, Sx, ALPHA, BETA, 
			  rowBeginBB, columnsBB, valuesBB,
			  rowBeginIB, columnsIB, valuesIB,
			  rowBeginBI, columnsBI, valuesBI, 
			  interiorSolver);
  }

  void applyBoundaryOperator(SX* x,
			     SX* Sx,
			     SX ALPHA,
			     SX BETA,
			     const LO* rowBeginBB,
			     const LO* columnsBB,
			     const SX* valuesBB,
			     const LO* rowBeginIB,
			     const LO* columnsIB,
			     const SX* valuesIB,
			     const LO* rowBeginBI,
			     const LO* columnsBI,
			     const SX* valuesBI,
			     SolverBase<SX>* & interiorSolver)
  {
    LO numInteriorDofs = m_interiorDofs.size();
    SX* rhsSolver = &m_work1[0];
    SX* solSolver = &m_work2[0];
    for (LO i=0; i<numInteriorDofs; i++) {
      rhsSolver[i] = 0;
      for (LO j=rowBeginIB[i]; j<rowBeginIB[i+1]; j++) {
	LO col = columnsIB[j];
	rhsSolver[i] -= valuesIB[j]*x[col];
      }
    }
    interiorSolver->Solve(1, rhsSolver, solSolver);
    for (LO i=0; i<m_numBoundaryDofs; i++) {
      Sx[i] *= BETA;
      for (LO j=rowBeginBB[i]; j<rowBeginBB[i+1]; j++) {
	LO col = columnsBB[j];
	Sx[i] += ALPHA*valuesBB[j]*x[col];
      }
      for (LO j=rowBeginBI[i]; j<rowBeginBI[i+1]; j++) {
	LO col = columnsBI[j];
	Sx[i] += ALPHA*valuesBI[j]*solSolver[col];
      }
    }
  }

  void makeInteriorAdjustments(SX* rhs,
			       SX* sol)
  {
    makeInteriorAdjustments(rhs, sol, !m_interfacePreconditioner);
  }

  void makeInteriorAdjustments(SX* rhs,
			       SX* sol,
			       const bool notInterfacePreconditioner)
  {
    // On input, rhs has load vector for interior unknowns
    // On output, sol has solution vector for interior unknowns
    // On output, the first m_numBoundaryDofs entries of rhs contain associated
    //  adjustments to boundary forces. If m_interfacePreconditioner is false,
    //  then the next rhs entries are adjustments to interior residuals
    SolverBase<SX>* interiorSolver = m_interiorSolver;
    interiorSolver->Solve(1, rhs, sol);
    // boundary adjustments to rhs
    for (LO i=0; i<m_numBoundaryDofs; i++) {
      SX sum = 0;
      for (LO j=m_rowBeginBI[i]; j<m_rowBeginBI[i+1]; j++) {
	const LO col = m_columnsBI[j];
	sum += m_valuesBI[j]*sol[col];
      }
      rhs[i] = sum;
    }
    // interior adjustments to rhs
    if (notInterfacePreconditioner) {
      const LO numInteriorDofs = m_interiorDofs.size();
      for (LO i=0; i<numInteriorDofs; i++) {
	SX sum = 0;
	for (LO j=m_rowBeginII[i]; j<m_rowBeginII[i+1]; j++) {
	  const LO col = m_columnsII[j];
	  sum += m_valuesII[j]*sol[col];
	}
	rhs[i+m_numBoundaryDofs] = sum;      
      }
    }
  }
  
  void applyWeights(const SX* g,
		    const bool applyTranspose,
		    SX* gScaled,
		    bool restrictToBoundary)
  {
    BDDC_TEST_FOR_EXCEPTION(m_rowBeginWeight.size() != m_equivClassDofs.size(),
			    std::runtime_error, "m_rowBeginWeight size error");
    LO numEquiv = m_rowBeginWeight.size();
    if (m_equivToBoundaryMap.size() == 0) initializeEquivBoundaryMaps();
    if (restrictToBoundary == true) {
      for (LO i=0; i<m_numBoundaryDofs; i++) gScaled[i] = 0;
      for (LO i=0; i<numEquiv; i++) {
	for (size_t j=0; j<m_equivClassDofs[i].size(); j++) {
	  LO rowB = m_equivToBoundaryMap[i][j];
	  for (LO k=m_rowBeginWeight[i][j]; k<m_rowBeginWeight[i][j+1]; k++) {
	    LO colEquiv = m_columnsWeight[i][k];
	    LO colB = m_equivToBoundaryMap[i][colEquiv];
	    if (applyTranspose == false) {
	      gScaled[rowB] += m_valuesWeight[i][k]*g[colB];
	    }
	    else {
	      gScaled[colB] += m_valuesWeight[i][k]*g[rowB];
	    }
	  }
	}
      }
    }
    else {
      for (LO i=0; i<m_numDofs; i++) gScaled[i] = g[i];
      for (LO i=0; i<m_numBoundaryDofs; i++) {
	LO rowBB = m_boundaryDofs[i];
	gScaled[rowBB] = 0;
      }
      for (LO i=0; i<numEquiv; i++) {
	for (size_t j=0; j<m_equivClassDofs[i].size(); j++) {
	  LO rowB = m_equivToBoundaryMap[i][j];
	  LO rowBB = m_boundaryDofs[rowB];
	  for (LO k=m_rowBeginWeight[i][j]; k<m_rowBeginWeight[i][j+1]; k++) {
	    LO colEquiv = m_columnsWeight[i][k];
	    LO colB = m_equivToBoundaryMap[i][colEquiv];
	    LO colBB = m_boundaryDofs[colB];
	    if (applyTranspose == false) {
	      gScaled[rowBB] += m_valuesWeight[i][k]*g[colBB];
	    }
	    else {
	      gScaled[colBB] += m_valuesWeight[i][k]*g[rowBB];
	    }
	  }
	}
      }
    }
  }

  void factorInteriorMatrix()
  {
    std::string solverName = m_Parameters.get("Dirichlet Solver", "SuperLU");
    std::string matrixName("Dirichlet");
    factorMatrix(m_interiorDofs, m_rowBeginII.data(), m_columnsII.data(), 
		 m_valuesII.data(), solverName, m_interiorSolver, matrixName);
    // free memory as appropriate
    /*
    if (m_interfacePreconditioner == true) {
      m_rowBeginII.resize(0);
      m_rowBeginII.shrink_to_fit();
      m_columnsII.resize(0);
      m_columnsII.shrink_to_fit();
      m_valuesII.resize(0);
      m_valuesII.shrink_to_fit();
    }
    */
  }

 private: // member data
  LO m_numDofs, m_numNodes;
  const LO *m_nodeBegin, *m_localDofs, *m_rowBegin, *m_columns;
  const SX *m_values;
  const SM *m_xCoord{nullptr}, *m_yCoord{nullptr}, *m_zCoord{nullptr};
  LO m_numBoundaryDofs, m_numAuxConstraints;
  const LO *m_boundaryDofs;
  bool m_interfacePreconditioner, m_useBlockPreconditioner;
  enum ProblemType m_problemType;
  int m_subNumber{0};
  std::vector<LO> m_interiorDofs;
  SolverBase<SX> *m_interiorSolver{nullptr}, *m_NeumannSolver{nullptr};
  SolverBase<SX> *m_interiorSolverOp{nullptr}, *m_subSolver{nullptr};
  std::vector<SX> m_work1, m_work2, m_work3, m_workNC1, m_workNC2, 
    m_Phi, m_PhiB, m_subVec1, m_subVec2;
  std::vector<LO> m_rowBeginIB, m_columnsIB, m_rowBeginBI, m_columnsBI,
    m_rowBeginBB, m_columnsBB, m_rowBeginII, m_columnsII, 
    m_pivotRows, m_remainRows, m_equivConstraintsBegin;
  std::vector< std::vector<LO> > m_equivConstraintsLocalDofs;
  std::vector<SX> m_valuesIB, m_valuesBI, m_valuesBB, m_valuesII;
  std::vector< std::vector<LO> >  m_equivClassDofs, m_pivotRowsLocal;
  std::vector< std::vector<LO> > m_rowBeginWeight, m_columnsWeight,
    m_equivToBoundaryMap;
  std::vector< std::vector<SX> > m_valuesWeight, m_equivConstraints;
  std::vector<SX> m_A11invCb1T, m_A11invA12;
  std::vector<SX> m_Sp1, m_Sp2, m_AhatInvCaT_X, m_AhatInvCaT_Lambda,
    m_Sp1NotFactored;
  std::vector<int> m_Sp1IPIV, m_Sp2IPIV;
  std::vector<LO> m_rowBeginBBOp, m_rowBeginIBOp, m_rowBeginBIOp,
    m_columnsBBOp, m_columnsIBOp, m_columnsBIOp, m_rowBeginIIOp,
    m_columnsIIOp;
  std::vector<SX> m_valuesBBOp, m_valuesIBOp, m_valuesBIOp, m_valuesIIOp;
  std::vector<LO> m_rowBeginBlock, m_columnsBlock, m_reorderRows;
  std::vector<SX> m_valuesBlock;
  Teuchos::ParameterList & m_Parameters;
  std::vector<SM> m_dofCoords;
  std::vector<LO> m_dofNodes;
  std::vector<double> m_xCoordAMG, m_yCoordAMG, m_zCoordAMG;
  bool m_reorderIsIdentity{false};

 private: // methods

  void determineReorderRows(const LO* nodes)
  {
    m_reorderIsIdentity = true;
    for (int i=1; i<m_numNodes; i++) {
      if (nodes[i-1] > nodes[i]) {
	m_reorderIsIdentity = false;
	break;
      }
    }
    if (m_reorderIsIdentity) return;
    m_reorderRows.resize(m_numDofs);
    LO numRows = 0;
    for (LO i=0; i<m_numNodes; i++) {
      LO node = nodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	m_reorderRows[numRows++] = j;
      }
    }
    m_work1.resize(m_numDofs);
    m_work2.resize(m_numDofs);
  }

  void setDofCoordsAndNodes(const LO* nodes)
  {
    m_dofCoords.resize(3*m_numDofs);
    m_dofNodes.resize(m_numDofs);
    for (LO i=0; i<m_numNodes; i++) {
      const LO node = nodes[i];
      for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	m_dofCoords[3*j+0] = m_xCoord[node];
	m_dofCoords[3*j+1] = m_yCoord[node];
	m_dofCoords[3*j+2] = m_zCoord[node];
	m_dofNodes[j] = i;
      }
    }
  }

  SM diagValue(LO row)
  {
    SM value(0);
    for (LO i=m_rowBegin[row]; i<m_rowBegin[row+1]; i++) {
      LO col = m_columns[i];
      if (col == row) value = std::abs(m_values[i]);
    }
    return value;
  }

  void loadG1G2(const int k,
		std::vector<SX> & g,
		LO numRowsEquiv, 
		std::vector<LO> & equivClassDofs, 
		LO numRows1,
		LO numRows2,
		SX* g1,
		SX* g2)
  {
    if (k < numRowsEquiv) {
      LO row = equivClassDofs[k];
      g[row] = 1;
      for (LO i=0; i<numRows1; i++) {
	g1[i] = g[m_remainRows[i]];
      }
      for (LO i=0; i<numRows2; i++) {
	g2[i] = g[m_pivotRows[i]];
      }
      g[row] = 0;
    }
    else {
      for (LO i=0; i<numRows1; i++) {
	g1[i] = 0;
      }
      LO row = (k - numRowsEquiv) + numRows2;
      g2[row] = 1;
    }
  }

  void extractSquareMatrix(const SX* Afull,
			   LO dimAfull,
			   SX* A, 
			   LO dimA,
			   const LO* indexMap)
  {
    for (LO i=0; i<dimAfull; i++) {
      LO row = indexMap[i];
      if (row != -1) {
	for (LO j=0; j<dimAfull; j++) {
	  LO col = indexMap[j];
	  if (col != -1) {
	    A[row+dimA*col] = Afull[i+j*dimAfull];
	  }
	}
      }
    }
  }

  SM getMaxAbsValue(const SX* values,
		    LO numValues)
  {
    SM maxAbsValue(0);
    for (LO i=0; i<numValues; i++) {
      SM absValue = std::abs(values[i]);
      if (absValue > maxAbsValue) maxAbsValue = absValue;
    }
    return maxAbsValue;
  }

  void setupSaddlePointProblem
    (std::vector<LO> & numEquivConstraints, 
     std::vector< std::vector<SX> > & equivConstraints)
  {
    LO numRows1 = m_remainRows.size();
    LO numRows2 = m_pivotRows.size(); // number of "base" constraints
    std::vector<LO> map1(m_numDofs, -1), map2(m_numDofs, -1);
    for (LO i=0; i<numRows1; i++) map1[m_remainRows[i]] = i;
    for (LO i=0; i<numRows2; i++) map2[m_pivotRows[i]] = i;
    std::vector<SX> A12, Cb1T;
    const LO *rowBegin, *columns;
    const SX *values;
    getSparseMatrix(rowBegin, columns, values);
    extractDenseMatrix(rowBegin, columns, values, m_remainRows, 
		       m_pivotRows, map2, A12);
    extractConstraintMatrix(numEquivConstraints, equivConstraints, numRows1,
			    numRows2, map1, Cb1T);
    m_A11invCb1T.resize(numRows1*numRows2);
    m_A11invA12.resize(numRows1*numRows2);
    for (LO i=0; i<numRows2; i++) {
      LO index = i*numRows1;
      m_NeumannSolver->Solve(1, &Cb1T[index], &m_A11invCb1T[index]);
      m_NeumannSolver->Solve(1, &A12[index],  &m_A11invA12[index]);
    }
    std::vector<SX> A22;
    extractDenseMatrix(rowBegin, columns, values, m_pivotRows,
		       m_pivotRows, map2, A22);
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(-1), BETA(1);
    if ((numRows2 > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, numRows2, 
		numRows1, ALPHA, &A12[0], numRows1, &m_A11invA12[0], numRows1, 
		BETA, &A22[0], numRows2);
    }
    std::vector<SX> Cb2T;
    extractConstraintMatrix(numEquivConstraints, equivConstraints, numRows2,
			    numRows2, map2, Cb2T);
    if ((numRows2 > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, numRows2, 
		numRows1, ALPHA, &A12[0], numRows1, &m_A11invCb1T[0], numRows1, 
		BETA, &Cb2T[0], numRows2);
    }
    std::vector<SX> Cb1A11invCb1T(numRows2*numRows2);
    ALPHA = 1; BETA = 0;
    if ((numRows2 > 0) && (numRows1 > 0)) {
      BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numRows2, numRows2, 
		numRows1, ALPHA, &Cb1T[0], numRows1, &m_A11invCb1T[0], numRows1, 
		BETA, &Cb1A11invCb1T[0], numRows2);
    }
    int N = 2*numRows2;
    m_work3.resize(N);
    m_Sp1.resize(N*N, 0);
    SX* A = &m_Sp1[0];
    for (int i=0; i<numRows2; i++) {
      for (int j=0; j<numRows2; j++) {
	A[i+N*(j+       0)] = A22[i+j*numRows2];
	A[i+N*(j+numRows2)] = Cb2T[i+j*numRows2];
	A[i+numRows2+N*(j+ 0      )] = 
	  UtilBDDC<SX,SM>::conj(Cb2T[j+i*numRows2]);
	A[i+numRows2+N*(j+numRows2)] = -Cb1A11invCb1T[i+j*numRows2];
      }
    }
    m_Sp1NotFactored = m_Sp1;
    m_Sp1IPIV.resize(N);
    int* IPIV = &m_Sp1IPIV[0];
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    if (N > 0) {
      LAPACK.GETRF(N, N, A, N, IPIV, &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRF error");
    }
  }

  void extractConstraintMatrix
      (std::vector<LO> & numEquivConstraints, 
       std::vector< std::vector<SX> > & equivConstraints, 
       LO numRows1,
       LO numRows2, 
       std::vector<LO> & rowMap,
       std::vector<SX> & C)
  {
    C.resize(numRows1*numRows2);
    size_t numEquiv = numEquivConstraints.size();
    LO count1 = 0;
    for (size_t i=0; i<numEquiv; i++) {
      LO count2 = 0;
      for (LO j=0; j<numEquivConstraints[i]; j++) {
	for (size_t k=0; k<m_equivClassDofs[i].size(); k++) {
	  LO row = rowMap[m_equivClassDofs[i][k]];
	  if (row != -1) {
	    C[row+count1*numRows1] = equivConstraints[i][count2];
	  }
	  count2++;
	}
	count1++;
      }
    }
    BDDC_TEST_FOR_EXCEPTION(count1 != numRows2, std::runtime_error, 
			    "count1 error");
  }

  void extractDenseMatrix(const LO* rowBeginA,
			  const LO* columnsA,
			  const SX* valuesA,
			  std::vector<LO> & rows1, 
			  std::vector<LO> & rows2,
			  std::vector<LO> & map2,
			  std::vector<SX> & A12)
  {
    LO numRows1 = rows1.size();
    LO numRows2 = rows2.size();
    A12.resize(numRows1*numRows2);
    for (LO i=0; i<numRows1; i++) {
      LO row = rows1[i];
      for (LO j=rowBeginA[row]; j<rowBeginA[row+1]; j++) {
	LO col = map2[columnsA[j]];
	if (col != -1) A12[i+col*numRows1] = valuesA[j];
      }
    }
  }

  void getNodalData(const std::vector<LO> & rows, 
		    std::vector<int> & nodeBegin,
		    std::vector<int> & localDofs,
		    std::vector<double> & xCoords, 
		    std::vector<double> & yCoords, 
		    std::vector<double> & zCoords)
  {
    const int numRows = rows.size();
    localDofs.resize(numRows);
    std::vector<int> nodeMap(m_numNodes, -1), count(m_numNodes, 0);
    // first, initialize localDofs and determine active nodes
    int numActiveNodes(0);
    for (int i=0; i<numRows; i++) {
      const int row = rows[i];
      localDofs[i] = m_localDofs[row];
      const int node = m_dofNodes[row];
      if (nodeMap[node] == -1) {
	nodeMap[node] = numActiveNodes;
	count[numActiveNodes]++;
	numActiveNodes++;
      }
      else {
	count[nodeMap[node]]++;
      }
    }
    nodeBegin.resize(numActiveNodes+1, 0);
    for (int i=0; i<numActiveNodes; i++) {
      nodeBegin[i+1] = nodeBegin[i] + count[i];
      count[i] = 0;
    }
    xCoords.resize(numActiveNodes);
    yCoords.resize(numActiveNodes);
    zCoords.resize(numActiveNodes);
    numActiveNodes = 0;
    for (int i=0; i<numRows; i++) {
      const int row = rows[i];
      const int node = nodeMap[m_dofNodes[row]];
      if (count[node] == 0) {
	xCoords[numActiveNodes] = m_dofCoords[3*row+0];
	yCoords[numActiveNodes] = m_dofCoords[3*row+1];
	zCoords[numActiveNodes] = m_dofCoords[3*row+2];
	numActiveNodes++;
      }
      count[node]++;
    }
  }

  void factorMatrix(const std::vector<LO> & rows, 
		    const LO* rowBegin, 
		    const LO* columns, 
		    const SX* values,
		    const std::string & solverName,
		    SolverBase<SX>* & solver,
		    const std::string & matrixName)
  {
    LO numRows = rows.size();
    if (matrixName == "Dohrmann") {
      std::string fileName = matrixName + std::to_string(m_subNumber) + ".dat";
      UtilBDDC<SX,SM>::printSparseMatrix(numRows, rowBegin, columns,
					 values, fileName.c_str());
    }
    std::vector<int> nodeBegin, localDofs;
    if ((solverName == "MueLu") || (solverName == "NodalAMG")) {
      getNodalData(rows, nodeBegin, localDofs, 
		   m_xCoordAMG, m_yCoordAMG, m_zCoordAMG);
      int numNodes = nodeBegin.size() - 1;
      m_Parameters.set("numNodes", numNodes);
      m_Parameters.set("nodeBegin", nodeBegin.data());
      m_Parameters.set("localDofs", localDofs.data());
      if (solverName == "NodalAMG") {
	m_Parameters.set("Need to Copy Matrix", 1);
      }
      m_Parameters.set("xCoords", m_xCoordAMG.data());
      m_Parameters.set("yCoords", m_yCoordAMG.data());
      m_Parameters.set("zCoords", m_zCoordAMG.data());
    }
    SolverFactory<SX> Factory;
    m_Parameters.set("Solver", solverName);
    solver = Factory.Generate(numRows,
			      const_cast<int*>(rowBegin),
			      const_cast<int*>(columns),
			      const_cast<SX*>(values),
			      m_Parameters);
    solver->Initialize();
  }
  
  void determineLocalPivotRows(LO numConstraints,
			       LO numEquivDofs, 
			       const std::vector<SX> & equivConstraints,
			       std::vector<LO> & localPivotRows) const
  {
    std::vector<SX> normalizedConstraints(numEquivDofs*numConstraints, 0);
    for (LO i=0; i<numConstraints; i++) {
      SM maxAbsValue(0);
      for (LO j=0; j<numEquivDofs; j++) {
	LO index = j + i*numEquivDofs;
	SM absValue = std::abs(equivConstraints[index]);
	if (absValue > maxAbsValue) maxAbsValue = absValue;
      }
      BDDC_TEST_FOR_EXCEPTION(maxAbsValue <= 0, std::runtime_error, 
			      "maxAbsValue must be greater than zero");
      for (LO j=0; j<numEquivDofs; j++) {
	LO index = j + i*numEquivDofs;
	normalizedConstraints[index] = 
	  equivConstraints[index]/maxAbsValue;
      }
    }
    // partial pivoting to identify pivot rows (constraints assumed
    // to be linearly independent)
    for (LO i=0; i<numConstraints; i++) {
      SM maxValue(0);
      LO pivot = -1;
      for (LO j=0; j<numEquivDofs; j++) {
	LO index = j + i*numEquivDofs;
	SM absValue = std::abs(normalizedConstraints[index]);
	if (absValue > maxValue) {
	  maxValue = absValue;
	  pivot = j;
	}
      }
      BDDC_TEST_FOR_EXCEPTION(pivot == -1, std::runtime_error, "invalid pivot");
      for (LO k=i+1; k<numConstraints; k++) {
	SX num = normalizedConstraints[pivot + k*numEquivDofs];
	SX den = normalizedConstraints[pivot + i*numEquivDofs];
	SX scaleFactor = num/den;
	for (LO j=0; j<numEquivDofs; j++) {
	  normalizedConstraints[j+k*numEquivDofs] -=
	    scaleFactor*normalizedConstraints[j+i*numEquivDofs];
	}
      }
      localPivotRows.push_back(pivot);
    }
  }

  void findPivotRows(LO numConstraints, 
		     std::vector<SX> & equivConstraints,
		     std::vector<LO> & equivClassDofs,
		     std::vector<LO> & localPivotRows,
		     std::vector<LO> & pivotRows) const
  {
    LO numEquivDofs = equivClassDofs.size();
    determineLocalPivotRows(numConstraints, numEquivDofs, equivConstraints,
			    localPivotRows);
    for (LO i=0; i<numConstraints; i++) {
      LO pivot = localPivotRows[i];
      pivotRows.push_back(equivClassDofs[pivot]);
    }
  }

  void extractMatrix(const LO* rowBegin, 
		     const LO* columns, 
		     const SX* values, 
		     const LO numRowsA,
		     const LO* rowsA, 
		     const LO* columnMapA,
		     std::vector<LO> & rowBeginA, 
		     std::vector<LO> & columnsA, 
		     std::vector<SX> & valuesA)
  {
    rowBeginA.resize(numRowsA+1, 0);
    LO numTerms(0);
    for (LO i=0; i<numRowsA; i++) {
      LO row = rowsA[i];
      for (LO j=rowBegin[row]; j<rowBegin[row+1]; j++) {
	LO col = columnMapA[columns[j]];
	if (col != -1) numTerms++;
      }
      rowBeginA[i+1] = numTerms;
    }
    columnsA.resize(numTerms);
    valuesA.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numRowsA; i++) {
      LO row = rowsA[i];
      for (LO j=rowBegin[row]; j<rowBegin[row+1]; j++) {
	LO col = columnMapA[columns[j]];
	if (col != -1) {
	  columnsA[numTerms] = col;
	  valuesA[numTerms] = values[j];
	  numTerms++;
	}
      }
    }
  }

  void checkForAppropriateBlocks()
  {
    LO numFlux(0); // localDofs > 6 are considered flux dofs for DPG methods
    for (LO i=0; i<m_numDofs; i++) {
      if (m_localDofs[i] > 6) {
	numFlux++;
      }
    }
    BDDC_TEST_FOR_EXCEPTION(numFlux <= 0, std::runtime_error, 
			    "numFlux must be positive");
  }
  
  void getBlockMatrix(std::vector<LO> & rowBeginBlock, 
		      std::vector<LO> & columnsBlock, 
		      std::vector<SX> & valuesBlock)
  {
    LO numTerms = 0;
    for (LO i=0; i<m_numDofs; i++) {
      LO row = i;
      for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	LO col = m_columns[j];
 	if (sameBlock(m_localDofs[row], m_localDofs[col])) {
	  numTerms++;
	}
      }
    }
    rowBeginBlock.resize(m_numDofs+1, 0);
    columnsBlock.resize(numTerms);
    valuesBlock.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<m_numDofs; i++) {
      LO row = i;
      for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	LO col = m_columns[j];
 	if (sameBlock(m_localDofs[row], m_localDofs[col])) {
	  columnsBlock[numTerms] = col;
	  valuesBlock[numTerms] = m_values[j];
	  numTerms++;
	}
      }
      rowBeginBlock[i+1] = numTerms;
    }
  }

  bool sameBlock(LO localDof1, 
		 LO localDof2)
  {
    bool isSameBlock = false;
    if ((localDof1 < 7) && (localDof2 < 7)) isSameBlock = true;
    if ((localDof1 > 6) && (localDof2 > 6)) isSameBlock = true;
    return isSameBlock;
  }

  void extractMatricesAndFactor()
  {
    BDDC_TEST_FOR_EXCEPTION(m_numDofs < m_numBoundaryDofs, std::runtime_error, 
			    "m_numDofs too small");
    LO numInteriorDofs = m_numDofs - m_numBoundaryDofs;
    m_interiorDofs.resize(numInteriorDofs);
    std::vector<LO> mapB(m_numDofs, -1), mapI(m_numDofs, -1);
    for (LO i=0; i<m_numBoundaryDofs; i++) {
      mapB[m_boundaryDofs[i]] = i;
    }
    numInteriorDofs = 0;
    for (LO i=0; i<m_numDofs; i++) {
      if (mapB[i] == -1) {
	mapI[i] = numInteriorDofs;
	m_interiorDofs[numInteriorDofs] = i;
	numInteriorDofs++;
      }
    }
    if (m_useBlockPreconditioner == true) {
      checkForAppropriateBlocks();
      getBlockMatrix(m_rowBeginBlock, m_columnsBlock, m_valuesBlock);
      const LO* rowBegin = &m_rowBeginBlock[0];
      const LO* columns = &m_columnsBlock[0];
      const SX* values = &m_valuesBlock[0];
      getBlockMatrices(rowBegin, columns, values, mapI, mapB,
		       m_rowBeginIB, m_columnsIB, m_valuesIB,
		       m_rowBeginBB, m_columnsBB, m_valuesBB,
		       m_rowBeginBI, m_columnsBI, m_valuesBI,
		       m_rowBeginII, m_columnsII, m_valuesII);
      getBlockMatrices(m_rowBegin, m_columns, m_values, mapI, mapB,
		       m_rowBeginIBOp, m_columnsIBOp, m_valuesIBOp,
		       m_rowBeginBBOp, m_columnsBBOp, m_valuesBBOp,
		       m_rowBeginBIOp, m_columnsBIOp, m_valuesBIOp,
		       m_rowBeginIIOp, m_columnsIIOp, m_valuesIIOp);
    }
    else {
      getBlockMatrices(m_rowBegin, m_columns, m_values, mapI, mapB,
		       m_rowBeginIB, m_columnsIB, m_valuesIB,
		       m_rowBeginBB, m_columnsBB, m_valuesBB,
		       m_rowBeginBI, m_columnsBI, m_valuesBI,
		       m_rowBeginII, m_columnsII, m_valuesII);
    }
  }
  
  void getBlockMatrices(const LO* rowBegin, 
			const LO* columns, 
			const SX* values,
			std::vector<LO> & mapI, 
			std::vector<LO> & mapB,
			std::vector<LO> & rowBeginIB, 
			std::vector<LO> & columnsIB, 
			std::vector<SX> & valuesIB,
			std::vector<LO> & rowBeginBB, 
			std::vector<LO> & columnsBB, 
			std::vector<SX> & valuesBB,
			std::vector<LO> & rowBeginBI, 
			std::vector<LO> & columnsBI, 
			std::vector<SX> & valuesBI,
			std::vector<LO> & rowBeginII, 
			std::vector<LO> & columnsII, 
			std::vector<SX> & valuesII)

  {
    LO numInteriorDofs = m_interiorDofs.size();
    // extract AIB
    extractMatrix(rowBegin, columns, values, numInteriorDofs,
		  &m_interiorDofs[0], &mapB[0], rowBeginIB, columnsIB, valuesIB);
    // extract ABB
    extractMatrix(rowBegin, columns, values, m_numBoundaryDofs,
		  m_boundaryDofs, &mapB[0], rowBeginBB, columnsBB, valuesBB);
    // extract ABI
    extractMatrix(rowBegin, columns, values, m_numBoundaryDofs,
		  m_boundaryDofs, &mapI[0], rowBeginBI, columnsBI, valuesBI);
    // extract AII
    extractMatrix(rowBegin, columns, values, numInteriorDofs,
		  &m_interiorDofs[0], &mapI[0], rowBeginII, columnsII, valuesII);
  }

  void initializeEquivBoundaryMaps()
  {
    LO numEquiv = m_equivClassDofs.size();
    m_equivToBoundaryMap.resize(numEquiv);
    std::vector<LO> allToBoundary(m_numDofs, -1);
    for (LO i=0; i<m_numBoundaryDofs; i++) {
      allToBoundary[m_boundaryDofs[i]] = i;
    }
    for (LO i=0; i<numEquiv; i++) {
      LO numEquivDofs = m_equivClassDofs[i].size();
      m_equivToBoundaryMap[i].resize(numEquivDofs);
      for (LO j=0; j<numEquivDofs; j++) {
	LO dof = m_equivClassDofs[i][j];
	m_equivToBoundaryMap[i][j] = allToBoundary[dof];
	BDDC_TEST_FOR_EXCEPTION(allToBoundary[dof] == -1, std::runtime_error, 
				"invalid allToBoundary[dof]");
      }
    }
  }

};

} // namespace bddc

#endif // BDDC_SUBDOMAIN_H
  
