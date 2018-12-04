
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

#ifndef BDDC_SOLVERMUELU_H
#define BDDC_SOLVERMUELU_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "shylu_SolverBaseBDDC.hpp"
#include "shylu_enumsBDDC.hpp"
#include "shylu_errorBDDC.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
template <class SX> class SolverMueLu : 
  public SolverBase<SX>
{
public:
  typedef int LO; // Local Ordinal
  typedef Tpetra::Map<>::global_ordinal_type GO; // Global Ordinal
  typedef Tpetra::Map<LO,GO>                                 Map;
  typedef Tpetra::MultiVector<double,LO,GO>                  MV;
  typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
  typedef Tpetra::CrsMatrix<double,LO,GO>                    CrsMatrix;
  typedef Tpetra::Operator<double,LO,GO>                     Operator;
  typedef MueLu::TpetraOperator<double,LO,GO>                MTOperator;

  SolverMueLu(int numRows,
	      int* rowBegin,
	      int* columns,
	      SX* values,
	      Teuchos::ParameterList & Parameters) :
  SolverBase<SX>(numRows, rowBegin, columns, values, Parameters),
    m_useNullSpaceCorrection(Parameters.get("nullSpaceCorrection", false))
  {
  }

  ~SolverMueLu()
  {
  }

  int Initialize()
  {
    RCP<Teuchos::ParameterList> paramListNull;
    RCP<Teuchos::ParameterList> paramList = 
      this->m_Parameters.get("MueLu Parameter List", paramListNull);

    // null space stuff
    const int numNodes = this->m_Parameters.get("numNodes", 0);
    int *nullPtrInt(nullptr);
    int* nodeBegin = this->m_Parameters.get("nodeBegin", nullPtrInt);
    int* localDofs = this->m_Parameters.get("localDofs", nullPtrInt);
    double *nullptrDouble(0);
    const double* xCoords = this->m_Parameters.get("xCoords", nullptrDouble);
    const double* yCoords = this->m_Parameters.get("yCoords", nullptrDouble);
    const double* zCoords = this->m_Parameters.get("zCoords", nullptrDouble);
    const int numRows = this->m_numRows;
    initializeNullSpaceCorrection(numNodes, nodeBegin, localDofs,
				  xCoords, yCoords, zCoords, paramList);
    RCP<const Teuchos::Comm<int> > Comm = 
      rcp( new Teuchos::MpiComm<int>(MPI_COMM_SELF) );
    initializeCoordinates(numNodes, xCoords, yCoords, zCoords, Comm);
    initializeNullspace(numRows, Comm);
    // MueLu preconditioner
    Teuchos::ArrayRCP<size_t> count(numRows, 0); 
    const LO* rowBegin = this->m_rowBegin;
    const LO* columns = this->m_columns;
    const double* values = reinterpret_cast<const double*>(this->m_values);
    for (LO i=0; i<numRows; i++) {
      count[i] = rowBegin[i+1] - rowBegin[i];
    }
    LO numTerms(0);
    RCP<const Map> rowMap = m_inNullspace->getMap();
    RCP<CrsMatrix> A = rcp( new CrsMatrix(rowMap, rowMap, count) );
    for (LO i=0; i<numRows; i++) {
      A->insertLocalValues
	(i, Teuchos::ArrayView<const LO>(&columns[numTerms], count[i]),
	 Teuchos::ArrayView<const double>(&values[numTerms], count[i]));
      numTerms += count[i];
    }
    A->fillComplete();
    /*
    const LO numNonzeros = A->getGlobalNumEntries();
    std::cout << "numNonzeros = " << numNonzeros << std::endl;
    for (LO i=0; i<numRows; i++) {
      std::cout << "row " << i << " numEntries = " 
		<< A->getNumEntriesInLocalRow(i) << std::endl;
    }
    */
    m_mueLuPreconditioner = 
      MueLu::CreateTpetraPreconditioner((RCP<Operator>)A, *paramList, 
					m_inCoords, m_inNullspace);
    m_rhsVec = rcp( new MV(rowMap, 1) );
    m_solVec = rcp( new MV(rowMap, 1) );
    return 0;
  }

  bool IsDirectSolver()
  {
    return false;
  }

  void MySolve(int NRHS,
	       double* Rhs, 
	       double* Sol)
  {
    int numRows = this->m_numRows;
    int numCols = m_IPIV.size();
    if (numRows == 0) return;
    Teuchos::ArrayRCP<double> rhsVals = m_rhsVec->getDataNonConst(0);
    Teuchos::ArrayRCP<double> solVals = m_solVec->getDataNonConst(0);
    // may make sense to optimize this later by using true multivectors
    for (int j=0; j<NRHS; j++) {
      SX* rhsCol = &Rhs[numRows*j];
      // apply null space pre-correction if requested
      if (m_useNullSpaceCorrection) {
	applyNullSpaceCorrection(numRows, numCols, rhsCol, m_sol1.data(), 
				 m_rhs2.data());
	for (int i=0; i<numRows; i++) rhsVals[i] = m_rhs2[i];
      }
      else {
	for (int i=0; i<numRows; i++) rhsVals[i] = rhsCol[i];
      }
      // apply preconditioner
      m_mueLuPreconditioner->apply(*m_rhsVec, *m_solVec);
      // apply null space post-correction if requested
      SX* solCol = &Sol[numRows*j];
      if (m_useNullSpaceCorrection) {
	applyNullSpaceCorrection2(numRows, numCols, solVals.getRawPtr(),
				  m_rhs2.data());
	// sum corrections
	for (int i=0; i<numRows; i++) {
	  solCol[i] = m_sol1[i] + solVals[i];
	}
      }
      else {
	for (int i=0; i<numRows; i++) solCol[i] = solVals[i];
      }
    }
  }

  void MySolve(int NRHS,
	       float* Rhs, 
	       float* Sol)
  {
    std::string msg("Error: SolverMueLu does not support float type");
    throw msg;
  }

  bool MyExactSolver() 
  {
    return false;
  }

private:
  RCP<MTOperator> m_mueLuPreconditioner;
  RCP<MV> m_rhsVec, m_solVec, m_inCoords, m_inNullspace;
  std::vector<SX> m_N, m_AN, m_NTAN, m_ATN, m_work, m_sol1, m_rhs2;
  std::vector<int> m_IPIV;
  const bool m_useNullSpaceCorrection{false};

  void initializeNullspace(const int numRows,
			   RCP<const Teuchos::Comm<int> > Comm)
  {
    const int numCols = m_IPIV.size();
    const int numTerms = m_N.size();
    BDDC_TEST_FOR_EXCEPTION(numTerms != numRows*numCols, std::runtime_error, 
			    "size error in initializeNullSpaceCorrection");
    RCP<const Map> rowMap = rcp( new Map(numRows, 0, Comm) );
    m_inNullspace = rcp( new MV(rowMap, numCols) );
    for (int j=0; j<numCols; j++) {
      const SX* vals = &m_N[numRows*j];
      Teuchos::ArrayRCP<double> valsVec = m_inNullspace->getDataNonConst(j);
      for (int i=0; i<numRows; i++) {
	valsVec[i] = vals[i];
      }
    }
  }

  void initializeCoordinates(const int numNodes, 
			     const double* xCoords, 
			     const double* yCoords, 
			     const double* zCoords,
			     RCP<const Teuchos::Comm<int> > & Comm)
  {
    RCP<const Map> nodeMap = rcp( new Map(numNodes, 0, Comm) );
    m_inCoords = rcp( new MV(nodeMap, 3) );
    Teuchos::ArrayRCP<double> xVals = m_inCoords->getDataNonConst(0);
    Teuchos::ArrayRCP<double> yVals = m_inCoords->getDataNonConst(1);
    Teuchos::ArrayRCP<double> zVals = m_inCoords->getDataNonConst(2);
    for (int i=0; i<numNodes; i++) {
      xVals[i] = xCoords[i];
      yVals[i] = yCoords[i];
      zVals[i] = zCoords[i];
    }
  }

  void applyNullSpaceCorrection(const LO numRows,
				const LO numCols,
				SX* rhsIn, 
				SX* solOut, 
				SX* rhsOut)
  {
    SX* NTrhsIn = m_work.data();
    Teuchos::BLAS<int, SX> BLAS;
    SX ALPHA(1), BETA(0);
    int INCX(1), INCY(1);
    BLAS.GEMV(Teuchos::CONJ_TRANS, numRows, numCols, ALPHA, m_N.data(), 
	      numRows, rhsIn, INCX, BETA, NTrhsIn, INCY);
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    LAPACK.GETRS('N', numCols, 1, m_NTAN.data(), numCols, m_IPIV.data(), 
		 NTrhsIn, numCols, &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
    memcpy(rhsOut, rhsIn, numRows*sizeof(SX));
    ALPHA = 1; BETA = 0;
    BLAS.GEMV(Teuchos::NO_TRANS, numRows, numCols, ALPHA, m_N.data(), 
	      numRows, NTrhsIn, INCX, BETA, solOut, INCY);
    ALPHA = -1; BETA = 1;
    BLAS.GEMV(Teuchos::NO_TRANS, numRows, numCols, ALPHA, m_AN.data(), 
	      numRows, NTrhsIn, INCX, BETA, rhsOut, INCY);
  }

  void applyNullSpaceCorrection2(const LO numRows,
				 const LO numCols,
				 SX* sol,
				 SX* rhs)
  {
    SX* NTrhs = m_work.data();
    Teuchos::BLAS<int, SX> BLAS;
    SX ALPHA(1), BETA(0);
    int INCX(1), INCY(1);
    BLAS.GEMV(Teuchos::CONJ_TRANS, numRows, numCols, ALPHA, m_N.data(), 
	      numRows, rhs, INCX, BETA, NTrhs, INCY);
    ALPHA = -1; BETA = 1;
    BLAS.GEMV(Teuchos::CONJ_TRANS, numRows, numCols, ALPHA, m_ATN.data(), 
	      numRows, sol, INCX, BETA, NTrhs, INCY);
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    LAPACK.GETRS('N', numCols, 1, m_NTAN.data(), numCols, m_IPIV.data(), 
		 NTrhs, numCols, &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
    ALPHA = 1; BETA = 1;
    BLAS.GEMV(Teuchos::NO_TRANS, numRows, numCols, ALPHA, m_N.data(), 
	      numRows, NTrhs, INCX, BETA, sol, INCY);
  }

  void initializeNullSpaceCorrection
    (const int numNodes,
     const int* nodeBegin,
     const int* localDofs,
     const double* xCoords, 
     const double* yCoords, 
     const double* zCoords,
     RCP<Teuchos::ParameterList> & paramList)
  {
    const std::string problemType = 
      paramList->get("problem: type", "Poisson-3D");
    LO numCols(0);
    const int numRows = nodeBegin[numNodes];
    if (problemType == "Elasticity-3D") {
      numCols = 6;
      m_N.resize(numRows*numCols, 0);
      const int numDofPerNode = 3;
      BDDC_TEST_FOR_EXCEPTION(numRows != 3*numNodes, std::runtime_error, 
			      "invalid numRows");
      // loop over nodes
      for (int i=0; i<numNodes; i++) {
	int localDofCount(0);
	for (int j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  if (localDofs[j] == 0) {
	    // x-translation
	    m_N[j + 0*numRows] =  1;
	    m_N[j + 4*numRows] =  zCoords[i];
	    m_N[j + 5*numRows] = -yCoords[i];
	    localDofCount++;
	  }
	  else if (localDofs[j] == 1) {
	    // y-translation
	    m_N[j + 1*numRows] =  1;
	    m_N[j + 5*numRows] =  xCoords[i];
	    m_N[j + 3*numRows] = -zCoords[i];
	    localDofCount++;
	  }
	  else if (localDofs[j] == 2) {
	    // z-translation
	    m_N[j + 2*numRows] =  1;
	    m_N[j + 3*numRows] =  yCoords[i];
	    m_N[j + 4*numRows] = -xCoords[i];
	    localDofCount++;
	  }
	}
	BDDC_TEST_FOR_EXCEPTION(localDofCount != numDofPerNode, 
				std::runtime_error, "invalid localDofCount");
      }
    }
    else if (problemType == "Elasticity-2D") {
      numCols = 3;
      m_N.resize(numRows*numCols, 0);
      const int numDofPerNode = 2;
      BDDC_TEST_FOR_EXCEPTION(numRows != 2*numNodes, std::runtime_error, 
			      "invalid numRows");
      for (int i=0; i<numNodes; i++) {
	int localDofCount(0);
	for (int j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  if (localDofs[j] == 0) {
	    // x-translation
	    m_N[j + 0*numRows] =  1;
	    m_N[j + 2*numRows] = -yCoords[i];
	    localDofCount++;
	  }
	  else if (localDofs[j] == 1) {
	    // y-translation
	    m_N[j + 1*numRows] =  1;
	    m_N[j + 2*numRows] =  xCoords[i];
	    localDofCount++;
	  }
	}
	BDDC_TEST_FOR_EXCEPTION(localDofCount != numDofPerNode, 
				std::runtime_error, "invalid localDofCount");
      }
    }
    else {
      numCols = 1;
      m_N.resize(numRows, 1);
    }
    m_work.resize(numCols);
    m_rhs2.resize(numRows);
    m_sol1.resize(numRows);
    const LO* rowBegin = this->m_rowBegin;
    const LO* columns = this->m_columns;
    const SX* values = this->m_values;
    m_AN.resize(numRows*numCols, 0);
    m_ATN.resize(numCols*numRows, 0);
    for (LO k=0; k<numCols; k++) {
      SX* Nk = &m_N[numRows*k];
      SX* ANk = &m_AN[numRows*k];
      SX* ATNk = &m_ATN[numRows*k];
      for (LO i=0; i<numRows; i++) {
	for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	  const LO col = columns[j];
	  ANk[i]    += values[j]*Nk[col];
	  // Note: values[j] = A(i,col) = AT(col,i)
	  ATNk[col] += values[j]*Nk[i];
	}
      }
    }
    m_NTAN.resize(numCols*numCols);
    Teuchos::BLAS<int, SX> BLAS;
    SX ALPHA(1), BETA(0);
    BLAS.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, numCols, numCols, 
	      numRows, ALPHA, m_N.data(), numRows, m_AN.data(), numRows, 
	      BETA, m_NTAN.data(), numCols);
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    m_IPIV.resize(numCols);
    LAPACK.GETRF(numCols, numCols, m_NTAN.data(), numCols, m_IPIV.data(), 
		 &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRF error");
  }
 };
  
} // namespace bddc

#endif // BDDC_SOLVERMUELU_H
  
