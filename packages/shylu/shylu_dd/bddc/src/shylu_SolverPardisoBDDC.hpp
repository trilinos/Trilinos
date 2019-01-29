
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

#ifndef BDDC_SOLVERPARDISO_H
#define BDDC_SOLVERPARDISO_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "mkl.h"
#include "mkl_pardiso.h"
#include "shylu_SolverBaseBDDC.hpp"
#include "shylu_UtilPardiso.hpp"

namespace bddc {

template <class SX> 
  class SolverPardiso : 
  public SolverBase<SX>
{
public:
  ~SolverPardiso()
  {
    delete [] m_rowBeginP;
    delete [] m_columnsP;
    delete [] m_valuesP;
    delete [] m_perm;
  }
  SolverPardiso(int numRows,
		int* rowBegin,
		int* columns,
		SX* values,
		Teuchos::ParameterList & Parameters) :
  SolverBase<SX>(numRows, rowBegin, columns, values, Parameters),
    m_matrixIsSymmetric(true),
    m_rowBeginP(0),
    m_columnsP(0),
    m_perm(0),
    m_matrixType(2),
    m_valuesP(0),
    m_messageLevel(0)
    {
    }
  
  int Initialize()
  {
    setMatrixType();
    int noMessage = 0;
    m_messageLevel = this->m_Parameters.get("Pardiso Message Level", 
					    noMessage);
    int error = InitializePardiso(this->m_numRows, 
				  this->m_rowBegin, 
				  this->m_columns, 
				  this->m_values);
    return error;
  }

  bool IsDirectSolver()
  {
    return true;
  }

  void MySolve(int NRHS,
	       SX* Rhs, 
	       SX* Sol)
  {
    int n = this->m_numRows;
    int one(1), error;
    int phase = 33; // solve phase
    pardiso((_MKL_DSS_HANDLE_t*)m_pt, &one, &one, &m_matrixType, &phase, 
	    &n, m_valuesP, m_rowBeginP, m_columnsP, m_perm, &NRHS, m_iparam, 
	    &m_messageLevel, Rhs, Sol, &error);
    BDDC_TEST_FOR_EXCEPTION(error != 0, std::runtime_error, 
			    "Pardiso solve error");
  }

  bool MyExactSolver() {
    return true;
  }

private:
  void setMatrixType()
  {
    // Note: only consider real matrices at this time
    std::string matrixName = 
      this->m_Parameters.get("Pardiso Matrix Type",
			     "Real And Symmetric Positive Definite");
    if (matrixName == "Real And Structurally Symmetric") {
      m_matrixType = 1;
      m_matrixIsSymmetric = false;
    }
    else if (matrixName == "Real And Symmetric Positive Definite") {
      m_matrixType = 2;
      m_matrixIsSymmetric = true;
    }
    else if (matrixName == "Real And Symmetric Indefinite") {
      m_matrixType = -2;
      m_matrixIsSymmetric = true;
    }
    else if (matrixName == "Real And Nonsymmetric") {
      m_matrixType = 11;
      m_matrixIsSymmetric = false;
    }
  }

  int InitializePardiso(int numRows,
			const int* rowBegin,
			const int* columns,
			const SX* values)
  {
    // get matrix in Pardiso format
    UtilPardiso<int,SX>::constructPardisoMatrix
      (numRows, rowBegin, columns, values, m_matrixIsSymmetric,
       m_rowBeginP, m_columnsP, m_valuesP);
    // Pardiso initialization
    int n = numRows;
    m_perm = new int[n];
    for (int i=0; i<64; i++) {
      m_pt[i] = 0;       // intialized to zero, don't change later
      m_iparam[i] = 0;   // to use default pardiso parameters
    }
    m_iparam[0] = 0; // use default Pardiso parameters
    int phase = 12;    // analysis and factorization phase
    int one(1), NRHS(1), error(0);
    SX RHS[1], SOL[1];
    double startWallTime = this->wall_time();
    pardiso((_MKL_DSS_HANDLE_t*)m_pt, &one, &one, &m_matrixType, &phase, &n, 
	    m_valuesP, m_rowBeginP, m_columnsP, m_perm, &NRHS, m_iparam, 
	    &m_messageLevel, RHS, SOL, &error);
    double endWallTime = this->wall_time();
    double deltaT = endWallTime - startWallTime;
    if (error != 0) std::cout << "> Pardiso error = " << error << std::endl;
    /*
      std::cout << "factoring problem matrix took " << deltaT << " seconds\n";
      std::cout << "> pardiso nnz " << iparam[17] << "\n";
    */
    return error;
  }

  bool m_matrixIsSymmetric;
  int *m_rowBeginP, *m_columnsP, *m_perm, m_matrixType;
  SX *m_valuesP;
  int m_messageLevel;
  long m_pt[64];
  MKL_INT m_iparam[64];
 };
  
} // namespace bddc

#endif // BDDC_SOLVERPARDISO_H
  
