
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

#ifndef SOLVERLAPACKBDDC_H
#define SOLVERLAPACKBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "shylu_SolverBaseBDDC.h"

namespace bddc {

template <class SX> 
  class SolverLAPACK : public SolverBase<SX>
{
public:
  ~SolverLAPACK()
  {
  }
  SolverLAPACK(int numRows,
	       int* rowBegin,
	       int* columns,
	       SX* values,
	       Teuchos::ParameterList & Parameters) :
  SolverBase<SX>(numRows, rowBegin, columns, values, Parameters)
    {
    }
  
  int Initialize()
  {
    int numRows = this->m_numRows;
    const int* rowBegin = this->m_rowBegin;
    const int* columns = this->m_columns;
    const SX* values = this->m_values;
    int INFO;
    m_factor.resize(numRows*numRows, 0); 
    for (int i=0; i<numRows; i++) {
      for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	int col = columns[j];
	m_factor[numRows*col+i] = values[j]; // column-major storage
      }
    }
    m_IPIV.resize(numRows);
    Teuchos::LAPACK<int, SX> LAPACK;
    LAPACK.GETRF(numRows, numRows, m_factor.data(), numRows, 
		 m_IPIV.data(), &INFO);
    assert (INFO == 0);
    return INFO;
  }

  bool IsDirectSolver()
  {
    return true;
  }

  void MySolve(int NRHS,
	       SX* Rhs, 
	       SX* Sol)
  {
    int numRows = this->m_numRows;
    memcpy(Sol, Rhs, numRows*NRHS*sizeof(SX));
    int INFO;
    char NOTRANS = 'N';
    Teuchos::LAPACK<int, SX> LAPACK;
    LAPACK.GETRS(NOTRANS, numRows, NRHS, m_factor.data(), numRows, 
		 m_IPIV.data(), Sol, numRows, &INFO);
    assert (INFO == 0);
  }

  bool MyExactSolver() {
    return true;
  }

private:
  std::vector<SX> m_factor;
  std::vector<int> m_IPIV;
 };
  
} // namespace bddc

#endif // SOLVERLAPACKBDDC_H
  
