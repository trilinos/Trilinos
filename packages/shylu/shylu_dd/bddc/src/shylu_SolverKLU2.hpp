
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

#ifndef BDDC_SOLVERKLU2_H
#define BDDC_SOLVERKLU2_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.hpp"
#include "klu2_defaults.hpp"
#include "klu2_analyze.hpp"
#include "klu2_factor.hpp"
#include "klu2_solve.hpp"
#include "klu2_free_symbolic.hpp"
#include "klu2_free_numeric.hpp"

namespace bddc {
  
  template <class SX> 
  class SolverKLU2 : public SolverBase<SX> 
  {
  public:

  SolverKLU2(int numRows,
	     int* rowBegin,
	     int* columns,
	     SX* values,
	     Teuchos::ParameterList & Parameters) :
    SolverBase<SX>(numRows, rowBegin, columns, values, Parameters)
  {
  }

  ~SolverKLU2() 
  {
    klu_free_symbolic<SX, int> (&m_Symbolic, &m_Common);
    klu_free_numeric<SX, int> (&m_Numeric, &m_Common);
  }

  int Initialize()
  {
    if (this->m_numRows == 0) return 0;
    int numRows = this->m_numRows;
    int* rowBegin = this->m_rowBegin;
    int* columns = this->m_columns;
    SX* values = this->m_values;
    klu_defaults<SX, int> (&m_Common);
    m_Symbolic = klu_analyze<SX, int> (numRows, rowBegin, columns, &m_Common);
    m_Numeric = klu_factor<SX, int> (rowBegin, columns, values, m_Symbolic, 
				     &m_Common);
    return 0;
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
    if (numRows == 0) return;
    memcpy(Sol, Rhs, numRows*NRHS*sizeof(SX));
    klu_solve<SX, int> (m_Symbolic, m_Numeric, numRows, NRHS, Sol, &m_Common);
  }
  
  bool MyExactSolver() 
  {
    return true;
  }
  
  private:
  klu_symbolic<SX, int> *m_Symbolic;
  klu_numeric<SX, int> *m_Numeric;
  klu_common<SX, int> m_Common;

  };
  
} // namespace bddc

#endif // BDDC_SOLVERKLU2_H
  
