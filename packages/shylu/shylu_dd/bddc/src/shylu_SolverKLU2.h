
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

#ifndef SOLVERKLU2BDDC_H
#define SOLVERKLU2BDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.h"
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
    SolverBase<SX>(numRows, rowBegin, columns, values, Parameters),
      m_Solver() 
  {
  }

  ~SolverKLU2() 
  {
  
  }

  int Initialize()
  {
    if (this->m_numRows == 0) return 0;
    int numRows = this m_numRows;
    int* rowBegin = this->m_rowBegin;
    int* columns = this->m_columns;
    SX* values = this->m_values;
    klu_symbolic<SX, int> *Symbolic ;
    klu_numeric<SX, int> *Numeric ;
    klu_common<SX, int> Common ;
    klu_defaults<SX, int> (&Common) ;
    Symbolic = klu_analyze<SX, int> (numRows, rowBegin, columns, &Common) ;
    Numeric = klu_factor<SX, int> (rowBegin, columns, values, Symbolic, 
				   &Common) ;
    /*
    klu_solve<double, int> (Symbolic, Numeric, 5, 1, b, &Common) ;
    klu_free_symbolic<double, int> (&Symbolic, &Common) ;
    klu_free_numeric<double, int> (&Numeric, &Common) ;
    for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
    */    
  }

  bool IsDirectSolver()
  {
    return true;
  }
  
  void MySolve(int NRHS,
	       SX* Rhs, 
	       SX* Sol)
  {
    if (this->m_numRows == 0) return;
    
    const int m = this->m_numRows;
    /*
    value_type_matrix x(Sol, m, NRHS);
    value_type_matrix b(Rhs, m, NRHS);
    
    if (static_cast<int>(m_TempRhs.extent(0)) < m || 
	static_cast<int>(m_TempRhs.extent(1)) < NRHS)
      m_TempRhs = value_type_matrix("temp rhs", m, NRHS);      
    
    m_Solver.solve(x, b, m_TempRhs);
    */
  }
  
  bool MyExactSolver() 
  {
    return true;
  }
  
  private:

  };
  
} // namespace bddc

#endif // SOLVERKLU2BDDC_H
  
