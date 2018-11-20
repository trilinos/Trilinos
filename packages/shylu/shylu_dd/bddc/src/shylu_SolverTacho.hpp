
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

#ifndef BDDC_SOLVERTACHO_H
#define BDDC_SOLVERTACHO_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.hpp"
#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

namespace bddc {
  
  template <class SX> 
  class SolverTacho : public SolverBase<SX> 
  {
  public:
    typedef Tacho::Solver<SX,Kokkos::DefaultHostExecutionSpace> solver_type;

    typedef Tacho::ordinal_type ordinal_type;
    typedef Tacho::size_type size_type;
    typedef typename solver_type::value_type value_type;

    typedef typename solver_type::ordinal_type_array ordinal_type_array;
    typedef typename solver_type::size_type_array size_type_array;
    typedef typename solver_type::value_type_array value_type_array;
    typedef typename solver_type::value_type_matrix value_type_matrix;

    SolverTacho(int numRows,
                int* rowBegin,
                int* columns,
                SX* values,
                Teuchos::ParameterList & Parameters) :
      SolverBase<SX>(numRows, rowBegin, columns, values, Parameters),
      m_Solver() 
    {
    }

    ~SolverTacho() 
    {
      m_Solver.release();
    }

    int Initialize()
    {
      if (this->m_numRows == 0) return 0;

      // from parameterlist (maybe later)
      //      const int verbose = 1;
      const int verbose = 0;

      //  -- by default it sets 4, but usually the solver performs better with 8
      const int max_num_superblocks = 8; 

      m_Solver.setVerbose(verbose);
      m_Solver.setMaxNumberOfSuperblocks(max_num_superblocks);
      
      const int m = this->m_numRows;
      size_type_array    ap((size_type*)   this->m_rowBegin, m+1);
      ordinal_type_array aj((ordinal_type*)this->m_columns,  ap(m));
      value_type_array   ax((value_type*)  this->m_values,   ap(m));
      
      m_Solver.analyze(m, ap, aj);
      m_Solver.factorize(ax);

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
      if (this->m_numRows == 0) return;

      const int m = this->m_numRows;

      value_type_matrix x(Sol, m, NRHS);
      value_type_matrix b(Rhs, m, NRHS);
      
      if (static_cast<int>(m_TempRhs.extent(0)) < m || 
          static_cast<int>(m_TempRhs.extent(1)) < NRHS)
        m_TempRhs = value_type_matrix("temp rhs", m, NRHS);      
      
      m_Solver.solve(x, b, m_TempRhs);
    }

    bool MyExactSolver() 
    {
      return true;
    }

  private:
    solver_type m_Solver;
    value_type_matrix m_TempRhs;

  };
  
} // namespace bddc

#endif // BDDC_SOLVERTACHO_H
  
