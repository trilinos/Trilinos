
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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "ShyLU_DDBDDC_config.h"

#if defined(HAVE_SHYLU_DDBDDC_SUPERLU)

#include "shylu_SolverBaseBDDC.hpp"
#include "shylu_errorBDDC.hpp"
#include "shylu_SolverSuperLU.hpp"
#include "shylu_SuperLU.hpp"

namespace bddc {

template <class SX> 
SolverSuperLU<SX>::SolverSuperLU(int numRows,
				 int* rowBegin,
				 int* columns,
				 SX* values,
				 Teuchos::ParameterList & Parameters) :
  SolverBase<SX>(numRows, rowBegin, columns, values, Parameters)
{
  m_solver = std::unique_ptr<SuperLU>(new SuperLU());
}

template <class SX> 
int SolverSuperLU<SX>::Initialize()
{
  int numRows = this->m_numRows;
  if (numRows == 0) return 0;
  int* rowBegin = this->m_rowBegin;
  int* columns = this->m_columns;
  SX* values =  this->m_values;
  int err = m_solver->initialize(numRows, rowBegin, columns, values);
  return err;
}

template <class SX> 
bool SolverSuperLU<SX>::IsDirectSolver()
{
  return true;
}

template <class SX> 
void SolverSuperLU<SX>::MySolve(int NRHS,
				SX* Rhs, 
				SX* Sol)
{
  int numRows = this->m_numRows;
  if (numRows == 0) return;
  m_solver->solve(NRHS, Rhs, Sol);
}

template <class SX> 
bool SolverSuperLU<SX>::MyExactSolver() 
  {
    return true;
  }

template class SolverSuperLU<double>;  
template class SolverSuperLU<float>;  

} // namespace bddc

#endif // HAVE_SHYLU_DDBDDC_SUPERLU

  
