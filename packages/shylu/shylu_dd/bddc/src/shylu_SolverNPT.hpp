
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

#ifndef BDDC_SOLVERNPT_H
#define BDDC_SOLVERNPT_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.hpp"
#include "SolverNPT.h"

namespace bddc {
  
template <class SX, class SM, class LO, class GO> class SolverNPT : 
  public SolverBase<SX,SM,LO,GO>
{
public:
  SolverNPT(LO numRows,
	    LO* rowBegin,
	    LO* columns,
	    SX* values,
	    Teuchos::ParameterList & Parameters) :
  SolverBase<SX,SM,LO,GO>(numRows, rowBegin, columns, values, Parameters)
  {
    // NPT solver parameters
    std::vector<int> iparamsNPT(3);
    std::string astring = Parameters.get("Metis Option", "Node");
    if (astring == "Edge") iparamsNPT[npt::NPT_METIS_INDEX] = 
			     npt::NPT_METIS_EDGEND;
    else iparamsNPT[npt::NPT_METIS_INDEX] = npt::NPT_METIS_NODEND;
    iparamsNPT[npt::NPT_SCALE_INDEX] = Parameters.get("NPT Scale Option", 0);
    astring = Parameters.get("NPT Solve Option", "NPT");
    if (astring == "HTS") iparamsNPT[npt::NPT_SOLVE_INDEX] = 
			    npt::NPT_HTS_SOLVER;
    else iparamsNPT[npt::NPT_SOLVE_INDEX] = npt::NPT_NPT_SOLVER;
    m_solverNPT = new npt::SolverNPT<SX, SM, LO, GO>(iparamsNPT);
  }

  ~SolverNPT()
  {
    delete m_solverNPT;
  }

  int Initialize()
  {
    if (this->m_numRows == 0) return 0;
    std::string matrixTypeString = 
      this->m_Parameters.get("MatrixType", "Symmetric");
    if (matrixTypeString != "Symmetric") {
      std::cerr << "Error: SolverNPT only supports symmetric matrices\n";
      return 1;
    }
    int rval = m_solverNPT->Initialize(this->m_numRows, 
				       this->m_rowBegin, 
				       this->m_columns, 
				       this->m_values);
    return rval;
  }

  bool IsDirectSolver()
  {
    return true;
  }

  void MySolve(int NRHS,
	       SX* Rhs, 
	       SX* Sol)
  {
    m_solverNPT->MySolve(NRHS, Rhs, Sol);
  }

  bool MyExactSolver() 
  {
    return true;
  }

 private: // functions

 private: // data
  npt::SolverNPT<SX, SM, LO, GO>* m_solverNPT;

 };
  
} // namespace bddc

#endif // BDDC_SOLVERNPT_H
  
