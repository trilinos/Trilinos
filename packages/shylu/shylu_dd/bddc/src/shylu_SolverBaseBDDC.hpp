
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

#ifndef BDDC_SOLVERBASE_H
#define BDDC_SOLVERBASE_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#include "Teuchos_ParameterList.hpp"  

namespace bddc {

enum MatrixTypeSolver{
  SymmetricMatrix = 0,
  NonSymmetricMatrix = 1
};

template <class SX> class SolverBase 
{
 protected: // variables
  int m_numRows{0};
  int *m_rowBegin{nullptr}, *m_columns{nullptr};
  SX *m_values{nullptr};
  Teuchos::ParameterList & m_Parameters;
  MPI_Comm* m_pComm;
  
 public: // functions
  SolverBase() { };
  
  SolverBase(int numRows,
	     int* rowBegin,
	     int* columns,
	     SX* values,
	     Teuchos::ParameterList & Parameters,
	     MPI_Comm* pComm = 0) :
  m_numRows(numRows),
    m_rowBegin(rowBegin),
    m_columns(columns),
    m_values(values),
    m_Parameters(Parameters),
    m_pComm(pComm)
  {
  }
  
  int getNumRows() const {
    return m_numRows;
  }
  virtual ~SolverBase() { };
  virtual int Initialize()=0;
  virtual bool IsDirectSolver()=0;
  virtual const int* getPermutation() const
  {
    return 0;
  }
  void Solve(int numRhs,
	     SX* rhs,
	     SX* sol)
  {
    MySolve(numRhs, rhs, sol);
  }
  
  bool exactSolver() {
    return MyExactSolver();
  }
  
 protected: // functions
  
  inline double wall_time () {
    // taken from Andrew Bradley's code
    static const double us = 1.0e6;
    timeval t;
    gettimeofday(&t, 0);
    return (t.tv_sec*us + t.tv_usec)/us;
  }

  virtual void MySolve(int NRHS,
		       SX* rhs,
		       SX* sol)=0;

  virtual bool MyExactSolver()=0;

};
 
} // namespace bddc

#endif // BDDC_SOLVERBASE_H
