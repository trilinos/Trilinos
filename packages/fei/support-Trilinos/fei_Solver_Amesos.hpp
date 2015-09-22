/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_Solver_Amesos_h_
#define _fei_Solver_Amesos_h_


#include <fei_trilinos_macros.hpp>

#ifdef HAVE_FEI_AMESOS

#include <fei_Solver.hpp>

namespace Teuchos {
  class ParameterList;
}
class Amesos;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_VbrMatrix;

class Solver_Amesos : public fei::Solver {
 public:
  Solver_Amesos();
  virtual ~Solver_Amesos();

  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    const fei::ParameterSet& parameterSet,
	    int& iterationsTaken,
	    int& status);

  Teuchos::ParameterList& get_ParameterList();

 private:
  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    int numParams,
	    const char* const* solverParams,
	    int& iterationsTaken,
	    int& status);

  int parseParameters(int numParams,
		      const char*const* params);

  int solve_private(Epetra_CrsMatrix* A,
		    Epetra_MultiVector* x,
		    Epetra_MultiVector* b,
		    fei::Matrix* preconditioningMatrix,
		    int numParams,
		    const char* const* solverParams,
		    int& iterationsTaken,
		    int& status);

  int solve_private(Epetra_VbrMatrix* A,
		    Epetra_MultiVector* x,
		    Epetra_MultiVector* b,
		    fei::Matrix* preconditioningMatrix,
		    int numParams,
		    const char* const* solverParams,
		    int& iterationsTaken,
		    int& status);

 private:
  double tolerance_;
  int maxIters_;
  Amesos* amesos_factory_;
  Amesos_BaseSolver* amesos_solver_;
  Epetra_LinearProblem* epetra_linearproblem_;
  Teuchos::ParameterList* paramlist_;
}; //class Solver_Amesos

#endif

#endif
