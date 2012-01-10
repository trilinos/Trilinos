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


#include <fei_macros.hpp>

#include <fei_Solver.hpp>

#include <fei_Matrix_Impl.hpp>
#include <fei_MatrixReducer.hpp>
#include <snl_fei_LinearSystem_FEData.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_utils.hpp>

#undef fei_file
#define fei_file "fei_Solver.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
int fei_Solver_solve(fei::LinearSystem* linearSystem,
			   fei::Matrix* preconditioningMatrix,
			   int numParams,
			   const char* const* solverParams,
			   int& iterationsTaken,
			   int& status)
{
  fei::SharedPtr<fei::Matrix> matrix = linearSystem->getMatrix();
  fei::Matrix_Impl<LinearSystemCore>* lscmatrix =
    dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matrix.get());

  fei::MatrixReducer* matred = dynamic_cast<fei::MatrixReducer*>(matrix.get());
  if (matred != NULL) {
    lscmatrix = dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matred->getTargetMatrix().get());
  }

  if (lscmatrix != NULL) {
    fei::SharedPtr<LinearSystemCore> linSysCore = lscmatrix->getMatrix();

    char** params = const_cast<char**>(solverParams);
    CHK_ERR( linSysCore->parameters(numParams, params) );

    CHK_ERR( linSysCore->launchSolver(status, iterationsTaken) );

    return(0);
  }

  snl_fei::LinearSystem_FEData* fedlinsys =
    dynamic_cast<snl_fei::LinearSystem_FEData*>(linearSystem);
  if (fedlinsys != NULL) {
    fei::SharedPtr<FiniteElementData> fedata = fedlinsys->getFiniteElementData();

    CHK_ERR( fedata->launchSolver(status, iterationsTaken) );

    return(0);
  }

  ERReturn(-1);
}

//----------------------------------------------------------------------------
int fei::Solver::solve(fei::LinearSystem* linearSystem,
			   fei::Matrix* preconditioningMatrix,
			   const fei::ParameterSet& parameterSet,
			   int& iterationsTaken,
			   int& status)
{
  int numParams = 0;
  const char** paramStrings = NULL;
  std::vector<std::string> stdstrings;
  fei::utils::convert_ParameterSet_to_strings(&parameterSet, stdstrings);
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

  int err = fei_Solver_solve(linearSystem, preconditioningMatrix,
		  numParams, paramStrings,
		  iterationsTaken, status);

  delete [] paramStrings;

  return(err);
}

