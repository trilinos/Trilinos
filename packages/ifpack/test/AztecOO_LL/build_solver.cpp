/*
//@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "AztecOO.h"
#include "Ifpack.h"

#include "ParameterHelper.hpp"

Teuchos::RCP<AztecOO>
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Epetra_LinearProblem> problem)
{
  Teuchos::ParameterList azparams;
  if (test_params.isSublist("AztecOO")) {
    azparams = test_params.sublist("AztecOO");
  }

  Teuchos::RCP<AztecOO> solver = Teuchos::rcp(new AztecOO);

  solver->SetProblem(*problem);

  solver->SetParameters(azparams);

  std::string ifpack_precond("not specified");
  helper::GetParameter(test_params, "Ifpack_Preconditioner", ifpack_precond);
  if (ifpack_precond != "not specified") {
    Ifpack factory;
    Ifpack_Preconditioner* precond = factory.Create(ifpack_precond, problem->GetMatrix());

    if (test_params.isSublist("Ifpack")) {
      Teuchos::ParameterList& ifparams = test_params.sublist("Ifpack");
      precond->SetParameters(ifparams);
    }
    Teuchos::Time prec_time("precond");
    prec_time.start();
    precond->Initialize();
    precond->Compute();
    prec_time.stop();
    int my_proc = problem->GetMatrix()->RowMatrixRowMap().Comm().MyPID();
    if (my_proc == 0) {
      std::cout << "Time to initialize/compute preconditioner: " << prec_time.totalElapsedTime() << "s" << std::endl;
    }
    
    solver->SetPrecOperator(precond);
  }

  return solver;
}
#endif

