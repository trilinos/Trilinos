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

