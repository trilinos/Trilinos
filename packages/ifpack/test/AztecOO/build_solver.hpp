#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

Teuchos::RCP<AztecOO>
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Epetra_LinearProblem> problem);

#endif

