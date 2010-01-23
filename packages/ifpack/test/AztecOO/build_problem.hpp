#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Comm.h"

#include "Epetra_LinearProblem.h"

Teuchos::RCP<Epetra_LinearProblem>
build_problem_mm(Teuchos::ParameterList& test_params, const Epetra_CrsMatrix* A, Epetra_MultiVector* b);

Teuchos::RCP< Epetra_LinearProblem >
build_problem(Teuchos::ParameterList& test_params, const Epetra_Comm& comm);


#endif

