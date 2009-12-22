#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include <string>
#include <sstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

#include "ParameterHelper.hpp"
#include "read_matrix.hpp"

Teuchos::RCP<Epetra_LinearProblem> build_problem_mm(Teuchos::ParameterList& test_params, Epetra_CrsMatrix* A)
{
  const Epetra_Map& rowmap = A->RowMap();

  Epetra_MultiVector* x = new Epetra_MultiVector(rowmap, 1);
  Epetra_MultiVector* b = new Epetra_MultiVector(rowmap, 1);
  x->PutScalar(1);

  A->Apply(*x, *b);
  x->PutScalar(0);

  Teuchos::RCP<Epetra_LinearProblem> problem = Teuchos::rcp(new Epetra_LinearProblem(A,x,b));

  return problem;
}

Teuchos::RCP< Epetra_LinearProblem >
build_problem(Teuchos::ParameterList& test_params,
              const Epetra_Comm& comm)
{
  Epetra_CrsMatrix* A;

  std::string mm_file("not specified");
  helper::GetParameter(test_params, "mm_file", mm_file);
  std::string hb_file("not specified");
  helper::GetParameter(test_params, "hb_file", hb_file);

  if (mm_file != "not specified") {
    if (comm.MyPID() == 0) {
      std::cout << "Matrix-Market file: " << mm_file << std::endl;
    }
    A = read_matrix_mm(mm_file, comm);
  }
  else if (hb_file != "not specified") {
    throw std::runtime_error("Harwell-Boeing not yet supported by test driver.");
    std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
  }
  else {
    throw std::runtime_error("No matrix file specified.");
  }

  return build_problem_mm(test_params, A);
}


#endif

