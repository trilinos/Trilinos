#ifndef MLAPI_KRYLOV
#define MLAPI_KRYLOV

#include "ml_config.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_EpetraPreconditioner.h"
#include "ml_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace MLAPI {

void Krylov(const Operator& A, const DoubleVector& LHS,
            const DoubleVector& RHS,
            const Preconditioner& Prec,
            Teuchos::ParameterList& List)
{

  Epetra_LinearProblem Problem;
  ML_Operator* A_ML = A.GetOperator();
  ML_Epetra::RowMatrix A_Epetra(A_ML,&GetEpetraComm());

  Epetra_Vector LHS_Epetra(View,A_Epetra.OperatorDomainMap(),
                           (double*)&(LHS(0)));
  Epetra_Vector RHS_Epetra(View,A_Epetra.OperatorRangeMap(),
                           (double*)&(RHS(0)));

  // FIXME: this works only for Epetra-based operators
  Problem.SetOperator((Epetra_RowMatrix*)A.GetOperator()->data);
  //Problem.SetOperator(&A_Epetra);
  Problem.SetLHS(&LHS_Epetra);
  Problem.SetRHS(&RHS_Epetra);

  AztecOO solver(Problem);

  EpetraPreconditioner Prec_Epetra(A_Epetra.OperatorDomainMap(),Prec);

  solver.SetPrecOperator(&Prec_Epetra);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 16);
  solver.Iterate(500, 1e-5);

}

} // namespace MLAPI

#endif // ifdef MLAPI_KRYLOV
