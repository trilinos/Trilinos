#ifndef MLAPI_KRYLOV
#define MLAPI_KRYLOV

#include "ml_config.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_EpetraPreconditioner.h"
#include "MLAPI_DataBase.h"
#include "ml_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace MLAPI {

/*!
\file MLAPI_Krylov

\brief Simple wrapper to use MLAPI::Preconditioner with AztecOO

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

void Krylov(const Operator& A, const MultiVector& LHS,
            const MultiVector& RHS,
            const Preconditioner& Prec, KrylovDataBase& KDB)
{

  Epetra_LinearProblem Problem;
  ML_Operator* A_ML = A.GetML_Operator();
  ML_Epetra::RowMatrix A_Epetra(A_ML,&GetEpetra_Comm());

  Epetra_Vector LHS_Epetra(View,A_Epetra.OperatorDomainMap(),
                           (double*)&(LHS(0)));
  Epetra_Vector RHS_Epetra(View,A_Epetra.OperatorRangeMap(),
                           (double*)&(RHS(0)));

  // FIXME: this works only for Epetra-based operators
  Problem.SetOperator((Epetra_RowMatrix*)A.GetML_Operator()->data);
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
