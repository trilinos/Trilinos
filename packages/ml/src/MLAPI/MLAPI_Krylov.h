#ifndef MLAPI_KRYLOV
#define MLAPI_KRYLOV

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "MLAPI_Error.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_EpetraBaseOperator.h"
#include "ml_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace MLAPI {

/*!
\file MLAPI_Krylov

\brief Simple wrapper to use MLAPI::BaseOperator's with AztecOO

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

void Krylov(const Operator& A, const MultiVector& LHS,
            const MultiVector& RHS, const BaseOperator& Prec, 
            Teuchos::ParameterList& List);

} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif // ifdef MLAPI_KRYLOV
