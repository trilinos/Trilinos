#ifndef MLAPI_KRYLOV
#define MLAPI_KRYLOV

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

namespace Teuchos {
  class List;
}

namespace MLAPI {

class Operator;
class BaseOperator;
class MultiVector;

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
