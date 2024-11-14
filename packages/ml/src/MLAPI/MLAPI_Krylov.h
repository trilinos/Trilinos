#ifndef MLAPI_KRYLOV
#define MLAPI_KRYLOV

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Krylov.h

\brief MLAPI interface to AztecOO's solvers.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"

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

#endif // ifdef MLAPI_KRYLOV
