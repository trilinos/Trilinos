#ifndef MLAPI_AGGREGATION_H
#define MLAPI_AGGREGATION_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Workspace.h"

/*!
\file MLAPI_Aggregation.h

\brief Functions to create aggregation-based prolongator operators.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

using namespace std;

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

//! Builds the tentative prolongator using aggregation.
void GetPtent(const Operator& A, Teuchos::ParameterList& List,
              const MultiVector& ThisNS, 
              Operator& Ptent, MultiVector& NextNS);

//! Builds the tentative prolongator with default null space.
void GetPtent(const Operator& A, Teuchos::ParameterList& List, Operator& Ptent);

} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif // MLAPI_AGGREGATION_H
