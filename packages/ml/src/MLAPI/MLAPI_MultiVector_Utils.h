#ifndef MLAPI_DOUBLEVECTOR_UTILS_H
#define MLAPI_DOUBLEVECTOR_UTILS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Operator_Utils.h

\brief Suite of utilities for MLAPI::Operator objects.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "MLAPI_Error.h"
#include "MLAPI_MultiVector.h"

namespace MLAPI {

/*!
\file MLAPI_MultiVector_Utils.h

\brief Utilities for MultiVector's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Creates a new vector, x, such that x = y.
MultiVector Duplicate(const MultiVector& y);

//! Creates a new vector, x, such that x = y(:,v)
MultiVector Duplicate(const MultiVector& y, const int v);

//! Extracts a component from a vector.
MultiVector Extract(const MultiVector& y, const int v);

//! Redistributes the entry of a vector as a multivector.
MultiVector Redistribute(const MultiVector& y, const int NumEquations);

} // namespace MLAPI

#endif
