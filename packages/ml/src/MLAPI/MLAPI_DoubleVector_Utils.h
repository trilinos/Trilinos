#ifndef MLAPI_DOUBLEVECTOR_UTILS_H
#define MLAPI_DOUBLEVECTOR_UTILS_H

#include "MLAPI_DoubleVector.h"

namespace MLAPI {

/*!
\file MLAPI_DoubleVector_Utils.h

\brief Utilities for DoubleVector's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Creates a new vector, x, such that x = y.
DoubleVector Duplicate(const DoubleVector& y)
{
  DoubleVector x(y.VectorSpace());
  x.Update(y);
  return(x);
}

} // namespace MLAPI

#endif
