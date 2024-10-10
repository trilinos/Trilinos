#ifndef MLAPI_PRECONDITIONER_H
#define MLAPI_PRECONDITIONER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_BaseOperator.h

\brief Base MLAPI operator.

\author Marzio Sala, SNL 9214.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "MLAPI_BaseObject.h"

namespace MLAPI {

class MultiVector;
class Space;

/*!
\class BaseOperator

\brief Base class for all MLAPI objects.

\author Marzio Sala, SNL 9214.

\date Last modified on Feb-05.

*/

class BaseOperator : public BaseObject {

public:

  //! Virtual destructor.
  virtual ~BaseOperator() {}

  //! Applies the operator to \c X, using \c Y as starting solution. Returns the solution in \c Y.
  virtual int Apply(const MultiVector& LHS, MultiVector& RHS) const = 0;

  //! Returns a copy of the domain space of \c this object.
  virtual const Space GetOperatorDomainSpace() const = 0;

  //! Returns a copy of the range space of \c this object.
  virtual const Space GetOperatorRangeSpace() const = 0;

};
} // namespace MLAPI

#endif
