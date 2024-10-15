#ifndef ML_BASELINEARCOMBINATION_H
#define ML_BASELINEARCOMBINATION_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_BaseLinearCombination.h

\brief Base class for all operator overloading related operations.

\author Marzio Sala, SNL 9214.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

namespace MLAPI {

class Space;
class BaseOperator;
class MultiVector;

class BaseLinearCombination
{
public:
  virtual ~BaseLinearCombination() {};

  //! Returns the vector space of the underlying object.
  virtual const Space GetVectorSpace() const = 0;
  // Computes v += <operations>
  virtual void Update(MultiVector& v) const = 0;
  // Computes v = <operations>
  virtual void Set(MultiVector& v) const = 0;
};

} // namespace MLAPI

#endif // ifdef ML_BASELINEARCOMBINATION_H
