#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "MLAPI_LinearCombinations.h"

namespace MLAPI {

/*!
\file MLAPI_Expressions.h

\brief Overloaded operators for MultiVector's, Operator's, and InverseOpereator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

// =============================================== //
// OPERATOR + BETWEEN VECTORS AND LINEAR OPERATORS //
// =============================================== //

// ======================================================================
//! Creates a new MultiVector, defined as x + y
// ======================================================================
#ifndef MLAPI_LC
MultiVectorCombination operator+(const MultiVector& x, const MultiVector& y);

LinearCombinationAdd operator+(const BaseLinearCombination& left,
                               const BaseLinearCombination& right);

LinearCombinationMixed operator+(const BaseLinearCombination& left,
                                 const MultiVector& right);

LinearCombinationMixed operator+(const MultiVector& left,
                                 const BaseLinearCombination& right);

// ======================================================================
//! Creates a new MultiVector, defined as alpha * x + beta * y
// ======================================================================

MultiVectorCombination operator+(const MultiVectorScaled& left,
                                 const MultiVectorScaled& right);

// ======================================================================
//! Creates a new MultiVector, defined as alpha * x + A * y
// ======================================================================

Residual operator+(const MultiVectorScaled& left,
                   const BaseOperatorTimesMultiVector& right);

// ======================================================================
//! Creates a new MultiVector, defined as x + A * y
// ======================================================================

Residual operator+(const MultiVector& left,
                   const BaseOperatorTimesMultiVector& right);

// ======================================================================
//! Creates a new MultiVector, defined as alpha * x + beta * y
// ======================================================================

MultiVectorCombination operator+(const MultiVectorScaled& left,
                                 const MultiVectorScaled& right);

// ======================================================================
//! Creates a new MultiVector, defined as alpha * x + A * y
// ======================================================================

Residual operator+(const MultiVectorScaled& left,
                   const BaseOperatorTimesMultiVector& right);

// ======================================================================
//! Creates a new MultiVector, defined as x + A * y
// ======================================================================

Residual operator+(const MultiVector& left,
                   const BaseOperatorTimesMultiVector& right);

// =============================================== //
// OPERATOR - BETWEEN VECTORS AND LINEAR OPERATORS //
// =============================================== //

// ======================================================================
//! Creates a new MultiVector, defined as x - y
// ======================================================================

MultiVectorCombination operator-(const MultiVector& x, const MultiVector& y);

LinearCombinationAdd operator-(const BaseLinearCombination& left,
                                const BaseLinearCombination& right);

LinearCombinationMixed operator-(const BaseLinearCombination& left,
                                const MultiVector& right);

LinearCombinationMixed operator-(const MultiVector& left,
                                const BaseLinearCombination& right);

// ======================================================================
//! Creates a new MultiVector, defined as x - A * y
// ======================================================================

Residual operator-(const MultiVector& left,
                   const BaseOperatorTimesMultiVector& right);
#else
inline
MultiVector operator+(const MultiVector& x, const MultiVector& y)
{
  MultiVector res(x.GetVectorSpace(), y.GetNumVectors());
  res.Update(1.0, x, 1.0, y);
  return(res);
}

inline
MultiVector operator-(const MultiVector& x, const MultiVector& y)
{
  MultiVector res(x.GetVectorSpace(), y.GetNumVectors());
  res.Update(1.0, x, -1.0, y);
  return(res);
}

inline
MultiVector operator*(const Operator& A, const MultiVector& y)
{
  MultiVector res(A.GetOperatorRangeSpace(), y.GetNumVectors());
  A.Apply(y, res);
  return(res);
}
#endif

// ======================================================================
//! Creates a new MultiVector, defined as x + alpha
// ======================================================================

MultiVector operator+(const MultiVector& x, const double alpha);

// ======================================================================
//! Creates a new MultiVector, defined as alpha + x
// ======================================================================

MultiVector operator+(const double alpha, const MultiVector& x);

// ======================================================================
//! Creates a new MultiVector, defined as x - alpha
// ======================================================================

MultiVector operator-(const MultiVector& x, const double alpha);

// ======================================================================
//! Creates a new MultiVector, defined as alpha - y
// ======================================================================

MultiVector operator-(const double alpha, const MultiVector& x);

#if 0
// ======================================================================
//! Adds a constant to a vector
// ======================================================================

MultiVector operator+= (const double alpha);

// ======================================================================
//! Subtracts a constant to a vector
// ======================================================================

MultiVector operator-= (const double alpha);
#endif

// ======================================================================
//! Creates a new Operator, defined as A + B
// ======================================================================

Operator operator+(const Operator& A, const Operator& B);

// ======================================================================
//! Creates a new Operator, defined as A - B
// ======================================================================

Operator operator-(const Operator& A, const Operator& B);

// ======================================================================
//! Creates a new Operator, defined as A * B
// ======================================================================

Operator operator*(const Operator& A, const Operator& B);

// ======================================================================
//! Creates a new Operator, defined as A * alpha
// ======================================================================

Operator operator*(const Operator& A, const double alpha);

// ======================================================================
//! Creates a new Operator, defined as alpha * A
// ======================================================================

Operator operator*(const double alpha, const Operator& A);

// ======================================================================
//! Creates a new Operator, defined as A / alpha
// ======================================================================

inline Operator operator/(const Operator& A, const double alpha)
{
  return(A * (1.0 / alpha));
}

// ======================================================================
//! Creates a new MultiVector, defined as x * alpha
// ======================================================================

MultiVector operator*(const MultiVector& x, const double alpha);

inline MultiVector operator*(const double alpha, const MultiVector&x)
{
  return(x * alpha);
}

// ======================================================================
//! Creates a new MultiVector y, such that y = x / alpha
// ======================================================================

MultiVector operator/(const MultiVector& x, const double alpha);

#ifndef MLAPI_LC
// ======================================================================
//! Creates a new MultiVector y, such that y = A * x.
// ======================================================================

BaseOperatorTimesMultiVector operator*(const BaseOperator& A, const MultiVector& x);

// ======================================================================
//! Creates a new MultiVector y, such that y = A * x (x is a BaseLinearCombination)
// ======================================================================

BaseOperatorTimesMultiVector operator*(const BaseOperator& A,
                                       const BaseLinearCombination& x);
#else
inline
MultiVector operator*(const BaseOperator& A, const MultiVector& x)
{
  MultiVector res(A.GetOperatorRangeSpace(), x.GetNumVectors());
  A.Apply(x, res);
  return(res);
}
#endif

// ======================================================================
//! Computes the dot product between the first vector in x and y
// ======================================================================

double operator* (const MultiVector& x, const MultiVector& y);

#ifndef MLAPI_LC
double operator* (const MultiVector& x, const BaseLinearCombination& y);

double operator* (const BaseLinearCombination& x, const MultiVector& y);

double operator* (const BaseLinearCombination& x, const BaseLinearCombination& y);
#endif

} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
