#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#include "ml_include.h"
#include "ml_epetra.h"
#include <iostream>
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"

namespace MLAPI {

/*!
\file MLAPI_Expressions.h

\brief Overloaded operators for DoubleVector's, Operator's, and InverseOpereator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Creates a new DoubleVector, defined as x + y
DoubleVector operator+(const DoubleVector& x, const DoubleVector& y)
{
  if (x.VectorSpace() != y.VectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  DoubleVector res(x.VectorSpace());
  res.Update(1.0, x, 1.0, y);
  return(res);
}

//! Creates a new DoubleVector, defined as x - y
DoubleVector operator-(const DoubleVector& x, const DoubleVector& y)
{
  if (x.VectorSpace() != y.VectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  DoubleVector res(x.VectorSpace());
  res.Update(1.0, x, -1.0, y);
  return(res);
}

//! Creates a new Operator, defined as A + B
Operator operator+(const Operator& A, const Operator& B)
{
  if (A.DomainSpace() != B.DomainSpace() ||
      A.RangeSpace() != B.RangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetData(),B.GetData(),ML_AplusB,MatrixType,1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

//! Creates a new Operator, defined as A - B
Operator operator-(const Operator& A, const Operator& B)
{
  if (A.DomainSpace() != B.DomainSpace() ||
      A.RangeSpace() != B.RangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetData(),B.GetData(),ML_AplusB,MatrixType,-1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

//! Creates a new Operator, defined as A * B
Operator operator*(const Operator& A, const Operator& B)
{
  if (A.DomainSpace() != B.RangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AtimesB = ML_Operator_Create(GetML_Comm());
  ML_2matmult(A.GetData(), B.GetData(), ML_AtimesB, MatrixType);
  Operator AtimesB(B.DomainSpace(),A.RangeSpace(), ML_AtimesB,true);
  return(AtimesB);
}

//! Creates a new DoubleVector, defined as x * alpha
DoubleVector
operator*(const DoubleVector& x, const double alpha) 
{
  DoubleVector y = Duplicate(x);
  y.Scale(alpha);
  return(y);
}

//! Creates a new DoubleVector y, such that y = x / alpha
DoubleVector
operator/(const DoubleVector& x, const double alpha) 
{
  if (alpha == 0.0)
    ML_THROW("Division by 0.0", -1);

  DoubleVector y = Duplicate(x);
  y.Scale(1.0 / alpha);
  return(y);
}

//! Creates a new DoubleVector y, such that y = A * x.
DoubleVector
operator*(const InverseOperator& A, const DoubleVector& x) 
{
  DoubleVector y(A.RangeSpace(), true);
  A.ApplyInverse(x,y);
  return(y);
}

//! Creates a new DoubleVector y, such that y = A * x.
DoubleVector
operator* (const Operator& A, const DoubleVector& x)
{
  DoubleVector y(A.RangeSpace(), true);
  A.Apply(x,y);
  return(y);
}

//! computes the dot product between x and y
double
operator* (const DoubleVector& x, const DoubleVector& y)
{
  return(x.DotProduct(y));
}

} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
