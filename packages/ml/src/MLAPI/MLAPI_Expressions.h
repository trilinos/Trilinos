#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#include "ml_include.h"
#include "ml_epetra.h"
#include <iostream>
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"

namespace MLAPI {

/*!
\file MLAPI_Expressions.h

\brief Overloaded operators for MultiVector's, Operator's, and InverseOpereator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Creates a new MultiVector, defined as x + y
MultiVector operator+(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  MultiVector res(x.GetVectorSpace());
  res.Update(1.0, x, 1.0, y);
  return(res);
}

//! Creates a new MultiVector, defined as x - y
MultiVector operator-(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  MultiVector res(x.GetVectorSpace());
  res.Update(1.0, x, -1.0, y);
  return(res);
}

//! Creates a new Operator, defined as A + B
Operator operator+(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetDomainSpace() ||
      A.GetRangeSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetML_Operator(),B.GetML_Operator(),ML_AplusB,MatrixType,1.);
  Operator AplusB(A.GetDomainSpace(),A.GetRangeSpace(), ML_AplusB,true);
  return(AplusB);
}

//! Creates a new Operator, defined as A - B
Operator operator-(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetDomainSpace() ||
      A.GetRangeSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetML_Operator(),B.GetML_Operator(),ML_AplusB,MatrixType,-1.);
  Operator AplusB(A.GetDomainSpace(),A.GetRangeSpace(), ML_AplusB,true);
  return(AplusB);
}

//! Creates a new Operator, defined as A * B
Operator operator*(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AtimesB = ML_Operator_Create(GetML_Comm());
  ML_2matmult(A.GetML_Operator(), B.GetML_Operator(), ML_AtimesB, MatrixType);
  Operator AtimesB(B.GetDomainSpace(),A.GetRangeSpace(), ML_AtimesB,true);
  return(AtimesB);
}

//! Creates a new MultiVector, defined as x * alpha
MultiVector
operator*(const MultiVector& x, const double alpha) 
{
  MultiVector y = Duplicate(x);
  y.Scale(alpha);
  return(y);
}

//! Creates a new MultiVector y, such that y = x / alpha
MultiVector
operator/(const MultiVector& x, const double alpha) 
{
  if (alpha == 0.0)
    ML_THROW("Division by 0.0", -1);

  MultiVector y = Duplicate(x);
  y.Scale(1.0 / alpha);
  return(y);
}

//! Creates a new MultiVector y, such that y = A * x.
MultiVector
operator*(const InverseOperator& A, const MultiVector& x) 
{
  MultiVector y(A.GetRangeSpace(), true);
  A.ApplyInverse(x,y);
  return(y);
}

//! Creates a new MultiVector y, such that y = A * x.
MultiVector
operator* (const Operator& A, const MultiVector& x)
{
  MultiVector y(A.GetRangeSpace(), true);
  A.Apply(x,y);
  return(y);
}

//! computes the dot product between x and y
double
operator* (const MultiVector& x, const MultiVector& y)
{
  return(x.DotProduct(y));
}

} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
