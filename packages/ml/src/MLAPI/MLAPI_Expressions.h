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

// x + y
DoubleVector operator+(const DoubleVector& x, const DoubleVector& y)
{
  assert (x.VectorSpace() == y.VectorSpace());
  DoubleVector res(x.VectorSpace());
  for (int i = 0 ; i < x.MyLength() ; ++i)
    res(i) = x(i) + y(i);
  return(res);
}

// x - y
DoubleVector operator-(const DoubleVector& x, const DoubleVector& y)
{
  assert (x.VectorSpace() == y.VectorSpace());
  DoubleVector res(x.VectorSpace());
  for (int i = 0 ; i < x.MyLength() ; ++i)
    res(i) = x(i) - y(i);
  return(res);
}

// A + B
Operator operator+(const Operator& A, const Operator& B)
{
  assert (A.DomainSpace() == B.DomainSpace());
  assert (A.RangeSpace() == B.RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetOperator(),ML_AplusB,MatrixType,1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// A - B
Operator operator-(const Operator& A, const Operator& B)
{
  assert (A.DomainSpace() == B.DomainSpace());
  assert (A.RangeSpace() == B.RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetOperator(),ML_AplusB,MatrixType,-1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// A * B
Operator operator*(const Operator& A, const Operator& B)
{
  // FIXME check on spaces
  ML_Operator* ML_AtimesB = ML_Operator_Create(GetMLComm());
  ML_2matmult(A.GetOperator(), B.GetOperator(), ML_AtimesB, MatrixType);
  Operator AtimesB(B.DomainSpace(),A.RangeSpace(), ML_AtimesB,true);
  return(AtimesB);
}

#if 0
// A + beta B
Operator operator+(const Operator& A, const BaseObjectMult<double,Operator>& B)
{
  assert (A.DomainSpace() == B.GetRight().DomainSpace());
  assert (A.RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetRight().GetOperator(),
                  ML_AplusB,MatrixType, B.GetLeft());
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// A - beta B
Operator operator-(const Operator& A, const BaseObjectMult<double,Operator>& B)
{
  assert (A.DomainSpace() == B.GetRight().DomainSpace());
  assert (A.RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetRight().GetOperator(),
                  ML_AplusB,MatrixType, - B.GetLeft());
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// alpha A + beta B
Operator operator+(const BaseObjectMult<double,Operator>& A, 
                    const BaseObjectMult<double,Operator>& B)
{
  assert (A.GetRight().DomainSpace() == B.GetRight().DomainSpace());
  assert (A.GetRight().RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add2(A.GetRight().GetOperator(),B.GetRight().GetOperator(),
                   ML_AplusB,MatrixType, A.GetLeft(),B.GetLeft());
  Operator AplusB(A.GetRight().DomainSpace(), A.GetRight().RangeSpace(),
                  ML_AplusB,true);
  return(AplusB);
}

// alpha A - beta B
Operator operator-(const BaseObjectMult<double,Operator>& A, 
                    const BaseObjectMult<double,Operator>& B)
{
  assert (A.GetRight().DomainSpace() == B.GetRight().DomainSpace());
  assert (A.GetRight().RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add2(A.GetRight().GetOperator(),B.GetRight().GetOperator(),
                   ML_AplusB,MatrixType, A.GetLeft(),-B.GetLeft());
  Operator AplusB(A.GetRight().DomainSpace(), A.GetRight().RangeSpace(),
                  ML_AplusB,true);
  return(AplusB);
}
#endif

// LHS = S / rhs
DoubleVector
operator/(const InverseOperator& A, const DoubleVector& RHS) {
  // FIXME: check on spaces
  DoubleVector LHS(RHS.VectorSpace());
  LHS = 0.0;
  A.ApplyInverse(RHS,LHS);
  return(LHS);
}

// LHS = A * RHS
DoubleVector
operator* (const Operator& A, const DoubleVector& RHS)
{
  DoubleVector LHS(RHS.VectorSpace());
  A.Apply(RHS,LHS);
  return(LHS);
}

// x * y
double
operator* (const DoubleVector& Left, const DoubleVector& Right)
{
  return(Left.DotProduct(Right));
}

} // namespace MLAPI
#endif // if ML_EXPRESSIONS_H
