/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)
#include "MLAPI_Error.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiVector_Utils.h"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_LinearCombinations.h"

namespace MLAPI {

// ====================================================================== 
LinearCombinationMixed operator+(const BaseLinearCombination& left, 
                                  const MultiVector& right)
{
  return(LinearCombinationMixed(left, right, 1.0));
}
// ====================================================================== 
LinearCombinationMixed operator+(const MultiVector& left, 
                                        const BaseLinearCombination& right)
{
  cout << "HERE" << endl;
  return(LinearCombinationMixed(right, left, 1.0));
}

// ====================================================================== 
LinearCombinationMixed operator-(const BaseLinearCombination& left, 
                                        const MultiVector& right)
{
  return(LinearCombinationMixed(left, right, -1.0));
}
// ====================================================================== 
LinearCombinationMixed operator-(const MultiVector& left, 
                                        const BaseLinearCombination& right)
{
  return(LinearCombinationMixed(right, left, -1.0));
}

// ====================================================================== 
LinearCombinationAdd operator+(const BaseLinearCombination& left, 
                                       const BaseLinearCombination& right)
{
  return(LinearCombinationAdd(left, right));
}

// ====================================================================== 
LinearCombinationAdd operator-(const BaseLinearCombination& left, 
                                       const BaseLinearCombination& right)
{
  return(LinearCombinationAdd(left, LinearCombinationScaled(right, -1.0)));
}

// ====================================================================== 
MultiVectorCombination operator+(const MultiVectorScaled& left, 
                                 const MultiVectorScaled& right)
{
  return(MultiVectorCombination(left.GetScalar(),
                                left.GetMultiVector(),
                                right.GetScalar(),
                                right.GetMultiVector()));
}

// ====================================================================== 
Residual operator+(const MultiVectorScaled& left, 
                   const BaseOperatorTimesMultiVector& right)
{
  return(Residual(left.GetScalar(), left.GetMultiVector(),
                  1.0, right.GetBaseOperator(), right.GetMultiVector()));
}

// ====================================================================== 
Residual operator+(const MultiVector& left, 
                   const BaseOperatorTimesMultiVector& right)
{
  return(Residual(1.0, left,
                  1.0, right.GetBaseOperator(), right.GetMultiVector()));
}

// ====================================================================== 
Residual operator-(const MultiVector& left, 
                   const BaseOperatorTimesMultiVector& right)
{
  return(Residual(1.0, left,
                  -1.0, right.GetBaseOperator(), right.GetMultiVector()));
}

// ====================================================================== 
MultiVectorCombination operator+(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  if (x.GetNumVectors() != y.GetNumVectors())
    ML_THROW("Number of vectors differs, " +
             GetString(x.GetNumVectors()) + " vs. " +
             GetString(y.GetNumVectors()), -1);

  return(MultiVectorCombination(1.0, x, 1.0, y));
}

// ====================================================================== 
MultiVectorCombination operator-(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  if (x.GetNumVectors() != y.GetNumVectors())
    ML_THROW("Number of vectors differs, " +
             GetString(x.GetNumVectors()) + " vs. " +
             GetString(y.GetNumVectors()), -1);

  return(MultiVectorCombination(1.0, x, -1.0, y));
}

// ====================================================================== 
MultiVector operator+(const MultiVector& x, const double alpha)
{
  MultiVector res(x.GetVectorSpace(), x.GetNumVectors());
  res = alpha;
  res.Update(1.0, x, 1.0);
  return(res);
}

// ====================================================================== 
MultiVector operator-(const MultiVector& x, const double alpha)
{
  MultiVector res(x.GetVectorSpace(), x.GetNumVectors());
  res = -alpha;
  res.Update(1.0, x, 1.0);
  return(res);
}

// ====================================================================== 
MultiVector operator+(const double alpha, const MultiVector& x)
{
  return(x + alpha);
}

// ====================================================================== 
MultiVector operator-(const double alpha, const MultiVector& x)
{
  return(x - alpha);
}

#if 0
// ====================================================================== 
MultiVector operator+= (const double alpha)
{
  return(x + alpha);
}

// ====================================================================== 
MultiVector operator-= (const double alpha)
{
  return(x - alpha);
}
#endif

// ====================================================================== 
Operator operator+(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetDomainSpace() ||
      A.GetRangeSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetML_Operator(),B.GetML_Operator(),ML_AplusB,
                  GetMatrixType(),1.);
  Operator AplusB(A.GetDomainSpace(),A.GetRangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// ====================================================================== 
Operator operator-(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetDomainSpace() ||
      A.GetRangeSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AplusB = ML_Operator_Create(GetML_Comm());
  ML_Operator_Add(A.GetML_Operator(),B.GetML_Operator(),ML_AplusB,
                  GetMatrixType(),-1.);
  Operator AplusB(A.GetDomainSpace(),A.GetRangeSpace(), ML_AplusB,true);
  return(AplusB);
}

// ====================================================================== 
Operator operator*(const Operator& A, const Operator& B)
{
  if (A.GetDomainSpace() != B.GetRangeSpace())
    ML_THROW("DomainSpace's or RangeSpace's are not compatible",-1);

  ML_Operator* ML_AtimesB = ML_Operator_Create(GetML_Comm());
  ML_2matmult(A.GetML_Operator(), B.GetML_Operator(), ML_AtimesB, 
              GetMatrixType());
  Operator AtimesB(B.GetDomainSpace(),A.GetRangeSpace(), ML_AtimesB,true);
  return(AtimesB);
}

// ====================================================================== 
Operator operator*(const Operator& A, const double alpha) 
{
  return(GetScaledOperator(A, alpha));
}

// ====================================================================== 
Operator operator*(const double alpha, const Operator& A) 
{
  return(GetScaledOperator(A, alpha));
}

// ====================================================================== 
MultiVector
operator*(const MultiVector& x, const double alpha) 
{
  MultiVector y = Duplicate(x);
  for (int v = 0 ; v < x.GetNumVectors() ; ++v)
    y.Scale(alpha, v);

  return(y);
}

// ====================================================================== 
MultiVector
operator/(const MultiVector& x, const double alpha) 
{
  if (alpha == 0.0)
    ML_THROW("Division by 0.0", -1);

  MultiVector y = Duplicate(x);
  for (int v = 0 ; v < x.GetNumVectors() ; ++v)
    y.Scale(1.0 / alpha, v);
  return(y);
}

// ====================================================================== 
BaseOperatorTimesMultiVector
operator*(const BaseOperator& A, const MultiVector& x) 
{
  return(BaseOperatorTimesMultiVector(A, x));
}

// ====================================================================== 
BaseOperatorTimesMultiVector
operator*(const BaseOperator& A, const BaseLinearCombination& LC) 
{
  MultiVector v(LC.GetVectorSpace());
  LC.Set(v);
  return(BaseOperatorTimesMultiVector(A, v));
}

// ====================================================================== 
double
operator* (const MultiVector& x, const MultiVector& y)
{
  if ((x.GetNumVectors() != 1) || (y.GetNumVectors() != 1))
    ML_THROW("operator* between MultiVectors requires NumVectors == 1", -1);

  return(x.DotProduct(y,0));
}

// ====================================================================== 
double
operator* (const MultiVector& y, const BaseLinearCombination& x)
{
  return(x * y);
}
  
// ====================================================================== 
double
operator* (const BaseLinearCombination& x, const MultiVector& y)
{
  MultiVector v(y.GetVectorSpace());
  x.Set(v);

  return(v.DotProduct(y,0));
}

// ====================================================================== 
double
operator* (const BaseLinearCombination& x, const BaseLinearCombination& y)
{
  MultiVector v(y.GetVectorSpace());
  MultiVector w(y.GetVectorSpace());
  x.Set(v);
  y.Set(w);

  return(v.DotProduct(w,0));
}

} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
