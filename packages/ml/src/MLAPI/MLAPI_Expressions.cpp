#include "ml_include.h"
#if defined(HAVE_ML_MLAPI)
#include "MLAPI_Error.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_Expressions.h"

namespace MLAPI {

// ====================================================================== 
MultiVector operator+(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  if (x.GetNumVectors() != y.GetNumVectors())
    ML_THROW("Number of vectors differs, " +
             GetString(x.GetNumVectors()) + " vs. " +
             GetString(y.GetNumVectors()), -1);

  MultiVector res(x.GetVectorSpace());
  res.Update(1.0, x, 1.0, y);
  return(res);
}

// ====================================================================== 
MultiVector operator-(const MultiVector& x, const MultiVector& y)
{
  if (x.GetVectorSpace() != y.GetVectorSpace())
    ML_THROW("VectorSpace's are not compatible",-1);

  if (x.GetNumVectors() != y.GetNumVectors())
    ML_THROW("Number of vectors differs, " +
             GetString(x.GetNumVectors()) + " vs. " +
             GetString(y.GetNumVectors()), -1);

  MultiVector res(x.GetVectorSpace());
  res.Update(1.0, x, -1.0, y);
  return(res);
}

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
  y.Scale(alpha);
  return(y);
}

// ====================================================================== 
MultiVector
operator/(const MultiVector& x, const double alpha) 
{
  if (alpha == 0.0)
    ML_THROW("Division by 0.0", -1);

  MultiVector y = Duplicate(x);
  y.Scale(1.0 / alpha);
  return(y);
}

// ====================================================================== 
MultiVector
operator*(const BaseOperator& A, const MultiVector& x) 
{
  MultiVector y(A.GetOperatorRangeSpace(), true);
  A.Apply(x,y);
  return(y);
}

// ====================================================================== 
double
operator* (const MultiVector& x, const MultiVector& y)
{
  if ((x.GetNumVectors() != 1) || (y.GetNumVectors() != 1))
    ML_THROW("operator* between MultiVectors requires NumVectors == 1", -1);

  return(x.DotProduct(y,0));
}

} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
