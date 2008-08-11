/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)
#include "MLAPI_Error.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiVector_Utils.h"

namespace MLAPI {

// ======================================================================
MultiVector Duplicate(const MultiVector& y)
{
  MultiVector x(y.GetVectorSpace(), y.GetNumVectors());
  x.Update(y);
  return(x);
}

// ======================================================================
MultiVector Duplicate(const MultiVector& y, const int v)
{
  if ((v < 0) || v >= y.GetNumVectors())
    ML_THROW("Wrong input parameter v (" +
             GetString(v) + ")", -1);
      
  // FIXME: use Extract
  MultiVector x(y.GetVectorSpace(), 1);
  for (int i = 0 ; i < x.GetMyLength() ; ++i)
    x(i) = y(i,v);

  return(x);
}

// ======================================================================
MultiVector Extract(const MultiVector& y, const int v)
{
  if ((v < 0) || v >= y.GetNumVectors())
    ML_THROW("Wrong input parameter v (" +
             GetString(v) + ")", -1);
      
  MultiVector x(y.GetVectorSpace(), y.GetRCPValues(v));

  return(x);
}

// ======================================================================
MultiVector Redistribute(const MultiVector& y, const int NumEquations)
{
  StackPush();

  if (y.GetMyLength() % NumEquations)
    ML_THROW("NumEquations does not divide MyLength()", -1);

  if (y.GetNumVectors() != 1)
    ML_THROW("Redistribute() works with single vectors only", -1);

  Space NewSpace(y.GetMyLength() / NumEquations);

  MultiVector y2(NewSpace, NumEquations);

  for (int i = 0 ; i < y2.GetMyLength() ; ++i)
    for (int j = 0 ; j < NumEquations ; ++j)
      y2(i, j) = y(j + NumEquations * i);

  StackPop();

  return(y2);
}

} // namespace MLAPI

#endif
