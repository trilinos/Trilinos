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

} // namespace MLAPI

#endif
