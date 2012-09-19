#include <stk_mms/stk_mms.h>

#include <assert.h>

namespace stk
{
namespace mms
{
  // first partial derivatives for first order FADs
  double DX(const FAD_Type &f) { return f.dx(_X);}
  double DY(const FAD_Type &f) { return f.dx(_Y);}
  double DZ(const FAD_Type &f) { return f.dx(_Z);}
  double DT(const FAD_Type &f) { return f.dx(_T);}

  double D(const FAD_Type &f, const int i) { 
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }

  FAD_Type D(const FAD2_Type &f, const int i) { 
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }

  // first order partial derivatives for second order FADs
  FAD_Type DX2(const FAD2_Type &f) {return f.dx(_X);}
  FAD_Type DY2(const FAD2_Type &f) {return f.dx(_Y);}
  FAD_Type DZ2(const FAD2_Type &f) {return f.dx(_Z);}
  FAD_Type DT2(const FAD2_Type &f) {return f.dx(_T);}

  FAD_Type D2(const FAD2_Type &f, const int i) {
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }
}
}
