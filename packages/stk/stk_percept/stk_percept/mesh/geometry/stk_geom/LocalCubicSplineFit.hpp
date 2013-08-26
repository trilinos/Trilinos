#ifndef LocalCubicSplineFit_hpp
#define LocalCubicSplineFit_hpp

#include <stk_percept/mesh/geometry/stk_geom/BSplineFit.hpp>

namespace stk {
  namespace geom {

    /// see P&T The Nurbs Book, 2nd Edition, section 9.3.4
    class LocalCubicSplineFit : public BSplineFit
    {
    public:
      /// create an OpenNURBS curve that fits the given input points
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T);
    };

  }
}

#endif
