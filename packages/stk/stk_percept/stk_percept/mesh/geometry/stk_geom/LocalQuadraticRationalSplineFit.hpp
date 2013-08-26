#ifndef LocalQuadraticRationalSplineFit_hpp
#define LocalQuadraticRationalSplineFit_hpp

#include <stk_percept/mesh/geometry/stk_geom/BSplineFit.hpp>

namespace stk {
  namespace geom {

    class LocalQuadraticRationalSplineFit : public BSplineFit
    {
    public:
      // computes m_CV, m_U
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T);
    };

  }
}

#endif
