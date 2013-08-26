#ifndef BSplineFit_hpp
#define BSplineFit_hpp

#include <stk_percept/mesh/geometry/stk_geom/SplineFit.hpp>

namespace stk {
  namespace geom {

    /// using Piegl and Tiller "The NURBS Book, 2nd Edition" notation and equation numbers

    /// base class for "local" spline fit methods creating a B-spline
    /// (note: local methods depend only on a few points around each point
    /// and don't require a tridiagonal solve)
    class BSplineFit : public SplineFit
    {
    public:
      BSplineFit() : m_n(0), m_isRational(false), m_curve(0) {}

      /// create an OpenNURBS curve that fits the given input points
      virtual ON_Curve * fit(Vectors2D& input);

      // computes m_CV, m_U
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T) = 0;

      // given m_CV, m_U, gives the ON_NurbsCurve
      ON_Curve * create();
      void print();

      // control vertices, knots
      Vectors2D m_CV;
      std::vector<double> m_U;
      int m_n;
      bool m_isRational;
      ON_NurbsCurve *m_curve;
    };
  }
}
#endif
