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
      BSplineFit() : m_n(0), m_isRational(false),
                     m_isPeriodic(false), m_isClosed(false), m_curve(0) {}

      void setIsPeriodic(bool per) { m_isPeriodic=per; }
      bool getIsPeriodic() { return m_isPeriodic; }

      bool getIsClosed() { return m_isClosed; }

      /// create an OpenNURBS curve that fits the given input points
      virtual ON_Curve * fit(Vectors2D& input);
      void print();

      const Vectors2D& Q() { return m_Q; }

      const std::vector<double>& U() { return m_U; }
      const Vectors2D& CV() { return m_CV; }

    protected:
      // computes m_CV, m_U
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T) = 0;

      // given m_CV, m_U, gives the ON_NurbsCurve
      ON_Curve * create();

    protected:
      // control vertices, knots
      Vectors2D m_CV;
      Vectors2D m_Q;
      std::vector<double> m_U;
      int m_n;
      bool m_isRational;
      bool m_isPeriodic;
      bool m_isClosed;
      ON_NurbsCurve *m_curve;
    };
  }
}
#endif
