#ifndef REGULARIZED_DISTORTION_METRIC
#define REGULARIZED_DISTORTION_METRIC

#include "distortion_metric.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// the regularization of the distortion metric used by Escobar for mesh
// untangling
// Note: the regularization is only active for invalid elements
class RegularizedDistortionMetric : public DistortionMetric
{
  public:
    explicit RegularizedDistortionMetric(const double delta = 0)
      : m_delta(delta)
    {}

    void set_delta(const double delta) { m_delta = delta * delta; }

    double get_value(const TriPts& pts, utils::impl::Mat2x2<double> w) override
    {
      /*
      std::cout << "pt1 = " << pts[0] << std::endl;
      std::cout << "pt2 = " << pts[1] << std::endl;
      std::cout << "pt3 = " << pts[2] << std::endl;
      std::cout << "W = " << W << std::endl;
      */
      auto a = compute_w(pts); // calculation of A is the same as W, just different
                               // data
      // std::cout << "A = " << A << std::endl;
      inverse2x2(w); // W is now Winv
      auto s = a * w;
      // std::cout << "S = " << S << std::endl;
      auto num = norm_f(s);
      num      = num * num;

      auto den = det2x2(s);
      // std::cout << "num = " << num << std::endl;
      // std::cout << "den = " << den << std::endl;
      double delta2Val = 0;
      if (den < m_delta)
        delta2Val = m_delta * (m_delta - den);
      auto den2 = (den + std::sqrt(den * den + 4 * delta2Val));

      return num / den2;
    }

    void get_deriv(const TriPts& pts, utils::impl::Mat2x2<double> w, TriPts& derivs, const double qBar) override
    {
      auto a = compute_w(pts); // calculation of A is the same as W, just different
                               // data
      inverse2x2(w);           // W is now Winv
      auto s    = a * w;
      auto num  = norm_f(s);
      auto num2 = num * num;

      auto den         = det2x2(s);
      double delta2Val = 0;
      if (den < m_delta)
        delta2Val = m_delta * (m_delta - den);
      auto den2 = (den + std::sqrt(den * den + 4 * delta2Val));

      // auto q =  num2/den2;
      //----------------------------------
      //  reverse
      auto num2Bar = qBar / den2;
      auto den2Bar = -num2 / (den2 * den2) * qBar;

      auto denBar       = den2Bar + den2Bar * den / std::sqrt(den * den + 4 * delta2Val);
      auto delta2ValBar = 2 * den2Bar / std::sqrt(den * den + 4 * delta2Val);
      if (den < m_delta)
        denBar += -m_delta * delta2ValBar;

      utils::impl::Mat2x2<double> sBar;
      det2x2_rev(s, sBar, denBar);

      auto numBar = 2 * num * num2Bar;
      utils::impl::Mat2x2<double> sTmpBar;
      norm_f_rev(s, sTmpBar, numBar);
      sBar += sTmpBar;

      transpose(w); // this overwrites W, but thats ok because we don't need
                    // it again
      auto aBar = sBar * w;

      compute_w_rev(pts, derivs, aBar);
    }

    // pts_dot is the vector the Hessian will be multiplied against
    // pts_bar_dot is the result of H * pts_dot
    // the outer array corresponds the the data in pts, the inner array
    // is the dual part
    void get_value_rev_dot(const TriPts& pts, utils::impl::Mat2x2<double> w, const PtsDot& ptsDot, const double qBar,
                           const DotVec& qBarDot, PtsDot& ptsBarDot) override
    {
      using T = double;
      // const double q_bar = 1;

      auto a = compute_w(pts); // calculation of A is the same as W, just different
                               // data
      auto aDot = compute_w_dot(pts, ptsDot);
      inverse2x2(w); // W is now Winv
      auto s    = a * w;
      auto sDot = matmat_dot(aDot, w);

      auto num    = norm_f(s);
      auto numDot = norm_f_dot(s, sDot);

      auto num2 = num * num;
      DotVec num2Dot;
      for (unsigned int i = 0; i < DOT_VEC_LEN; ++i)
        num2Dot[i] = 2 * num * numDot[i];

      auto den            = det2x2(s);
      auto denDot         = det2x2_dot(s, sDot);
      double delta2Val    = 0;
      DotVec delta2ValDot = {0, 0};
      if (den < m_delta)
      {
        delta2Val = m_delta * (m_delta - den);
        for (int i = 0; i < DOT_VEC_LEN; ++i)
          delta2ValDot[i] = -m_delta * denDot[i];
      }
      auto valTmp = std::sqrt(den * den + 4 * delta2Val);
      auto den2   = (den + valTmp);
      std::array<T, DOT_VEC_LEN> den2Dot;
      for (int i = 0; i < DOT_VEC_LEN; ++i)
        den2Dot[i] = denDot[i] + den * denDot[i] / valTmp + 2 * delta2ValDot[i] / valTmp;

      // auto q =  num2/den2;
      //----------------------------------
      //  reverse
      auto num2Bar = qBar / den2;
      auto den2Bar = -num2 / (den2 * den2) * qBar;
      std::array<T, DOT_VEC_LEN> num2BarDot, den2BarDot;
      for (int i = 0; i < DOT_VEC_LEN; ++i)
      {
        num2BarDot[i] = -qBar * den2Dot[i] / (den2 * den2) + qBarDot[i] / den2;
        den2BarDot[i] = qBar * (-num2Dot[i] / (den2 * den2) + 2 * num2 * den2Dot[i] / (den2 * den2 * den2)) -
                        num2 * qBarDot[i] / (den2 * den2);
      }

      auto denBar       = den2Bar + den2Bar * den / valTmp;
      auto delta2ValBar = 2 * den2Bar / valTmp;

      std::array<T, DOT_VEC_LEN> denBarDot, delta2ValBarDot;
      auto denTmp = std::pow(valTmp, 3);
      for (int i = 0; i < DOT_VEC_LEN; ++i)
      {
        denBarDot[i] =
            den2BarDot[i] + den2BarDot[i] * den / valTmp +
            den2Bar * (denDot[i] / valTmp - den * den * denDot[i] / denTmp + -2 * den * delta2ValDot[i] / denTmp);
        delta2ValBarDot[i] = 2 * den2BarDot[i] / valTmp + -2 * den2Bar * den * denDot[i] / denTmp +
                             -4 * den2Bar * delta2ValDot[i] / denTmp;
      }

      if (den < m_delta)
      {
        denBar += -m_delta * delta2ValBar;
        for (int i = 0; i < DOT_VEC_LEN; ++i)
          denBarDot[i] += -m_delta * delta2ValBarDot[i];
      }

      utils::impl::Mat2x2<double> sBar;
      det2x2_rev(s, sBar, denBar);

      utils::impl::Mat2x2<DotVec> sBarDot;
      det2x2_rev_dot(s, sDot, denBar, denBarDot, sBarDot);

      auto numBar = 2 * num * num2Bar;
      DotVec numBarDot;
      for (int i = 0; i < DOT_VEC_LEN; ++i)
        numBarDot[i] = 2 * num2Bar * numDot[i] + 2 * num * num2BarDot[i];

      utils::impl::Mat2x2<double> sTmpBar;
      norm_f_rev(s, sTmpBar, numBar);
      sBar += sTmpBar;

      norm_f_rev_dot(s, sDot, numBar, numBarDot, sBarDot);

      transpose(w); // this overwrites W, but thats ok because we don't need
                    // it again
      // auto A_bar = S_bar * W;
      auto aBarDot = matmat_dot(sBarDot, w);

      // computeW_rev(pts, pts_bar, A_bar);
      compute_w_rev_dot(pts, aBarDot, ptsBarDot);
    }

  private:
    double m_delta;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
