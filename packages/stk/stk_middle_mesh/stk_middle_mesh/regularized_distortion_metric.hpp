#ifndef REGULARIZED_DISTORTION_METRIC
#define REGULARIZED_DISTORTION_METRIC

#include "distortion_metric.hpp"
#include <iomanip>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// the regularization of the distortion metric used by Escobar for mesh
// untangling
// Note: the regularization is only active for invalid elements
template <typename T>
class RegularizedDistortionMetric : public DistortionMetric<T>
{
  public:
    explicit RegularizedDistortionMetric(const double gamma = 0.1, double denominatorMin=1)
      : m_gamma(gamma)
      , m_denominatorMin(denominatorMin)
    {}

    template <typename T2>
    RegularizedDistortionMetric(const RegularizedDistortionMetric<T2>& other) :
      m_gamma(other.m_gamma),
      m_denominatorMin(other.m_denominatorMin)
    {}

    static constexpr int DOT_VEC_LEN = DistortionMetric<T>::DOT_VEC_LEN;

    template <typename T2>
    using DotVec = typename DistortionMetric<T>::template DotVec<T2>;

    template <typename T2>
    using PtsDot = typename DistortionMetric<T>::template PtsDot<T2>;

    void set_min_denominator(double denominatorMin) override { m_denominatorMin = denominatorMin; }

    //void set_delta(const double delta) { m_delta = delta * delta; }

    T get_value(const TriPts<T>& pts, utils::impl::Mat2x2<double> w) override;

    // uses reverse-mode algorithmic differentiation to c
    // computes the derivative of the distorion metric wrt each input coordinates
    void get_deriv(const TriPts<T>& pts, utils::impl::Mat2x2<double> w, TriPts<T>& derivs, const T qBare) override;

    // pts_dot is the vector the Hessian will be multiplied against
    // pts_bar_dot is the result of H * pts_dot
    // the outer array corresponds the the data in pts, the inner array
    // is the dual part
    void get_value_rev_dot(const TriPts<T>& pts, utils::impl::Mat2x2<double> w, const PtsDot<T>& ptsDot,
                           const T qIBar, const DotVec<T>& qIBarDot, PtsDot<T>& ptsBarDot) override;

    double compute_denominator(const TriPts<double>& pts, utils::impl::Mat2x2<double> w) override
    {
      auto a = this->compute_w(pts);      
      inverse2x2(w); // W is now Winv
      auto s = a * w;   
      auto den = det2x2(s);

      return den;
    }

  private:
    double m_gamma;
    double m_denominatorMin;

    template <typename T2>
    friend class RegularizedDistortionMetric;
};


template <typename T>
T RegularizedDistortionMetric<T>::get_value(const TriPts<T>& pts, utils::impl::Mat2x2<double> w) 
{
  auto a = this->compute_w(pts); // calculation of A is the same as W, just different
                            // data
  inverse2x2(w); // W is now Winv

  auto s = a * w;
  auto num = norm_f(s);
  num      = num * num;

  auto den = det2x2(s);

  double delta2Val = 0;
  if (m_denominatorMin < m_gamma)
  {
    delta2Val = m_gamma*(m_gamma - m_denominatorMin);
  }

  auto den2 = (den + std::sqrt(den * den + 4 * delta2Val));

  return num / den2;
}

template <typename T>
void RegularizedDistortionMetric<T>::get_deriv(const TriPts<T>& pts, utils::impl::Mat2x2<double> w, TriPts<T>& derivs, const T qBar)
{
  auto a = this->compute_w(pts); // calculation of A is the same as W, just different
                                 // data
  inverse2x2(w);           // W is now Winv

  auto s    = a * w;
  auto num  = norm_f(s);
  auto num2 = num * num;
  

  auto den         = det2x2(s);

  double delta2Val = 0;
  if (m_denominatorMin < m_gamma)
  {
    delta2Val = m_gamma*(m_gamma - m_denominatorMin);
  }  
  auto den2 = (den + std::sqrt(den * den + 4 * delta2Val));
  
  //----------------------------------
  //  reverse
  
  auto num2Bar = qBar / den2;
  auto den2Bar = -num2 / (den2 * den2) * qBar;

  auto denBar       = den2Bar + den2Bar * den / std::sqrt(den * den + 4 * delta2Val);

  utils::impl::Mat2x2<T> sBar;
  det2x2_rev(s, sBar, denBar);
  
  auto numBar = 2 * num * num2Bar;
  
  utils::impl::Mat2x2<T> sTmpBar;
  norm_f_rev(s, sTmpBar, numBar);
  sBar += sTmpBar;

  transpose(w); // this overwrites W, but thats ok because we don't need
                // it again
  auto aBar = sBar * w;

  this->compute_w_rev(pts, derivs, aBar);
}

template <typename T>
void RegularizedDistortionMetric<T>::get_value_rev_dot(const TriPts<T>& pts, utils::impl::Mat2x2<double> w, 
                                                       const PtsDot<T>& ptsDot, const T qBar,
                                                       const DotVec<T>& qBarDot, PtsDot<T>& ptsBarDot)
{
  // const double q_bar = 1;

  auto a = this->compute_w(pts); // calculation of A is the same as W, just different
                            // data
  auto aDot = this->compute_w_dot(pts, ptsDot);
  inverse2x2(w); // W is now Winv
  auto s    = a * w;
  auto sDot = matmat_dot(aDot, w);

  auto num    = norm_f(s);
  auto numDot = norm_f_dot(s, sDot);

  auto num2 = num * num;
  DotVec<T> num2Dot;
  for (unsigned int i = 0; i < DOT_VEC_LEN; ++i)
    num2Dot[i] = 2 * num * numDot[i];

  auto den            = det2x2(s);
  auto denDot         = det2x2_dot(s, sDot);

  double delta2Val = 0;
  if (m_denominatorMin < m_gamma)
  {
    delta2Val = m_gamma*(m_gamma - m_denominatorMin);
  }

  auto valTmp = std::sqrt(den * den + 4 * delta2Val);

  DotVec<T> valTmpDot;
  for (int i=0; i < DOT_VEC_LEN; ++i)
    valTmpDot[i] = den * denDot[i]/valTmp;

  auto den2   = (den + valTmp);
  std::array<T, DOT_VEC_LEN> den2Dot;
  for (int i = 0; i < DOT_VEC_LEN; ++i)
    den2Dot[i] = denDot[i] + valTmpDot[i];

  // auto q =  num2/den2;
  //----------------------------------
  //  reverse
  auto num2Bar = qBar / den2;
  auto den2Bar = -num2 / (den2 * den2) * qBar;
  std::array<T, DOT_VEC_LEN> num2BarDot, den2BarDot;
  for (int i = 0; i < DOT_VEC_LEN; ++i)
  {
    num2BarDot[i] = -qBar * den2Dot[i] / (den2 * den2) + qBarDot[i] / den2;
    den2BarDot[i] = -qBar * num2Dot[i] / (den2 * den2) + 2 * qBar * num2 * den2Dot[i] / (den2 * den2 * den2) -
                    num2 * qBarDot[i] / (den2 * den2);                     
  }

  auto denBar       = den2Bar + den2Bar * den / valTmp;

  std::array<T, DOT_VEC_LEN> denBarDot;      
  for (int i = 0; i < DOT_VEC_LEN; ++i)
  {
    T t1 = den2Bar * den;
    T t1Dot = den2BarDot[i] * den + den2Bar * denDot[i];
    denBarDot[i] = den2BarDot[i] + (t1Dot * valTmp - t1 * valTmpDot[i])/(valTmp*valTmp);
  }

  utils::impl::Mat2x2<T> sBar;
  det2x2_rev(s, sBar, denBar);

  utils::impl::Mat2x2<DotVec<T>> sBarDot;
  det2x2_rev_dot(s, sDot, denBar, denBarDot, sBarDot);

  auto numBar = 2 * num * num2Bar;
  DotVec<T> numBarDot;
  for (int i = 0; i < DOT_VEC_LEN; ++i)
  {
    numBarDot[i] = 2 * num2Bar * numDot[i] + 2 * num * num2BarDot[i];
  }

  utils::impl::Mat2x2<T> sTmpBar;
  norm_f_rev(s, sTmpBar, numBar);
  sBar += sTmpBar;

  norm_f_rev_dot(s, sDot, numBar, numBarDot, sBarDot);

  transpose(w); // this overwrites W, but thats ok because we don't need
                // it again
  // auto A_bar = S_bar * W;

  auto aBarDot = matmat_dot(sBarDot, w);

  // computeW_rev(pts, pts_bar, A_bar);
  this->compute_w_rev_dot(pts, aBarDot, ptsBarDot);   
}


} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
