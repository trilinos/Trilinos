#ifndef DISTORTION_METRIC_H
#define DISTORTION_METRIC_H

#include "mat2x2.h"
#include "projection.h"
#include <array>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

using TriPts = std::array<utils::Point, 3>; // points that define a triangle;

// returns a value in the range [1, infty), where 1 is the ideal
// element shape and infty is an inverted elemented
class DistortionMetric
{
  public:
    virtual ~DistortionMetric() {}

    // represents the derivatives of a single value
    static constexpr int DOT_VEC_LEN = 2;
    using DotVec                     = std::array<double, DOT_VEC_LEN>;

    // vector for holding the derivatives of all entries in TriPts
    using PtsDot = std::array<DotVec, 6>;

    // computes the distortion metric
    // W is the matrix describing the ideal element shape (see Knupp)
    virtual double get_value(const TriPts& pts, utils::impl::Mat2x2<double> w) = 0;

    // uses reverse-mode algorithmic differentiation to c
    // computes the derivative of the distorion metric wrt each input coordinates
    virtual void get_deriv(const TriPts& pts, utils::impl::Mat2x2<double> w, TriPts& derivs, const double qBar) = 0;

    // pts_dot is the vector the Hessian will be multiplied against
    // pts_bar_dot is the result of H * pts_dot
    // the outer array corresponds the the data in pts, the inner array
    // is the dual part
    virtual void get_value_rev_dot(const TriPts& pts, utils::impl::Mat2x2<double> w, const PtsDot& ptsDot,
                                   const double qIBar, const DotVec& qIBarDot, PtsDot& ptsBarDot) = 0;

    // computes the W matrix that specifies the given element as the ideal
    // element
    static utils::impl::Mat2x2<double> compute_w(const TriPts& pts)
    {
      return utils::impl::Mat2x2<double>{pts[1].x - pts[0].x, pts[2].x - pts[0].x, pts[1].y - pts[0].y,
                                         pts[2].y - pts[0].y};
    }

    // forward mode differentiation of computeW
    static utils::impl::Mat2x2<DotVec> compute_w_dot(const TriPts& pts, const PtsDot& ptsDot)
    {
      // return utils::impl::Mat2x2<double>{pts[1].x - pts[0].x, pts[2].x - pts[0].x,
      //                       pts[1].y - pts[0].y, pts[2].y - pts[0].y};
      utils::impl::Mat2x2<DotVec> wDot;
      auto& v1 = wDot(0, 0);
      auto& v2 = wDot(0, 1);
      auto& v3 = wDot(1, 0);
      auto& v4 = wDot(1, 1);

      v1.fill(0);
      v2.fill(0);
      v3.fill(0);
      v4.fill(0);

      for (int i = 0; i < DOT_VEC_LEN; ++i)
      {
        v1[i] = ptsDot[2][i] - ptsDot[0][i];
        v2[i] = ptsDot[4][i] - ptsDot[0][i];
        v3[i] = ptsDot[3][i] - ptsDot[1][i];
        v4[i] = ptsDot[5][i] - ptsDot[1][i];
      }

      return wDot;
    }

    // reverse-mode algorithmic differentiation version of computeW
    // pts_bar gets overwritten with the result
    static void compute_w_rev(const TriPts& pts, TriPts& ptsBar, const utils::impl::Mat2x2<double>& wBar)
    {
      // return utils::impl::Mat2x2<double>{pts[1].x - pts[0].x, pts[2].x - pts[0].x,
      //                       pts[1].y - pts[0].y, pts[2].y - pts[0].y};
      ptsBar[1].x = wBar(0, 0);
      ptsBar[0].x = -wBar(0, 0);
      ptsBar[2].x = wBar(0, 1);
      ptsBar[0].x += -wBar(0, 1);
      ptsBar[1].y = wBar(1, 0);
      ptsBar[0].y = -wBar(1, 0);
      ptsBar[2].y = wBar(1, 1);
      ptsBar[0].y += -wBar(1, 1);
    }

    // forward over reverse mode algorithmic differentiation
    static void compute_w_rev_dot(const TriPts& pts, const utils::impl::Mat2x2<DotVec> wBarDot, PtsDot& ptsBarDot)
    {
      // return utils::impl::Mat2x2<double>{pts[1].x - pts[0].x, pts[2].x - pts[0].x,
      //                       pts[1].y - pts[0].y, pts[2].y - pts[0].y};
      // pts_bar[1].x = W_bar(0, 0); pts_bar[0].x  = -W_bar(0, 0);
      // pts_bar[2].x = W_bar(0, 1); pts_bar[0].x += -W_bar(0, 1);
      // pts_bar[1].y = W_bar(1, 0); pts_bar[0].y  = -W_bar(1, 0);
      // pts_bar[2].y = W_bar(1, 1); pts_bar[0].y += -W_bar(1, 1);

      for (int i = 0; i < DOT_VEC_LEN; ++i)
      {
        ptsBarDot[2][i] = wBarDot(0, 0)[i];
        ptsBarDot[0][i] = -wBarDot(0, 0)[i];
        ptsBarDot[4][i] = wBarDot(0, 1)[i];
        ptsBarDot[0][i] += -wBarDot(0, 1)[i];

        ptsBarDot[3][i] = wBarDot(1, 0)[i];
        ptsBarDot[1][i] = -wBarDot(1, 0)[i];
        ptsBarDot[5][i] = wBarDot(1, 1)[i];
        ptsBarDot[1][i] += -wBarDot(1, 1)[i];
      }
    }
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
