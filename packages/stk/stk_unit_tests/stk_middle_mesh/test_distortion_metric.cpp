#include "gtest/gtest.h"

#include "stk_middle_mesh/distortion_metric.hpp"
#include "stk_middle_mesh/regularized_distortion_metric.hpp"

#include <fstream>

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh::impl;

TriPts<double> convert_pts(std::array<double, 6> pts)
{
  return {utils::Point(pts[0], pts[1]), utils::Point(pts[2], pts[3]), utils::Point(pts[4], pts[5])};
}

std::array<double, 6> convert_pts(const TriPts<double>& pts)
{
  return {pts[0].x, pts[0].y, pts[1].x, pts[1].y, pts[2].x, pts[2].y};
}

utils::impl::Mat2x2<double> compute_w_fd(std::array<double, 6> pts, int idx)
{
  auto ptsAos = convert_pts(pts);
  auto w      = DistortionMetric<double>::compute_w(ptsAos);

  double eps = 1e-7;
  pts[idx] += eps;
  auto wPert = DistortionMetric<double>::compute_w(convert_pts(pts));
  auto dWDpt = (wPert - w) / eps;
  pts[idx] -= eps;

  return dWDpt;
}

TEST(DistortionMetric, ComputeW)
{
  TriPts<double> pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric<double>::compute_w(pts);

  EXPECT_FLOAT_EQ(w(0, 0), 1.0);
  EXPECT_FLOAT_EQ(w(0, 1), 0.5);
  EXPECT_FLOAT_EQ(w(1, 0), 0.0);
  EXPECT_FLOAT_EQ(w(1, 1), std::sqrt(3) / 2.0);

  TriPts<double> ptsBar;
  utils::impl::Mat2x2<double> wBar{1, 2, 3, 4};
  DistortionMetric<double>::compute_w_rev(pts, ptsBar, wBar);

  // compute finite difference
  std::vector<utils::impl::Mat2x2<double>> dWDptsFd(6);
  std::array<double, 6> ptsArr{pts[0].x, pts[0].y, pts[1].x, pts[1].y, pts[2].x, pts[2].y};
  for (int i = 0; i < 6; ++i)
    dWDptsFd[i] = compute_w_fd(ptsArr, i);

  std::array<double, 6> ptsBarFd;
  for (int i = 0; i < 6; ++i)
  {
    ptsBarFd[i] = 0;
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        ptsBarFd[i] += dWDptsFd[i](j, k) * wBar(j, k);
  }

  EXPECT_FLOAT_EQ(ptsBar[0].x, ptsBarFd[0]);
  EXPECT_FLOAT_EQ(ptsBar[0].y, ptsBarFd[1]);
  EXPECT_FLOAT_EQ(ptsBar[1].x, ptsBarFd[2]);
  EXPECT_FLOAT_EQ(ptsBar[1].y, ptsBarFd[3]);
  EXPECT_FLOAT_EQ(ptsBar[2].x, ptsBarFd[4]);
  EXPECT_FLOAT_EQ(ptsBar[2].y, ptsBarFd[5]);
}

TEST(DistortionMetric, ComputeW_dot)
{
  TriPts<double> pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric<double>::compute_w(pts);

  std::array<std::array<double, 2>, 6> ptsDot;
  for (int i = 0; i < 6; ++i)
    ptsDot[i] = std::array<double, 2>{2.0 * i, 2.0 * i + 1};

  auto wDot = DistortionMetric<double>::compute_w_dot(pts, ptsDot);

  // finite difference
  double eps = 1e-7;
  for (int i = 0; i < 2; ++i)
  {
    auto ptsArr = convert_pts(pts);
    for (int j = 0; j < 6; ++j)
      ptsArr[j] += ptsDot[j][i] * eps;

    auto wp = DistortionMetric<double>::compute_w(convert_pts(ptsArr));

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        EXPECT_NEAR((wp(j, k) - w(j, k)) / eps, wDot(j, k)[i], 1e-6);

    // pts_arr gets re-created every iteration, no need to reset it
  }
}

TEST(DistortionMetric, ComputeW_rev_dot)
{
  TriPts<double> pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};

  utils::impl::Mat2x2<DistortionMetric<double>::DotVec<double>> wBarDot;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
    {
      int offset    = 2 * i + j;
      wBarDot(i, j) = DistortionMetric<double>::DotVec<double>{2.0 * offset + 1, 2.0 * offset + 1 + 1};
    }

  DistortionMetric<double>::PtsDot<double> ptsBarDot;
  DistortionMetric<double>::compute_w_rev_dot(pts, wBarDot, ptsBarDot);

  // finite difference
  TriPts<double> ptsBar;
  utils::impl::Mat2x2<double> wBar = {1, 2, 3, 4}; // this value is arbitrary because
                                                   // compute_w_rev_dot does not depend
                                                   // on it
  DistortionMetric<double>::compute_w_rev(pts, ptsBar, wBar);
  double eps = 1e-7;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        wBar(j, k) += wBarDot(j, k)[i] * eps;

    TriPts<double> ptspBar;
    DistortionMetric<double>::compute_w_rev(pts, ptspBar, wBar);

    for (int j = 0; j < 3; ++j)
    {
      EXPECT_NEAR((ptspBar[j].x - ptsBar[j].x) / eps, ptsBarDot[2 * j][i], 1e-6);
      EXPECT_NEAR((ptspBar[j].y - ptsBar[j].y) / eps, ptsBarDot[2 * j + 1][i], 1e-6);
    }

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        wBar(j, k) -= wBarDot(j, k)[i] * eps;
  }
}

void compute_distortion_fd(DistortionMetric<double>& metric, TriPts<double>& pts, TriPts<double>& ptsBar, TriPts<double>& ptsBarFd)
{
  double qBar = 0.816497, eps = 1e-7;
  utils::impl::Mat2x2<double> w{0.5, 0, 0, 0.5};
  metric.get_deriv(pts, w, ptsBar, qBar);

  auto qVal = metric.get_value(pts, w);
  for (int i = 0; i < 3; ++i)
  {
    pts[i].x += eps;
    ptsBarFd[i].x = qBar * (metric.get_value(pts, w) - qVal) / eps;
    pts[i].x -= eps;

    pts[i].y += eps;
    ptsBarFd[i].y = qBar * (metric.get_value(pts, w) - qVal) / eps;
    pts[i].y -= eps;
  }
}

TEST(RegularizedDistortionMetric, Gamma0)
{
  RegularizedDistortionMetric<double> metric(0);
  TriPts<double> ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric<double>::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test scale invariance
  {
    TriPts<double> pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    EXPECT_FLOAT_EQ(metric.get_value(pts, w), 1.0);
  }

  // test that value increases for lower quality elements
  {
    TriPts<double> pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3))};
    EXPECT_GE(metric.get_value(pts, w), 1.0);
  }

  // test invalid element
  {
    TriPts<double> pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, -std::sqrt(3))};
    EXPECT_FALSE(std::isfinite(metric.get_value(pts, w)));
  }

  // test derivative
  {
    TriPts<double> pts{utils::Point(0.6, 0.6), utils::Point(1, 0.5), utils::Point(0.5, 1)};
    TriPts<double> ptsBar, ptsBarFd;
    compute_distortion_fd(metric, pts, ptsBar, ptsBarFd);
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(ptsBar[i].x, ptsBarFd[i].x, 5e-6);
      EXPECT_NEAR(ptsBar[i].y, ptsBarFd[i].y, 5e-6);
    }
  }
}

TEST(RegularizedDistortionMetric, Gamma01)
{
  RegularizedDistortionMetric<double> metric(0.01);
  TriPts<double> ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric<double>::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test derivative
  {
    TriPts<double> pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    TriPts<double> ptsBar, ptsBarFd;
    compute_distortion_fd(metric, pts, ptsBar, ptsBarFd);
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(ptsBar[i].x, ptsBarFd[i].x, 1e-6);
      EXPECT_NEAR(ptsBar[i].y, ptsBarFd[i].y, 1e-6);
    }
  }

  // test derivative on invalid element
  {
    TriPts<double> pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    pts[1].y = -pts[1].y;
    TriPts<double> ptsBar, ptsBarFd;
    compute_distortion_fd(metric, pts, ptsBar, ptsBarFd);
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(ptsBar[i].x, ptsBarFd[i].x, 1e-6);
      EXPECT_NEAR(ptsBar[i].y, ptsBarFd[i].y, 1e-6);
    }
  }
}

TEST(RegularizedDistortionMetric, HessianProduct)
{
  using PtsDot = DistortionMetric<double>::PtsDot<double>;
  using DotVec = DistortionMetric<double>::DotVec<double>;

  RegularizedDistortionMetric<double> metric(0.2);
  TriPts<double> ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric<double>::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test derivative
  {
    TriPts<double> pts{ptsIdeal[0], ptsIdeal[1], ptsIdeal[2]};
    // make element invalid
    pts[2].y = -pts[2].y;

    double minDenominator = metric.compute_denominator(pts, w);
    metric.set_min_denominator(minDenominator);

    PtsDot ptsDot, ptsBarDot;
    double qBar    = 2;
    DotVec qBarDot = {13.0, 14.0};
    for (int i = 0; i < 6; ++i)
      ptsDot[i] = DistortionMetric<double>::DotVec<double>{2.0 * i + 1, 2.0 * i + 2};

    metric.get_value_rev_dot(pts, w, ptsDot, qBar, qBarDot, ptsBarDot);

    // finite difference get_deriv
    TriPts<double> derivs0;
    double eps = 1e-7;
    metric.get_deriv(pts, w, derivs0, qBar);
    auto derivs0Arr = convert_pts(derivs0);
    for (int i = 0; i < 2; ++i)
    {
      auto ptsArr = convert_pts(pts);
      for (int j = 0; j < 6; ++j)
        ptsArr[j] += eps * ptsDot[j][i];
      qBar += eps * qBarDot[i];

      TriPts<double> derivs;
      metric.get_deriv(convert_pts(ptsArr), w, derivs, qBar);

      auto derivsArr = convert_pts(derivs);
      for (int j = 0; j < 6; ++j)
      {
        double valFd = (derivsArr[j] - derivs0Arr[j]) / eps;        
        double scale = std::min(std::abs(valFd), std::abs(ptsBarDot[j][i]));
        EXPECT_NEAR(valFd, ptsBarDot[j][i], 5e-6 * scale);
      }

      qBar -= eps * qBarDot[i];
    }
  }
}


TEST(RegularizedDistortionMetric, Regression)
{
  TriPts<double> pts = { utils::Point(0.875, 0.875, 0), utils::Point(1, 0.75, 0), utils::Point(0.75, 1, 0) };
  utils::impl::Mat2x2<double> W = {0.25, 0, 0, 0.25};
  double gamma = 1e-1;
  int localidx = 0;

  RegularizedDistortionMetric<double> metric(gamma);
  metric.set_min_denominator(metric.compute_denominator(pts, W));

  TriPts<double> derivs;
  metric.get_deriv(pts, W, derivs, 1);

  EXPECT_GT(std::abs(derivs[localidx][0]), 0.01);
  EXPECT_GT(std::abs(derivs[localidx][1]), 0.01);
}

TEST(RegularizedDistortionMetric, Regression2)
{
  using Complex = std::complex<double>;

  TriPts<double> pts = { utils::Point(1.456521739130435, 0, 0), 
                         utils::Point(1.288338373852929, 0.2913272354863573,       0),
                         utils::Point(1.402509962006577, 0.3929648623374049, 0)};

  utils::impl::Mat2x2<double> w({-0.012145808108736, -0.05401177712385752, 0.4046951567355364, 0.3929648623374049});

  double gamma = 1e-3;
  int localidx = 0;

  RegularizedDistortionMetric<double> metric(gamma);
  metric.set_min_denominator(metric.compute_denominator(pts, w));
  RegularizedDistortionMetric<Complex> metricC(metric);


  TriPts<double> derivs;
  metric.get_deriv(pts, w, derivs, 1);

  utils::Point grad = derivs[localidx];

  // complex step
  TriPts<Complex> ptsC;
  ptsC[0] = pts[0];
  ptsC[1] = pts[1];
  ptsC[2] = pts[2];

  double h = 1e-40;
  Complex pert(0, h);
  ptsC[localidx].x += pert;
  Complex qx = metricC.get_value(ptsC, w);
  ptsC[localidx].x -= pert;

  ptsC[localidx].y += pert;
  Complex qy = metricC.get_value(ptsC, w);
  ptsC[localidx].y -= pert;

  utils::Point gradC(qx.imag()/h, qy.imag()/h);

  for (int i=0; i < 2; ++i)
  {
    double minMagnitude = std::min(std::abs(grad[i]), std::abs(gradC[i]));
    double errRelative = std::abs(grad[i] - gradC[i])/minMagnitude;
    double errAbs = std::abs(grad[i] - gradC[i]);
    double err = std::min(errAbs, errRelative); 
    EXPECT_LT(err, 1e-12);
  } 
}

TEST(RegularizedDistortionMetric, Regression3)
{
  using Complex = std::complex<double>;
  using DotVec = DistortionMetric<double>::DotVec<double>;
  using PtsDot = DistortionMetric<double>::PtsDot<double>;

  TriPts<double> pts = { utils::Point(1.456521739130435, 0, 0), 
                         utils::Point(1.288338373852929, 0.2913272354863573, 0      ), 
                         utils::Point(1.402509962006577, 0.3929648623374049, 0)};

  utils::impl::Mat2x2<double> w({-0.012145808108736, -0.05401177712385752, 0.4046951567355364, 0.3929648623374049});
  double qIBar = 0.2041241452319302;
  DotVec qIBarDot = {6.661338147750939e-15, -0.5597855687575399};

  double gamma = 1e-3;
  int localidx = 0;

  RegularizedDistortionMetric<double> metric(gamma);
  metric.set_min_denominator(metric.compute_denominator(pts, w));
  RegularizedDistortionMetric<Complex> metricC(metric);


  PtsDot ptsBarDot;
  PtsDot ptsIDot;
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 2; ++k)
      ptsIDot[j][k] = 0;

  ptsIDot[2 * localidx][0]     = 1;
  ptsIDot[2 * localidx + 1][1] = 1;

  metric.get_value_rev_dot(pts, w, ptsIDot, qIBar, qIBarDot, ptsBarDot);

  utils::impl::Mat2x2<double> hess;
  hess(0, 0) = ptsBarDot[2 * localidx][0];
  hess(0, 1) = ptsBarDot[2 * localidx][1];
  hess(1, 0) = ptsBarDot[2 * localidx + 1][0];
  hess(1, 1) = ptsBarDot[2 * localidx + 1][1];


  // complex step
  TriPts<Complex> ptsC;
  ptsC[0] = pts[0];
  ptsC[1] = pts[1];
  ptsC[2] = pts[2];

  double h = 1e-40;
  Complex pert(0, h);

  Complex qBarx(qIBar, h*qIBarDot[0]);
  TriPts<Complex> derivsx;
  ptsC[localidx].x += pert;
  metricC.get_deriv(ptsC, w, derivsx, qBarx);
  ptsC[localidx].x -= pert;

  ptsC[localidx].y += pert;
  Complex qBary(qIBar, h*qIBarDot[1]);
  TriPts<Complex> derivsy;
  metricC.get_deriv(ptsC, w, derivsy, qBary);
  ptsC[localidx].y -= pert;

  utils::impl::Mat2x2<double> hessC;
  hessC(0, 0) = derivsx[localidx].x.imag()/h;
  hessC(1, 0) = derivsx[localidx].y.imag()/h;
  hessC(0, 1) = derivsy[localidx].x.imag()/h;
  hessC(1, 1) = derivsy[localidx].y.imag()/h;  

  for (int i=0; i < 2; ++i)
    for (int j=0; j < 2; ++j)
    {
      double minMagnitude = std::min(std::abs(hess(i, j)), std::abs(hessC(i, j)));
      double errRelative = std::abs(hess(i, j) - hessC(i, j))/minMagnitude;
      double errAbs = std::abs(hess(i, j) - hessC(i, j));
      double err = std::min(errAbs, errRelative); 
      EXPECT_LT(err, 1e-12);
    } 
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
