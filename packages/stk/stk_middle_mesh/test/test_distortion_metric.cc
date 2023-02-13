#include "gtest/gtest.h"

#include "distortion_metric.h"
#include "regularized_distortion_metric.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh::impl;

TriPts convert_pts(std::array<double, 6> pts)
{
  return {utils::Point(pts[0], pts[1]), utils::Point(pts[2], pts[3]), utils::Point(pts[4], pts[5])};
}

std::array<double, 6> convert_pts(const TriPts& pts)
{
  return {pts[0].x, pts[0].y, pts[1].x, pts[1].y, pts[2].x, pts[2].y};
}

utils::impl::Mat2x2<double> compute_w_fd(std::array<double, 6> pts, int idx)
{
  auto ptsAos = convert_pts(pts);
  auto w      = DistortionMetric::compute_w(ptsAos);

  double eps = 1e-7;
  pts[idx] += eps;
  auto wPert = DistortionMetric::compute_w(convert_pts(pts));
  auto dWDpt = (wPert - w) / eps;
  pts[idx] -= eps;

  return dWDpt;
}

TEST(DistortionMetric, ComputeW)
{
  TriPts pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric::compute_w(pts);

  EXPECT_FLOAT_EQ(w(0, 0), 1.0);
  EXPECT_FLOAT_EQ(w(0, 1), 0.5);
  EXPECT_FLOAT_EQ(w(1, 0), 0.0);
  EXPECT_FLOAT_EQ(w(1, 1), std::sqrt(3) / 2.0);

  TriPts ptsBar;
  utils::impl::Mat2x2<double> wBar{1, 2, 3, 4};
  DistortionMetric::compute_w_rev(pts, ptsBar, wBar);

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
  TriPts pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric::compute_w(pts);

  std::array<std::array<double, 2>, 6> ptsDot;
  for (int i = 0; i < 6; ++i)
    ptsDot[i] = std::array<double, 2>{2.0 * i, 2.0 * i + 1};

  auto wDot = DistortionMetric::compute_w_dot(pts, ptsDot);

  // finite difference
  double eps = 1e-7;
  for (int i = 0; i < 2; ++i)
  {
    auto ptsArr = convert_pts(pts);
    for (int j = 0; j < 6; ++j)
      ptsArr[j] += ptsDot[j][i] * eps;

    auto wp = DistortionMetric::compute_w(convert_pts(ptsArr));

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        EXPECT_NEAR((wp(j, k) - w(j, k)) / eps, wDot(j, k)[i], 1e-6);

    // pts_arr gets re-created every iteration, no need to reset it
  }
}

TEST(DistortionMetric, ComputeW_rev_dot)
{
  TriPts pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};

  utils::impl::Mat2x2<DistortionMetric::DotVec> wBarDot;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
    {
      int offset    = 2 * i + j;
      wBarDot(i, j) = DistortionMetric::DotVec{2.0 * offset + 1, 2.0 * offset + 1 + 1};
    }

  DistortionMetric::PtsDot ptsBarDot;
  DistortionMetric::compute_w_rev_dot(pts, wBarDot, ptsBarDot);

  // finite difference
  TriPts ptsBar;
  utils::impl::Mat2x2<double> wBar = {1, 2, 3, 4}; // this value is arbitrary because
                                                   // compute_w_rev_dot does not depend
                                                   // on it
  DistortionMetric::compute_w_rev(pts, ptsBar, wBar);
  double eps = 1e-7;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        wBar(j, k) += wBarDot(j, k)[i] * eps;

    TriPts ptspBar;
    DistortionMetric::compute_w_rev(pts, ptspBar, wBar);

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

void compute_distortion_fd(DistortionMetric& metric, TriPts& pts, TriPts& ptsBar, TriPts& ptsBarFd)
{
  double qBar = 0.816497, eps = 1e-7;
  // TriPts pts_ideal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3)/2.0)};
  // auto W = DistortionMetric::compute_w(pts_ideal);
  utils::impl::Mat2x2<double> w{0.5, 0, 0, 0.5};
  metric.get_deriv(pts, w, ptsBar, qBar);

  auto qVal = metric.get_value(pts, w);
  std::cout << "q_var = " << qVal << std::endl;
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

TEST(RegularizedDistortionMetric, Delta0)
{
  RegularizedDistortionMetric metric(0);
  TriPts ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test scale invariance
  {
    TriPts pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    EXPECT_FLOAT_EQ(metric.get_value(pts, w), 1.0);
  }

  // test that value increases for lower quality elements
  {
    TriPts pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3))};
    EXPECT_GE(metric.get_value(pts, w), 1.0);
  }

  // test invalid element
  {
    TriPts pts{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, -std::sqrt(3))};
    EXPECT_FALSE(std::isfinite(metric.get_value(pts, w)));
  }

  // test derivative
  {
    // TriPts pts{pts_ideal[0]*2, pts_ideal[1]*2 + utils::Point(2, 0), pts_ideal[2]*2};
    TriPts pts{utils::Point(0.6, 0.6), utils::Point(1, 0.5), utils::Point(0.5, 1)};
    TriPts ptsBar, ptsBarFd;
    compute_distortion_fd(metric, pts, ptsBar, ptsBarFd);
    for (int i = 0; i < 3; ++i)
    {
      std::cout << "pts_bar_fd = " << ptsBarFd[i] << std::endl;
      std::cout << "pts_bar    = " << ptsBar[i] << std::endl;
      EXPECT_NEAR(ptsBar[i].x, ptsBarFd[i].x, 5e-6);
      EXPECT_NEAR(ptsBar[i].y, ptsBarFd[i].y, 5e-6);
    }
  }
}

TEST(RegularizedDistortionMetric, Delta01)
{
  RegularizedDistortionMetric metric(0.01);
  TriPts ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test derivative
  {
    TriPts pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    TriPts ptsBar, ptsBarFd;
    compute_distortion_fd(metric, pts, ptsBar, ptsBarFd);
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(ptsBar[i].x, ptsBarFd[i].x, 1e-6);
      EXPECT_NEAR(ptsBar[i].y, ptsBarFd[i].y, 1e-6);
    }
  }

  // test derivative on invalid element
  {
    TriPts pts{ptsIdeal[0] * 2, ptsIdeal[1] * 2, ptsIdeal[2] * 2};
    pts[1].y = -pts[1].y;
    TriPts ptsBar, ptsBarFd;
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
  using PtsDot = DistortionMetric::PtsDot;
  using DotVec = DistortionMetric::DotVec;

  RegularizedDistortionMetric metric(1.0);
  TriPts ptsIdeal{utils::Point(0, 0), utils::Point(1, 0), utils::Point(0.5, std::sqrt(3) / 2.0)};
  auto w = DistortionMetric::compute_w(ptsIdeal);

  EXPECT_FLOAT_EQ(metric.get_value(ptsIdeal, w), 1.0);

  // test derivative
  {
    // TriPts pts{pts_ideal[0]*2, pts_ideal[1]*2 + utils::Point(2, 0), pts_ideal[2]*2};
    TriPts pts{ptsIdeal[0], ptsIdeal[1], ptsIdeal[2]};
    // make element invalid
    pts[2].y = -pts[2].y;

    PtsDot ptsDot, ptsBarDot;
    double qBar    = 2;
    DotVec qBarDot = {13.0, 14.0};
    for (int i = 0; i < 6; ++i)
      ptsDot[i] = DistortionMetric::DotVec{2.0 * i + 1, 2.0 * i + 2};

    metric.get_value_rev_dot(pts, w, ptsDot, qBar, qBarDot, ptsBarDot);

    // finite difference get_deriv
    TriPts derivs0;
    double eps = 1e-7;
    metric.get_deriv(pts, w, derivs0, qBar);
    auto derivs0Arr = convert_pts(derivs0);
    for (int i = 0; i < 2; ++i)
    {
      auto ptsArr = convert_pts(pts);
      for (int j = 0; j < 6; ++j)
        ptsArr[j] += eps * ptsDot[j][i];
      qBar += eps * qBarDot[i];

      TriPts derivs;
      metric.get_deriv(convert_pts(ptsArr), w, derivs, qBar);

      auto derivsArr = convert_pts(derivs);
      for (int j = 0; j < 6; ++j)
      {
        std::cout << "val_fd      = " << (derivsArr[j] - derivs0Arr[j]) / eps << std::endl;
        std::cout << "pts_bar_dot = " << ptsBarDot[j][i] << std::endl;

        std::cout << "diff        = " << (derivsArr[j] - derivs0Arr[j]) / eps - ptsBarDot[j][i] << std::endl;
        EXPECT_NEAR((derivsArr[j] - derivs0Arr[j]) / eps, ptsBarDot[j][i], 2e-5);
      }

      qBar -= eps * qBarDot[i];
    }
  }
}

TEST(DistortionMetric, Regression)
{
  auto pt1 = utils::Point(-1.254914017606327, -0.716672555794041);
  auto pt2 = utils::Point(-1.213525491562421, -0.8816778784387096);
  auto pt3 = utils::Point(-0.4635254915624213, -1.42658477444273);
  auto pt4 = utils::Point(-0.590411130252923, -1.398520034038108);

  auto w1   = DistortionMetric::compute_w({pt1, pt2, pt4});
  auto det1 = det2x2(w1);

  auto w2   = DistortionMetric::compute_w({pt2, pt3, pt3});
  auto det2 = det2x2(w2);

  auto w3   = DistortionMetric::compute_w({pt1, pt2, pt3});
  auto det3 = det2x2(w3);

  auto w4   = DistortionMetric::compute_w({pt1, pt3, pt4});
  auto det4 = det2x2(w4);

  std::cout << "det1 = " << det1 << std::endl;
  std::cout << "det2 = " << det2 << std::endl;
  std::cout << "det3 = " << det3 << std::endl;
  std::cout << "det4 = " << det4 << std::endl;
  std::cout << "W2 = " << w2 << std::endl;
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
