// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MetricTensorUtils.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace krino {
namespace metric {

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

constexpr double tol = 1.e-6;

template <unsigned DIM>
static void expect_tensor(const minitensor::Tensor<double, DIM> & gold, const minitensor::Tensor<double, DIM> & A)
{
  for (unsigned i=0; i<DIM; ++i)
    for (unsigned j=0; j<DIM; ++j)
      EXPECT_NEAR(gold(i,j), A(i,j), tol);
}

// A rotated anisotropic metric with eigenvalues 3, 0.75 (and 2 in 3D).
static Tensor2d rotated_metric_2d(const double theta)
{
  const Tensor2d rotation(std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta));
  return rotation * Tensor2d(3.0, 0.0, 0.0, 0.75) * minitensor::transpose(rotation);
}

static Tensor rotated_metric(const double theta)
{
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const Tensor R = Tensor(c, -s, 0., s, c, 0., 0., 0., 1.) * Tensor(1., 0., 0., 0., c, -s, 0., s, c);
  return R * Tensor(3.0, 0., 0., 0., 0.75, 0., 0., 0., 2.0) * minitensor::transpose(R);
}

TEST(MetricTensorUtils, metricMapping2d_factorsMetricsOfAnyOrientation)
{
  // Covers both orderings of the diagonal entries and a diagonal metric.
  for (const double theta : {0.0, 0.3, 1.2, 1.5707963267948966, 2.0, 2.8})
  {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Tensor2d rotation(c, -s, s, c);
    const Tensor2d diagonal(3.0, 0.0, 0.0, 0.75);
    const Tensor2d metric = rotation * diagonal * minitensor::transpose(rotation);
    const Tensor2d metricMapping = metric::compute_metric_mapping_2d(metric);

    expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
    EXPECT_NEAR(1.5, minitensor::det(metricMapping), tol);
  }
}

TEST(MetricTensorUtils, metricMapping_isProperForMetricsOfAnyOrientation)
{
  // A metric mapping with det(F_M) < 0 would reflect the ideal element, so that elements of positive
  // orientation would appear inverted.
  for (const double theta : {0.0, 0.3, 1.2, 1.5707963267948966, 2.0, 2.8})
  {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Tensor rotation(c, -s, 0., s, c, 0., 0., 0., 1.);
    const Tensor tilt(1., 0., 0., 0., c, -s, 0., s, c);
    const Tensor diagonal(3.0, 0., 0., 0., 0.75, 0., 0., 0., 2.0);
    const Tensor R = rotation * tilt;
    const Tensor metric = R * diagonal * minitensor::transpose(R);
    const Tensor metricMapping = metric::compute_metric_mapping(metric);

    expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
    EXPECT_NEAR(std::sqrt(4.5), minitensor::det(metricMapping), tol);
  }
}

TEST(MetricTensorUtils, symmetricEigenDecomposition_reconstructsTensorOfAnyOrientation)
{
  for (const double theta : {0.0, 0.3, 1.2, 1.5707963267948966, 2.0, 2.8})
  {
    const Tensor2d metric2d = rotated_metric_2d(theta);
    const auto & [V2, D2] = symmetric_eigen_decomposition(metric2d);
    expect_tensor(metric2d, Tensor2d(V2 * D2 * minitensor::transpose(V2)));
    expect_tensor(minitensor::identity<double, 2>(), Tensor2d(minitensor::transpose(V2) * V2));

    const Tensor metric = rotated_metric(theta);
    const auto & [V, D] = symmetric_eigen_decomposition(metric);
    expect_tensor(metric, Tensor(V * D * minitensor::transpose(V)));
    expect_tensor(minitensor::identity<double, 3>(), Tensor(minitensor::transpose(V) * V));
  }
}

TEST(MetricTensorUtils, logThenExp_givesBackTheTensor)
{
  const auto log = [](const double x) { return std::log(x); };
  const auto exp = [](const double x) { return std::exp(x); };

  const Tensor2d metric2d = rotated_metric_2d(0.7);
  expect_tensor(metric2d, apply_to_symmetric_tensor_eigenvalues<2>(apply_to_symmetric_tensor_eigenvalues<2>(metric2d, log), exp));

  const Tensor metric = rotated_metric(0.7);
  expect_tensor(metric, apply_to_symmetric_tensor_eigenvalues<3>(apply_to_symmetric_tensor_eigenvalues<3>(metric, log), exp));
}

TEST(MetricTensorUtils, logEuclideanMean_ofMetricAndItsInverseIsIdentity)
{
  const auto metric2d = to_array(rotated_metric_2d(0.4));
  const auto inverse2d = to_array(Tensor2d(minitensor::inverse(rotated_metric_2d(0.4))));
  expect_tensor(minitensor::identity<double, 2>(), log_euclidean_mean_of_metrics<2>({metric2d.data(), inverse2d.data()}));

  const auto metric = to_array(rotated_metric(0.4));
  const auto inverse = to_array(Tensor(minitensor::inverse(rotated_metric(0.4))));
  expect_tensor(minitensor::identity<double, 3>(), log_euclidean_mean_of_metrics<3>({metric.data(), inverse.data()}));
}

TEST(MetricTensorUtils, logEuclideanMean_ofIsotropicMetricsIsMetricForGeometricMeanSize)
{
  const auto small = to_array(metric_tensor_for_size(1.));
  const auto large = to_array(metric_tensor_for_size(4.));
  expect_tensor(metric_tensor_for_size(2.), log_euclidean_mean_of_metrics<3>({small.data(), large.data()}));
}

TEST(MetricTensorUtils, computeLogEuclideanMean_rowMajorDataMatchesTensorVersion)
{
  const auto metric2d = to_array(rotated_metric_2d(0.4));
  const auto other2d = to_array(Tensor2d(2.*rotated_metric_2d(1.3)));
  const std::vector<double> mean2d = compute_log_euclidean_mean_of_metrics(2, {metric2d.data(), other2d.data()});
  const auto gold2d = to_array(log_euclidean_mean_of_metrics<2>({metric2d.data(), other2d.data()}));
  ASSERT_EQ(4u, mean2d.size());
  for (unsigned i=0; i<4; ++i)
    EXPECT_NEAR(gold2d[i], mean2d[i], tol);

  const auto metric = to_array(rotated_metric(0.4));
  const auto other = to_array(Tensor(2.*rotated_metric(1.3)));
  const std::vector<double> mean = compute_log_euclidean_mean_of_metrics(3, {metric.data(), other.data()});
  const auto gold = to_array(log_euclidean_mean_of_metrics<3>({metric.data(), other.data()}));
  ASSERT_EQ(9u, mean.size());
  for (unsigned i=0; i<9; ++i)
    EXPECT_NEAR(gold[i], mean[i], tol);
}

TEST(MetricTensorUtils, isSymmetricPositiveDefiniteMetric)
{
  const auto metric = to_array(rotated_metric(0.9));
  EXPECT_TRUE(is_symmetric_positive_definite_metric(3, metric.data()));
  EXPECT_TRUE(is_symmetric_positive_definite_metric(2, to_array(rotated_metric_2d(0.9)).data()));

  auto roundoff = metric;
  roundoff[1] *= 1. + 1.e-14; // M_01
  EXPECT_TRUE(is_symmetric_positive_definite_metric(3, roundoff.data()));

  auto nonSymmetric = metric;
  nonSymmetric[1] += 0.5; // M_01, with M_10 unchanged
  EXPECT_FALSE(is_symmetric_positive_definite_metric(3, nonSymmetric.data()));

  const std::array<double,4> zeroEigenvalue{{1., 0., 0., 0.}};
  EXPECT_FALSE(is_symmetric_positive_definite_metric(2, zeroEigenvalue.data()));
  const std::array<double,9> negative{{-1., 0., 0., 0., -1., 0., 0., 0., -1.}};
  EXPECT_FALSE(is_symmetric_positive_definite_metric(3, negative.data()));
}

// Diagonal metric with sizes hx, hy (and hz in 3D), row-major.
static std::vector<double> diagonal_metric(const std::vector<double> & sizes)
{
  const size_t dim = sizes.size();
  std::vector<double> metric(dim*dim, 0.);
  for (size_t i=0; i<dim; ++i)
    metric[i*dim+i] = 1./(sizes[i]*sizes[i]);
  return metric;
}

TEST(MetricTensorUtils, relativeLengthInMetric_isotropicIsLengthOverSize)
{
  const std::vector<double> metric = diagonal_metric({0.5, 0.5, 0.5});
  const stk::math::Vector3d edge(1., 2., 2.); // length 3
  EXPECT_NEAR(6., compute_relative_length_in_metric(3, metric.data(), edge), tol);
}

TEST(MetricTensorUtils, relativeLengthInMetric_anisotropicUsesSizeInEdgeDirection)
{
  const std::vector<double> metric = diagonal_metric({0.1, 2.0, 1.0});
  EXPECT_NEAR(30., compute_relative_length_in_metric(3, metric.data(), stk::math::Vector3d(3., 0., 0.)), tol);
  EXPECT_NEAR(1.5, compute_relative_length_in_metric(3, metric.data(), stk::math::Vector3d(0., 3., 0.)), tol);

  const std::vector<double> metric2d = diagonal_metric({0.1, 2.0});
  // In 2D the z component of the edge is ignored.
  EXPECT_NEAR(30., compute_relative_length_in_metric(2, metric2d.data(), stk::math::Vector3d(3., 0., 7.)), tol);
  EXPECT_NEAR(1.5, compute_relative_length_in_metric(2, metric2d.data(), stk::math::Vector3d(0., 3., 7.)), tol);
}

TEST(MetricTensorUtils, relativeLengthInMetric_rotatedMetricGivesSameLengthForRotatedEdge)
{
  // rotated_metric_2d has eigenvalues 3 and 0.75 along the rotated axes.
  const double theta = 0.6;
  const auto metric = to_array(rotated_metric_2d(theta));
  const stk::math::Vector3d alongFirstAxis(std::cos(theta), std::sin(theta), 0.);
  const stk::math::Vector3d alongSecondAxis(-std::sin(theta), std::cos(theta), 0.);
  EXPECT_NEAR(std::sqrt(3.), compute_relative_length_in_metric(2, metric.data(), alongFirstAxis), tol);
  EXPECT_NEAR(std::sqrt(0.75), compute_relative_length_in_metric(2, metric.data(), alongSecondAxis), tol);
}

TEST(MetricTensorUtils, logarithmicMean)
{
  EXPECT_NEAR(1./std::log(2.), logarithmic_mean(2., 1.), tol);
  EXPECT_DOUBLE_EQ(logarithmic_mean(2., 1.), logarithmic_mean(1., 2.));
  EXPECT_DOUBLE_EQ(1.7, logarithmic_mean(1.7, 1.7));
  EXPECT_DOUBLE_EQ(0., logarithmic_mean(0., 0.));

  // Between the geometric and arithmetic means.
  const double a = 3.;
  const double b = 0.5;
  EXPECT_LT(std::sqrt(a*b), logarithmic_mean(a, b));
  EXPECT_GT(0.5*(a+b), logarithmic_mean(a, b));

  // Nearly continuous across the switch to the arithmetic mean: near the switch, (a - b)/ln(a/b) differs
  // from (a + b)/2 by about (a - b)^2/(12 a), here ~1e-7.
  const double justOutside = 1. + 1.001e-3;
  EXPECT_NEAR(0.5*(1.+justOutside), logarithmic_mean(1., justOutside), 1.e-7);
}

TEST(MetricTensorUtils, relativeEdgeLength_sameMetricAtBothNodesIsLengthInThatMetric)
{
  const auto metric = to_array(rotated_metric(0.8));
  const stk::math::Vector3d x0(0.1, -0.2, 0.3);
  const stk::math::Vector3d x1(1.4, 0.5, -0.6);
  EXPECT_NEAR(compute_relative_length_in_metric(3, metric.data(), x1-x0),
      compute_relative_edge_length(3, x0, metric.data(), x1, metric.data()), tol);
}

TEST(MetricTensorUtils, relativeEdgeLength_exactForGeometricSizeVariation)
{
  // Isotropic sizes h0 at x0 and h1 at x1: with h(t) = h0^(1-t) h1^t along the edge, the relative length is
  // the integral of |e|/h(t) over t in [0,1], which is |e| (1/h0 - 1/h1)/ln(h1/h0).
  const double h0 = 0.1;
  const double h1 = 0.8;
  const stk::math::Vector3d x0(0., 0., 0.);
  const stk::math::Vector3d x1(0.3, 0.4, 0.); // length 0.5
  const std::vector<double> metric0 = diagonal_metric({h0, h0, h0});
  const std::vector<double> metric1 = diagonal_metric({h1, h1, h1});
  const double length = compute_relative_edge_length(3, x0, metric0.data(), x1, metric1.data());

  EXPECT_NEAR(0.5*(1./h0 - 1./h1)/std::log(h1/h0), length, tol);

  const int numSteps = 100000;
  double integral = 0.;
  for (int i=0; i<numSteps; ++i)
  {
    const double t = (i+0.5)/numSteps;
    integral += 0.5/(std::pow(h0, 1.-t)*std::pow(h1, t))/numSteps;
  }
  EXPECT_NEAR(integral, length, 1.e-6);
}

TEST(MetricTensorUtils, relativeEdgeLength_symmetricTranslationInvariantAndLinearInEdge)
{
  const auto metric0 = to_array(rotated_metric(0.3));
  const auto metric1 = to_array(Tensor(4.*rotated_metric(1.1)));
  const stk::math::Vector3d x0(0.2, 0.1, -0.4);
  const stk::math::Vector3d x1(0.9, -0.3, 0.5);
  const double length = compute_relative_edge_length(3, x0, metric0.data(), x1, metric1.data());

  EXPECT_NEAR(length, compute_relative_edge_length(3, x1, metric1.data(), x0, metric0.data()), tol);

  const stk::math::Vector3d shift(5., -3., 2.);
  EXPECT_NEAR(length, compute_relative_edge_length(3, x0+shift, metric0.data(), x1+shift, metric1.data()), tol);

  EXPECT_NEAR(2.*length, compute_relative_edge_length(3, x0, metric0.data(), x0+2.*(x1-x0), metric1.data()), tol);
}

TEST(MetricTensorUtils, relativeEdgeLength_2d)
{
  const std::vector<double> metric0 = diagonal_metric({0.2, 1.0});
  const std::vector<double> metric1 = diagonal_metric({0.4, 1.0});
  const stk::math::Vector3d x0(0., 0., 0.);
  const stk::math::Vector3d x1(1., 0., 0.);
  // Lengths 5 and 2.5 along x.
  EXPECT_NEAR(2.5/std::log(2.), compute_relative_edge_length(2, x0, metric0.data(), x1, metric1.data()), tol);
}

} // namespace metric
} // namespace krino
