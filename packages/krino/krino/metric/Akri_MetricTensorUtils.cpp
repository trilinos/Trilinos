// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MetricTensorUtils.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

namespace krino {
namespace metric {

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

std::array<double,9> to_array(const minitensor::Tensor<double, 3> & tensor)
{
  std::array<double,9> arr;
  copy_tensor_to_data(tensor, arr.data());
  return arr;
}

std::array<double,4> to_array(const minitensor::Tensor<double, 2> & tensor)
{
  std::array<double,4> arr;
  copy_tensor_to_data(tensor, arr.data());
  return arr;
}

// WORKAROUND for minitensor::eig_sym: its 2x2 specialization returns eigenvectors that do not
// diagonalize the tensor, so 2D tensors use this closed form instead.  The eigenvectors of a symmetric
// 2x2 tensor A are the columns of the rotation V by angle theta, with tan(2 theta) = 2 A01 / (A00 - A11).
// When a version of minitensor with a correct eig_sym is available, remove this function and call
// minitensor::eig_sym directly in symmetric_eigen_decomposition below.
static std::pair<Tensor2d, Tensor2d> eig_sym_2x2_workaround(const Tensor2d & A)
{
  const double offDiag = 0.5 * (A(0,1) + A(1,0));
  const double theta = 0.5 * std::atan2(2.0 * offDiag, A(0,0) - A(1,1));
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const double eigenvalue0 = c*c*A(0,0) + 2.0*c*s*offDiag + s*s*A(1,1);
  const double eigenvalue1 = s*s*A(0,0) - 2.0*c*s*offDiag + c*c*A(1,1);
  return {Tensor2d(c, -s, s, c), Tensor2d(eigenvalue0, 0., 0., eigenvalue1)};
}

template<minitensor::Index N>
std::pair<minitensor::Tensor<double, N>, minitensor::Tensor<double, N>> symmetric_eigen_decomposition(const minitensor::Tensor<double, N> & A)
{
  if constexpr (N == 2)
    return eig_sym_2x2_workaround(A);
  else
    return minitensor::eig_sym(A);
}

template<minitensor::Index N>
static minitensor::Tensor<double, N> compute_metric_mapping_of_dimension(const minitensor::Tensor<double, N> & metric)
{
  // eig_sym returns (eigenvectors V, eigenvalues D) with metric = V D V^T, so Q = sqrt(D) V^T gives Q^T Q = metric.
  // Flip an eigenvector if needed so that det(Q) > 0; otherwise Q would reflect the ideal element.
  const auto & [eigenvectors, eigenvalues] = symmetric_eigen_decomposition(metric);
  minitensor::Tensor<double, N> metricMapping = minitensor::sqrt(eigenvalues) * minitensor::transpose(eigenvectors);
  if (minitensor::det(metricMapping) < 0.)
    for (minitensor::Index j = 0; j < N; ++j)
      metricMapping(0, j) = -metricMapping(0, j);
  return metricMapping;
}

Tensor compute_metric_mapping(
    const Tensor & metric)
{
  return compute_metric_mapping_of_dimension(metric);
}

template<minitensor::Index N>
minitensor::Tensor<double, N> apply_to_symmetric_tensor_eigenvalues(const minitensor::Tensor<double, N> & A, const std::function<double(double)> & fn)
{
  const auto & [V, D] = symmetric_eigen_decomposition(A);
  minitensor::Tensor<double, N> fnD(minitensor::Filler::ZEROS);
  for (minitensor::Index i=0; i<N; ++i)
    fnD(i,i) = fn(D(i,i));
  return V * fnD * minitensor::transpose(V);
}

template<minitensor::Index N>
minitensor::Tensor<double, N> log_euclidean_mean_of_metrics(const std::vector<const double *> & metrics)
{
  STK_ThrowRequire(!metrics.empty());
  minitensor::Tensor<double, N> sumLogMetric(minitensor::Filler::ZEROS);
  for (const double * metric : metrics)
    sumLogMetric += apply_to_symmetric_tensor_eigenvalues(minitensor::Tensor<double, N>(metric), [](const double x) { return std::log(x); });
  return apply_to_symmetric_tensor_eigenvalues<N>((1./metrics.size()) * sumLogMetric, [](const double x) { return std::exp(x); });
}

template<minitensor::Index N>
static std::vector<double> compute_log_euclidean_mean_of_metrics_of_dimension(const std::vector<const double *> & metrics)
{
  std::vector<double> meanData(N*N);
  copy_tensor_to_data(log_euclidean_mean_of_metrics<N>(metrics), meanData.data());
  return meanData;
}

std::vector<double> compute_log_euclidean_mean_of_metrics(const unsigned spatialDim, const std::vector<const double *> & metrics)
{
  if (spatialDim == 2)
    return compute_log_euclidean_mean_of_metrics_of_dimension<2>(metrics);
  return compute_log_euclidean_mean_of_metrics_of_dimension<3>(metrics);
}

template<minitensor::Index N>
static std::vector<double> compute_metric_mapping_data_from_nodal_metrics_of_dimension(const std::vector<const double *> & nodalMetrics)
{
  const minitensor::Tensor<double, N> meanMetric = log_euclidean_mean_of_metrics<N>(nodalMetrics);

  std::vector<double> metricMappingData(N*N);
  copy_tensor_to_data(compute_metric_mapping_of_dimension(meanMetric), metricMappingData.data());
  return metricMappingData;
}

std::vector<double> compute_metric_mapping_data_from_nodal_metrics(
    const unsigned spatialDim,
    const std::vector<const double *> & nodalMetrics)
{
  if (spatialDim == 2)
    return compute_metric_mapping_data_from_nodal_metrics_of_dimension<2>(nodalMetrics);
  return compute_metric_mapping_data_from_nodal_metrics_of_dimension<3>(nodalMetrics);
}

template<minitensor::Index N>
static bool is_symmetric_positive_definite_metric_of_dimension(const double * metric)
{
  double maxAbsComponent = 0.;
  for (minitensor::Index i=0; i<N*N; ++i)
    maxAbsComponent = std::max(maxAbsComponent, std::abs(metric[i]));
  const double symmetryTol = 1.e-12*maxAbsComponent;
  for (minitensor::Index i=0; i<N; ++i)
    for (minitensor::Index j=i+1; j<N; ++j)
      if (!(std::abs(metric[i*N+j] - metric[j*N+i]) <= symmetryTol))
        return false;

  const auto & [V, D] = symmetric_eigen_decomposition(minitensor::Tensor<double, N>(metric));
  for (minitensor::Index i=0; i<N; ++i)
    if (!(D(i,i) > 0.))
      return false;
  return true;
}

bool is_symmetric_positive_definite_metric(
    const unsigned spatialDim,
    const double * metric)
{
  if (spatialDim == 2)
    return is_symmetric_positive_definite_metric_of_dimension<2>(metric);
  return is_symmetric_positive_definite_metric_of_dimension<3>(metric);
}

Tensor2d compute_metric_mapping_2d(
    const Tensor2d & metric)
{
  return compute_metric_mapping_of_dimension(metric);
}

Tensor metric_tensor_for_size(
    const double h)
{
  Tensor metric(minitensor::Filler::ZEROS);
  const double c = 1./(h*h);
  metric(0,0) = c;
  metric(1,1) = c;
  metric(2,2) = c;
  return metric;
}

Tensor2d metric_tensor_for_size_2d(
    const double h)
{
  Tensor2d metric(minitensor::Filler::ZEROS);
  const double c = 1./(h*h);
  metric(0,0) = c;
  metric(1,1) = c;
  return metric;
}

double compute_relative_length_in_metric(
    const unsigned spatialDim,
    const double * metric,
    const stk::math::Vector3d & edge)
{
  double lengthSquared = 0.;
  for (unsigned i=0; i<spatialDim; ++i)
    for (unsigned j=0; j<spatialDim; ++j)
      lengthSquared += edge[i] * metric[i*spatialDim+j] * edge[j];
  return std::sqrt(lengthSquared);
}

double logarithmic_mean(const double a, const double b)
{
  // Same threshold as ngs, so the two agree on which edges are too long.
  constexpr double tol = 1.e-3;
  if (std::abs(a-b) < tol)
    return 0.5*(a+b);
  return (a-b) / std::log(a/b);
}

double compute_relative_edge_length(
    const unsigned spatialDim,
    const stk::math::Vector3d & x0,
    const double * metric0,
    const stk::math::Vector3d & x1,
    const double * metric1)
{
  const stk::math::Vector3d edge = x1 - x0;
  return logarithmic_mean(compute_relative_length_in_metric(spatialDim, metric0, edge),
      compute_relative_length_in_metric(spatialDim, metric1, edge));
}

template std::pair<Tensor2d, Tensor2d> symmetric_eigen_decomposition<2>(const Tensor2d &);
template std::pair<Tensor, Tensor> symmetric_eigen_decomposition<3>(const Tensor &);
template Tensor2d apply_to_symmetric_tensor_eigenvalues<2>(const Tensor2d &, const std::function<double(double)> &);
template Tensor apply_to_symmetric_tensor_eigenvalues<3>(const Tensor &, const std::function<double(double)> &);
template Tensor2d log_euclidean_mean_of_metrics<2>(const std::vector<const double *> &);
template Tensor log_euclidean_mean_of_metrics<3>(const std::vector<const double *> &);

} // namespace metric
} // namespace krino
