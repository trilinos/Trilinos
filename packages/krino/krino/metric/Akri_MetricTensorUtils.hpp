// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_METRIC_AKRI_METRICTENSORUTILS_HPP_
#define KRINO_KRINO_METRIC_AKRI_METRICTENSORUTILS_HPP_

#include <MiniTensor.h>
#include <stk_math/StkVector.hpp>

#include <array>
#include <functional>
#include <utility>
#include <vector>

namespace krino {
namespace metric {

// Metric tensors M (symmetric positive definite), with M = F_M^T F_M for a metric mapping F_M.

template<minitensor::Index N>
void copy_tensor_to_data(const minitensor::Tensor<double, N> & tensor, double * tensorData)
{
  for (minitensor::Index i = 0; i < N; ++i)
    for (minitensor::Index j = 0; j < N; ++j)
      tensorData[i*N+j] = tensor(i,j);
}

std::array<double,9> to_array(const minitensor::Tensor<double, 3> & tensor);
std::array<double,4> to_array(const minitensor::Tensor<double, 2> & tensor);

// Metric mapping F_M with F_M^T F_M = metric and det(F_M) > 0.
minitensor::Tensor<double, 3> compute_metric_mapping(const minitensor::Tensor<double, 3> & metric);
minitensor::Tensor<double, 2> compute_metric_mapping_2d(const minitensor::Tensor<double, 2> & metric);

// Eigen-decomposition (eigenvectors V, eigenvalues D) of a symmetric tensor A = V D V^T, as
// minitensor::eig_sym, but correct for 2x2 tensors (see the workaround in the source).
template<minitensor::Index N>
std::pair<minitensor::Tensor<double, N>, minitensor::Tensor<double, N>> symmetric_eigen_decomposition(const minitensor::Tensor<double, N> & A);

// V fn(D) V^T for a symmetric tensor A = V D V^T, for example its logarithm or exponential.
template<minitensor::Index N>
minitensor::Tensor<double, N> apply_to_symmetric_tensor_eigenvalues(const minitensor::Tensor<double, N> & A, const std::function<double(double)> & fn);

// Log-Euclidean mean exp(mean(log(M_i))) of N*N row-major metrics M_i (the same average that ngs uses to
// interpolate metrics).
template<minitensor::Index N>
minitensor::Tensor<double, N> log_euclidean_mean_of_metrics(const std::vector<const double *> & metrics);

// As above, for spatialDim*spatialDim row-major metrics, returning the mean as row-major data.
std::vector<double> compute_log_euclidean_mean_of_metrics(const unsigned spatialDim, const std::vector<const double *> & metrics);

// Metric mapping F_M of an element, with M = F_M^T F_M the log-Euclidean mean of the metrics M_n at its
// nodes.  Each nodal metric, and the metric mapping, are spatialDim*spatialDim row-major matrices.
std::vector<double> compute_metric_mapping_data_from_nodal_metrics(
    const unsigned spatialDim,
    const std::vector<const double *> & nodalMetrics);

// Whether a spatialDim*spatialDim row-major matrix is a valid metric: symmetric (to a relative
// tolerance of the largest component) and positive definite.
bool is_symmetric_positive_definite_metric(
    const unsigned spatialDim,
    const double * metric);

// Isotropic metric I/h^2 for size h.
minitensor::Tensor<double, 3> metric_tensor_for_size(const double h);
minitensor::Tensor<double, 2> metric_tensor_for_size_2d(const double h);

// Length of the vector edge in the metric M (a spatialDim*spatialDim row-major matrix), sqrt(edge^T M edge):
// the length of edge relative to the size the metric specifies in its direction.  In 2D only the first two
// components of edge are used.
double compute_relative_length_in_metric(
    const unsigned spatialDim,
    const double * metric,
    const stk::math::Vector3d & edge);

// The logarithmic mean (a - b)/ln(a/b) of two non-negative lengths, or their mean (a + b)/2 when they differ
// by less than 1e-3 (absolute, in relative length units), as ngs computes it.
double logarithmic_mean(const double a, const double b);

// Relative length of the edge from x0 to x1, with the metrics metric0 at x0 and metric1 at x1 (row-major):
// the logarithmic mean of its lengths in the two metrics.  This is the integral of the relative length along
// the edge when the size varies geometrically from x0 to x1, and matches how ngs (e.g. improve_mesh) measures
// edges in a nodal metric.
double compute_relative_edge_length(
    const unsigned spatialDim,
    const stk::math::Vector3d & x0,
    const double * metric0,
    const stk::math::Vector3d & x1,
    const double * metric1);

} // namespace metric
} // namespace krino

#endif /* KRINO_KRINO_METRIC_AKRI_METRICTENSORUTILS_HPP_ */
