// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "Legendre.hpp"

#include <cmath>
#include <algorithm>
#include <cassert>                   // for assert
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

double legendre(int order, double x) {
  if (order == 0) return 1.0;
  if (order == 1) return x;
  double Pn_2 = 1.0;
  double Pn_1 = x;
  double Pn = 0.0;
  for (int k = 2; k <= order; ++k) {
    Pn = ((2.0 * k - 1.0) * x * Pn_1 - (k - 1.0) * Pn_2) / k;
    Pn_2 = Pn_1;
    Pn_1 = Pn;
  }
  return Pn;
}

double legendre_derivative(int order, double x) { return order * (x * legendre(order, x) - legendre(order - 1, x)) / (x * x - 1); }

std::vector<double> gauss_weights(int num_intg_points, const std::vector<double>& roots) {
  std::vector<double> weights(num_intg_points);
  for (int i = 0; i < num_intg_points; ++i) {
    double x = roots[i];
    weights[i] =
        2.0 / ((1.0 - x * x) * legendre_derivative(num_intg_points, x) * legendre_derivative(num_intg_points, x));
  }
  return weights;
}

std::vector<double> gauss_abscissas(int num_intg_points) {
  std::vector<double> roots(num_intg_points);
  for (int i = 0; i < num_intg_points; ++i) {
    double x = std::cos(M_PI * (i + 0.75) / (num_intg_points + 0.5));
    double x0 = x;
    do {
      x0 = x;
      x -= legendre(num_intg_points, x) / legendre_derivative(num_intg_points, x);
    } while (std::abs(x - x0) > 1e-15);
    roots[i] = x;
  }
  std::sort(roots.begin(), roots.end());
  return roots;
}

void gauss_legendre_1D(unsigned q, std::vector<double>& points, std::vector<double>& weights)
{
  assert(q > 0);

  points = gauss_abscissas(q);
  weights = gauss_weights(q, points);
}

void gauss_legendre_2D(unsigned q, std::vector<double>& points, std::vector<double>& weights)
{
  assert(q > 0);

  unsigned numIntgPoints(q * q);

  if(points.size() < 2*numIntgPoints) {
    points.resize(2*numIntgPoints);
  }

  if(weights.size() < numIntgPoints) {
    weights.resize(numIntgPoints);
  }

  std::vector<double> g_abscissas = gauss_abscissas(q);
  std::vector<double> g_weights = gauss_weights(q, g_abscissas);

  unsigned l = 0;
  unsigned m = 0;

  for(unsigned i = 0; i < q; i++) {
    for(unsigned j = 0; j < q; j++) {

      points[l++] = g_abscissas[i];
      points[l++] = g_abscissas[j];

      weights[m++] = g_weights[i] * g_weights[j];
    }
  }
}

void gauss_legendre_3D(unsigned q, std::vector<double>& points, std::vector<double>& weights)
{
  assert(q > 0);

  unsigned numIntgPoints(q * q * q);

  if(points.size() < 3*numIntgPoints) {
    points.resize(3*numIntgPoints);
  }

  if(weights.size() < numIntgPoints) {
    weights.resize(numIntgPoints);
  }

  std::vector<double> g_abscissas = gauss_abscissas(q);
  std::vector<double> g_weights = gauss_weights(q, g_abscissas);

  unsigned l = 0;
  unsigned m = 0;

  for(unsigned i = 0; i < q; i++) {
    for(unsigned j = 0; j < q; j++) {
      for(unsigned k = 0; k < q; k++) {

        points[l++] = g_abscissas[i];
        points[l++] = g_abscissas[j];
        points[l++] = g_abscissas[k];

        weights[m++] = g_weights[i] * g_weights[j] * g_weights[k];
      }
    }
  }
}

}
}
