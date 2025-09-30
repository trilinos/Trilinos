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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTHEX8_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTHEX8_HPP_

#include "Legendre.hpp"
#include "MockMasterElement.hpp"

#include <math.h>  // for sqrt
#include <stddef.h>
#include <algorithm>  // for sort, max, min
#include <array>
#include <cmath>
#include <cstddef>  // for size_t
#include <cstdint>  // for int64_t, uint64_t
#include <iomanip>
#include <iostream>
#include <limits>  // for numeric_limits
#include <memory>  // for __shared_ptr_ac...
#include <sstream>
#include <stdexcept>  // for logic_error
#include <string>     // for string, basic_s...
#include <typeinfo>   // for type_info
#include <utility>    // for move, pair
#include <vector>     // for vector, swap

namespace stk {
namespace unit_test_util {

class Hex8
{
 public:
  static double invSqrt(double x)
  {
    // use the bit-shifting magic of the fast inverse square root method for 64-bit numbers

    union {
      double f;
      std::int64_t i;
    } conv;

    double x2 = 0.5 * x;
    conv.f = x;
    conv.i = 0x5fe6eb50c7b537a9 - (conv.i >> 1);
    conv.f = conv.f * (1.5 - (x2 * conv.f * conv.f));
    return conv.f;
  }

  template <typename T, std::size_t N>
  static T vector_norm(std::array<T, N>& x)
  {
    T norm_sq = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
      norm_sq += x[i] * x[i];
    }
    return norm_sq;
  }

  static bool within_tol(const double& value, const double& tolerance) { return (std::fabs(value) < tolerance); }

  static double parametric_distance(const std::array<double, 3>& x)
  {
    std::array<double, 3> y = {{std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2])}};

    double d = 0.0;
    for (int i = 0; i < 3; ++i) {
      if (d < y[i]) {
        d = y[i];
      }
    }
    return d;
  }

  static double is_in_element(const double* elem_nodal_coor,  // (8,3)
                              const double* point_coor,       // (3)
                              double* par_coor)
  {
    const double isInElemConverged = 1.0e-16;
    // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)
    double x[] = {0., 0.125 * (elem_nodal_coor[1] - elem_nodal_coor[0]),
        0.125 * (elem_nodal_coor[2] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[3] - elem_nodal_coor[0]),
        0.125 * (elem_nodal_coor[4] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[5] - elem_nodal_coor[0]),
        0.125 * (elem_nodal_coor[6] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[7] - elem_nodal_coor[0])};
    double y[] = {0., 0.125 * (elem_nodal_coor[9] - elem_nodal_coor[8]),
        0.125 * (elem_nodal_coor[10] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[11] - elem_nodal_coor[8]),
        0.125 * (elem_nodal_coor[12] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[13] - elem_nodal_coor[8]),
        0.125 * (elem_nodal_coor[14] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[15] - elem_nodal_coor[8])};
    double z[] = {0., 0.125 * (elem_nodal_coor[17] - elem_nodal_coor[16]),
        0.125 * (elem_nodal_coor[18] - elem_nodal_coor[16]), 0.125 * (elem_nodal_coor[19] - elem_nodal_coor[16]),
        0.125 * (elem_nodal_coor[20] - elem_nodal_coor[16]), 0.125 * (elem_nodal_coor[21] - elem_nodal_coor[16]),
        0.125 * (elem_nodal_coor[22] - elem_nodal_coor[16]), 0.125 * (elem_nodal_coor[23] - elem_nodal_coor[16])};

    // (xp,yp,zp) is the point at which we're searching for (xi,eta,zeta)
    // (must translate this also)
    double xp = point_coor[0] - elem_nodal_coor[0];
    double yp = point_coor[1] - elem_nodal_coor[8];
    double zp = point_coor[2] - elem_nodal_coor[16];

    // Newton-Raphson iteration for (xi,eta,zeta)
    double j[9];
    double f[3];
    double shapefct[8];
    double xinew = 0.5;  // initial guess
    double etanew = 0.5;
    double zetanew = 0.5;
    double xicur = xinew;
    double etacur = etanew;
    double zetacur = zetanew;
    std::array<double, 3> xidiff = {{1.0, 1.0, 1.0}};
    unsigned i = 1;
    const unsigned MAX_NR_ITER = 100;

    double xp8 = 0.125 * xp;
    double yp8 = 0.125 * yp;
    double zp8 = 0.125 * zp;

    constexpr std::array<std::array<int, 4>, 5> t2n = {
        {{0, 1, 3, 4}, {4, 1, 5, 6}, {7, 3, 6, 4}, {6, 1, 2, 3}, {6, 1, 3, 4}}};
    constexpr std::array<std::array<int, 12>, 5> tmat = {{{2, 0, 0, -1, 0, 2, 0, -1, 0, 0, 2, -1},
        {2, 2, 2, -1, 0, 0, 2, -1, -2, 0, 0, 1}, {0, 2, 0, -1, 0, 0, -2, 1, -2, 0, 0, 1},
        {0, 0, -2, 1, -2, 0, 0, 1, -2, -2, -2, 1}, {0, -2, -2, 1, -2, 0, -2, 1, -2, -2, 0, 1}}};

    // Break the hex into five tets, and search inside each
    bool found = false;
    for (int tindex = 0; tindex < 5; tindex++) {
      double a11 = x[t2n[tindex][1]] - x[t2n[tindex][0]];
      double a21 = y[t2n[tindex][1]] - y[t2n[tindex][0]];
      double a31 = z[t2n[tindex][1]] - z[t2n[tindex][0]];
      double a12 = x[t2n[tindex][2]] - x[t2n[tindex][0]];
      double a22 = y[t2n[tindex][2]] - y[t2n[tindex][0]];
      double a32 = z[t2n[tindex][2]] - z[t2n[tindex][0]];
      double a13 = x[t2n[tindex][3]] - x[t2n[tindex][0]];
      double a23 = y[t2n[tindex][3]] - y[t2n[tindex][0]];
      double a33 = z[t2n[tindex][3]] - z[t2n[tindex][0]];
      double f1 = xp8 - x[t2n[tindex][0]];
      double f2 = yp8 - y[t2n[tindex][0]];
      double f3 = zp8 - z[t2n[tindex][0]];
      double oden =
          1.0 / (a31 * (a13 * a22 - a12 * a23) + a32 * (a11 * a23 - a13 * a21) + a33 * (a12 * a21 - a11 * a22));
      double myxi = (f1 * (a23 * a32 - a22 * a33) + f2 * (a12 * a33 - a13 * a32) + f3 * (a13 * a22 - a12 * a23)) * oden;
      double myeta =
          -(f1 * (a23 * a31 - a21 * a33) + f2 * (a11 * a33 - a13 * a31) + f3 * (a13 * a21 - a11 * a23)) * oden;
      double myzeta =
          (f1 * (a22 * a31 - a21 * a32) + f2 * (a11 * a32 - a12 * a31) + f3 * (a12 * a21 - a11 * a22)) * oden;

      if (myxi >= 0 && myeta >= 0 && myzeta >= 0 && myzeta <= 1.0 - myxi - myeta) {
        xicur = tmat[tindex][0] * myxi + tmat[tindex][1] * myeta + tmat[tindex][2] * myzeta + tmat[tindex][3];
        etacur = tmat[tindex][4] * myxi + tmat[tindex][5] * myeta + tmat[tindex][6] * myzeta + tmat[tindex][7];
        zetacur = tmat[tindex][8] * myxi + tmat[tindex][9] * myeta + tmat[tindex][10] * myzeta + tmat[tindex][11];
        found = true;
        break;
      }
    }

    // If the point is not found inside any of the tetrahedra, fall back to IDW
    if (!found) {
      double w0 = invSqrt((xp8 - x[0]) * (xp8 - x[0]) + (yp8 - y[0]) * (yp8 - y[0]) + (zp8 - z[0]) * (zp8 - z[0]));
      double w1 = invSqrt((xp8 - x[1]) * (xp8 - x[1]) + (yp8 - y[1]) * (yp8 - y[1]) + (zp8 - z[1]) * (zp8 - z[1]));
      double w2 = invSqrt((xp8 - x[2]) * (xp8 - x[2]) + (yp8 - y[2]) * (yp8 - y[2]) + (zp8 - z[2]) * (zp8 - z[2]));
      double w3 = invSqrt((xp8 - x[3]) * (xp8 - x[3]) + (yp8 - y[3]) * (yp8 - y[3]) + (zp8 - z[3]) * (zp8 - z[3]));
      double w4 = invSqrt((xp8 - x[4]) * (xp8 - x[4]) + (yp8 - y[4]) * (yp8 - y[4]) + (zp8 - z[4]) * (zp8 - z[4]));
      double w5 = invSqrt((xp8 - x[5]) * (xp8 - x[5]) + (yp8 - y[5]) * (yp8 - y[5]) + (zp8 - z[5]) * (zp8 - z[5]));
      double w6 = invSqrt((xp8 - x[6]) * (xp8 - x[6]) + (yp8 - y[6]) * (yp8 - y[6]) + (zp8 - z[6]) * (zp8 - z[6]));
      double w7 = invSqrt((xp8 - x[7]) * (xp8 - x[7]) + (yp8 - y[7]) * (yp8 - y[7]) + (zp8 - z[7]) * (zp8 - z[7]));

      double wt = 1.0 / (w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7);
      double p6m0 = w6 - w0;
      double p7m1 = w7 - w1;
      double p2m4 = w2 - w4;
      double p5m3 = w5 - w3;
      xicur = (p6m0 - p7m1 + p2m4 + p5m3) * wt;
      etacur = (p6m0 + p7m1 + p2m4 - p5m3) * wt;
      zetacur = (p6m0 + p7m1 - p2m4 + p5m3) * wt;
    }

    // Constants for the iteration
    double x3mx2 = x[3] - x[2];
    double x4mx5 = x[4] - x[5];
    double x7mx6 = x[7] - x[6];
    double x1mx2 = x[1] - x[2];
    double x4mx7 = x[4] - x[7];
    double x5mx6 = x[5] - x[6];
    double x1mx5 = x[1] - x[5];
    double x2mx6 = x[2] - x[6];
    double x3mx7 = x[3] - x[7];

    double y3my2 = y[3] - y[2];
    double y4my5 = y[4] - y[5];
    double y7my6 = y[7] - y[6];
    double y1my2 = y[1] - y[2];
    double y4my7 = y[4] - y[7];
    double y5my6 = y[5] - y[6];
    double y1my5 = y[1] - y[5];
    double y2my6 = y[2] - y[6];
    double y3my7 = y[3] - y[7];

    double z3mz2 = z[3] - z[2];
    double z4mz5 = z[4] - z[5];
    double z7mz6 = z[7] - z[6];
    double z1mz2 = z[1] - z[2];
    double z4mz7 = z[4] - z[7];
    double z5mz6 = z[5] - z[6];
    double z1mz5 = z[1] - z[5];
    double z2mz6 = z[2] - z[6];
    double z3mz7 = z[3] - z[7];

    // Actual NR iteration
    do {
      double one_minu_xi = 1.0 - xicur;
      double one_plus_xi = 1.0 + xicur;
      double one_minu_eta = 1.0 - etacur;
      double one_plus_eta = 1.0 + etacur;
      double one_minu_zeta = 1.0 - zetacur;
      double one_plus_zeta = 1.0 + zetacur;

      double memz = one_minu_eta * one_minu_zeta;
      double mepz = one_minu_eta * one_plus_zeta;
      double pepz = one_plus_eta * one_plus_zeta;
      double pemz = one_plus_eta * one_minu_zeta;

      double mxmz = one_minu_xi * one_minu_zeta;
      double mxpz = one_minu_xi * one_plus_zeta;
      double pxpz = one_plus_xi * one_plus_zeta;
      double pxmz = one_plus_xi * one_minu_zeta;

      double mxme = one_minu_xi * one_minu_eta;
      double mxpe = one_minu_xi * one_plus_eta;
      double pxpe = one_plus_xi * one_plus_eta;
      double pxme = one_plus_xi * one_minu_eta;

      j[0] = -memz * x[1] + pemz * x3mx2 + mepz * x4mx5 + pepz * x7mx6;
      j[1] = pxmz * x1mx2 - mxmz * x[3] + mxpz * x4mx7 + pxpz * x5mx6;
      j[2] = pxme * x1mx5 + pxpe * x2mx6 + mxpe * x3mx7 - mxme * x[4];
      j[3] = -memz * y[1] + pemz * y3my2 + mepz * y4my5 + pepz * y7my6;
      j[4] = pxmz * y1my2 - mxmz * y[3] + mxpz * y4my7 + pxpz * y5my6;
      j[5] = pxme * y1my5 + pxpe * y2my6 + mxpe * y3my7 - mxme * y[4];
      j[6] = -memz * z[1] + pemz * z3mz2 + mepz * z4mz5 + pepz * z7mz6;
      j[7] = pxmz * z1mz2 - mxmz * z[3] + mxpz * z4mz7 + pxpz * z5mz6;
      j[8] = pxme * z1mz5 + pxpe * z2mz6 + mxpe * z3mz7 - mxme * z[4];

      double jdet = -(j[2] * j[4] * j[6]) + j[1] * j[5] * j[6] + j[2] * j[3] * j[7] - j[0] * j[5] * j[7] -
                    j[1] * j[3] * j[8] + j[0] * j[4] * j[8];
      double odet = 1.0 / jdet;

      if (!jdet) {
        i = MAX_NR_ITER;
        break;
      }

      shapefct[0] = mxme * one_minu_zeta;
      shapefct[1] = pxme * one_minu_zeta;
      shapefct[2] = pxpe * one_minu_zeta;
      shapefct[3] = mxpe * one_minu_zeta;
      shapefct[4] = mxme * one_plus_zeta;
      shapefct[5] = pxme * one_plus_zeta;
      shapefct[6] = pxpe * one_plus_zeta;
      shapefct[7] = mxpe * one_plus_zeta;

      f[0] = xp - shapefct[1] * x[1] - shapefct[2] * x[2] - shapefct[3] * x[3] - shapefct[4] * x[4] -
             shapefct[5] * x[5] - shapefct[6] * x[6] - shapefct[7] * x[7];
      f[1] = yp - shapefct[1] * y[1] - shapefct[2] * y[2] - shapefct[3] * y[3] - shapefct[4] * y[4] -
             shapefct[5] * y[5] - shapefct[6] * y[6] - shapefct[7] * y[7];
      f[2] = zp - shapefct[1] * z[1] - shapefct[2] * z[2] - shapefct[3] * z[3] - shapefct[4] * z[4] -
             shapefct[5] * z[5] - shapefct[6] * z[6] - shapefct[7] * z[7];

      double relax = 1.0;
      xinew = xicur + relax *
                          (f[2] * (j[2] * j[4] - j[1] * j[5]) + f[1] * (j[1] * j[8] - j[2] * j[7]) +
                              f[0] * (j[5] * j[7] - j[4] * j[8])) *
                          odet;
      etanew = etacur + relax *
                            (f[2] * (-j[2] * j[3] + j[0] * j[5]) + f[1] * (j[2] * j[6] - j[0] * j[8]) +
                                f[0] * (j[3] * j[8] - j[5] * j[6])) *
                            odet;
      zetanew = zetacur + relax *
                              (f[2] * (j[1] * j[3] - j[0] * j[4]) + f[1] * (j[0] * j[7] - j[1] * j[6]) +
                                  f[0] * (j[4] * j[6] - j[3] * j[7])) *
                              odet;

      xidiff[0] = xinew - xicur;
      xidiff[1] = etanew - etacur;
      xidiff[2] = zetanew - zetacur;
      xicur = xinew;
      etacur = etanew;
      zetacur = zetanew;
    } while (!within_tol(vector_norm(xidiff), isInElemConverged) && ++i < MAX_NR_ITER);

    par_coor[0] = par_coor[1] = par_coor[2] = std::numeric_limits<double>::max();
    double dist = std::numeric_limits<double>::max();

    if (i < MAX_NR_ITER) {
      par_coor[0] = xinew;
      par_coor[1] = etanew;
      par_coor[2] = zetanew;

      std::array<double, 3> xtmp = {{par_coor[0], par_coor[1], par_coor[2]}};
      dist = parametric_distance(xtmp);
    }

    return dist;
  }

  static const std::vector<double>& coordinate_center()
  {
    static const std::vector<double> C(3, 0.);
    return C;
  }

  static void interpolate_point(const double* par_coord, // (3)
                                const int& ncomp_field,
                                const double* field,  // (8,ncomp_field)
                                double* result)       // (ncomp_field)
  {
    // 'field' is a flat array of dimension (8,ncomp_field) (Fortran ordering);
    double xi = par_coord[0];
    double eta = par_coord[1];
    double zeta = par_coord[2];

    for(int i = 0; i < ncomp_field; i++) {
      // Base 'field array' index for ith component
      int b = 8 * i;

      result[i] = 0.125 * (1.0 - eta) * (1.0 - xi) * (1.0 - zeta) * field[b + 0] +
                  0.125 * (1.0 - eta) * (1.0 + xi) * (1.0 - zeta) * field[b + 1] +
                  0.125 * (1.0 + eta) * (1.0 + xi) * (1.0 - zeta) * field[b + 2] +
                  0.125 * (1.0 + eta) * (1.0 - xi) * (1.0 - zeta) * field[b + 3] +
                  0.125 * (1.0 - eta) * (1.0 - xi) * (1.0 + zeta) * field[b + 4] +
                  0.125 * (1.0 - eta) * (1.0 + xi) * (1.0 + zeta) * field[b + 5] +
                  0.125 * (1.0 + eta) * (1.0 + xi) * (1.0 + zeta) * field[b + 6] +
                  0.125 * (1.0 + eta) * (1.0 - xi) * (1.0 + zeta) * field[b + 7];
    }
  }
};

/**
 * A 3D Gauss-Legendre quadrature rule (traditionally
 * called the Gauss quadrature) of arbitrary order q x q x q, on the
 * interval [-1,1] x [-1,1] x [-1,1].
 *
 */
class Hex8GaussQuadrature : public GaussQuadrature {
public:
  Hex8GaussQuadrature(unsigned q)
  {
    m_order = q;
    m_numIntgPoints = q * q * q;
    m_numParametricCoordinates = 3;

    // initialize the points and weights
    gauss_legendre_3D(q, m_intgLocations, m_intgWeights);
  }

  ~Hex8GaussQuadrature() = default;
};

class MasterElementHex8 : public MasterElement {
 public:

  MasterElementHex8(const unsigned integrationOrder)
  : MasterElement(stk::topology::HEX_8)
  {
    m_name = "MasterElementHex8";
    m_integrationOrder = get_integration_order(integrationOrder);
    m_quadrature = std::make_shared<Hex8GaussQuadrature>(m_integrationOrder);
  }

  MasterElementHex8()
  : MasterElement(stk::topology::HEX_8)
  {
    m_name = "MasterElementHex8";
    m_integrationOrder = get_integration_order(0);
    m_quadrature = std::make_shared<Hex8GaussQuadrature>(m_integrationOrder);
  }

  ~MasterElementHex8() override = default;

  const std::vector<double>& coordinate_center() const override { return Hex8::coordinate_center(); }

  double is_in_element(const double* elem_nodal_coor,
                       const double* point_coor,
                       double* par_coor) const override
  {
    return Hex8::is_in_element(elem_nodal_coor, point_coor, par_coor);
  }

  void interpolate_point(const double* par_coord,
                         const int& ncomp_field,
                         const double* field,
                         double* result) const override
  {
    Hex8::interpolate_point(par_coord, ncomp_field, field, result);
  }
};

}
}

#endif
