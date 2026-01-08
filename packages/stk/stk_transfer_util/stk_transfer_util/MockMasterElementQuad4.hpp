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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTQUAD4_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTQUAD4_HPP_

#include "stk_transfer_util/MockMasterElement.hpp"

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
namespace transfer_util {

class Quad4
{
 public:

  static bool within_tol(const double& value, const double& tolerance)
  {
    return ( std::fabs(value) < tolerance );
  }

  template<typename T, std::size_t N>
  static T vector_norm(std::array<T,N> &x)
  {
    T norm_sq = 0.0;
    for (std::size_t i = 0; i < N; ++i ) {
      norm_sq += x[i] * x[i];
    }
    return norm_sq;
  }

  static double vector_norm(const double* vect, int N)
  {
    double norm_sq = 0.0;
    for (int i = 0; i < N; ++i ) {
      norm_sq += vect[i] * vect[i];
    }
    return norm_sq;
  }

  static double linear_quad_parametric_distance(const std::array<double, 3>& x)
  {
    const double ELEM_THICK = 0.01;
    std::array<double, 3> y = { { std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) } };
    double d = y[0];
    if(d < y[1]) d = y[1];
    if(ELEM_THICK < y[2] && d < 1 + y[2]) d = 1 + y[2];
    return d;
  }

  static double compute_linear_quad_cpp_distance(const double* normal, const double* solcur, const double* deltasol,
                                                 double* par_coor)
  {
    // Rescale the distance vector by the length of the (non-unit) normal vector,
    // which was used above in the NR iteration.
    const double area = std::sqrt(vector_norm(normal, 3));
    const double length = std::sqrt(area);
    const double par_coor_2 = (solcur[2] + deltasol[2]) * length;

    std::array<double, 3> xtmp = { { par_coor[0], par_coor[1], par_coor_2 } };
    return linear_quad_parametric_distance(xtmp);
  }

  static void non_unit_face_normal(const double* par_coord,       // (2)
                                   const double* elem_nodal_coor, // (4,3)
                                   double* normal_vector)         // (3)
  {
    double xi  = par_coord[0];
    double eta = par_coord[1];

    // Translate element so that node 0 is at (x,y,z) = (0,0,0)
    double x[3] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                    elem_nodal_coor[2] - elem_nodal_coor[0],
                    elem_nodal_coor[3] - elem_nodal_coor[0] };

    double y[3] = { elem_nodal_coor[5] - elem_nodal_coor[4],
                    elem_nodal_coor[6] - elem_nodal_coor[4],
                    elem_nodal_coor[7] - elem_nodal_coor[4] };

    double z[3] = { elem_nodal_coor[9] - elem_nodal_coor[8],
                    elem_nodal_coor[10] - elem_nodal_coor[8],
                    elem_nodal_coor[11] - elem_nodal_coor[8] };

    // Mathematica-generated and simplified code for the normal vector

    double n0 = 0.125* (xi * y[2] * z[0] + y[0] * z[1] + xi * y[0] * z[1] - y[2] * z[1] - xi * y[0] * z[2] +
                        y[1] * (-((1.0 + xi) * z[0]) + (1.0 + eta) * z[2]) +
                        eta * (y[2] * z[0] - y[2] * z[1] - y[0] * z[2]));

    double n1 = 0.125* (-(xi * x[2] * z[0]) - x[0] * z[1] - xi * x[0] * z[1] + x[2] * z[1] + xi * x[0] * z[2] +
                        x[1] * ((1.0 + xi) * z[0] - (1.0 + eta) * z[2]) +
                        eta * (-(x[2] * z[0]) + x[2] * z[1] + x[0] * z[2]));

    double n2 = 0.125* (xi * x[2] * y[0] + x[0] * y[1] + xi * x[0] * y[1] - x[2] * y[1] - xi * x[0] * y[2] +
                        x[1] * (-((1.0 + xi) * y[0]) + (1.0 + eta) * y[2]) +
                        eta * (x[2] * y[0] - x[2] * y[1] - x[0] * y[2]));

    normal_vector[0] = n0;
    normal_vector[1] = n1;
    normal_vector[2] = n2;
  }

  static void unit_face_normal(
      const double* par_coord,
      const double* elem_nodal_coor, // (4,3)
      double* normal_vector)         // (3)
  {
    non_unit_face_normal(par_coord, elem_nodal_coor, normal_vector);

    double nlen = std::sqrt(vector_norm(normal_vector, 3));

    normal_vector[0] /= nlen;
    normal_vector[1] /= nlen;
    normal_vector[2] /= nlen;
  }

  static void interpolate_point(const double* par_coord, // (2)
                                const int& ncomp_field,
                                const double* field,  // (4,ncomp_field)
                                double* result)       // (ncomp_field)
  {
    // 'field' is a flat array of dimension (4,ncomp_field) (Fortran ordering);

    double xi = par_coord[0];
    double eta = par_coord[1];

    for(int i = 0; i < ncomp_field; i++) {
      // Base 'field array' index for ith component
      int b = 4 * i;

      result[i] = 0.25 * ((1.0 - eta) * (1.0 - xi) * field[b + 0] +
                          (1.0 - eta) * (1.0 + xi) * field[b + 1] +
                          (1.0 + eta) * (1.0 + xi) * field[b + 2] +
                          (1.0 + eta) * (1.0 - xi) * field[b + 3]);
    }
  }

  static double is_in_element(const double* elem_nodal_coor, // (4,3)
                              const double* point_coor,      // (3)
                              double* par_coor)
  {
    const double isInElemConverged = 1.0e-16;
    // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

    double x[3] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                    elem_nodal_coor[2] - elem_nodal_coor[0],
                    elem_nodal_coor[3] - elem_nodal_coor[0] };

    double y[3] = { elem_nodal_coor[5] - elem_nodal_coor[4],
                    elem_nodal_coor[6] - elem_nodal_coor[4],
                    elem_nodal_coor[7] - elem_nodal_coor[4] };

    double z[3] = { elem_nodal_coor[9] - elem_nodal_coor[8],
                    elem_nodal_coor[10] - elem_nodal_coor[8],
                    elem_nodal_coor[11] - elem_nodal_coor[8] };

    // (xp,yp,zp) is the point at which we're searching for (xi,eta,d)
    // (must translate this also)
    // d = (scaled) distance in (x,y,z) space from point (xp,yp,zp) to the
    //     surface defined by the face element (the distance is scaled by
    //     the length of the non-unit normal vector; rescaling of d is done
    //     following the NR iteration below).

    double xp = point_coor[0] - elem_nodal_coor[0];
    double yp = point_coor[1] - elem_nodal_coor[4];
    double zp = point_coor[2] - elem_nodal_coor[8];

    // Newton-Raphson iteration for (xi,eta,d)

    double j[9];
    double gn[3];
    double xcur[3];   // current (x,y,z) point on element surface
    double normal[3]; // (non-unit) normal computed at xcur

    // Solution vector solcur[3] = {xi,eta,d}
    double solcur[3] = { -0.5, -0.5, -0.5 }; // initial guess
    std::array<double,3> deltasol = {{ 1.0, 1.0, 1.0 }};

    unsigned i = 0;
    const unsigned MAX_NR_ITER = 100;
    do {
      // Update guess vector
      solcur[0] += deltasol[0];
      solcur[1] += deltasol[1];
      solcur[2] += deltasol[2];

      interpolate_point(solcur, 3, elem_nodal_coor, xcur);

      // Translate xcur ((x,y,z) point corresponding
      // to current (xi,eta) guess)

      xcur[0] -= elem_nodal_coor[0];
      xcur[1] -= elem_nodal_coor[4];
      xcur[2] -= elem_nodal_coor[8];

      non_unit_face_normal(solcur, elem_nodal_coor, normal);

      gn[0] = xcur[0] - xp + solcur[2] * normal[0];
      gn[1] = xcur[1] - yp + solcur[2] * normal[1];
      gn[2] = xcur[2] - zp + solcur[2] * normal[2];

      // Mathematica-generated code for the jacobian

      j[0] = 0.125* (-2.0* (-1.0 + solcur[1]) * x[0] + (2.0* (1.0 +
             solcur[1]) * (x[1] - x[2]) + solcur[2] * (-(y[1] * z[0]) + y[2] * z[0] + y[0] * z[1] - y[0] * z[2])));

      j[1] = 0.125* (-2.0* (1.0+ solcur[0]) * x[0] + 2.0* (1.0 + solcur[0]) * x[1] - 2.0* (-1.0 + solcur[0]) *
                                      x[2] + (solcur[2] * (y[2] * (z[0] - z[1]) + (-y[0] + y[1]) * z[2])));

      j[2] = normal[0];

      j[3] = 0.125* (-2.0* (-1.0 + solcur[1]) * y[0] + (2.0* (1.0 +
             solcur[1]) * (y[1] - y[2]) + solcur[2] * (x[1] * z[0] - x[2] * z[0] - x[0] * z[1] + x[0] * z[2])));

      j[4] = 0.125* (-2.0* (1.0+ solcur[0]) * y[0] +
                                  2.0* (1.0 + solcur[0]) * y[1] -
                                  2.0* (-1.0 + solcur[0]) * y[2] +
                                  (solcur[2] * (x[2] * (-z[0] + z[1]) + (x[0] - x[1]) * z[2])));

      j[5] = normal[1];

      j[6] = 0.125* ((solcur[2] * (-(x[1] * y[0]) + x[2] * y[0] + x[0] * y[1] - x[0] * y[2])) - 2.0* ((-1.0 + solcur[1]) * z[0] - (1.0 + solcur[1]) * (z[1] - z[2])));

      j[7] = 0.125* ((solcur[2] * (x[2] * (y[0] - y[1]) + (-x[0] + x[1]) * y[2])) - 2.0* (1.0 + solcur[0]) * z[0] + 2.0* (1.0 + solcur[0]) *
                                      z[1] - 2.0* (-1.0 + solcur[0]) * z[2]);

      j[8] = normal[2];

      double jdet = -(j[2] * j[4] * j[6]) + j[1] * j[5] * j[6] + j[2] * j[3] * j[7] - j[0] * j[5] * j[7] - j[1] * j[3] * j[8] +
             j[0] * j[4] * j[8];

      // Solve linear system (j*deltasol = -gn) for deltasol at step n+1

      deltasol[0] = (gn[2] * (j[2] * j[4] - j[1] * j[5]) + gn[1] * (-(j[2] * j[7]) + j[1] * j[8]) +
                     gn[0] * (j[5] * j[7] - j[4] * j[8])) /
                    jdet;
      deltasol[1] = (gn[2] * (-(j[2] * j[3]) + j[0] * j[5]) + gn[1] * (j[2] * j[6] - j[0] * j[8]) +
                     gn[0] * (-(j[5] * j[6]) + j[3] * j[8])) /
                    jdet;
      deltasol[2] = (gn[2] * (j[1] * j[3] - j[0] * j[4]) + gn[1] * (-(j[1] * j[6]) + j[0] * j[7]) +
                     gn[0] * (j[4] * j[6] - j[3] * j[7])) /
                    jdet;

    } while(!within_tol(vector_norm(deltasol), isInElemConverged) && ++i < MAX_NR_ITER);

    // Fill in solution vector; only include the distance (in the third
    // solution vector slot) if npar_coord = 3 (this is how the user
    // requests it)

    par_coor[0] = par_coor[1] = std::numeric_limits<double>::max();

    if(i < MAX_NR_ITER) {
      par_coor[0] = solcur[0] + deltasol[0];
      par_coor[1] = solcur[1] + deltasol[1];
    }

    double dist = compute_linear_quad_cpp_distance(normal, solcur, deltasol.data(), par_coor);

    return dist;
  }

  static const std::vector<double>& coordinate_center()
  {
    static const std::vector<double> C(2, 0.);
    return C;
  }
};


/**
 * A 2D Gauss-Legendre quadrature rule (traditionally
 * called the Gauss quadrature) of arbitrary order q x q, on the
 * interval [-1,1] x [-1,1].
 *
 */
class Quad4GaussQuadrature : public GaussQuadrature {
public:
  Quad4GaussQuadrature(unsigned q)
  {
    m_order = q;
    m_numIntgPoints = q * q;
    m_numParametricCoordinates = 2;

    // initialize the points and weights
    gauss_legendre_2D(q, m_intgLocations, m_intgWeights);
  }

  ~Quad4GaussQuadrature() = default;
};

class MasterElementQuad4 : public MasterElement {
 public:

  MasterElementQuad4(const unsigned integrationOrder)
  : MasterElement(stk::topology::QUAD_4)
  {
     m_name = "MasterElementQuad4";
     m_integrationOrder = get_integration_order(integrationOrder);
     m_quadrature = std::make_shared<Quad4GaussQuadrature>(m_integrationOrder);
  }

  MasterElementQuad4()
  : MasterElement(stk::topology::QUAD_4)
  {
    m_name = "MasterElementQuad4";
    m_integrationOrder = get_integration_order(0);
    m_quadrature = std::make_shared<Quad4GaussQuadrature>(m_integrationOrder);
  }

  ~MasterElementQuad4() override = default;

  const std::vector<double>& coordinate_center() const override { return Quad4::coordinate_center(); }

  double is_in_element(const double* elem_nodal_coor,
                       const double* point_coor,
                       double* par_coor) const override
  {
    return Quad4::is_in_element(elem_nodal_coor, point_coor, par_coor);
  }

  void interpolate_point(const double* par_coord,
                         const int& ncomp_field,
                         const double* field,
                         double* result) const override
  {
    Quad4::interpolate_point(par_coord, ncomp_field, field, result);
  }
};

}
}


#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTQUAD4_HPP_ */
