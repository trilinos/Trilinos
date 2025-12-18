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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTLINE2_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTLINE2_HPP_

#include "stk_transfer_util/MockMasterElement.hpp"
#include <stk_util/util/ReportHandler.hpp>  // for eval_test_condition, STK_...

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

class Line2
{
 public:

  static double parametric_distance(const std::array<double, 2>& x)
  {
    const double ELEMENT_THICKNESS = 0.01;
    double dist = std::fabs(x[0]);
    if(ELEMENT_THICKNESS < x[1] && dist < 1 + x[1]) dist = 1 + x[1];
    return dist;
  }

  static double is_in_element(const double* elem_nodal_coor, // (3,3)
                              const double* point_coor,      // (3)
                              double* par_coor)
  {
    // elem_nodal_coor has the endpoints of the line
    // segment defining this element.  Set the first
    // endpoint to zero.  This means subtrace the
    // first endpoint from the second.
    const double X1 = elem_nodal_coor[1] - elem_nodal_coor[0];
    const double X2 = elem_nodal_coor[3] - elem_nodal_coor[2];
    const double X3 = elem_nodal_coor[5] - elem_nodal_coor[4];

    // Now subtract the first endpoint from the target point
    const double P1 = point_coor[0] - elem_nodal_coor[0];
    const double P2 = point_coor[1] - elem_nodal_coor[2];
    const double P3 = point_coor[2] - elem_nodal_coor[4];

    // Now find the projection along the line of the point
    // This is the parametric coordinate in range (0,1)
    const double norm2 = X1 * X1 + X2 * X2 + X3 * X3;
    STK_ThrowAssert(norm2 != 0.0);

    const double xi = (P1 * X1 + P2 * X2 + P3 * X3) / norm2;
    // rescale to (-1,1)
    par_coor[0] = 2.0 * xi - 1.0;

    // d = ||P x X||/||X||
    // gives the normalized distance from the point to the element.

    const double alpha = std::sqrt((P2 * X1 - P1 * X2) * (P2 * X1 - P1 * X2) + (-P3 * X1 + P1 * X3) * (-P3 * X1 + P1 * X3) +
                                   (P3 * X2 - P2 * X3) * (P3 * X2 - P2 * X3)) / norm2;

    std::array x = { par_coor[0], alpha };
    const double dist = parametric_distance(x);

    return dist;
  }

  static void interpolate_point(const double* par_coord, // (1)
                                const int& ncomp_field,
                                const double* field,     // (2,ncomp_field)
                                double* result)          // (ncomp_field)
  {
    // 'field' is a flat array of dimension (2,ncomp_field) (Fortran ordering);
    double one = 1.0;

    double s = par_coord[0];

    for(int i = 0; i < ncomp_field; i++) {
      int b = 2 * i; // Base 'field array' index for ith component

      result[i] = 0.5 * (one - s) * field[b + 0] + 0.5 * (one + s) * field[b + 1];
    }
  }

  static const std::vector<double>& coordinate_center()
  {
    static const std::vector<double> C(1, 0.);
    return C;
  }
};

/**
 * A 1D Gauss-Legendre quadrature rule (traditionally
 * called the Gauss quadrature) of arbitrary order q, on the
 * interval [-1,1].
 *
 */
class Line2GaussQuadrature : public GaussQuadrature {
public:
  Line2GaussQuadrature(unsigned q)
  {
    m_order = q;
    m_numIntgPoints = q;
    m_numParametricCoordinates = 1;

    // initialize the points and weights
    gauss_legendre_1D(q, m_intgLocations, m_intgWeights);
  }

  ~Line2GaussQuadrature() = default;
};

class MasterElementLine2 : public MasterElement {
 public:

  MasterElementLine2(const unsigned integrationOrder)
    : MasterElement(stk::topology::LINE_2)
  {
    m_name = "MasterElementLine2";
    m_integrationOrder = get_integration_order(integrationOrder);
    m_quadrature = std::make_shared<Line2GaussQuadrature>(m_integrationOrder);
  }

  MasterElementLine2()
  : MasterElement(stk::topology::LINE_2)
  {
    m_name = "MasterElementLine2";
    m_integrationOrder = get_integration_order(0);
    m_quadrature = std::make_shared<Line2GaussQuadrature>(m_integrationOrder);
  }

  ~MasterElementLine2() override = default;

  const std::vector<double>& coordinate_center() const override { return Line2::coordinate_center(); }

  double is_in_element(const double* elem_nodal_coor,
                       const double* point_coor,
                       double* par_coor) const override
  {
    return Line2::is_in_element(elem_nodal_coor, point_coor, par_coor);
  }

  void interpolate_point(const double* par_coord,
                         const int& ncomp_field,
                         const double* field,
                         double* result) const override
  {
    Line2::interpolate_point(par_coord, ncomp_field, field, result);
  }
};

}
}

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKELEMENTLINE2_HPP_ */
