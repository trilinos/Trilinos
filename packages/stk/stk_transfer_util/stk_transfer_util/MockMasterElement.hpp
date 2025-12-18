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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKMASTERELEMENT_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKMASTERELEMENT_HPP_

#include "stk_topology/topology.hpp"              // for topology, topol...
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

/**
 * A Gauss-Legendre quadrature rule (traditionally
 * called the Gauss quadrature) of arbitrary order q on the
 * interval [-1,1] x [-1,1] x [-1,1].
 *
 */
class GaussQuadrature {
public:
  GaussQuadrature() { }

  virtual ~GaussQuadrature() {};

  unsigned order() const { return m_order; }
  unsigned num_intg_pts() const { return m_numIntgPoints; }
  unsigned num_parametric_coordinates() const { return m_numParametricCoordinates; }
  const std::vector<double>& intg_pt_weights() const { return m_intgWeights; }
  const std::vector<double>& intg_pt_locations() const { return m_intgLocations; }

protected:
  unsigned m_order{0};
  unsigned m_numIntgPoints{0};
  unsigned m_numParametricCoordinates{0};

  std::vector<double> m_intgWeights; // array holding values of the weights
  std::vector<double> m_intgLocations; // array holding locations of the pts
};

class MasterElement {
 public:
  //: The constructor takes a name and a mesh object topology;
  //: these are to be supplied by the derived element class.
  MasterElement(const stk::topology topology)
  : m_topology(topology)
  {

  }

  virtual ~MasterElement() {};

  // Copy and assignment are not allowed
  MasterElement(const MasterElement&) = delete;
  MasterElement& operator=(const MasterElement&) = delete;

  static unsigned get_integration_order(const unsigned integrationOrder)
  {
    return integrationOrder == 0 ? 2 : integrationOrder;
  }

  //: Returns the name of the element
  const std::string& name() const { return m_name; }

  //: Returns the underlying topology of the master element.
  stk::topology get_topology() const { return m_topology; }

  //: Return number of integration points/stations
  virtual unsigned num_intg_pts() const
  {
    STK_ThrowRequire(m_quadrature);
    return m_quadrature->num_intg_pts();
  }

  //: Return the number of nodes
  virtual unsigned num_nodes() const { return m_topology.num_nodes(); }

  //: Return parametric locations of integration points/stations
  virtual void intg_pt_locations(std::vector<double>& gaussPoints) const
  {
    STK_ThrowRequire(m_quadrature);
    gaussPoints = m_quadrature->intg_pt_locations();
  }

  virtual unsigned num_parametric_coordinates() const
  {
    STK_ThrowRequire(m_quadrature);
    return m_quadrature->num_parametric_coordinates();
  }

  /* Dimensions for the following methods:                         */
  /*   npar_coord  number of parametric coordinates                */
  /*   npts        number of locations to interpolate per element  */
  /*   ncomp       number of components in interpolated field      */
  /*   nelem       number of elements                              */

  //: Return the logical center coordinate of the element.
  //  This is not necessarily the centroid.  The center coordinate
  //  is used with the scalar return value for the is_in_element()
  //  function to scale a parametric coordinate to the parametric
  //  surface of the master element.  Given the center coordinate
  //  C and the parametric distance d for the parametric point X then:
  //            Y = (X-C)/d + C
  //  will define a point Y on the element parametric surface.
  virtual const std::vector<double>& coordinate_center() const = 0;

  //: Find the element parametric coordinates of a given point in space.
  //  The return value is the distance from the center coordinate
  //  See the coordinate_center() function.  For the return value d:
  //    0 <=  d < 1 : the point is inside the element
  //    1 ==  d     : the point is on the element surface
  //    1 <   d     : the point is outside the element

  virtual double is_in_element(const double* elem_nodal_coor, // (npe,ndim)
                               const double* point_coor,      // (ndim)
                               double* par_coor) const = 0;   // (ncoord)

  //: interpolate the given nodal field onto the point inside the element,
  //: given the point's parametric coordinates.
  //  Assume the minimum of the master element is (xi,eta,zeta) = (-1,-1,-1);
  //  Assume the maximum of the master element is (xi,eta,zeta) = (+1,+1,+1);

  virtual void interpolate_point(const double* par_coord,   // (npar_coord)
                                 const int& ncomp_field,
                                 const double* field,       // (num_nodes,ncomp_field)
                                 double* result) const = 0; // (ncomp_field)

 protected:
  // All elements have a pointer to a topology object
  stk::topology m_topology;
  std::string m_name;
  unsigned m_integrationOrder{2};

  std::shared_ptr<GaussQuadrature> m_quadrature;
};

}
}

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_MOCKMASTERELEMENT_HPP_ */
