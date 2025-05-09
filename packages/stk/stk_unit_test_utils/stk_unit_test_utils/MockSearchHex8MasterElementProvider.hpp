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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MOCKSEARCHHEX8MASTERELEMENTPROVIDER_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MOCKSEARCHHEX8MASTERELEMENTPROVIDER_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "Legendre.hpp"
#include "MockSearchHex8Element.hpp"
#include "stk_search/Box.hpp"                         // for Box
#include "stk_search/Point.hpp"                       // for Point
#include "stk_search/Sphere.hpp"                      // for Sphere
#include "stk_util/parallel/Parallel.hpp"             // for ParallelMachine
#include "stk_mesh/base/Types.hpp"                    // for EntityId, PartV...
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_search/IdentProc.hpp"                   // for IdentProc
#include "stk_topology/topology.hpp"                  // for topology
#include "stk_search_util/MasterElementProvider.hpp"
#include <memory>                                     // for shared_ptr
#include <set>                                        // for set
#include <string>                                     // for string, basic_s...
#include <utility>                                    // for pair
#include <vector>                                     // for vector
#include <cassert>                   // for assert
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

/**
 * A 3D Gauss-Legendre quadrature rule (traditionally
 * called the Gauss quadrature) of arbitrary order q x q x q, on the
 * interval [-1,1] x [-1,1] x [-1,1].
 *
 */
class Hex8GaussQuadrature {
public:
  Hex8GaussQuadrature(unsigned q)
  : m_order(q)
  , m_numIntgPoints(q * q * q)
  , m_numParametricCoordinates(3)
  , m_intgWeights(m_numIntgPoints)
  , m_intgLocations(static_cast<size_t>(m_numParametricCoordinates) * m_numIntgPoints)
  {
    // initialize the points and weights
    gauss_legendre(q, m_intgLocations, m_intgWeights);
  }

  ~Hex8GaussQuadrature() = default;

  unsigned order() const { return m_order; }
  unsigned num_intg_pts() const { return m_numIntgPoints; }
  unsigned num_parametric_coordinates() const { return m_numParametricCoordinates; }
  const std::vector<double>& intg_pt_weights() const { return m_intgWeights; }
  const std::vector<double>& intg_pt_locations() const { return m_intgLocations; }

private:
  unsigned m_order;
  unsigned m_numIntgPoints; // (num_intg_pts)
  unsigned m_numParametricCoordinates;

  std::vector<double> m_intgWeights; // array holding values of the weights
  std::vector<double> m_intgLocations; // array holding locations of the pts
};

class Hex8MasterElementProvider : public stk::search::MasterElementProviderInterface {
 public:
  Hex8MasterElementProvider(const unsigned integrationOrder)
  : m_integrationOrder(get_integration_order(integrationOrder))
  , m_quadrature(m_integrationOrder)
  {

  }

  Hex8MasterElementProvider()
  : m_integrationOrder(get_integration_order(0))
  , m_quadrature(m_integrationOrder)
  {

  }

  unsigned num_integration_points(const stk::search::MasterElementTopology& meTopo) override;

  void integration_points(const stk::search::MasterElementTopology& meTopo, std::vector<double>& gaussPoints) override;

  //: Interpolate the given nodal field onto the point inside the element,
  //: given the point's parametric coordinates.
  //  Assume the minimum of the master element is (xi,eta,zeta) = (-1,-1,-1);
  //  Assume the maximum of the master element is (xi,eta,zeta) = (+1,+1,+1);
  void evaluate_field(const stk::search::MasterElementTopology& meTopo,
                      const std::vector<double>& paramCoords,
                      const unsigned numFieldComponents,
                      const std::vector<double>& fieldData,
                      std::vector<double>& result) override;

  /*--------------------------------------------------------------------*/
  /* Arrays are contiguous blocks of memory laid out in fortran order,  */
  /* i.e. first index cycles fastest and last index cycles slowest when */
  /* array is accessed in memory order i.e transposed C order           */
  /*                                                                    */
  /* Dimensions:  ncoord number of parametric coordinates per           */
  /*                     interpolation location                         */
  /*              npts   number of locations to interpolate per element */
  /*              ncomp  number of components in interpolated field     */
  /*                            npe        number nodes per element     */
  /*              nelem  number of elements                             */
  /*--------------------------------------------------------------------*/
  void nodal_field_data(const stk::mesh::Entity entity,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        unsigned& numNodes,
                        std::vector<double>& fieldData) override;

  void nodal_field_data(const stk::mesh::EntityKey key,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        unsigned& numNodes,
                        std::vector<double>& fieldData) override;

  void nodal_field_data(const std::vector<stk::mesh::Entity>& nodes,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        std::vector<double>& fieldData) override;

  void nodal_field_data(const std::vector<stk::mesh::EntityKey>& nodeKeys,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        std::vector<double>& fieldData) override;

  //: Find the element parametric coordinates of a given point in space.
  //  Assume the minimum of the master element is (xi,eta,zeta) = (-1,-1,-1);
  //  Assume the maximum of the master element is (xi,eta,zeta) = (+1,+1,+1);
  //  The return value is the distance from the center coordinate
  //  See the coordinate_center() function.  For the return value d:
  //    0 <=  d < 1 : the point is inside the element
  //    1 ==  d     : the point is one the element surface
  //    1 <   d     : the point is outside the element
  void find_parametric_coordinates(const stk::search::MasterElementTopology& meTopo,
                                   const unsigned numCoordComponents,
                                   const std::vector<double>& elementNodeCoords,
                                   const std::vector<double>& inputCoords,
                                   std::vector<double>& paramCoords,
                                   double& paramDistance) override;

  unsigned num_parametric_coordinates(const stk::search::MasterElementTopology& meTopo) override;

  //: Return the logical center coordinate of the element.
  //  This is not necessarily the centroid.  The center coordinate
  //  is used with the scalar return value for the is_in_element()
  //  function to scale a parametric coordinate to the parametric
  //  surface of the master element.  Given the center coordinate
  //  C and the parametric distance d for the parametric point X then:
  //            Y = (X-C)/d + C
  //  will define a point Y on the element parametric surface.
  void coordinate_center(const stk::search::MasterElementTopology& meTopo, std::vector<double>& center) override;

  ~Hex8MasterElementProvider() = default;

 private:
  unsigned m_integrationOrder{2};
  Hex8GaussQuadrature m_quadrature;

  void check_consistent_topology(const stk::search::MasterElementTopology& meTopo);

  unsigned get_integration_order(const unsigned integrationOrder);
};

}
}

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MOCKSEARCHHEX8MASTERELEMENTPROVIDER_HPP_ */
