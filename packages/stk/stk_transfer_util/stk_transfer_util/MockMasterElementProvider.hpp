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
#include "stk_transfer_util/Legendre.hpp"
#include "stk_transfer_util/MockMasterElementLine2.hpp"
#include "stk_transfer_util/MockMasterElementQuad4.hpp"
#include "stk_transfer_util/MockMasterElementHex8.hpp"
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
namespace transfer_util {

class MasterElementProvider : public stk::search::MasterElementProviderInterface {
 public:
  MasterElementProvider(const unsigned integrationOrder)
  : m_integrationOrder(MasterElement::get_integration_order(integrationOrder))
  {
    m_line2MasterElement = std::make_shared<MasterElementLine2>(m_integrationOrder);
    m_quad4MasterElement = std::make_shared<MasterElementQuad4>(m_integrationOrder);
    m_hex8MasterElement  = std::make_shared<MasterElementHex8>(m_integrationOrder);
  }

  MasterElementProvider()
  : m_integrationOrder(MasterElement::get_integration_order(0))
  {
    m_line2MasterElement = std::make_shared<MasterElementLine2>(m_integrationOrder);
    m_quad4MasterElement = std::make_shared<MasterElementQuad4>(m_integrationOrder);
    m_hex8MasterElement  = std::make_shared<MasterElementHex8>(m_integrationOrder);
  }

  unsigned num_integration_points(const stk::search::SearchTopology& meTopo) const override;

  void integration_points(const stk::search::SearchTopology& meTopo, std::vector<double>& gaussPoints) const override;

  //: Interpolate the given nodal field onto the point inside the element,
  //: given the point's parametric coordinates.
  //  Assume the minimum of the master element is (xi,eta,zeta) = (-1,-1,-1);
  //  Assume the maximum of the master element is (xi,eta,zeta) = (+1,+1,+1);
  void evaluate_field(const stk::search::SearchTopology& meTopo,
                      const std::vector<double>& paramCoords,
                      const unsigned numFieldComponents,
                      const std::vector<double>& fieldData,
                      std::vector<double>& result) const override;

  void evaluate_field(const stk::search::SearchTopology& meTopo,
                      const unsigned numEvalPoints,
                      const std::vector<double>& paramCoords,
                      const unsigned numFieldComponents,
                      const std::vector<double>& fieldData,
                      std::vector<double>& result) const override;

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
  void nodal_field_data(const stk::search::spmd::EntityKeyPair& key,
                        const stk::search::SearchField& field,
                        unsigned& numFieldComponents,
                        unsigned& numNodes,
                        std::vector<double>& fieldData) const override;

  void nodal_field_data(const std::vector<stk::search::spmd::EntityKeyPair>& nodeKeys,
                        const stk::search::SearchField& field,
                        unsigned& numFieldComponents,
                        std::vector<double>& fieldData) const override;

  //: Find the element parametric coordinates of a given point in space.
  //  Assume the minimum of the master element is (xi,eta,zeta) = (-1,-1,-1);
  //  Assume the maximum of the master element is (xi,eta,zeta) = (+1,+1,+1);
  //  The return value is the distance from the center coordinate
  //  See the coordinate_center() function.  For the return value d:
  //    0 <=  d < 1 : the point is inside the element
  //    1 ==  d     : the point is one the element surface
  //    1 <   d     : the point is outside the element
  void find_parametric_coordinates(const stk::search::SearchTopology& meTopo,
                                   const unsigned numCoordComponents,
                                   const std::vector<double>& elementNodeCoords,
                                   const std::vector<double>& inputCoords,
                                   std::vector<double>& paramCoords,
                                   double& paramDistance) const override;

  unsigned num_parametric_coordinates(const stk::search::SearchTopology& meTopo) const override;

  //: Return the logical center coordinate of the element.
  //  This is not necessarily the centroid.  The center coordinate
  //  is used with the scalar return value for the is_in_element()
  //  function to scale a parametric coordinate to the parametric
  //  surface of the master element.  Given the center coordinate
  //  C and the parametric distance d for the parametric point X then:
  //            Y = (X-C)/d + C
  //  will define a point Y on the element parametric surface.
  void coordinate_center(const stk::search::SearchTopology& meTopo, std::vector<double>& center) const override;

  ~MasterElementProvider() = default;

 private:
  unsigned m_integrationOrder{2};

  std::shared_ptr<MasterElementLine2> m_line2MasterElement;
  std::shared_ptr<MasterElementQuad4> m_quad4MasterElement;
  std::shared_ptr<MasterElementHex8>  m_hex8MasterElement;

  void check_consistent_topology(const stk::search::SearchTopology& meTopo) const;
  const MasterElement* get_master_element(const stk::search::SearchTopology& meTopo) const;
};

}
}

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MOCKSEARCHHEX8MASTERELEMENTPROVIDER_HPP_ */
