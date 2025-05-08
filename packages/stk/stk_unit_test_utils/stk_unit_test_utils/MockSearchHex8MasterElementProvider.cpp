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
#include "MockSearchHex8MasterElementProvider.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

  unsigned Hex8MasterElementProvider::num_integration_points(const stk::search::MasterElementTopology& meTopo)
  {
    check_consistent_topology(meTopo);
    return m_quadrature.num_intg_pts();
  }

  void Hex8MasterElementProvider::integration_points(const stk::search::MasterElementTopology& meTopo, std::vector<double>& gaussPoints)
  {
    check_consistent_topology(meTopo);
    const std::vector<double>& gp = m_quadrature.intg_pt_locations();
    gaussPoints.assign(gp.begin(), gp.end());
  }

  void Hex8MasterElementProvider::evaluate_field(const stk::search::MasterElementTopology& meTopo,
                      const std::vector<double>& paramCoords,
                      const unsigned numFieldComponents,
                      const std::vector<double>& fieldData,
                      std::vector<double>& result)
  {
    check_consistent_topology(meTopo);
    const stk::topology topo = meTopo.get_topology();
    const unsigned numParCoords = num_parametric_coordinates(meTopo);
    STK_ThrowRequireMsg(paramCoords.size() >= numParCoords,
                        "Insufficient length for parametric coordinates: " << paramCoords.size() << " Expected: " << numParCoords);

    if(numFieldComponents == 0) return;

    result.clear();
    result.assign(numFieldComponents, 0.0);

    const unsigned numNodes = topo.num_nodes();
    STK_ThrowRequireMsg(fieldData.size() >= numNodes*numFieldComponents,
                        "Insufficient length for fieldData: " << fieldData.size() << " Expected: " << numNodes*numFieldComponents);

    Hex8::interpolate_point(numParCoords, paramCoords.data(), numFieldComponents, fieldData.data(), result.data());
  }

  void Hex8MasterElementProvider::nodal_field_data(const stk::mesh::Entity entity,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        unsigned& numNodes,
                        std::vector<double>& fieldData)
  {
    // Extract transposed field data
    const stk::mesh::BulkData& bulk = field.get_mesh();

    // load nodal coordinates from entity
    stk::mesh::Entity const* nodes = bulk.begin_nodes(entity);

    numNodes = bulk.num_nodes(entity);
    if(numNodes > 0) {
      numFieldComponents = stk::mesh::field_extent0_per_entity(field, nodes[0]);
    } else {
      numFieldComponents = 0;
    }

    fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

    for(unsigned ni = 0u; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      const double* data = static_cast<const double *>(stk::mesh::field_data(field, node));
      for(unsigned j = 0; j < numFieldComponents; ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = data[j];
      }
    }
  }

  void Hex8MasterElementProvider::nodal_field_data(const stk::mesh::EntityKey key,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        unsigned& numNodes,
                        std::vector<double>& fieldData)
  {
    const stk::mesh::Entity entity = field.get_mesh().get_entity(key);
    nodal_field_data(entity, field, numFieldComponents, numNodes, fieldData);
  }

  void Hex8MasterElementProvider::nodal_field_data(const std::vector<stk::mesh::Entity>& nodes,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        std::vector<double>& fieldData)
  {
    // Extract transposed field data
    unsigned numNodes = nodes.size();
    numFieldComponents = numNodes > 0 ? stk::mesh::field_extent0_per_entity(field, nodes[0]) : 0;
    fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

    for(unsigned ni = 0u; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      unsigned numNodeFieldComponents = stk::mesh::field_extent0_per_entity(field, node);

      const double* data = static_cast<const double *>(stk::mesh::field_data(field, node));
      for(unsigned j = 0; j < numNodeFieldComponents; ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = data[j];
      }
    }
  }

  void Hex8MasterElementProvider::nodal_field_data(const std::vector<stk::mesh::EntityKey>& nodeKeys,
                        const stk::mesh::FieldBase& field,
                        unsigned& numFieldComponents,
                        std::vector<double>& fieldData)
  {
    stk::mesh::EntityVector nodes = stk::mesh::impl::convert_keys_to_entities(field.get_mesh(), nodeKeys);
    nodal_field_data(nodes, field, numFieldComponents, fieldData);
  }

  void Hex8MasterElementProvider::find_parametric_coordinates(const stk::search::MasterElementTopology& meTopo,
                                   const unsigned numCoordComponents,
                                   const std::vector<double>& elementNodeCoords,
                                   const std::vector<double>& inputCoords,
                                   std::vector<double>& paramCoords,
                                   double& paramDistance)
  {
    check_consistent_topology(meTopo);
    const stk::topology topo = meTopo.get_topology();
    const unsigned nDim = numCoordComponents;
    STK_ThrowRequireMsg(inputCoords.size() >= nDim,
                        "Insufficient length for input coordinates: " << inputCoords.size() << " Expected: " << nDim);

    const unsigned long numNodes = topo.num_nodes();
    STK_ThrowRequireMsg(elementNodeCoords.size() >= numNodes*nDim,
                        "Insufficient length for elementNodeCoords: " << elementNodeCoords.size() << " Expected: " <<
                        numNodes*nDim);

    const unsigned numParCoords = num_parametric_coordinates(meTopo);
    paramCoords.clear();
    paramCoords.assign(numParCoords, std::numeric_limits<double>::max());

    paramDistance = Hex8::is_in_element(elementNodeCoords.data(), inputCoords.data(), paramCoords.data());
  }

  unsigned Hex8MasterElementProvider::num_parametric_coordinates(const stk::search::MasterElementTopology& meTopo)
  {
    check_consistent_topology(meTopo);
    return m_quadrature.num_parametric_coordinates();
  }

  void Hex8MasterElementProvider::coordinate_center(const stk::search::MasterElementTopology& meTopo, std::vector<double>& center)
  {
    check_consistent_topology(meTopo);
    center = Hex8::coordinate_center();
  }

  void Hex8MasterElementProvider::check_consistent_topology(const stk::search::MasterElementTopology& meTopo)
  {
    // Only support Hex8
    const stk::topology topo = meTopo.get_topology();
    STK_ThrowRequireMsg(topo == stk::topology::HEX_8, "Invalid topology " << topo);
  }

  unsigned Hex8MasterElementProvider::get_integration_order(const unsigned integrationOrder)
  {
    return integrationOrder == 0 ? 2 : integrationOrder;
  }

}
}

