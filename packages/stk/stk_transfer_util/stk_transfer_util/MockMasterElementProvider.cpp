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
#include "stk_transfer_util/MockMasterElementProvider.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer_util {

  unsigned MasterElementProvider::num_integration_points(const stk::search::SearchTopology& meTopo) const
  {
    const MasterElement* me = get_master_element(meTopo);
    return me->num_intg_pts();
  }

  void MasterElementProvider::integration_points(const stk::search::SearchTopology& meTopo, std::vector<double>& gaussPoints) const
  {
    const MasterElement* me = get_master_element(meTopo);
    me->intg_pt_locations(gaussPoints);
  }

  void MasterElementProvider::evaluate_field(const stk::search::SearchTopology& meTopo,
                                             const unsigned numEvalPoints,
                                             const std::vector<double>& paramCoords,
                                             const unsigned numFieldComponents,
                                             const std::vector<double>& fieldData,
                                             std::vector<double>& result) const
  {
    check_consistent_topology(meTopo);
    const stk::topology topo = meTopo.get_topology();
    const unsigned numParCoords = num_parametric_coordinates(meTopo);
    STK_ThrowRequireMsg(paramCoords.size() >= numParCoords*numEvalPoints,
                        "Insufficient length for parametric coordinates: " << paramCoords.size() << " Expected: " << numParCoords*numEvalPoints);

    if(numFieldComponents == 0 || numEvalPoints == 0) return;

    result.clear();
    result.assign(numFieldComponents*numEvalPoints, 0.0);

    const unsigned numNodes = topo.num_nodes();
    STK_ThrowRequireMsg(fieldData.size() >= numNodes*numFieldComponents,
                        "Insufficient length for fieldData: " << fieldData.size() << " Expected: " << numNodes*numFieldComponents);

    const MasterElement* me = get_master_element(meTopo);
    for(unsigned i=0; i<numEvalPoints; i++) {
      me->interpolate_point(&paramCoords[i*numParCoords], numFieldComponents, fieldData.data(), &result[i*numFieldComponents]);
    }
  }

  void MasterElementProvider::evaluate_field(const stk::search::SearchTopology& meTopo,
                                             const std::vector<double>& paramCoords,
                                             const unsigned numFieldComponents,
                                             const std::vector<double>& fieldData,
                                             std::vector<double>& result) const
  {
    evaluate_field(meTopo, 1, paramCoords, numFieldComponents, fieldData, result);
  }

  void MasterElementProvider::nodal_field_data(const stk::search::spmd::EntityKeyPair& key,
                                               const stk::search::SearchField& meField,
                                               unsigned& numFieldComponents,
                                               unsigned& numNodes,
                                               std::vector<double>& fieldData) const
  {
    const stk::mesh::Entity entity = key;

    // Extract transposed field data
    const stk::mesh::FieldBase& field = *meField.get_field();
    const stk::mesh::BulkData& bulk = field.get_mesh();

    // load nodal coordinates from entity
    stk::mesh::Entity const* nodes = bulk.begin_nodes(entity);

    numNodes = bulk.num_nodes(entity);
    numFieldComponents = numNodes > 0 ? stk::mesh::field_extent0_per_entity(field, nodes[0]) : 0;

    fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(field,
      [&](auto& meFieldData) {
        for (unsigned ni = 0u; ni < numNodes; ++ni) {
          stk::mesh::Entity node = nodes[ni];

          if (field.defined_on(node)) {
            auto data = meFieldData.entity_values(node);
            for (stk::mesh::ComponentIdx j(0); j < static_cast<int>(numFieldComponents); ++j) {
              const auto offSet = ni + j * numNodes;
              fieldData[offSet] = meField.transform(data(j));
            }
          }
          else {
            for (unsigned j = 0; j < numFieldComponents; ++j) {
              const auto offSet = ni + j * numNodes;
              fieldData[offSet] = meField.default_value();
            }
          }
        }
      }
    );
  }

  void MasterElementProvider::nodal_field_data(const std::vector<stk::search::spmd::EntityKeyPair>& nodeKeys,
                                               const stk::search::SearchField& meField,
                                               unsigned& numFieldComponents,
                                               std::vector<double>& fieldData) const
  {
    // Extract transposed field data
    const stk::mesh::FieldBase& field = *meField.get_field();
    unsigned numNodes = nodeKeys.size();
    numFieldComponents = numNodes > 0 ? stk::mesh::field_extent0_per_entity(field, nodeKeys[0]) : 0;
    fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(field,
      [&](auto& meFieldData) {
        for (unsigned ni = 0u; ni < numNodes; ++ni) {
          stk::mesh::Entity node = nodeKeys[ni];

          unsigned numNodeFieldComponents = stk::mesh::field_extent0_per_entity(field, node);

          if (field.defined_on(node)) {
            auto data = meFieldData.entity_values(node);
            for (stk::mesh::ComponentIdx j = 0_comp; j < static_cast<int>(numNodeFieldComponents); ++j) {
              const auto offSet = ni + j * numNodes;
              fieldData[offSet] = meField.transform(data(j));
            }
          }
          else {
            for (unsigned j = 0; j < numNodeFieldComponents; ++j) {
              const auto offSet = ni + j * numNodes;
              fieldData[offSet] = meField.default_value();
            }
          }
        }
      }
    );
  }

  void MasterElementProvider::find_parametric_coordinates(const stk::search::SearchTopology& meTopo,
                                                              const unsigned numCoordComponents,
                                                              const std::vector<double>& elementNodeCoords,
                                                              const std::vector<double>& inputCoords,
                                                              std::vector<double>& paramCoords,
                                                              double& paramDistance) const
  {
    const MasterElement* me = get_master_element(meTopo);
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

    paramDistance = me->is_in_element(elementNodeCoords.data(), inputCoords.data(), paramCoords.data());
  }

  unsigned MasterElementProvider::num_parametric_coordinates(const stk::search::SearchTopology& meTopo) const
  {
    const MasterElement* me = get_master_element(meTopo);
    return me->num_parametric_coordinates();
  }

  void MasterElementProvider::coordinate_center(const stk::search::SearchTopology& meTopo, std::vector<double>& center) const
  {
    const MasterElement* me = get_master_element(meTopo);
    center = me->coordinate_center();
  }

  const MasterElement* MasterElementProvider::get_master_element(const stk::search::SearchTopology& meTopo) const
  {
    check_consistent_topology(meTopo);

    stk::topology topology = meTopo.get_topology();

    if(topology == stk::topology::HEX_8) {
      return m_hex8MasterElement.get();
    }
    if(topology == stk::topology::QUAD_4) {
      return m_quad4MasterElement.get();
    }
    if(topology == stk::topology::LINE_2) {
      return m_line2MasterElement.get();
    }

    return nullptr;
  }

  void MasterElementProvider::check_consistent_topology([[maybe_unused]] const stk::search::SearchTopology& meTopo) const
  {
    // Only support Hex8, Quad4 and Line2
    STK_ThrowRequireMsg(meTopo.get_topology() == stk::topology::HEX_8 ||
                       meTopo.get_topology() == stk::topology::QUAD_4 ||
                       meTopo.get_topology() == stk::topology::LINE_2,
                       "Invalid topology " << meTopo.get_topology());
  }
}
}

