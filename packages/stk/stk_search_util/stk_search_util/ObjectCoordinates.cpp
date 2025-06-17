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

#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_search/DistanceComparison.hpp"          // for distance_sq

namespace stk {
namespace search {

namespace impl {
void determine_node_centroid(const unsigned spatialDimension, stk::mesh::Entity node,
                             const stk::mesh::FieldBase& coordinateField, double* centroid)
{
  const stk::mesh::BulkData& bulkData = coordinateField.get_mesh();

  STK_ThrowRequireMsg(bulkData.entity_rank(node) == stk::topology::NODE_RANK, "Input entity must be a node");

  double* coor = static_cast<double*>(stk::mesh::field_data(coordinateField, node));
  STK_ThrowRequireMsg(coor != nullptr, "Input node " << bulkData.entity_key(node)
                                                     << " has no data for coordinate field: "
                                                     << coordinateField.name());
  for(unsigned i = 0; i < spatialDimension; ++i) {
    centroid[i] = coor[i];
  }
}

}

void determine_centroid(const unsigned spatialDimension, stk::mesh::Entity entity,
                        const stk::mesh::FieldBase& coordinateField, double* centroid)
{
  const stk::mesh::BulkData& bulkData = coordinateField.get_mesh();

  if(bulkData.entity_rank(entity) == stk::topology::NODE_RANK) {
    impl::determine_node_centroid(spatialDimension, entity, coordinateField, centroid);
    return;
  }

  const stk::mesh::Entity* const nodes = bulkData.begin_nodes(entity);
  const unsigned numNodes = bulkData.num_nodes(entity);
  STK_ThrowAssertMsg(numNodes != 0, "Input entity " << bulkData.entity_key(entity) << " has no connected nodes");

  for(unsigned i = 0; i < spatialDimension; ++i) {
    centroid[i] = 0.0;
  }

  for(unsigned iNode = 0; iNode < numNodes; ++iNode) {
    stk::mesh::Entity node = nodes[iNode];
    double* coor = static_cast<double*>(stk::mesh::field_data(coordinateField, node));
    STK_ThrowRequireMsg(coor != nullptr, "Node " << bulkData.entity_key(node)
                                                 << " connected to input entity "
                                                 << bulkData.entity_key(entity)
                                                 << " has no data for coordinate field: " << coordinateField.name());
    for(unsigned i = 0; i < spatialDimension; ++i) {
      centroid[i] += coor[i];
    }
  }
  for(unsigned i = 0; i < spatialDimension; ++i) {
    centroid[i] /= numNodes;
  }
}

void determine_centroid(const unsigned spatialDimension, stk::mesh::Entity element,
                        const stk::mesh::FieldBase& coordinateField, std::vector<double>& centroid)
{
  centroid.clear();
  centroid.resize(spatialDimension);

  determine_centroid(spatialDimension, element, coordinateField, centroid.data());
}

double distance_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                         const stk::mesh::FieldBase* coordinateField, const double* point)
{
  STK_ThrowRequireMsg(bulk.entity_rank(entity) > stk::topology::NODE_RANK,
                      "Invalid entity rank for object: " << bulk.entity_rank(entity));

  double minDistance = std::numeric_limits<double>::max();
  const unsigned nDim = coordinateField->mesh_meta_data().spatial_dimension();

  if(coordinateField->entity_rank() == stk::topology::NODE_RANK) {
    const stk::mesh::Entity* const nodes = bulk.begin_nodes(entity);
    const unsigned numNodes = bulk.num_nodes(entity);

    for(unsigned i = 0; i < numNodes; ++i) {
      double distance = 0.0;
      double* coordinates = static_cast<double*>(stk::mesh::field_data(*coordinateField, nodes[i]));

      for(unsigned j = 0; j < nDim; ++j) {
        const double t = point[j] - coordinates[j];
        distance += t * t;
      }
      if(distance < minDistance) minDistance = distance;
    }
    minDistance = std::sqrt(minDistance);
  }
  else if(coordinateField->entity_rank() == stk::topology::ELEM_RANK) {
    const double* coor = static_cast<const double*>(stk::mesh::field_data(*coordinateField, entity));
    minDistance = distance(nDim, coor, point);
  }

  return minDistance;
}

double distance_squared_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                                 const stk::mesh::FieldBase* coordinateField, const double* point)
{
  STK_ThrowRequireMsg(bulk.entity_rank(entity) > stk::topology::NODE_RANK,
                      "Invalid entity rank for object: " << bulk.entity_rank(entity));

  double minDistance = std::numeric_limits<double>::max();
  const unsigned nDim = coordinateField->mesh_meta_data().spatial_dimension();

  if(coordinateField->entity_rank() == stk::topology::NODE_RANK) {
    const stk::mesh::Entity* const nodes = bulk.begin_nodes(entity);
    const unsigned numNodes = bulk.num_nodes(entity);

    for(unsigned i = 0; i < numNodes; ++i) {
      double distance = 0.0;
      double* coordinates = static_cast<double*>(stk::mesh::field_data(*coordinateField, nodes[i]));

      for(unsigned j = 0; j < nDim; ++j) {
        const double t = point[j] - coordinates[j];
        distance += t * t;
      }
      if(distance < minDistance) minDistance = distance;
    }
  }
  else if(coordinateField->entity_rank() == stk::topology::ELEM_RANK) {
    const double* coor = static_cast<const double*>(stk::mesh::field_data(*coordinateField, entity));
    minDistance = distance_sq(nDim, coor, point);
  }

  return minDistance;
}

void determine_gauss_points(const stk::mesh::BulkData& recvBulk, stk::mesh::Entity element,
                            const MasterElementProviderInterface& masterElemProvider,
                            const stk::mesh::FieldBase& coordinateField, std::vector<double>& location)
{
  std::vector<double> gaussPoints;
  std::vector<double> fieldData;

  unsigned numFieldComponents;
  unsigned numNodes;

  SearchField masterElementCoordField(&coordinateField);
  SearchTopology topo(recvBulk.bucket(element).topology(), spmd::make_entity_key_pair(recvBulk,element));

  unsigned numGaussPoints = masterElemProvider.num_integration_points(topo);

  masterElemProvider.integration_points(topo, gaussPoints);
  masterElemProvider.nodal_field_data(topo.get_key(), masterElementCoordField, numFieldComponents, numNodes, fieldData);

  masterElemProvider.evaluate_field(topo, numGaussPoints, gaussPoints, numFieldComponents, fieldData, location);
}

double distance_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                         const stk::mesh::FieldBase* coordinateField, const std::vector<double>& point)
{
  return distance_from_nearest_entity_node(bulk, entity, coordinateField, point.data());
}

double distance_squared_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                                 const stk::mesh::FieldBase* coordinateField, const std::vector<double>& point)
{
  return distance_squared_from_nearest_entity_node(bulk, entity, coordinateField, point.data());
}

} // end namespace search
} // end namespace stk
