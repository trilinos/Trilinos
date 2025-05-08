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

/** Convert a single mesh object to a single point.
 * This is used in the case where a mesh object is
 * an element with many nodes.  Then something like
 * the element centroid can be used to define
 * a single coordinate point for it.
 * coor is the output coordinate and is assumed to
 * be large enough to hold a single coordinate.
 * Note the maximum needed in all cases is length 3.
 *
 * entity_coordinates is used to return the nodal coordinates
 * that compute_entity_centroid averages.  The return value is
 * the mesh objects from which the coordinates were obtained.
 * compute_entity_centroid returns a vector of vectors containing
 * the coordinates of the nodes that were used to compute the
 * centroid.
 */
std::vector<stk::mesh::Entity>
entity_coordinates(stk::mesh::Entity entity, const stk::mesh::FieldBase* stkField,
                   std::vector<std::vector<double> >& coordinates)
{
  coordinates.clear();
  std::vector<stk::mesh::Entity> meshNodes;
  const unsigned ndim = stkField->mesh_meta_data().spatial_dimension();

  const stk::mesh::BulkData& meshBulk = stkField->get_mesh();

  const stk::mesh::EntityRank fieldType = stkField->entity_rank();
  const stk::mesh::EntityRank objType = meshBulk.entity_rank(entity);
  const bool fieldExistsOnEntity = (fieldType == objType);

  if(objType == stk::topology::NODE_RANK || fieldExistsOnEntity) {
    const double* coor = static_cast<const double *>(stk::mesh::field_data(*stkField, entity));
    if(ndim) {
      STK_ThrowRequireMsg(coor, " Error: The field does not exist, NULL pointer found.\n"
                                << " Field Name:" << stkField->name() << "\n"
                                << " Mesh Object number:" << meshBulk.identifier(entity) << "\n");
      coordinates.emplace_back(coor, coor + ndim);
      meshNodes.push_back(entity);
    }
  }
  else {
    // Loop over node relations in mesh object
    stk::mesh::Entity const* nodes = meshBulk.begin_nodes(entity);
    for(int i = 0, ie = meshBulk.num_nodes(entity); i < ie; ++i) {
      stk::mesh::Entity node = nodes[i];
      double* coor = (double*)stk::mesh::field_data(*stkField, node);

      STK_ThrowRequireMsg(coor, " Error: The field does not exist, NULL pointer found.\n"
                                << " Field Name:" << stkField->name() << "\n"
                                << " Mesh Object number:" << meshBulk.identifier(entity) << "\n");

      coordinates.emplace_back(coor, coor + ndim);
      meshNodes.push_back(node);
    }
  }
  return meshNodes;
}

void
average_coordinates_for_centroid(const std::vector<std::vector<double> >& coordinates, std::vector<double>& centroid)
{
  const int ndim = coordinates.front().size();
  const int num_nodes = coordinates.size();

  centroid.assign(ndim, 0.0);

  for(int j = 0; j < num_nodes; ++j) {
    for(int i = 0; i < ndim; ++i) {
      centroid[i] += coordinates[j][i];
    }
  }
  if(1 != num_nodes) {
    for(int i = 0; i < ndim; ++i) {
      centroid[i] /= num_nodes;
    }
  }
}

std::vector<std::vector<double> >
compute_entity_centroid(stk::mesh::Entity entity, const stk::mesh::FieldBase& nodalCoorRef, std::vector<double>& centroid)
{
  std::vector<std::vector<double> > coordinates;
  entity_coordinates(entity, &nodalCoorRef, coordinates);

  average_coordinates_for_centroid(coordinates, centroid);
  return coordinates;
}


void determine_centroid(const unsigned spatialDimension, stk::mesh::Entity element,
                        const stk::mesh::FieldBase& coordinateField, double* centroid)
{
  const stk::mesh::BulkData& bulk_data = coordinateField.get_mesh();
  const stk::mesh::Entity* const nodes = bulk_data.begin_nodes(element);
  const unsigned numNodes = bulk_data.num_nodes(element);

  for(unsigned i = 0; i < spatialDimension; ++i) {
    centroid[i] = 0.0;
  }

  for(unsigned iNode = 0; iNode < numNodes; ++iNode) {
    stk::mesh::Entity node = nodes[iNode];
    double* coor = static_cast<double*>(stk::mesh::field_data(coordinateField, node));
    STK_ThrowRequireMsg(coor != nullptr, "Entity " << bulk_data.entity_key(node)
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

double distance_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                         const stk::mesh::FieldBase* coordinateField, const std::vector<double>& point)
{
  return distance_from_nearest_entity_node(bulk, entity, coordinateField, point.data());
}

} // end namespace search
} // end namespace stk
