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
#include "stk_search_util/PointProjection.hpp"
#include "stk_search_util/MeshUtility.hpp"

#include "stk_mesh/base/Bucket.hpp"                   // for Bucket
#include "stk_mesh/base/BulkData.hpp"                 // for BulkData
#include "stk_mesh/base/CompositeRank.hpp"            // for CompositeRank
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/FieldBase.hpp"                // for field_data, Fie...
#include "stk_mesh/base/MetaData.hpp"                 // for MetaData
#include "stk_mesh/base/Part.hpp"                     // for Part
#include "stk_util/util/ReportHandler.hpp"            // for eval_test_condi...
#include "stk_util/parallel/ParallelReduce.hpp"       // for all_reduce_max

#include <algorithm>                                  // for all_of, max, min
#include <array>                                      // for array
#include <cmath>                                      // for sqrt
#include <cstddef>                                    // for size_t
#include <limits>                                     // for numeric_limits
#include <memory>                                     // for __shared_ptr_ac...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

namespace impl {

void gather_nodal_coordinates(const ProjectionData& data, stk::mesh::Entity entity, std::vector<double>& coords)
{
  unsigned numFieldComponents;
  unsigned numNodes;

  spmd::EntityKeyPair key(entity, data.bulk.entity_key(entity));
  data.masterElemProvider->nodal_field_data(key, data.projectedCoordinateField, numFieldComponents, numNodes, coords);

  const unsigned spatialDimension = data.bulk.mesh_meta_data().spatial_dimension();
  STK_ThrowRequireMsg(numFieldComponents == spatialDimension, "Invalid coordinate field: " << data.coordinateField.name());
}

void gather_nodal_coordinates(const ProjectionData& data, const spmd::EntityKeyPair& key, std::vector<double>& coords)
{
  unsigned numFieldComponents;
  unsigned numNodes;

  data.masterElemProvider->nodal_field_data(key, data.projectedCoordinateField, numFieldComponents, numNodes, coords);

  const unsigned spatialDimension = data.bulk.mesh_meta_data().spatial_dimension();
  STK_ThrowRequireMsg(numFieldComponents == spatialDimension, "Invalid coordinate field: " << data.coordinateField.name());
}

void gather_nodal_coordinates(const ProjectionData& data, stk::mesh::EntityVector& sideNodes, std::vector<double>& coords)
{
  unsigned numFieldComponents;
  unsigned numNodes = sideNodes.size();

  data.scratchKeySpace.resize(numNodes);
  for(unsigned i=0; i<numNodes; i++) {
    data.scratchKeySpace[i] = spmd::EntityKeyPair(sideNodes[i], data.bulk.entity_key(sideNodes[i]));
  }
  data.masterElemProvider->nodal_field_data(data.scratchKeySpace, data.projectedCoordinateField, numFieldComponents, coords);

  const unsigned spatialDimension = data.bulk.mesh_meta_data().spatial_dimension();
  STK_ThrowRequireMsg(numFieldComponents == spatialDimension, "Invalid coordinate field: " << data.coordinateField.name());
}

bool fill_projected_side_location(const ProjectionData& data, stk::mesh::Entity side,
                                  stk::mesh::EntityVector& sideNodes, stk::topology sideTopology,
                                  std::vector<double>& evaluatedLocation)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const spmd::EntityKeyPair key(side, bulk.entity_key(side));
  const unsigned spatialDimension = bulk.mesh_meta_data().spatial_dimension();

  const unsigned numParametricCoordinates = get_number_of_parametric_coordinates(sideTopology);
  std::vector<double> parametricCoordinates(numParametricCoordinates);
  std::vector<double> sideCoordinates;

  gather_nodal_coordinates(data, sideNodes, sideCoordinates);

  double parametricDistance;
  auto meTopo = SearchTopology(sideTopology, key, bulk.bucket_ptr(side));
  data.masterElemProvider->find_parametric_coordinates(meTopo, spatialDimension, sideCoordinates,
                                                       data.evalPoint, parametricCoordinates, parametricDistance);

  if(parametricDistance == std::numeric_limits<double>::max()) {
   // Convergence issues
    evaluatedLocation = data.evalPoint;
    return false;
  }

  data.masterElemProvider->evaluate_field(meTopo, parametricCoordinates, spatialDimension, sideCoordinates, evaluatedLocation);

  data.masterElemProvider->find_parametric_coordinates(meTopo, spatialDimension, sideCoordinates,
                                                       evaluatedLocation, parametricCoordinates, parametricDistance);

  if(1 < parametricDistance) {
    std::vector<double> center;
    data.masterElemProvider->coordinate_center(meTopo, center);

    for(size_t j = 0; j < numParametricCoordinates; ++j) {
      parametricCoordinates[j] = ((parametricCoordinates[j] - center[j]) / parametricDistance) + center[j];
    }

    data.masterElemProvider->evaluate_field(meTopo, parametricCoordinates, spatialDimension, sideCoordinates, evaluatedLocation);
  }

  return true;
}

bool fill_projected_side_location(const ProjectionData& data, stk::mesh::Entity side, std::vector<double>& evaluatedLocation)
{
  const stk::mesh::BulkData& send_mesh = data.bulk;
  stk::topology sideTopology = send_mesh.bucket(side).topology();
  stk::mesh::EntityVector sideNodes(send_mesh.begin_nodes(side), send_mesh.end_nodes(side));
  return fill_projected_side_location(data, side, sideNodes, sideTopology, evaluatedLocation);
}

stk::mesh::EntityVector
get_nodes_from_side_ordinal(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                            stk::mesh::ConnectivityOrdinal ordinal, stk::mesh::EntityRank subRank)
{
  stk::topology objTopology = bulk.bucket(entity).topology();
  stk::topology subTopology = objTopology.sub_topology(subRank, ordinal);

  stk::mesh::EntityVector nodes(subTopology.num_nodes());
  objTopology.sub_topology_nodes(bulk.begin_nodes(entity), subRank, ordinal, nodes.data());
  return nodes;
}
}

double fill_projected_object_location(const ProjectionData& data, stk::mesh::Entity entity,
                                      const std::vector<double>& location,
                                      std::vector<double>& evaluatedParametricCoordinates)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const unsigned spatialDimension = bulk.mesh_meta_data().spatial_dimension();

  stk::mesh::Bucket& bucket = bulk.bucket(entity);
  stk::topology objTopology = bucket.topology();
  const spmd::EntityKeyPair key(entity, bulk.entity_key(entity));
  std::vector<double> coordinates;

  impl::gather_nodal_coordinates(data, entity, coordinates);

  double parametricDistance;
  auto meTopo = SearchTopology(objTopology, key, &bucket);
  data.masterElemProvider->find_parametric_coordinates(meTopo, spatialDimension, coordinates,
                                                       location, evaluatedParametricCoordinates, parametricDistance);

  return parametricDistance;
}

void project_to_entity(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const int spatialDimension = meta.spatial_dimension();
  std::vector<double> location(spatialDimension, 0);
  std::vector<double> parametricCoordinates;

  impl::fill_projected_side_location(data, entity, location);
  double parametricDistance = fill_projected_object_location(data, entity, location, parametricCoordinates);

  const double projectedDistanceSquared =
      distance_sq(data.evalPoint.size(), location.data(), data.evalPoint.data());

  result.parametricCoords.swap(parametricCoordinates);
  result.geometricDistanceSquared = projectedDistanceSquared;
  result.parametricDistance = parametricDistance;
  result.doneProjection = true;
}

void project_to_closest_side(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityRank entityRank = bulk.entity_rank(entity);
  stk::topology objTopology = bulk.bucket(entity).topology();
  stk::mesh::EntityRank subTopologyRank =
      (entityRank == stk::topology::ELEM_RANK) ? meta.side_rank() : stk::topology::EDGE_RANK;

  STK_ThrowRequireMsg(entityRank == stk::topology::ELEM_RANK || entityRank == meta.side_rank(),
                  "Send object rank: " << entityRank);

  result.geometricDistanceSquared = std::numeric_limits<double>::max();
  result.parametricDistance = std::numeric_limits<double>::max();
  result.doneProjection = false;

  if(objTopology.num_sub_topology(subTopologyRank) == 0) {
    project_to_entity(data, entity, result);
    return;
  }

  double closestProjectedDistanceSquared = std::numeric_limits<double>::max();

  const int spatialDimension = meta.spatial_dimension();
  stk::mesh::EntityVector sideNodes;

  std::vector<double> location(spatialDimension, 0);
  std::vector<double> parametricCoordinates;

  for(unsigned i = 0; i < objTopology.num_sub_topology(subTopologyRank); i++) {
    stk::mesh::ConnectivityOrdinal ordinal = i;
    sideNodes = impl::get_nodes_from_side_ordinal(bulk, entity, ordinal, subTopologyRank);
    stk::topology sideTopology = objTopology.sub_topology(subTopologyRank, ordinal);

    if(impl::fill_projected_side_location(data, stk::mesh::Entity(), sideNodes, sideTopology, location)) {
      double distance = fill_projected_object_location(data, entity, location, parametricCoordinates);

      const double projectedDistanceSquared = distance_sq(data.evalPoint.size(),
                                                          location.data(),
                                                          data.evalPoint.data());
      if(projectedDistanceSquared < closestProjectedDistanceSquared) {
        closestProjectedDistanceSquared = projectedDistanceSquared;
        result.parametricCoords.swap(parametricCoordinates);
        result.geometricDistanceSquared = closestProjectedDistanceSquared;
        result.parametricDistance = distance;
        result.doneProjection = true;
      }
    }
  }
}

void project_to_closest_face(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result)
{
  double closestProjectedDistanceSquared = std::numeric_limits<double>::max();
  const stk::mesh::BulkData& bulk = data.bulk;
  const int spatialDimension = bulk.mesh_meta_data().spatial_dimension();
  const stk::mesh::Entity* faces = bulk.num_faces(entity) ? bulk.begin_faces(entity) : bulk.begin_edges(entity);
  const unsigned numFaces = bulk.num_faces(entity) ? bulk.num_faces(entity) : bulk.num_edges(entity);
  std::vector<double> location(spatialDimension, 0);
  std::vector<double> parametricCoordinates;
  result.doneProjection = false;

  for(unsigned i = 0; i < numFaces; i++) {
    impl::fill_projected_side_location(data, faces[i], location);

    double parametricDistance = fill_projected_object_location(data, entity, location, parametricCoordinates);

    const double projectedDistanceSquared = distance_sq(data.evalPoint.size(),
                                                        location.data(),
                                                        data.evalPoint.data());
    if(projectedDistanceSquared < closestProjectedDistanceSquared) {
      closestProjectedDistanceSquared = projectedDistanceSquared;
      result.parametricCoords.swap(parametricCoordinates);
      result.geometricDistanceSquared = closestProjectedDistanceSquared;
      result.parametricDistance = parametricDistance;
      result.doneProjection = true;
    }
  }
}

void truncate_to_entity(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const int spatialDimension = bulk.mesh_meta_data().spatial_dimension();
  std::vector<double> location(spatialDimension, 0);
  std::vector<double> parametricCoordinates;
  result.doneProjection = false;

  impl::fill_projected_side_location(data, entity, location);

  double parametricDistance = fill_projected_object_location(data, entity, location, parametricCoordinates);

  const double truncatedDistanceSquared = distance_sq(data.evalPoint.size(),
                                                      location.data(),
                                                      data.evalPoint.data());
  result.parametricCoords.swap(parametricCoordinates);
  result.geometricDistanceSquared = truncatedDistanceSquared;
  result.parametricDistance = parametricDistance;
  result.doneProjection = true;
}

double compute_parametric_distance(const ProjectionData& data, stk::mesh::Entity entity,
                                   const std::vector<double>& parametricCoordinates,
                                   std::vector<double>& evaluatedLocation)
{
  const stk::mesh::BulkData& bulk = data.bulk;
  const unsigned spatialDimension = bulk.mesh_meta_data().spatial_dimension();
  stk::mesh::Bucket& bucket = bulk.bucket(entity);
  stk::topology objTopology = bucket.topology();

  const unsigned numParametricCoordinates = get_number_of_parametric_coordinates(objTopology);

  STK_ThrowRequireMsg(numParametricCoordinates <= parametricCoordinates.size(),
                      "Size of input parametric coordinate does not match topology = " << objTopology << ".");

  const spmd::EntityKeyPair key(entity, bulk.entity_key(entity));
  std::vector<double> coordinates;

  impl::gather_nodal_coordinates(data, entity, coordinates);

  evaluatedLocation.resize(spatialDimension);
  auto meTopo = SearchTopology(objTopology, key, &bucket);
  data.masterElemProvider->evaluate_field(meTopo, parametricCoordinates,
                                          spatialDimension, coordinates, evaluatedLocation);

  std::vector<double> interpolatedParametricCoordinates(numParametricCoordinates);
  double parametricDistance;

  data.masterElemProvider->find_parametric_coordinates(meTopo, spatialDimension,
                                                       coordinates, evaluatedLocation,
                                                       interpolatedParametricCoordinates,
                                                       parametricDistance);

  return parametricDistance;
}

} // end namespace search
} // end namespace stk

