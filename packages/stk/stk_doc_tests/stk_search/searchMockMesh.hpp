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

#ifndef STK_DOC_TEST_SEARCH_MOCK_MESH_HPP
#define STK_DOC_TEST_SEARCH_MOCK_MESH_HPP

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

#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_mesh/base/Entity.hpp"           // for Entity
#include "stk_mesh/base/Part.hpp"             // for Part
#include "stk_mesh/base/Selector.hpp"         // for Selector, opera...
#include "stk_mesh/base/Types.hpp"            // for EntityRank, Ent...
#include "stk_search/Box.hpp"                 // for Box
#include "stk_search/DistanceComparison.hpp"  // for stk_distance
#include "stk_search/FilterCoarseSearch.hpp"
#include "stk_search/IdentProc.hpp"  // for IdentProc
#include "stk_search/Point.hpp"      // for Point
#include "stk_search/SearchInterface.hpp"
#include "stk_search/Sphere.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"  // for compute_entity_centroid
#include "stk_topology/topology.hpp"              // for topology, topol...
#include "stk_util/parallel/Parallel.hpp"         // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"        // for ThrowRequireMsg
#include "searchHex8.hpp"

namespace doc_test
{
class Hex8SourceMesh;
class SinglePointMesh;
}  // namespace doc_test

namespace stk
{
namespace search
{
template <>
struct MeshTraits<doc_test::Hex8SourceMesh> {
  using Entity = stk::mesh::Entity;
  using EntityVec = std::vector<Entity>;
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using BoundingBox = std::pair<Box, EntityProc>;
};
}  // namespace search
}  // namespace stk

namespace stk
{
namespace search
{
//BEGINMeshTrait
template <>
struct MeshTraits<doc_test::SinglePointMesh> {
  using Entity = int;
  using EntityVec = std::vector<Entity>;
  using EntityKey = int;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Sphere, EntityProc>;
};
//ENDMeshTrait
}  // namespace search
}  // namespace stk

namespace doc_test
{

class Hex8SourceMesh : public stk::search::SourceMeshInterface<Hex8SourceMesh>
{
 public:
  Hex8SourceMesh(stk::mesh::BulkData& bulkData,
      const stk::mesh::PartVector& sendParts,
      const stk::ParallelMachine comm,
      const double parametricTolerance)
      : m_meta(bulkData.mesh_meta_data()),
        m_bulk(bulkData),
        m_coordinateField(bulkData.mesh_meta_data().coordinate_field()),
        m_parts(sendParts),
        m_comm(comm),
        m_parametricTolerance(parametricTolerance) 
  {
    for (const stk::mesh::Part* part : sendParts) {
      STK_ThrowRequireMsg(
          part->primary_entity_rank() == stk::topology::ELEM_RANK, "All source parts must be {ELEM_RANK}");
    }
  }

  stk::ParallelMachine comm() const { return m_comm; }

  std::string name() const { return "Hex8SourceMesh"; }

  //BEGINSource_bounding_boxes
  void bounding_boxes(std::vector<BoundingBox>& boxes) const
  {
    Point min_corner, max_corner;

    stk::mesh::Selector selector = stk::mesh::selectUnion(m_parts);
    stk::mesh::BucketVector const& buckets = m_bulk.get_buckets(stk::topology::ELEM_RANK, selector);

    for (auto&& ib : buckets) {
      stk::mesh::Bucket& b = *ib;

      for (auto elem : b) {
        fill_bounding_box(elem, min_corner, max_corner);

        EntityProc theIdent(m_bulk.entity_key(elem), m_bulk.parallel_rank());
        BoundingBox theBox(Box(min_corner, max_corner), theIdent);
        boxes.push_back(theBox);
      }
    }
    std::sort(boxes.begin(), boxes.end(),
        [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
  }
  //ENDSource_bounding_boxes

  //BEGINSource_find_parametric_coords
  void find_parametric_coords(const EntityKey k,
      const double* toCoords,
      std::vector<double>& parametricCoords,
      double& parametricDistance,
      bool& isWithinParametricTolerance) const
  {
    stk::mesh::Entity elem = m_bulk.get_entity(k);
    stk::topology topology = m_bulk.bucket(elem).topology();
    STK_ThrowRequireMsg(topology == stk::topology::HEX_8, "Invalid topology: " << topology.name());

    // load nodal coordinates from element
    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    const auto numNodes = m_bulk.num_nodes(elem);
    unsigned nDim = m_meta.spatial_dimension();

    std::vector<double> transposedElementCoords(nDim * numNodes);

    for (auto ni = 0u; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      const double* fromCoords = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, node));
      for (unsigned j = 0; j < nDim; ++j) {
        const auto offSet = ni + j * numNodes;
        transposedElementCoords[offSet] = fromCoords[j];
      }
    }

    parametricCoords.assign(3, std::numeric_limits<double>::max());
    parametricDistance = Hex8::is_in_element(transposedElementCoords.data(), toCoords, parametricCoords.data());

    isWithinParametricTolerance = parametricDistance <= (1 + m_parametricTolerance);
  }
  //ENDSource_find_parametric_coords

  bool modify_search_outside_parametric_tolerance(const EntityKey k,
      const double* toCoords,
      std::vector<double>& parametricCoords,
      double& geometricDistanceSquared,
      bool& isWithinGeometricTolerance) const
  {
    return false;
  }

  //BEGINSource_get_distance_from_nearest_node
  double get_distance_from_nearest_node(const EntityKey k, const double* point) const
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);

    STK_ThrowRequireMsg(
        m_bulk.entity_rank(e) == stk::topology::ELEM_RANK, "Invalid entity rank for object: " << m_bulk.entity_rank(e));

    double minDistance = std::numeric_limits<double>::max();
    const unsigned nDim = m_meta.spatial_dimension();

    const stk::mesh::Entity* const nodes = m_bulk.begin_nodes(e);
    const int num_nodes = m_bulk.num_nodes(e);

    for (int i = 0; i < num_nodes; ++i) {
      double d = 0.0;
      double* node_coordinates = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, nodes[i]));

      for (unsigned j = 0; j < nDim; ++j) {
        const double t = point[j] - node_coordinates[j];
        d += t * t;
      }
      if (d < minDistance) minDistance = d;
    }

    minDistance = std::sqrt(minDistance);
    return minDistance;
  }
  //ENDSource_get_distance_from_nearest_node

  double get_closest_geometric_distance_squared(const EntityKey k, const double* toCoords) const
  {
    double distance = get_distance_from_nearest_node(k, toCoords);
    return distance * distance;
  }

  double get_distance_from_centroid(const EntityKey k, const double* toCoords) const
  {
    double distanceSquared = get_distance_squared_from_centroid(k, toCoords);
    return std::sqrt(distanceSquared);
  }

  double get_distance_squared_from_centroid(const EntityKey k, const double* toCoords) const
  {
    std::vector<double> centroidVec;
    centroid(k, centroidVec);

    const unsigned nDim = m_meta.spatial_dimension();
    return stk::search::distance_sq(nDim, centroidVec.data(), toCoords);
  }

  void centroid(const EntityKey k, std::vector<double>& centroidVec) const
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);
    stk::search::compute_entity_centroid(e, *m_coordinateField, centroidVec);
  }

  const double* coord(const EntityKey k) const
  {
    centroid(k, m_coordVector);
    return m_coordVector.data();
  }

 protected:
  stk::mesh::MetaData& m_meta;
  stk::mesh::BulkData& m_bulk;
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

 private:
  stk::mesh::PartVector m_parts;
  const stk::ParallelMachine m_comm;
  const double m_parametricTolerance;

  mutable std::vector<double> m_coordVector;

  Hex8SourceMesh(const Hex8SourceMesh&) = delete;
  const Hex8SourceMesh& operator()(const Hex8SourceMesh&) = delete;

  void fill_bounding_box(
      stk::mesh::Entity elem, stk::search::Point<double>& min_corner, stk::search::Point<double>& max_corner) const
  {
    const unsigned nDim = m_meta.spatial_dimension();

    STK_ThrowRequireMsg(m_bulk.is_valid(elem), "Invalid entity: " << m_bulk.entity_key(elem));

    for (unsigned j = 0; j < nDim; ++j) {
      min_corner[j] = std::numeric_limits<double>::max();
      max_corner[j] = -std::numeric_limits<double>::max();
    }

    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    int numNodes = m_bulk.num_nodes(elem);
    for (int ni = 0; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      double* coords = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, node));

      for (unsigned j = 0; j < nDim; ++j) {
        min_corner[j] = std::min(min_corner[j], coords[j]);
        max_corner[j] = std::max(max_corner[j], coords[j]);
      }
    }
  }
};


class SinglePointMesh : public stk::search::DestinationMeshInterface<SinglePointMesh>
{
 public:
  SinglePointMesh(const stk::ParallelMachine comm, double x, double y, double z, double paramTol, double geomTol)
      : m_comm(comm), m_parametricTolerance(paramTol), m_geometricTolerance(geomTol)
  {
    m_coords[0] = x;
    m_coords[1] = y;
    m_coords[2] = z;
  }

  stk::ParallelMachine comm() const { return m_comm; };

  std::string name() const { return "SinglePointMesh"; }

  //BEGINDestination_bounding_boxes
  void bounding_boxes(std::vector<BoundingBox>& v) const
  {
    Point center(m_coords[0], m_coords[1], m_coords[2]);

    EntityKey key = 1;
    EntityProc theIdent(key, stk::parallel_machine_rank(m_comm));
    BoundingBox theBox(Sphere(center, m_geometricTolerance), theIdent);
    v.push_back(theBox);
  }
  //ENDDestination_bounding_boxes

  const double* coord(const EntityKey k) const { return m_coords; }
  double get_search_tolerance() const { return m_geometricTolerance; }
  double get_parametric_tolerance() const { return m_parametricTolerance; }

  void centroid(const EntityKey k, std::vector<double>& centroidVec) const
  {
    centroidVec.assign(m_coords, m_coords + 3);
  }
  double get_distance_from_nearest_node(const EntityKey k, const double* toCoords) const
  {
    return stk::search::distance(3, m_coords, toCoords);
  }

 private:
  const stk::ParallelMachine m_comm;
  double m_coords[3];
  double m_parametricTolerance = 0.00001;
  double m_geometricTolerance = 0.1;
};

}  // namespace doc_test
#endif
