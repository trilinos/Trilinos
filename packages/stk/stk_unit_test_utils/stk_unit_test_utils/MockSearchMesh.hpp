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

#include "stk_transfer_util/MockMasterElementHex8.hpp"
#include "stk_transfer_util/MockMasterElementLine2.hpp"
#include "stk_transfer_util/MockMasterElementQuad4.hpp"
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
#include "stk_search_util/PointEvaluator.hpp"   // for EvaluatePointsInte...
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_topology/topology.hpp"              // for topology, topol...
#include "stk_util/parallel/Parallel.hpp"         // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"        // for ThrowRequireMsg

namespace stk {
namespace unit_test_util {
class Hex8SendMesh;
class Hex8RecvMesh;
class SinglePointMesh;
}
}

namespace stk
{
namespace search
{
template <>
struct MeshTraits<stk::unit_test_util::Hex8SendMesh> {
  using Entity = stk::mesh::Entity;
  using EntityVec = std::vector<Entity>;
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Box, EntityProc>;
};
}  // namespace search
}  // namespace stk

namespace stk
{
namespace search
{
template <>
struct MeshTraits<stk::unit_test_util::Hex8RecvMesh> {
  using Entity = stk::mesh::Entity;
  using EntityVec = std::vector<Entity>;
  using EntityKey = std::pair<stk::mesh::EntityKey, int>;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Sphere, EntityProc>;
};
}  // namespace search
}  // namespace stk

namespace stk
{
namespace search
{
//BEGINMeshTrait
template <>
struct MeshTraits<stk::unit_test_util::SinglePointMesh> {
  using Entity = int;
  using EntityVec = std::vector<Entity>;
  using EntityKey = int;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Sphere, EntityProc>;
};
//ENDMeshTrait
}  // namespace search
}  // namespace stk

namespace stk {
namespace unit_test_util {

class Hex8SendMesh : public stk::search::SourceMeshInterface<Hex8SendMesh>
{
 public:
  Hex8SendMesh(stk::mesh::BulkData& bulkData,
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

  stk::ParallelMachine comm() const override { return m_comm; }

  std::string name() const override { return m_name; }

  void set_name(const std::string& meshName) override { m_name = meshName; }

  std::vector<std::string> get_part_membership(const EntityKey& k) const override
  {
    return stk::search::get_part_membership(m_bulk, k, m_parts);
  }

  void bounding_boxes(std::vector<BoundingBox>& boxes, [[maybe_unused]] bool /*includeGhosts*/=false) const override
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

  void find_parametric_coords(const EntityKey& k,
      const std::vector<double>& toCoords,
      std::vector<double>& parametricCoords,
      double& parametricDistance,
      bool& isWithinParametricTolerance) const override
  {
    stk::mesh::Entity elem = m_bulk.get_entity(k);
    stk::topology topology = m_bulk.bucket(elem).topology();
    STK_ThrowRequireMsg(topology == stk::topology::HEX_8, "Invalid topology: " << topology.name());

    // load nodal coordinates from element
    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    const auto numNodes = m_bulk.num_nodes(elem);
    unsigned nDim = m_meta.spatial_dimension();
    std::vector<double> transposedElementCoords(nDim * numNodes);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*m_coordinateField,
      [&](auto& coordFieldData) {
        for (auto ni = 0u; ni < numNodes; ++ni) {
          stk::mesh::Entity node = nodes[ni];

          auto fromCoords = coordFieldData.entity_values(node);
          for (stk::mesh::ComponentIdx j : fromCoords.components()) {
            const auto offSet = ni + j * numNodes;
            transposedElementCoords[offSet] = fromCoords(j);
          }
        }
      }
    );

    parametricCoords.assign(3, std::numeric_limits<double>::max());
    parametricDistance = stk::transfer_util::Hex8::is_in_element(transposedElementCoords.data(), toCoords.data(), parametricCoords.data());

    isWithinParametricTolerance = parametricDistance <= (1 + m_parametricTolerance);
  }

  bool modify_search_outside_parametric_tolerance(const EntityKey& /*k*/,
      const std::vector<double>& /*toCoords*/,
      std::vector<double>& /*parametricCoords*/,
      double& /*geometricDistanceSquared*/,
      bool& /*isWithinGeometricTolerance*/) const override
  {
    return false;
  }

  //BEGINSource_get_distance_from_nearest_node
  double get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& point) const override
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);

    STK_ThrowRequireMsg(
        m_bulk.entity_rank(e) == stk::topology::ELEM_RANK, "Invalid entity rank for object: " << m_bulk.entity_rank(e));

    double minDistance = std::numeric_limits<double>::max();
    const stk::mesh::Entity* const nodes = m_bulk.begin_nodes(e);
    const int num_nodes = m_bulk.num_nodes(e);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*m_coordinateField,
      [&](auto& coordFieldData) {
        for (int i = 0; i < num_nodes; ++i) {
          double d = 0.0;
          auto node_coordinates = coordFieldData.entity_values(nodes[i]);

          for (stk::mesh::ComponentIdx j : node_coordinates.components()) {
            const double t = point[j] - node_coordinates(j);
            d += t * t;
          }
          if (d < minDistance) minDistance = d;
        }
      }
    );

    minDistance = std::sqrt(minDistance);
    return minDistance;
  }
  //ENDSource_get_distance_from_nearest_node

  double get_closest_geometric_distance_squared(const EntityKey& k, const std::vector<double>& toCoords) const override
  {
    double distance = get_distance_from_nearest_node(k, toCoords);
    return distance * distance;
  }

  double get_distance_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const override
  {
    double distanceSquared = get_distance_squared_from_centroid(k, toCoords);
    return std::sqrt(distanceSquared);
  }

  double get_distance_squared_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const override
  {
    std::vector<double> centroidVec;
    centroid(k, centroidVec);

    const unsigned nDim = m_meta.spatial_dimension();
    return stk::search::distance_sq(nDim, centroidVec.data(), toCoords.data());
  }

  void centroid(const EntityKey& k, std::vector<double>& centroidVec) const override
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);
    const unsigned ndim = m_coordinateField->mesh_meta_data().spatial_dimension();
    stk::search::determine_centroid(ndim, e, *m_coordinateField, centroidVec);
  }

  void coordinates(const EntityKey& k, std::vector<double>& coords) const override
  {
    centroid(k, coords);
  }

  void initialize() override { }

  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const override
  {
    return stk::search::ObjectOutsideDomainPolicy::IGNORE;
  }

  void update_ghosting(const EntityProcVec& /*entity_keys*/, const std::string& /*suffix*/ = "") override { }

  void update_ghosted_key(EntityKey& /*k*/) override { }

  void post_mesh_modification_event() override { }

  void destroy_ghosting() override { }

 protected:
  stk::mesh::MetaData& m_meta;
  stk::mesh::BulkData& m_bulk;
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

 private:
  stk::mesh::PartVector m_parts;
  const stk::ParallelMachine m_comm;
  const double m_parametricTolerance;

  std::string m_name{"Hex8SourceMesh"};

  Hex8SendMesh(const Hex8SendMesh&) = delete;
  const Hex8SendMesh& operator()(const Hex8SendMesh&) = delete;

  void fill_bounding_box(stk::mesh::Entity elem, stk::search::Point<double>& min_corner,
                         stk::search::Point<double>& max_corner) const
  {
    const unsigned nDim = m_meta.spatial_dimension();

    STK_ThrowRequireMsg(m_bulk.is_valid(elem), "Invalid entity: " << m_bulk.entity_key(elem));

    for (unsigned j = 0; j < nDim; ++j) {
      min_corner[j] = std::numeric_limits<double>::max();
      max_corner[j] = std::numeric_limits<double>::lowest();
    }

    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    int numNodes = m_bulk.num_nodes(elem);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*m_coordinateField,
      [&](auto& coordsFieldData) {
        for (int ni = 0; ni < numNodes; ++ni) {
          stk::mesh::Entity node = nodes[ni];

          auto coords = coordsFieldData.entity_values(node);

          //for (unsigned j = 0; j < nDim; ++j) {
          for (stk::mesh::ComponentIdx j : coords.components()) {
            min_corner[j] = std::min(min_corner[j], coords(j));
            max_corner[j] = std::max(max_corner[j], coords(j));
          }
        }
      }
    );
  }
};

class Hex8RecvMesh : public stk::search::DestinationMeshInterface<Hex8RecvMesh>
{
 public:

  Hex8RecvMesh(stk::mesh::BulkData& bulkData,
               const stk::mesh::FieldBase* coordinateField,
               const stk::mesh::PartVector& recvParts,
               const stk::ParallelMachine recvComm,
               std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator,
               double parametricTolerance, double geometricTolerance)
  : m_bulk(bulkData)
  , m_meta(bulkData.mesh_meta_data())
  , m_coordinateField(coordinateField)
  , m_parts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_pointEvaluator(pointEvaluator)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
  {

  }

  Hex8RecvMesh(stk::mesh::BulkData& bulkData,
               const stk::mesh::FieldBase* coordinateField,
               const stk::mesh::PartVector& recvParts,
               const stk::mesh::Selector& activeSelector,
               const stk::ParallelMachine recvComm,
               std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator,
               double parametricTolerance, double geometricTolerance)
  : m_bulk(bulkData)
  , m_meta(bulkData.mesh_meta_data())
  , m_coordinateField(coordinateField)
  , m_parts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(activeSelector)
  , m_pointEvaluator(pointEvaluator)
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
  {

  }

  Hex8RecvMesh(stk::mesh::BulkData& bulkData,
               const stk::mesh::PartVector& recvParts,
               const stk::ParallelMachine recvComm,
               const double parametricTolerance,
               double geometricTolerance)
  : m_bulk(bulkData)
  , m_meta(bulkData.mesh_meta_data())
  , m_coordinateField(bulkData.mesh_meta_data().coordinate_field())
  , m_parts(recvParts)
  , m_comm(recvComm)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_pointEvaluator(std::make_shared<stk::search::CentroidEvaluator>(m_bulk, m_coordinateField))
  , m_parametricTolerance(parametricTolerance)
  , m_searchTolerance(geometricTolerance)
  {
    consistency_check();
  }

  stk::ParallelMachine comm() const override { return m_comm; }

  std::string name() const override { return m_name; }

  void set_name(const std::string& meshName) override { m_name = meshName; }

  //BEGINDestination_bounding_boxes
  void bounding_boxes(std::vector<BoundingBox>& v_box) const override
  {
    const unsigned nDim = m_meta.spatial_dimension();

    stk::mesh::Selector selector = stk::search::get_objects_selector(m_meta, m_parts, &m_activeSelector);
    stk::mesh::BucketVector const& buckets = m_bulk.get_buckets(stk::topology::ELEM_RANK, selector);

    std::vector<double> coords;

    for(auto&& ib : buckets) {
      stk::mesh::Bucket& b = *ib;

      for(auto entity : b) {

        stk::mesh::EntityKey eKey = m_bulk.entity_key(entity);
        stk::search::spmd::EntityKeyPair spmdKey(stk::search::spmd::make_entity_key_pair(m_bulk, entity));

        stk::topology topo = b.topology();
        size_t numPoints = m_pointEvaluator->num_points(spmdKey, topo);

        for(size_t p = 0; p < numPoints; ++p) {
          m_pointEvaluator->coordinates(spmdKey, p, coords);

          Point center;
          if(nDim == 2) {
            center = Point(coords[0], coords[1]);
          }
          else {
            center = Point(coords[0], coords[1], coords[2]);
          }

          EntityKey key(eKey, p);
          EntityProc theIdent(key, m_bulk.parallel_rank());
          BoundingBox theBox(Sphere(center, m_searchTolerance), theIdent);

          v_box.push_back(theBox);
        }
      }
    }
    std::sort(v_box.begin(), v_box.end(),
              [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
  }
  //ENDDestination_bounding_boxes

  //BEGINDestination_get_distance_from_nearest_node
  double get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& point) const override
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k.first);

    STK_ThrowRequireMsg(
        m_bulk.entity_rank(e) == stk::topology::ELEM_RANK, "Invalid entity rank for object: " << m_bulk.entity_rank(e));

    double minDistance = std::numeric_limits<double>::max();
    const stk::mesh::Entity* const nodes = m_bulk.begin_nodes(e);
    const int num_nodes = m_bulk.num_nodes(e);

    stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(*m_coordinateField,
      [&](auto& coordsFieldData) {
        for (int i = 0; i < num_nodes; ++i) {
          double d = 0.0;
          auto node_coordinates = coordsFieldData.entity_values(nodes[i]);

          for (stk::mesh::ComponentIdx j : node_coordinates.components()) {
            const double t = point[j] - node_coordinates(j);
            d += t * t;
          }
          if (d < minDistance) minDistance = d;
        }
      }
    );

    minDistance = std::sqrt(minDistance);
    return minDistance;
  }
  //ENDDestination_get_distance_from_nearest_node

  void centroid(const EntityKey& k, std::vector<double>& centroidVec) const override
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k.first);
    const unsigned ndim = m_coordinateField->mesh_meta_data().spatial_dimension();
    stk::search::determine_centroid(ndim, e, *m_coordinateField, centroidVec);
  }

  void coordinates(const EntityKey& k, std::vector<double>& coords) const override
  {
    size_t p = k.second;
    stk::search::spmd::EntityKeyPair spmdKey(stk::search::spmd::make_entity_key_pair(m_bulk, k.first));
    m_pointEvaluator->coordinates(spmdKey, p, coords);
  }

  double get_parametric_tolerance() const final { return m_parametricTolerance; }

  double get_search_tolerance() const final { return m_searchTolerance; }

 protected:
  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

 private:
  stk::mesh::PartVector m_parts;
  const stk::ParallelMachine m_comm;
  const stk::mesh::Selector m_activeSelector;

  std::shared_ptr<stk::search::PointEvaluatorInterface> m_pointEvaluator;

  const double m_parametricTolerance;
  const double m_searchTolerance;

  std::string m_name{"Hex8DestinationMesh"};

  Hex8RecvMesh(const Hex8RecvMesh&) = delete;
  const Hex8RecvMesh& operator()(const Hex8RecvMesh&) = delete;

  void consistency_check() const
  {
    for (const stk::mesh::Part* part : m_parts) {
      STK_ThrowRequireMsg(
          part->primary_entity_rank() == stk::topology::ELEM_RANK, "All source parts must be {ELEM_RANK}");
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

  stk::ParallelMachine comm() const override { return m_comm; };

  std::string name() const override { return m_name; }

  void set_name(const std::string& meshName) override { m_name = meshName; }

  //BEGINDestination_bounding_boxes
  void bounding_boxes(std::vector<BoundingBox>& v) const override
  {
    Point center(m_coords[0], m_coords[1], m_coords[2]);

    EntityKey key = 1;
    EntityProc theIdent(key, stk::parallel_machine_rank(m_comm));
    BoundingBox theBox(Sphere(center, m_geometricTolerance), theIdent);
    v.push_back(theBox);
  }
  //ENDDestination_bounding_boxes

  void coordinates(const EntityKey& /*k*/, std::vector<double>& coords) const override
  {
    coords.assign(m_coords, m_coords + 3);
  }
  double get_search_tolerance() const override { return m_geometricTolerance; }
  double get_parametric_tolerance() const override { return m_parametricTolerance; }

  void centroid(const EntityKey& /*k*/, std::vector<double>& centroidVec) const override
  {
    centroidVec.assign(m_coords, m_coords + 3);
  }
  double get_distance_from_nearest_node(const EntityKey& /*k*/, const std::vector<double>& toCoords) const override
  {
    return stk::search::distance(3, m_coords, toCoords.data());
  }

  void initialize() override { }

 private:
  const stk::ParallelMachine m_comm;
  double m_coords[3];
  double m_parametricTolerance = 0.00001;
  double m_geometricTolerance = 0.1;

  std::string m_name{"SinglePointMesh"};
};

}}  // namespace stk::unit_test_util
#endif
