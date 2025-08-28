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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_NODESENDMESH_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_NODESENDMESH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/EntityKey.hpp"                // for EntityKey
#include "stk_mesh/base/Field.hpp"                    // for Field
#include "stk_mesh/base/Types.hpp"                    // for PartVector, Ent...
#include "stk_search/Sphere.hpp"                      // for Sphere
#include "stk_search/IdentProc.hpp"                   // for IdentProc
#include "stk_search/FilterCoarseSearch.hpp"          // for ObjectOutsideDo...
#include "stk_search/Point.hpp"                       // for Point
#include "stk_search_util/FindParametricCoordinates.hpp"  // for FindParametr...
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_search/SearchInterface.hpp"
#include "stk_util/parallel/Parallel.hpp"             // for ParallelMachine

#include <cstddef>                                    // for size_t
#include <memory>                                     // for shared_ptr
#include <set>                                        // for set
#include <string>                                     // for string
#include <utility>                                    // for pair
#include <vector>                                     // for vector

namespace stk::mesh { class BulkData; }
namespace stk::mesh { class FieldBase; }
namespace stk::mesh { class Ghosting; }
namespace stk::mesh { class MetaData; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace search { namespace spmd { class NodeSendMesh; } } }

namespace stk {
namespace search {

template <>
struct MeshTraits<spmd::NodeSendMesh> {
  using Entity = stk::mesh::Entity;
  using EntityVec = std::vector<Entity>;
  using EntityKey = stk::search::spmd::EntityKeyPair;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Sphere, EntityProc>;
};

namespace spmd {

using NodeSendMeshSearchBaseClass = stk::search::SourceMeshInterface<NodeSendMesh>;

class NodeSendMesh : public NodeSendMeshSearchBaseClass  {
 public:
  using BaseClass = NodeSendMeshSearchBaseClass;

  NodeSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
               const stk::mesh::PartVector& sendParts, const stk::ParallelMachine sendComm,
               const double coarseSearchTolerance);

  NodeSendMesh(stk::mesh::BulkData* sendBulk, const stk::mesh::FieldBase* coordinateField,
               const stk::mesh::PartVector& sendParts, const stk::mesh::Selector& activeSelector,
               const stk::ParallelMachine sendComm, const double coarseSearchTolerance);

  virtual ~NodeSendMesh() = default;

  stk::ParallelMachine comm() const final { return m_comm; }

  void bounding_boxes(std::vector<BoundingBox>& v, bool includeGhosts=false) const override;

  virtual void find_parametric_coords(const EntityKey& k, const std::vector<double>& toCoords,
                                      std::vector<double>& parametricCoords,
                                      double& parametricDistance,
                                      bool& isWithinParametricTolerance) const override;

  virtual bool modify_search_outside_parametric_tolerance(const EntityKey& k, const std::vector<double>& toCoords,
                                                          std::vector<double>& parametricCoords,
                                                          double& geometricDistanceSquared,
                                                          bool& isWithinGeometricTolerance) const override;

  virtual double get_closest_geometric_distance_squared(const EntityKey& k,
                                                        const std::vector<double>& toCoords) const override;

  virtual double get_distance_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const override;

  virtual double get_distance_squared_from_centroid(const EntityKey& k, const std::vector<double>& toCoords) const override;

  virtual double get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const override;

  virtual void centroid(const EntityKey& k, std::vector<double>& centroid) const override;

  void coordinates(const EntityKey& k, std::vector<double>& coords) const override;

  virtual std::string name() const override { return m_name; }

  virtual void set_name(const std::string& meshName) override { m_name = meshName; }

  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const override { return m_extrapolateOption; }

  virtual void update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix = "") override;

  virtual void update_ghosted_key(EntityKey& k) override;

  virtual void initialize() override;

  virtual void post_mesh_modification_event() override;


  virtual void destroy_ghosting() override;

  stk::mesh::EntityId id(const EntityKey& k) const;

  std::vector<std::string> get_part_membership(const EntityKey& k) const override;

  bool is_valid(const EntityKey& e) const;

  void set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy option) { m_extrapolateOption = option; }

  void set_mesh_modified(bool modified) { m_meshModified = modified; }
  bool is_mesh_modified() const { return m_meshModified; }

  stk::mesh::Ghosting* get_ghosting() { return m_ghosting; }

  virtual std::string get_inspector_delimiter() const { return get_time_stamp(); }

  const stk::mesh::BulkData* get_bulk() const { return m_bulk; }
  const stk::mesh::MetaData* get_meta() const { return m_meta; }

 protected:
  stk::mesh::BulkData* m_bulk{nullptr};
  stk::mesh::MetaData* m_meta{nullptr};
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

  stk::mesh::PartVector m_meshParts;
  const stk::ParallelMachine m_comm;
  const stk::mesh::Selector m_activeSelector;

  bool m_meshModified{false};
  stk::mesh::Ghosting* m_ghosting{nullptr};
  stk::mesh::EntityProcVec m_ghostingMap;

  std::shared_ptr<FindParametricCoordsInterface> m_findParametricCoords;

  mutable std::vector<double> m_coordVector;
  size_t m_syncCount{0};

  const double m_coarseSearchTol;

  bool m_isInitialized{false};
  stk::search::ObjectOutsideDomainPolicy m_extrapolateOption{stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG};

  void consistency_check();

 private:
  std::string m_name{"<UNKNOWN SEND NODE SEARCH MESH>"};

  NodeSendMesh(const NodeSendMesh&) = delete;
  const NodeSendMesh& operator()(const NodeSendMesh&) = delete;
};

} // namespace spmd
} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_NODESENDMESH_HPP_ */
