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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_ELEMRECVMESH_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_ELEMRECVMESH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/Entity.hpp"             // for Entity
#include "stk_mesh/base/EntityKey.hpp"          // for EntityKey
#include "stk_mesh/base/Field.hpp"              // for Field
#include "stk_mesh/base/Part.hpp"               // for Part
#include "stk_mesh/base/Selector.hpp"           // for Selector
#include "stk_mesh/base/Types.hpp"              // for EntityRank, PartVector
#include "stk_search_util/CachedEntity.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/PointEvaluator.hpp"   // for EvaluatePointsInte...
#include "stk_search_util/spmd/ElementRecvMesh.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_transfer/TransferInterface.hpp"
#include "stk_transfer_util/PointInterpolation.hpp"
#include "stk_util/parallel/Parallel.hpp"       // for ParallelMachine

#include <memory>                               // for shared_ptr
#include <set>                                  // for set
#include <string>                               // for string
#include <utility>                              // for pair
#include <vector>                               // for vector

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace transfer { namespace spmd { class ElementRecvMesh; } } }

namespace stk {
namespace search {

template <>
struct MeshTraits<stk::transfer::spmd::ElementRecvMesh> {
  using Entity        = stk::search::spmd::ElementRecvMesh::Entity;
  using EntityVec     = stk::search::spmd::ElementRecvMesh::EntityVec;
  using EntityKey     = stk::search::spmd::ElementRecvMesh::EntityKey;
  using EntityProc    = stk::search::spmd::ElementRecvMesh::EntityProc;
  using EntityProcVec = stk::search::spmd::ElementRecvMesh::EntityProcVec;
  using Point         = stk::search::spmd::ElementRecvMesh::Point;
  using Box           = stk::search::spmd::ElementRecvMesh::Box;
  using Sphere        = stk::search::spmd::ElementRecvMesh::Sphere;
  using BoundingBox   = stk::search::spmd::ElementRecvMesh::BoundingBox;
};
}
}

namespace stk {
namespace transfer {

template <>
struct MeshTraits<stk::transfer::spmd::ElementRecvMesh> {
  using Entity        = stk::search::spmd::ElementRecvMesh::Entity;
  using EntityVec     = stk::search::spmd::ElementRecvMesh::EntityVec;
  using EntityKey     = stk::search::spmd::ElementRecvMesh::EntityKey;
  using EntityProc    = stk::search::spmd::ElementRecvMesh::EntityProc;
  using EntityProcVec = stk::search::spmd::ElementRecvMesh::EntityProcVec;
  using Point         = stk::search::spmd::ElementRecvMesh::Point;
  using Box           = stk::search::spmd::ElementRecvMesh::Box;
  using Sphere        = stk::search::spmd::ElementRecvMesh::Sphere;
  using BoundingBox   = stk::search::spmd::ElementRecvMesh::BoundingBox;
};

namespace spmd {

using ElementRecvMeshSearchBaseClass   = stk::search::DestinationMeshInterface<ElementRecvMesh>;
using ElementRecvMeshTransferBaseClass = stk::transfer::DestinationMeshInterface<ElementRecvMesh>;

class ElementRecvMesh : public ElementRecvMeshSearchBaseClass,
                        public ElementRecvMeshTransferBaseClass {
 public:
  using SearchBaseClass   = ElementRecvMeshSearchBaseClass;
  using TransferBaseClass = ElementRecvMeshTransferBaseClass;

  // Explicit enumeration due to multiple inheritance
  using Entity        = TransferBaseClass::Entity;
  using EntityVec     = TransferBaseClass::EntityVec;
  using EntityKey     = TransferBaseClass::EntityKey;
  using EntityProc    = TransferBaseClass::EntityProc;
  using EntityProcVec = TransferBaseClass::EntityProcVec;
  using Point         = TransferBaseClass::Point;
  using Box           = TransferBaseClass::Box;
  using Sphere        = TransferBaseClass::Sphere;
  using BoundingBox   = TransferBaseClass::BoundingBox;

  ElementRecvMesh(stk::mesh::BulkData* recvBulk,
                  const stk::mesh::FieldBase* coordinateField,
                  const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                  const stk::mesh::EntityRank transferEntityRank,
                  const stk::mesh::PartVector& recvParts,
                  const stk::ParallelMachine recvComm,
                  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator,
                  double parametricTolerance, double geometricTolerance);

  ElementRecvMesh(stk::mesh::BulkData* recvBulk,
                  const stk::mesh::FieldBase* coordinateField,
                  const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                  const stk::mesh::EntityRank transferEntityRank,
                  const stk::mesh::PartVector& recvParts,
                  const stk::mesh::Selector& activeSelector,
                  const stk::ParallelMachine recvComm,
                  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator,
                  double parametricTolerance, double geometricTolerance);

  virtual ~ElementRecvMesh() = default;

  stk::ParallelMachine comm() const final { return m_comm; }

  double get_parametric_tolerance() const final { return m_parametricTolerance; }

  double get_search_tolerance() const final { return m_searchTolerance; }

  virtual void bounding_boxes(std::vector<BoundingBox>& v) const override;

  virtual void coordinates(const EntityKey& k, std::vector<double>& coords) const override;

  virtual std::string name() const override { return m_name; }

  virtual void set_name(const std::string& meshName) override { m_name = meshName; }

  virtual double get_distance_from_nearest_node(const EntityKey& k, const std::vector<double>& toCoords) const override;

  virtual void centroid(const EntityKey& k, std::vector<double>& centroid) const override;

  virtual void initialize() override;


  virtual double* value(const EntityKey& k, const unsigned fieldIndex) const override;

  virtual unsigned value_size(const EntityKey& k, const unsigned fieldIndex) const override;

  virtual unsigned num_values(const EntityKey& e) const override;

  virtual unsigned max_num_values() const override;

  virtual unsigned value_key(const EntityKey& k, const unsigned fieldIndex) const override;

  virtual void update_values() override;

  virtual unsigned get_index(const unsigned i) const override;


  stk::mesh::EntityId id(const EntityKey& k) const;

  TransferBaseClass::Entity entity(const EntityKey& k) const;

  void fill_entity_keys(const stk::mesh::EntityKeyVector& rangeEntities, std::vector<EntityKey>& elementEntityKeys);

  void initialize(stk::mesh::BulkData* recvBulk);

  virtual std::string get_inspector_delimiter() const { return stk::search::get_time_stamp(); }

  unsigned spatial_dimension() const;

  const stk::mesh::BulkData* get_bulk() const { return m_bulk; }
  const stk::mesh::MetaData* get_meta() const { return m_meta; }

 protected:
  stk::mesh::BulkData* m_bulk{nullptr};
  stk::mesh::MetaData* m_meta{nullptr};
  const stk::mesh::FieldBase* m_coordinateField{nullptr};
  const stk::mesh::EntityRank m_transferEntityRank;

  std::vector<stk::transfer::FieldSpec> m_fieldSpecs;
  std::vector<IndexedField> m_fieldVec;

  stk::mesh::PartVector m_meshParts;
  const stk::ParallelMachine m_comm;
  const stk::mesh::Selector m_activeSelector;
  bool m_meshModified{false};

  std::shared_ptr<stk::search::PointEvaluatorInterface> m_pointEvaluator;

  const double m_parametricTolerance;
  const double m_searchTolerance;

  stk::search::spmd::ElementRecvMesh m_searchMesh;

  bool m_isInitialized{false};

  void consistency_check();

 private:
  std::string m_name{"<UNKNOWN RECV ELEMENT TRANSFER MESH>"};

  ElementRecvMesh(const ElementRecvMesh&) = delete;
  const ElementRecvMesh& operator()(const ElementRecvMesh&) = delete;
};

} // namespace spmd
} // namespace transfer
} // namespace stk


#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_ELEMRECVMESH_HPP_ */
