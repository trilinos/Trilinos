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

#ifndef NgpDeviceCommMapUpdateTestUtils_hpp
#define NgpDeviceCommMapUpdateTestUtils_hpp

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <stk_ngp_test/ngp_test.hpp>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "../UnitTestUtils.hpp"
#include "ngp/UnitTestNgpMeshModificationUtils.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/CreateFaces.hpp"
#include "stk_mesh/base/DestroyRelations.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/util/ReportHandler.hpp"

#ifdef STK_USE_DEVICE_MESH

namespace {

using DeviceEntityView = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
using DevicePartOrdinalView = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;

template <typename ViewType, typename ValueType>
ViewType make_device_view(const std::string& name, const std::vector<ValueType>& values)
{
  ViewType view(Kokkos::view_alloc(Kokkos::WithoutInitializing, name), values.size());
  fill_views(view, values);
  return view;
}

inline stk::mesh::EntityKey::entity_key_t entity_key_value(stk::mesh::EntityRank rank, stk::mesh::EntityId id)
{
  return stk::mesh::EntityKey(rank, id).m_value;
}

inline std::vector<stk::mesh::EntityKey::entity_key_t> entity_key_values(stk::mesh::EntityRank rank, const std::vector<stk::mesh::EntityId>& ids)
{
  std::vector<stk::mesh::EntityKey::entity_key_t> keyValues;
  keyValues.reserve(ids.size());
  for (stk::mesh::EntityId id : ids) {
    keyValues.push_back(entity_key_value(rank, id));
  }
  return keyValues;
}

inline size_t device_comm_map_size(const stk::mesh::NgpMesh& ngpMesh, stk::mesh::EntityRank rank, int proc,
                                   bool includeGhosts)
{
  Kokkos::View<size_t, stk::ngp::MemSpace> size("commMapSize");
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const int /*i*/) {
                         size() = ngpMesh.volatile_fast_shared_comm_map(rank, proc, includeGhosts).extent(0);
                       });
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), size)();
}

// The FastMeshIndex in each slot of the device comm map for (rank, proc), shared and ghosted, in slot order.
inline std::vector<stk::mesh::FastMeshIndex> device_comm_map_indices(const stk::mesh::NgpMesh& ngpMesh,
                                                                     stk::mesh::EntityRank rank, int proc)
{
  const size_t size = device_comm_map_size(ngpMesh, rank, proc, true);
  Kokkos::View<stk::mesh::FastMeshIndex*, stk::ngp::MemSpace> indices("commMapIndices", size);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, static_cast<int>(size)),
                       KOKKOS_LAMBDA(const int i) {
                         indices(i) = ngpMesh.volatile_fast_shared_comm_map(rank, proc, true)(i);
                       });
  auto hostIndices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), indices);
  return std::vector<stk::mesh::FastMeshIndex>(hostIndices.data(), hostIndices.data() + size);
}

// The key of the entity each slot of the device comm map for (rank, proc) points at, in slot order.  A slot
// whose FastMeshIndex is stale points at some other entity, so it shows up here as the wrong key.
inline std::vector<stk::mesh::EntityKey::entity_key_t> device_comm_map_keys(const stk::mesh::NgpMesh& ngpMesh,
                                                                            stk::mesh::EntityRank rank, int proc,
                                                                            bool includeGhosts)
{
  const size_t size = device_comm_map_size(ngpMesh, rank, proc, includeGhosts);
  Kokkos::View<stk::mesh::EntityKey::entity_key_t*, stk::ngp::MemSpace> keyValues("commMapKeyValues", size);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, static_cast<int>(size)),
                       KOKKOS_LAMBDA(const int i) {
                         const stk::mesh::FastMeshIndex fmi =
                             ngpMesh.volatile_fast_shared_comm_map(rank, proc, includeGhosts)(i);
                         keyValues(i) = ngpMesh.entity_key(ngpMesh.get_entity(rank, fmi)).m_value;
                       });
  auto hostKeyValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keyValues);
  return std::vector<stk::mesh::EntityKey::entity_key_t>(hostKeyValues.data(), hostKeyValues.data() + size);
}

inline void write_device_field_values(const stk::mesh::NgpMesh& ngpMesh, const stk::mesh::Field<double>& field,
                               const std::vector<stk::mesh::Entity>& entities,
                               const std::vector<double>& values)
{
  auto deviceEntities = make_device_view<DeviceEntityView>("fieldWriteEntities", entities);
  auto deviceValues =
      make_device_view<Kokkos::View<double*, stk::ngp::MemSpace>>("fieldWriteValues", values);

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, static_cast<int>(entities.size())),
                       KOKKOS_LAMBDA(const int i) {
                         auto entityValues =
                             fieldData.entity_values(ngpMesh.device_mesh_index(deviceEntities(i)));
                         entityValues(0_comp) = deviceValues(i);
                       });
  Kokkos::fence();
}

inline std::vector<double> read_device_field_values(const stk::mesh::NgpMesh& ngpMesh,
                                             const stk::mesh::Field<double>& field,
                                             const std::vector<stk::mesh::Entity>& entities)
{
  const size_t numEntities = entities.size();
  auto deviceEntities = make_device_view<DeviceEntityView>("fieldReadEntities", entities);
  Kokkos::View<double*, stk::ngp::MemSpace> values("fieldReadValues", numEntities);

  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, static_cast<int>(numEntities)),
                       KOKKOS_LAMBDA(const int i) {
                         auto entityValues =
                             fieldData.entity_values(ngpMesh.device_mesh_index(deviceEntities(i)));
                         values(i) = entityValues(0_comp);
                       });

  auto hostValues = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
  return std::vector<double>(hostValues.data(), hostValues.data() + numEntities);
}

const std::string twoHexMesh = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                               "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";

const std::string threeHexMesh = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                                 "2,3,HEX_8,9,10,11,12,13,14,15,16,block_1";

// A text mesh of one column of 'numLayers' stacked hexes per proc (1 to 4 procs), proc p owning column p.  The
// columns sit in the four quadrants around a vertical axis, so every column has the axis as one vertical edge.
// Viewed from above, one plane of nodes:
//
//          j=0     j=1     j=2
//    i=0    o-------o-------o
//           |  P0   |  P1   |
//    i=1    o-------A-------o
//           |  P3   |  P2   |
//    i=2    o-------o-------o
//
// A is the axis.  Its nodes are 1..numLayers+1, bottom to top; the other nodes are numbered after them, eight
// per plane.  Column p's elements are p*numLayers+1 .. (p+1)*numLayers, bottom to top.
//
// So every proc shares the axis nodes, and once edges are created, the axis edges.  Procs in neighboring
// quadrants (P0-P1, P1-P2, P2-P3, P3-P0) also share the nodes on the side between them, and once faces are
// created, one vertical face per layer.  Procs in opposite quadrants (P0-P2, P1-P3) share only the axis.
inline std::string axis_hex_mesh(int numProcs, int numLayers)
{
  STK_ThrowRequireMsg(numProcs >= 1 && numProcs <= 4,
                      "axis_hex_mesh has one column per quadrant, so it needs 1 to 4 procs, not " << numProcs);
  // The id of the node at plane z (0 at the bottom) and grid point (i, j) in the picture above.
  const int numAxisNodes = numLayers + 1;
  auto nodeId = [&](int z, int i, int j) {
    if (i == 1 && j == 1) { return 1 + z; }
    const int k = i * 3 + j;
    return numAxisNodes + 1 + z * 8 + (k < 4 ? k : k - 1);
  };
  // Column p occupies quadrant (qi, qj) of the grid.
  const int quadrants[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};

  std::ostringstream os;
  int elemId = 1;
  for (int p = 0; p < numProcs; ++p) {
    const int qi = quadrants[p][0];
    const int qj = quadrants[p][1];
    for (int z = 0; z < numLayers; ++z) {
      if (elemId > 1) { os << "\n"; }
      os << p << "," << elemId++ << ",HEX_8";
      for (int dz = 0; dz < 2; ++dz) {
        os << "," << nodeId(z + dz, qi, qj) << "," << nodeId(z + dz, qi + 1, qj) << ","
           << nodeId(z + dz, qi + 1, qj + 1) << "," << nodeId(z + dz, qi, qj + 1);
      }
      os << ",block_1";
    }
  }
  return os.str();
}

class NgpDeviceCommMapUpdate : public ::ngp_testing::Test
{
public:
  void create_bulk_data(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    bulkPtr = stk::mesh::MeshBuilder(stk::parallel_machine_world())
                  .set_spatial_dimension(3)
                  .set_aura_option(auraOption)
                  .set_symmetric_ghost_info(useSymmetricGhostInfo)
                  .create();
    meta = &bulkPtr->mesh_meta_data();

    newPart = &meta->declare_part("newNodePart");
    splitPart = &meta->declare_part("splitNodePart");
    nodeField = &meta->declare_field<double>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*nodeField, meta->universal_part(), nullptr);
    edgeField = &meta->declare_field<double>(stk::topology::EDGE_RANK, "edgeField");
    stk::mesh::put_field_on_mesh(*edgeField, meta->universal_part(), nullptr);
    faceField = &meta->declare_field<double>(stk::topology::FACE_RANK, "faceField");
    stk::mesh::put_field_on_mesh(*faceField, meta->universal_part(), nullptr);
  }

  void setup_hex_mesh(const std::string& meshDesc,
                      stk::mesh::BulkData::AutomaticAuraOption auraOption =
                          stk::mesh::BulkData::NO_AUTO_AURA)
  {
    create_bulk_data(auraOption);
    stk::unit_test_util::setup_text_mesh(*bulkPtr, meshDesc);
  }

  void setup_hex_mesh_with_faces(const std::string& meshDesc,
                                 stk::mesh::BulkData::AutomaticAuraOption auraOption =
                                     stk::mesh::BulkData::NO_AUTO_AURA)
  {
    setup_hex_mesh(meshDesc, auraOption);
    stk::mesh::create_faces(*bulkPtr);
  }

  // An axis_hex_mesh with one column per proc in this run, with edges and faces created.
  void setup_axis_hex_mesh(int numLayers)
  {
    create_bulk_data(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::unit_test_util::setup_text_mesh(*bulkPtr, axis_hex_mesh(bulkPtr->parallel_size(), numLayers));
    stk::mesh::create_edges(*bulkPtr);
    stk::mesh::create_faces(*bulkPtr);
  }

  // Records the keys of the host-built comm map, which check_device_comm_map_unchanged compares the device
  // map against.  Call it after the mesh's sharing and ghosting are set up and before any device modification.
  void capture_nodes_and_gold_values()
  {
    capture_host_comm_map_keys();

    stk::mesh::get_entities(*bulkPtr, stk::topology::NODE_RANK, allNodes);

    nodeValues.clear();
    copyCounts.clear();
    std::vector<int> commProcs;
    for (stk::mesh::Entity node : allNodes) {
      bulkPtr->comm_procs(node, commProcs);
      nodeValues.push_back(static_cast<double>(bulkPtr->identifier(node)));
      copyCounts.push_back(commProcs.size() + 1);
    }
  }

  void capture_gold_values(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& entities,
                           std::vector<double>& values, std::vector<size_t>& counts)
  {
    stk::mesh::get_entities(*bulkPtr, rank, entities);

    values.clear();
    counts.clear();
    std::vector<int> commProcs;
    for (stk::mesh::Entity entity : entities) {
      bulkPtr->comm_procs(entity, commProcs);
      values.push_back(static_cast<double>(bulkPtr->identifier(entity)));
      counts.push_back(commProcs.size() + 1);
    }
  }

  void capture_faces_and_gold_values() { capture_gold_values(stk::topology::FACE_RANK, allFaces, faceValues, faceCopyCounts); }

  std::vector<stk::mesh::Entity> get_owned_shared_nodes() const
  {
    stk::mesh::Selector ownedAndShared = meta->locally_owned_part() & meta->globally_shared_part();
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(*bulkPtr, stk::topology::NODE_RANK, ownedAndShared, nodes, true);
    return nodes;
  }

  void change_parts_on_device(const std::vector<stk::mesh::Entity>& entities,
                              const std::vector<stk::mesh::PartOrdinal>& addPartOrdinals)
  {
    auto entityView = make_device_view<DeviceEntityView>("entities", entities);
    auto addPartView = make_device_view<DevicePartOrdinalView>("addPartOrdinals", addPartOrdinals);
    DevicePartOrdinalView removePartView(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                            "removePartOrdinals"), 0);

    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    ngpMesh.batch_change_entity_parts(entityView, addPartView, removePartView);
  }

  void destroy_on_device(const std::vector<stk::mesh::Entity>& entities)
  {
    auto entityView = make_device_view<DeviceEntityView>("entitiesToDestroy", entities);
    Kokkos::View<bool*, stk::ngp::MemSpace> wasDestroyed(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "wasDestroyed"), entities.size());

    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    ngpMesh.batch_destroy_entities(entityView, wasDestroyed);

    auto hostWasDestroyed = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), wasDestroyed);
    for (size_t i = 0; i < entities.size(); ++i) {
      EXPECT_TRUE(hostWasDestroyed(i)) << "entity " << i << " was not destroyable";
    }
  }

  void declare_nodes_on_device(const std::vector<unsigned>& newNodeIds)
  {
    batch_declare_entities_on_device_for_test(*bulkPtr, stk::topology::NODE_RANK, newNodeIds, {});
  }

  std::vector<stk::mesh::Entity> get_owned_elements() const
  {
    std::vector<stk::mesh::Entity> elems;
    stk::mesh::get_entities(*bulkPtr, stk::topology::ELEM_RANK, meta->locally_owned_part(), elems, true);
    return elems;
  }

  std::vector<stk::mesh::EntityId> get_shared_ids(stk::mesh::EntityRank rank) const
  {
    std::vector<stk::mesh::Entity> entities;
    stk::mesh::get_entities(*bulkPtr, rank, meta->globally_shared_part(), entities, true);
    std::vector<stk::mesh::EntityId> ids;
    for (stk::mesh::Entity entity : entities) { ids.push_back(bulkPtr->identifier(entity)); }
    return ids;
  }

  std::vector<stk::mesh::EntityId> get_shared_face_ids() const { return get_shared_ids(stk::topology::FACE_RANK); }

  // The relations of 'fromEntities' to 'connectedRank', in the flattened form batch_declare_relations takes,
  // so that they can be put back exactly as they were -- ordinals and permutations included.
  struct RelationsToDeclare
  {
    std::vector<stk::mesh::Entity> fromEntities;
    std::vector<unsigned> offsets;
    std::vector<stk::mesh::Entity> toEntities;
    std::vector<stk::mesh::RelationIdentifier> ordinals;
    std::vector<stk::mesh::Permutation> permutations;
  };

  RelationsToDeclare capture_relations(const std::vector<stk::mesh::Entity>& fromEntities,
                                       stk::mesh::EntityRank connectedRank) const
  {
    RelationsToDeclare relations;
    relations.fromEntities = fromEntities;
    relations.offsets.push_back(0);
    for (stk::mesh::Entity from : fromEntities) {
      const unsigned numConn = bulkPtr->num_connectivity(from, connectedRank);
      const stk::mesh::Entity* conn = bulkPtr->begin(from, connectedRank);
      const stk::mesh::ConnectivityOrdinal* ords = bulkPtr->begin_ordinals(from, connectedRank);
      const stk::mesh::Permutation* perms = bulkPtr->begin_permutations(from, connectedRank);
      for (unsigned j = 0; j < numConn; ++j) {
        relations.toEntities.push_back(conn[j]);
        relations.ordinals.push_back(ords[j]);
        relations.permutations.push_back(perms[j]);
      }
      relations.offsets.push_back(static_cast<unsigned>(relations.toEntities.size()));
    }
    if (fromEntities.empty()) { relations.offsets.clear(); }
    return relations;
  }

  void declare_relations_on_device(const std::vector<stk::mesh::Entity>& fromEntities,
                                   const std::vector<std::vector<stk::mesh::Entity>>& downwardEntities)
  {
    std::vector<unsigned> offsets{0};
    std::vector<stk::mesh::Entity> toEntities;
    std::vector<stk::mesh::RelationIdentifier> ordinals;
    for (const std::vector<stk::mesh::Entity>& conn : downwardEntities) {
      for (size_t j = 0; j < conn.size(); ++j) {
        toEntities.push_back(conn[j]);
        ordinals.push_back(static_cast<stk::mesh::RelationIdentifier>(j));
      }
      offsets.push_back(static_cast<unsigned>(toEntities.size()));
    }
    if (fromEntities.empty()) { offsets.clear(); }

    auto fromView = make_device_view<DeviceEntityView>("relationFromEntities", fromEntities);
    auto offsetsView = make_device_view<Kokkos::View<unsigned*, stk::ngp::MemSpace>>("relationOffsets", offsets);
    auto toView = make_device_view<DeviceEntityView>("relationToEntities", toEntities);
    auto ordinalsView =
        make_device_view<Kokkos::View<stk::mesh::RelationIdentifier*, stk::ngp::MemSpace>>("relationOrdinals",
                                                                                         ordinals);

    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    ngpMesh.batch_declare_relations(fromView, offsetsView, toView, ordinalsView);
  }

  void declare_relations_on_device(const RelationsToDeclare& relations)
  {
    auto fromView = make_device_view<DeviceEntityView>("relationFromEntities", relations.fromEntities);
    auto offsetsView =
        make_device_view<Kokkos::View<unsigned*, stk::ngp::MemSpace>>("relationOffsets", relations.offsets);
    auto toView = make_device_view<DeviceEntityView>("relationToEntities", relations.toEntities);
    auto ordinalsView = make_device_view<Kokkos::View<stk::mesh::RelationIdentifier*, stk::ngp::MemSpace>>(
        "relationOrdinals", relations.ordinals);
    auto permutationsView = make_device_view<Kokkos::View<stk::mesh::Permutation*, stk::ngp::MemSpace>>(
        "relationPermutations", relations.permutations);

    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    ngpMesh.batch_declare_relations(fromView, offsetsView, toView, ordinalsView, permutationsView);
  }

  void destroy_relations_on_device(const std::vector<stk::mesh::Entity>& entities, stk::mesh::EntityRank connectedRank)
  {
    auto entityView = make_device_view<DeviceEntityView>("entitiesToDetach", entities);
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    ngpMesh.batch_destroy_relations(entityView, connectedRank);
  }

  void check_comm_map_keys(stk::mesh::EntityRank rank, int proc,
                           const std::vector<stk::mesh::EntityId>& goldSharedIds,
                           const std::vector<stk::mesh::EntityId>& goldGhostIds)
  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);

    auto goldKeyValues = entity_key_values(rank, goldSharedIds);
    EXPECT_EQ(goldKeyValues, device_comm_map_keys(ngpMesh, rank, proc, false))
        << "P" << bulkPtr->parallel_rank() << " rank " << rank << " comm map for proc " << proc
        << ": wrong shared entries, or the right entries in the wrong order";

    const auto goldGhostKeyValues = entity_key_values(rank, goldGhostIds);
    goldKeyValues.insert(goldKeyValues.end(), goldGhostKeyValues.begin(), goldGhostKeyValues.end());
    EXPECT_EQ(goldKeyValues, device_comm_map_keys(ngpMesh, rank, proc, true))
        << "P" << bulkPtr->parallel_rank() << " rank " << rank << " comm map for proc " << proc
        << ": wrong entries, or the right entries in the wrong order";
  }

  void check_node_comm_map_keys(int proc, const std::vector<stk::mesh::EntityId>& goldSharedIds,
                                const std::vector<stk::mesh::EntityId>& goldGhostIds)
  {
    check_comm_map_keys(stk::topology::NODE_RANK, proc, goldSharedIds, goldGhostIds);
  }

  // Keys of the host-built comm map for each (rank, proc, includeGhosts), in slot order.
  using CommMapKeys = std::map<std::tuple<stk::mesh::EntityRank, int, bool>,
                               std::vector<stk::mesh::EntityKey::entity_key_t>>;

  void capture_host_comm_map_keys()
  {
    hostCommMapKeys.clear();
    const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meta->entity_rank_count());
    for (stk::mesh::EntityRank rank = stk::topology::BEGIN_RANK; rank < endRank; ++rank) {
      const stk::mesh::BucketVector& buckets = bulkPtr->buckets(rank);
      for (int proc = 0; proc < bulkPtr->parallel_size(); ++proc) {
        if (proc == bulkPtr->parallel_rank()) { continue; }
        for (bool includeGhosts : {false, true}) {
          auto hostCommMap =
              bulkPtr->template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(rank, proc, includeGhosts);
          std::vector<stk::mesh::EntityKey::entity_key_t>& keyValues = hostCommMapKeys[{rank, proc, includeGhosts}];
          for (size_t i = 0; i < hostCommMap.extent(0); ++i) {
            const stk::mesh::Entity entity = (*buckets[hostCommMap(i).bucket_id])[hostCommMap(i).bucket_ord];
            keyValues.push_back(bulkPtr->entity_key(entity).m_value);
          }
        }
      }
    }
    hostCommMapKeysCaptured = true;
  }

  // No device modification in these tests changes who communicates what, so the device comm map must still
  // hold the host map's entities, slot for slot.  The keys are read through each slot's current
  // FastMeshIndex, so this also catches an index the refresh missed.
  void check_device_comm_map_unchanged()
  {
    ASSERT_TRUE(hostCommMapKeysCaptured) << "capture_nodes_and_gold_values() was not called";
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    for (const auto& [rankProcGhosts, goldKeyValues] : hostCommMapKeys) {
      const auto& [rank, proc, includeGhosts] = rankProcGhosts;
      EXPECT_EQ(goldKeyValues, device_comm_map_keys(ngpMesh, rank, proc, includeGhosts))
          << "P" << bulkPtr->parallel_rank() << " rank " << rank << " proc " << proc
          << (includeGhosts ? " (shared and ghosted)" : " (shared)") << ": device comm map differs from host's";
    }
  }

  void check_parallel_sum_on_device(stk::mesh::Field<double>& field, const std::vector<stk::mesh::Entity>& entities,
                                    const std::vector<double>& values, const std::vector<size_t>& counts)
  {
    auto& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
    write_device_field_values(ngpMesh, field, entities, values);

    stk::mesh::parallel_sum_including_ghosts<stk::ngp::DeviceSpace>(*bulkPtr, {&field});

    const std::vector<double> summedValues = read_device_field_values(ngpMesh, field, entities);
    for (size_t i = 0; i < entities.size(); ++i) {
      EXPECT_DOUBLE_EQ(counts[i] * values[i], summedValues[i])
          << "P" << bulkPtr->parallel_rank() << " rank " << field.entity_rank() << " entity with id " << values[i]
          << ": wrong sum, so the comm map is inconsistent across ranks";
    }
  }

  void check_parallel_sum_on_device() {
    check_parallel_sum_on_device(*nodeField, allNodes, nodeValues, copyCounts);
  }

  void check_face_parallel_sum_on_device()
  {
    check_parallel_sum_on_device(*faceField, allFaces, faceValues, faceCopyCounts);
  }

  bool useSymmetricGhostInfo = true;

  CommMapKeys hostCommMapKeys;
  bool hostCommMapKeysCaptured = false;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr;
  stk::mesh::MetaData* meta = nullptr;
  stk::mesh::Part* newPart = nullptr;
  stk::mesh::Part* splitPart = nullptr;
  stk::mesh::Field<double>* nodeField = nullptr;
  stk::mesh::Field<double>* edgeField = nullptr;
  stk::mesh::Field<double>* faceField = nullptr;
  std::vector<stk::mesh::Entity> allNodes;
  std::vector<double> nodeValues;
  std::vector<size_t> copyCounts;
  std::vector<stk::mesh::Entity> allFaces;
  std::vector<double> faceValues;
  std::vector<size_t> faceCopyCounts;
  std::vector<stk::mesh::Entity> allEdges;
  std::vector<double> edgeValues;
  std::vector<size_t> edgeCopyCounts;
};

}  // namespace

#endif  // STK_USE_DEVICE_MESH

#endif  // NgpDeviceCommMapUpdateTestUtils_hpp
