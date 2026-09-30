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

#include <gtest/gtest.h>
#include <algorithm>
#include <iterator>
#include <map>
#include "Kokkos_Core.hpp"
#include "ngp/NgpUnitTestUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SkinMesh.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"

#ifndef UnitTestNgpMeshModificationUtils_hpp
#define UnitTestNgpMeshModificationUtils_hpp

#ifdef STK_USE_DEVICE_MESH
namespace
{
using ngp_unit_test_utils::check_bucket_layout;
using ngp_unit_test_utils::check_entity_parts_on_device;

using DeviceEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::MemSpace>;
using DevicePartOrdinalsType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::MemSpace>;

using HostEntitiesType = Kokkos::View<stk::mesh::Entity*, stk::ngp::HostExecSpace>;
using HostPartOrdinalsType = Kokkos::View<stk::mesh::PartOrdinal*, stk::ngp::HostExecSpace>;

template <typename MeshType>
void confirm_host_mesh_is_not_synchronized_from_device(const MeshType& ngpMesh)
{
  if constexpr (std::is_same_v<MeshType, stk::mesh::DeviceMesh>) {
    EXPECT_TRUE(ngpMesh.needs_update_bulk_data());
  }
  else {
    EXPECT_FALSE(ngpMesh.needs_update_bulk_data());  // If host build, HostMesh can't ever be stale
  }
}

template <typename MeshType>
void confirm_host_mesh_is_synchronized_from_device(const MeshType& ngpMesh)
{
  EXPECT_FALSE(ngpMesh.needs_update_bulk_data());
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_connectivity(DeviceMeshType& ngpMesh, stk::mesh::EntityRank entityRank, const EntitiesViewType& entities,
                            stk::mesh::EntityRank connRank, unsigned expectedNumConnected)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      auto connectedEntities = ngpMesh.get_connected_entities(entityRank, fastMeshIndex, connRank);
      NGP_EXPECT_EQ(expectedNumConnected, connectedEntities.size());
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType>
void check_part_is_on_device(DeviceMeshType& ngpMesh, stk::mesh::Part& part, stk::topology::rank_t rank) {
  const auto partOrdinal = part.mesh_meta_data_ordinal();
  bool is_present = false;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& local_result) {
    local_result = ngpMesh.get_device_bucket_repository().get_part_rank(partOrdinal) == rank;
  }, is_present);
  EXPECT_TRUE(is_present);
}


template <typename HostMeshType, typename EntitiesViewType>
void check_host_connectivity(HostMeshType& ngpMesh, stk::mesh::EntityRank entityRank, const EntitiesViewType& entities,
                            stk::mesh::EntityRank connRank, unsigned expectedNumConnected)
{
  for (unsigned idx = 0; idx < entities.extent(0); ++idx) {
    auto entity = entities(idx);
    auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
    auto connectedEntities = ngpMesh.get_connected_entities(entityRank, fastMeshIndex, connRank);
    EXPECT_EQ(expectedNumConnected, connectedEntities.size());
  }
}

template <typename ViewDataType, typename... ViewArgs>
void fill_views(const Kokkos::View<ViewDataType*, ViewArgs...>& view,
                const std::vector<ViewDataType>& data)
{
  STK_ThrowRequire(view.extent(0) == data.size());
  using ViewType = Kokkos::View<ViewDataType*, ViewArgs...>;
  using ConstHostViewType = typename ViewType::host_mirror_type::const_type;
  auto host_entities = ConstHostViewType(data.data(), view.extent(0));
  Kokkos::deep_copy(view, host_entities);
}

class NgpMeshMod : public ::ngp_testing::Test
{
public:
  NgpMeshMod()
  {
  }

  void build_empty_mesh(unsigned initialBucketCapacity, unsigned maximumBucketCapacity)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
  }

  void commit_meta_data()
  {
    m_meta->commit();
  }

protected:
  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData * m_meta;
};

class NgpBatchDeclareDestroyEntities : public NgpMeshMod
{
public:
  NgpBatchDeclareDestroyEntities() {}
};

inline stk::mesh::Entity create_node(stk::mesh::BulkData& bulk, stk::mesh::EntityId nodeId,
                              const stk::mesh::PartVector& initialParts = stk::mesh::PartVector())
{
  bulk.modification_begin();
  stk::mesh::Entity newNode = bulk.declare_node(nodeId, initialParts);
  bulk.modification_end();

  return newNode;
}

class NgpBatchDestroyConnectivities : public NgpMeshMod
{
public:
  NgpBatchDestroyConnectivities() {};

  void create_one_hex_mesh() {
    stk::io::fill_mesh("generated:1x1x1", *m_bulk);
    stk::mesh::Part& nodePart = m_meta->declare_part_with_topology("nodePart", stk::topology::NODE);
    stk::mesh::Part& quadPart = m_meta->declare_part_with_topology("quadPart", stk::topology::QUAD_4);
    stk::mesh::Part& edgePart = m_meta->declare_part_with_topology("edgePart", stk::topology::LINE_2);
    stk::mesh::Part& hexPart = m_meta->declare_part_with_topology("hexPart", stk::topology::HEX_8);

    constexpr auto connectFacesToEdges = true;
    stk::mesh::create_edges(*m_bulk, m_meta->universal_part(), &edgePart);
    stk::mesh::create_all_sides(*m_bulk, m_meta->universal_part(), {&quadPart}, connectFacesToEdges);

    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::ELEMENT_RANK, hexes);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::FACE_RANK, quads);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::EDGE_RANK, edges);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::NODE_RANK, nodes);

    m_bulk->modification_begin();
    m_bulk->change_entity_parts(nodes, stk::mesh::PartVector{&nodePart});
    m_bulk->change_entity_parts(edges, stk::mesh::PartVector{&edgePart});
    m_bulk->change_entity_parts(quads, stk::mesh::PartVector{&quadPart});
    m_bulk->change_entity_parts(hexes, stk::mesh::PartVector{&hexPart});
    m_bulk->modification_end();
  }

  void create_two_hex_mesh() {
    stk::io::fill_mesh("generated:1x1x2", *m_bulk);
    stk::mesh::Part& nodePart = m_meta->declare_part_with_topology("nodePart", stk::topology::NODE);
    stk::mesh::Part& quadPart = m_meta->declare_part_with_topology("quadPart", stk::topology::QUAD_4);
    stk::mesh::Part& edgePart = m_meta->declare_part_with_topology("edgePart", stk::topology::LINE_2);
    stk::mesh::Part& hexPart = m_meta->declare_part_with_topology("hexPart", stk::topology::HEX_8);

    constexpr auto connectFacesToEdges = true;
    stk::mesh::create_edges(*m_bulk, m_meta->universal_part(), &edgePart);
    stk::mesh::create_all_sides(*m_bulk, m_meta->universal_part(), {&quadPart}, connectFacesToEdges);

    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::ELEMENT_RANK, hexes);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::FACE_RANK, quads);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::EDGE_RANK, edges);
    stk::mesh::get_entities(*m_bulk, stk::mesh::EntityRank::NODE_RANK, nodes);

    m_bulk->modification_begin();
    m_bulk->change_entity_parts(nodes, stk::mesh::PartVector{&nodePart});
    m_bulk->change_entity_parts(edges, stk::mesh::PartVector{&edgePart});
    m_bulk->change_entity_parts(quads, stk::mesh::PartVector{&quadPart});
    m_bulk->change_entity_parts(hexes, stk::mesh::PartVector{&hexPart});
    m_bulk->modification_end();
  }

  std::vector<stk::mesh::Entity> nodes;
  std::vector<stk::mesh::Entity> edges;
  std::vector<stk::mesh::Entity> quads;
  std::vector<stk::mesh::Entity> hexes;

};

inline DevicePartOrdinalsType create_device_part_ordinal(stk::mesh::PartVector const& vector)
{
  HostPartOrdinalsType hostPartOrdinals("", vector.size());
  for (unsigned i = 0; i < vector.size(); ++i) {
    hostPartOrdinals(i) = vector[i]->mesh_meta_data_ordinal();
  }
  DevicePartOrdinalsType devicePartOrdinals = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace{}, hostPartOrdinals);
  Kokkos::fence();
  return devicePartOrdinals;
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_part_ordinal_match(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank, EntitiesViewType const& entities, DevicePartOrdinalsType const& expectedPartOrdinals)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);

      auto bucket = deviceBucketRepo.get_bucket(rank, fastMeshIndex.bucket_id);
      auto& bucketPartOrdinals = bucket->get_part_ordinals();

      for (int i = expectedPartOrdinals.extent(0)-1, j = bucketPartOrdinals.extent(0)-1; i >= 0; --i, --j) {
        NGP_EXPECT_EQ(bucketPartOrdinals(j), expectedPartOrdinals(i));
      }

      auto partition = deviceBucketRepo.get_partition(rank, bucket->partition_id());
      auto partitionPartOrdinals = partition->superset_part_ordinals();

      for (int i = expectedPartOrdinals.extent(0)-1, j = partitionPartOrdinals.extent(0)-1; i >= 0; --i, --j) {
        NGP_EXPECT_EQ(partitionPartOrdinals(j), expectedPartOrdinals(i));
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_has_parts(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank, EntitiesViewType const& entities, DevicePartOrdinalsType const& expectedPartOrdinals)
{
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);

      auto& bucket = ngpMesh.get_bucket(rank, fastMeshIndex.bucket_id);

      for (unsigned i = 0; i<expectedPartOrdinals.extent(0); ++i) {
        NGP_EXPECT_TRUE(bucket.member(expectedPartOrdinals(i)));
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType>
void check_device_mesh_indices(DeviceMeshType& ngpMesh, stk::mesh::EntityRank rank)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(deviceBucketRepo.num_buckets(rank),
    KOKKOS_LAMBDA(const int idx) {
      auto& buckets = deviceBucketRepo.m_buckets[rank];
      auto& bucket = buckets[idx];

      if (!bucket.is_active()) { return; }

      for (unsigned i = 0; i < bucket.size(); ++i) {
        auto entity = bucket[i];

        if (!entity.is_local_offset_valid()) { continue; }

        auto fastMeshIndex = ngpMesh.fast_mesh_index(entity);
        NGP_EXPECT_EQ(fastMeshIndex.bucket_id, bucket.bucket_id());
        NGP_EXPECT_EQ(fastMeshIndex.bucket_ord, i);
      }
    }
  );
  Kokkos::fence();
}

template<typename EntitiesHostViewType>
void init_host_field_data(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& field, EntitiesHostViewType entities)
{
  auto fieldData = field.data<stk::mesh::ReadWrite,stk::ngp::HostSpace>();
  for(unsigned idx=0; idx<entities.extent(0); ++idx) {
    auto entity = entities(idx);
    stk::mesh::EntityId id = mesh.identifier(entity);
    auto fieldEntityValues = fieldData.entity_values(entity);
    for(stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
      fieldEntityValues(i) = id;
    }
  }
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_field_data_is_id(DeviceMeshType& ngpMesh,
    stk::mesh::EntityRank /*rank*/,
    stk::mesh::Field<double>& field,
    EntitiesViewType const& entities)
{
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      stk::mesh::EntityId id = ngpMesh.identifier(entity);
      auto fieldEntityValues = fieldData.entity_values(fastMeshIndex);

      for (stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
        stk::mesh::EntityId fieldValue = static_cast<stk::mesh::EntityId>(fieldEntityValues(i));
        NGP_EXPECT_EQ(id, fieldValue);
      }
    }
  );
  Kokkos::fence();
}

template <typename DeviceMeshType, typename EntitiesViewType>
void check_device_entity_field_data(DeviceMeshType& ngpMesh,
    stk::mesh::EntityRank /*rank*/,
    stk::mesh::Field<double>& field,
    EntitiesViewType const& entities,
    double expectedValue)
{
  constexpr double tol = 1.e-12;
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(entities.extent(0),
    KOKKOS_LAMBDA(const int idx) {
      auto entity = entities(idx);
      auto fastMeshIndex = ngpMesh.device_mesh_index(entity);
      auto fieldEntityValues = fieldData.entity_values(fastMeshIndex);

      for (stk::mesh::ComponentIdx i : fieldEntityValues.components()) {
        double fieldValue = fieldEntityValues(i);
        NGP_EXPECT_NEAR(expectedValue, fieldValue, tol);
      }
    }
  );
  Kokkos::fence();
}

}  // namespace
#endif

#endif
