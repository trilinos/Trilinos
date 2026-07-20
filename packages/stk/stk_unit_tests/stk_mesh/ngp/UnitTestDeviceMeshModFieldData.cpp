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

#ifdef STK_USE_DEVICE_MESH

#include <gtest/gtest.h>
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_topology/topology.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/NgpFieldBLAS.hpp"
#include "stk_unit_test_utils/DeviceBucketTestUtils.hpp"
#include <stk_ngp_test/ngp_test.hpp>

namespace {

class DeviceMeshModFieldTester : public ::ngp_testing::Test
{
public:
  DeviceMeshModFieldTester() :
    meta(3),
    bulk(meta,
         stk::parallel_machine_world(),
         stk::mesh::BulkData::AUTO_AURA,
         false,
         std::unique_ptr<stk::mesh::FieldDataManager>(),
         maximumBucketCapacity,
         maximumBucketCapacity),
    block2(&meta.declare_part("block_2", stk::topology::ELEM_RANK)),
    testField(meta.declare_field<double>(stk::topology::ELEM_RANK, "test_elem_field1"))
  {
    stk::mesh::put_field_on_mesh(testField, meta.universal_part(), &initOffset);
  }

  void setup_mesh(std::string const& meshDesc)
  {
    stk::io::fill_mesh(meshDesc, bulk);
    block1 = bulk.mesh_meta_data().get_part("block_1");

    stk::mesh::EntityVector elems;
    bulk.get_entities(stk::topology::ELEM_RANK, meta.universal_part(), elems);
    modify_values_in_test_field(elems, initOffset);
  }

  void change_entity_parts(stk::mesh::EntityVector const& entities, stk::mesh::PartVector const& addParts, stk::mesh::PartVector const& removeParts)
  {
    Kokkos::View<stk::mesh::Entity*> entityView("entities", entities.size());
    Kokkos::View<stk::mesh::PartOrdinal*> addPartsView("addPartOrds", addParts.size());
    Kokkos::View<stk::mesh::PartOrdinal*> removePartsView("removePartOrds", removeParts.size());

    fill_view(addParts, addPartsView);
    fill_view(removeParts, removePartsView);
    fill_view(entities, entityView);

    stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    deviceMesh.batch_change_entity_parts(entityView, addPartsView, removePartsView);
  }

  template <typename VectorType>
  void modify_values_in_test_field(VectorType const& modifiedEntities, double initOffset)
  {
    auto fieldData = testField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

    Kokkos::View<stk::mesh::Entity*> entityView("entities", modifiedEntities.size());
    fill_view(modifiedEntities, entityView);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entity = entityView(i);
          auto entityValues = fieldData.entity_values(entity);

          for (int j = 0; j < entityValues.num_components(); ++j) {
            entityValues(stk::mesh::ComponentIdx{j}) = initOffset + entity.local_offset();
          }
        }
      }
    );
  }

  template <typename VectorType>
  void check_field_data_on_device(VectorType const& modifiedEntities, double expectedOffset)
  {
    auto fieldData = testField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    EXPECT_TRUE(testField.has_device_data());

    Kokkos::View<stk::mesh::Entity*> entityView("entities", modifiedEntities.size());
    fill_view(modifiedEntities, entityView);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entity = entityView(i);
          auto entityValues = fieldData.entity_values(entity);

          for (int j = 0; j < entityValues.num_components(); ++j) {
            NGP_EXPECT_EQ(expectedOffset + entity.local_offset(), entityValues(stk::mesh::ComponentIdx{j}));
          }
        }
      }
    );
  }

  template <typename VectorType, typename ViewType>
  void fill_view(VectorType const& vector, ViewType view)
  {
    auto hostView = Kokkos::create_mirror_view_and_copy(stk::ngp::HostMemSpace{}, view);

    for (unsigned i = 0; i < hostView.extent(0); ++i) {
      if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vector)>::value_type, stk::mesh::Part*>) {
        hostView(i) = vector[i]->mesh_meta_data_ordinal();
      } else {
        hostView(i) = vector[i];
      }
    }
    Kokkos::deep_copy(view, hostView);
  }

  stk::mesh::MetaData meta;
  stk::unit_test_util::BulkDataTester bulk;
  stk::mesh::Part* block1;
  stk::mesh::Part* block2;
  stk::mesh::Field<double>& testField;

  static constexpr unsigned maximumBucketCapacity = 2;
  static constexpr double initOffset = 4;
};

TEST_F(DeviceMeshModFieldTester, check_device_field_no_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem = (*buckets[0])[0];
  change_entity_parts(stk::mesh::EntityVector{elem}, stk::mesh::PartVector{}, stk::mesh::PartVector{});
  check_field_data_on_device(stk::mesh::EntityVector{elem}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_change_field_value_change_before_skipped_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  double newOffset = 3;
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem = (*buckets[0])[0];
  modify_values_in_test_field(stk::mesh::EntityVector{elem}, newOffset);
  change_entity_parts(stk::mesh::EntityVector{elem}, stk::mesh::PartVector{}, stk::mesh::PartVector{});
  check_field_data_on_device(stk::mesh::EntityVector{elem}, newOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_change_field_value_change_after_skipped_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem = (*buckets[0])[0];
  change_entity_parts(stk::mesh::EntityVector{elem}, stk::mesh::PartVector{}, stk::mesh::PartVector{});

  double newOffset = 3;
  modify_values_in_test_field(stk::mesh::EntityVector{elem}, newOffset);
  check_field_data_on_device(stk::mesh::EntityVector{elem}, newOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_after_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem = (*buckets[0])[0];
  change_entity_parts(stk::mesh::EntityVector{elem}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_after_parts_change_move_one_element)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem = (*buckets[0])[0];
  change_entity_parts(stk::mesh::EntityVector{elem}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_after_parts_change_move_all_elements)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];
  auto elem2 = (*buckets[0])[1];
  change_entity_parts(stk::mesh::EntityVector{elem1, elem2}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem1, elem2}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_after_parts_change_move_all_elements_to_new_partition)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];
  auto elem2 = (*buckets[0])[1];
  change_entity_parts(stk::mesh::EntityVector{elem1, elem2}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem1, elem2}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_after_parts_change_move_partial_elements)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x4");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];
  auto elem2 = (*buckets[0])[1];
  change_entity_parts(stk::mesh::EntityVector{elem1, elem2}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem1, elem2}, initOffset);

  auto elem3 = (*buckets[1])[0];
  auto elem4 = (*buckets[1])[1];
  check_field_data_on_device(stk::mesh::EntityVector{elem3, elem4}, initOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_update_field_data_before_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  auto newOffset = initOffset + 10;
  modify_values_in_test_field(stk::mesh::EntityVector{elem1}, newOffset);

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem1}, newOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_update_field_data_after_parts_change)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});

  auto newOffset = initOffset + 10;
  modify_values_in_test_field(stk::mesh::EntityVector{elem1}, newOffset);
  check_field_data_on_device(stk::mesh::EntityVector{elem1}, newOffset);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_parts_change_to_new_part)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  setup_mesh("generated:1x1x2");

  meta.declare_part("block_3", stk::topology::ELEM_RANK);
  auto addBlock = meta.get_part("block_3");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});
  check_field_data_on_device(stk::mesh::EntityVector{elem1}, initOffset);
}

}

#endif
