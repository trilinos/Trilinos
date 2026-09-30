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
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/BucketConnectivity.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/base/NgpUtils.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"
#include "stk_mesh/baseImpl/DeviceBucketRepository.hpp"
#include "stk_topology/topology.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/FieldBLAS.hpp"
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
  void modify_values_in_test_field(VectorType const& modifiedEntities, double valueOffset)
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
            entityValues(stk::mesh::ComponentIdx{j}) = valueOffset + entity.local_offset();
          }
        }
      }
    );
  }

  // Without STK_USE_DEVICE_MESH, DeviceSpace aliases HostSpace and no separate device allocation is
  // ever created, so has_device_data() is legitimately false.
  void expect_has_device_data([[maybe_unused]] stk::mesh::FieldBase const& field)
  {
#ifdef STK_USE_DEVICE_MESH
    EXPECT_TRUE(field.has_device_data());
#endif
  }

  template <typename VectorType>
  void check_field_data_on_device(VectorType const& modifiedEntities, double expectedOffset)
  {
    auto fieldData = testField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    expect_has_device_data(testField);

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

  template <typename VectorType>
  void check_constant_field_value_on_device(stk::mesh::Field<double>& field, VectorType const& entities,
                                            double expectedValue)
  {
    auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    expect_has_device_data(field);

    Kokkos::View<stk::mesh::Entity*> entityView("entities", entities.size());
    fill_view(entities, entityView);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entity = entityView(i);
          auto entityValues = fieldData.entity_values(entity);

          NGP_EXPECT_TRUE(entityValues.is_field_defined());
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            NGP_EXPECT_EQ(expectedValue, entityValues(component));
          }
        }
      }
    );
  }

  void check_constant_field_value_on_host(stk::mesh::Field<double>& field,
                                          stk::mesh::EntityVector const& entities,
                                          double expectedValue)
  {
    auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

    for (stk::mesh::Entity entity : entities) {
      auto entityValues = fieldData.entity_values(entity);

      EXPECT_TRUE(entityValues.is_field_defined());
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        EXPECT_DOUBLE_EQ(expectedValue, entityValues(component));
      }
    }
  }

  template <typename VectorType>
  void check_field_component_on_device(stk::mesh::Field<double>& field, VectorType const& entities,
                                       int component, double expectedValue)
  {
    auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    expect_has_device_data(field);

    Kokkos::View<stk::mesh::Entity*> entityView("entities", entities.size());
    fill_view(entities, entityView);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entity = entityView(i);
          auto entityValues = fieldData.entity_values(entity);

          NGP_EXPECT_TRUE(entityValues.is_field_defined());
          NGP_EXPECT_EQ(expectedValue, entityValues(stk::mesh::ComponentIdx{component}));
        }
      }
    );
  }

  Kokkos::View<double*> make_device_view(std::vector<double> const& values)
  {
    Kokkos::View<double*> valuesView("values", values.size());
    auto valuesHost = Kokkos::create_mirror_view(valuesView);

    for (unsigned i = 0; i < values.size(); ++i) {
      valuesHost(i) = values[i];
    }
    Kokkos::deep_copy(valuesView, valuesHost);
    return valuesView;
  }

  // The helpers above share one value across every component.  These walk the slots in ascending
  // (copy, component) order against an explicit list of values, so a transposed or mis-strided copy
  // is visible.
  template <typename FieldType>
  void set_field_values_on_device(FieldType& field, stk::mesh::EntityVector const& entities,
                                  std::vector<double> const& values)
  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

    Kokkos::View<stk::mesh::Entity*> entityView("entities", entities.size());
    fill_view(entities, entityView);
    Kokkos::View<double*> valuesView = make_device_view(values);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entityValues = fieldData.entity_values(entityView(i));

          for (stk::mesh::CopyIdx copy : entityValues.copies()) {
            for (stk::mesh::ComponentIdx component : entityValues.components()) {
              entityValues(copy, component) = valuesView(copy()*entityValues.num_components() + component());
            }
          }
        }
      }
    );
  }

  template <typename FieldType>
  void check_field_values_on_device(FieldType& field, stk::mesh::EntityVector const& entities,
                                    std::vector<double> const& expectedValues)
  {
    auto fieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    expect_has_device_data(field);

    Kokkos::View<stk::mesh::Entity*> entityView("entities", entities.size());
    fill_view(entities, entityView);
    Kokkos::View<double*> expectedView = make_device_view(expectedValues);

    Kokkos::parallel_for(1,
      KOKKOS_LAMBDA(const int) {
        for (unsigned i = 0; i < entityView.extent(0); ++i) {
          auto entityValues = fieldData.entity_values(entityView(i));

          NGP_EXPECT_TRUE(entityValues.is_field_defined());
          NGP_EXPECT_EQ(static_cast<int>(expectedView.extent(0)),
                        entityValues.num_copies()*entityValues.num_components());
          for (stk::mesh::CopyIdx copy : entityValues.copies()) {
            for (stk::mesh::ComponentIdx component : entityValues.components()) {
              NGP_EXPECT_EQ(expectedView(copy()*entityValues.num_components() + component()),
                            entityValues(copy, component));
            }
          }
        }
      }
    );
  }

  template <typename FieldType>
  void set_field_values_on_host(FieldType& field, stk::mesh::EntityVector const& entities,
                                std::vector<double> const& values)
  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::HostSpace>();

    for (stk::mesh::Entity entity : entities) {
      auto entityValues = fieldData.entity_values(entity);

      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = values[copy()*entityValues.num_components() + component()];
        }
      }
    }
  }

  template <typename FieldType>
  void check_field_values_on_host(FieldType& field, stk::mesh::EntityVector const& entities,
                                  std::vector<double> const& expectedValues)
  {
    auto fieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::HostSpace>();

    for (stk::mesh::Entity entity : entities) {
      auto entityValues = fieldData.entity_values(entity);

      EXPECT_TRUE(entityValues.is_field_defined());
      ASSERT_EQ(static_cast<int>(expectedValues.size()),
                entityValues.num_copies()*entityValues.num_components());
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          EXPECT_DOUBLE_EQ(expectedValues[copy()*entityValues.num_components() + component()],
                           entityValues(copy, component));
        }
      }
    }
  }

  stk::mesh::FastMeshIndex host_fast_mesh_index(stk::mesh::Entity entity)
  {
    const stk::mesh::MeshIndex& meshIndex = bulk.mesh_index(entity);
    return stk::mesh::FastMeshIndex{meshIndex.bucket->bucket_id(), meshIndex.bucket_ordinal};
  }

  stk::mesh::Entity declare_entity_on_device(stk::topology::rank_t rank, unsigned entityId,
                                             stk::mesh::PartVector const& addParts)
  {
    Kokkos::View<unsigned*> entityIds("entityIds", 1);
    auto entityIds_host = Kokkos::create_mirror_view(entityIds);
    entityIds_host(0) = entityId;
    Kokkos::deep_copy(entityIds, entityIds_host);

    Kokkos::View<stk::mesh::PartOrdinal*> addPartsView("addPartOrds", addParts.size());
    fill_view(addParts, addPartsView);

    Kokkos::View<stk::mesh::Entity*> requestedEntities("requestedEntities", 1);

    stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    deviceMesh.batch_declare_entities(rank, entityIds, addPartsView, requestedEntities);

    auto requested_host = Kokkos::create_mirror_view(requestedEntities);
    Kokkos::deep_copy(requested_host, requestedEntities);
    return requested_host(0);
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

TEST_F(DeviceMeshModFieldTester, check_device_field_initialized_when_moved_to_part_with_new_field)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  const double partFieldInitValue = 99.0;
  auto& partField = meta.declare_field<double>(stk::topology::ELEM_RANK, "part_only_field");
  stk::mesh::put_field_on_mesh(partField, *block2, &partFieldInitValue);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});

  check_constant_field_value_on_device(partField, stk::mesh::EntityVector{elem1}, partFieldInitValue);
}

TEST_F(DeviceMeshModFieldTester, check_device_field_zeroed_when_moved_to_part_with_new_uninitialized_field)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  auto& partField = meta.declare_field<double>(stk::topology::ELEM_RANK, "part_only_field_no_init");
  stk::mesh::put_field_on_mesh(partField, *block2, nullptr);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});

  check_constant_field_value_on_device(partField, stk::mesh::EntityVector{elem1}, 0.0);
}

TEST_F(DeviceMeshModFieldTester, check_host_field_initialized_after_sync_to_host_when_moved_to_part_with_new_field)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  const double partFieldInitValue = 99.0;
  auto& partField = meta.declare_field<double>(stk::topology::ELEM_RANK, "part_only_field");
  stk::mesh::put_field_on_mesh(partField, *block2, &partFieldInitValue);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock}, stk::mesh::PartVector{removeBlock});

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  deviceMesh.update_bulk_data();

  check_constant_field_value_on_host(partField, stk::mesh::EntityVector{elem1}, partFieldInitValue);
}

// A newly-declared Entity has no prior location, so its field data must come from the Field's
// registered initial value.  Declare a node, not an element: modification_end() requires every
// EDGE_RANK-through-ELEMENT_RANK entity to have connected nodes, but exempts nodes themselves.
TEST_F(DeviceMeshModFieldTester, check_device_field_initialized_when_entity_created)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  const double nodeFieldInitValue = 42.0;
  auto& nodeField = meta.declare_field<double>(stk::topology::NODE_RANK, "new_node_field");
  auto& nodePart = meta.declare_part_with_topology("new_node_part", stk::topology::NODE);
  stk::mesh::put_field_on_mesh(nodeField, nodePart, &nodeFieldInitValue);

  setup_mesh("generated:1x1x2");

  auto newNode = declare_entity_on_device(stk::topology::NODE_RANK, 99u, stk::mesh::PartVector{&nodePart});

  check_constant_field_value_on_device(nodeField, stk::mesh::EntityVector{newNode}, nodeFieldInitValue);
}

// Growing a Field's per-entity size across a part move is rejected, matching the host invariant in
// BulkData::copy_entity_fields_callback().  The check is an assert, so it compiles out under NDEBUG.
#ifndef NDEBUG

TEST_F(DeviceMeshModFieldTester, check_field_grow_on_move_is_rejected)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  auto& growField = meta.declare_field<double>(stk::topology::ELEM_RANK, "grow_field");
  const double growFieldInitVals[3] = {7.0, 7.0, 7.0};
  stk::mesh::put_field_on_mesh(growField, *meta.get_part("block_1"), 1, growFieldInitVals);  // source: 1 comp
  stk::mesh::put_field_on_mesh(growField, *block2, 3, growFieldInitVals);                    // dest:   3 comp

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  EXPECT_ANY_THROW(change_entity_parts(stk::mesh::EntityVector{elem1},
                                       stk::mesh::PartVector{addBlock},
                                       stk::mesh::PartVector{removeBlock}));
}

// With numCopies > 1 a clamped copy misplaces data rather than merely under-filling it: growing
// 2 components x 2 copies to 3 x 2, the leading 4 scalars put source (copy 1, component 0) into
// destination (copy 0, component 2), since components are the inner extent.
TEST_F(DeviceMeshModFieldTester, check_field_grow_on_move_is_rejected_with_multiple_copies)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  auto& growField = meta.declare_field<double>(stk::topology::ELEM_RANK, "grow_multi_copy_field");
  const double growFieldInitVals[6] = {7.0, 7.0, 7.0, 7.0, 7.0, 7.0};
  const unsigned numCopies = 2;
  stk::mesh::put_field_on_mesh(growField, *meta.get_part("block_1"), 2, numCopies, growFieldInitVals);
  stk::mesh::put_field_on_mesh(growField, *block2, 3, numCopies, growFieldInitVals);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  EXPECT_ANY_THROW(change_entity_parts(stk::mesh::EntityVector{elem1},
                                       stk::mesh::PartVector{addBlock},
                                       stk::mesh::PartVector{removeBlock}));
}

#endif  // NDEBUG

// Write a distinct value into every slot and check the exact (copy, component) ordering on device and,
// after a sync, on host (where the layout may be Layout::Right while device is always Layout::Left).
TEST_F(DeviceMeshModFieldTester, check_device_field_component_order_preserved_across_move)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  auto& multiField = meta.declare_field<double>(stk::topology::ELEM_RANK, "multi_comp_copy_field");
  const double multiFieldInitVals[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  const unsigned numComponents = 3;
  const unsigned numCopies = 2;
  stk::mesh::put_field_on_mesh(multiField, *meta.get_part("block_1"), numComponents, numCopies, multiFieldInitVals);
  stk::mesh::put_field_on_mesh(multiField, *block2, numComponents, numCopies, multiFieldInitVals);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  const std::vector<double> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  set_field_values_on_device(multiField, stk::mesh::EntityVector{elem1}, values);

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock},
                      stk::mesh::PartVector{removeBlock});

  check_field_values_on_device(multiField, stk::mesh::EntityVector{elem1}, values);

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  deviceMesh.update_bulk_data();

  check_field_values_on_host(multiField, stk::mesh::EntityVector{elem1}, values);
}

// A Field newly-defined on the destination part is filled from its registered initial value.  Use a
// multi-component Field, so each init value must land in its own component.
TEST_F(DeviceMeshModFieldTester, check_device_field_per_component_init_values_when_moved_to_part_with_new_field)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  auto& partField = meta.declare_field<double>(stk::topology::ELEM_RANK, "part_only_multi_comp_field");
  const double partFieldInitVals[3] = {5.0, 6.0, 7.0};
  const unsigned numComponents = 3;
  stk::mesh::put_field_on_mesh(partField, *block2, numComponents, partFieldInitVals);

  setup_mesh("generated:1x1x2");

  auto addBlock = meta.get_part("block_2");
  auto removeBlock = meta.get_part("block_1");
  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];

  change_entity_parts(stk::mesh::EntityVector{elem1}, stk::mesh::PartVector{addBlock},
                      stk::mesh::PartVector{removeBlock});

  const std::vector<double> expectedValues{5.0, 6.0, 7.0};
  check_field_values_on_device(partField, stk::mesh::EntityVector{elem1}, expectedValues);

  stk::mesh::NgpMesh& deviceMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  deviceMesh.update_bulk_data();

  check_field_values_on_host(partField, stk::mesh::EntityVector{elem1}, expectedValues);
}

// copy_entity_field_bytes<HostSpace>() dispatches on Field::host_data_layout(), which DefaultHostLayout
// would otherwise pin to one branch.  Declare both layouts explicitly to cover both.
TEST_F(DeviceMeshModFieldTester, check_host_field_bytes_copied_for_both_host_layouts)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  auto& leftField = meta.declare_field<double, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "left_field");
  auto& rightField = meta.declare_field<double, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "right_field");
  const double layoutFieldInitVals[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  const unsigned numComponents = 3;
  const unsigned numCopies = 2;
  stk::mesh::put_field_on_mesh(leftField, meta.universal_part(), numComponents, numCopies, layoutFieldInitVals);
  stk::mesh::put_field_on_mesh(rightField, meta.universal_part(), numComponents, numCopies, layoutFieldInitVals);

  ASSERT_EQ(stk::mesh::Layout::Left, leftField.host_data_layout());
  ASSERT_EQ(stk::mesh::Layout::Right, rightField.host_data_layout());

  setup_mesh("generated:1x1x2");

  auto& buckets = bulk.buckets(stk::topology::ELEM_RANK);
  auto elem1 = (*buckets[0])[0];
  auto elem2 = (*buckets[0])[1];

  const std::vector<double> srcValues{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const std::vector<double> destValues{11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
  set_field_values_on_host(leftField, stk::mesh::EntityVector{elem1}, srcValues);
  set_field_values_on_host(rightField, stk::mesh::EntityVector{elem1}, srcValues);
  set_field_values_on_host(leftField, stk::mesh::EntityVector{elem2}, destValues);
  set_field_values_on_host(rightField, stk::mesh::EntityVector{elem2}, destValues);

  std::vector<stk::mesh::FieldBase*> fields{&leftField, &rightField};
  stk::mesh::copy_entity_field_bytes<stk::ngp::HostSpace>(fields, host_fast_mesh_index(elem1),
                                                          host_fast_mesh_index(elem2));

  check_field_values_on_host(leftField, stk::mesh::EntityVector{elem2}, srcValues);
  check_field_values_on_host(rightField, stk::mesh::EntityVector{elem2}, srcValues);
  check_field_values_on_host(leftField, stk::mesh::EntityVector{elem1}, srcValues);
  check_field_values_on_host(rightField, stk::mesh::EntityVector{elem1}, srcValues);
}

}
