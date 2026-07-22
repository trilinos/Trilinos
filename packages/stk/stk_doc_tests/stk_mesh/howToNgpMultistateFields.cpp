#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_io/FillMesh.hpp>

namespace {

void set_field_on_host(const stk::mesh::BulkData& stkMesh,
                         stk::mesh::Field<double>& stkField,
                         double fieldVal)
{
  stk::mesh::Selector select(stkField);

  auto stkFieldData = stkField.data<stk::mesh::ReadWrite>();

  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), select,
    [&](const stk::mesh::BulkData& /*bulk*/, stk::mesh::Entity entity) {
      auto entityValues = stkFieldData.entity_values(entity);
      entityValues() = fieldVal;
    });
}

void set_field_on_device(const stk::mesh::NgpMesh& ngpMesh,
                         stk::mesh::Field<double>& stkField,
                         double fieldVal)
{
  stk::mesh::Selector select(stkField);

  auto stkFieldData = stkField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stkField.entity_rank(), select,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = stkFieldData.entity_values(entity);
      entityValues() = fieldVal;
    });
}

template<typename FieldDataType>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         stk::mesh::Selector & select,
                         stk::mesh::EntityRank rank,
                         FieldDataType & stkFieldData,
                         double expectedFieldValue)
{
  stk::mesh::for_each_entity_run(bulk, rank, select,
    [&](const stk::mesh::BulkData& /*mesh*/, stk::mesh::Entity entity) {
      auto entityValues = stkFieldData.entity_values(entity);
      EXPECT_NEAR(entityValues(), expectedFieldValue, 1.e-9);
    });
}

template<typename FieldDataType>
void check_field_on_device(const stk::mesh::NgpMesh& ngpMesh,
                           stk::mesh::Selector & select,
                           stk::mesh::EntityRank rank,
                           FieldDataType & stkFieldData,
                           double expectedFieldValue)
{
  stk::mesh::for_each_entity_run(ngpMesh, rank, select,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = stkFieldData.entity_values(entity);
      NGP_EXPECT_NEAR(entityValues(), expectedFieldValue, 1.e-9);
    });
}

NGP_TEST(NgpMultistateField, setOnHost_swap_checkOnDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                     .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  constexpr unsigned numStates = 2;
  stk::mesh::Field<double>& stkFieldNew = meta.declare_field<double>(stk::topology::ELEM_RANK,
                                                                     "myElemField", numStates);
  stk::mesh::put_field_on_mesh(stkFieldNew, meta.universal_part(), 1, nullptr);

  stk::io::fill_mesh("generated:1x1x1", *bulkPtr);

  EXPECT_EQ(stk::mesh::StateNew, stkFieldNew.state());
  stk::mesh::Field<double>& stkFieldOld = stkFieldNew.field_of_state(stk::mesh::StateOld);
  EXPECT_EQ(stk::mesh::StateOld, stkFieldOld.state());
  stk::mesh::Selector selectFieldOld(stkFieldOld);
  stk::mesh::Selector selectFieldNew(stkFieldNew);
  stk::mesh::EntityRank entityRankFieldOld = stkFieldOld.entity_rank();
  stk::mesh::EntityRank entityRankFieldNew = stkFieldNew.entity_rank();

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;
  set_field_on_host(*bulkPtr, stkFieldOld, oldValue);
  set_field_on_host(*bulkPtr, stkFieldNew, newValue);

  bulkPtr->update_field_data_states();

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);

  auto stkFieldOldData = stkFieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  auto stkFieldNewData = stkFieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData, newValue);
  check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData, oldValue);
}

NGP_TEST(NgpMultistateField, setOnHost_swap_preExistingNgpFieldsNeedSync)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                     .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  constexpr unsigned numStates = 2;
  stk::mesh::Field<double>& stkFieldNew = meta.declare_field<double>(stk::topology::ELEM_RANK,
                                                                     "myElemField", numStates);
  stk::mesh::put_field_on_mesh(stkFieldNew, meta.universal_part(), 1, nullptr);

  stk::io::fill_mesh("generated:1x1x1", *bulkPtr);

  //BEGINNgpMultiStateField
  EXPECT_EQ(stk::mesh::StateNew, stkFieldNew.state());
  stk::mesh::Field<double>& stkFieldOld = stkFieldNew.field_of_state(stk::mesh::StateOld);
  EXPECT_EQ(stk::mesh::StateOld, stkFieldOld.state());
  stk::mesh::Selector selectFieldOld(stkFieldOld);
  stk::mesh::Selector selectFieldNew(stkFieldNew);
  stk::mesh::EntityRank entityRankFieldOld = stkFieldOld.entity_rank();
  stk::mesh::EntityRank entityRankFieldNew = stkFieldNew.entity_rank();

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;
  set_field_on_host(*bulkPtr, stkFieldOld, oldValue);
  set_field_on_host(*bulkPtr, stkFieldNew, newValue);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);

  {
    auto stkFieldOldData = stkFieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto stkFieldNewData = stkFieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData, oldValue);
    check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData, newValue);
  }

  {
    bulkPtr->update_field_data_states();
    auto stkFieldOldData = stkFieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto stkFieldNewData = stkFieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

#ifdef STK_USE_DEVICE_MESH
    check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData, oldValue);
    check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData, newValue);
#else
    check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData, newValue);
    check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData, oldValue);
#endif
  }

  {
    stkFieldOld.synchronize<stk::mesh::ReadWrite>();
    stkFieldNew.synchronize<stk::mesh::ReadWrite>();
    auto stkFieldOldData = stkFieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto stkFieldNewData = stkFieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData, newValue);
    check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData, oldValue);
  }
  //ENDNgpMultiStateField
}

NGP_TEST(NgpMultistateField, setOnDevice_swap_checkOnDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                     .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  constexpr unsigned numStates = 2;
  stk::mesh::Field<double>& stkFieldNew = meta.declare_field<double>(stk::topology::ELEM_RANK,
                                                                     "myElemField", numStates);
  stk::mesh::put_field_on_mesh(stkFieldNew, meta.universal_part(), 1, nullptr);
  stk::io::fill_mesh("generated:1x1x1", *bulkPtr);

  EXPECT_EQ(stk::mesh::StateNew, stkFieldNew.state());
  stk::mesh::Field<double>& stkFieldOld = stkFieldNew.field_of_state(stk::mesh::StateOld);
  EXPECT_EQ(stk::mesh::StateOld, stkFieldOld.state());
  stk::mesh::Selector selectFieldOld(stkFieldOld);
  stk::mesh::Selector selectFieldNew(stkFieldNew);
  stk::mesh::EntityRank entityRankFieldOld = stkFieldOld.entity_rank();
  stk::mesh::EntityRank entityRankFieldNew = stkFieldNew.entity_rank();

  //BEGINNgpMultiStateFieldRotateOnDevice
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;

  set_field_on_device(ngpMesh, stkFieldOld, oldValue);
  set_field_on_device(ngpMesh, stkFieldNew, newValue);

  auto stkFieldOldData = stkFieldOld.data();
  auto stkFieldNewData = stkFieldNew.data();

  check_field_on_host(*bulkPtr, selectFieldOld, entityRankFieldOld, stkFieldOldData, oldValue);
  check_field_on_host(*bulkPtr, selectFieldNew, entityRankFieldNew, stkFieldNewData, newValue);

  const bool rotateNgpFieldViews = true;
  bulkPtr->update_field_data_states(rotateNgpFieldViews);

  {
    auto stkFieldOldData2 = stkFieldOld.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    auto stkFieldNewData2 = stkFieldNew.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

    check_field_on_device(ngpMesh, selectFieldOld, entityRankFieldOld, stkFieldOldData2, newValue);
    check_field_on_device(ngpMesh, selectFieldNew, entityRankFieldNew, stkFieldNewData2, oldValue);
  }
  //ENDNgpMultiStateFieldRotateOnDevice
}
}
