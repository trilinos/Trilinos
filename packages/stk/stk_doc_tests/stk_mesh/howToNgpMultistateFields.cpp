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

  stkField.sync_to_host();

  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), select,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity entity) {
      *stk::mesh::field_data(stkField, entity) = fieldVal;
    });

  stkField.modify_on_host();
}

void set_field_on_device(const stk::mesh::NgpMesh& ngpMesh,
                         stk::mesh::NgpField<double>& ngpField,
                         double fieldVal)
{
  stk::mesh::Selector select(*ngpField.get_field_base());

  ngpField.sync_to_device();

  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), select,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
      ngpField(entityIndex, 0) = fieldVal;
    });

  ngpField.modify_on_device();
}

void check_field_on_host(const stk::mesh::BulkData & bulk,
                         stk::mesh::Field<double> & stkField,
                         double expectedFieldValue)
{
  stk::mesh::Selector select(stkField);
  stk::mesh::for_each_entity_run(bulk, stkField.entity_rank(), select,
    [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity entity) {
      EXPECT_NEAR(*stk::mesh::field_data(stkField, entity), expectedFieldValue, 1.e-9);
    });
}

void check_field_on_device(const stk::mesh::NgpMesh& ngpMesh,
                           stk::mesh::NgpField<double> & ngpField,
                           double expectedFieldValue)
{
  stk::mesh::Selector select(*ngpField.get_field_base());
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), select,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
      NGP_EXPECT_NEAR(ngpField(entityIndex, 0), expectedFieldValue, 1.e-9);
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

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;
  set_field_on_host(*bulkPtr, stkFieldOld, oldValue);
  set_field_on_host(*bulkPtr, stkFieldNew, newValue);

  bulkPtr->update_field_data_states();

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(stkFieldOld);
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(stkFieldNew);

  check_field_on_device(ngpMesh, ngpFieldOld, newValue);
  check_field_on_device(ngpMesh, ngpFieldNew, oldValue);
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

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;
  set_field_on_host(*bulkPtr, stkFieldOld, oldValue);
  set_field_on_host(*bulkPtr, stkFieldNew, newValue);

  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(stkFieldOld);
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(stkFieldNew);

  check_field_on_device(ngpMesh, ngpFieldOld, oldValue);
  check_field_on_device(ngpMesh, ngpFieldNew, newValue);

  stk::mesh::sync_to_host_and_mark_modified(meta);
  bulkPtr->update_field_data_states();

#ifdef STK_USE_DEVICE_MESH
  check_field_on_device(ngpMesh, ngpFieldOld, oldValue);
  check_field_on_device(ngpMesh, ngpFieldNew, newValue);
#else
  check_field_on_device(ngpMesh, ngpFieldOld, newValue);
  check_field_on_device(ngpMesh, ngpFieldNew, oldValue);
#endif

  ngpFieldOld.sync_to_device();
  ngpFieldNew.sync_to_device();

  check_field_on_device(ngpMesh, ngpFieldOld, newValue);
  check_field_on_device(ngpMesh, ngpFieldNew, oldValue);
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

  //BEGINNgpMultiStateFieldRotateOnDevice
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulkPtr);
  stk::mesh::NgpField<double>& ngpFieldOld = stk::mesh::get_updated_ngp_field<double>(stkFieldOld);
  stk::mesh::NgpField<double>& ngpFieldNew = stk::mesh::get_updated_ngp_field<double>(stkFieldNew);

  constexpr double oldValue = 1.0;
  constexpr double newValue = 2.0;

  set_field_on_device(ngpMesh, ngpFieldOld, oldValue);
  set_field_on_device(ngpMesh, ngpFieldNew, newValue);
#ifdef STK_USE_DEVICE_MESH
  check_field_on_host(*bulkPtr, stkFieldOld, 0.0);
  check_field_on_host(*bulkPtr, stkFieldNew, 0.0);
#else
  check_field_on_host(*bulkPtr, stkFieldOld, oldValue);
  check_field_on_host(*bulkPtr, stkFieldNew, newValue);
#endif

  const bool rotateNgpFieldViews = true;
  bulkPtr->update_field_data_states(rotateNgpFieldViews);

  check_field_on_device(ngpMesh, ngpFieldOld, newValue);
  check_field_on_device(ngpMesh, ngpFieldNew, oldValue);
  //ENDNgpMultiStateFieldRotateOnDevice
}

}
