#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_io/FillMesh.hpp>
#include "NgpUnitTestUtils.hpp"

namespace {

using IntDualViewType = Kokkos::DualView<int*, stk::ngp::ExecSpace>;

void test_view_of_fields(const stk::mesh::BulkData& bulk,
                         stk::mesh::Field<double>& field1,
                         stk::mesh::Field<double>& field2)
{
  using FieldViewType = Kokkos::View<stk::mesh::NgpField<double>*,stk::ngp::MemSpace>;

  FieldViewType fields(Kokkos::ViewAllocateWithoutInitializing("fields"),2);
  FieldViewType::HostMirror hostFields = Kokkos::create_mirror_view(fields);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 2),
                       KOKKOS_LAMBDA(const unsigned& i)
                       {
                         new (&fields(i)) stk::mesh::NgpField<double>();
                       });

  hostFields(0) = stk::mesh::get_updated_ngp_field<double>(field1);
  hostFields(1) = stk::mesh::get_updated_ngp_field<double>(field2);

  Kokkos::deep_copy(fields, hostFields);

  unsigned numResults = 2;
  IntDualViewType result = ngp_unit_test_utils::create_dualview<IntDualViewType>("result",numResults);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 2),
                       KOKKOS_LAMBDA(const unsigned& i)
                       {
                         result.d_view(i) = fields(i).get_ordinal() == i ? 1 : 0;
                       });

  result.modify<IntDualViewType::execution_space>();
  result.sync<IntDualViewType::host_mirror_space>();

  EXPECT_EQ(1, result.h_view(0));
  EXPECT_EQ(1, result.h_view(1));
}

TEST(UnitTestNgp, viewOfFields)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  auto &field1 = meta.declare_field<double>(stk::topology::NODE_RANK, "myField1");
  auto &field2 = meta.declare_field<double>(stk::topology::NODE_RANK, "myField2");
  stk::mesh::put_field_on_mesh(field1, meta.universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(field2, meta.universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x4|sideset:zZ", *bulk);

  test_view_of_fields(*bulk, field1, field2);
}

}
