#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpFieldManager.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/stk_config.h>

namespace ngp_field_test {

class NgpFieldFixture : public stk::unit_test_util::MeshFixture {};

void modify_field_on_device(stk::mesh::NgpMesh &mesh,
                            const stk::mesh::NgpFieldManager& fieldManager,
                            const unsigned fieldOrdinal,
                            const int multiplier)
{
  const stk::mesh::BulkData& bulk = fieldManager.get_bulk();
  stk::mesh::NgpField<int> & ngpField = fieldManager.get_field<int>(fieldOrdinal);
  ngpField.sync_to_device();

  const int component = 0;
  stk::mesh::for_each_entity_run(mesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(),
                                 KOKKOS_LAMBDA(stk::mesh::NgpMesh::MeshIndex entity)
                                 {
                                   ngpField.get(entity, component) *= multiplier;
                                 });
  ngpField.modify_on_device();
}

void move_data_between_fields_on_host(const stk::mesh::NgpFieldManager& fieldManager,
                                      const stk::mesh::Field<int>& source,
                                      stk::mesh::Field<int>& dest)
{
  const stk::mesh::BucketVector& buckets = fieldManager.get_bulk().buckets(stk::topology::ELEMENT_RANK);
  stk::mesh::NgpField<int>& ngpSource = fieldManager.get_field<int>(source.mesh_meta_data_ordinal());
  ngpSource.sync_to_host();

  for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
  {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];

    int* sourceData = static_cast<int*>(stk::mesh::field_data(source, bucket));
    int* destData   = static_cast<int*>(stk::mesh::field_data(dest, bucket));
    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
    {
      *destData = *sourceData;
    }
  }

  stk::mesh::NgpField<int>& ngpDest = fieldManager.get_field<int>(dest.mesh_meta_data_ordinal());
  ngpDest.modify_on_host();
}

void test_field_values_on_device(const stk::mesh::NgpMesh &mesh,
                                 const stk::mesh::NgpFieldManager &fieldManager,
                                 const unsigned field1Ordinal,
                                 const unsigned field2Ordinal,
                                 const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField1 = fieldManager.get_field<int>(field1Ordinal);
  ngpField1.sync_to_device();

  stk::mesh::NgpField<int> & ngpField2 = fieldManager.get_field<int>(field2Ordinal);
  ngpField2.sync_to_device();

  const int component = 0;
  stk::mesh::Part& locallyOwned = fieldManager.get_bulk().mesh_meta_data().locally_owned_part();
  stk::mesh::for_each_entity_run(mesh, stk::topology::ELEM_RANK, locallyOwned,
                                 KOKKOS_LAMBDA(stk::mesh::NgpMesh::MeshIndex entity)
                                 {
                                   NGP_ThrowRequire(ngpField1.get(entity, component) == expectedFieldValue);
                                   NGP_ThrowRequire(ngpField2.get(entity, component) == expectedFieldValue);
                                 });
}

void test_field_values_on_host(const stk::mesh::NgpFieldManager &fieldManager,
                               const stk::mesh::Field<int>& field1,
                               const stk::mesh::Field<int>& field2,
                               const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField1 = fieldManager.get_field<int>(field1.mesh_meta_data_ordinal());
  ngpField1.sync_to_host();

  stk::mesh::NgpField<int> & ngpField2 = fieldManager.get_field<int>(field2.mesh_meta_data_ordinal());
  ngpField2.sync_to_host();

  const stk::mesh::BucketVector& buckets = fieldManager.get_bulk().buckets(stk::topology::ELEMENT_RANK);
  for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
  {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];

    int* field1Data = static_cast<int*>(stk::mesh::field_data(field1, bucket));
    int* field2Data = static_cast<int*>(stk::mesh::field_data(field2, bucket));
    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
    {
      EXPECT_EQ(expectedFieldValue, *field1Data);
      EXPECT_EQ(expectedFieldValue, *field2Data);
    }
  }
}

TEST_F(NgpFieldFixture, TestAriaAlgorithm)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  unsigned numStates = 1;

  const int init = 1;
  stk::mesh::Field<int> &f1 = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, "f1", numStates);
  stk::mesh::put_field_on_mesh(f1, get_meta().universal_part(), &init);

  stk::mesh::Field<int> &f2 = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, "f2", numStates);
  stk::mesh::put_field_on_mesh(f2, get_meta().universal_part(), &init);

  setup_mesh("generated:1x1x1", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());

  int multiplier = 2;
  modify_field_on_device(ngpMesh, fieldManager, f1.mesh_meta_data_ordinal(), multiplier);

  move_data_between_fields_on_host(fieldManager, f1, f2);

  test_field_values_on_device(ngpMesh,
                              fieldManager,
                              f1.mesh_meta_data_ordinal(),
                              f2.mesh_meta_data_ordinal(),
                              multiplier);

  test_field_values_on_host(fieldManager, f1, f2, multiplier);
}
}
