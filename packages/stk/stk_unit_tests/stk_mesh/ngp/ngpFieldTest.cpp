#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
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
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_util/stk_config.h>

namespace ngp_field_test {

class NgpFieldFixture : public stk::unit_test_util::MeshFixture
{
public:
  template <typename T>
  stk::mesh::Field<T> & create_field(stk::topology::rank_t rank, const std::string & name)
  {
    unsigned numStates = 1;
    const T init = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, name, numStates);
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    return field;
  }
};

void move_data_between_fields_on_host(const stk::mesh::BulkData & bulk,
                                      const stk::mesh::Field<int>& source,
                                      stk::mesh::Field<int>& dest)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEMENT_RANK);
  stk::mesh::NgpField<int>& ngpSource = stk::mesh::get_updated_ngp_field<int>(source);
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

  stk::mesh::NgpField<int>& ngpDest = stk::mesh::get_updated_ngp_field<int>(dest);
  ngpDest.modify_on_host();
}

void test_field_values_on_device(stk::mesh::BulkData &bulk,
                                 const stk::mesh::Field<int> & stkField1,
                                 const stk::mesh::Field<int> & stkField2,
                                 const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField1 = stk::mesh::get_updated_ngp_field<int>(stkField1);
  stk::mesh::NgpField<int> & ngpField2 = stk::mesh::get_updated_ngp_field<int>(stkField2);
  ngpField1.sync_to_device();
  ngpField2.sync_to_device();

  const int component = 0;
  stk::mesh::Part& locallyOwned = bulk.mesh_meta_data().locally_owned_part();
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, locallyOwned,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   NGP_ThrowRequire(ngpField1(entity, component) == expectedFieldValue);
                                   NGP_ThrowRequire(ngpField2(entity, component) == expectedFieldValue);
                                 });
}

void test_field_values_on_host(const stk::mesh::BulkData &bulk,
                               const stk::mesh::Field<int>& stkField1,
                               const stk::mesh::Field<int>& stkField2,
                               const int expectedFieldValue)
{
  stk::mesh::NgpField<int> & ngpField1 = stk::mesh::get_updated_ngp_field<int>(stkField1);
  stk::mesh::NgpField<int> & ngpField2 = stk::mesh::get_updated_ngp_field<int>(stkField2);
  ngpField1.sync_to_host();
  ngpField2.sync_to_host();

  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEMENT_RANK);
  for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
  {
    const stk::mesh::Bucket &bucket = *buckets[iBucket];

    int* field1Data = reinterpret_cast<int*>(stk::mesh::field_data(stkField1, bucket));
    int* field2Data = reinterpret_cast<int*>(stk::mesh::field_data(stkField2, bucket));
    for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
    {
      EXPECT_EQ(expectedFieldValue, *field1Data);
      EXPECT_EQ(expectedFieldValue, *field2Data);
    }
  }
}

template <typename T>
void initialize_ngp_field(stk::mesh::Field<T> & stkField)
{
  stk::mesh::get_updated_ngp_field<T>(stkField);
}

template <typename T>
void modify_field_on_host(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & field, int multiplier)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(field.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(field, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      fieldData[iEntity] *= multiplier;
    }
  }
  field.modify_on_host();
}

void modify_mesh_part_membership(stk::mesh::BulkData & bulk, stk::mesh::Part & part)
{
  stk::mesh::PartVector addParts(1, &part);
  stk::mesh::EntityVector allElements;

  bulk.get_entities(stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), allElements);

  bulk.modification_begin();
  bulk.change_entity_parts(allElements[0], addParts);
  bulk.modification_end();
}

template <typename T>
void sync_field_to_device(stk::mesh::Field<T> & field)
{
  field.sync_to_device();
}

template <typename T>
void modify_field_on_device(stk::mesh::BulkData & bulk, stk::mesh::Field<T> & stkField, int multiplier)
{
  const int component = 0;
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, meta.locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
                                   ngpField(entity, component) *= multiplier;
                                 });
  stkField.modify_on_device();
}

template <typename T>
void sync_field_to_host(stk::mesh::Field<T> & stkField)
{
  stkField.sync_to_host();
}

template <typename T>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         const stk::mesh::Field<T>& stkField,
                         int expectedValue)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stkField.entity_rank());
  for (stk::mesh::Bucket * bucket : buckets) {
    T * fieldData = stk::mesh::field_data(stkField, *bucket);
    for(size_t iEntity = 0; iEntity < bucket->size(); ++iEntity) {
      EXPECT_EQ(fieldData[iEntity], expectedValue);
    }
  }
}

TEST_F(NgpFieldFixture, TestAriaAlgorithm)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & f1 = create_field<int>(stk::topology::ELEM_RANK, "f1");
  stk::mesh::Field<int> & f2 = create_field<int>(stk::topology::ELEM_RANK, "f2");

  setup_mesh("generated:1x1x1", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::NgpMesh ngpMesh = get_bulk().get_updated_ngp_mesh();

  int multiplier = 2;
  modify_field_on_device(get_bulk(), f1, multiplier);

  move_data_between_fields_on_host(get_bulk(), f1, f2);

  test_field_values_on_device(get_bulk(), f1, f2, multiplier);

  test_field_values_on_host(get_bulk(), f1, f2, multiplier);
}

TEST_F(NgpFieldFixture, GetNgpField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int> & ngpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  EXPECT_EQ(ngpIntField.get_ordinal(), stkIntField.mesh_meta_data_ordinal());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::DeviceField<int>&>::value));
#else
  EXPECT_TRUE((std::is_same<decltype(ngpIntField), stk::mesh::HostField<int>&>::value));
#endif
}

TEST_F(NgpFieldFixture, GetMultipleNgpFields)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField1 = create_field<int>(stk::topology::ELEM_RANK, "intField1");
  stk::mesh::Field<int> & stkIntField2 = create_field<int>(stk::topology::ELEM_RANK, "intField2");
  stk::mesh::Field<double> & stkDoubleField1 = create_field<double>(stk::topology::ELEM_RANK, "doubleField1");
  stk::mesh::Field<double> & stkDoubleField2 = create_field<double>(stk::topology::ELEM_RANK, "doubleField2");

  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::NgpField<int> & ngpIntField1 = stk::mesh::get_updated_ngp_field<int>(stkIntField1);
  stk::mesh::NgpField<int> & ngpIntField2 = stk::mesh::get_updated_ngp_field<int>(stkIntField2);
  stk::mesh::NgpField<double> & ngpDoubleField1 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField1);
  stk::mesh::NgpField<double> & ngpDoubleField2 = stk::mesh::get_updated_ngp_field<double>(stkDoubleField2);

  EXPECT_EQ(ngpIntField1.get_ordinal(), stkIntField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpIntField2.get_ordinal(), stkIntField2.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField1.get_ordinal(), stkDoubleField1.mesh_meta_data_ordinal());
  EXPECT_EQ(ngpDoubleField2.get_ordinal(), stkDoubleField2.mesh_meta_data_ordinal());
}

TEST_F(NgpFieldFixture, ModifyAndSync)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

#ifdef STK_USE_DEVICE_MESH
  size_t expectedSyncsToDevice = 2;
  size_t expectedSyncsToHost = 1;
#else
  size_t expectedSyncsToDevice = 0;
  size_t expectedSyncsToHost = 0;
#endif

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnHost)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_host(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);
  check_field_on_host(get_bulk(), stkIntField, multiplier);

  sync_field_to_device(stkIntField);
  modify_field_on_device(get_bulk(), stkIntField, multiplier);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier*multiplier);

#ifdef STK_USE_DEVICE_MESH
  size_t expectedSyncsToDevice = 2;
  size_t expectedSyncsToHost = 1;
#else
  size_t expectedSyncsToDevice = 0;
  size_t expectedSyncsToHost = 0;
#endif

  EXPECT_EQ(expectedSyncsToDevice, stkIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, stkIntField.num_syncs_to_host());

  stk::mesh::NgpField<int>& deviceNgpIntField = stk::mesh::get_updated_ngp_field<int>(stkIntField);
  stk::mesh::HostField<int> hostNgpIntField(get_bulk(), stkIntField);

  EXPECT_EQ(expectedSyncsToDevice, deviceNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToDevice, hostNgpIntField.num_syncs_to_device());
  EXPECT_EQ(expectedSyncsToHost, deviceNgpIntField.num_syncs_to_host());
  EXPECT_EQ(expectedSyncsToHost, hostNgpIntField.num_syncs_to_host());
}

TEST_F(NgpFieldFixture, UpdateNgpFieldAfterMeshMod_WithMostCurrentDataOnDevice)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  stk::mesh::Field<int> & stkIntField = create_field<int>(stk::topology::ELEM_RANK, "intField");
  stk::mesh::Part & dummyPart = get_meta().declare_part("DummyPart", stk::topology::ELEM_RANK);
  setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_ngp_field(stkIntField);

  int multiplier = 2;
  modify_field_on_device(get_bulk(), stkIntField, multiplier);
  modify_mesh_part_membership(get_bulk(), dummyPart);

  sync_field_to_host(stkIntField);
  check_field_on_host(get_bulk(), stkIntField, multiplier);
}

}
