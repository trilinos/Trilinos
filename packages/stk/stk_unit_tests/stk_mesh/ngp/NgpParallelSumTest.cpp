#include <stk_ngp_test/ngp_test.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpFieldManager.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace  {

class NgpParallelSum : public stk::unit_test_util::MeshFixture
{
protected:
  NgpParallelSum() : stk::unit_test_util::MeshFixture(3)
  {
  }

  void initialize_shared_values(stk::mesh::FieldBase & userField, stk::mesh::FieldBase & goldValues,
                                bool leaveGoldValuesNotSummed=false)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().globally_shared_part());
    std::vector<int> sharingProcs;
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & node : *bucket) {
        stk::mesh::EntityKey nodeKey = get_bulk().entity_key(node);
        double id = static_cast<double>(get_bulk().identifier(node));
        double * gold = static_cast<double*>(stk::mesh::field_data(goldValues, node));
        *gold = id;

        double * user = static_cast<double*>(stk::mesh::field_data(userField, node));
        get_bulk().comm_procs(nodeKey, sharingProcs);
        const size_t numSharers = sharingProcs.size() + 1;
        *user = id / numSharers;
        if (leaveGoldValuesNotSummed) {
          *gold = *user;
        }
      }
    }
  }

};

class NgpCopyOwnedToShared : public stk::unit_test_util::MeshFixture
{
protected:
  NgpCopyOwnedToShared() : stk::unit_test_util::MeshFixture(3)
  {
  }

  void initialize_owned_shared_values(stk::mesh::FieldBase & userField, stk::mesh::FieldBase & goldValues,
                                      bool leaveGoldValuesZero=false)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().globally_shared_part());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & node : *bucket) {
        std::vector<int> sharingProcs;
        double id = static_cast<double>(get_bulk().identifier(node));
        double * gold = static_cast<double*>(stk::mesh::field_data(goldValues, node));
        *gold = id;

        double * user = static_cast<double*>(stk::mesh::field_data(userField, node));
        *user = bucket->owned() ? id : 0;
        if (leaveGoldValuesZero) {
          *gold = *user;
        }
      }
    }
  }

};

class NgpCommunicateFieldData : public stk::unit_test_util::MeshFixture
{
protected:
  NgpCommunicateFieldData() : stk::unit_test_util::MeshFixture(3)
  {
  }

  void initialize_owned_ghosted_values(stk::mesh::FieldBase & userField, stk::mesh::FieldBase & goldValues,
                                       bool leaveGoldValuesZero=false)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, (get_meta().locally_owned_part() & !get_meta().globally_shared_part()) | get_meta().aura_part());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & node : *bucket) {
        std::vector<int> sharingProcs;
        double id = static_cast<double>(get_bulk().identifier(node));
        double * gold = static_cast<double*>(stk::mesh::field_data(goldValues, node));
        *gold = id;

        double * user = static_cast<double*>(stk::mesh::field_data(userField, node));
        *user = bucket->owned() ? id : 0;
        if (leaveGoldValuesZero) {
          *gold = *user;
        }
      }
    }
  }

};

template <typename T>
void check_field_on_device(stk::mesh::NgpMesh &mesh,
                           stk::mesh::NgpField<T> & userField,
                           stk::mesh::NgpField<T> & goldValues)
{
  stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK, mesh.get_bulk_on_host().mesh_meta_data().universal_part(), KOKKOS_LAMBDA(stk::mesh::NgpMesh::MeshIndex entity)
  {
                                   NGP_EXPECT_NEAR(userField.get(entity, 0), goldValues.get(entity, 0), 1.e-12);
                                 });
}

NGP_TEST_F(NgpParallelSum, simpleVersion)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_shared_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  stk::mesh::parallel_sum<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}

#ifdef KOKKOS_ENABLE_CUDA
NGP_TEST_F(NgpParallelSum, simpleVersion_noSyncToDeviceAfterwards)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  const bool leaveSharedGoldValuesNotSummed = true;
  initialize_shared_values(userField, goldValues, leaveSharedGoldValuesNotSummed);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  const bool finalSyncToDevice = false;
  stk::mesh::parallel_sum<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField}, finalSyncToDevice);

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}
#endif

NGP_TEST_F(NgpCopyOwnedToShared, simpleVersion)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_owned_shared_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  stk::mesh::copy_owned_to_shared<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}

#ifdef KOKKOS_ENABLE_CUDA
NGP_TEST_F(NgpCopyOwnedToShared, simpleVersion_noSyncToDeviceAfterwards)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  const bool leaveSharedGoldValuesZero = true;
  initialize_owned_shared_values(userField, goldValues, leaveSharedGoldValuesZero);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  const bool finalSyncToDevice = false;
  stk::mesh::copy_owned_to_shared<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField}, finalSyncToDevice);

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}
#endif

NGP_TEST_F(NgpCommunicateFieldData, simpleVersion_takesGhosting)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_owned_ghosted_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  stk::mesh::communicate_field_data<double>(*get_bulk().ghostings()[stk::mesh::BulkData::AURA], std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}

#ifdef KOKKOS_ENABLE_CUDA
NGP_TEST_F(NgpCommunicateFieldData, simpleVersion_takesGhosting_noSyncToDeviceAfterwards)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  const bool leaveRecvGhostGoldValuesZero = true;
  initialize_owned_ghosted_values(userField, goldValues, leaveRecvGhostGoldValuesZero);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  const bool finalSyncToDevice = false;
  stk::mesh::communicate_field_data<double>(*get_bulk().ghostings()[stk::mesh::BulkData::AURA], std::vector<stk::mesh::NgpField<double>*>{&deviceUserField}, finalSyncToDevice);

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}
#endif

NGP_TEST_F(NgpCommunicateFieldData, simpleVersion_takesBulkData)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_owned_ghosted_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  stk::mesh::communicate_field_data<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}

#ifdef KOKKOS_ENABLE_CUDA
NGP_TEST_F(NgpCommunicateFieldData, simpleVersion_takesBulkData_noSyncToDeviceAfterwards)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  const bool leaveRecvGhostGoldValuesZero = true;
  initialize_owned_ghosted_values(userField, goldValues, leaveRecvGhostGoldValuesZero);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  const bool finalSyncToDevice = false;
  stk::mesh::communicate_field_data<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField}, finalSyncToDevice);

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}
#endif

NGP_TEST_F(NgpParallelSum, DeviceMPIVersion)
{
  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_shared_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  stk::mesh::parallel_sum_device_mpi<double>(ngpMesh, std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});

  check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
}

NGP_TEST_F(NgpParallelSum, Performance)
{
  const std::string serialMeshName = "serialParallelSumMesh.g";
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
  {
    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_SELF);

    std::string meshSpecDefault = "10x10x10";
    std::string meshSpec = stk::unit_test_util::get_command_line_option("-m", meshSpecDefault);

    stk::io::fill_mesh("generated:" + meshSpec, bulk);
    stk::io::write_mesh(serialMeshName, bulk);
  }

  const double initValue = 0.0;
  const int numStates = 1;
  stk::mesh::Field<double> & userField  = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "userField", numStates);
  stk::mesh::Field<double> & goldValues = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "goldValues", numStates);
  stk::mesh::put_field_on_mesh(userField, get_meta().universal_part(), &initValue);
  stk::mesh::put_field_on_mesh(goldValues, get_meta().universal_part(), &initValue);

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::fill_mesh_with_auto_decomp(serialMeshName, get_bulk());

  initialize_shared_values(userField, goldValues);

  stk::mesh::NgpMesh ngpMesh(get_bulk());
  stk::mesh::NgpFieldManager fieldManager(get_bulk());
  stk::mesh::NgpField<double> & deviceUserField = fieldManager.get_field<double>(userField.mesh_meta_data_ordinal());
  stk::mesh::NgpField<double> & deviceGoldValues = fieldManager.get_field<double>(goldValues.mesh_meta_data_ordinal());

  const bool useSimpleDefault = true;
  bool useSimple = stk::unit_test_util::get_command_line_option("-s", useSimpleDefault);

  const int numIterationsDefault = 1;
  int numIterations = stk::unit_test_util::get_command_line_option("-n", numIterationsDefault);

  for (int i = 0; i < numIterations; ++i) {
    if (useSimple) {
      const double startTime = stk::wall_time();
      stk::mesh::parallel_sum<double>(get_bulk(), std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});
      const double stopTime = stk::wall_time();
      const double localTime = stopTime - startTime;
      double globalTime = 0;
      stk::all_reduce_max(MPI_COMM_WORLD, &localTime, &globalTime, 1);

      if (get_bulk().parallel_rank() == 0) {
        std::cout << "Time for simple parallel_sum(): " << globalTime << " s" << std::endl;
      }
    }
    else {
      const double startTime = stk::wall_time();
      stk::mesh::parallel_sum_device_mpi<double>(ngpMesh, std::vector<stk::mesh::NgpField<double>*>{&deviceUserField});
      const double stopTime = stk::wall_time();
      const double localTime = stopTime - startTime;
      double globalTime = 0;
      stk::all_reduce_max(MPI_COMM_WORLD, &localTime, &globalTime, 1);

      if (get_bulk().parallel_rank() == 0) {
        std::cout << "Time for NGP-aware parallel_sum(): " << globalTime << " s" << std::endl;
      }
    }
  }

  if (numIterations == 1) {
    check_field_on_device<double>(ngpMesh, deviceUserField, deviceGoldValues);
  }

  unlink(serialMeshName.c_str());
}

}
