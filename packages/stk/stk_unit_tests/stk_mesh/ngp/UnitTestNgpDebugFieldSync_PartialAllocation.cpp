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
#include "NgpDebugFieldSync_Fixtures.hpp"

namespace {

class NgpDebugFieldSync_PartialAllocation : public NgpDebugFieldSyncFixture
{
public:
  template <typename T>
  void write_scalar_field_on_device_using_entity_field_data(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::mesh::Selector fieldSelector(stkField);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, fieldSelector,
                                   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
#if defined(DEVICE_USE_LOCATION_BUILTINS)
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity);
#else
                                     stk::mesh::EntityFieldData<double> vals = ngpField(entity, __FILE__, __LINE__);
#endif
                                     vals[component] = value;
                                   });
  }

  template <typename T>
  void write_scalar_field_on_device_using_mesh_index(stk::mesh::Field<T> & stkField, T value)
  {
    const int component = 0;
    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
    stk::mesh::NgpField<T, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<T, NgpDebugger>(stkField);
    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stkField.entity_rank(), stkField);
    stk::mesh::EntityRank rank = ngpField.get_rank();

    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(unsigned ) {
                           for (unsigned i = 0; i < bucketIds.size(); ++i) {
                             const stk::mesh::NgpMesh::BucketType & bucket = ngpMesh.get_bucket(rank, bucketIds.device_get(i));
                             for (unsigned j = 0; j < bucket.size(); ++j) {
                               stk::mesh::NgpMesh::MeshIndex index{bucket.bucket_id(), static_cast<unsigned>(j)};
                               ngpField(index, component) = value;
                             }
                           }
                         });
  }
};

TEST_F(NgpDebugFieldSync_PartialAllocation, FirstBlock_ScalarAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntityFieldData_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntityFieldData_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_field_on_device_using_entity_field_data(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingMeshIndex_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_field_on_device_using_mesh_index(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingMeshIndex_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_field_on_device_using_mesh_index(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketId_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketId_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketIdAndOrdinal_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketIdAndOrdinal_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketIdAndOrdinalAndSize_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucketIdAndOrdinalAndSize_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_bucket_id_and_ordinal_and_size(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, FirstBlock_ScalarAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_DeviceEntityFieldDataAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device_using_entity_field_data(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_DeviceEntityFieldDataAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device_using_entity_field_data(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_DeviceMeshIndexAccess_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device_using_mesh_index(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_DeviceMeshIndexAccess_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device_using_mesh_index(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_bucket(stkField);
  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarWriteOnHost_ProperlyMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarWriteOnHost_ProperlyMarkAsModified_ClearHostSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_host_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarWriteOnDevice_ProperlyMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarWriteOnDevice_ProperlyMarkAsModified_ClearDeviceSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_device_sync_state();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarDeviceSetAll_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarDeviceSetAll_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarDeviceSetAll_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarDeviceSetAll_MissingAllModifySyncCallsToHost_AccessOnHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarDeviceSetAll_MissingAllModifySyncCallsToHost_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkField, 3.14);

  device_field_set_all(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=40");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_ProperlyMarkAsModified_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_entity(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.clear_sync_state();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  write_scalar_field_on_host_using_entity(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=40.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(stkField, 3.14+timeStep);
    stkField.modify_on_host();
    stkField.sync_to_device();
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_host_using_bucket(stkField, 3.14+timeStep);
    read_scalar_field_on_device(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);

  read_scalar_field_on_device(stkField);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  write_scalar_field_on_host_using_bucket(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_host_using_bucket(stkField, 3.14);
  write_scalar_field_on_host_using_bucket(stkField, 2.18);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=40.000000");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleTimestep_ProperlyMarkAsModified_ClearSyncState_AccessOnHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.clear_sync_state();
    read_scalar_field_on_host_using_entity(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_entity(stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleTimestep_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    stkField.modify_on_device();
    stkField.sync_to_host();
    read_scalar_field_on_host_using_bucket(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleTimestep_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  const size_t numTimeSteps = 2;
  for (size_t timeStep = 0; timeStep < numTimeSteps; ++timeStep) {
    write_scalar_field_on_device(stkField, 3.14+timeStep);
    read_scalar_field_on_host_using_bucket(stkField);
  }

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleAccesses_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);
  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleStaleAccesses_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  read_scalar_field_on_host_using_bucket(stkField);
  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleWrites_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingBucket_MultipleWrites_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  write_scalar_field_on_device(stkField, 2.18);

  read_scalar_field_on_host_using_bucket(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=30");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, EmptyField_MeshModification_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  get_meta().declare_field<double>(stk::topology::ELEM_RANK, "doubleScalarField", 1);
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  get_bulk().modification_begin();
  get_bulk().modification_end();

  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}


TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_MoveBucketOffField_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_MeshModification_Batch_Empty_MoveBucketOffField_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  batch_modify_element_part_membership({{3, "Part3", "Part2"}});
  get_bulk().modification_begin();
  get_bulk().modification_end();

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialAllocation, SecondBlock_ScalarAccessUsingEntity_DuringMeshModification_MoveBucketOffField_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 1}, {"Part2", 1}, {"Part3", 3}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership_with_scalar_field_write_using_entity({{2, "Part1", "Part2"}, {3, "Part1", "Part3"}}, stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

}
