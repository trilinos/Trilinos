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
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace {

class NgpDebugFieldSync_PartialAllocation : public NgpDebugFieldSyncFixture
{
public:

};

TEST_F(NgpDebugFieldSync_PartialAllocation, FirstBlock_ScalarAccessUsingEntity_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
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






}
