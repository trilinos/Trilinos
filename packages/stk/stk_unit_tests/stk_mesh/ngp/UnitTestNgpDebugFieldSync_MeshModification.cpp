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
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/NgpFieldSyncDebugger.hpp>
#include <stk_util/stk_config.h>
#include "NgpDebugFieldSync_Fixtures.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <numeric>

namespace {

class NgpDebugFieldSync_MeshModification : public NgpDebugFieldSyncFixture {};

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_ProperlySyncToDevice_WithSimultaneousHostField_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::HostField<double> hostField(get_bulk(), stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_MissingAllModifySyncCallsToDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, UsingBucketIdAndOrdinal_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, UsingBucketIdAndOrdinal_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, UsingBucketIdAndOrdinal_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_bucket_id_and_ordinal(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_AccessOnDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_MissingDeviceFieldUpdate_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ModifyBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpFieldCopy, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_StaleDeviceFieldCopy_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ModifyBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  modify_element_part_membership({{2, "Part2", "Part1"}});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpFieldCopy, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_StaleDeviceFieldCopy_ClearSyncState_AccessOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> ngpFieldCopy = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});

  // The device Field is currently out-of-date, so our debugging code on the host side needs to not
  // mysteriously seg-fault before the user does the read on the Device side, where they will get
  // a useful warning.  Do a host-side write to confirm that we skip over dangerous code properly.
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_old_scalar_field_on_device(stkField, ngpFieldCopy);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element({{3, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element({2});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  delete_element({2});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 3, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  stkField.clear_sync_state();
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ModifyOnHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  delete_element({2});
  delete_element({3});

  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 4, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  modify_element_part_membership({{2, "Part2", "Part1"}});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  create_element({{3, "Part1"}}, stkField);
  create_element({{4, "Part1"}}, stkField);

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoConsecutiveMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();

  delete_element({2});
  delete_element({3});

  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_ProperlySyncToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.sync_to_device();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_ProperlySyncToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.modify_on_host();
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  stkField.clear_sync_state();
  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  modify_element_part_membership({{3, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  modify_element_part_membership({{2, "Part2", "Part1"}});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();
  const stk::mesh::EntityId maxIdToRead = 1;  // Avoid memory corruption due to accessing old Field after new bucket allocation

  create_element({{3, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  create_element({{4, "Part1"}}, stkField);
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField, maxIdToRead);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_MissingDeviceFieldUpdate_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::NgpField<double, NgpDebugger> & ngpField = stk::mesh::get_updated_ngp_field<double, NgpDebugger>(stkField);

  testing::internal::CaptureStdout();

  delete_element({2});
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  read_old_scalar_field_on_device(stkField, ngpField);

  delete_element({3});
  write_scalar_field_on_host_using_entity(stkField, 2.18);
  read_old_scalar_field_on_device(stkField, ngpField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Accessing un-updated Field doubleScalarField on Device after mesh modification");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_ProperlySyncToHost_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.sync_to_host();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_ModifyOnDevice_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.modify_on_device();
  stkField.clear_sync_state();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 8, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 7, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 5, "WARNING: Lost Device values for Field doubleScalarField due to a mesh modification before a sync to Host");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_ChangeBucket_ChangeBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  modify_element_part_membership({{3, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  modify_element_part_membership({{2, "Part2", "Part1"}});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_CreateBucket_CreateBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 1}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  create_element({{3, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  create_element({{4, "Part1"}}, stkField);
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MeshModification, TwoMods_DeleteBucket_DeleteBucket_MissingAllModifySyncCallsToHost_ClearSyncState_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  const unsigned bucketCapacity = 1;
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 3}, {"Part2", 1}}, bucketCapacity);
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");

  testing::internal::CaptureStdout();

  write_scalar_field_on_device(stkField, 3.14);
  stkField.clear_sync_state();
  delete_element({2});
  read_scalar_field_on_host_using_entity(stkField);

  write_scalar_field_on_device(stkField, 2.18);
  stkField.clear_sync_state();
  delete_element({3});
  read_scalar_field_on_host_using_entity(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

}
