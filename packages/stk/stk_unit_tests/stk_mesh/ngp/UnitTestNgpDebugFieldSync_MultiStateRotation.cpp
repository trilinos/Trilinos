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

class NgpDebugFieldSync_MultiStateRotation : public NgpDebugFieldSyncFixture
{
public:
  template <typename T>
  stk::mesh::Field<T>& get_initialized_field_state(const std::string& fieldName, const stk::mesh::FieldState& state)
  {
    stk::mesh::Field<T>& stkMultiStateField = *static_cast<stk::mesh::Field<T>*>(get_meta().get_field(stk::topology::ELEM_RANK, fieldName));
    stk::mesh::Field<T>& stkField = stkMultiStateField.field_of_state(state);
    fill_initial_field<T>(stkField);
    initialize_ngp_field<T>(stkField);
    return stkField;
  }

  void perform_field_state_rotation_correctly()
  {
    stk::mesh::sync_to_host_and_mark_modified(get_meta());
    get_bulk().update_field_data_states();
  }

  void perform_field_state_rotation_without_sync_to_host_mark_modified()
  {
    get_bulk().update_field_data_states();
  }
};

TEST_F(NgpDebugFieldSync_MultiStateRotation, HostToDevice_HappyPath)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkFieldOld, 3.14);
  write_scalar_field_on_host_using_entity(stkFieldNew, 6.28);

  perform_field_state_rotation_correctly();
  stkFieldOld.sync_to_device();
  stkFieldNew.sync_to_device();

  read_scalar_field_on_device(stkFieldOld);
  read_scalar_field_on_device(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, HostToDevice_noSync_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkFieldOld, 3.14);
  write_scalar_field_on_host_using_entity(stkFieldNew, 6.28);

  perform_field_state_rotation_correctly();

  read_scalar_field_on_device(stkFieldOld);
  read_scalar_field_on_device(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField_STKFS_OLD[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, HostToDevice_SyncStateOld_WarningForStateNew)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkFieldOld, 3.14);
  write_scalar_field_on_host_using_entity(stkFieldNew, 6.28);

  perform_field_state_rotation_correctly();
  stkFieldOld.sync_to_device();

  read_scalar_field_on_device(stkFieldOld);
  read_scalar_field_on_device(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, HostToDevice_noSyncToHostMarkModified_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkFieldOld, 3.14);
  write_scalar_field_on_host_using_entity(stkFieldNew, 6.28);

  perform_field_state_rotation_without_sync_to_host_mark_modified();
  stkFieldOld.sync_to_device();
  stkFieldNew.sync_to_device();

  read_scalar_field_on_device(stkFieldOld);
  read_scalar_field_on_device(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField_STKFS_OLD[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, HostToDevice_MeshMod_HappyPath)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  get_bulk().modification_begin();
  get_bulk().modification_end();

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkFieldOld, 3.14);
  write_scalar_field_on_host_using_entity(stkFieldNew, 6.28);

  perform_field_state_rotation_correctly();
  stkFieldOld.sync_to_device();
  stkFieldNew.sync_to_device();

  read_scalar_field_on_device(stkFieldOld);
  read_scalar_field_on_device(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, DeviceToHost_HappyPath)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkFieldOld, 3.14);
  write_scalar_field_on_device(stkFieldNew, 6.28);
  stkFieldOld.modify_on_device();
  stkFieldNew.modify_on_device();

  perform_field_state_rotation_correctly();

  read_scalar_field_on_host_using_entity(stkFieldOld);
  read_scalar_field_on_host_using_entity(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, DeviceToHost_noModifyOnDevice_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkFieldOld, 3.14);
  write_scalar_field_on_device(stkFieldNew, 6.28);

  perform_field_state_rotation_correctly();

  read_scalar_field_on_host_using_entity(stkFieldOld);
  read_scalar_field_on_host_using_entity(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField_STKFS_OLD[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField_STKFS_OLD[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, DeviceToHost_ModifyOnDeviceStateNew_WarnForStateOld)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkFieldOld, 3.14);
  write_scalar_field_on_device(stkFieldNew, 6.28);
  stkFieldNew.modify_on_device();

  perform_field_state_rotation_correctly();

  read_scalar_field_on_host_using_entity(stkFieldOld);
  read_scalar_field_on_host_using_entity(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 2, "WARNING: Accessing stale data on Host for Field doubleScalarField_STKFS_OLD[0]=6.28");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_MultiStateRotation, DeviceToHost_noSyncToHostMarkModified_Warning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  unsigned numStates = 2;
  declare_scalar_field<double>("doubleScalarField", {"Part1"}, numStates);
  build_mesh({{"Part1", 2}});

  stk::mesh::Field<double>& stkFieldOld = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateOld);
  stk::mesh::Field<double>& stkFieldNew = get_initialized_field_state<double>("doubleScalarField", stk::mesh::StateNew);

  testing::internal::CaptureStdout();
  write_scalar_field_on_device(stkFieldOld, 3.14);
  write_scalar_field_on_device(stkFieldNew, 6.28);
  stkFieldOld.modify_on_device();
  stkFieldNew.modify_on_device();

  perform_field_state_rotation_without_sync_to_host_mark_modified();
  stkFieldOld.sync_to_device();
  stkFieldNew.sync_to_device();

  read_scalar_field_on_host_using_entity(stkFieldOld);
  read_scalar_field_on_host_using_entity(stkFieldNew);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField_STKFS_OLD[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField_STKFS_OLD[0]=20");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=10");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Host for Field doubleScalarField[0]=20");
  check_no_warnings(stdoutString);
}

}
