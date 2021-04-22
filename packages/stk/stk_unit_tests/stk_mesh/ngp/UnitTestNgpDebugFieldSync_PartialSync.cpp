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
//#include "stk_mesh/base/Selector.hpp"

namespace {

class NgpDebugFieldSync_PartialSync : public NgpDebugFieldSyncFixture
{
};


TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncSelector_ReadAll_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_ModifySelectorClear_ReadAll_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.clear_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_ModifySelectorClearHost_ReadAll_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.clear_host_sync_state();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_ModifySelector_ReadAll_WarnAll)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteSelector_SyncSelector_ReadAll_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Vector_WriteSelector_SyncSelector_ReadAll_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  unsigned numComponents = 3;
  declare_vector_field<double>("doubleVectorField", numComponents, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double,stk::mesh::Cartesian> & stkField = initialized_vector_field<double>("doubleVectorField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_vector_field_on_device<double>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, ScalarPartial_WriteSelector_SyncSelector_ReadAll_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part3");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteSelector_SyncNothing_ReadAll_WarnSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, selector, 3.14);

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteOutsideSelector_SyncSelector_ReadAll_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, !selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Vector_WriteOutsideSelector_SyncSelector_ReadAll_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  unsigned numComponents = 3;
  declare_vector_field<double>("doubleVectorField", numComponents, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double,stk::mesh::Cartesian> & stkField = initialized_vector_field<double>("doubleVectorField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, !selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_vector_field_on_device<double>(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, ScalarPartial_WriteOutsideSelector_SyncSelector_ReadAll_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part3");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, !selector, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncSelector_ReadSelector_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField, selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Vector_WriteAll_SyncSelector_ReadSelector_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  unsigned numComponents = 3;
  declare_vector_field<double>("doubleVectorField", numComponents, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double,stk::mesh::Cartesian> & stkField = initialized_vector_field<double>("doubleVectorField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_vector_field_on_device<double>(stkField, selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, ScalarPartial_WriteAll_SyncSelector_ReadSelector_NoWarning)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part3");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField, selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncNothing_ReadSelector_WarnSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);

  read_scalar_field_on_device(stkField, selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncSelector_ReadOutsideSelector_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField, !selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=20.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, ScalarPartial_WriteAll_SyncSelector_ReadOutsideSelector_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector = *get_meta().get_part("Part3");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField, !selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=40.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Vector_WriteAll_SyncSelector_ReadOutsideSelector_WarnOutsideSelector)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2"});
  unsigned numComponents = 3;
  declare_vector_field<double>("doubleVectorField", numComponents, {"Part1", "Part2"});
  build_mesh({{"Part1", 2}, {"Part2", 2}});
  stk::mesh::Field<double,stk::mesh::Cartesian> & stkField = initialized_vector_field<double>("doubleVectorField");
  stk::mesh::Selector selector = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_vector_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector);
  stkField.sync_to_device();

  read_vector_field_on_device<double>(stkField, !selector);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=10.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=11.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=12.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[0]=20.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[1]=21.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleVectorField[2]=22.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncMultipleSelectors_ReadAll_WarnOutsideSelectors)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector1 = *get_meta().get_part("Part1");
  stk::mesh::Selector selector3 = *get_meta().get_part("Part3");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector1);
  stkField.modify_on_host(selector3);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=30.000000");
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=40.000000");
  check_no_warnings(stdoutString);
}

TEST_F(NgpDebugFieldSync_PartialSync, Scalar_WriteAll_SyncOverlappingSelectors_ReadAll_WarnOutsideSelectors)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  create_parts({"Part1", "Part2", "Part3"});
  declare_scalar_field<double>("doubleScalarField", {"Part1", "Part2", "Part3"});
  build_mesh({{"Part1", 2}, {"Part2", 2}, {"Part3", 2}});
  stk::mesh::Field<double> & stkField = initialized_field<double>("doubleScalarField");
  stk::mesh::Selector selector1 = *get_meta().get_part("Part1");
  stk::mesh::Selector selector2 = *get_meta().get_part("Part2");

  testing::internal::CaptureStdout();
  write_scalar_field_on_host_using_entity(stkField, 3.14);
  stkField.modify_on_host(selector1);
  stkField.modify_on_host(selector1 | selector2);
  stkField.sync_to_device();

  read_scalar_field_on_device(stkField);

  std::string stdoutString = testing::internal::GetCapturedStdout();
  extract_warning(stdoutString, 1, "WARNING: Accessing stale data on Device for Field doubleScalarField[0]=50.000000");
  check_no_warnings(stdoutString);
}


}
