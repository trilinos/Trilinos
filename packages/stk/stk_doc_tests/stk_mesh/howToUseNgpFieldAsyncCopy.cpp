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
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpExecutionSpace.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <Kokkos_Core.hpp>

namespace {

template <typename DoubleDataType, typename IntDataType>
void check_field_data_on_device(stk::mesh::BulkData& bulk, DoubleDataType& doubleFieldData, 
                                IntDataType& intFieldData, double expectedDoubleValue, int expectedIntValue)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::Selector selector = bulk.mesh_meta_data().universal_part();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, selector,
                                KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                {
                                  auto doubleFieldValues = doubleFieldData.entity_values(elem);
                                  NGP_EXPECT_EQ(expectedDoubleValue, doubleFieldValues());

                                  auto intFieldValues = intFieldData.entity_values(elem);
                                  NGP_EXPECT_EQ(expectedIntValue, intFieldValues());
                                });
}

TEST(stkMeshHowTo, ngpFieldAsyncCopy)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) > 1) { GTEST_SKIP();}

  using DoubleField = stk::mesh::Field<double>;
  using IntField = stk::mesh::Field<int>;

  const unsigned spatialDimension = 3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;

  unsigned numStates = 1;
  DoubleField& doubleField = meta.declare_field<double>(stk::topology::ELEM_RANK, "doubleField", numStates);
  IntField& intField = meta.declare_field<int>(stk::topology::ELEM_RANK, "intField", numStates);
  
  double initialDoubleFieldValue = 1.0;
  double modifiedDoubleFieldValue = initialDoubleFieldValue*2;
  int initialIntFieldValue = 2;
  int modifiedIntFieldValue = initialIntFieldValue*2;

  stk::mesh::put_field_on_entire_mesh_with_initial_value(doubleField, &initialDoubleFieldValue);
  stk::mesh::put_field_on_entire_mesh_with_initial_value(intField, &initialIntFieldValue);
  stk::io::fill_mesh("generated:1x1x1", bulk);

  stk::mesh::ExecSpaceWrapper<> execSpaceWithStream1 = stk::mesh::get_execution_space_with_stream();
  stk::mesh::ExecSpaceWrapper<> execSpaceWithStream2 = stk::mesh::get_execution_space_with_stream();

  {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
    auto doubleFieldData = doubleField.data<stk::mesh::ReadWrite>(execSpaceWithStream1.get_execution_space());
    auto doubleData = doubleFieldData.entity_values(elem);
    doubleData() = initialDoubleFieldValue*2;

    auto intFieldData = intField.data<stk::mesh::ReadWrite>(execSpaceWithStream2.get_execution_space());
    auto intData = intFieldData.entity_values(elem);
    intData() = initialIntFieldValue*2;
  }

  auto doubleFieldData = doubleField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(execSpaceWithStream1.get_execution_space());
  auto intFieldData = intField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(execSpaceWithStream2.get_execution_space());

  check_field_data_on_device(bulk, doubleFieldData, intFieldData, modifiedDoubleFieldValue, modifiedIntFieldValue);
}

}
