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
#include <unistd.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_io/MeshField.hpp"

namespace
{

void write_five_steps_to_file(stk::mesh::BulkData& bulkData,
                              stk::mesh::FieldBase& nodeField,
                              double& time,
                              const std::string& ouputName,
                              stk::io::DatabasePurpose purpose)
{
  stk::io::StkMeshIoBroker stkIo(bulkData.parallel());
  stkIo.set_bulk_data(bulkData);
  size_t outputFileIndex = stkIo.create_output_mesh(ouputName, purpose);
  stkIo.add_field(outputFileIndex, nodeField);
  stkIo.write_output_mesh(outputFileIndex);

  const int numSteps = 5;
  for(int i = 0; i < numSteps; ++i) {
    stkIo.process_output_request(outputFileIndex, time);
    time += 1.0;
  }
}

void expect_ten_steps_in_file(const std::string& ouputName, MPI_Comm comm)
{
  stk::io::StkMeshIoBroker stkIo(comm);
  stkIo.add_mesh_database(ouputName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  EXPECT_EQ(10, stkIo.get_num_time_steps());
}

TEST(StkIoHowToAppend, toResultsFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                  .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::FieldBase& nodeField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodalField");
  stk::mesh::put_field_on_mesh(nodeField, meta.universal_part(), nullptr);
  stk::io::fill_mesh("generated:8x8x8|sideset:xX|nodeset:xX", *bulkPtr);

  std::string ouputName = "output.exo";
  double time = 1.0;
  write_five_steps_to_file(*bulkPtr, nodeField, time, ouputName, stk::io::WRITE_RESULTS);
  write_five_steps_to_file(*bulkPtr, nodeField, time, ouputName, stk::io::APPEND_RESULTS);
  expect_ten_steps_in_file(ouputName, bulkPtr->parallel());
  unlink(ouputName.c_str());
}

}
