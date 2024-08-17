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

#include <stddef.h>
#include <unistd.h>
#include <gtest/gtest.h>
#include <string>
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker


namespace {

size_t open_results_file(const std::string& resultsFileName,
                         stk::io::StkMeshIoBroker& stkIo,
                         std::shared_ptr<stk::mesh::BulkData> bulkPtr,
                         stk::mesh::Field<double>& elemField)
{
  size_t resultsFileIndex = stkIo.create_output_mesh(resultsFileName, stk::io::WRITE_RESULTS);
  stkIo.set_bulk_data(bulkPtr);
  stkIo.add_field(resultsFileIndex, elemField);
  stkIo.write_output_mesh(resultsFileIndex);
  return resultsFileIndex;
}

void remove_entity_from_mesh(stk::mesh::BulkData& bulk,
                             stk::mesh::EntityRank rank,
                             stk::mesh::EntityId entityId)
{
  bulk.modification_begin();
  stk::mesh::Entity entity = bulk.get_entity(rank, entityId);
  EXPECT_TRUE(bulk.is_valid(entity));
  bulk.destroy_entity(entity);
  bulk.modification_end();
}

void test_results_output(MPI_Comm comm,
                         const std::string& fileName,
                         int goldNumSteps)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                               .set_spatial_dimension(3).create();
  stk::io::StkMeshIoBroker stkIo;
  stk::io::fill_mesh_preexisting(stkIo, fileName, *bulkPtr);
  EXPECT_EQ(goldNumSteps, stkIo.get_num_time_steps());
  unlink(fileName.c_str());
}

TEST(StkIoResultsOutputMeshMod, writeResultsElemDelete)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                               .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Field<double>& elemField = meta.declare_field<double>(stk::topology::ELEM_RANK, "myElementField");
  stk::mesh::put_field_on_mesh(elemField, meta.universal_part(), 1, nullptr);

  const std::string meshDesc("generated:1x1x2|sideset:x");
  stk::io::fill_mesh(meshDesc, *bulkPtr);

  const std::string resultsFileName("results.e");
  const std::string resultsFile2Name("results.e-s0002");

  stk::io::StkMeshIoBroker stkIo(comm);
  size_t resultsFileIndex = open_results_file(resultsFileName, stkIo, bulkPtr, elemField);
  
  constexpr int nSteps = 5;
  for(int i = 1; i <= nSteps; ++i) {
    if (i == 3) {
      stk::mesh::EntityId elemId = 2;
      remove_entity_from_mesh(*bulkPtr, stk::topology::ELEM_RANK, elemId);

      stkIo.close_output_mesh(resultsFileIndex);
      resultsFileIndex = open_results_file(resultsFile2Name, stkIo, bulkPtr, elemField);
    }

    stkIo.process_output_request(resultsFileIndex, (double)(i-1));
  }

  stkIo.close_output_mesh(resultsFileIndex);

  constexpr int nSteps_before_meshMod = 2;
  constexpr int nSteps_after_meshMod = 3;
  test_results_output(comm, resultsFileName, nSteps_before_meshMod);
  test_results_output(comm, resultsFile2Name, nSteps_after_meshMod);
}

}

