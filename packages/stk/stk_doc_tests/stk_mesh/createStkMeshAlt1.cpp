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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc

namespace
{

void create_example_exodus_file(MPI_Comm communicator, const std::string & exodusFileName);

//-BEGIN
TEST(StkMeshHowTo, CreateStkMesh)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  const std::string exodusFileName = "example.exo";

  create_example_exodus_file(communicator, exodusFileName);
  // Creation of STK Mesh objects.
  // MetaData creates the universal_part, locally-owned part, and globally shared part.
  std::shared_ptr<stk::mesh::BulkData> stkMeshBulkDataPtr = stk::mesh::MeshBuilder(communicator).create();
  stk::mesh::MetaData& stkMeshMetaData = stkMeshBulkDataPtr->mesh_meta_data();

  // Read the mesh data from the Exodus file and populate an STK Mesh.
  // The order of the following lines in {} are important
  {
    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    // Provide STK Mesh object to be populated
    exodusFileReader.set_bulk_data(*stkMeshBulkDataPtr);

    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

    // Populate the MetaData which has the descriptions of the Parts and Fields.
    exodusFileReader.create_input_mesh();

    // Populate entities in STK Mesh from Exodus file
    exodusFileReader.populate_bulk_data();
  }

  // Verify that the STK Mesh has 512 elements.
  stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, *stkMeshBulkDataPtr, entityCounts);
  EXPECT_EQ(512u, entityCounts[stk::topology::ELEMENT_RANK]);
  unlink(exodusFileName.c_str());
}
//-END

void create_example_exodus_file(MPI_Comm communicator, const std::string & exodusFileName)
{
  // ============================================================
  //+ INITIALIZATION:
  //+ Create a mesh
  stk::io::StkMeshIoBroker stkIo(communicator);

  const std::string generatedFileName = "generated:8x8x8";
  size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.set_active_mesh(index);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  size_t fh = stkIo.create_output_mesh(exodusFileName, stk::io::WRITE_RESULTS);
  stkIo.write_output_mesh(fh);
}
} // namespace
