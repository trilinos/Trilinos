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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <stddef.h>                     // for size_t, NULL
#include <unistd.h>                     // for unlink
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_util/parallel/Parallel.hpp"
namespace stk { namespace mesh { class Part; } }
namespace {

TEST(StkMeshIoBrokerHowTo, readMesh)
{
  std::string mesh_file_name = "input_mesh_example.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  {
    // ============================================================
    //BeginBasicReadWrite
    std::shared_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(communicator).create();

    //+ Create a basic mesh with a hex block, 3 shell blocks, 3 nodesets, and 3 sidesets.
    const std::string generatedFileName = "generated:8x8x8|shell:xyz|nodeset:xyz|sideset:XYZ";
    stk::io::fill_mesh(generatedFileName, *stkMesh);
    stk::io::write_mesh(mesh_file_name, *stkMesh);
    //EndBasicReadWrite
  }

  {
    //-BEGIN
    // ============================================================
    //+ EXAMPLE:
    //+ Read mesh data from the specified file.
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(mesh_file_name, stk::io::READ_MESH);

    //+ Creates meta data; creates parts
    stkIo.create_input_mesh();

    //+ Modifications to the meta data (such as creating extra parts and fields)
    //+ is generally be done here.

    //+ Commit the meta data and create the bulk data.
    //+ Populate the bulk data with data from the mesh file.
    stkIo.populate_bulk_data();

    // ============================================================
    //+ VERIFICATION
    //+ In this case we know that the mesh (specified above) contains
    //+ 4 element blocks, 3 nodesets, and 3 sidesets
    //+ There should be a STK Mesh Part for each of those.
    std::shared_ptr<const stk::mesh::MetaData> meta = stkIo.meta_data_ptr();
    EXPECT_NE(nullptr, meta->get_part("block_1"));
    EXPECT_NE(nullptr, meta->get_part("block_2"));
    EXPECT_NE(nullptr, meta->get_part("block_3"));
    EXPECT_NE(nullptr, meta->get_part("block_4"));

    EXPECT_NE(nullptr, meta->get_part("nodelist_1"));
    EXPECT_NE(nullptr, meta->get_part("nodelist_2"));
    EXPECT_NE(nullptr, meta->get_part("nodelist_3"));

    EXPECT_NE(nullptr, meta->get_part("surface_1"));
    EXPECT_NE(nullptr, meta->get_part("surface_2"));
    EXPECT_NE(nullptr, meta->get_part("surface_3"));

    std::shared_ptr<const stk::mesh::BulkData> bulk = stkIo.bulk_data_ptr();
    stk::mesh::EntityVector shellElems;
    stk::mesh::get_entities(*bulk, stk::topology::ELEM_RANK, *meta->get_part("block_2"), shellElems);
    EXPECT_EQ(64u, shellElems.size());
    //-END
  }
  // ============================================================
  // Cleanup
  unlink(mesh_file_name.c_str());
}
}
