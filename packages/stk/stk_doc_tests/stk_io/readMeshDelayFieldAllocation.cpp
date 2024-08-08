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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
namespace stk { namespace mesh { class Part; } }
namespace {

TEST(StkMeshIoBrokerHowTo, readMeshDelayFieldAllocation)
{
  std::string mesh_name = "input_mesh_example.e";
  MPI_Comm communicator = MPI_COMM_WORLD;

  {
    // ============================================================
    //+ INITIALIZATION:
    //+ Create a basic mesh with a hex block, 3 shell blocks, 3 nodesets, and 3 sidesets.
    stk::io::StkMeshIoBroker stkIo(communicator);

    const std::string generatedFileName = "generated:8x8x8|shell:xyz|nodeset:xyz|sideset:XYZ";
    size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);

    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t fh = stkIo.create_output_mesh(mesh_name, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);
  }

  {
    //-BEGIN
    // ============================================================
    //+ EXAMPLE:
    //+ Read mesh data from the specified file.
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);

    //+ Creates meta data; creates parts
    stkIo.create_input_mesh();

    //+ Any modifications to the meta data must be done here.
    //+ This includes declaring fields.

    //+ Commit the meta data and create the bulk data.
    //+ populate the bulk data with data from the mesh file.
    stkIo.populate_mesh();

    //+ Application would call mesh modification here.
    //+ for example, create_edges().

    //+ Mesh modifications complete, allocate field data.
    stkIo.populate_field_data();

    //-END
    // ============================================================
    //+ VERIFICATION
    //+ There should be:
    //+ 4 parts corresponding to the 1 hex block and 3 shell blocks
    stk::mesh::MetaData &meta = stkIo.meta_data();
    stk::mesh::Part *invalid = NULL;
    EXPECT_NE(invalid, meta.get_part("block_1"));
    EXPECT_NE(invalid, meta.get_part("block_2"));
    EXPECT_NE(invalid, meta.get_part("block_3"));
    EXPECT_NE(invalid, meta.get_part("block_4"));

    //+ 3 parts corresponding to the 3 nodesets.
    EXPECT_NE(invalid, meta.get_part("nodelist_1"));
    EXPECT_NE(invalid, meta.get_part("nodelist_2"));
    EXPECT_NE(invalid, meta.get_part("nodelist_3"));

    //+ 3 parts corresponding to the 3 sidesets.
    EXPECT_NE(invalid, meta.get_part("surface_1"));
    EXPECT_NE(invalid, meta.get_part("surface_2"));
    EXPECT_NE(invalid, meta.get_part("surface_3"));

  }
  // ============================================================
  // Cleanup
  unlink(mesh_name.c_str());
}
}
