// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector
namespace stk { namespace mesh { class BulkData; } }

namespace
{
  //-BEGIN
  TEST(StkMeshHowTo, CreateSelectedFacesHex)
  {
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) != 1) { return; }
    stk::io::StkMeshIoBroker stkIo(communicator);

    // Generate a mesh containing 1 hex part and 6 shell parts
    const std::string generatedFileName = "generated:8x8x8|shell:xyzXYZ";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();
    const stk::mesh::PartVector &all_parts = stkIo.meta_data().get_mesh_parts();

    // ============================================================
    //+ EXAMPLE
    //+ Create a selector containing just the shell parts.
    stk::mesh::Selector shell_subset;
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];
      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
	shell_subset |= *part;
      }
    }

    //+ Create the faces on just the selected shell parts.
    stk::mesh::create_faces(stkIo.bulk_data(), shell_subset);

    // ==================================================
    // VERIFICATION
    stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
    EXPECT_EQ( 896u, entityCounts[stk::topology::ELEMENT_RANK]);
    EXPECT_EQ( 768u, entityCounts[stk::topology::FACE_RANK]);

    // Edges are not generated, only faces.
    EXPECT_EQ(0u,   entityCounts[stk::topology::EDGE_RANK]);
  }
  //-END
}
