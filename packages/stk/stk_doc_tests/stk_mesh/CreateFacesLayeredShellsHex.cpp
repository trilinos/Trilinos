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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
namespace stk { namespace mesh { class BulkData; } }

namespace
{
//-BEGIN
TEST(StkMeshHowTo, CreateFacesLayeredShellsHex)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  stk::io::StkMeshIoBroker stkIo(communicator);

  // Generate a mesh containing 1 hex part and 12 shell parts
  // Shells are layered 2 deep.
  const std::string generatedFileName = "generated:8x8x8|shell:xxyyzzXYZXYZ";
  stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();

  // ============================================================
  //+ EXAMPLE
  //+ Create the faces
  stk::mesh::create_faces(stkIo.bulk_data());

  // ==================================================
  // VERIFICATION
  stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
  EXPECT_EQ(1280u, entityCounts[stk::topology::ELEMENT_RANK]);
  //+ The shell faces are the same as the boundary hex faces
  EXPECT_EQ(2112u, entityCounts[stk::topology::FACE_RANK]);

  // Edges are not generated, only faces.
  EXPECT_EQ(0u,   entityCounts[stk::topology::EDGE_RANK]);
}
//-END
}
