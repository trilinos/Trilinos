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
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/Part.hpp"       // for Part

namespace
{
//-BEGIN
TEST(StkMeshHowTo, CreateSelectedFacesHex)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();

  // Generate a mesh containing 1 hex part and 6 shell parts
  const std::string generatedFileName = "generated:8x8x8|shell:xyzXYZ";
  stk::io::fill_mesh(generatedFileName, *bulkPtr);

  // ============================================================
  //+ EXAMPLE
  //+ Create a selector containing just the shell parts.
  stk::mesh::Selector shell_subset = bulkPtr->mesh_meta_data().get_topology_root_part(stk::topology::SHELL_QUAD_4);

  //+ Create the faces on just the selected shell parts.
  stk::mesh::create_all_sides(*bulkPtr, shell_subset);

  // ==================================================
  // VERIFICATION
  stk::mesh::Selector allEntities = bulkPtr->mesh_meta_data().universal_part();
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(allEntities, *bulkPtr, entityCounts);
  EXPECT_EQ( 896u, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ( 768u, entityCounts[stk::topology::FACE_RANK]);

  // Edges are not generated, only faces.
  EXPECT_EQ(0u,   entityCounts[stk::topology::EDGE_RANK]);
}
//-END
}
