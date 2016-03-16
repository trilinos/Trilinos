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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/SkinBoundary.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
namespace stk { namespace mesh { class BulkData; } }

namespace
{
//BEGIN_CREATE_EXPOSED_BOUNDARY
TEST(StkMeshHowTo, SkinExposedHex)
{
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) != 1) { return; }
    stk::io::StkMeshIoBroker stkIo(communicator);

    const std::string generatedFileName = "generated:1x1x1";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    // ============================================================
    //+ EXAMPLE
    //+ Skin the mesh and create the exposed boundary sides..
    stk::mesh::MetaData &metaData = stkIo.meta_data();
    stk::mesh::BulkData &bulkData = stkIo.bulk_data();
    stk::mesh::Selector allEntities = metaData.universal_part();
    stk::mesh::Part &skinPart = metaData.declare_part("skin", metaData.side_rank());
    stk::io::put_io_part_attribute(skinPart);

    stk::mesh::create_exposed_block_boundary_sides(bulkData, allEntities, {&skinPart});

    // ==================================================
    // VERIFICATION
    EXPECT_TRUE(stk::mesh::check_exposed_block_boundary_sides(bulkData, allEntities, skinPart));
    stk::mesh::Selector skin(skinPart & metaData.locally_owned_part());
    unsigned numSkinnedSides = stk::mesh::count_selected_entities(skin, bulkData.buckets(metaData.side_rank()));
    EXPECT_EQ(6u, numSkinnedSides) << "in part " << skinPart.name();
}
//END_CREATE_EXPOSED_BOUNDARY

//BEGIN_CREATE_INTERIOR_BOUNDARY
TEST(StkMeshHowTo, SkinInteriorHex)
{
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) != 1) { return; }
    stk::io::StkMeshIoBroker stkIo(communicator);

    const std::string generatedFileName = "generated:1x1x2";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    // ============================================================
    //+ EXAMPLE
    //+ Skin the mesh and create the exposed boundary sides..
    stk::mesh::MetaData &metaData = stkIo.meta_data();
    stk::mesh::BulkData &bulkData = stkIo.bulk_data();
    stk::mesh::Selector allEntities = metaData.universal_part();
    stk::mesh::Part &skinPart = metaData.declare_part("skin", metaData.side_rank());
    stk::io::put_io_part_attribute(skinPart);

    stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2u);
    stk::mesh::Part *block_1 = metaData.get_part("block_1");

    bulkData.modification_begin();
    stk::mesh::Part &block_2 = metaData.declare_part("block_2", stk::topology::ELEM_RANK);
    stk::io::put_io_part_attribute(block_2);
    bulkData.change_entity_parts(elem2, {&block_2}, {block_1});
    bulkData.modification_end();

    stk::mesh::create_interior_block_boundary_sides(bulkData, allEntities, {&skinPart});

    // ==================================================
    // VERIFICATION
    EXPECT_TRUE(stk::mesh::check_interior_block_boundary_sides(bulkData, allEntities, skinPart));
    stk::mesh::Selector skin(skinPart & metaData.locally_owned_part());
    unsigned numSkinnedSides = stk::mesh::count_selected_entities(skin, bulkData.buckets(metaData.side_rank()));
    EXPECT_EQ(1u, numSkinnedSides) << "in part " << skinPart.name();
}
//END_CREATE_INTERIOR_BOUNDARY
}
