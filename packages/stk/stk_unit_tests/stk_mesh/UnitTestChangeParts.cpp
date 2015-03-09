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

#include <gtest/gtest.h>
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "UnitTestUtils.hpp"

namespace
{

TEST(UnitTestChangeParts, test_throw_on_internal_part_change)
{
    const int spatialDim = 3;
    stk::mesh::MetaData metaData(spatialDim);
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::mesh::BulkData bulkData(metaData, communicator);

    std::string generatedMeshSpec = "generated:1x1x4";
    stk::mesh::unit_test::readMesh(generatedMeshSpec, bulkData, communicator);

    stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, 1);

    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts;

    addParts.push_back(&metaData.locally_owned_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    addParts.clear();
    addParts.push_back(&metaData.globally_shared_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    addParts.clear();
    removeParts.push_back(&metaData.locally_owned_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    removeParts.clear();
    removeParts.push_back(&metaData.globally_shared_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);
}

}
