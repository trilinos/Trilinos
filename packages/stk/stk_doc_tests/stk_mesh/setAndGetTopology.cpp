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
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
namespace stk { namespace mesh { class Part; } }

namespace {
//-BEGIN
TEST(stkMeshHowTo, setAndGetTopology)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    //There are two methods of attaching topology to a part:

    //method 1: declare a part with a specified topology
    stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);

    stk::mesh::Part &hexPart = metaData.declare_part("hexElementPart");
    //method 2: set a topology on an existing part:
    stk::mesh::set_topology(hexPart, stk::topology::HEX_8);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId elem1Id = 1, elem2Id = 2;
    stk::mesh::Entity elem1=mesh.declare_entity(stk::topology::ELEMENT_RANK, elem1Id, tetPart);
    stk::mesh::Entity elem2=mesh.declare_entity(stk::topology::ELEMENT_RANK, elem2Id, hexPart);

    //No common nodes between elements - so topologically disconnected elem1 and elem2
    for(unsigned node_ord = 0 ; node_ord < 4; ++node_ord)
    {
      stk::mesh::Entity new_node = mesh.declare_entity(stk::topology::NODE_RANK, node_ord+100*elem1Id);
      mesh.declare_relation( elem1 , new_node , node_ord);
    }

    for(unsigned node_ord = 0 ; node_ord < 8; ++node_ord)
    {
      stk::mesh::Entity new_node2 = mesh.declare_entity(stk::topology::NODE_RANK, node_ord+100*elem2Id);
      mesh.declare_relation( elem2 , new_node2 , node_ord);
    }

    mesh.modification_end();

    stk::topology elem1_topology = mesh.bucket(elem1).topology();
    stk::topology elem2_topology = mesh.bucket(elem2).topology();

    EXPECT_EQ(stk::topology::TET_4, elem1_topology);
    EXPECT_EQ(stk::topology::HEX_8, elem2_topology);
}
//-END
}
