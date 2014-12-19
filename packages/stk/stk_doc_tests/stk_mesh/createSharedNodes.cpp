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
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_mesh/base/Comm.hpp"      // for comm_mesh_counts
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {
//BEGIN
TEST(stkMeshHowTo, createSharedNodes)
{
    const unsigned spatialDimension = 2;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());
    stk::mesh::Part &tri_part = metaData.declare_part_with_topology("tri_part", stk::topology::TRIANGLE_3_2D);
    metaData.commit();

    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    if (mesh.parallel_size() != 2) {
      return; //this test only runs on 2 procs
    }
    mesh.modification_begin();

    //  proc 0   proc 1
    //         2   2
    //elem 1  /|   |\ elem 2
    //       1 |   | 4
    //        \|   |/
    //         3   3
 
    const unsigned nodesPerElem = 3;
    stk::mesh::EntityId elemIds[] = {1, 2};//one elemId for each proc
    stk::mesh::EntityId elemNodeIds[][nodesPerElem] = { {1, 3, 2}, {4, 2, 3} };
    const int myproc = mesh.parallel_rank();
    int otherproc = 1;
    if (myproc == 1) otherproc = 0;

    stk::mesh::Entity elem = mesh.declare_entity(stk::topology::ELEM_RANK, elemIds[myproc], tri_part);
    stk::mesh::Entity elemNodes[nodesPerElem];
    elemNodes[0] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[myproc][0]);
    elemNodes[1] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[myproc][1]);
    elemNodes[2] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[myproc][2]);

    mesh.declare_relation(elem, elemNodes[0], 0);
    mesh.declare_relation(elem, elemNodes[1], 1);
    mesh.declare_relation(elem, elemNodes[2], 2);

    mesh.add_node_sharing(elemNodes[1], otherproc);
    mesh.add_node_sharing(elemNodes[2], otherproc);

    mesh.modification_end();

    //now verify there are 4 nodes globally, since nodes 2 and 3 are shared.
    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(mesh, entity_counts);
    const size_t expectedTotalNumNodes = 4;
    EXPECT_EQ(expectedTotalNumNodes, entity_counts[stk::topology::NODE_RANK]);
}
//END

//BEGIN_INDEP
TEST(stkMeshHowTo, createIndependentSharedNodes)
{
    const unsigned spatialDimension = 2;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());
    metaData.commit();

    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    if (mesh.parallel_size() != 2) {
      return; //this test only runs on 2 procs
    }
    mesh.modification_begin();

    //   proc 0  |   proc 1
    //         2 | 2
    //           |
    //       1   |   4
    //           |
    //         3 | 3

    const unsigned nodesPerProc = 3;
    stk::mesh::EntityId nodeIds[][nodesPerProc] = { {1, 3, 2}, {4, 2, 3} };
    const int myproc = mesh.parallel_rank();
    int otherproc = 1;
    if (myproc == 1) otherproc = 0;

    stk::mesh::Entity nodes[nodesPerProc];
    nodes[0] = mesh.declare_entity(stk::topology::NODE_RANK, nodeIds[myproc][0]);
    nodes[1] = mesh.declare_entity(stk::topology::NODE_RANK, nodeIds[myproc][1]);
    nodes[2] = mesh.declare_entity(stk::topology::NODE_RANK, nodeIds[myproc][2]);

    mesh.add_node_sharing(nodes[1], otherproc);
    mesh.add_node_sharing(nodes[2], otherproc);

    mesh.modification_end();

    //now verify there are 4 nodes globally, since nodes 2 and 3 are shared.
    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(mesh, entity_counts);
    const size_t expectedTotalNumNodes = 4;
    EXPECT_EQ(expectedTotalNumNodes, entity_counts[stk::topology::NODE_RANK]);
}
//END_INDEP

}
