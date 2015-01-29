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
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <algorithm>                    // for sort
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include <stk_mesh/base/FEMHelpers.hpp>
#include <unit_tests/BulkDataTester.hpp>

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

namespace
{

stk::mesh::Part& setupDavidNobleTestCase(stk::mesh::BulkData& bulk)
{
    //
    //        5____1  1  1____3
    //        |   /  /|\  \   |
    //        |E3/  / | \  \E1|
    //        | /  /E4|E2\  \ |
    //        |/  /___|___\  \|
    //        6   6   4   2   2
    //
    //        P2     P1      P0
    //

    stk::mesh::MetaData& meta = bulk.mesh_meta_data();

    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::TRIANGLE_3_2D);
    stk::mesh::Part& nonConformalPart = meta.declare_part("noconform", stk::topology::ELEMENT_RANK);

    meta.commit();

    bulk.modification_begin();

    const int nodesPerElem = 3;
    stk::mesh::EntityId elem1_nodes[nodesPerElem] = {1, 2, 3}; // 1
    stk::mesh::EntityId elem2_nodes[nodesPerElem] = {1, 4, 2}; // 2
    stk::mesh::EntityId elem3_nodes[nodesPerElem] = {6, 1, 5}; // 3
    stk::mesh::EntityId elem4_nodes[nodesPerElem] = {6, 4, 1}; // 4

    stk::mesh::EntityId elemId1 = 1; // p0
    stk::mesh::EntityId elemId2 = 2; // p1
    stk::mesh::EntityId elemId3 = 3; // p2
    stk::mesh::EntityId elemId4 = 4; // p1

    if(bulk.parallel_rank() == 0)
    {
        stk::mesh::declare_element(bulk, block_1, elemId1, elem1_nodes);
        stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
        stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
        bulk.add_node_sharing(node1, 1);
        bulk.add_node_sharing(node1, 2);
        bulk.add_node_sharing(node2, 1);
    }
    else if(bulk.parallel_rank() == 1)
    {
        stk::mesh::declare_element(bulk, block_1, elemId2, elem2_nodes);
        stk::mesh::declare_element(bulk, block_1, elemId4, elem4_nodes);

        stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
        stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
        stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
        bulk.add_node_sharing(node1, 2);
        bulk.add_node_sharing(node6, 2);

        bulk.add_node_sharing(node1, 0);
        bulk.add_node_sharing(node2, 0);
    }
    else
    {
        stk::mesh::declare_element(bulk, block_1, elemId3, elem3_nodes);
        stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
        stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
        bulk.add_node_sharing(node1, 0);
        bulk.add_node_sharing(node1, 1);
        bulk.add_node_sharing(node6, 1);
    }

    bulk.modification_end();

    return nonConformalPart;
}

bool isEntityInPart(stk::mesh::BulkData &bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id, const stk::mesh::Part &part)
{
    stk::mesh::Entity entity = bulk.get_entity(rank, id);
    const stk::mesh::PartVector &partsNodes6 = bulk.bucket(entity).supersets();
    bool isInPart = false;
    for(size_t i = 0; i < partsNodes6.size(); ++i)
    {
        if(partsNodes6[i] == &part)
        {
            isInPart = true;
            break;
        }
    }
    return isInPart;
}

TEST(BulkDataTest, testRemovingPartsOnNodeSharedWithOneProcAndAuraToAnotherProc)
{
    //unit test for ticket #12837

    stk::mesh::MetaData meta_data(2);
    stk::mesh::BulkData bulk(meta_data, MPI_COMM_WORLD);

    int num_procs = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if(num_procs == 3)
    {
        stk::mesh::Part& nonConformalPart = setupDavidNobleTestCase(bulk);

        {
            bulk.modification_begin();
            stk::mesh::PartVector add_parts;
            stk::mesh::PartVector rm_parts;
            add_parts.push_back(&nonConformalPart);
            if(bulk.parallel_rank() == 2)
            {
                stk::mesh::Entity element_3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
                bulk.change_entity_parts(element_3, add_parts, rm_parts);
            }
            bulk.modification_end();

            EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));
        }

        {
            bulk.modification_begin();
            if(bulk.parallel_rank() == 2)
            {
                stk::mesh::PartVector add_parts;
                stk::mesh::PartVector rm_parts;
                rm_parts.push_back(&nonConformalPart);
                stk::mesh::Entity element_3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
                bulk.change_entity_parts(element_3, add_parts, rm_parts);
            }

            EXPECT_TRUE(isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));

            bulk.modification_end();

            EXPECT_TRUE(!isEntityInPart(bulk, stk::topology::NODE_RANK, 6, nonConformalPart));

            stk::mesh::Entity node6 = bulk.get_entity(stk::topology::NODE_RANK, 6);
            if(bulk.parallel_rank() == 0)
            {
                EXPECT_TRUE(bulk.bucket(node6).in_aura());
            }
            else if(bulk.parallel_rank() == 1)
            {
                EXPECT_TRUE(bulk.bucket(node6).owned());
                EXPECT_TRUE(bulk.bucket(node6).shared());
            }
            else
            {
                EXPECT_TRUE(!bulk.bucket(node6).owned());
                EXPECT_TRUE(bulk.bucket(node6).shared());
            }
        }

    }
}

} // empty namespace

