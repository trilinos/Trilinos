// Copyright (c) 2014, Sandia Corporation.
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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <iostream>                     // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <string>                       // for string
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "UnitTestCEOCommonUtils.hpp"
#include "UnitTestCEO2Elem.hpp"
#include "UnitTestCEO3Elem.hpp"
#include "UnitTestCEO4ElemEdge.hpp"
#include "UnitTestCEO4ElemRotate.hpp"
#include "UnitTestCEO8Elem.hpp"
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for EntityProcVec, EntityProc, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Ghosting; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { namespace fixtures { class RingFixture; } } }

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace
{
//==============================================================================

TEST(CEO, change_entity_owner_2Elem2ProcMove)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(pm);
    const int p_size = stk::parallel_machine_size(pm);

    if(p_size != 2)
    {
        return;
    }

    const int spatial_dimension = 2;
    stk::mesh::MetaData meta(spatial_dimension);
    stk::unit_test_util::BulkDataTester bulk(meta, pm);

    stk::mesh::EntityVector elems;
    CEOUtils::fillMeshfor2Elem2ProcMoveAndTest(bulk, meta, elems);

    stk::mesh::EntityProcVec entity_procs;
    if(p_rank == 0)
    {
        entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
        entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
        entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
        entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
    }

    bulk.modification_begin("testing CEO");
    bulk.my_internal_change_entity_owner(entity_procs);

    CEOUtils::checkStatesAfterCEO_2Elem2ProcMove(bulk);
}

void test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);
    int prank = stk::parallel_machine_rank(communicator);
    if(psize == 3)
    { // Skip unless we're on 2 processors

        const int spatialDim = 3;
        std::vector<std::string> rankNames;
        rankNames.push_back("node");
        rankNames.push_back("edge");
        rankNames.push_back("face");
        rankNames.push_back("elem");
        rankNames.push_back("comst");

        stk::mesh::MetaData stkMeshMetaData(spatialDim, rankNames);
        //stk::mesh::Part &part = stkMeshMetaData.declare_part("constraints", stk::topology::CONSTRAINT_RANK);

        stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, communicator, autoAuraOption);
        const std::string generatedMeshSpecification = "generated:1x1x6";

        // STK IO module will be described in separate chapter.
        // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
        // The order of the following lines in {} are important
        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);

            // Inform STK IO which STK Mesh objects to populate later
            exodusFileReader.set_bulk_data(stkMeshBulkData);

            exodusFileReader.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);

            // Populate the MetaData which has the descriptions of the Parts and Fields.
            exodusFileReader.create_input_mesh();

            // Populate entities in STK Mesh from Exodus file
            exodusFileReader.populate_bulk_data();
        }

        /////////////////////////
        stkMeshBulkData.modification_begin();
        std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;
        stk::mesh::Ghosting &ghosting = stkMeshBulkData.create_ghosting("Ghost Node 1");
        if(prank == 0)
        {
            int proc2 = 2;
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            ghostingStruct.push_back(std::make_pair(node1, proc2));
        }
        stkMeshBulkData.change_ghosting(ghosting, ghostingStruct);
        stkMeshBulkData.modification_end();

        if (prank==2)
        {
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            ASSERT_TRUE(stkMeshBulkData.is_valid(node1));
            EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(node1));
        }
        /////////////////////////

        stkMeshBulkData.modification_begin();

        if (prank==2)
        {
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            stk::mesh::Entity constraint = stkMeshBulkData.declare_constraint(1);
            stkMeshBulkData.declare_relation(constraint,node1,0);
        }
        stk::mesh::fixup_ghosted_to_shared_nodes(stkMeshBulkData);
        stkMeshBulkData.modification_end();

        if(prank==2)
        {
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            EXPECT_TRUE(stkMeshBulkData.bucket(node1).shared());
        }

        if(prank==0 || prank==2)
        {
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            stk::mesh::EntityKey key = stkMeshBulkData.entity_key(node1);
            int otherProc = 2-prank;
            EXPECT_TRUE(stkMeshBulkData.in_ghost(ghosting, key, otherProc));
        }

        /////////////////////////

        stk::mesh::EntityProcVec entity_procs;

        if(prank == 0)
        {
            int proc2 = 2;
            stk::mesh::Entity elem1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
            entity_procs.push_back(stk::mesh::EntityProc(elem1, proc2));
            stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
            entity_procs.push_back(stk::mesh::EntityProc(node1, proc2));
        }

        stkMeshBulkData.change_entity_owner(entity_procs);

        if (prank == 2)
        {
            stk::mesh::Entity elem1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
            EXPECT_TRUE(stkMeshBulkData.parallel_owner_rank(elem1) == 2);
        }
    }
}

TEST(CEO, change_entity_owner_3Elem3Proc_WithCustomGhosts_WithAura)
{
    test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::AUTO_AURA);
}

TEST(CEO, change_entity_owner_3Elem3Proc_WithCustomGhosts_WithoutAura)
{
    test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::NO_AUTO_AURA);
}

}
