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
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>  // for test_shift_ring
#include <unit_tests/Setup8Quad4ProcMesh.hpp>
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for stk::mesh::EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include "UnitTestCEOCommonUtils.hpp"

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
using stk::mesh::BaseEntityRank;
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

void update_mesh_based_on_changed_sharing_info(BulkDataTester& bulk, std::vector<stk::mesh::Entity>& modifiedEntities)
{
    bulk.my_resolve_ownership_of_modified_entities(modifiedEntities);
    bulk.my_move_entities_to_proper_part_ownership(modifiedEntities);
    bulk.my_update_comm_list_based_on_changes_in_comm_map();
    bulk.my_update_comm_list(modifiedEntities);
}

TEST(CEOME, change_entity_owner_2Elem2ProcMove)
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
    BulkDataTester bulk(meta, pm);

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

    bulk.change_entity_owner_exp(entity_procs);

    CEOUtils::checkStatesAfterCEO_2Elem2ProcMove(bulk);

    ////////////////////////////////////////////////////////////////////////////

    //   id/owner_proc
    //
    //   1/0---4/0---5/0      1/0---4/1---5/1
    //    |     |     |        |     |     |
    //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
    //    |     |     |        |     |     |
    //   2/0---3/0---6/0      2/0---3/0---6/1

    std::vector<bool> nodeSharing;
    if(p_rank == 0)
    {
        bool isNodeShared[] = {false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + 6);
    }
    else
    {
        bool isNodeShared[] = {false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + 6);
    }

    int otherProc = 1;
    if(p_rank == 1)
    {
        otherProc = 0;
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(bulk, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, i + 1);
            CEOUtils::addSharingInfo(bulk, node, stk::mesh::BulkData::SHARED, otherProc);
            modifiedEntities.push_back(node);
        }
    }

    update_mesh_based_on_changed_sharing_info(bulk, modifiedEntities);

    bulk.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    bulk.modification_begin();
    bulk.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_2Elem2ProcMove(bulk);
}

TEST(CEOME, change_entity_owner_2Elem2ProcFlip)
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
    BulkDataTester mesh(meta, pm);

    CEOUtils::fillMeshfor2Elem2ProcFlipAndTest(mesh, meta);

    //okay now flip
    //    1/0---4/0---5/1          1/1---4/0---5/0
    //     |     |     |            |     |     |
    //     | 1/0 | 2/1 |       =>   | 1/1 | 2/0 |
    //     |     |     |            |     |     |
    //    2/0---3/0---6/1          2/1---3/0---6/0

    stk::mesh::EntityProcVec entity_procs_flip;
    if(p_rank == 0)
    {
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEMENT_RANK, 1), 1));
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 1), 1));
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 2), 1));
    }
    else
    {
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEMENT_RANK, 2), 0));
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 5), 0));
        entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 6), 0));
    }
    mesh.change_entity_owner_exp(entity_procs_flip);

    ////////////////////////////////////////////////////////////////////////////

    std::vector<bool> nodeSharing;
    if(p_rank == 0)
    {
        bool isNodeShared[] = {false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + 6);
    }
    else
    {
        bool isNodeShared[] = {false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + 6);
    }

    int otherProc = 1;
    if(p_rank == 1)
    {
        otherProc = 0;
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);
            CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
            modifiedEntities.push_back(node);
        }
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////
    CEOUtils::checkStatesAfterCEOME_2Elem2ProcFlip(mesh);
}

TEST(CEOME, change_entity_owner_3Elem2ProcMoveRight)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    BulkDataTester mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    if(p_size != 2)
    {
        return;
    }

    EntityVector nodes;
    EntityVector elements;

    CEOUtils::fillMeshfor3Elem2ProcMoveRightAndTest(mesh, meta_data, nodes, elements);

    std::vector<EntityProc> change;
    if(p_rank == 0)
    {
        change.push_back(EntityProc(elements[1], 1));
        change.push_back(EntityProc(nodes[4], 1));
        change.push_back(EntityProc(nodes[5], 1));
    }

    mesh.change_entity_owner_exp(change);

    CEOUtils::checkStatesAfterCEO_3Elem2ProcMoveRight(mesh);
    ////////////////////////////////////////////////////////////////////////////

    //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    size_t numNodes = 8;
    std::vector<bool> nodeSharing;

    if(p_rank == 0)
    {
        bool isNodeShared[] = {false, false, true, true, false, false, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
    }
    else
    {
        bool isNodeShared[] = {false, false, true, true, false, false, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    int otherProc = 1;
    if(p_rank == 1)
    {
        otherProc = 0;
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);

        if(CEOUtils::isEntityInSharingCommMap(mesh, node) && !nodeSharing[i])
        {
            modifiedEntities.push_back(node);
        }

        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
            modifiedEntities.push_back(node);
        }
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_3Elem2ProcMoveRight(mesh);
}

TEST(CEOME, change_entity_owner_3Elem2ProcMoveLeft)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    BulkDataTester mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    if(p_size != 2)
    {
        return;
    }

    EntityVector nodes;
    EntityVector elements;

    CEOUtils::fillMeshfor3Elem2ProcMoveLeftAndTest(mesh, meta_data, nodes, elements);

    std::vector<EntityProc> change;
    if(p_rank == 1)
    {
        change.push_back(EntityProc(elements[0], 0));
        change.push_back(EntityProc(nodes[2], 0));
        change.push_back(EntityProc(nodes[3], 0));
    }

    mesh.change_entity_owner_exp(change);

    CEOUtils::checkStatesAfterCEO_3Elem2ProcMoveLeft(mesh);

    ////////////////////////////////////////////////////////////////////////////

    //   id/owner_proc
    //
    //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

    size_t numNodes = 8;
    std::vector<bool> nodeSharing;

    if(p_rank == 0)
    {
        bool isNodeShared[] = {false, false, false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
    }
    else
    {
        bool isNodeShared[] = {false, false, false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    int otherProc = 1;
    if(p_rank == 1)
    {
        otherProc = 0;
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);

        if(CEOUtils::isEntityInSharingCommMap(mesh, node) && !nodeSharing[i])
        {
            modifiedEntities.push_back(node);
        }

        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, otherProc);
            modifiedEntities.push_back(node);
        }
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////
    CEOUtils::checkStatesAfterCEOME_3Elem2ProcMoveLeft(mesh);
}

TEST(CEOME, change_entity_owner_4Elem4ProcEdge)
{
    // This unit-test is designed to test the conditions that results that
    // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
    // it will test the changing-of-ownership of a shared edge to a proc that
    // either ghosted it or did not know about it.
    //
    //         id/proc                             id/proc
    //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
    //        |      |     |    {|     |          |      |     |    {|     |
    //        | 1/0  | 2/1 | 3/2{| 4/3 |          | 1/0  | 2/1 | 3/0{| 4/3 |
    //        |      |     |    {|     |          |      |     |    {|     |
    //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
    //  this edge moves to p0 --^
    //  element 3 moves to proc 0.
    //  nodes 7&8 move to proc 0.
    //  proc 2 forgets everything.
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. We will take the edge shared by the last
    // two (rightmost) elements and change the ownership to proc 0.

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    BulkDataTester mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    if(p_size != 4)
    {
        return;
    }

    stk::mesh::EntityKey elem_key_chg_own;
    CEOUtils::fillMeshfor4Elem4ProcEdgeAndTest(mesh, meta_data, elem_key_chg_own);

    std::vector<EntityProc> change;
    if(p_rank == 2)
    {
        // Change ownership of changing elem and all entities in it's closure that
        // we own to proc 0.

        Entity changing_elem = mesh.get_entity(elem_key_chg_own);
        ASSERT_TRUE( mesh.is_valid(changing_elem));
        EntityProc eproc(changing_elem, 0 /*new owner*/);
        change.push_back(eproc);

        const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
                for(stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
        {
            stk::mesh::Entity const *to_i = mesh.begin(changing_elem, irank);
            stk::mesh::Entity const *to_e = mesh.end(changing_elem, irank);
            for (; to_i != to_e; ++to_i)
            {
                if (mesh.parallel_owner_rank(*to_i) == p_rank)
                {
                    EntityProc eproc_new(*to_i, 0 /*new owner*/);
                    change.push_back(eproc_new);
                }
            }
        }
    }

    stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_TRUE(mesh.is_valid(node)) << " is not valid on processor " << p_rank;

    mesh.change_entity_owner_exp(change);

    ////////////////////////////////////////////////////////////////////////////

    //
    //         id/proc                             id/proc
    //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
    //        |      |     |    {|     |          |      |     |    {|     |
    //        | 1/0  | 2/1 | 3/2{| 4/3 |          | 1/0  | 2/1 | 3/0{| 4/3 |
    //        |      |     |    {|     |          |      |     |    {|     |
    //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3

    size_t numNodes = 10;
    std::vector<bool> nodeSharing;
    std::vector<int> nodeSharingProcs;

    size_t numEdges = 1;
    std::vector<bool> edgeSharing;
    std::vector<int> edgeSharingProcs;

    std::vector<stk::mesh::Entity> modifiedEntities;

    if(p_rank == 0)
    {
        bool isNodeShared[] = {false, false, true, true, true, true, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
        int procIdsNode[] = {-1, -1, 1, 1, 1, 1, 3, 3, -1, -1};
        nodeSharingProcs.assign(procIdsNode, procIdsNode + numNodes);
        bool isEdgeShared[] = {true};
        int procIdsEdge[] = {3};
        edgeSharing.assign(isEdgeShared, isEdgeShared + numEdges);
        edgeSharingProcs.assign(procIdsEdge, procIdsEdge + numEdges);
    }
    else if(p_rank == 1)
    {
        bool isNodeShared[] = {false, false, true, true, true, true, false, false, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
        int procIdsNode[] = {-1, -1, 0, 0, 0, 0, -1, -1, -1, -1};
        nodeSharingProcs.assign(procIdsNode, procIdsNode + numNodes);
        bool isEdgeShared[] = {false};
        int procIdsEdge[] = {-1};
        edgeSharing.assign(isEdgeShared, isEdgeShared + numEdges);
        edgeSharingProcs.assign(procIdsEdge, procIdsEdge + numEdges);
    }
    else if(p_rank == 2)
    {
        bool isNodeShared[] = {false, false, false, false, false, false, false, false, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
        int procIdsNode[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        nodeSharingProcs.assign(procIdsNode, procIdsNode + numNodes);
        bool isEdgeShared[] = {false};
        int procIdsEdge[] = {-1};
        edgeSharing.assign(isEdgeShared, isEdgeShared + numEdges);
        edgeSharingProcs.assign(procIdsEdge, procIdsEdge + numEdges);
    }
    else
    {
        bool isNodeShared[] = {false, false, false, false, false, false, true, true, false, false};
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);
        int procIdsNode[] = {-1, -1, -1, -1, -1, -1, 0, 0, -1, -1};
        nodeSharingProcs.assign(procIdsNode, procIdsNode + numNodes);
        bool isEdgeShared[] = {true};
        int procIdsEdge[] = {0};
        edgeSharing.assign(isEdgeShared, isEdgeShared + numEdges);
        edgeSharingProcs.assign(procIdsEdge, procIdsEdge + numEdges);
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);
        if(CEOUtils::isEntityInSharingCommMap(mesh, node) && !nodeSharing[i])
        {
            modifiedEntities.push_back(node);
        }

        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeSharingProcs[i]);
            modifiedEntities.push_back(node);
        }
    }

    for(size_t i = 0; i < edgeSharing.size(); i++)
    {
        stk::mesh::Entity edge = mesh.get_entity(stk::topology::EDGE_RANK, i + 1);

        if(CEOUtils::isEntityInSharingCommMap(mesh, edge) && !edgeSharing[i])
        {
            modifiedEntities.push_back(edge);
        }

        stk::mesh::EntityKey key(stk::topology::EDGE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(edgeSharing[i])
        {
            CEOUtils::addSharingInfo(mesh, edge, stk::mesh::BulkData::SHARED, edgeSharingProcs[i]);
            modifiedEntities.push_back(edge);
        }
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////
    CEOUtils::checkStatesAfterCEOME_4Elem4ProcEdge(mesh);
}

TEST(CEOME, change_entity_owner_8Elem4ProcMoveTop)
{
    //
    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
    //
    // This test moves ownership of elements 6 and 7 (as well as their locally-owned
    // nodes) to procs 3 and 0, respectively.
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if(numProcs != 4)
    {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    BulkDataTester mesh(meta, pm);
    const int p_rank = mesh.parallel_rank();

    CEOUtils::fillMeshfor8Elem4ProcMoveTopAndTest(mesh, meta);

    std::vector<stk::mesh::EntityProc> entities_to_move;
    if(mesh.parallel_rank() == 1)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 6);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    if(mesh.parallel_rank() == 2)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 7);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }

    mesh.change_entity_owner_exp(entities_to_move);

    ////////////////////////////////////////////////////////////////////////////

    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3

    size_t numNodes = 15;
    std::vector<bool> nodeSharing;
    std::vector<std::vector<int> > nodeSharingProcs(15);

    std::vector<stk::mesh::Entity> modifiedEntities;

    if(p_rank == 0)
    {
        bool isNodeShared[] = {
                false, true, false, false, false,
                false, true, true, true, false,
                false, true, true, true, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // node 11: nothing
        // nodes 12, 13, 14
        {
            int procsForNodeSharing[] = {3};
            nodeSharingProcs[12 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[13 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[14 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 15
        // node 6
        // node 7
        {
            int procsForNodeSharing[] = {1, 3};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {1, 2, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {2, 3};
            nodeSharingProcs[9 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 10
        // node 1
        // node 2
        {
            int procsForNodeSharing[] = {1};
            nodeSharingProcs[2 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // nodes 3, 4, 5
    }
    else if(p_rank == 1)
    {
        bool isNodeShared[] = {
                false, true, true, false, false,
                false, true, true, false, false,
                false, false, false, false, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 11-15: nothing
        // nodes 6, 9, 10: nothing
        // node 7
        {
            int procsForNodeSharing[] = {0, 3};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {0, 2, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // nodes 1, 4, 5: nothing
        // node 2
        {
            int procsForNodeSharing[] = {0};
            nodeSharingProcs[2 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 3
        {
            int procsForNodeSharing[] = {2};
            nodeSharingProcs[3 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }
    else if(p_rank == 2)
    {
        bool isNodeShared[] = {
                false, false, true, true, false,
                false, false, true, true, false,
                false, false, false, false, false
        };

        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 11-15: nothing
        // nodes 6, 7, 10: nothing
        // node 8
        {
            int procsForNodeSharing[] = {0, 1, 3};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {0, 3};
            nodeSharingProcs[9 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // nodes 1, 2, 5: nothing
        // node 3
        {
            int procsForNodeSharing[] = {1};
            nodeSharingProcs[3 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 4
        {
            int procsForNodeSharing[] = {3};
            nodeSharingProcs[4 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }
    else
    {
        bool isNodeShared[] = {
                false, false, false, true, false,
                false, true, true, true, false,
                false, true, true, true, false
        };
        nodeSharing.assign(isNodeShared, isNodeShared + numNodes);

        // nodes 1, 2, 3, 5, 6, 10, 11, 15: nothing
        // nodes 12, 13, 14
        {
            int procsForNodeSharing[] = {0};
            nodeSharingProcs[12 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[13 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
            nodeSharingProcs[14 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
        // node 7
        {
            int procsForNodeSharing[] = {0, 1};
            nodeSharingProcs[7 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 8
        {
            int procsForNodeSharing[] = {0, 1, 2};
            nodeSharingProcs[8 - 1].assign(procsForNodeSharing, procsForNodeSharing + 3);
        }
        // node 9
        {
            int procsForNodeSharing[] = {0, 2};
            nodeSharingProcs[9 - 1].assign(procsForNodeSharing, procsForNodeSharing + 2);
        }
        // node 4
        {
            int procsForNodeSharing[] = {2};
            nodeSharingProcs[4 - 1].assign(procsForNodeSharing, procsForNodeSharing + 1);
        }
    }

    EXPECT_EQ(numNodes, nodeSharing.size());

    for(size_t i = 0; i < nodeSharing.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i + 1);
        if(CEOUtils::isEntityInSharingCommMap(mesh, node) && !nodeSharing[i])
        {
            modifiedEntities.push_back(node);
        }

        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

        if(nodeSharing[i])
        {
            for(size_t j = 0; j < nodeSharingProcs[i].size(); j++)
            {
                CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeSharingProcs[i][j]);
            }
            modifiedEntities.push_back(node);
        }
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////
    CEOUtils::checkStatesAfterCEOME_8Elem4ProcMoveTop(mesh);
}

TEST(CEOME, change_entity_owner_4Elem4ProcRotate)
{
    //
    //     id/proc                id/proc
    //      7/3---8/2---9/2        7/2---8/1---9/1
    //       |     |     |          |     |     |
    //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
    //       |     |     |          |     |     |
    //      4/0---5/0---6/1  -->   4/3---5/3---6/0
    //       |     |     |          |     |     |
    //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
    //       |     |     |          |     |     |
    //      1/0---2/0---3/1        1/3---2/3---3/0

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if(numProcs != 4)
    {
        return;
    }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    BulkDataTester mesh(meta, pm);
    const int p_rank = mesh.parallel_rank();
    CEOUtils::fillMeshfor4Elem4ProcRotateAndTest(mesh, meta);

    std::vector<stk::mesh::EntityProc> entities_to_move;
    if(p_rank == 0)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 1)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 2)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
        int dest_proc = 1;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 3)
    {
        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 4);
        int dest_proc = 2;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }

    mesh.change_entity_owner_exp(entities_to_move);

    std::vector<std::pair<int, int> > entities;

    int numNodes = 9;
    if(p_rank == 0)
    {
        entities.push_back(std::make_pair(2, 3));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 2));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(6, 1));

    }
    else if(p_rank == 1)
    {
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 2));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(6, 0));
        entities.push_back(std::make_pair(8, 2));
    }
    else if(p_rank == 2)
    {
        entities.push_back(std::make_pair(4, 3));
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 3));
        entities.push_back(std::make_pair(8, 1));
    }
    else
    {
        entities.push_back(std::make_pair(2, 0));
        entities.push_back(std::make_pair(4, 2));
        entities.push_back(std::make_pair(5, 0));
        entities.push_back(std::make_pair(5, 1));
        entities.push_back(std::make_pair(5, 2));
    }

    std::vector<stk::mesh::Entity> modifiedEntities;
    std::set<int> ids;
    for (size_t i=0;i<entities.size();i++)
    {
        ids.insert(entities[i].first);
    }

    for(int  i = 0; i < numNodes; i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
        if ( CEOUtils::isEntityInSharingCommMap(mesh, node ) && ids.find(i+1) == ids.end() )
        {
            modifiedEntities.push_back(node);
        }
    }

    for(int i = 0; i < numNodes; i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
    }

    for(size_t i = 0; i < entities.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, entities[i].first);

        CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, entities[i].second);
        modifiedEntities.push_back(node);
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_4Elem4ProcRotate(mesh, meta);
}

TEST(CEOME, change_entity_owner_3Elem4Proc1Edge3D)
{
    //  ID.proc
    //                    15.2--------16.2                      15.1--------16.1
    //                     /|          /|                        /|          /|
    //                    / |         / |                       / |         / |
    //                  7.2---------8.2 |                     7.1---------8.1 |
    //                   |  |  3.2   |  |                      |  |  3.1   |  |
    //                   |  |        |  |                      |  |        |  |
    //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
    //         /|        | *|        | /|   -->      /|        | *|        | /|
    //        / |        |* |        |/ |           / |        |* |        |/ |
    //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
    //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
    //       |  |        |  |        |  |          |  |        |  |        |  |
    //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
    //       | /         | /         | /           | /         | /         | /
    //       |/          |/          |/            |/          |/          |/
    //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
    //
    //      (*)edge: 1.0                          (*)edge: 1.1

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if(numProcs != 4)
    {
        return;
    }

    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    BulkDataTester mesh(meta, pm);
    const int p_rank = mesh.parallel_rank();
    CEOUtils::fillMeshfor3Elem4Proc1Edge3DAndTest(mesh, meta);

    std::vector<stk::mesh::EntityProc> entities_to_move;
    if(p_rank == 0)
    {
        Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
        int dest_proc = 3;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);

        elem = mesh.get_entity(stk::topology::EDGE_RANK, 1);
        dest_proc = 1;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    }
    else if(p_rank == 1)
    {
        Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
        int dest_proc = 0;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }
    else if(p_rank == 2)
    {
        Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
        int dest_proc = 1;
        entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
        CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
    }

    mesh.change_entity_owner_exp(entities_to_move);

//    CEOUtils::checkStatesAfterCEO_3Elem4Proc1Edge3D(mesh);

    ////////////////////////////////////////////////////////////////////////////

    std::vector<std::pair<int, int> > nodeEntities;
    std::vector<std::pair<int, int> > edgeEntities;

    int numNodes = 16;
    int numEdges = 1;

    if(p_rank == 0)
    {
        nodeEntities.push_back(std::make_pair(2, 3));
        nodeEntities.push_back(std::make_pair(10, 3));
        nodeEntities.push_back(std::make_pair(5, 1));
        nodeEntities.push_back(std::make_pair(5, 3));
        nodeEntities.push_back(std::make_pair(13, 1));
        nodeEntities.push_back(std::make_pair(13, 3));
        nodeEntities.push_back(std::make_pair(6, 1));
        nodeEntities.push_back(std::make_pair(14, 1));

        edgeEntities.push_back(std::make_pair(1, 1));
        edgeEntities.push_back(std::make_pair(1, 3));
    }
    else if(p_rank == 1)
    {
        nodeEntities.push_back(std::make_pair(5, 0));
        nodeEntities.push_back(std::make_pair(5, 3));
        nodeEntities.push_back(std::make_pair(13, 0));
        nodeEntities.push_back(std::make_pair(13, 3));
        nodeEntities.push_back(std::make_pair(6, 0));
        nodeEntities.push_back(std::make_pair(14, 0));

        edgeEntities.push_back(std::make_pair(1, 0));
        edgeEntities.push_back(std::make_pair(1, 3));
    }
    else if(p_rank == 2)
    {

    }
    else
    {
        nodeEntities.push_back(std::make_pair(2, 0));
        nodeEntities.push_back(std::make_pair(10, 0));
        nodeEntities.push_back(std::make_pair(5, 0));
        nodeEntities.push_back(std::make_pair(5, 1));
        nodeEntities.push_back(std::make_pair(13, 0));
        nodeEntities.push_back(std::make_pair(13, 1));

        edgeEntities.push_back(std::make_pair(1, 0));
        edgeEntities.push_back(std::make_pair(1, 1));
    }

    std::vector<stk::mesh::Entity> modifiedEntities;

    std::set<int> ids;
    for (size_t i=0;i<nodeEntities.size();i++)
    {
        ids.insert(nodeEntities[i].first);
    }

    for(int  i = 0; i < numNodes; i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, i+1);
        if ( CEOUtils::isEntityInSharingCommMap(mesh, node) && ids.find(i+1) == ids.end() )
        {
            modifiedEntities.push_back(node);
        }
    }

    for(int i = 0; i < numNodes; i++)
    {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, i + 1);
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
    }

    stk::mesh::EntityKey key(stk::topology::EDGE_RANK, 1);
    CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);

    for(size_t i = 0; i < nodeEntities.size(); i++)
    {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeEntities[i].first);

        CEOUtils::addSharingInfo(mesh, node, stk::mesh::BulkData::SHARED, nodeEntities[i].second);
        modifiedEntities.push_back(node);
    }

    ///////////////

    std::set<int> edgeIds;
    for (size_t i=0;i<edgeEntities.size();i++)
    {
        edgeIds.insert(edgeEntities[i].first);
    }

    for(int  i = 0; i < numEdges; i++)
    {
        stk::mesh::Entity edge = mesh.get_entity(stk::topology::EDGE_RANK, i+1);
        if ( CEOUtils::isEntityInSharingCommMap(mesh, edge) && edgeIds.find(i+1) == edgeIds.end() )
        {
            modifiedEntities.push_back(edge);
        }
    }

    for(size_t i = 0; i < edgeEntities.size(); i++)
    {
        stk::mesh::Entity edge = mesh.get_entity(stk::topology::EDGE_RANK, edgeEntities[i].first);

        CEOUtils::addSharingInfo(mesh, edge, stk::mesh::BulkData::SHARED, edgeEntities[i].second);
        modifiedEntities.push_back(edge);
    }

    update_mesh_based_on_changed_sharing_info(mesh, modifiedEntities);

    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    mesh.modification_begin();
    mesh.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_3Elem4Proc1Edge3D(mesh);
}

}

