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
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>

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

void printMemoryStats(MPI_Comm comm)
{
    size_t maxHwm = 0, minHwm = 0, avgHwm = 0;
    stk::get_memory_high_water_mark_across_processors(comm, maxHwm, minHwm, avgHwm);

    int proc=-1;
    MPI_Comm_rank(comm, &proc);

    int numProcs=0;
    MPI_Comm_size(comm, &numProcs);

    if (proc == 0)
    {
      std::ostringstream os;
      const double bytesInMegabyte = 1024*1024;
      os << std::setw(6) << std::fixed << std::setprecision(1) << "Max HWM: "<<double(maxHwm)/double(bytesInMegabyte)
        <<", Min HWM: "<<double(minHwm)/double(bytesInMegabyte)<<", Avg HWM: "<<avgHwm/bytesInMegabyte<<std::endl;

      std::cerr << os.str();
    }
}

void printPeformanceStats(double elapsedTime, MPI_Comm comm)
{
    size_t maxHwm = 0, minHwm = 0, avgHwm = 0;
    stk::get_memory_high_water_mark_across_processors(comm, maxHwm, minHwm, avgHwm);

    int proc=-1;
    MPI_Comm_rank(comm, &proc);

    int numProcs=0;
    MPI_Comm_size(comm, &numProcs);

    double minTime = 0, maxTime = 0, avgTime = 0;
    MPI_Allreduce(&elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&elapsedTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, comm);
    double elapsedTimeDivided = elapsedTime/numProcs;
    MPI_Allreduce(&elapsedTimeDivided, &avgTime, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (proc == 0)
    {
      std::ostringstream os;
      const double bytesInMegabyte = 1024*1024;
      os << "Max time: "  << maxTime << ", Min time: " << minTime << ", Avg time: " << avgTime << std::endl;
      os << std::setw(6) << std::fixed << std::setprecision(1) << "Max HWM: "<<double(maxHwm)/double(bytesInMegabyte)
        <<", Min HWM: "<<double(minHwm)/double(bytesInMegabyte)<<", Avg HWM: "<<avgHwm/bytesInMegabyte<<std::endl;
      std::cerr << os.str();
    }
}

namespace
{
//==============================================================================

void update_comm_maps_based_on_modified_entities_in_sharing(BulkDataTester& bulk, std::vector<stk::mesh::Entity>& modifiedEntities)
{
    bulk.my_resolve_ownership_of_modified_entities(modifiedEntities);
    bulk.my_move_entities_to_proper_part_ownership(modifiedEntities);
    bulk.my_update_comm_list_based_on_changes_in_comm_map();
    bulk.my_update_comm_list(modifiedEntities);
}

void erase_all_sharing_for_valid_entities_on_comm_map(BulkDataTester &mesh)
{
    for(size_t i=0; i<mesh.comm_list().size(); ++i)
    {
        stk::mesh::EntityKey key = mesh.comm_list()[i].key;
        stk::mesh::Entity entity = mesh.get_entity(key);
        if(!mesh.is_valid(entity))
        {
            CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
        }
    }
}

void fill_entities_that_have_lost_sharing_info(BulkDataTester &mesh, const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
        const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating)
{
    std::set<stk::mesh::EntityKey> keysThatNeedToHaveCorrectSharingInfo;
    for (size_t i=0;i<sharedEntities.size();i++)
    {
        keysThatNeedToHaveCorrectSharingInfo.insert(sharedEntities[i].first);
    }

    for (size_t i=0;i<entitiesThatUsedToHaveSharingInfoBeforeCEO.size();i++)
    {
        stk::mesh::EntityKey entityKey = mesh.entity_key(entitiesThatUsedToHaveSharingInfoBeforeCEO[i]);
        const bool keyNotInSharedKeysList = keysThatNeedToHaveCorrectSharingInfo.find(entityKey) == keysThatNeedToHaveCorrectSharingInfo.end();
        if ( keyNotInSharedKeysList )
        {
            modifiedEntitiesForWhichCommMapsNeedUpdating.push_back(entitiesThatUsedToHaveSharingInfoBeforeCEO[i]);
            CEOUtils::eraseSharingInfoUsingKey(mesh, entityKey, stk::mesh::BulkData::SHARED);
        }
    }
}

void fillModifiedEntitiesListBasedOnSharingInfo(BulkDataTester &mesh, const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
        const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating)
{
    erase_all_sharing_for_valid_entities_on_comm_map(mesh);
    fill_entities_that_have_lost_sharing_info(mesh, sharedEntities, entitiesThatUsedToHaveSharingInfoBeforeCEO, modifiedEntitiesForWhichCommMapsNeedUpdating);

    for (size_t i=0;i<sharedEntities.size();i++)
    {
        stk::mesh::EntityKey key = sharedEntities[i].first;
        CEOUtils::eraseSharingInfoUsingKey(mesh, key, stk::mesh::BulkData::SHARED);
    }

    for(size_t i = 0; i < sharedEntities.size(); i++)
    {
        stk::mesh::Entity entity = mesh.get_entity(sharedEntities[i].first);
        CEOUtils::addSharingInfo(mesh, entity, stk::mesh::BulkData::SHARED, sharedEntities[i].second);
        modifiedEntitiesForWhichCommMapsNeedUpdating.push_back(entity);
    }
}

void getEntitiesThatHaveSharing(BulkDataTester &bulk, std::vector<stk::mesh::Entity> &entitiesThatHaveSharingInfo,
        stk::mesh::NodeToDependentProcessorsMap &entityKeySharing)
{
    int myProcId = bulk.parallel_rank();
    stk::CommAll commStage1(bulk.parallel());

    size_t numEntitiesThatHaveSharingInfo = 0;

    for(int phase = 0; phase < 2; ++phase)
    {
        if(phase == 1)
        {
            entitiesThatHaveSharingInfo.resize(numEntitiesThatHaveSharingInfo);
            numEntitiesThatHaveSharingInfo = 0;
        }

        for(stk::mesh::EntityRank irank = stk::topology::NODE_RANK; irank <= stk::topology::FACE_RANK; ++irank)
        {
            for(stk::mesh::impl::EntityRepository::const_iterator i = bulk.begin_entities(irank); i != bulk.end_entities(irank); ++i)
            {
                stk::mesh::Entity entity = i->second;

                //TODD: entities after mesh creation are in created state. Could "skip" this state? Or different algorithm?
                if(bulk.state(entity) != stk::mesh::Unchanged)
                {
                    if(phase == 0 && CEOUtils::isEntityInSharingCommMap(bulk, entity))
                    {
                        numEntitiesThatHaveSharingInfo++;
                    }
                    else if(phase == 1 && CEOUtils::isEntityInSharingCommMap(bulk, entity))
                    {
                        entitiesThatHaveSharingInfo[numEntitiesThatHaveSharingInfo] = entity;
                        numEntitiesThatHaveSharingInfo++;
                    }

                    int procThatOwnsEntity = bulk.parallel_owner_rank(entity);
                    bool anotherProcOwnsThisEntity = procThatOwnsEntity != myProcId;
                    if(anotherProcOwnsThisEntity)
                    {
                        //TODD: should not need to send myProcId, because buffer is based on proc id
                        stk::mesh::EntityKey entityKey = bulk.entity_key(entity);
                        commStage1.send_buffer(procThatOwnsEntity).pack < stk::mesh::EntityKey > (entityKey);
                    }
                }
            }
        }

        if(phase == 0)
        {
            commStage1.allocate_buffers(bulk.parallel_size() / 4, 0);
        }
        else
        {
            commStage1.communicate();
        }
    }

    for(int procIndex = 0; procIndex < bulk.parallel_size(); procIndex++)
    {
        if(myProcId == procIndex) continue;
        stk::CommBuffer & dataFromAnotherProc = commStage1.recv_buffer(procIndex);
        EntityKey key;
        int sharingProc = procIndex;
        while(dataFromAnotherProc.remaining())
        {
            dataFromAnotherProc.unpack<stk::mesh::EntityKey>(key);
            Entity entity = bulk.get_entity(key);
            ThrowRequireMsg(bulk.is_valid(entity) && bulk.parallel_owner_rank(entity) == myProcId, "Entitykey " << key << " is not owned by receiving processor " << myProcId);
            entityKeySharing[key].insert(sharingProc);
        }
    }
}

//TODD: NodeToDependentProcessMap should be renamed to something like EntityToDepen...

void getLocallyModifiedSharedEntities(BulkDataTester &bulk, stk::mesh::NodeToDependentProcessorsMap &entityKeySharing, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities)
{
    int myProcId = bulk.parallel_rank();

    stk::mesh::NodeToDependentProcessorsMap::iterator iter = entityKeySharing.begin();
    for(; iter != entityKeySharing.end(); iter++)
    {
        std::vector<int> sharingProcs(iter->second.begin(), iter->second.end());
        iter->second.insert(myProcId);
        for(size_t j = 0; j < sharingProcs.size(); j++)
        {
            sharedEntities.push_back(std::pair<stk::mesh::EntityKey, int>(iter->first, sharingProcs[j]));
        }
    }

    stk::CommAll commStage2(bulk.parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        iter = entityKeySharing.begin();
        for(; iter != entityKeySharing.end(); iter++)
        {
            std::vector<int> sharingProcs(iter->second.begin(), iter->second.end());
            for(size_t j = 0; j < sharingProcs.size(); j++)
            {
                if(sharingProcs[j] == myProcId)continue;
                commStage2.send_buffer(sharingProcs[j]).pack<stk::mesh::EntityKey>(iter->first);
                commStage2.send_buffer(sharingProcs[j]).pack<size_t>(sharingProcs.size());
                for(size_t k = 0; k < sharingProcs.size(); k++)
                {
                    commStage2.send_buffer(sharingProcs[j]).pack<int>(sharingProcs[k]);
                }
            }
        }

        if(phase == 0)
        { // allocation phase
            commStage2.allocate_buffers(bulk.parallel_size() / 4, 0);
        }
        else
        { // communication phase
            commStage2.communicate();
        }
    }

    for(int procIndex = 0; procIndex < bulk.parallel_size(); procIndex++)
    {
        if(myProcId == procIndex) continue;
        stk::CommBuffer & dataFromAnotherProc = commStage2.recv_buffer(procIndex);
        while(dataFromAnotherProc.remaining())
        {
            EntityKey key;
            size_t numSharingProcs = 0;
            dataFromAnotherProc.unpack<stk::mesh::EntityKey>(key);
            dataFromAnotherProc.unpack<size_t>(numSharingProcs);
            for(size_t j = 0; j < numSharingProcs; j++)
            {
                int sharingProc = -1;
                dataFromAnotherProc.unpack<int>(sharingProc);
                if(sharingProc
             != myProcId)
                {
                    sharedEntities.push_back(std::pair<stk::mesh::EntityKey, int>(key, sharingProc));
                }
            }
        }
    }
}

void makeThingsConsistent(BulkDataTester &bulk)
{
    printMemoryStats(bulk.parallel());
    double startTime = stk::wall_time();

    std::vector<stk::mesh::Entity> entitiesThatHaveSharingInfo;
    stk::mesh::NodeToDependentProcessorsMap ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns;

    getEntitiesThatHaveSharing(bulk, entitiesThatHaveSharingInfo, ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns);

    std::vector<std::pair<stk::mesh::EntityKey, int> > sharedEntities;
    getLocallyModifiedSharedEntities(bulk, ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns, sharedEntities);

    std::vector<stk::mesh::Entity> modifiedEntities;
    fillModifiedEntitiesListBasedOnSharingInfo(bulk, sharedEntities, entitiesThatHaveSharingInfo, modifiedEntities);

    update_comm_maps_based_on_modified_entities_in_sharing(bulk, modifiedEntities);

    bulk.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, bulk.parallel());

    bulk.modification_begin();
    bulk.internal_modification_end_for_change_entity_owner_exp(true, stk::mesh::BulkData::MOD_END_SORT);
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

    makeThingsConsistent(bulk);

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

    CEOUtils::checkStatesAfterCEO_2Elem2ProcFlip(mesh);

    ////////////////////////////////////////////////////////////////////////////

    makeThingsConsistent(mesh);

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

    makeThingsConsistent(mesh);

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

    makeThingsConsistent(mesh);

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
    //        |      |     |    ||     |          |      |     |    ||     |
    //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
    //        |      |     |    ||     |          |      |     |    ||     |
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

    if (p_rank == 1 || p_rank == 2) {
        // Fill elem_field data on nodes 5 and 6 with data.
        Entity node5 = mesh.get_entity(stk::topology::NODE_RANK,5);
        Entity node6 = mesh.get_entity(stk::topology::NODE_RANK,6);
        stk::mesh::FieldBase* elem_field = meta_data.get_field(stk::topology::NODE_RANK, "elem_field");
        double * elem_field_data_node5 = static_cast<double*>(stk::mesh::field_data(*elem_field,node5));
        *elem_field_data_node5 = 5.0;
        double * elem_field_data_node6 = static_cast<double*>(stk::mesh::field_data(*elem_field,node6));
        *elem_field_data_node6 = 6.0;
    }

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

    CEOUtils::checkStatesAfterCEO_4Elem4ProcEdge(mesh);

    if (p_rank == 1 || p_rank == 0) {
        Entity node5 = mesh.get_entity(stk::topology::NODE_RANK,5);
        Entity node6 = mesh.get_entity(stk::topology::NODE_RANK,6);
        stk::mesh::FieldBase* elem_field = meta_data.get_field(stk::topology::NODE_RANK, "elem_field");
        double * elem_field_data_node5 = static_cast<double*>(stk::mesh::field_data(*elem_field,node5));
        EXPECT_EQ( 5.0, *elem_field_data_node5 );
        double * elem_field_data_node6 = static_cast<double*>(stk::mesh::field_data(*elem_field,node6));
        EXPECT_EQ( 6.0, *elem_field_data_node6 );
    }

    ////////////////////////////////////////////////////////////////////////////

    makeThingsConsistent(mesh);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_4Elem4ProcEdge(mesh);

    if (p_rank == 1 || p_rank == 0) {
        Entity node5 = mesh.get_entity(stk::topology::NODE_RANK,5);
        Entity node6 = mesh.get_entity(stk::topology::NODE_RANK,6);
        stk::mesh::FieldBase* elem_field = meta_data.get_field(stk::topology::NODE_RANK, "elem_field");
        double * elem_field_data_node5 = static_cast<double*>(stk::mesh::field_data(*elem_field,node5));
        EXPECT_EQ( 5.0, *elem_field_data_node5 );
        double * elem_field_data_node6 = static_cast<double*>(stk::mesh::field_data(*elem_field,node6));
        EXPECT_EQ( 6.0, *elem_field_data_node6 );
    }

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

    CEOUtils::checkStatesAfterCEO_8Elem4ProcMoveTop(mesh);

    ////////////////////////////////////////////////////////////////////////////

    makeThingsConsistent(mesh);

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

    CEOUtils::checkStatesAfterCEO_4Elem4ProcRotate(mesh, meta);

    ////////////////////////////////////////////////////////////////////////////

    makeThingsConsistent(mesh);

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

    CEOUtils::checkStatesAfterCEO_3Elem4Proc1Edge3D(mesh);

    ////////////////////////////////////////////////////////////////////////////

    makeThingsConsistent(mesh);

    ////////////////////////////////////////////////////////////////////////////

    CEOUtils::checkStatesAfterCEOME_3Elem4Proc1Edge3D(mesh);
}

}

