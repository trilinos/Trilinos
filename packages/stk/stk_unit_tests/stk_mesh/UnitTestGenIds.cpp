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
//       from this software getNumIdsPerProcwithout specific prior written permission.
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
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/MPI.hpp>
#include <fstream>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>
#include <stk_util/stk_config.h>
#include <stk_unit_test_utils/getOption.h>
#include <stdint.h>

//====================

#if defined(STK_HAS_MPI)  // this means MPI is available

namespace {

typedef int INTMPI;

class MpiInfo
{
public:
    MpiInfo(MPI_Comm comm) : mProcId(-1), mNumProcs(-1), mComm(comm)
    {
        MPI_Comm_rank(mComm, &mProcId);
        MPI_Comm_size(mComm, &mNumProcs);
    }
    ~MpiInfo() {}

    INTMPI getProcId() const { return mProcId; }
    INTMPI getNumProcs() const { return mNumProcs; }
    MPI_Comm getMpiComm() const { return mComm; }

private:
    INTMPI mProcId;
    INTMPI mNumProcs;
    MPI_Comm mComm;
};

void retrieveIds(const INTMPI root, uint64_t id, MPI_Comm comm, uint64_t numIdsToGetPerProc, std::vector<int>& areIdsBeingUsed);
bool sendIdToCheck(const INTMPI root, uint64_t id, MPI_Comm comm);

// not stkmesh
void respondToRootProcessorAboutIdsOwnedOnThisProc(const int root, const uint64_t maxId, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);
void receiveIdAndCheck(const int root, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);
void checkUniqueIds(const std::vector<uint64_t> &myIds, const std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo);
void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        std::vector<uint64_t> &sortedIds, const uint64_t maxId, const MpiInfo& mpiInfo);

// is stkmesh
void respondToRootProcessorAboutIdsOwnedOnThisProc(const int root, const uint64_t maxId, stk::mesh::BulkData& stkMeshBulkData, MPI_Comm comm);
void receiveIdAndCheck(const int root, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm comm);
void checkUniqueIds(stk::mesh::BulkData &stkMeshBulkData, const std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo);
void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        stk::mesh::BulkData& stkMeshBulkData, const uint64_t maxId, const MpiInfo& mpiInfo);

void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t> &uniqueIds);
void terminateIdRequestForThisProc(INTMPI root, MPI_Comm comm);
void getAvailableIds_exp(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo);
void getAvailableIds_exp(stk::mesh::BulkData &stkMeshBulkData, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo);

////////////////////////////////////////////////////////////////////

TEST(GeneratedIds, StkMeshApproach1)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    std::string exodusFileName = unitTestUtils::getOption("-i", "generated:10x10x10");
    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, mpiInfo.getMpiComm());

    // STK IO module will be described in separate chapter.
    // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
    // The order of the following lines in {} are important
    {
      stk::io::StkMeshIoBroker exodusFileReader(mpiInfo.getMpiComm());

      // Inform STK IO which STK Mesh objects to populate later
      exodusFileReader.set_bulk_data(stkMeshBulkData);

      exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

      // Populate the MetaData which has the descriptions of the Parts and Fields.
      exodusFileReader.create_input_mesh();

      // Populate entities in STK Mesh from Exodus file
      exodusFileReader.populate_bulk_data();
    }

    stk::mesh::Selector select_owned( stkMeshMetaData.locally_owned_part() );

    std::vector<unsigned> count;
    stk::mesh::count_entities( select_owned , stkMeshBulkData , count );
    size_t numIdsThisProc = count[stk::topology::NODE_RANK];

    std::vector<size_t> count1;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, count1);
    size_t totalIdsInUse = count1[stk::topology::NODE_RANK];

    std::vector<uint64_t> myIds(numIdsThisProc,0);

    const stk::mesh::BucketVector& buckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);

    size_t counter = 0;
    for (size_t i=0;i<buckets.size();i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if ( bucket.owned() )
        {
            for (size_t j=0;j<bucket.size();j++)
            {
                stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(bucket[j]);
                myIds[counter] = entityKey;
                counter++;
            }
        }
    }

    EXPECT_EQ(counter, numIdsThisProc);

    uint64_t startingIdToSearchForNewIds = 1;
    uint64_t numIdsNeeded = 100;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 72057594037927935;

    double startTime = stk::wall_time();

    getAvailableIds_exp(myIds, numIdsNeeded, idsObtained, startingIdToSearchForNewIds, maxId, mpiInfo);

    double endTime = stk::wall_time();

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(myIds, idsObtained, mpiInfo);
    writeIdsToFile("ids_", mpiInfo.getProcId(), idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

////////////////////////////////////////////////////////////////////

TEST(GeneratedIds, StkMeshApproach2)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    std::string exodusFileName = unitTestUtils::getOption("-i", "generated:10x10x10");
    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, mpiInfo.getMpiComm());

    // STK IO module will be described in separate chapter.
    // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
    // The order of the following lines in {} are important
    {
      stk::io::StkMeshIoBroker exodusFileReader(mpiInfo.getMpiComm());

      // Inform STK IO which STK Mesh objects to populate later
      exodusFileReader.set_bulk_data(stkMeshBulkData);

      exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

      // Populate the MetaData which has the descriptions of the Parts and Fields.
      exodusFileReader.create_input_mesh();

      // Populate entities in STK Mesh from Exodus file
      exodusFileReader.populate_bulk_data();
    }

    ////////////////////////////////////////////////////

    std::vector<size_t> count1;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, count1);
    size_t totalIdsInUse = count1[stk::topology::NODE_RANK];

    uint64_t startingIdToSearchForNewIds = 1;
    uint64_t numIdsNeeded = 100;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 72057594037927935;

    double startTime = stk::wall_time();

    getAvailableIds_exp(stkMeshBulkData, numIdsNeeded, idsObtained, startingIdToSearchForNewIds, maxId, mpiInfo);

    double endTime = stk::wall_time();

    ////////////////////////////////////////////////////

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(stkMeshBulkData, idsObtained, mpiInfo);
    writeIdsToFile("ids_", mpiInfo.getProcId(), idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

////////////////////////////////////////////////////////////////////

void checkUniqueIds(const std::vector<uint64_t> &myIds, const std::vector<uint64_t> &uniqueIds, const MpiInfo &mpiInfo)
{
    std::vector<uint64_t> sortedIds(uniqueIds.begin(), uniqueIds.end());
    std::sort(sortedIds.begin(), sortedIds.end());

    for (INTMPI i=0;i<mpiInfo.getNumProcs();i++)
    {
        if (mpiInfo.getNumProcs() == i)
        {
            for (size_t j=0;j<uniqueIds.size();j++)
            {
                ThrowRequireMsg(uniqueIds[j]>0, "Id generation error. Please contact sierra-help for support.");
                EXPECT_FALSE(std::binary_search(myIds.begin(), myIds.end(), uniqueIds[j]));
                EXPECT_EQ(true, sendIdToCheck(i, uniqueIds[j], mpiInfo.getMpiComm()));
            }
            terminateIdRequestForThisProc(i, mpiInfo.getMpiComm());
        }
        else
        {
            receiveIdAndCheck(i, sortedIds, mpiInfo.getMpiComm());
        }
    }
}

void checkUniqueIds(stk::mesh::BulkData &stkMeshBulkData, const std::vector<uint64_t> &uniqueIds, const MpiInfo &mpiInfo)
{
    for (INTMPI i=0;i<mpiInfo.getNumProcs();i++)
    {
        if (mpiInfo.getNumProcs() == i)
        {
            for (size_t j=0;j<uniqueIds.size();j++)
            {
                ThrowRequireMsg(uniqueIds[j]>0, "Id generation error. Please contact sierra-help for support.");

                stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, uniqueIds[j]);
                EXPECT_FALSE(stkMeshBulkData.is_valid(entity));
                EXPECT_EQ(true, sendIdToCheck(i, uniqueIds[j], mpiInfo.getMpiComm()));
            }
            terminateIdRequestForThisProc(i, mpiInfo.getMpiComm());
        }
        else
        {
            receiveIdAndCheck(i, stkMeshBulkData, mpiInfo.getMpiComm());
        }
    }
}

////////////////////////////////////////////////////////////////////

void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t> &uniqueIds)
{
#ifdef DEBUG_THIS_
    std::ostringstream os;
    os << "ids_" << myProcId << ".m";
    std::ofstream out(os.str().c_str());
//    for (size_t i=0;i<myIds.size();i++)
//    {
//        out << "myids(" << i+1 << ") = " << myIds[i] << ";" << std::endl;
//    }
    for (size_t i=0;i<uniqueIds.size();i++)
    {
        out << "unique_ids(" << i+1 << ") = " << uniqueIds[i] << ";" << std::endl;
    }
    out.close();
#endif
}

////////////////////////////////////////////////////////////////////

bool sendIdToCheck(const INTMPI root, uint64_t id, MPI_Comm comm)
{
    MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
    bool goodId = true;
    if ( id != 0 )
    {
        uint64_t good = 0;
        uint64_t received = 0;
        MPI_Reduce(&good, &received, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, root, comm);

        if ( received > 0 )
        {
            goodId = false;
        }
    }
    return goodId;
}

////////////////////////////////////////////////////////////////////

void receiveIdAndCheck(const int root, const std::vector<uint64_t> &idsInUse, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
        if ( id == 0) break;

        bool found = std::binary_search(idsInUse.begin(), idsInUse.end(), id);
        uint64_t result = 0;
        if ( found )
        {
            result = 1;
        }
        MPI_Reduce(&result, &result, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////
void receiveIdAndCheck(const int root, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
        if ( id == 0) break;

        stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, id);
        bool found = stkMeshBulkData.is_valid(entity);
        uint64_t result = 0;
        if ( found )
        {
            result = 1;
        }
        MPI_Reduce(&result, &result, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////

void respondToRootProcessorAboutIdsOwnedOnThisProc(const int root, const uint64_t maxId, const std::vector<uint64_t> &idsInUse, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
        if ( id == 0) break;
        uint64_t numIdsToGet=0;
        MPI_Bcast(&numIdsToGet, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);

        std::vector<int> areIdsBeingused(numIdsToGet,0);
        for (size_t i=0;i<areIdsBeingused.size();i++)
        {
            if ( std::binary_search(idsInUse.begin(), idsInUse.end(), id+i) )
            {
                areIdsBeingused[i] = 1;
            }
        }
        int *rbuff = 0;
        MPI_Reduce(&areIdsBeingused[0], rbuff, numIdsToGet, MPI_INT, MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////

void respondToRootProcessorAboutIdsOwnedOnThisProc(const int root, const uint64_t maxId, stk::mesh::BulkData& stkMeshBulkData, MPI_Comm comm)
{
    uint64_t id=0;

    while( true )
    {
        MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
        if ( id == 0) break;
        uint64_t numIdsToGet=0;
        MPI_Bcast(&numIdsToGet, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);

        std::vector<int> areIdsBeingused(numIdsToGet,0);
        for (size_t i=0;i<areIdsBeingused.size();i++)
        {
            stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, id+i);
            if ( stkMeshBulkData.is_valid(entity) )
            {
                areIdsBeingused[i] = 1;
            }
        }
        int *rbuff = 0;
        MPI_Reduce(&areIdsBeingused[0], rbuff, numIdsToGet, MPI_INT, MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////

void retrieveIds(const INTMPI root, uint64_t id, MPI_Comm comm, uint64_t numIdsToGetPerProc, std::vector<int>& areIdsBeingUsed)
{
    MPI_Bcast(&id, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
    if ( id != 0 )
    {
        MPI_Bcast(&numIdsToGetPerProc, 1, sierra::MPI::Datatype<uint64_t>::type(), root, comm);
        std::vector<uint64_t> zeroids(numIdsToGetPerProc,0);
        MPI_Reduce(&zeroids[0], &areIdsBeingUsed[0], numIdsToGetPerProc, MPI_INT, MPI_SUM, root, comm);
    }
}

////////////////////////////////////////////////////////////////////


void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        std::vector<uint64_t> &sortedIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    while ( startingIdToSearchForNewIds < maxId && idsObtained.size() < numIdsNeeded )
    {
        uint64_t requestNumIds = numIdsNeeded - idsObtained.size();
        std::vector<int> areIdsBeingUsed(scaleFactorForNumIds*requestNumIds,0);
        if ( !std::binary_search(sortedIds.begin(), sortedIds.end(), startingIdToSearchForNewIds) )
        {
            retrieveIds(root, startingIdToSearchForNewIds, mpiInfo.getMpiComm(), scaleFactorForNumIds*requestNumIds, areIdsBeingUsed);
            uint64_t numIdsChecked=0;
            for (size_t i=0;i<areIdsBeingUsed.size();i++)
            {
                numIdsChecked=i;
                if ( areIdsBeingUsed[i] == 0 )
                {
                    idsObtained.push_back(startingIdToSearchForNewIds+i);
                    if ( idsObtained.size() == numIdsNeeded ) break;
                }
            }
            startingIdToSearchForNewIds += numIdsChecked+1;
        }
        else
        {
            startingIdToSearchForNewIds++;
        }
    }
}

////////////////////////////////////////////////////////////////////


void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        stk::mesh::BulkData &stkMeshBulkData, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    while ( startingIdToSearchForNewIds < maxId && idsObtained.size() < numIdsNeeded )
    {
        uint64_t requestNumIds = numIdsNeeded - idsObtained.size();
        std::vector<int> areIdsBeingUsed(scaleFactorForNumIds*requestNumIds,0);
        stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, startingIdToSearchForNewIds);
        if ( !stkMeshBulkData.is_valid(entity) )
        {
            retrieveIds(root, startingIdToSearchForNewIds, mpiInfo.getMpiComm(), scaleFactorForNumIds*requestNumIds, areIdsBeingUsed);
            uint64_t numIdsChecked=0;
            for (size_t i=0;i<areIdsBeingUsed.size();i++)
            {
                numIdsChecked=i;
                if ( areIdsBeingUsed[i] == 0 )
                {
                    idsObtained.push_back(startingIdToSearchForNewIds+i);
                    if ( idsObtained.size() == numIdsNeeded ) break;
                }
            }
            startingIdToSearchForNewIds += numIdsChecked+1;
        }
        else
        {
            startingIdToSearchForNewIds++;
        }
    }
}

////////////////////////////////////////////////////////////////////

void terminateIdRequestForThisProc(INTMPI root, MPI_Comm comm)
{
    sendIdToCheck(root, 0, comm); // a zero terminates communication
}

void getAvailableIds_exp(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    INTMPI numprocs = mpiInfo.getNumProcs();
    std::vector<uint64_t> receivedInfo(numprocs,0);
    MPI_Allgather(&numIdsNeeded, 1, sierra::MPI::Datatype<uint64_t>::type(), &receivedInfo[0], 1, sierra::MPI::Datatype<uint64_t>::type(), mpiInfo.getMpiComm());

    std::vector<uint64_t> sortedIds(myIds.begin(), myIds.end());
    std::sort(sortedIds.begin(), sortedIds.end());

    uint64_t largestIdHere = sortedIds.back();
    uint64_t largestIdEverywhere = 0;

    MPI_Allreduce(&largestIdHere, &largestIdEverywhere, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_MAX, mpiInfo.getMpiComm());

    uint64_t totalNumberOfIdsNeeded = 0;
    uint64_t offsetId=0;
    for (INTMPI i=0;i<numprocs;i++)
    {
        totalNumberOfIdsNeeded += receivedInfo[i];
        if ( i < mpiInfo.getProcId() )
        {
            offsetId += receivedInfo[i];
        }
    }

    if ( maxId - totalNumberOfIdsNeeded > largestIdHere )
    {
        startingIdToSearchForNewIds = largestIdHere;
        idsObtained.resize(numIdsNeeded);
        for (size_t i=0;i<idsObtained.size();i++)
        {
            idsObtained[i] = largestIdEverywhere + offsetId + i + 1;
        }
    }
    else
    {
        int scaleFactorForNumIds = 1;

        for (INTMPI procIndex=0;procIndex<numprocs;procIndex++)
        {
            if ( receivedInfo[procIndex] != 0 )
            {
                if ( procIndex == mpiInfo.getProcId() )
                {
                    getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(procIndex, startingIdToSearchForNewIds, idsObtained, numIdsNeeded, scaleFactorForNumIds, sortedIds, maxId, mpiInfo);
                    ThrowRequireMsg(idsObtained.size()==numIdsNeeded, "Id generation error. Ran out of ids. Please contact sierra-help for support.");
                    terminateIdRequestForThisProc(procIndex, mpiInfo.getMpiComm());
                }
                else
                {
                    respondToRootProcessorAboutIdsOwnedOnThisProc(procIndex, maxId, sortedIds, mpiInfo.getMpiComm());
                }
                // updated starting id across all procs
                MPI_Bcast(&startingIdToSearchForNewIds, 1, sierra::MPI::Datatype<uint64_t>::type(), procIndex, mpiInfo.getMpiComm());
            }
        }
    }
}

void getAvailableIds_exp(stk::mesh::BulkData &stkMeshBulkData, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    INTMPI numprocs = mpiInfo.getNumProcs();
    std::vector<uint64_t> receivedInfo(numprocs,0);
    MPI_Allgather(&numIdsNeeded, 1, sierra::MPI::Datatype<uint64_t>::type(), &receivedInfo[0], 1, sierra::MPI::Datatype<uint64_t>::type(), mpiInfo.getMpiComm());

    stk::mesh::EntityId largestIdHere = 0;

    const stk::mesh::BucketVector& buckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);

    for (size_t i=0;i<buckets.size();i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if ( bucket.owned() )
        {
            for (size_t j=0;j<bucket.size();j++)
            {
                stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(bucket[j]);
                if ( entityKey.id() > largestIdHere )
                {
                    largestIdHere = entityKey.id();
                }
            }
        }
    }

    uint64_t largestIdEverywhere = 0;
    MPI_Allreduce(&largestIdHere, &largestIdEverywhere, 1, sierra::MPI::Datatype<uint64_t>::type(), MPI_MAX, mpiInfo.getMpiComm());

    uint64_t totalNumberOfIdsNeeded = 0;
    uint64_t offsetId=0;
    for (INTMPI i=0;i<numprocs;i++)
    {
        totalNumberOfIdsNeeded += receivedInfo[i];
        if ( i < mpiInfo.getProcId() )
        {
            offsetId += receivedInfo[i];
        }
    }

    if ( maxId - totalNumberOfIdsNeeded > largestIdHere )
    {
        startingIdToSearchForNewIds = largestIdHere;
        idsObtained.resize(numIdsNeeded);
        for (size_t i=0;i<idsObtained.size();i++)
        {
            idsObtained[i] = largestIdEverywhere + offsetId + i + 1;
        }
    }
    else
    {
        int scaleFactorForNumIds = 1;

        for (INTMPI procIndex=0;procIndex<numprocs;procIndex++)
        {
            if ( receivedInfo[procIndex] != 0 )
            {
                if ( procIndex == mpiInfo.getProcId() )
                {
                    getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(procIndex, startingIdToSearchForNewIds, idsObtained, numIdsNeeded, scaleFactorForNumIds, stkMeshBulkData, maxId, mpiInfo);
                    ThrowRequireMsg(idsObtained.size()==numIdsNeeded, "Id generation error. Ran out of ids. Please contact sierra-help for support.");
                    terminateIdRequestForThisProc(procIndex, mpiInfo.getMpiComm());
                }
                else
                {
                    respondToRootProcessorAboutIdsOwnedOnThisProc(procIndex, maxId, stkMeshBulkData, mpiInfo.getMpiComm());
                }
                // updated starting id across all procs
                MPI_Bcast(&startingIdToSearchForNewIds, 1, sierra::MPI::Datatype<uint64_t>::type(), procIndex, mpiInfo.getMpiComm());
            }
        }
    }
}

}

#endif // STK_HAS_MPI

