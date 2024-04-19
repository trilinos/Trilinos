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

#include "gtest/gtest.h"
#include "stk_util/environment/WallTime.hpp"                // for wall_time
#include "stk_util/parallel/DistributedIndex.hpp"           // for DistributedIndex::KeyTypeVector
#include "stk_util/parallel/GenerateParallelUniqueIDs.hpp"  // for generate_parallel_unique_ids
#include "stk_util/parallel/MPI.hpp"                        // for Datatype
#include "stk_util/parallel/Parallel.hpp"                   // for MPI_Bcast, MPI_Comm, MPI_Reduce
#include "stk_util/stk_config.h"                            // for STK_HAS_MPI
#include "stk_util/util/ReportHandler.hpp"                  // for ThrowRequireMsg
#include <cstddef>                                          // for size_t
#include <algorithm>                                        // for binary_search, sort, copy
#include <cstdint>                                          // for uint64_t
#include <iostream>                                         // for operator<<, basic_ostream, endl
#include <memory>                                           // for allocator_traits<>::value_type
#include <string>                                           // for string
#include <vector>                                           // for vector, vector<>::iterator

#if defined(STK_HAS_MPI)  // means that MPI is available

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

INTMPI whichProcOwnsId(const uint64_t maxId, const uint64_t id, INTMPI numProcs);
void distributeIds(std::vector<uint64_t> &myIds, const MpiInfo& mpiInfo);
uint64_t getNumIdsPerProc(const uint64_t maxId, const INTMPI numProcs);
void getIdUsageAcrossAllProcs(std::vector< std::vector<uint64_t> > &idsToComm, std::vector<uint64_t> &idsInUseAcrossAllProcsInMyRange, const MpiInfo& mpiInfo);

void retrieveIds(const INTMPI root, uint64_t id, MPI_Comm comm, uint64_t numIdsToGetPerProc, std::vector<int>& areIdsBeingUsed);
void respondToRootProcessorAboutIdsOwnedOnThisProc(const int root, const uint64_t maxId, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);

bool sendIdToCheck(const INTMPI root, uint64_t id, MPI_Comm comm);
void receiveIdAndCheck(const int root, const std::vector<uint64_t> &idsInUse, MPI_Comm comm);
void checkUniqueIds(const std::vector<uint64_t> &myIds, const std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo);
void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t>& myIds, const std::vector<uint64_t> &uniqueIds);
void generate_ids(const uint64_t maxId, const std::vector<uint64_t> &idsInUse, std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo);
void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        std::vector<uint64_t> &sortedIds, const uint64_t maxId, const MpiInfo& mpiInfo);
void terminateIdRequestForThisProc(INTMPI root, MPI_Comm comm);
void getAvailableIds(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo);
void getAvailableIds_exp(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo);

////////////////////////////////////////////////////////////////////


TEST(GeneratedIds, whichProcOwnsId)
{
    uint64_t maxId=10;
    INTMPI numProcs=3;

    uint64_t ids[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    INTMPI procOwner[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3};
    for (int i=0;i<10;i++)
    {
        EXPECT_EQ(procOwner[i], whichProcOwnsId(maxId, ids[i], numProcs));
    }
}

TEST(GeneratedIds, distributeIds)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    int numIdsPerProc = 10;
    std::vector<uint64_t> myIds(numIdsPerProc);
    distributeIds(myIds, mpiInfo);

    uint64_t offsetForProc = numIdsPerProc*mpiInfo.getProcId();
    for (size_t i=0;i<myIds.size();i++)
    {
        uint64_t oneBasedId = i+1;
        EXPECT_EQ(oneBasedId+offsetForProc, myIds[i]);
    }
}

TEST(GeneratedIds, getNumIdsPerProc)
{
    INTMPI numProcs = 10;
    uint64_t findIdsLessThan = 1000;

    uint64_t numIdsOneOneProc = getNumIdsPerProc(findIdsLessThan, numProcs);
    EXPECT_EQ(findIdsLessThan/numProcs, numIdsOneOneProc);
}

TEST(GeneratedIds, getIdUsageAcrossAllProcs)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() != 1 )
    {
        std::vector< std::vector<uint64_t> > idsToComm(mpiInfo.getNumProcs());
        INTMPI procToSendId = 0;
        uint64_t madeUpId = static_cast<uint64_t>(mpiInfo.getProcId());
        idsToComm[procToSendId].push_back(madeUpId);
        std::vector<uint64_t> idsInUseAcrossAllProcsInMyRange;
        getIdUsageAcrossAllProcs(idsToComm, idsInUseAcrossAllProcsInMyRange, mpiInfo);
        if ( mpiInfo.getProcId() == 0 )
        {
            for (uint64_t id=0;id<idsInUseAcrossAllProcsInMyRange.size();id++)
            {
                EXPECT_EQ(id, idsInUseAcrossAllProcsInMyRange[id]);
            }
        }
        else
        {
            EXPECT_EQ(0u, idsInUseAcrossAllProcsInMyRange.size());
        }
    }
}

TEST(GeneratedIds, findUniqueIdAcrossProcs)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    uint64_t totalIdsInUse = 362880;
    uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();

    std::vector<uint64_t> myIds(numIdsThisProc,0);
    distributeIds(myIds, mpiInfo);

    uint64_t startingIdToSearchForNewIds = 1;
    uint64_t numIdsNeeded = 100;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 10000000;

    double startTime = stk::wall_time();

    getAvailableIds(myIds, numIdsNeeded, idsObtained, startingIdToSearchForNewIds, maxId, mpiInfo);

    double endTime = stk::wall_time();

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(myIds, idsObtained, mpiInfo);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

TEST(GeneratedIds, findUniqueIdAcrossProcsApproach1)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    uint64_t numIdsThisProc = 50000;
    uint64_t totalIdsInUse = numIdsThisProc*mpiInfo.getNumProcs();

    std::vector<uint64_t> myIds(numIdsThisProc,0);
    distributeIds(myIds, mpiInfo);

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
    writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

TEST(GeneratedIds, findUniqueIdAcrossProcsVaryingNumIdsInUse)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    uint64_t numIdsThisProc = 50000;
    uint64_t totalIdsInUse = numIdsThisProc*mpiInfo.getNumProcs();

    std::vector<uint64_t> myIds(numIdsThisProc,0);
    distributeIds(myIds, mpiInfo);

    uint64_t startingIdToSearchForNewIds = 1;
    uint64_t numIdsNeeded = 1000;
    std::vector<uint64_t> idsObtained;
    uint64_t maxId = 100000000000;

    double startTime = stk::wall_time();

    getAvailableIds(myIds, numIdsNeeded, idsObtained, startingIdToSearchForNewIds, maxId, mpiInfo);

    double endTime = stk::wall_time();

    if ( mpiInfo.getProcId() == 0 )
    {
        std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
    }

    checkUniqueIds(myIds, idsObtained, mpiInfo);
    writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, idsObtained);

    for (size_t i=0;i<idsObtained.size();i++)
    {
        uint64_t goldId = totalIdsInUse + numIdsNeeded*mpiInfo.getProcId()+1 + i;
        EXPECT_EQ(goldId, idsObtained[i]);
    }
}

TEST(GeneratedIds, numIdsInUseVariesWithNumberProcessorsNumIdsNeededVariesWithNumberProcessors)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() < 10 )
    {
        uint64_t maxId = 10000000;
        uint64_t numIdsThisProc = 100000;

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t numIdsNeeded = 10000;
        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

// 2^56 = 72, 057, 594, 037, 927, 936 --> 72 quadrillion, 57 trillion, 594 billion, 37 million, 927 thousand, 936 ids!

TEST(GeneratedIds, numIdsInUseConstantNumIdsNeededConstant)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() <= 10 )
    {
        uint64_t totalIdsInUse = 362880;
        uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t maxId = 10000000;

        uint64_t numIdsNeeded = 10000;
        if ( mpiInfo.getProcId() != 0 )
        {
            numIdsNeeded = 0;
        }

        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

TEST(GeneratedIds, numIdsInUseConstantNumIdsNeededVariesWithNumberProcessors)
{
    MpiInfo mpiInfo(MPI_COMM_WORLD);

    if ( mpiInfo.getNumProcs() <= 10 )
    {
        uint64_t totalIdsInUse = 362880;
        uint64_t numIdsThisProc = totalIdsInUse/mpiInfo.getNumProcs();

        std::vector<uint64_t> myIds(numIdsThisProc,0);
        distributeIds(myIds, mpiInfo);

        uint64_t maxId = 10000000;

        uint64_t numIdsNeeded = 10000*mpiInfo.getNumProcs();
        if ( mpiInfo.getProcId() != 0 )
        {
            numIdsNeeded = 0;
        }

        std::vector<uint64_t> uniqueIds(numIdsNeeded,0);

        MPI_Barrier(mpiInfo.getMpiComm());

        double startTime = stk::wall_time();

        generate_ids(maxId, myIds, uniqueIds, mpiInfo);

        double endTime = stk::wall_time();

        if ( mpiInfo.getProcId() == 0 )
        {
            std::cerr << "Took " << endTime-startTime << " seconds." << std::endl;
        }

        EXPECT_EQ(numIdsNeeded, uniqueIds.size());
        checkUniqueIds(myIds, uniqueIds, mpiInfo);

        writeIdsToFile("ids_", mpiInfo.getProcId(), myIds, uniqueIds);
    }
}

TEST(GeneratedIds, distributed_index_vs_generate_parallel_unique_ids)
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  PDIndex::KeySpanVector partition_spans(1) ;

  partition_spans[0].first  = 1;
  partition_spans[0].second = 10000;

  PDIndex di( comm , partition_spans );

  size_t numIdsNeeded = 1;
  std::vector<size_t> requests( partition_spans.size() , numIdsNeeded );

  std::vector< PDIndex::KeyTypeVector > generated_ids_from_di ;
  di.generate_new_keys( requests , generated_ids_from_di );

  uint64_t maxAllowedId = 10000;
  std::vector<uint64_t> existingIds;
  std::vector<uint64_t> ids_from_gi = stk::generate_parallel_unique_ids(maxAllowedId, existingIds, numIdsNeeded, comm);

  ASSERT_EQ(ids_from_gi.size(), generated_ids_from_di[0].size());
}

////////////////////////////////////////////////////////////////////

void distributeIds(std::vector<uint64_t> &myIds, const MpiInfo &mpiInfo)
{
    uint64_t offset = myIds.size();

    for (uint64_t i=0;i<myIds.size();i++)
    {
        myIds[i] = offset*mpiInfo.getProcId() + i + 1;
    }
}

////////////////////////////////////////////////////////////////////

void getIdUsageAcrossAllProcs(std::vector< std::vector<uint64_t> > &idsToComm, std::vector<uint64_t> &idsInUseAcrossAllProcsInMyRange, const MpiInfo &mpiInfo)
{
    for (INTMPI i=0;i<mpiInfo.getNumProcs();i++)
    {
        if ( i == mpiInfo.getProcId() )
        {
            for (INTMPI j=0;j<mpiInfo.getNumProcs();j++)
            {
                if ( j != mpiInfo.getProcId() )
                {
                    uint64_t numItemsToComm = idsToComm[j].size();
                    MPI_Send(&numItemsToComm, 1, sierra::MPI::Datatype<uint64_t>::type(), j, mpiInfo.getNumProcs()*i+j, mpiInfo.getMpiComm());
                    if ( numItemsToComm > 0 )
                    {
                      MPI_Send(idsToComm[j].data(), numItemsToComm, sierra::MPI::Datatype<uint64_t>::type(), j,
                          mpiInfo.getNumProcs() * i + j, mpiInfo.getMpiComm());
                    }
                }
            }
        }
        else
        {
            uint64_t numItemsToReceive=0;
            MPI_Status status1;
            MPI_Recv(&numItemsToReceive, 1, sierra::MPI::Datatype<uint64_t>::type(), i, mpiInfo.getNumProcs()*i+mpiInfo.getProcId(), mpiInfo.getMpiComm(), &status1);
            if ( numItemsToReceive > 0 )
            {
                std::vector<uint64_t> idsFromOtherProc(numItemsToReceive,0);
                MPI_Request request;
                MPI_Irecv(idsFromOtherProc.data(), numItemsToReceive, sierra::MPI::Datatype<uint64_t>::type(), i,
                    mpiInfo.getNumProcs() * i + mpiInfo.getProcId(), mpiInfo.getMpiComm(), &request);
                MPI_Status status2;
                MPI_Wait(&request, &status2);
                idsInUseAcrossAllProcsInMyRange.insert(idsInUseAcrossAllProcsInMyRange.end(), idsFromOtherProc.begin(), idsFromOtherProc.end());
            }
        }
    }

    idsInUseAcrossAllProcsInMyRange.insert(idsInUseAcrossAllProcsInMyRange.end(), idsToComm[mpiInfo.getProcId()].begin(), idsToComm[mpiInfo.getProcId()].end());
    std::sort(idsInUseAcrossAllProcsInMyRange.begin(), idsInUseAcrossAllProcsInMyRange.end());
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
                STK_ThrowRequireMsg(uniqueIds[j]>0, "Id generation error. Please contact sierra-help for support.");
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

////////////////////////////////////////////////////////////////////

void writeIdsToFile(const std::string &filename, const INTMPI myProcId, const std::vector<uint64_t>& myIds, const std::vector<uint64_t> &uniqueIds)
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

uint64_t getNumIdsPerProc(const uint64_t maxId, const INTMPI numProcs)
{
    return maxId/numProcs;
}

////////////////////////////////////////////////////////////////////

INTMPI whichProcOwnsId(const uint64_t maxId, const uint64_t id, INTMPI numProcs)
{
    STK_ThrowRequireMsg(id>0, "Invalid Id. Contact sierra-help for support.");
    uint64_t numIdsPerProc = getNumIdsPerProc(maxId, numProcs);
    INTMPI procOwner = (id-1)/numIdsPerProc;
    return procOwner;
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
        MPI_Reduce(areIdsBeingused.data(), rbuff, numIdsToGet, MPI_INT, MPI_SUM, root, comm);
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
        MPI_Reduce(zeroids.data(), areIdsBeingUsed.data(), numIdsToGetPerProc, MPI_INT, MPI_SUM, root, comm);
    }
}


TEST(GeneratedIds, multiPerf) {
#ifdef COMMENT_OUT
  //
  //  Performance and scalability comparision of a number of id generation routines.
  //
  //  16384 ids per processor initially.
  //  Add 1000 ids in each of ten steps
  //
  //  Case 1:  Ids initially densly packed
  //  Case 2:  Ids initially spread
  //

  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  MpiInfo mpiInfo(MPI_COMM_WORLD);

  unsigned numInit = 32768*mpi_size;
  unsigned numNew  = 1000;


  {
    srand((unsigned)0);


    std::vector<unsigned> inUse32a;
    std::vector<unsigned> inUse32c;
    std::vector<uint64_t> inUse64a;
    std::vector<unsigned> new32a;
    std::vector<unsigned> new32c;
    std::vector<uint64_t> new64a;
    std::vector<unsigned> inUse32b;
    std::vector<unsigned> order32c;
    std::vector<uint64_t> inUse64b;
    std::vector<unsigned> new32b;
    std::vector<uint64_t> new64b;

    unsigned maxAllowableId32 = ~0U;
    uint64_t maxAllowableId64 = maxAllowableId32;

    for(unsigned i=1; i<numInit+1; ++i) {
      unsigned targetRank = rand()%mpi_size;
      if(targetRank == (unsigned)mpi_rank) {
        inUse32a.push_back(i);
        inUse64a.push_back(i);
      }
    }

    unsigned numToFill = numInit;
    unsigned firstIndex = 1;

    unsigned maxIndex = 0;

    while(numToFill > 0) {
      unsigned curNumToFill = numToFill/2;
      if(curNumToFill == 0) curNumToFill = 1;

      for(unsigned i=0; i<curNumToFill; ++i) {

        unsigned targetRank = rand()%mpi_size;

        maxIndex = std::max(maxIndex, firstIndex+1);

        if(targetRank == (unsigned)mpi_rank) {
          inUse32b.push_back(firstIndex + i);
          inUse32c.push_back(firstIndex + i);
          inUse64b.push_back(firstIndex + i);
        }
      }
      numToFill -= curNumToFill;
      firstIndex = firstIndex + (maxAllowableId32 - firstIndex)/1.8;
    }
    order32c.resize(numNew);
    for(unsigned i=0; i<numNew; ++i) {
      order32c[i] = mpi_rank*numNew + i;
    }


    double timeGenerateParallelUniqueA = 0.0;
    double timeGenerateA               = 0.0;
    double timeGenerateParallelUniqueB = 0.0;
    double timeGenerateParallelUniqueC = 0.0;
    double timeGenerateB               = 0.0;
    double timeDI                      = 0.0;
    double timeRED                         = 0.0;

    {
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();


    for(int iter=0; iter<512; ++iter) {
      unsigned globalSubDivCount[16];
      unsigned localSubDivCount[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
      int mpiResult = MPI_Scan(localSubDivCount, globalSubDivCount, 16, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    }
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeRED += (endTime-startTime);

    }

    for(int iter=0; iter<64; ++iter) {
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();
      new32a = stk::generate_parallel_unique_ids(maxAllowableId32, inUse32a, numNew, MPI_COMM_WORLD, inUse32a.back());
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeGenerateParallelUniqueA += (endTime-startTime);
      for(unsigned inew=0; inew<new32a.size(); ++inew) {
        inUse32a.push_back(new32a[inew]);
      }
    }

    for(int iter=0; iter<64; ++iter) {
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();
      new32b = stk::generate_parallel_unique_ids(maxAllowableId32, inUse32b, numNew, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeGenerateParallelUniqueB += (endTime-startTime);
      for(unsigned inew=0; inew<new32b.size(); ++inew) {
        inUse32b.push_back(new32b[inew]);
      }
    }

    for(int iter=0; iter<64; ++iter) {
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();
      new32c = stk::generate_parallel_consistent_ids(maxAllowableId32, inUse32c, order32c, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeGenerateParallelUniqueC += (endTime-startTime);
      for(unsigned inew=0; inew<new32c.size(); ++inew) {
        inUse32c.push_back(new32c[inew]);
        //order32c.push_back(mpi_size*numNew + iter*mpi_size*);
      }
    }

    if(mpi_rank == 0) {
      std::cout<<"CASE 1: TIME 32 A: "<<timeGenerateParallelUniqueA<<std::endl;
      std::cout<<"CASE 1: TIME 64 A: "<<timeGenerateA<<std::endl;
      std::cout<<"CASE 1: TIME 32 B: "<<timeGenerateParallelUniqueB<<std::endl;
      std::cout<<"CASE 1: TIME 32 C: "<<timeGenerateParallelUniqueC<<std::endl;
      std::cout<<"CASE 1: TIME 64 B: "<<timeGenerateB<<std::endl;
      std::cout<<"CASE 1: TIME DI: "<<timeDI<<std::endl;
      std::cout<<"CASE 1: TIME RED: "<<timeRED<<std::endl;
    }


    for(int iter=0; iter<64; ++iter) {
      new64a.resize(numNew);
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();
      generate_ids(maxAllowableId64, inUse64a, new64a, mpiInfo);
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeGenerateA += (endTime-startTime);
      for(unsigned inew=0; inew<new64a.size(); ++inew) {
        inUse64a.push_back(new64a[inew]);
      }
    }
    for(int iter=0; iter<64; ++iter) {
      new64b.resize(numNew);
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();
      generate_ids(maxAllowableId64, inUse64b, new64b, mpiInfo);
      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeGenerateB += (endTime-startTime);
      for(unsigned inew=0; inew<new64b.size(); ++inew) {
        inUse64b.push_back(new64b[inew]);
      }
    }


    if(mpi_rank == 0) {
      std::cout<<"CASE 1: TIME 32 A: "<<timeGenerateParallelUniqueA<<std::endl;
      std::cout<<"CASE 1: TIME 64 A: "<<timeGenerateA<<std::endl;
      std::cout<<"CASE 1: TIME 32 B: "<<timeGenerateParallelUniqueB<<std::endl;
      std::cout<<"CASE 1: TIME 64 B: "<<timeGenerateB<<std::endl;
      std::cout<<"CASE 1: TIME DI: "<<timeDI<<std::endl;
    }


    typedef stk::parallel::DistributedIndex PDIndex;
    PDIndex::KeySpanVector partition_spans;
    enum { test_spans_count = 100 };
    enum { test_spans_size  = 1000000 };

    partition_spans.resize( test_spans_count );

    for ( unsigned i = 0 ; i < test_spans_count ; ++i ) {
      partition_spans[i].first  = 1 + test_spans_size * i * 2 ;
      partition_spans[i].second = test_spans_size * ( i * 2 + 1 );
    }
    PDIndex di( MPI_COMM_WORLD , partition_spans );
    std::vector<size_t> requests( partition_spans.size() , 0u );
    std::vector< PDIndex::KeyTypeVector > generated_keys ;
    PDIndex::KeyProcVector sharing_of_local_keys ;

    PDIndex::KeyTypeVector keys_to_add ;
    PDIndex::KeyTypeVector keys_to_remove ;

    //------------------------------
    // Add 32768 odd keys per process
    // starting at the beginning of the partition.

    const size_t old_size_multiplier = 16 ;

    for ( size_t j = 0 ; j < partition_spans.size() ; ++j ) {
      PDIndex::KeyType key_first = partition_spans[j].first ;
      if ( 0 == key_first % 2 ) { ++key_first ; } // Make it odd
      key_first += old_size_multiplier * mpi_rank ;

      const size_t n = old_size_multiplier * 1024 ;
      for ( size_t i = 0 ; i < n ; ++i ) {
        PDIndex::KeyType key = key_first + 2 * i ;
        keys_to_add.push_back( key );
      }
    }

    di.update_keys( keys_to_add.begin(), keys_to_add.end() , keys_to_remove.begin(), keys_to_remove.end() );

    for(int iter=0; iter<64; ++iter) {
      const size_t gen_count = 10 ;
      for ( size_t i = 0 ; i < requests.size() ; ++i ) {
        if ( i % 2 ) {
          requests[i] = gen_count ;
        }
        else {
          requests[i] = 0 ;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      double startTime = stk::wall_time();

      di.generate_new_keys( requests , generated_keys );

      MPI_Barrier(MPI_COMM_WORLD);
      double endTime = stk::wall_time();
      timeDI += (endTime-startTime);
    }


    if(mpi_rank == 0) {
      std::cout<<"CASE 1: TIME 32 A: "<<timeGenerateParallelUniqueA<<std::endl;
      std::cout<<"CASE 1: TIME 64 A: "<<timeGenerateA<<std::endl;
      std::cout<<"CASE 1: TIME 32 B: "<<timeGenerateParallelUniqueB<<std::endl;
      std::cout<<"CASE 1: TIME 64 B: "<<timeGenerateB<<std::endl;
      std::cout<<"CASE 1: TIME DI: "<<timeDI<<std::endl;
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////

void generate_ids(const uint64_t maxId, const std::vector<uint64_t> &idsInUse, std::vector<uint64_t> &uniqueIds, const MpiInfo& mpiInfo)
{
    std::vector< std::vector<uint64_t> > idsToComm(mpiInfo.getNumProcs());

    for (size_t i=0;i<idsInUse.size();i++)
    {
        INTMPI procOwner = whichProcOwnsId(maxId, idsInUse[i], mpiInfo.getNumProcs());
        STK_ThrowRequireMsg(static_cast<int>(procOwner)<mpiInfo.getNumProcs(), "Id generation error. Please contact sierra-help. procOwner = " << procOwner << ", maxId = " << maxId << ", and id = " << idsInUse[i] << std::endl);
        idsToComm[procOwner].push_back(idsInUse[i]);
    }

    std::vector<uint64_t> idsInUseAcrossAllProcsInMyRange;
    getIdUsageAcrossAllProcs(idsToComm, idsInUseAcrossAllProcsInMyRange, mpiInfo);

    std::vector<uint64_t> uniqueIdsFound;

    uint64_t myIndexStart = getNumIdsPerProc(maxId, mpiInfo.getNumProcs())*mpiInfo.getProcId()+1;
    uint64_t myIndexEnd = getNumIdsPerProc(maxId, mpiInfo.getNumProcs())*(mpiInfo.getProcId()+1)+1;

    size_t numIdsNeeded = uniqueIds.size();
    size_t numIdsAvailableThisProc =  getNumIdsPerProc(maxId, mpiInfo.getNumProcs()) - idsInUseAcrossAllProcsInMyRange.size();
    STK_ThrowRequireMsg(numIdsNeeded <= numIdsAvailableThisProc, "Not enough unique ids available (Id generation error). Plrease contact sierra-help for support. Number of ids needed: "
            << numIdsNeeded << " and num available ids: " << numIdsAvailableThisProc);

    if ( !uniqueIds.empty() )
    {
        uniqueIdsFound.reserve(numIdsNeeded);
        for (uint64_t i=myIndexStart;i<myIndexEnd;i++)
        {
            if ( !std::binary_search(idsInUseAcrossAllProcsInMyRange.begin(), idsInUseAcrossAllProcsInMyRange.end(), i) )
            {
                uniqueIdsFound.push_back(i);
                if ( uniqueIdsFound.size() == uniqueIds.size() )
                {
                    break;
                }
            }
        }
    }

    STK_ThrowRequireMsg(uniqueIdsFound.size() == uniqueIds.size(), "Id generation error. Could not obtain needed ids. Please contact sierra-help for support.");
    std::copy(uniqueIdsFound.begin(), uniqueIdsFound.end(), uniqueIds.begin());
}

////////////////////////////////////////////////////////////////////

void getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(INTMPI root, uint64_t &startingIdToSearchForNewIds, std::vector<uint64_t> &idsObtained, uint64_t numIdsNeeded, int scaleFactorForNumIds,
        std::vector<uint64_t> &sortedIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    while ( startingIdToSearchForNewIds < maxId && idsObtained.size() < numIdsNeeded )
    {
        int requestNumIds = numIdsNeeded - idsObtained.size();
        std::vector<int> areIdsBeingUsed(scaleFactorForNumIds*requestNumIds,0);
        if ( !std::binary_search(sortedIds.begin(), sortedIds.end(), startingIdToSearchForNewIds) )
        {
            retrieveIds(root, startingIdToSearchForNewIds, mpiInfo.getMpiComm(), scaleFactorForNumIds*requestNumIds, areIdsBeingUsed);
            int numIdsChecked=0;
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

////////////////////////////////////////////////////////////////////

void getAvailableIds(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    std::vector<uint64_t> receivedInfo(mpiInfo.getNumProcs(),0);
    MPI_Allgather(&numIdsNeeded, 1, sierra::MPI::Datatype<uint64_t>::type(), receivedInfo.data(), 1,
        sierra::MPI::Datatype<uint64_t>::type(), mpiInfo.getMpiComm());

    std::vector<uint64_t> sortedIds(myIds.begin(), myIds.end());
    std::sort(sortedIds.begin(), sortedIds.end());

    int scaleFactorForNumIds = 10;

    for (INTMPI procIndex=0;procIndex<mpiInfo.getNumProcs();procIndex++)
    {
        if ( receivedInfo[procIndex] != 0 )
        {
            if ( procIndex == mpiInfo.getProcId() )
            {
                getBatchesOfIdsFromOtherProcessorsUntilRequestOnThisProcIsFulfilled(procIndex, startingIdToSearchForNewIds, idsObtained, numIdsNeeded, scaleFactorForNumIds, sortedIds, maxId, mpiInfo);
                STK_ThrowRequireMsg(idsObtained.size()==numIdsNeeded, "Id generation error. Ran out of ids. Please contact sierra-help for support.");
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

void getAvailableIds_exp(const std::vector<uint64_t> &myIds, uint64_t numIdsNeeded, std::vector<uint64_t> &idsObtained, uint64_t &startingIdToSearchForNewIds, const uint64_t maxId, const MpiInfo& mpiInfo)
{
    INTMPI numprocs = mpiInfo.getNumProcs();
    std::vector<uint64_t> receivedInfo(numprocs,0);
    MPI_Allgather(&numIdsNeeded, 1, sierra::MPI::Datatype<uint64_t>::type(), receivedInfo.data(), 1,
        sierra::MPI::Datatype<uint64_t>::type(), mpiInfo.getMpiComm());

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
                    STK_ThrowRequireMsg(idsObtained.size()==numIdsNeeded, "Id generation error. Ran out of ids. Please contact sierra-help for support.");
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

}

#endif // STK_HAS_MPI

