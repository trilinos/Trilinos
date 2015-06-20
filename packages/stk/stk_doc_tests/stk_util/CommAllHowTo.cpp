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

#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll
#include <stk_util/parallel/MPI.hpp>
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include <stk_util/stk_config.h>
#include <limits>

#if defined ( STK_HAS_MPI )

//DocTest1
TEST(ParallelComm, HowToCommunicateOneValue)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);

    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    double sendSomeNumber = 100-myProcId;

    for(int phase = 0; phase < 2; ++phase)
    {
        for (int proc=0;proc<numProcs;proc++)
        {
            if ( proc != myProcId )
            {
                stk::CommBuffer& proc_buff = comm_all.send_buffer(proc);
                proc_buff.pack<double>(sendSomeNumber);
            }
        }
        if(phase == 0)
        {
            // NOTE: The value passed into allocate_buffers determines whether a sparse or dense communication will occur
            comm_all.allocate_buffers(numProcs / 4);
        }
        else
        {
            comm_all.communicate();
        }
    }


    for (int proc=0;proc<numProcs;proc++)
    {
        if ( proc != myProcId )
        {
            stk::CommBuffer& dataReceived = comm_all.recv_buffer(proc);
            double val = -1;
            dataReceived.unpack(val);
            EXPECT_EQ(100-proc, val);
        }
    }
}
//DocTest2
TEST(ParallelComm, HowToCommunicateAnArbitraryNumberOfValues)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);

    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    double sendSomeNumber = 100-myProcId;

    for(int phase = 0; phase < 2; ++phase)
    {
        for (int proc=0;proc<numProcs;proc++)
        {
            if ( proc != myProcId )
            {
                stk::CommBuffer& proc_buff = comm_all.send_buffer(proc);
                for (int i=0;i<myProcId;i++)
                {
                    proc_buff.pack<double>(sendSomeNumber+i);
                }
            }
        }
        if(phase == 0)
        {
            // NOTE: The value passed into allocate_buffers determines whether a sparse or dense communication will occur
            comm_all.allocate_buffers(numProcs / 4);
        }
        else
        {
            comm_all.communicate();
        }
    }


    for (int procFromWhichDataIsReceived=0;procFromWhichDataIsReceived<numProcs;procFromWhichDataIsReceived++)
    {
        if ( procFromWhichDataIsReceived != myProcId )
        {
            stk::CommBuffer& dataReceived = comm_all.recv_buffer(procFromWhichDataIsReceived);
            int numItemsReceived = 0;
            while ( dataReceived.remaining() )
            {
                double val = -1;
                dataReceived.unpack(val);
                EXPECT_EQ(100-procFromWhichDataIsReceived+numItemsReceived, val);
                numItemsReceived++;
            }
            int goldNumItemsReceived = procFromWhichDataIsReceived;
            EXPECT_EQ(goldNumItemsReceived, numItemsReceived);
        }
    }
}
//EndDocTest
#endif
