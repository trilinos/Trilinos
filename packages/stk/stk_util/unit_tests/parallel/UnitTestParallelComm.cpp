/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include <stk_util/stk_config.h>
#if defined ( STK_HAS_MPI )
#  include <mpi.h>                        // for MPI_Comm


TEST(ParallelComm, CommAllDestructor)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD ;
    int mpi_size = stk::parallel_machine_size(comm);
    const unsigned zero = 0 ;
    std::vector<unsigned> msg_size( mpi_size , zero );
    const unsigned * const s_size = & msg_size[0] ;
    {
        stk::CommAll sparse;
        // This will allocate both m_send and m_recv regardless of mpi_size
        sparse.allocate_symmetric_buffers( comm, s_size );
        stk::CommBuffer * recv_buffer = &sparse.recv_buffer(0);
        stk::CommBuffer * send_buffer = &sparse.send_buffer(0);
        ASSERT_TRUE( recv_buffer != send_buffer );
    }
    // This should not produce a memory leak.
}

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
