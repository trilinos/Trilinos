/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(ParallelComm, CommAllDestructor) {
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
        STKUNIT_ASSERT( recv_buffer != send_buffer );
    }
    // This should not produce a memory leak.
}
