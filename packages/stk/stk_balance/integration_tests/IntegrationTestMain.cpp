/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for TestInfo, etc
#include <mpi.h>                        // for MPI_Comm_rank, MPI_Finalize, etc
#include <stdarg.h>                     // for va_end, va_list, va_start
#include <stdio.h>                      // for printf, vprintf, fflush, etc
#include "gtest/gtest-test-part.h"      // for TestPartResult
#include "stk_unit_test_utils/ParallelGtestOutput.hpp"

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
    int procId = 0;

#if defined(STK_BUILT_IN_SIERRA) // this means MPI is available
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
#endif

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    stk::unit_test_util::create_parallel_output(procId);

    int returnVal = RUN_ALL_TESTS();

#if defined(STK_BUILT_IN_SIERRA)
    MPI_Finalize();
#endif

    return returnVal;
}
