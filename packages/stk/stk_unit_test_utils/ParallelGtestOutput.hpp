#ifndef STK_PARALLEL_GTEST_OUTPUT_H
#define STK_PARALLEL_GTEST_OUTPUT_H

#include <mpi.h>

namespace stk
{
namespace unit_test_util
{

void create_parallel_output(int procId);
void create_parallel_output_with_comm(int procId, MPI_Comm comm);

}
}


#endif
