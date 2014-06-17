#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#if defined ( STK_HAS_MPI )
#  include <mpi.h>                        // for MPI_Comm
#endif

int main(int argc, char **argv)
{
#if defined ( STK_HAS_MPI )
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int errorCode = RUN_ALL_TESTS();
    
#if defined ( STK_HAS_MPI )
    MPI_Finalize();
#endif
    return errorCode;
}
