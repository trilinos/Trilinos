// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include <mpi.h>          // for MPI_Finalize, MPI_Init
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

int main(int argc, char **argv)
{
    MPI_Init( &argc , &argv );
    testing::InitGoogleTest(&argc, argv);
    int returnVal = RUN_ALL_TESTS();
    MPI_Finalize();

    return returnVal;
}


