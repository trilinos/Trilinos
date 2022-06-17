#include <gtest/gtest.h>
#include <Compadre_KokkosParser.hpp>
#include "unittests/test_XYZ.hpp"
#include "unittests/test_NeighborLists.hpp"
#include "unittests/test_PointCloudSearch.hpp"
#include "unittests/test_LinearAlgebra.hpp"
#include "unittests/test_Targets.hpp"
#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

#define ASSERT_NO_DEATH(statement) \
ASSERT_EXIT({{ statement } ::exit(EXIT_SUCCESS); }, ::testing::ExitedWithCode(0), "")

//
// KokkosParser tests go here because they modify the Kokkos backend
//
TEST (KokkosInitialize, NoArgsGiven) { 
    Kokkos::InitArguments args;
    ASSERT_NO_DEATH({
            // default constructor is hidden for KokkosParser
            // but still visible from this test
            auto kp = Compadre::KokkosParser(false);
            kp.finalize();
    });
}
TEST (KokkosInitialize, NoCommandLineArgsGiven) { 
    ASSERT_NO_DEATH({
            std::vector<std::string> arguments = {"--kokkos-threads=4"};
            auto kp = KokkosParser(arguments);
            kp.finalize();
    });
}


// this provides main(),
// but all other tests come from ./unittests/*.cpp
int main(int argc, char **argv) {

    // initializes MPI (if available) with command line arguments given
    #ifdef COMPADRE_USE_MPI
    MPI_Init(&argc, &argv);
    #endif

    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "Kokkos*";
    int sig = RUN_ALL_TESTS();

    // initializes kokkos
    auto kp = KokkosParser(argc, argv, true);

    // execute all tests
    ::testing::GTEST_FLAG(filter) = "-Kokkos*";
    sig += RUN_ALL_TESTS();

    // finalize Kokkos and MPI (if available)
    kp.finalize();

    // finialize MPI (if available)
    #ifdef COMPADRE_USE_MPI
    MPI_Finalize();
    #endif
    return sig; 
}
