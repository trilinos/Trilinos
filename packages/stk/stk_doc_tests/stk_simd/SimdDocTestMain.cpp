#include <gtest/gtest.h>                // for InitGoogleTest, etc
#include <Kokkos_Core.hpp>

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
    Kokkos::initialize(argc, argv);

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    int returnVal = RUN_ALL_TESTS();

    Kokkos::finalize_all();

    return returnVal;
}
