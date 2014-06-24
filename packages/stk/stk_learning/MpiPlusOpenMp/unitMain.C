#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

int main(int argc, char **argv)
{
    stk::parallel_machine_init(&argc, &argv);

    testing::InitGoogleTest(&argc, argv);
    int errorCode = RUN_ALL_TESTS();
    
    stk::parallel_machine_finalize();
    return errorCode;
}
