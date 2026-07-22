#include "gtest/gtest.h"
#include "stk_unit_test_utils/algorithmTimer.hpp"  // for RunInfo, print_run_info, time_algorithm
#include "stk_unit_test_utils/getOption.h"         // for get_command_line_option
#include "stk_util/parallel/Parallel.hpp"          // for parallel_machine_size, MPI_COMM_WORLD
#include <cstddef>                                 // for size_t
#include <iostream>                                // for cerr

TEST(RunTimer, runs)
{
    const size_t numWork = stk::unit_test_util::get_command_line_option<size_t>("-work", 10000);
    const double tolerance = stk::unit_test_util::get_command_line_option<double>("-tol", 1e-6);
    const size_t minRuns = stk::unit_test_util::get_command_line_option<double>("-min", 10000);
    ASSERT_TRUE(tolerance > 0);

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    unitTestUtils::RunInfo runInfo = unitTestUtils::time_algorithm(tolerance, minRuns, MPI_COMM_WORLD,
        [numWork]()
        {
            double sum = 0;
            for(size_t i=0; i<numWork; i++)
            {
                sum++;
            }
            EXPECT_EQ(static_cast<double>(numWork), sum);
        }
    );
    EXPECT_TRUE(runInfo.mean > 0.0);
    EXPECT_TRUE(runInfo.numRuns >= minRuns);

    unitTestUtils::print_run_info(std::cerr,
                                  "test",
                                  stk::parallel_machine_size(comm),
                                  runInfo);
}
