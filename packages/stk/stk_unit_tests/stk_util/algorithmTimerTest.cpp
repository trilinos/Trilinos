#include <stk_unit_test_utils/algorithmTimer.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <gtest/gtest.h>                // for TEST, AssertHelper, etc


TEST(RunTimer, runs)
{
    const size_t numWork = unitTestUtils::get_command_line_option<size_t>("-work", "10000");
    const double tolerance = unitTestUtils::get_command_line_option<double>("-tol", "1e-6");
    const size_t minRuns = unitTestUtils::get_command_line_option<double>("-min", "10000");
    ASSERT_TRUE(tolerance > 0);

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
    EXPECT_TRUE(runInfo.numRuns > 0);

    std::cerr << "num runs = " << runInfo.numRuns << std::endl;
    std::cerr << "mean = " << runInfo.mean << std::endl;
}
