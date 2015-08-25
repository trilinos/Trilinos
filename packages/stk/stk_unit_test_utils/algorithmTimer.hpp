#ifndef UNITTESTUTILS_ALGORITHM_TIMER
#define UNITTESTUTILS_ALGORITHM_TIMER

#include <stk_util/environment/CPUTime.hpp>
#include <cmath>
#include <stddef.h>

namespace unitTestUtils
{

struct RunInfo
{
    RunInfo(double t, size_t num) : mean(t), numRuns(num) {}
    double mean;
    size_t numRuns;
};

template <typename ALGORITHM>
void initial_run_to_allocate_memory_and_stuff(const ALGORITHM &alg)
{
    alg();
}

template <typename ALGORITHM>
RunInfo time_algorithm(const double tolerance, const ALGORITHM &alg)
{
    initial_run_to_allocate_memory_and_stuff(alg);

    double mean = 0.0;
    double lastMean = 0.0;
    double duration = 0.0;
    size_t numRuns = 0;
    bool done = false;
    while(!done)
    {
        double startTime = stk::cpu_time();
        alg();
        duration += stk::cpu_time() - startTime;

        mean = duration / numRuns;
        double meanDiff = std::abs(mean - lastMean);
        if(numRuns > 1000 && duration > 0.0 && meanDiff < tolerance)
        {
            done = true;
        }
        lastMean = mean;
        numRuns++;
    }
    return RunInfo(mean, numRuns);
}

}

#endif
