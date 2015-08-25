#ifndef UNITTESTUTILS_ALGORITHM_TIMER
#define UNITTESTUTILS_ALGORITHM_TIMER

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <iostream>
#include <cmath>
#include <stddef.h>

namespace unitTestUtils
{

struct RunInfo
{
    RunInfo() : mean(0.0), max(0.0), min(1e20), numRuns(0) {}
    RunInfo(double t, double x, double n, size_t num) : mean(t), max(x), min(n), numRuns(num) {}
    double mean;
    double max;
    double min;
    size_t numRuns;
};

template <typename ALGORITHM>
void initial_run_to_allocate_memory_and_stuff(const ALGORITHM &alg)
{
    alg();
}

template <typename ALGORITHM>
RunInfo time_algorithm(const double tolerance, const size_t minRuns, stk::ParallelMachine comm, const ALGORITHM &alg)
{
    initial_run_to_allocate_memory_and_stuff(alg);

    RunInfo runInfo;

    double lastMean = 0.0;
    double duration = 0.0;
    int done = false;
    while(!done)
    {
        double startTime = stk::cpu_time();
        alg();
        double timeForOneRun = stk::cpu_time() - startTime;
        duration += timeForOneRun;
        runInfo.numRuns++;

        if(timeForOneRun > runInfo.max)
        {
            runInfo.max = timeForOneRun;
        }
        if(timeForOneRun < runInfo.min)
        {
            runInfo.min = timeForOneRun;
        }

        runInfo.mean = duration / runInfo.numRuns;
        double meanDiff = std::abs(runInfo.mean - lastMean);
        lastMean = runInfo.mean;
        if(runInfo.numRuns >= minRuns && duration > 0.0 && meanDiff < tolerance)
        {
            done = true;
        }

        int sendDone = done ? 1 : 0;
        int recvDone = 0;
        MPI_Allreduce(&sendDone, &recvDone, 1, MPI_INTEGER, MPI_MIN, comm);
        done = recvDone;
    }
    return runInfo;
}

}

#endif
