#ifndef UNITTESTUTILS_ALGORITHM_TIMER
#define UNITTESTUTILS_ALGORITHM_TIMER

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <iostream>
#include <iomanip>
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

int convert_bool_to_int(bool done)
{
    return done ? 1 : 0;
}

void update_max(double timeForOneRun, RunInfo &runInfo)
{
    if(timeForOneRun > runInfo.max)
    {
        runInfo.max = timeForOneRun;
    }
}
void update_min(double timeForOneRun, RunInfo &runInfo)
{
    if(timeForOneRun < runInfo.min)
    {
        runInfo.min = timeForOneRun;
    }
}

void update_run_info(double timeForOneRun, double duration, RunInfo &runInfo)
{
    runInfo.numRuns++;
    runInfo.mean = duration / runInfo.numRuns;

    update_max(timeForOneRun, runInfo);
    update_min(timeForOneRun, runInfo);
}

template <typename ALGORITHM>
double get_time_to_run_algorithm(const ALGORITHM &alg)
{
    double startTime = stk::cpu_time();
    alg();
    return stk::cpu_time() - startTime;
}

bool has_mean_converged_within_tolerance(const RunInfo &runInfo, const double lastMean, const size_t minRuns, const double tolerance)
{
    bool done = false;
    double meanDiff = std::abs(runInfo.mean - lastMean);
    if(runInfo.numRuns >= minRuns && runInfo.mean > 0.0 && meanDiff < tolerance)
    {
        done = true;
    }
    return done;
}

bool has_converged_on_all_procs(bool done, stk::ParallelMachine comm)
{
    int sendDone = convert_bool_to_int(done);
    int recvDone = 0;
    MPI_Allreduce(&sendDone, &recvDone, 1, MPI_INTEGER, MPI_MIN, comm);
    return recvDone == 1;
}

template <typename ALGORITHM>
RunInfo time_algorithm(const double tolerance, const size_t minRuns, stk::ParallelMachine comm, const ALGORITHM &alg)
{
    initial_run_to_allocate_memory_and_stuff(alg);

    RunInfo runInfo;

    double lastMean = 0.0;
    double duration = 0.0;
    bool isConverged = false;
    while(!isConverged)
    {
        double timeForOneRun = get_time_to_run_algorithm(alg);

        duration += timeForOneRun;
        update_run_info(timeForOneRun, duration, runInfo);

        bool isConvergedOnThisProc = has_mean_converged_within_tolerance(runInfo, lastMean, minRuns, tolerance);
        isConverged = has_converged_on_all_procs(isConvergedOnThisProc, comm);

        lastMean = runInfo.mean;
    }
    return runInfo;
}

void print_run_info(std::ostream &os, const std::string &tag, const int numProcs, const RunInfo &runInfo)
{
    os      << tag << " "
            << std::setw(10) << runInfo.numRuns << " runs"
            << " over " << numProcs << " procs"
            << " mean was "<< std::setw(15) << runInfo.mean
            << " with (min " << runInfo.min << ", max " << runInfo.max << ")" << std::endl;
}

}

#endif
