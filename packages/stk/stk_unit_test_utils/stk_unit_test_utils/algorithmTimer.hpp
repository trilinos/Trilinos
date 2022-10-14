// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef UNITTESTUTILS_ALGORITHM_TIMER
#define UNITTESTUTILS_ALGORITHM_TIMER

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stddef.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

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
    double startTime = stk::wall_time();
    alg();
    return stk::wall_time() - startTime;
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
    return stk::is_true_on_all_procs(comm, done);
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

        duration += std::abs(timeForOneRun);
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
