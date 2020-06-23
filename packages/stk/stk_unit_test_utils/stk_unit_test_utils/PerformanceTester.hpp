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

#ifndef UNITTEST_PERFORMANCETESTER_HPP
#define UNITTEST_PERFORMANCETESTER_HPP

#include <stk_mesh/base/Comm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>

namespace stk
{
namespace unit_test_util
{

inline void print_output_for_pass_fail_test(double duration, unsigned iterCount, MPI_Comm communicator)
{
    std::ofstream out("forPassFailScript.log");
    double maxTime = stk::get_max_time_across_procs(duration, communicator);
    double maxHwmInMB = stk::get_max_hwm_across_procs(communicator) / (1024.0 * 1024.0);
    stk::print_stats_for_performance_compare(out, maxTime, maxHwmInMB, iterCount, communicator);
}

inline void print_output_for_graph_generation(double duration, const stk::diag::Timer &rootTimer, MPI_Comm communicator)
{
    std::ofstream out("forGraphs.log");
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(out, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint, communicator);
    stk::parallel_print_time_without_output_and_hwm(communicator, duration, out);
}

class PerformanceTester
{
public:
    void run_performance_test()
    {
        time_algorithm();
        generate_output();
    }

protected:
    PerformanceTester(MPI_Comm comm) :
            duration(0.0),
            enabledTimerSet(CHILDMASK1),
            rootTimer(createRootTimer("totalTestRuntime", enabledTimerSet)),
            childTimer("timed algorithm", CHILDMASK1, rootTimer),
            communicator(comm)
    {
        rootTimer.start();
    }

    virtual ~PerformanceTester()
    {
        stk::diag::deleteRootTimer(rootTimer);
    }

    virtual void run_algorithm_to_time() = 0;
    virtual size_t get_value_to_output_as_iteration_count() = 0;

    double duration;

private:
    const int CHILDMASK1 = 1;
    stk::diag::TimerSet enabledTimerSet;
    stk::diag::Timer rootTimer;
    stk::diag::Timer childTimer;
    MPI_Comm communicator;

    void time_algorithm()
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(childTimer, communicator);
        double startTime = stk::wall_time();
        run_algorithm_to_time();
        duration += stk::wall_time() - startTime;
    }

    void generate_output()
    {
        print_output_for_pass_fail_test(duration, get_value_to_output_as_iteration_count(), communicator);
        print_output_for_graph_generation(duration, rootTimer, communicator);
    }
};

}
}
#endif
