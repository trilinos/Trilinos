#ifndef UNITTEST_PERFORMANCETESTER_HPP
#define UNITTEST_PERFORMANCETESTER_HPP

#include <stk_mesh/base/Comm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>

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
            enabledTimerSet(CHILDMASK1),
            rootTimer(createRootTimer("totalTestRuntime", enabledTimerSet)),
            childTimer("timed algorithm", CHILDMASK1, rootTimer),
            communicator(comm),
            duration(0.0)
    {
        rootTimer.start();
    }

    virtual ~PerformanceTester()
    {
        stk::diag::deleteRootTimer(rootTimer);
    }

    virtual void run_algorithm_to_time() = 0;
    virtual size_t get_value_to_output_as_iteration_count() = 0;

protected:
    const int CHILDMASK1 = 1;
    stk::diag::TimerSet enabledTimerSet;
    stk::diag::Timer rootTimer;
    stk::diag::Timer childTimer;
    MPI_Comm communicator;

    double duration;

private:
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
