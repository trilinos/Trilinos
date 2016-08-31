
#ifndef IMBALANCETIMER_HPP_
#define IMBALANCETIMER_HPP_

namespace stk { namespace util {

inline
void calculate_max_ave_and_imbalance(const std::vector<float> &allTimes, int &procWithMaxTime, float &maxTime, float &averageTime, float &imbalance)
{
    for(size_t i = 0; i < allTimes.size(); ++i)
    {
        averageTime += allTimes[i];
        if(maxTime < allTimes[i])
        {
            maxTime = allTimes[i];
            procWithMaxTime = i;
        }
    }

    averageTime /= allTimes.size();
    imbalance = maxTime/averageTime;
}

struct ImbalanceTimings
{
    ImbalanceTimings(MPI_Comm communicator, float time)
    : myProc(-1), imbalance(0), averageTime(0), maxTime(0), procWithMaxTime(-1)
    {
        int numProcs = 0;
        MPI_Comm_size(communicator, &numProcs);
        MPI_Comm_rank(communicator, &myProc);

        allTimes.resize(numProcs, 0);

        MPI_Allgather(&time, 1, MPI_FLOAT, allTimes.data(), 1, MPI_FLOAT, communicator);

        calculate_max_ave_and_imbalance(allTimes, procWithMaxTime, maxTime, averageTime, imbalance);
    }

    float get_imbalance() const { return imbalance; }
    int get_proc_with_max_time() const { return procWithMaxTime; }
    float get_max_time() const { return maxTime; }
    float get_average_time() const { return averageTime; }

    void print_table_timings(std::ostream& os) const
    {
        if (myProc == 0)
        {
            std::ostringstream out;
            out << "imbalance: " << imbalance << std::endl;
            out << "Proc\tTime" << std::endl;
            for(size_t i=0;i<allTimes.size();++i)
                out << i << "\t" << allTimes[i] << std::endl;
            os << out.str();
        }
    }

private:
    int myProc;
    float imbalance;
    float averageTime;
    float maxTime;
    int procWithMaxTime;
    std::vector<float> allTimes;
};

template<class FUNCTION>
ImbalanceTimings measure_imbalance(MPI_Comm communicator, FUNCTION do_work)
{
    enum {CHILDMASK1 = 1};
    stk::diag::TimerSet enabledTimerSet(CHILDMASK1);
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
    rootTimer.start();
    stk::diag::Timer childTimer1("childTimer1", CHILDMASK1, rootTimer);

    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(childTimer1, communicator);
        do_work();
    }

    float time = childTimer1.getMetric<stk::diag::WallTime>().getAccumulatedLap(false);
    ImbalanceTimings timings(communicator, time);
    stk::diag::deleteRootTimer(rootTimer);

    return timings;
}

}}

#endif /* IMBALANCETIMER_HPP_ */
