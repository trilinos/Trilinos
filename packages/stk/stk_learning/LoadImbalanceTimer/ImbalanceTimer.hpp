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

#ifndef IMBALANCETIMER_HPP_
#define IMBALANCETIMER_HPP_

#include <stk_util/diag/Timer.hpp>

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
