#include <gtest/gtest.h>

#include "mpi.h"
#include <stk_util/diag/Timer.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_unit_test_utils/getOption.h>
#include "ImbalanceTimer.hpp"

const float workTimeProc0 = 5;
const float workTimeProcOthers = 1;

void send_data_to_other_procs(const std::vector<double>& data, stk::CommSparse& comm)
{
    int myprocId = comm.parallel_rank();
    for(int i=0;i<comm.parallel_size();++i)
    {
        if(i!=myprocId)
        {
            stk::pack_vector_to_proc(comm, data, i);
        }
    }
}

double calculate_sum(size_t numbersToAdd)
{
    double sum = 0;
    for(size_t i=0;i<numbersToAdd;++i)
        sum += i;
    return sum;
}
void do_sleep(int proc, size_t numDoubles)
{
    stk::CommSparse comm(MPI_COMM_WORLD);

    if(proc==0)
        ::usleep(workTimeProc0*1e6);
    else
        calculate_sum(workTimeProcOthers*1e6);

    std::vector<double> dataSent(numDoubles,1);
    std::vector<std::vector<double> > dataRecvd(comm.parallel_size());
    stk::pack_and_communicate(comm, [&comm, &dataSent]() { send_data_to_other_procs(dataSent, comm); });
    stk::unpack_communications(comm,
                               [&comm, &dataRecvd](int procId)
                               {
                                    stk::unpack_vector_from_proc(comm, dataRecvd[procId], procId);
                               });
}


TEST(LoadImbalanceTimer, unevenWork_imbalanceIsMeasured)
{
    int procId = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    size_t numDoubles = stk::unit_test_util::get_command_line_option<size_t>("-n", 1000u);

    stk::util::ImbalanceTimings imbalanceTimings = stk::util::measure_imbalance(MPI_COMM_WORLD, [=](){
        do_sleep(procId, numDoubles);
    });

    int num_procs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    std::vector<float> timingsGold(num_procs, workTimeProcOthers);
    timingsGold[0] = workTimeProc0;

    int goldProcWithMaxTime = -1;
    float goldMaxTime = 0, goldAverageTime = 0, goldImbalance = 0;
    stk::util::calculate_max_ave_and_imbalance(timingsGold, goldProcWithMaxTime, goldMaxTime, goldAverageTime, goldImbalance);


    EXPECT_NEAR(goldImbalance, imbalanceTimings.get_imbalance(), 0.1);
    EXPECT_NEAR(goldMaxTime, imbalanceTimings.get_max_time(), 0.1);
    EXPECT_EQ(goldProcWithMaxTime, imbalanceTimings.get_proc_with_max_time());
    EXPECT_NEAR(goldAverageTime, imbalanceTimings.get_average_time(), 0.1);

    if(procId == 0)
        imbalanceTimings.print_table_timings(std::cerr);
}

