#include "gtest/gtest.h"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_performance_tests/stk_mesh/timer.hpp"
#include <vector>


TEST(ParallelDataExchangePerf, parallel_data_exchange_t_dense)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int commSize = stk::parallel_machine_size(comm);
  int dataSize = 1024 * 1024;
  int numRuns  = 100;

  std::vector<std::vector<int>> sendData(commSize);
  std::vector<std::vector<int>> recvData(commSize);

  for (int i=0; i < commSize; ++i)
  {
    sendData[i].resize(dataSize);
    recvData[i].resize(dataSize);
    for (int j=0; j < dataSize; ++j)
      sendData[i][j] = j;
  }

  stk::performance_tests::Timer timer(comm);
  MPI_Barrier(comm);
  timer.start_timing();  

  for (int i=0; i < numRuns; ++i)
  {
    stk::parallel_data_exchange_t(sendData, recvData, comm);
  }

  MPI_Barrier(comm);
  timer.update_timing();
  timer.print_timing(numRuns);

  // Intel MPI
  // np=1: [0.078, 0.110]
  // np=2: [0.70, 0.95]
  // np=4: [1.84, 2.17]
  // np=8: [4.64, 4.87]

  // OpenMPI
  // np=1: [0.09, 0.167]
  // np=2: [0.71, 0.85]
  // np=4: [1.85, 2.09]
  // np=8: [4.25, 4.46]
}