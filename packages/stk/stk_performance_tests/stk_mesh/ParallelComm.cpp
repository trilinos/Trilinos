#include "gtest/gtest.h"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_unit_test_utils/timer.hpp"
#include <vector>


TEST(ParallelDataExchangePerf, parallel_data_exchange_t_dense)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int commSize = stk::parallel_machine_size(comm);
  int dataSize = 1024 * 1024;
  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS  = 100;

  std::vector<std::vector<int>> sendData(commSize);
  std::vector<std::vector<int>> recvData(commSize);

  for (int i=0; i < commSize; ++i)
  {
    sendData[i].resize(dataSize);
    recvData[i].resize(dataSize);
    for (int j=0; j < dataSize; ++j)
      sendData[i][j] = j;
  }

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();

    for (int i = 0; i < NUM_ITERS; ++i)
    {
      stk::parallel_data_exchange_t(sendData, recvData, comm);
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

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
