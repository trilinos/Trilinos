#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"
#include <stk_unit_test_utils/timer.hpp>
#include <stk_unit_test_utils/getOption.h>
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingPrepost.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingProbe.hpp"
#include <fstream>
#include <chrono>
#include <thread>

namespace {

template <typename Exchanger>
class ExchangerUnknownPatternTester : public ::testing::Test
{
  public:
    void setupSendBufs(int numDestProcs, int nvalsPerProc)
    {
      int myrank = stk::parallel_machine_rank(comm);
      int commSize = stk::parallel_machine_size(comm);
      
      STK_ThrowRequireMsg(numDestProcs <= commSize, "number of dest procs cannot be greater than comm size");
      
      int firstDestProc = myrank - numDestProcs/2 + 1 - (numDestProcs % 2);
      int lastDestProc = myrank + numDestProcs/2;
      for (int proc=firstDestProc; proc < lastDestProc; ++proc)
      {
        int proc_wrapped = (proc + commSize) % commSize;
        EXPECT_TRUE(proc_wrapped >= 0 && proc_wrapped < commSize);
        destProcs.push_back(proc_wrapped);
      }
      
      sendBufs.resize(commSize);
      recvBufs.resize(commSize);
      for (int proc : destProcs)
        sendBufs[proc].resize(nvalsPerProc, 42);
    }
    
    stk::ParallelMachine comm = stk::parallel_machine_world();
    std::vector<int> destProcs;
    std::vector<std::vector<int>> sendBufs;
    std::vector<std::vector<int>> recvBufs;
};

using ExchangerTypes = ::testing::Types<stk::DataExchangeUnknownPatternNonBlocking,
                                        stk::DataExchangeUnknownPatternNonBlockingProbe,
                                        stk::DataExchangeUnknownPatternNonBlockingPrepost>;

class NameGenerator
{
  public:
    template <typename T>
    static std::string GetName(int)
    {
      if constexpr (std::is_same_v<T, stk::DataExchangeUnknownPatternNonBlockingProbe>)
      {
        return "DataExchangeUnknownPatternNonBlockingProbe";
      } else if constexpr (std::is_same_v<T, stk::DataExchangeUnknownPatternNonBlockingPrepost>)
      {
        return "DataExchangeUnknownPatternNonBlockingPrepost";        
      } else
      {
        return "OtherDataExchanger";
      }
    }
}; 


}

TYPED_TEST_SUITE(ExchangerUnknownPatternTester, ExchangerTypes, NameGenerator);


TYPED_TEST(ExchangerUnknownPatternTester, Perf)
{
  using ExchangerType = TypeParam;
  const unsigned NUM_BATCHES = 5;
  const unsigned NUM_RUNS = stk::unit_test_util::get_command_line_option("-r", 5);
  const unsigned NUM_DEST_PROCS = stk::unit_test_util::get_command_line_option("-d", 8);
  const unsigned NUM_VALUES_PER_PROC = stk::unit_test_util::get_command_line_option("-v", 20);
  const double SLEEP_FAC = stk::unit_test_util::get_command_line_option("-s", 0.0);
  
  if (stk::parallel_machine_rank(this->comm) == 0)
  {
    std::cout << "using " << NameGenerator::GetName<ExchangerType>(0) << " for " << NUM_RUNS << " round of communciation, sending " << NUM_VALUES_PER_PROC*sizeof(int)/(1024*1024)
              << " Megabytes of data to " << NUM_DEST_PROCS << " procs" << std::endl;
  }

  std::chrono::duration<double> sleep_time(SLEEP_FAC * stk::parallel_machine_rank(this->comm));
  stk::unit_test_util::BatchTimer batchTimer(this->comm);
  batchTimer.initialize_batch_timer();
  
  for(unsigned b=0; b<NUM_BATCHES; ++b)
  {
    batchTimer.start_batch_timer();
    ExchangerType exchanger(this->comm);
    this->setupSendBufs(NUM_DEST_PROCS, NUM_VALUES_PER_PROC);

    for(unsigned r=0; r<NUM_RUNS; ++r)
    {
      exchanger.start_nonblocking(this->sendBufs, this->recvBufs);
      std::this_thread::sleep_for(sleep_time);
      exchanger.post_nonblocking_receives(this->recvBufs);
      exchanger.complete_receives(this->recvBufs, [](int /*rank*/, const std::vector<int>& /*buf*/) {});
      exchanger.complete_sends();
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_RUNS);
  
  if (stk::parallel_machine_rank(this->comm) == 0)
  {
    std::ofstream of("unknown_pattern_results.txt", std::ios::app);
    of << std::setprecision(16);
    of << stk::parallel_machine_size(this->comm) << " " << NUM_DEST_PROCS << " " << NUM_VALUES_PER_PROC*sizeof(int)
       << " " << SLEEP_FAC << " " << NameGenerator::GetName<ExchangerType>(0) << " " << batchTimer.get_min_batch_time() << "\n";
    of.close();
  }
  
  stk::parallel_machine_barrier(this->comm);
}

