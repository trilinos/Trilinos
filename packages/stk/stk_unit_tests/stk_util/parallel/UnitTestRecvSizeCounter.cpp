#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ReceiveSizeCounter.hpp"

#ifdef STK_HAS_MPI

TEST(RecvSizeCounter, AllToAllNonblocking)
{
  int commRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int commSize = stk::parallel_machine_size(stk::parallel_machine_world());

  std::vector<std::vector<double>> send_bufs(commSize);
  for (int recvRank=0; recvRank < commSize; ++recvRank)
  {
    send_bufs[recvRank].resize(commRank + recvRank);
  }

  stk::ReceiveSizeCounter counter(stk::parallel_machine_world());
  counter.start(send_bufs, false);
  while (!counter.is_complete()) {}

  const std::vector<stk::ReceiveSizeCounter::ULL>& recvSizes = counter.get_receive_sizes();

  for (int senderRank=0; senderRank < commSize; ++senderRank)
  {
    EXPECT_EQ(recvSizes[senderRank], size_t(senderRank + commRank));
  }
}

TEST(RecvSizeCounter, AllToAllBlocking)
{
  int commRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int commSize = stk::parallel_machine_size(stk::parallel_machine_world());

  std::vector<std::vector<double>> send_bufs(commSize);
  for (int recvRank=0; recvRank < commSize; ++recvRank)
  {
    send_bufs[recvRank].resize(commRank + recvRank);
  }

  stk::ReceiveSizeCounter counter(stk::parallel_machine_world());
  counter.start(send_bufs, true);
  EXPECT_TRUE(counter.is_complete());

  const std::vector<stk::ReceiveSizeCounter::ULL>& recvSizes = counter.get_receive_sizes();

  for (int senderRank=0; senderRank < commSize; ++senderRank)
  {
    EXPECT_EQ(recvSizes[senderRank], size_t(senderRank + commRank));
  }
}

TEST(RecvSizeCounter, AllToAllNonblockingRepeat)
{
  int commRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int commSize = stk::parallel_machine_size(stk::parallel_machine_world());
  stk::ReceiveSizeCounter counter(stk::parallel_machine_world());

  for (int round=0; round < 5; ++round)
  {

    std::vector<std::vector<double>> send_bufs(commSize);
    for (int recvRank=0; recvRank < commSize; ++recvRank)
    {
      send_bufs[recvRank].resize(commRank + recvRank + round);
    }

    counter.start(send_bufs, false);
    while (!counter.is_complete()) {}

    const std::vector<stk::ReceiveSizeCounter::ULL>& recvSizes = counter.get_receive_sizes();

    for (int senderRank=0; senderRank < commSize; ++senderRank)
    {
      EXPECT_EQ(recvSizes[senderRank], size_t(senderRank + commRank + round));
    }
  }
}

TEST(RecvSizeCounter, NeighorNonblocking)
{
  int commRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int commSize = stk::parallel_machine_size(stk::parallel_machine_world());

  int destRank = (commRank + 1) % commSize;
  int sendRank = (commRank - 1 + commSize) % commSize;

  std::vector<std::vector<double>> send_bufs(commSize);
  send_bufs[destRank].resize(commRank);

  stk::ReceiveSizeCounter counter(stk::parallel_machine_world());
  counter.start(send_bufs, false);
  while (!counter.is_complete()) {}

  const std::vector<stk::ReceiveSizeCounter::ULL>& recvSizes = counter.get_receive_sizes();

  for (int rank=0; rank < commSize; ++rank)
  {
    if (rank == sendRank)
    {
      EXPECT_EQ(recvSizes[rank], size_t(sendRank));
    } else
    {
      EXPECT_EQ(recvSizes[rank], 0U);
    }
  }
}

TEST(RecvSizeCounter, NeighorBlocking)
{
  int commRank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int commSize = stk::parallel_machine_size(stk::parallel_machine_world());

  int destRank = (commRank + 1) % commSize;
  int sendRank = (commRank - 1 + commSize) % commSize;

  std::vector<std::vector<double>> send_bufs(commSize);
  send_bufs[destRank].resize(commRank);

  stk::ReceiveSizeCounter counter(stk::parallel_machine_world());
  counter.start(send_bufs, true);
  EXPECT_TRUE(counter.is_complete());

  const std::vector<stk::ReceiveSizeCounter::ULL>& recvSizes = counter.get_receive_sizes();

  for (int rank=0; rank < commSize; ++rank)
  {
    if (rank == sendRank)
    {
      EXPECT_EQ(recvSizes[rank], size_t(sendRank));
    } else
    {
      EXPECT_EQ(recvSizes[rank], 0U);
    }
  }
}

#endif