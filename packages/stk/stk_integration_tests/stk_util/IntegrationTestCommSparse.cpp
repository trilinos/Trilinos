#include "gtest/gtest.h"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternUserDataNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlocking.hpp"

namespace {
using DataType = long long;
// make message large enough to not fit in any int (signed or unsigned)
constexpr size_t approxMsgSizeInBytes = (size_t(1) << 32) + 2*sizeof(DataType);  

constexpr size_t numElements = approxMsgSizeInBytes / sizeof(DataType);
constexpr size_t msgSizeInBytes = numElements * sizeof(DataType);
}

TEST(CommSparse, LargerThan2GB)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  stk::CommSparse comm(stk::parallel_machine_world());
  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  constexpr int senderRank = 0;
  constexpr int receiverRank = 1;
  
  for (int phase=0; phase < 2; ++phase)
  {
    if (myrank == senderRank)
    {
      for (size_t i=0; i < numElements; ++i)
      {
        comm.send_buffer(receiverRank).pack(DataType(i));
      }
    }
    
    if (phase == 0)
    {
      comm.allocate_buffers();
      if (myrank == senderRank)
      {
        EXPECT_EQ(comm.send_buffer(receiverRank).capacity(), msgSizeInBytes);
      }
    } else
    {
      comm.communicate();
    }
  }
  
  if (myrank == 1)
  {
    EXPECT_EQ(comm.recv_buffer(senderRank).capacity(), msgSizeInBytes);
    for (size_t i=0; i < numElements; ++i)
    {
      DataType val;
      comm.recv_buffer(senderRank).unpack(val);
      EXPECT_EQ(val, DataType(i));
    }
  }
}

TEST(DataExchangeUnknownPatternNonBlocking, LargerThan2GB)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  stk::DataExchangeUnknownPatternNonBlocking comm(stk::parallel_machine_world());
  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  constexpr int senderRank = 0;
  constexpr int receiverRank = 1;
  
  for (int repeat=0; repeat < 3; ++repeat)
  {
    std::vector<std::vector<DataType>> sendBufs(2);
    std::vector<std::vector<DataType>> recvBufs(2);
    if (myrank == senderRank)
    {
      sendBufs[receiverRank].resize(numElements);
      for (size_t i=0; i < numElements; ++i)
      {
        sendBufs[receiverRank][i] = i;
      }
    }
    
    comm.start_nonblocking(sendBufs, recvBufs);
    comm.post_nonblocking_receives(recvBufs);
    
    auto unpacker = [&](int rank, const std::vector<DataType>& recv_buf)
    {
      EXPECT_EQ(rank, senderRank);
      EXPECT_EQ(recv_buf.size(), numElements);
      for (size_t i=0; i < numElements; ++i)
      {
        EXPECT_EQ(recv_buf[i], DataType(i));
      }
    };
    
    comm.complete_receives(recvBufs, unpacker);
    comm.complete_sends();
  }
}

TEST(DataExchangeKnownPatternNonBlocking, LargerThan2GB)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  stk::DataExchangeKnownPatternNonBlocking comm(stk::parallel_machine_world());
  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  constexpr int senderRank = 0;
  constexpr int receiverRank = 1;
  
  for (int repeat = 0; repeat < 3; ++repeat)
  {
    std::vector<std::vector<DataType>> sendBufs(2);
    std::vector<std::vector<DataType>> recvBufs(2);
    if (myrank == senderRank)
    {
      sendBufs[receiverRank].reserve(numElements);
      for (size_t i=0; i < numElements; ++i)
      {
        sendBufs[receiverRank].push_back(i);
      }
    } else
    {
      recvBufs[senderRank].resize(numElements);
    }
    
    comm.start_nonblocking(sendBufs, recvBufs);
    
    auto unpacker = [&](int rank, const std::vector<DataType>& recv_buf)
    {
      EXPECT_EQ(rank, senderRank);
      EXPECT_EQ(recv_buf.size(), numElements);
      for (size_t i=0; i < numElements; ++i)
      {
        EXPECT_EQ(recv_buf[i], DataType(i));
      }
    };
    
    comm.complete_receives(recvBufs, unpacker);
    comm.complete_sends();
  }
}

TEST(DataExchangeKnownPatternUserDataNonBlocking, LargerThan2GB)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  stk::DataExchangeKnownPatternUserDataNonBlocking comm(stk::parallel_machine_world());
  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  constexpr int senderRank = 0;
  constexpr int receiverRank = 1;
  
  for (int repeat = 0; repeat < 3; ++repeat)
  {
    std::vector<DataType> allSendBufs(numElements), allRecvBufs(numElements);
    std::vector<stk::PointerAndSize> sendBufPointers, recvBufPointers;
    std::vector<int> sendRanks, recvRanks;
    if (myrank == senderRank)
    {
      for (size_t i=0; i < numElements; ++i)
      {
        allSendBufs[i] = i;
      }
      sendBufPointers.emplace_back(reinterpret_cast<unsigned char*>(allSendBufs.data()),
                                   numElements * sizeof(DataType));
      sendRanks.push_back(receiverRank);
    } else
    {
      recvBufPointers.emplace_back(reinterpret_cast<unsigned char*>(allRecvBufs.data()),
                                   numElements * sizeof(DataType));
      recvRanks.push_back(senderRank);      
    }
    
    comm.start_nonblocking(sendBufPointers, sendRanks, recvBufPointers, recvRanks);
    
    auto unpacker = [&](int rank, const stk::PointerAndSize& recv_buf)
    {
      EXPECT_EQ(rank, senderRank);
      EXPECT_EQ(recv_buf.size, numElements * sizeof(DataType));
      const DataType* ptr = reinterpret_cast<const DataType*>(recv_buf.ptr);
      for (size_t i=0; i < numElements; ++i)
      {
        EXPECT_EQ(ptr[i], DataType(i));
        EXPECT_EQ(allRecvBufs[i], DataType(i));
      }
    };
    
    comm.complete_receives(recvBufPointers, recvRanks, unpacker);
    comm.complete_sends();
  }
}

