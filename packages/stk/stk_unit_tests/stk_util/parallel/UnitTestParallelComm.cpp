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
// 

#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"      // for parallel_machine_rank, parallel_machine_size

#ifdef STK_HAS_MPI

#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternBlockingBuffer.hpp"

#include <vector>
#include <algorithm>
#include <fstream>


namespace {

using DataExchangeBlocking              = stk::DataExchangeUnknownPatternBlocking;
using DataExchangeNonBlocking           = stk::DataExchangeUnknownPatternNonBlocking;
template <typename T>
using DataExchangeBlockingBuffer        = stk::DataExchangeUnknownPatternBlockingBuffer<T>;
template <typename T>
using DataExchangeNonBlockingBuffer     = stk::DataExchangeUnknownPatternNonBlockingBuffer<T>;
using DataExchangeBlockingCommBuffer    = stk::DataExchangeUnknownPatternBlockingCommBuffer;
using DataExchangeNonBlockingCommBuffer = stk::DataExchangeUnknownPatternNonBlockingCommBuffer;


template <typename T>
std::vector<int> get_ranks(std::vector< std::vector<T> >& lists)
{
  std::vector<int> ranks;
  for (unsigned int i=0; i < lists.size(); ++i) {
    if (lists[i].size() > 0) {
      ranks.push_back(i);
    }
  }

  return ranks;
}

template <typename T>
std::vector<std::vector<T>> copySendBufs(const stk::ManagedBufferBase<T>& exchanger)
{
  int commSize = stk::parallel_machine_size(exchanger.get_comm());
  std::vector< std::vector<T> > sendBufs(commSize);

  for (int i=0; i < commSize; ++i)
  {
    auto& sendBuf = exchanger.get_send_buf(i);
    sendBufs[i].assign(sendBuf.begin(), sendBuf.end());
  }

  return sendBufs;
}


template <typename T>
std::vector<std::vector<T>> copySendBufs(stk::ManagedCommBufferBase& exchanger)
{
  int commSize = stk::parallel_machine_size(exchanger.get_comm());
  std::vector< std::vector<T> > sendBufs(commSize);

  for (int i=0; i < commSize; ++i)
  {
    stk::CommBuffer& sendBuf = exchanger.get_send_buf(i);
    sendBuf.reset();
    size_t nvals = sendBuf.capacity() / sizeof(T);
    for (size_t j=0; j < nvals; ++j)
    {
      T val;
      sendBuf.unpack(val);
      sendBufs[i].push_back(val);
    }
  }

  return sendBufs;
}


template <typename T>
std::vector<std::vector<T>> copyRecvBufs(const stk::ManagedBufferBase<T>& exchanger)
{
  int commSize = stk::parallel_machine_size(exchanger.get_comm());
  std::vector< std::vector<T> > recvBufs(commSize);

  for (int i=0; i < commSize; ++i)
  {
    auto& recvBuf = exchanger.get_recv_buf(i);
    recvBufs[i].assign(recvBuf.begin(), recvBuf.end());
  }

  return recvBufs;
}


template <typename T>
std::vector<std::vector<T>> copyRecvBufs(stk::ManagedCommBufferBase& exchanger)
{
  int commSize = stk::parallel_machine_size(exchanger.get_comm());
  std::vector< std::vector<T> > recvBufs(commSize);

  for (int i=0; i < commSize; ++i)
  {
    stk::CommBuffer& recvBuf = exchanger.get_recv_buf(i);
    size_t nvals = recvBuf.capacity() / sizeof(T);
    for (size_t j=0; j < nvals; ++j)
    {
      T val;
      recvBuf.unpack(val);
      recvBufs[i].push_back(val);
    }
  }

  return recvBufs;
}

template <typename T>
class ParallelCommTester
{
  public:
    using ElementType = T;

    ParallelCommTester(bool split_last = false)
    {
      
      myrank    = stk::parallel_machine_rank(MPI_COMM_WORLD);
      commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
      if (split_last) {
        color = myrank ==  commSize -1 ? 0 : 1;
        int key = myrank;
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);

        myrank    = stk::parallel_machine_rank(comm);
        commSize = stk::parallel_machine_size(comm);
      } else {
        comm = MPI_COMM_WORLD;
        color = 0;
      }
    }

    virtual ~ParallelCommTester() {}

    int get_value(const int rankSrc, const int rankDest, const int idx)
    {
        return rankSrc + 2*rankDest + 3*idx + m_offset;
    }

    int get_size(const int rankSrc, const int rankDest)
    {
        return m_nvalsBase + rankSrc + 2*rankDest;
    }

    virtual void set_offset(int val) { m_offset = val; }

    std::vector< std::vector<T> >& getSendLists() { return sendLists; }

    std::vector< std::vector<T> >& getReceiveLists() { return recvLists; }

    MPI_Comm comm;
    int color = -1;
    int myrank = -1;
    int commSize = -1;

  protected:
    std::vector< std::vector<T> > sendLists;
    std::vector< std::vector<T> > recvLists;

  private:
    const int m_nvalsBase = 512;
    int m_offset = 0;

};

// test communication when all processes send to all other processes
template <typename T>
class DenseParallelCommTesterBase : public ParallelCommTester<T>
{
  public:
    using ParallelCommTester<T>::sendLists;
    using ParallelCommTester<T>::recvLists;
    using ParallelCommTester<T>::comm;
    using ParallelCommTester<T>::myrank;
    using ParallelCommTester<T>::commSize;

    explicit DenseParallelCommTesterBase(bool split_last = false)
    {
      recvLists.resize(commSize);
      sendLists.resize(commSize);

      set_send_buffers_values();
    }

    virtual ~DenseParallelCommTesterBase() {}

    void set_offset(int val) override
    {
      ParallelCommTester<T>::set_offset(val);
      set_send_buffers_values();
    }

    void test_results(std::vector< std::vector<T> >& recvLists)
    {
      for (int src=0; src < commSize; ++src) {
        test_recv_vals(recvLists[src], src);
      }
    }

    void test_recv_vals(const std::vector<T>& recvList, int src)
    {
      size_t sizeEx = this->get_size(src, myrank);
      EXPECT_EQ(recvList.size(), sizeEx);

      if (recvList.size() == sizeEx) {
        for (size_t j=0; j < sizeEx; ++j) {
          EXPECT_EQ(recvList[j], this->get_value(src, myrank, j));
        }
      }
    }

    void test_send_ranks(std::vector< std::vector<T> >& sendLists)
    {
      test_ranks_inner(sendLists);
    }


    void test_recv_ranks(std::vector<std::vector<T>>& recvLists)
    {
      test_ranks_inner(recvLists);
    }

  private:
    void test_ranks_inner(std::vector< std::vector<T> >& buffers)
    {
      auto ranks = get_ranks(buffers);
      int len = ranks.size();
      EXPECT_EQ(len, commSize);

      for (int i=0; i < len; ++i) {
        EXPECT_EQ(ranks[i], i);
      }
    }

    void set_send_buffers_values()
    {
      for (int dest=0; dest < commSize; ++dest) {
        int size_d = this->get_size(myrank, dest);
        sendLists[dest].resize(size_d);
        for (int j=0; j < size_d; ++j) {
          sendLists[dest][j] = this->get_value(myrank, dest, j);
        }
      }
    }
};


template <typename T>
class DenseParallelCommTester : public ::testing::Test,
                                public DenseParallelCommTesterBase<T>
{
};

using DenseParallelCommTesterInt = DenseParallelCommTester<int>;
using DenseParallelCommTesterDouble = DenseParallelCommTester<double>;

// test when process i only sends to i + 1 and i + 2 (with wrap around)
template <typename T>
class NeighborParallelCommTesterBase : public ParallelCommTester<T>
{
  public:
    using ParallelCommTester<T>::sendLists;
    using ParallelCommTester<T>::recvLists;
    using ParallelCommTester<T>::comm;
    using ParallelCommTester<T>::myrank;
    using ParallelCommTester<T>::commSize;

    NeighborParallelCommTesterBase()
    {
      recvLists.resize(commSize);
      sendLists.resize(commSize);

      set_send_buffers_values();
    }

    virtual ~NeighborParallelCommTesterBase() {}

    void set_offset(int val) override
    {
      ParallelCommTester<T>::set_offset(val);
      set_send_buffers_values();
    }

    void test_results(std::vector< std::vector<T> >& recvLists)
    {
      for (int src=0; src < commSize; ++src) {
        test_recv_vals(recvLists[src], src);
      }
    }

    void test_recv_vals(const std::vector<T>& recvVals, int src)
    {
      int src1 = (myrank - 1 + commSize) % commSize;
      int src2 = (myrank - 2 + commSize) % commSize;
      size_t sizeEx = 0;
      if (src == src1 || src == src2) {
        sizeEx = this->get_size(src, myrank);
      }

      EXPECT_EQ(recvVals.size(), sizeEx);
      for (size_t j=0; j < sizeEx; ++j) {
        EXPECT_EQ(recvVals[j], this->get_value(src, myrank, j));
      }
    }

    void test_send_ranks(std::vector<std::vector<T>>& sendLists)
    {

      std::vector<int> sendRanks = get_ranks(sendLists);

      int len = sendRanks.size();
      int dest1 = (myrank + 1) % commSize;
      int dest2 = (myrank + 2) % commSize;
      int nOthers = 2;
      if (myrank == dest1) {
        nOthers = 1;
      }

      EXPECT_EQ(len, nOthers);

      std::vector<int> sendRanksEx{dest1, dest2};
      std::sort(sendRanksEx.begin(), sendRanksEx.end());
      for (int i=0; i < nOthers; ++i) {
        EXPECT_EQ(sendRanks[i], sendRanksEx[i]);
      }
    }

    void test_recv_ranks(std::vector<std::vector<T>>& recvLists)
    {
      auto recvRanks = get_ranks(recvLists);

      int len = recvRanks.size();
      int dest1 = (myrank - 1 + commSize) % commSize;
      int dest2 = (myrank - 2 + commSize) % commSize;
      int nOthers = 2;
      if (dest1 == myrank && dest2 == myrank) {
        nOthers = 1;
      }

      EXPECT_EQ(len, nOthers);

      std::vector<int> recvRanksEx{dest1, dest2};
      std::sort(recvRanksEx.begin(), recvRanksEx.end());
      for (int i=0; i < nOthers; ++i) {
        EXPECT_EQ(recvRanks[i], recvRanksEx[i]);
      }
    }

  private:
    void set_send_buffers_values()
    {
      int dest1 = (myrank + 1) % commSize;
      int dest2 = (myrank + 2) % commSize;
      int size_d1 = this->get_size(myrank, dest1);
      int size_d2 = this->get_size(myrank, dest2);
      sendLists[dest1].resize(size_d1);
      sendLists[dest2].resize(size_d2);
      for (int j=0; j < size_d1; ++j) {
        sendLists[dest1][j] = this->get_value(myrank, dest1, j);
      }

      for (int j=0; j < size_d2; ++j) {
        sendLists[dest2][j] = this->get_value(myrank, dest2, j);
      }
    }

};

template <typename T>
class NeighborParallelCommTester : public ::testing::Test,
                                   public NeighborParallelCommTesterBase<T>
{
};

using NeighborParallelCommTesterInt = NeighborParallelCommTester<int>;
using NeighborParallelCommTesterDouble = NeighborParallelCommTester<double>;


}

//-----------------------------------------------------------------------------
// Test parallel_data_exchange_t

TEST_F(DenseParallelCommTesterInt, FuncBlocking)
{
  for (int i=0; i < 100; ++i) {
    stk::parallel_data_exchange_t(getSendLists(), getReceiveLists(), comm);
    test_results(getReceiveLists());
  }
}


TEST_F(DenseParallelCommTesterDouble, FuncBlocking)
{
  for (int i=0; i < 100; ++i) {
    stk::parallel_data_exchange_t(getSendLists(), getReceiveLists(), comm);
    test_results(getReceiveLists());
  }
}


TEST_F(NeighborParallelCommTesterInt, FuncBlocking)
{
  stk::parallel_data_exchange_t(getSendLists(), getReceiveLists(), comm);
  test_results(getReceiveLists());
}


TEST_F(NeighborParallelCommTesterDouble, FuncBlocking)
{
  stk::parallel_data_exchange_t(getSendLists(), getReceiveLists(), comm);
  test_results(getReceiveLists());
}

//-----------------------------------------------------------------------------
// test ParallelDataExchange class

TEST_F(DenseParallelCommTesterInt, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    exchanger.execute(getSendLists(), getReceiveLists());
    test_results(getReceiveLists());
    test_send_ranks(getSendLists());
    test_recv_ranks(getReceiveLists());
  }

}

TEST_F(DenseParallelCommTesterInt, ClassBlockingBuffer)
{
  DataExchangeBlockingBuffer<int> exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int i=0; i < commSize; ++i)
    {
      auto& buf_in = getSendLists()[i];
      auto& buf_out = exchanger.get_send_buf(i);
      buf_out.assign(buf_in.begin(), buf_in.end());
    }

    exchanger.execute();
    
    std::vector<std::vector<int>> send_bufs = copySendBufs(exchanger);
    std::vector<std::vector<int>> recv_bufs = copyRecvBufs(exchanger);

    test_results(recv_bufs);
    test_send_ranks(send_bufs);
    test_recv_ranks(recv_bufs);

    exchanger.clear_send_bufs();
    exchanger.clear_recv_bufs();
  }
}

TEST_F(DenseParallelCommTesterInt, ClassBlockingCommBuffer)
{
  DataExchangeBlockingCommBuffer exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int phase=0; phase < 2; ++phase) {
      for (int i=0; i < commSize; ++i) {
        auto& send_buf = exchanger.get_send_buf(i);
        for (auto& val : getSendLists()[i]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger.allocate_send_buffers();
      }
    }

    exchanger.execute();
    
    std::vector<std::vector<int>> send_bufs = copySendBufs<int>(exchanger);
    std::vector<std::vector<int>> recv_bufs = copyRecvBufs<int>(exchanger);

    test_results(recv_bufs);
    test_send_ranks(send_bufs);
    test_recv_ranks(recv_bufs);

    exchanger.clear_send_bufs();
    exchanger.clear_recv_bufs();
  }
}


TEST_F(DenseParallelCommTesterDouble, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    exchanger.execute(getSendLists(), getReceiveLists());
    test_results(getReceiveLists());
    test_send_ranks(getSendLists());
    test_recv_ranks(getReceiveLists());
  }
}


TEST_F(NeighborParallelCommTesterInt, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);
  for (int i=0; i < 100; ++i) {
    set_offset(i);
    exchanger.execute(getSendLists(), getReceiveLists());
    test_results(getReceiveLists());
    test_send_ranks(getSendLists());
    test_recv_ranks(getReceiveLists());
  }

}

TEST_F(NeighborParallelCommTesterInt, ClassBlockingBuffer)
{
  DataExchangeBlockingBuffer<int> exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int i=0; i < commSize; ++i)
    {
      auto& buf_in = getSendLists()[i];
      auto& buf_out = exchanger.get_send_buf(i);
      buf_out.assign(buf_in.begin(), buf_in.end());
    }

    exchanger.execute();
    
    std::vector<std::vector<int>> send_bufs = copySendBufs(exchanger);
    std::vector<std::vector<int>> recv_bufs = copyRecvBufs(exchanger);

    test_results(recv_bufs);
    test_send_ranks(send_bufs);
    test_recv_ranks(recv_bufs);

    exchanger.clear_send_bufs();
    exchanger.clear_recv_bufs();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassBlockingCommBuffer)
{
  DataExchangeBlockingCommBuffer exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int phase=0; phase < 2; ++phase) {
      for (int i=0; i < commSize; ++i) {
        auto& send_buf = exchanger.get_send_buf(i);
        for (auto& val : getSendLists()[i]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger.allocate_send_buffers();
      }
    }

    exchanger.execute();
    
    std::vector<std::vector<int>> send_bufs = copySendBufs<int>(exchanger);
    std::vector<std::vector<int>> recv_bufs = copyRecvBufs<int>(exchanger);

    test_results(recv_bufs);
    test_send_ranks(send_bufs);
    test_recv_ranks(recv_bufs);

    exchanger.clear_send_bufs();
    exchanger.clear_recv_bufs();
  }
}


TEST_F(NeighborParallelCommTesterDouble, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);
  for (int i=0; i < 100; ++i) {
    set_offset(i);
    exchanger.execute(getSendLists(), getReceiveLists());
    test_results(getReceiveLists());
    test_send_ranks(getSendLists());
    test_recv_ranks(getReceiveLists());
  }
}


TEST_F(DenseParallelCommTesterInt, ClassNonBlocking)
{
  DataExchangeNonBlocking exchanger1(comm);
  DataExchangeNonBlocking exchanger2(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());
  EXPECT_FALSE(exchanger2.are_recvs_in_progress());
  EXPECT_FALSE(exchanger2.are_sends_in_progress());

  auto& sendLists = getSendLists();
  auto& recvLists = getReceiveLists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    exchanger1.start_nonblocking(sendLists, recvLists);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_TRUE(exchanger2.are_recvs_in_progress());
    EXPECT_TRUE(exchanger2.are_sends_in_progress());
    exchanger1.post_nonblocking_receives(recvLists);
    exchanger2.post_nonblocking_receives(recvLists2);
    test_send_ranks(sendLists);
    test_recv_ranks(recvLists);

    test_send_ranks(sendLists2);
    test_recv_ranks(recvLists2);

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvLists, f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    test_results(recvLists);

    exchanger2.complete_receives(recvLists2, f);
    EXPECT_FALSE(exchanger2.are_recvs_in_progress());
    EXPECT_TRUE(exchanger2.are_sends_in_progress());
    test_results(recvLists2);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());

    exchanger2.complete_sends();
    EXPECT_FALSE(exchanger2.are_sends_in_progress());
    EXPECT_FALSE(exchanger2.are_recvs_in_progress());
  }
}


TEST_F(DenseParallelCommTesterDouble, ClassNonBlocking)
{
  DataExchangeNonBlocking exchanger1(comm);
  DataExchangeNonBlocking exchanger2(comm);

  auto& sendLists = getSendLists();
  auto& recvLists = getReceiveLists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    exchanger1.start_nonblocking(sendLists, recvLists);
    exchanger2.start_nonblocking(sendLists2, recvLists2);
    exchanger1.post_nonblocking_receives(recvLists);
    exchanger2.post_nonblocking_receives(recvLists2);
    test_send_ranks(sendLists);
    test_recv_ranks(recvLists);

    test_send_ranks(sendLists2);
    test_recv_ranks(recvLists2);

    auto f = [&](int rank, std::vector<double>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvLists, f);
    test_results(recvLists);

    exchanger2.complete_receives(recvLists2, f);
    test_results(recvLists2);

    exchanger1.complete_sends();
    exchanger2.complete_sends();
  }
}


TEST_F(DenseParallelCommTesterInt, ClassNonBlockingBuffer)
{
  DataExchangeNonBlockingBuffer<int> exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int i=0; i < commSize; ++i)
    {
      auto& send_buf = exchanger1.get_send_buf(i);
      send_buf.assign(getSendLists()[i].begin(), getSendLists()[i].end());
    }

    exchanger1.start_nonblocking();
    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int rank, std::vector<int>&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    std::vector< std::vector<int> >recvLists = copyRecvBufs(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendLists = copySendBufs(exchanger1);
    test_send_ranks(sendLists);

    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
  }
}

TEST_F(DenseParallelCommTesterInt, ClassNonBlockingCommBuffer)
{
  DataExchangeNonBlockingCommBuffer exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int phase=0; phase < 2; ++phase) {
      for (int i=0; i < commSize; ++i) {
        auto& send_buf = exchanger1.get_send_buf(i);
        for (auto& val : getSendLists()[i]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger1.allocate_send_buffers();
      }
    }

    exchanger1.start_nonblocking();
    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int rank, stk::CommBuffer&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    recvLists = copyRecvBufs<int>(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendLists = copySendBufs<int>(exchanger1);
    test_send_ranks(sendLists);

    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
    exchanger1.clear_send_bufs();
    exchanger1.clear_recv_bufs();
  }
}


TEST_F(NeighborParallelCommTesterInt, ClassNonBlocking)
{
  DataExchangeNonBlocking exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  auto& sendLists = getSendLists();
  auto& recvLists = getReceiveLists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    exchanger1.start_nonblocking(sendLists, recvLists);

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    exchanger1.post_nonblocking_receives(recvLists);
    test_send_ranks(sendLists);
    test_recv_ranks(recvLists);

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvLists, f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    test_results(recvLists);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
  }
}


TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingBuffer)
{
  DataExchangeNonBlockingBuffer<int> exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int i=0; i < commSize; ++i)
    {
      auto& send_buf = exchanger1.get_send_buf(i);
      send_buf.assign(getSendLists()[i].begin(), getSendLists()[i].end());
    }

    exchanger1.start_nonblocking();

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());

    exchanger1.post_nonblocking_receives();

    auto f = [](int rank, std::vector<int>&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    recvLists = copyRecvBufs(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
    std::vector< std::vector<int> > sendLists = copySendBufs(exchanger1);
    test_send_ranks(sendLists);
  }
}


TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingCommBuffer)
{
  DataExchangeNonBlockingCommBuffer exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int phase=0; phase < 2; ++phase) {
      for (int i=0; i < commSize; ++i) {
        auto& send_buf = exchanger1.get_send_buf(i);
        for (auto& val : getSendLists()[i]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger1.allocate_send_buffers();
      }
    }

    exchanger1.start_nonblocking();
    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int rank, stk::CommBuffer&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    recvLists = copyRecvBufs<int>(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendLists = copySendBufs<int>(exchanger1);
    test_send_ranks(sendLists);

    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
    exchanger1.clear_send_bufs();
    exchanger1.clear_recv_bufs();
  }
}


TEST(DenseParallelCommTester, SplitComm)
{
  DenseParallelCommTesterBase<int> tester1(false);
  DenseParallelCommTesterBase<int> tester2(true);

  for (int i=0; i < 100; ++i) {
    tester1.set_offset(i);
    tester2.set_offset(i);
    stk::parallel_data_exchange_t(tester1.getSendLists(), tester1.getReceiveLists(), tester1.comm);
    if (tester2.color == 0) {
      stk::parallel_data_exchange_t(tester2.getSendLists(), tester2.getReceiveLists(), tester2.comm);
    }

    tester1.test_results(tester1.getReceiveLists());
    tester2.test_results(tester2.getReceiveLists());
  }
}

#endif

