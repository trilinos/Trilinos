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
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingProbe.hpp"
#include "stk_util/parallel/Parallel.hpp"      // for parallel_machine_rank, parallel_machine_size

#ifdef STK_HAS_MPI

#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingPrepost.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternUserDataNonBlocking.hpp"


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

using DataExchangeNonBlockingKnown = stk::DataExchangeKnownPatternNonBlocking;
template <typename T>
using DataExchangeNonBlockingBufferKnown = stk::DataExchangeKnownPatternNonBlockingBuffer<T>;
using DataExchangeNonBlockingCommBufferKnown = stk::DataExchangeKnownPatternNonBlockingCommBuffer;
using DataExchangeNonBlockingPtr = stk::DataExchangeKnownPatternUserDataNonBlocking;


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
std::vector<std::vector<T>> copy_send_bufs(const stk::ManagedBufferBase<T>& exchanger)
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
std::vector<std::vector<T>> copy_send_bufs(stk::ManagedCommBufferBase& exchanger)
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
std::vector<std::vector<T>> copy_recv_bufs(const stk::ManagedBufferBase<T>& exchanger)
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
std::vector<std::vector<T>> copy_recv_bufs(stk::ManagedCommBufferBase& exchanger)
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

    std::vector< std::vector<T> >& get_send_lists() { return sendLists; }

    std::vector< std::vector<T> >& get_receive_lists() { return recvLists; }

    virtual int get_num_sends() = 0;

    virtual int get_num_recvs() = 0;

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
      if (split_last) { /*silly statement to get past build warning*/ }

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

    void set_recv_buffer_sizes(std::vector< std::vector<T> >& rcvLists)
    {
      for (int src=0; src < commSize; ++src) {
        rcvLists[src].resize(this->get_size(src, myrank));
      }
    }

    virtual int get_num_sends() override { return commSize; }

    virtual int get_num_recvs() override { return commSize; }

    void test_results(std::vector< std::vector<T> >& rcvLists)
    {
      for (int src=0; src < commSize; ++src) {
        test_recv_vals(rcvLists[src], src);
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

    void test_send_ranks(std::vector< std::vector<T> >& sndLists)
    {
      test_ranks_inner(sndLists);
    }


    void test_recv_ranks(std::vector<std::vector<T>>& recvListsArg)
    {
      test_ranks_inner(recvListsArg);
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

template <typename ExchangerTypeT>
class DenseParallelCommTesterIntT : public ::testing::Test,
                                    public DenseParallelCommTesterBase<int>
{
  using ExchangerType = ExchangerTypeT;
};

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

    void set_recv_buffer_sizes(std::vector< std::vector<T> >& rcvLists)
    {
      int src1 = (myrank - 1 + commSize) % commSize;
      int src2 = (myrank - 2 + commSize) % commSize;
      rcvLists[src1].resize(this->get_size(src1, myrank));
      rcvLists[src2].resize(this->get_size(src2, myrank));
    }

    virtual int get_num_sends() override { return std::min(2, commSize); }

    virtual int get_num_recvs() override { return std::min(2, commSize); }
    
    std::vector<int> get_dest_ranks()
    {
      int dest1 = (myrank + 1) % commSize;
      int dest2 = (myrank + 2) % commSize;
      if (dest1 == dest2)
      {
        return {dest1};
      } else
      {
        return {dest1, dest2};
      }
    }
  
  std::vector<int> get_src_ranks()
  {
    int src1 = (myrank - 1 + commSize) % commSize;
    int src2 = (myrank - 2 + commSize) % commSize;
    
    if (src1 == src2)
    {
      return {src1};
    } else
    {
      return {src1, src2};
    }
  }

    void test_results(std::vector< std::vector<T> >& rcvLists)
    {
      for (int src=0; src < commSize; ++src) {
        test_recv_vals(rcvLists[src], src);
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

    void test_send_ranks(std::vector<std::vector<T>>& sndLists)
    {

      std::vector<int> sendRanks = get_ranks(sndLists);

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

    void test_recv_ranks(std::vector<std::vector<T>>& rcvLists)
    {
      auto recvRanks = get_ranks(rcvLists);

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

TEST(UnknownPatternExchanger, StringToEnum)
{
  EXPECT_EQ(stk::stringToUnknownPatternExchangerEnum("Default"), stk::UnknownPatternExchanger::Default);
  EXPECT_EQ(stk::stringToUnknownPatternExchangerEnum("Probe"),   stk::UnknownPatternExchanger::Probe);
  EXPECT_EQ(stk::stringToUnknownPatternExchangerEnum("Prepost"), stk::UnknownPatternExchanger::Prepost);
  EXPECT_ANY_THROW(stk::stringToUnknownPatternExchangerEnum("foo"));
}

//-----------------------------------------------------------------------------
// Test parallel_data_exchange_t

TEST_F(DenseParallelCommTesterInt, FuncBlocking)
{
  for (int i=0; i < 100; ++i) {
    stk::parallel_data_exchange_t(get_send_lists(), get_receive_lists(), comm);
    test_results(get_receive_lists());
  }
}


TEST_F(DenseParallelCommTesterDouble, FuncBlocking)
{
  for (int i=0; i < 100; ++i) {
    stk::parallel_data_exchange_t(get_send_lists(), get_receive_lists(), comm);
    test_results(get_receive_lists());
  }
}


TEST_F(NeighborParallelCommTesterInt, FuncBlocking)
{
  stk::parallel_data_exchange_t(get_send_lists(), get_receive_lists(), comm);
  test_results(get_receive_lists());
}


TEST_F(NeighborParallelCommTesterDouble, FuncBlocking)
{
  stk::parallel_data_exchange_t(get_send_lists(), get_receive_lists(), comm);
  test_results(get_receive_lists());
}

//-----------------------------------------------------------------------------
// test parallel_data_exchange_nonsym_known_sizes_t
TEST(NonsymKnownsizes, nominal_correct)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs != 4) { GTEST_SKIP(); }

  std::vector<int> sendOffsets, recvOffsets;
  std::vector<double> sendData, recvData, expectedRecvData;

  const int myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);

  if (myProc == 2) {
    sendOffsets = { 0, 0, 0, 0, 0 };
    recvOffsets = { 0, 2, 4, 4, 6 }; //recving from procs 0, 1, 3
    expectedRecvData = { 42.0, 99.0, 42.0, 99.0, 42.0, 99.0 };
    recvData.resize(6);
  }
  else {
    sendOffsets = { 0, 0, 0, 2, 2 }; //sending to proc 2
    sendData = { 42.0, 99.0 };
    recvOffsets = { 0, 0, 0, 0, 0 };
  }

  const bool checkInput = true;
  EXPECT_NO_THROW(
  stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendData.data(),
                                                   recvOffsets.data(), recvData.data(),
                                                   MPI_COMM_WORLD, checkInput));

  EXPECT_EQ(recvData.size(), expectedRecvData.size());
  if (myProc == 2) {
    EXPECT_EQ(6u, recvData.size());
  }

  for(unsigned i=0; i<recvData.size(); ++i) {
    EXPECT_NEAR(recvData[i], expectedRecvData[i], 1.e-8);
  }
}

TEST(NonsymKnownsizes, inconsistent_send_input)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs != 4) { GTEST_SKIP(); }

  std::vector<int> sendOffsets, recvOffsets;
  std::vector<double> sendData, recvData, expectedRecvData;

  const int myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);

  if (myProc == 2) {
    sendOffsets = { 0, 0, 0, 0, 0 };
    recvOffsets = { 0, 2, 4, 4, 6 }; //recving from procs 0, 1, 3
    expectedRecvData = { 42.0, 99.0, 42.0, 99.0, 42.0, 99.0 };
    recvData.resize(6);
  }
  else if (myProc == 0) {
    sendOffsets = { 0, 0, 1, 3, 3 }; //sending to procs 1 and 2
    sendData = { 42.0, 99.0 };
    recvOffsets = { 0, 0, 0, 0, 0 };
  }
  else {
    sendOffsets = { 0, 0, 0, 2, 2 };
    sendData = { 42.0, 99.0 };
    recvOffsets = { 0, 0, 0, 0, 0 };
  }

  const bool checkInput = true;
  EXPECT_ANY_THROW(
  stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendData.data(),
                                                   recvOffsets.data(), recvData.data(),
                                                   MPI_COMM_WORLD, checkInput));

}

TEST(NonsymKnownsizes, inconsistent_recv_input)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs != 4) { GTEST_SKIP(); }

  std::vector<int> sendOffsets, recvOffsets;
  std::vector<double> sendData, recvData, expectedRecvData;

  const int myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);

  if (myProc == 2) {
    sendOffsets = { 0, 0, 0, 0, 0 };
    recvOffsets = { 0, 0, 2, 2, 4 }; // only recving from procs 1 and 3
    expectedRecvData = { 42.0, 99.0, 42.0, 99.0, 42.0, 99.0 };
    recvData.resize(6);
  }
  else {
    sendOffsets = { 0, 0, 0, 2, 2 }; //sending to proc 2
    sendData = { 42.0, 99.0 };
    recvOffsets = { 0, 0, 0, 0, 0 };
  }

  const bool checkInput = true;
  EXPECT_ANY_THROW(
  stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendData.data(),
                                                   recvOffsets.data(), recvData.data(),
                                                   MPI_COMM_WORLD, checkInput));

}

TEST(NonsymKnownsizes, inconsistent_msg_sizing)
{
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numProcs != 4) { GTEST_SKIP(); }

  std::vector<int> sendOffsets, recvOffsets;
  std::vector<double> sendData, recvData, expectedRecvData;

  const int myProc = stk::parallel_machine_rank(MPI_COMM_WORLD);

  if (myProc == 2) {
    sendOffsets = { 0, 0, 0, 0, 0 };
    recvOffsets = { 0, 2, 4, 4, 6 }; //recving 2 values each from procs 0, 1, 3
    expectedRecvData = { 42.0, 99.0, 42.0, 99.0, 42.0, 99.0 };
    recvData.resize(6);
  }
  else {
    sendOffsets = { 0, 0, 0, 1, 1 }; //sending 1 value to proc 2
    sendData = { 42.0, 99.0 };
    recvOffsets = { 0, 0, 0, 0, 0 };
  }

  const bool checkInput = true;
  EXPECT_ANY_THROW(
  stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), sendData.data(),
                                                   recvOffsets.data(), recvData.data(),
                                                   MPI_COMM_WORLD, checkInput));
}

//-----------------------------------------------------------------------------
// test ParallelDataExchange class

TEST_F(DenseParallelCommTesterInt, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    if (i % 2 == 0) {
      exchanger.execute(get_send_lists(), get_receive_lists());
    } else {
      exchanger.execute(get_send_lists(), get_receive_lists(), get_num_recvs());
    }
    test_results(get_receive_lists());
    test_send_ranks(get_send_lists());
    test_recv_ranks(get_receive_lists());
  }

}

TEST_F(DenseParallelCommTesterInt, ClassBlockingBuffer)
{
  DataExchangeBlockingBuffer<int> exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int j=0; j < commSize; ++j)
    {
      auto& buf_in = get_send_lists()[j];
      auto& buf_out = exchanger.get_send_buf(j);
      buf_out.assign(buf_in.begin(), buf_in.end());
    }

    if (i % 2 == 0) {
      exchanger.execute();
    } else {
      exchanger.execute(get_num_recvs());
    }
    
    
    std::vector<std::vector<int>> send_bufs = copy_send_bufs(exchanger);
    std::vector<std::vector<int>> recv_bufs = copy_recv_bufs(exchanger);

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
      for (int j=0; j < commSize; ++j) {
        auto& send_buf = exchanger.get_send_buf(j);
        for (auto& val : get_send_lists()[j]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger.allocate_send_buffers();
      }
    }

    if (i % 2 == 0) {
      exchanger.execute();
    } else {
      exchanger.execute(get_num_recvs());
    }
    
    std::vector<std::vector<int>> send_bufs = copy_send_bufs<int>(exchanger);
    std::vector<std::vector<int>> recv_bufs = copy_recv_bufs<int>(exchanger);

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
    if (i % 2 == 0) {
      exchanger.execute(get_send_lists(), get_receive_lists());
    } else {
      exchanger.execute(get_send_lists(), get_receive_lists(), get_num_recvs());
    }

    test_results(get_receive_lists());
    test_send_ranks(get_send_lists());
    test_recv_ranks(get_receive_lists());
  }
}


TEST_F(NeighborParallelCommTesterInt, ClassBlocking)
{
  DataExchangeBlocking exchanger(comm);
  for (int i=0; i < 100; ++i) {
    set_offset(i);
    if (i % 2 == 0) {
      exchanger.execute(get_send_lists(), get_receive_lists());
    } else {
      exchanger.execute(get_send_lists(), get_receive_lists(), get_num_recvs());
    }
    test_results(get_receive_lists());
    test_send_ranks(get_send_lists());
    test_recv_ranks(get_receive_lists());
  }

}

TEST_F(NeighborParallelCommTesterInt, ClassBlockingBuffer)
{
  DataExchangeBlockingBuffer<int> exchanger(comm);

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int j=0; j < commSize; ++j)
    {
      auto& buf_in = get_send_lists()[j];
      auto& buf_out = exchanger.get_send_buf(j);
      buf_out.assign(buf_in.begin(), buf_in.end());
    }

    if (i % 2 == 0) {
      exchanger.execute();
    } else {
      exchanger.execute(get_num_recvs());
    }
    
    std::vector<std::vector<int>> send_bufs = copy_send_bufs(exchanger);
    std::vector<std::vector<int>> recv_bufs = copy_recv_bufs(exchanger);

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
      for (int j=0; j < commSize; ++j) {
        auto& send_buf = exchanger.get_send_buf(j);
        for (auto& val : get_send_lists()[j]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger.allocate_send_buffers();
      }
    }

    if (i % 2 == 0) {
      exchanger.execute();
    } else {
      exchanger.execute(get_num_recvs());
    }
    
    std::vector<std::vector<int>> send_bufs = copy_send_bufs<int>(exchanger);
    std::vector<std::vector<int>> recv_bufs = copy_recv_bufs<int>(exchanger);

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

    if (i % 2 == 0) {
      exchanger.execute(get_send_lists(), get_receive_lists());
    } else {
      exchanger.execute(get_send_lists(), get_receive_lists(), get_num_recvs());

    }
    test_results(get_receive_lists());
    test_send_ranks(get_send_lists());
    test_recv_ranks(get_receive_lists());
  }
}

namespace {
  using ExchangerTypes = ::testing::Types<DataExchangeNonBlocking,
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

TYPED_TEST_SUITE(DenseParallelCommTesterIntT, ExchangerTypes, NameGenerator);

TYPED_TEST(DenseParallelCommTesterIntT, ClassNonBlocking)
{
  using ExchType = TypeParam;
  using ElemType = int;
  
  ExchType exchanger1(this->comm);
  ExchType exchanger2(this->comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());
  EXPECT_FALSE(exchanger2.are_recvs_in_progress());
  EXPECT_FALSE(exchanger2.are_sends_in_progress());

  auto& sendListsTest = this->get_send_lists();
  auto& recvListsTest = this->get_receive_lists();

  for (int i=0; i < 100; ++i) {
    this->set_offset(i);
    
    std::vector< std::vector<ElemType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElemType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvListsTest);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_TRUE(exchanger2.are_recvs_in_progress());
    EXPECT_TRUE(exchanger2.are_sends_in_progress());
    exchanger1.post_nonblocking_receives(recvListsTest);
    exchanger2.post_nonblocking_receives(recvLists2);
    
    this->test_send_ranks(sendListsTest);
    this->test_recv_ranks(recvListsTest);

    this->test_send_ranks(sendLists2);
    this->test_recv_ranks(recvLists2);

    auto f = [&](int rank, std::vector<int>& buf) { this->test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvListsTest, f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    this->test_results(recvListsTest);

    exchanger2.complete_receives(recvLists2, f);
    EXPECT_FALSE(exchanger2.are_recvs_in_progress());
    EXPECT_TRUE(exchanger2.are_sends_in_progress());
    this->test_results(recvLists2);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());

    exchanger2.complete_sends();
    EXPECT_FALSE(exchanger2.are_sends_in_progress());
    EXPECT_FALSE(exchanger2.are_recvs_in_progress());
  }
}

/*
TEST_F(DenseParallelCommTesterInt, ClassNonBlockingPrepost)
{
  stk::DataExchangeUnknownPatternNonBlockingPrepost exchanger1(comm);
  stk::DataExchangeUnknownPatternNonBlockingPrepost exchanger2(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());
  EXPECT_FALSE(exchanger2.are_recvs_in_progress());
  EXPECT_FALSE(exchanger2.are_sends_in_progress());

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElementType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvLists);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_TRUE(exchanger2.are_recvs_in_progress());
    EXPECT_TRUE(exchanger2.are_sends_in_progress());
    exchanger1.post_nonblocking_receives(recvListsTest);
    exchanger2.post_nonblocking_receives(recvLists2);
    test_send_ranks(sendListsTest);
    test_recv_ranks(recvListsTest);

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
*/

TEST_F(DenseParallelCommTesterInt, ClassNonBlockingKnownPattern)
{
  DataExchangeNonBlockingKnown exchanger1(comm);
  DataExchangeNonBlockingKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElementType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvLists);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvListsTest, f);
    test_results(recvListsTest);

    exchanger2.complete_receives(recvLists2, f);
    test_results(recvLists2);

    exchanger1.complete_sends();
    exchanger2.complete_sends();
  }
}

TEST_F(DenseParallelCommTesterInt, ClassNonBlockingUserDataKnownPattern)
{
  DataExchangeNonBlockingPtr exchanger(comm);
  
  for (int i=0; i < 1; ++i)
  {
    set_offset(i);
    std::vector<int> allSendBufs, allRecvBufs;
    std::vector<stk::PointerAndSize> sendBufPointers, recvBufPointers;
    std::vector<int> sendRanks, recvRanks;
    
    size_t totalSendSize = 0;
    size_t totalRecvSize = 0;
    for (int rank=0; rank < commSize; ++rank)
    {
      totalSendSize += get_size(myrank, rank);
      totalRecvSize += get_size(rank, myrank);      
    }
    
    allSendBufs.reserve(totalSendSize);
    allRecvBufs.resize(totalRecvSize);
    
    size_t valuesSent = 0;
    size_t valuesReceived = 0;
    for (int rank=0; rank < commSize; ++rank)
    {
      int sendCount = get_size(myrank, rank);
      int recvCount = get_size(rank, myrank);
      
      for (int j=0; j < sendCount; ++j)
      {
        allSendBufs.push_back(get_value(myrank, rank, j));
      }
      auto sendBufBegin = reinterpret_cast<unsigned char*>(allSendBufs.data() + valuesSent);
      sendBufPointers.emplace_back(sendBufBegin, sendCount * sizeof(int));
      sendRanks.push_back(rank);
      valuesSent += sendCount;
      
      auto recvBufBegin = reinterpret_cast<unsigned char*>(allRecvBufs.data() + valuesReceived);
      recvRanks.push_back(rank);
      recvBufPointers.emplace_back(recvBufBegin, recvCount * sizeof (int));
      valuesReceived += recvCount;
    }
    
    exchanger.start_nonblocking(sendBufPointers, sendRanks,
                                recvBufPointers, recvRanks);
    
    std::vector<bool> senderRankSeen(commSize, false);
    auto unpacker = [&](int senderRank, const stk::PointerAndSize& buf)
    {
      size_t offset = 0;
      for (int rank=0; rank < senderRank; ++rank)
      {
        offset += get_size(rank, myrank);
      }
      
      int recvCount  = get_size(senderRank, myrank);
      const int* ptr = reinterpret_cast<const int*>(buf.ptr);      
      EXPECT_EQ(buf.size, recvCount * sizeof(int));
      EXPECT_EQ(ptr, allRecvBufs.data() + offset);
      
      for (int k=0; k < recvCount; ++k)
      {
        int expectedValue = get_value(senderRank, myrank, k);
        EXPECT_EQ(ptr[k],                  expectedValue);
        EXPECT_EQ(allRecvBufs[k + offset], expectedValue);
      }
      senderRankSeen[senderRank] = true;    
    };
    
    exchanger.complete_receives(recvBufPointers, recvRanks, unpacker);
    for (bool val : senderRankSeen)
    {
      EXPECT_TRUE(val);
    }
    
    exchanger.complete_sends();    
  }
}

TEST_F(DenseParallelCommTesterInt, ClassNonBlockingBufferKnownPattern)
{
  DataExchangeNonBlockingBufferKnown<int> exchanger1(comm);
  DataExchangeNonBlockingBufferKnown<int> exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int rank=0; rank < commSize; ++rank) {
      exchanger1.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger2.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger1.get_recv_buf(rank).resize(recvListsTest[rank].size());
      exchanger2.get_recv_buf(rank).resize(recvListsTest[rank].size());
    }

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(f);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_recv_bufs();
  }
}

TEST_F(DenseParallelCommTesterInt, ClassNonBlockingCommBufferKnownPattern)
{
  DataExchangeNonBlockingCommBufferKnown exchanger1(comm);
  DataExchangeNonBlockingCommBufferKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int phase=0; phase < 2; ++phase)
    {
      for (int rank=0; rank < commSize; ++rank) {
        for (auto& val : sendListsTest[rank]) {
          exchanger1.get_send_buf(rank).pack(val);
          exchanger2.get_send_buf(rank).pack(val);
        }

        exchanger1.set_recv_buffer_size(rank, sizeof(int)*recvListsTest[rank].size());
        exchanger2.set_recv_buffer_size(rank, sizeof(int)*recvListsTest[rank].size());
      }

      if (phase == 0)
      {
        exchanger1.allocate_send_buffers();
        exchanger2.allocate_send_buffers();
      }
    }

    exchanger1.allocate_recv_buffers();
    exchanger2.allocate_recv_buffers();

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int /*rank*/, stk::CommBuffer& /*buf*/) {};

    exchanger1.complete_receives(f);
    auto recvListsCopy = copy_recv_bufs<int>(exchanger1);
    test_results(recvListsCopy);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    recvListsCopy = copy_recv_bufs<int>(exchanger2);
    test_results(recvListsCopy);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_send_bufs();
  }
}


TEST_F(DenseParallelCommTesterDouble, ClassNonBlocking)
{
  DataExchangeNonBlocking exchanger1(comm);
  DataExchangeNonBlocking exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElementType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvListsTest);
    exchanger2.start_nonblocking(sendLists2, recvLists2);
    exchanger1.post_nonblocking_receives(recvListsTest);
    exchanger2.post_nonblocking_receives(recvLists2);
    test_send_ranks(sendListsTest);
    test_recv_ranks(recvListsTest);

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

    for (int j=0; j < commSize; ++j)
    {
      auto& send_buf = exchanger1.get_send_buf(j);
      send_buf.assign(get_send_lists()[j].begin(), get_send_lists()[j].end());
    }

    if (i % 2 == 0) {
      exchanger1.start_nonblocking();
    } else {
      exchanger1.start_nonblocking(get_num_recvs());
    }

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int , std::vector<int>&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    std::vector< std::vector<int> >recvListsTest = copy_recv_bufs(exchanger1);
    test_recv_ranks(recvListsTest);
    test_results(recvListsTest);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendListsTest = copy_send_bufs(exchanger1);
    test_send_ranks(sendListsTest);

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
      for (int j=0; j < commSize; ++j) {
        auto& send_buf = exchanger1.get_send_buf(j);
        for (auto& val : get_send_lists()[j]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger1.allocate_send_buffers();
      }
    }

    if (i % 2 == 0) {
      exchanger1.start_nonblocking();
    } else { 
      exchanger1.start_nonblocking(get_num_recvs());
    }

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int , stk::CommBuffer&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    recvLists = copy_recv_bufs<int>(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendListsTest = copy_send_bufs<int>(exchanger1);
    test_send_ranks(sendListsTest);

    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
    exchanger1.clear_send_bufs();
    exchanger1.clear_recv_bufs();
  }
}

TEST_F(DenseParallelCommTesterDouble, ClassNonBlockingKnownPattern)
{
  DataExchangeNonBlockingKnown exchanger1(comm);
  DataExchangeNonBlockingKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElementType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvListsTest);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    auto f = [&](int rank, std::vector<double>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvListsTest, f);
    test_results(recvListsTest);

    exchanger2.complete_receives(recvLists2, f);
    test_results(recvLists2);

    exchanger1.complete_sends();
    exchanger2.complete_sends();
  }
}

TEST_F(DenseParallelCommTesterDouble, ClassNonBlockingUserDataKnownPattern)
{
  DataExchangeNonBlockingPtr exchanger(comm);
  
  for (int i=0; i < 1; ++i)
  {
    set_offset(i);
    std::vector<double> allSendBufs, allRecvBufs;
    std::vector<stk::PointerAndSize> sendBufPointers, recvBufPointers;
    std::vector<int> sendRanks, recvRanks;
    
    size_t totalSendSize = 0;
    size_t totalRecvSize = 0;
    for (int rank=0; rank < commSize; ++rank)
    {
      totalSendSize += get_size(myrank, rank);
      totalRecvSize += get_size(rank, myrank);      
    }
    
    allSendBufs.reserve(totalSendSize);
    allRecvBufs.resize(totalRecvSize);
    
    size_t valuesSent = 0;
    size_t valuesReceived = 0;
    for (int rank=0; rank < commSize; ++rank)
    {
      int sendCount = get_size(myrank, rank);
      int recvCount = get_size(rank, myrank);
      
      for (int j=0; j < sendCount; ++j)
      {
        allSendBufs.push_back(get_value(myrank, rank, j));
      }
      auto sendBufBegin = reinterpret_cast<unsigned char*>(allSendBufs.data() + valuesSent);
      sendBufPointers.emplace_back(sendBufBegin, sendCount * sizeof(double));
      sendRanks.push_back(rank);
      valuesSent += sendCount;
      
      auto recvBufBegin = reinterpret_cast<unsigned char*>(allRecvBufs.data() + valuesReceived);
      recvRanks.push_back(rank);
      recvBufPointers.emplace_back(recvBufBegin, recvCount * sizeof(double));
      valuesReceived += recvCount;
    }
    
    exchanger.start_nonblocking(sendBufPointers, sendRanks,
                                recvBufPointers, recvRanks);
    
    std::vector<bool> senderRankSeen(commSize, false);
    auto unpacker = [&](int senderRank, const stk::PointerAndSize& buf)
    {
      size_t offset = 0;
      for (int rank=0; rank < senderRank; ++rank)
      {
        offset += get_size(rank, myrank);
      }
      
      int recvCount  = get_size(senderRank, myrank);
      const double* ptr = reinterpret_cast<const double*>(buf.ptr);      
      EXPECT_EQ(buf.size, recvCount * sizeof(double));
      EXPECT_EQ(ptr, allRecvBufs.data() + offset);
      
      for (int k=0; k < recvCount; ++k)
      {
        int expectedValue = get_value(senderRank, myrank, k);
        EXPECT_EQ(ptr[k],                  expectedValue);
        EXPECT_EQ(allRecvBufs[k + offset], expectedValue);
      }
      senderRankSeen[senderRank] = true;    
    };
    
    exchanger.complete_receives(recvBufPointers, recvRanks, unpacker);
    for (bool val : senderRankSeen)
    {
      EXPECT_TRUE(val);
    }
    
    exchanger.complete_sends();    
  }
}


TEST_F(DenseParallelCommTesterDouble, ClassNonBlockingBufferKnownPattern)
{
  DataExchangeNonBlockingBufferKnown<double> exchanger1(comm);
  DataExchangeNonBlockingBufferKnown<double> exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int rank=0; rank < commSize; ++rank) {
      exchanger1.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger2.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger1.get_recv_buf(rank).resize(recvListsTest[rank].size());
      exchanger2.get_recv_buf(rank).resize(recvListsTest[rank].size());
    }

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int rank, std::vector<double>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(f);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_recv_bufs();
  }
}

TEST_F(DenseParallelCommTesterDouble, ClassNonBlockingCommBufferKnownPattern)
{
  DataExchangeNonBlockingCommBufferKnown exchanger1(comm);
  DataExchangeNonBlockingCommBufferKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int phase=0; phase < 2; ++phase)
    {
      for (int rank=0; rank < commSize; ++rank) {
        for (auto& val : sendListsTest[rank]) {
          exchanger1.get_send_buf(rank).pack(val);
          exchanger2.get_send_buf(rank).pack(val);
        }

        exchanger1.set_recv_buffer_size(rank, sizeof(double)*recvListsTest[rank].size());
        exchanger2.set_recv_buffer_size(rank, sizeof(double)*recvListsTest[rank].size());
      }

      if (phase == 0)
      {
        exchanger1.allocate_send_buffers();
        exchanger2.allocate_send_buffers();
      }
    }

    exchanger1.allocate_recv_buffers();
    exchanger2.allocate_recv_buffers();

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int , stk::CommBuffer& ) {};

    exchanger1.complete_receives(f);
    auto recvListsCopy = copy_recv_bufs<double>(exchanger1);
    test_results(recvListsCopy);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    recvListsCopy = copy_recv_bufs<double>(exchanger2);
    test_results(recvListsCopy);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_send_bufs();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlocking)
{
  DataExchangeNonBlocking exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    if (i % 2 == 0) {
      exchanger1.start_nonblocking(sendListsTest, recvListsTest);
    } else {
      exchanger1.start_nonblocking(sendListsTest, recvListsTest, get_num_recvs());
    }

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    exchanger1.post_nonblocking_receives(recvListsTest);
    test_send_ranks(sendListsTest);
    test_recv_ranks(recvListsTest);

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvListsTest, f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    test_results(recvListsTest);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingKnownPattern)
{
  DataExchangeNonBlockingKnown exchanger1(comm);
  DataExchangeNonBlockingKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvLists);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    std::vector< std::vector<ElementType> > sendLists2 = sendListsTest;
    std::vector< std::vector<ElementType> > recvLists2 = recvListsTest;

    exchanger1.start_nonblocking(sendListsTest, recvListsTest);
    exchanger2.start_nonblocking(sendLists2, recvLists2);

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(recvListsTest, f);
    test_results(recvListsTest);

    exchanger2.complete_receives(recvLists2, f);
    test_results(recvLists2);

    exchanger1.complete_sends();
    exchanger2.complete_sends();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingUserDataKnownPattern)
{
  DataExchangeNonBlockingPtr exchanger(comm);
  
  for (int i=0; i < 1; ++i)
  {
    set_offset(i);
    std::vector<int> allSendBufs, allRecvBufs;
    std::vector<stk::PointerAndSize> sendBufPointers, recvBufPointers;
    std::vector<int> sendRanks = get_dest_ranks();
    std::vector<int> recvRanks = get_src_ranks();
    
    size_t totalSendSize = 0;
    size_t totalRecvSize = 0;
    for (int rank : sendRanks)
    {
      totalSendSize += get_size(myrank, rank);
    }
    
    for (int rank : recvRanks)
    {
      totalRecvSize += get_size(rank, myrank);
    }
    
    allSendBufs.reserve(totalSendSize);
    allRecvBufs.resize(totalRecvSize);
    
    size_t valuesSent = 0;
    size_t valuesReceived = 0;
    for (int rank : sendRanks)
    {
      int sendCount = get_size(myrank, rank);
      for (int j=0; j < sendCount; ++j)
      {
        allSendBufs.push_back(get_value(myrank, rank, j));
      }
      auto sendBufBegin = reinterpret_cast<unsigned char*>(allSendBufs.data() + valuesSent);
      sendBufPointers.emplace_back(sendBufBegin, sendCount * sizeof(int));
      valuesSent += sendCount;      
    }
    
    for (int rank : recvRanks)
    {
      int recvCount = get_size(rank, myrank);      
      auto recvBufBegin = reinterpret_cast<unsigned char*>(allRecvBufs.data() + valuesReceived);
      recvBufPointers.emplace_back(recvBufBegin, recvCount * sizeof (int));
      valuesReceived += recvCount;
    }
    
    exchanger.start_nonblocking(sendBufPointers, sendRanks,
                                recvBufPointers, recvRanks);
    
    auto unpacker = [&](int senderRank, const stk::PointerAndSize& buf)
    {
      size_t offset = 0;
      auto senderRankIt = std::find(recvRanks.begin(), recvRanks.end(), senderRank);
      EXPECT_NE(senderRankIt, recvRanks.end());
      for (auto it = recvRanks.begin(); it != senderRankIt; ++it)
      {
        offset += get_size(*it, myrank);
      }
      
      int recvCount  = get_size(senderRank, myrank);
      const int* ptr = reinterpret_cast<const int*>(buf.ptr);      
      EXPECT_EQ(buf.size, recvCount * sizeof(int));
      EXPECT_EQ(ptr, allRecvBufs.data() + offset);
      
      for (int k=0; k < recvCount; ++k)
      {
        int expectedValue = get_value(senderRank, myrank, k);
        EXPECT_EQ(ptr[k],                  expectedValue);
        EXPECT_EQ(allRecvBufs[k + offset], expectedValue);
      }
    };
    
    exchanger.complete_receives(recvBufPointers, recvRanks, unpacker);
    exchanger.complete_sends();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingBufferKnownPattern)
{
  DataExchangeNonBlockingBufferKnown<int> exchanger1(comm);
  DataExchangeNonBlockingBufferKnown<int> exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvListsTest);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int rank=0; rank < commSize; ++rank) {
      exchanger1.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger2.get_send_buf(rank).assign(sendListsTest[rank].begin(), sendListsTest[rank].end());
      exchanger1.get_recv_buf(rank).resize(recvListsTest[rank].size());
      exchanger2.get_recv_buf(rank).resize(recvListsTest[rank].size());
    }

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int rank, std::vector<int>& buf) { test_recv_vals(buf, rank); };

    exchanger1.complete_receives(f);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_recv_bufs();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingCommBufferKnownPattern)
{
  DataExchangeNonBlockingCommBufferKnown exchanger1(comm);
  DataExchangeNonBlockingCommBufferKnown exchanger2(comm);

  auto& sendListsTest = get_send_lists();
  auto& recvListsTest = get_receive_lists();
  set_recv_buffer_sizes(recvListsTest);

  for (int i=0; i < 100; ++i) {
    set_offset(i);
    //std::vector< std::vector<ElementType> > sendLists2 = sendLists;
    //std::vector< std::vector<ElementType> > recvLists2 = recvLists;

    for (int phase=0; phase < 2; ++phase)
    {
      for (int rank=0; rank < commSize; ++rank) {
        for (auto& val : sendListsTest[rank]) {
          exchanger1.get_send_buf(rank).pack(val);
          exchanger2.get_send_buf(rank).pack(val);
        }

        exchanger1.set_recv_buffer_size(rank, sizeof(int)*recvListsTest[rank].size());
        exchanger2.set_recv_buffer_size(rank, sizeof(int)*recvListsTest[rank].size());
      }

      if (phase == 0)
      {
        exchanger1.allocate_send_buffers();
        exchanger2.allocate_send_buffers();
      }
    }

    exchanger1.allocate_recv_buffers();
    exchanger2.allocate_recv_buffers();

    exchanger1.start_nonblocking();
    exchanger2.start_nonblocking();

    auto f = [&](int , stk::CommBuffer& ) {};

    exchanger1.complete_receives(f);
    auto recvListsCopy = copy_recv_bufs<int>(exchanger1);
    test_results(recvListsCopy);
    exchanger1.clear_recv_bufs();

    exchanger2.complete_receives(f);
    recvListsCopy = copy_recv_bufs<int>(exchanger2);
    test_results(recvListsCopy);
    exchanger2.clear_recv_bufs();

    exchanger1.complete_sends();
    exchanger2.complete_sends();

    exchanger1.clear_send_bufs();
    exchanger2.clear_send_bufs();
  }
}

TEST_F(NeighborParallelCommTesterInt, ClassNonBlockingBuffer)
{
  DataExchangeNonBlockingBuffer<int> exchanger1(comm);

  EXPECT_FALSE(exchanger1.are_recvs_in_progress());
  EXPECT_FALSE(exchanger1.are_sends_in_progress());

  for (int i=0; i < 100; ++i) {
    set_offset(i);

    for (int j=0; j < commSize; ++j)
    {
      auto& send_buf = exchanger1.get_send_buf(j);
      send_buf.assign(get_send_lists()[j].begin(), get_send_lists()[j].end());
    }

    if (i % 2 == 0) {
      exchanger1.start_nonblocking();
    } else {
      exchanger1.start_nonblocking(get_num_recvs());
    }

    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());

    exchanger1.post_nonblocking_receives();

    auto f = [](int , std::vector<int>&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    recvLists = copy_recv_bufs(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_FALSE(exchanger1.are_sends_in_progress());
    std::vector< std::vector<int> > sendListsTest = copy_send_bufs(exchanger1);
    test_send_ranks(sendListsTest);
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
      for (int j=0; j < commSize; ++j) {
        auto& send_buf = exchanger1.get_send_buf(j);
        for (auto& val : get_send_lists()[j]) {
          send_buf.pack(val);
        }
      }
      if (phase == 0) {
        exchanger1.allocate_send_buffers();
      }
    }

    if (i % 2 == 0) {
      exchanger1.start_nonblocking();
    } else {
      exchanger1.start_nonblocking(get_num_recvs());
    }
    EXPECT_TRUE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    EXPECT_ANY_THROW(exchanger1.get_recv_buf(0));

    exchanger1.post_nonblocking_receives();

    auto f = [](int , stk::CommBuffer&) {};

    exchanger1.complete_receives(f);
    EXPECT_FALSE(exchanger1.are_recvs_in_progress());
    EXPECT_TRUE(exchanger1.are_sends_in_progress());
    EXPECT_ANY_THROW(exchanger1.get_send_buf(0));
    recvLists = copy_recv_bufs<int>(exchanger1);
    test_recv_ranks(recvLists);
    test_results(recvLists);

    exchanger1.complete_sends();
    std::vector< std::vector<int> > sendListsTest = copy_send_bufs<int>(exchanger1);
    test_send_ranks(sendListsTest);

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
    stk::parallel_data_exchange_t(tester1.get_send_lists(), tester1.get_receive_lists(), tester1.comm);
    if (tester2.color == 0) {
      stk::parallel_data_exchange_t(tester2.get_send_lists(), tester2.get_receive_lists(), tester2.comm);
    }

    tester1.test_results(tester1.get_receive_lists());
    tester2.test_results(tester2.get_receive_lists());
  }
}

namespace {
void send_to_self(stk::DataExchangeKnownPatternNonBlockingCommBuffer& exchanger, int size)
{
  using T = int;
  int myrank = stk::parallel_machine_rank(MPI_COMM_WORLD);

  for (int phase=0; phase < 2; ++phase)
  {
    for (int i=0; i < size; ++i)
    {
      exchanger.get_send_buf(myrank).pack<T>(i);
    }

    if (phase == 0)
      exchanger.allocate_send_buffers();
  }
  exchanger.set_recv_buffer_size(myrank, size * sizeof(T));
  exchanger.allocate_recv_buffers();

  exchanger.start_nonblocking();

  auto unpacker = [&](int rank, stk::CommBuffer& buf)
  {
    EXPECT_EQ(rank, myrank);
    EXPECT_EQ(buf.remaining(), ptrdiff_t(size * sizeof(T)));
    for (int i=0; i < size; ++i)
    {
      T val;
      buf.unpack(val);
      EXPECT_EQ(val, i);
    }
  };

  exchanger.complete_receives(unpacker);
  exchanger.complete_sends();
}
}

TEST(ManagedCommBufferBase, IncreasingSize)
{
  stk::DataExchangeKnownPatternNonBlockingCommBuffer exchanger(MPI_COMM_WORLD);
  send_to_self(exchanger, 16);

  exchanger.clear_send_bufs();
  exchanger.clear_recv_bufs();

  send_to_self(exchanger, 32);
}

#endif

