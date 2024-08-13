#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_search/DeviceMPIUtils.hpp"
#include "Kokkos_Sort.hpp"

namespace {
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using DeviceBuffers = stk::search::impl::DeviceMPIBuffers<int, ExecutionSpace>;
using DeviceBufferAppender = stk::search::impl::DeviceMPIBufferAppender<int, ExecutionSpace>;
using DeviceDataExchanger = stk::search::impl::DeviceDataExchangeUnknownPattern<int, ExecutionSpace>;

template <typename ViewType>
void set_on_device(ViewType viewDevice, const std::vector<typename ViewType::value_type>& vals)
{
  STK_ThrowRequire(viewDevice.extent(0) == vals.size());
  auto viewHost = Kokkos::create_mirror_view(viewDevice);
  for (size_t i=0; i < vals.size(); ++i)
  {
    viewHost(i) = vals[i];
  }

  Kokkos::deep_copy(viewDevice, viewHost);
}

}

TEST(DeviceMPIBuffers, ConstructorSizeOnly)
{
  DeviceBuffers buffers(3);
  EXPECT_EQ(buffers.bufferSizes.extent(0), 3u);
  EXPECT_EQ(buffers.buffers.extent(0), 0u);
}


TEST(DeviceMPIBuffers, ConstructorSizesAndBuffersGiven)
{
  DeviceBuffers::BufferSizesView bufferSizes("buffer_sizes", 3);
  set_on_device(bufferSizes, {5, 2, 4});

  DeviceBuffers::BufferView bufferData("buffers", 5+2+4);
  DeviceBuffers buffers(bufferSizes, bufferData);

  EXPECT_EQ(buffers.bufferSizes.extent(0), 3u);
  EXPECT_EQ(buffers.buffers.extent(0), 5u + 2u + 4u);
}

namespace {
void test_device_buffers_set_data()
{
  DeviceBuffers::BufferSizesView bufferSizes("buffer_sizes", 3);
  set_on_device(bufferSizes, {5, 2, 4});

  int totalSize = 5 + 2 + 4;
  DeviceBuffers::BufferView bufferData("buffers", totalSize);
  DeviceBuffers buffers(bufferSizes, bufferData);

  Kokkos::RangePolicy<ExecutionSpace> policy(0, totalSize);
  auto setFunc = KOKKOS_LAMBDA(int i)
  {
    buffers.buffers(i) = i;
  };

  Kokkos::parallel_for("set_buffer_vals", policy, setFunc);

  auto buffersHost = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, buffers.buffers);
  for (int i=0; i < totalSize; ++i)
  {
    EXPECT_EQ(buffersHost(i), i);
  }
}
}

TEST(DeviceMPIBuffers, SetData)
{
  test_device_buffers_set_data();
}

namespace {
template <typename ExecutionSpace>
struct DeviceAppenderUser2Ranks
{
  DeviceAppenderUser2Ranks(int numEntriesOnRank1, ExecutionSpace execSpace) :
    m_numEntriesOnRank1(numEntriesOnRank1),
    appender(2, execSpace)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    int rank = i < m_numEntriesOnRank1 ? 0 : 1;
    appender.push_back(rank, i);
  }

  int m_numEntriesOnRank1;
  DeviceBufferAppender appender;
};
}

TEST(DeviceMPIBufferAppender, push_back)
{
  Kokkos::DefaultHostExecutionSpace hostSpace{};
  ExecutionSpace execSpace{};
  int numRanks = 2;
  Kokkos::View<int*, ExecutionSpace> bufferExpectedSizes("expected_buffer_sizes", numRanks);
  set_on_device(bufferExpectedSizes, {3, 4});
  auto bufferExpectedSizesHost = Kokkos::create_mirror_view_and_copy(hostSpace, bufferExpectedSizes);
  int totalSize = bufferExpectedSizesHost(0) + bufferExpectedSizesHost(1);

  Kokkos::RangePolicy<> policy(execSpace, 0, totalSize);
  DeviceAppenderUser2Ranks setValuesFunc(bufferExpectedSizesHost(0), execSpace);

  Kokkos::parallel_for("size_buffers", policy, setValuesFunc);
  execSpace.fence();
  auto bufferSizesHost = Kokkos::create_mirror_view_and_copy(hostSpace, setValuesFunc.appender.getBuffers().bufferSizes);
  EXPECT_EQ(bufferSizesHost(0), bufferExpectedSizesHost(0));
  EXPECT_EQ(bufferSizesHost(1), bufferExpectedSizesHost(1));
  EXPECT_EQ(setValuesFunc.appender.getBuffers().buffers.extent(0), 0u);
  setValuesFunc.appender.allocate_buffers();
  EXPECT_EQ(setValuesFunc.appender.getBuffers().buffers.extent(0), size_t(totalSize));

  Kokkos::parallel_for("fill_buffers", policy, setValuesFunc);
  execSpace.fence();
  DeviceBuffers deviceBuffers = setValuesFunc.appender.getBuffers();

  auto buffersHost = Kokkos::create_mirror_view_and_copy(hostSpace, deviceBuffers.buffers);
  auto rank0Buffer = Kokkos::subview(buffersHost, std::pair<int, int>(0,                          bufferExpectedSizesHost(0)));
  auto rank1Buffer = Kokkos::subview(buffersHost, std::pair<int, int>(bufferExpectedSizesHost(0), buffersHost.extent(0)));
  // the order the values get put into the buffer are undefined, sort them for testing
  Kokkos::sort(rank0Buffer);
  Kokkos::sort(rank1Buffer);

  for (size_t i=0; i < rank0Buffer.extent(0); ++i)
  {
    EXPECT_EQ(rank0Buffer(i), int(i));
  }

  for (size_t i=0; i < rank1Buffer.extent(0); ++i)
  {
    EXPECT_EQ(rank1Buffer(i), int(i + bufferExpectedSizesHost(0)));
  }
}

void test_device_data_exchange_all_to_all(int numValuesPerRank)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  HostSpace hostSpace{};
  ExecutionSpace execSpace{};
  MPI_Comm comm = stk::parallel_machine_world();
  int commRank = stk::parallel_machine_rank(comm);
  int commSize = stk::parallel_machine_size(comm);

  DeviceBuffers sendBuffers(commSize);
  Kokkos::resize(sendBuffers.buffers, commSize*numValuesPerRank);

  Kokkos::RangePolicy<> setValuesPolicy(execSpace, 0, commSize*numValuesPerRank);
  auto setValuesFunc = KOKKOS_LAMBDA(int i) { sendBuffers.buffers(i) = i + numValuesPerRank * commSize * commRank; };
  Kokkos::parallel_for("set_send_values", setValuesPolicy, setValuesFunc);

  Kokkos::RangePolicy<> setSizesPolicy(execSpace, 0, commSize);
  auto setSizesFunc = KOKKOS_LAMBDA(int i) { sendBuffers.bufferSizes(i) = numValuesPerRank; };
  Kokkos::parallel_for("set_send_sizes", setSizesPolicy, setSizesFunc);

  DeviceDataExchanger exchanger(sendBuffers, execSpace, comm);
  DeviceBuffers recvBuffers = exchanger.communicate();

  auto recvSizesHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.bufferSizes);
  auto recvBuffersHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.buffers);

  EXPECT_EQ(recvSizesHost.extent(0), size_t(commSize));
  EXPECT_EQ(recvBuffersHost.extent(0), size_t(commSize * numValuesPerRank));

  int idx = 0;
  for (int senderRank=0; senderRank < commSize; ++senderRank)
  {
    EXPECT_EQ(recvSizesHost(senderRank), numValuesPerRank);
    int offset = numValuesPerRank * commSize * senderRank + commRank * numValuesPerRank;
    for (int i=0; i < numValuesPerRank; ++i)
    {
      EXPECT_EQ(recvBuffersHost(idx++), i + offset);
    }
  }
}

void test_device_data_exchange_all_to_all_skip_sender(int numValuesPerRank, int senderRankToSkip)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  HostSpace hostSpace{};
  ExecutionSpace execSpace{};
  MPI_Comm comm = stk::parallel_machine_world();
  int commRank = stk::parallel_machine_rank(comm);
  int commSize = stk::parallel_machine_size(comm);

  DeviceBuffers sendBuffers(commSize);
  Kokkos::resize(sendBuffers.buffers, commSize*numValuesPerRank);

  if (commRank != senderRankToSkip)
  {
    Kokkos::RangePolicy<> setValuesPolicy(execSpace, 0, commSize*numValuesPerRank);
    auto setValuesFunc = KOKKOS_LAMBDA(int i) { sendBuffers.buffers(i) = i + numValuesPerRank * commSize * commRank; };
    Kokkos::parallel_for("set_send_values", setValuesPolicy, setValuesFunc);

    Kokkos::RangePolicy<> setSizesPolicy(execSpace, 0, commSize);
    auto setSizesFunc = KOKKOS_LAMBDA(int i) { sendBuffers.bufferSizes(i) = numValuesPerRank; };
    Kokkos::parallel_for("set_send_sizes", setSizesPolicy, setSizesFunc);
  }

  DeviceDataExchanger exchanger(sendBuffers, execSpace, comm);
  DeviceBuffers recvBuffers = exchanger.communicate();

  auto recvSizesHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.bufferSizes);
  auto recvBuffersHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.buffers);

  EXPECT_EQ(recvSizesHost.extent(0), size_t(commSize));
  EXPECT_EQ(recvBuffersHost.extent(0), size_t((commSize-1) * numValuesPerRank));

  int idx = 0;
  for (int senderRank=0; senderRank < commSize; ++senderRank)
  {
    if (senderRank == senderRankToSkip)
    {
      EXPECT_EQ(recvSizesHost(senderRank), 0);
      continue;
    }

    EXPECT_EQ(recvSizesHost(senderRank), numValuesPerRank);
    int offset = numValuesPerRank * commSize * senderRank + commRank * numValuesPerRank;
    for (int i=0; i < numValuesPerRank; ++i)
    {
      EXPECT_EQ(recvBuffersHost(idx++), i + offset);
    }
  }
}


void test_device_data_exchange_all_to_all_skip_receiver(int numValuesPerRank, int receiverToSkip)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  HostSpace hostSpace{};
  ExecutionSpace execSpace{};
  MPI_Comm comm = stk::parallel_machine_world();
  int commRank = stk::parallel_machine_rank(comm);
  int commSize = stk::parallel_machine_size(comm);

  DeviceBuffers sendBuffers(commSize);
  Kokkos::resize(sendBuffers.buffers, commSize*numValuesPerRank);

  Kokkos::RangePolicy<> setValuesPolicy(execSpace, 0, (commSize - 1)*numValuesPerRank);
  auto setValuesFunc = KOKKOS_LAMBDA(int i) { sendBuffers.buffers(i) = i + numValuesPerRank * commSize * commRank; };
  Kokkos::parallel_for("set_send_values", setValuesPolicy, setValuesFunc);

  Kokkos::RangePolicy<> setSizesPolicy(execSpace, 0, commSize);
  auto setSizesFunc = KOKKOS_LAMBDA(int rank) { sendBuffers.bufferSizes(rank) = rank == receiverToSkip ? 0 : numValuesPerRank;};
  Kokkos::parallel_for("set_send_sizes", setSizesPolicy, setSizesFunc);

  DeviceDataExchanger exchanger(sendBuffers, execSpace, comm);
  DeviceBuffers recvBuffers = exchanger.communicate();

  auto recvSizesHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.bufferSizes);
  auto recvBuffersHost = Kokkos::create_mirror_view_and_copy(hostSpace, recvBuffers.buffers);

  EXPECT_EQ(recvSizesHost.extent(0), size_t(commSize));

  if (commRank == receiverToSkip)
  {
    for (int senderRank = 0; senderRank < commSize; ++senderRank)
    {
      EXPECT_EQ(recvSizesHost(senderRank), 0);
    }
    EXPECT_EQ(recvBuffersHost.extent(0), 0U);

    return;
  } else
  {
    EXPECT_EQ(recvBuffersHost.extent(0), size_t(commSize * numValuesPerRank));

    int idx = 0;
    for (int senderRank=0; senderRank < commSize; ++senderRank)
    {
      EXPECT_EQ(recvSizesHost(senderRank), numValuesPerRank);
      int offset = numValuesPerRank * commSize * senderRank;
      if (commRank < receiverToSkip)
      {
        offset += commRank * numValuesPerRank;
      } else
      {
        offset += (commRank - 1) * numValuesPerRank;
      }
      for (int i=0; i < numValuesPerRank; ++i)
      {
        EXPECT_EQ(recvBuffersHost(idx++), i + offset);
      }
    }
  }
}


TEST(DeviceDataExchangeUnknownPattern, AllToAllSmallMsg)
{
  test_device_data_exchange_all_to_all(4);
}

TEST(DeviceDataExchangeUnknownPattern, AllToAllLargeMsg)
{
  test_device_data_exchange_all_to_all(20 * 1024 * 1024 / sizeof(int));
}

TEST(DeviceDataExchangeUnknownPattern, AllToAllSkipSender)
{
  test_device_data_exchange_all_to_all_skip_sender(4, 0);
}

TEST(DeviceDataExchangeUnknownPattern, AllToAllSkipReceiver)
{
  test_device_data_exchange_all_to_all_skip_receiver(4, 0);
}