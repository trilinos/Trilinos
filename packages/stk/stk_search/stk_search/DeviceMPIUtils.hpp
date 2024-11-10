#ifndef STK_SEARCH_DEVICE_MPI_UTILS_HPP
#define STK_SEARCH_DEVICE_MPI_UTILS_HPP

#include <numeric>
#include "Kokkos_Core.hpp"
#include "stk_util/parallel/DeviceAwareMPI.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_util/parallel/DeviceAwareMPI.hpp"

namespace stk::search {

namespace impl {

// throws an error if the view data layout is different than a C array (or std::vector).
// This can happen if Kokkos puts padding in, for example for memory alignment reasons.
// If no error is thrown, the view can be passed into MPI as a send or receive buffer.
template <typename ViewType>
void check_view_c_layout(ViewType view)
{
  static_assert(ViewType::rank() == 1);
  if (view.extent(0) == 0)
  {
    return;
  }

  size_t viewSizeBytes = (&(view(view.extent(0) - 1)) - &(view(0)) + 1) * sizeof(typename ViewType::value_type);
  size_t cArraySize = sizeof(typename ViewType::value_type) * view.extent(0);
  STK_ThrowRequireMsg(viewSizeBytes == cArraySize, "The view must not have any padding");
}

template <typename BufferSizesView, typename BuffersView>
void check_buffer_sizes(const BufferSizesView bufferSizes, const BuffersView buffers)
{
  using ExecutionSpace = typename BufferSizesView::execution_space;
  Kokkos::RangePolicy<ExecutionSpace> policy(0, bufferSizes.extent(0));
  auto func = KOKKOS_LAMBDA(int i, int& lsum)
  {
    lsum += bufferSizes(i);
  };

  int totalSize = 0;
  Kokkos::parallel_reduce("check DeviceMPIBuffers size", policy, func, totalSize);
  STK_ThrowRequireMsg(size_t(totalSize) == buffers.extent(0), "the buffer view must be the the total size of the buffers");
}


template <typename T, typename ExecutionSpace>
struct DeviceMPIBuffers
{
  using BufferSizesView = Kokkos::View<int*, ExecutionSpace>;
  using BufferView      = Kokkos::View<T*, ExecutionSpace>;
  using ValueType       = T;

  DeviceMPIBuffers() :
    bufferSizes("mpi_buffer_sizes", 0),
    buffers("mpi_buffers", 0)
  {}

  DeviceMPIBuffers(int numRanks) :
    bufferSizes("mpi_buffer_sizes", numRanks),
    buffers("mpi_buffers", 0)
  {}

  DeviceMPIBuffers(BufferSizesView _bufferSizes, BufferView _buffers) :
    bufferSizes(_bufferSizes),
    buffers(_buffers)
  {
#ifndef NDEBUG
    check_buffer_sizes(bufferSizes, buffers);
#endif
  }

  BufferSizesView bufferSizes;
  BufferView buffers;
};

template <typename T, typename ExecutionSpace>
class DeviceMPIBufferAppender
{
  public:
    using DeviceBuffers = DeviceMPIBuffers<T, ExecutionSpace>;

    DeviceMPIBufferAppender() = default;

    DeviceMPIBufferAppender(int numBuffers, ExecutionSpace execSpace) :
      m_deviceBuffers(numBuffers),
      m_bufferIdxs("mpi_buffer_current_idx", numBuffers),
      m_execSpace(execSpace),
      m_areBuffersSized("areBuffersSized")
    {
      Kokkos::deep_copy(m_areBuffersSized, false);
    }

    KOKKOS_INLINE_FUNCTION
    void push_back(int rank, const T& val) const
    {
      if (!(m_areBuffersSized()))
      {
        Kokkos::atomic_add(&m_deviceBuffers.bufferSizes(rank), 1);
      } else
      {
        int idx = Kokkos::atomic_fetch_add(&m_bufferIdxs(rank), 1);
#ifndef NDEBUG
        int offset = 0;
        for (int i=0; i < rank; ++i)
        {
          offset += m_deviceBuffers.bufferSizes(i);
        }
        KOKKOS_ASSERT(idx >= offset && idx < (m_deviceBuffers.bufferSizes(rank) + offset))
#endif
        m_deviceBuffers.buffers(idx) = val;
      }
    }

    void allocate_buffers()
    {
      bool areBuffersSized;
      Kokkos::deep_copy(areBuffersSized, m_areBuffersSized);
      STK_ThrowRequireMsg(!areBuffersSized, "DeviceMPIBufferAppender can only be used once");

      Kokkos::View<int, ExecutionSpace> totalSize("total_size");

      Kokkos::RangePolicy<ExecutionSpace> policy(m_execSpace, 0, m_bufferIdxs.extent(0));
      int lastIdx             = m_deviceBuffers.bufferSizes.extent(0) - 1;
      auto deviceBuffersLocal = m_deviceBuffers;
      auto bufferIdxsLocal    = m_bufferIdxs;
      auto func = KOKKOS_LAMBDA(int i, int& update, bool final)
      {
        if (final)
        {
          bufferIdxsLocal(i) = update;
          if (i == lastIdx)
          {
            totalSize() = update + deviceBuffersLocal.bufferSizes(lastIdx);
          }
        }

        update += deviceBuffersLocal.bufferSizes(i);
      };

      Kokkos::parallel_scan("set_buffer_offsets", policy, func);
      m_execSpace.fence();

      int totalSizeHost;
      Kokkos::deep_copy(totalSizeHost, totalSize);
      Kokkos::resize(Kokkos::WithoutInitializing, m_deviceBuffers.buffers, totalSizeHost);
      Kokkos::deep_copy(m_areBuffersSized, true);
    }

    DeviceBuffers getBuffers() const { return m_deviceBuffers; }

  private:
    DeviceBuffers m_deviceBuffers;
    Kokkos::View<int*, ExecutionSpace> m_bufferIdxs;
    ExecutionSpace m_execSpace;

    Kokkos::View<bool, ExecutionSpace> m_areBuffersSized;
};

template <typename T, typename ExecutionSpace>
class DeviceDataExchangeUnknownPattern
{
  public:
    using DeviceBuffers = impl::DeviceMPIBuffers<T, ExecutionSpace>;
    using HostMirrorSpace = typename DeviceBuffers::BufferSizesView::host_mirror_space;

    using BufferSizesHostView = typename DeviceBuffers::BufferSizesView::HostMirror;
    using BufferHostView      = typename DeviceBuffers::BufferView::HostMirror;

    // Note: running the ConeCrush and Jenga performance tests on Vortex showed device MPI
    //       was never faster than host MPI and sometimes ~5% slower.  Setting useDeviceMPI
    //       to false for now, maybe things will be better on ATS-4
    explicit DeviceDataExchangeUnknownPattern(const DeviceBuffers& sendBuffers, ExecutionSpace execSpace,
                                              MPI_Comm comm, bool useDeviceMPI=false) :
      m_sendBuffers(sendBuffers),
      m_recvBuffers(stk::parallel_machine_size(comm)),
      m_sendBufferSizesHost(Kokkos::create_mirror_view_and_copy(HostMirrorSpace{}, sendBuffers.bufferSizes)),
      m_sendBufferHost(),
      m_comm(comm),
      m_useDeviceMPI(useDeviceMPI)
    {
      if (useDeviceMPI && !stk::have_device_aware_mpi())
      {
        STK_ThrowRequireMsg(false, "STK did detect device-aware MPI on this system");
      }

      if (!useDeviceMPI)
      {
        m_sendBufferHost = Kokkos::create_mirror_view_and_copy(HostMirrorSpace{}, sendBuffers.buffers);
      }
    }

    DeviceBuffers communicate()
    {
      int commSize = stk::parallel_machine_size(m_comm);
      MPITag tag = get_mpi_tag_manager().get_tag(m_comm);

      Kokkos::Profiling::pushRegion("Compute Recv lists");
      std::vector<int> recvSizes = getRecvCounts();
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("Create recv buffers");
      createReceiveBuffers(recvSizes);
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("Post sends and receives");
      std::vector<MPI_Request> sendRequests(commSize, MPI_REQUEST_NULL);
      std::vector<MPI_Request> recvRequests(commSize, MPI_REQUEST_NULL);

      postReceives(recvSizes, recvRequests, tag);
      postSends(sendRequests, tag);
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("Wait for MPI completion");
      MPI_Waitall(commSize, recvRequests.data(), MPI_STATUSES_IGNORE);
      if (!m_useDeviceMPI)
      {
        Kokkos::deep_copy(m_execSpace, m_recvBuffers.buffers, m_recvBufferHost);
      }
      MPI_Waitall(commSize, sendRequests.data(), MPI_STATUSES_IGNORE);
      Kokkos::Profiling::popRegion();

      return m_recvBuffers;
    }

  private:
    void createReceiveBuffers(const std::vector<int>& recvSizes)
    {
      using BufferSizesVectorWrap = Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      BufferSizesVectorWrap recvBufferSizesHost(Kokkos::view_wrap(recvSizes.data()), recvSizes.size());
      Kokkos::deep_copy(m_recvBuffers.bufferSizes, recvBufferSizesHost);

      int totalSize = std::accumulate(recvSizes.begin(), recvSizes.end(), int(0));
      Kokkos::resize(Kokkos::WithoutInitializing, m_recvBuffers.buffers, totalSize);
      if (!m_useDeviceMPI)
      {
        m_recvBufferHost = Kokkos::create_mirror_view(m_recvBuffers.buffers);
      }
    }

    int postReceives(const std::vector<int>& recvSizes, std::vector<MPI_Request>& recvRequests, stk::MPITag tag)
    {
      size_t recvBufOffset = 0;
      int numRecvs = 0;
      T* recvBufferData = m_useDeviceMPI ? m_recvBuffers.buffers.data() : m_recvBufferHost.data();
      for (int rank=0; rank < int(recvRequests.size()); ++rank)
      {
        if (recvSizes[rank] > 0)
        {
          MPI_Irecv(recvBufferData + recvBufOffset, recvSizes[rank]*sizeof(T), MPI_CHAR, rank, tag, m_comm, &recvRequests[rank]);
          numRecvs++;
        }

        recvBufOffset += recvSizes[rank];
      }

      return numRecvs;
    }

    std::vector<int> getRecvCounts()
    {
      auto tag = stk::get_mpi_tag_manager().get_tag(m_comm);
      check_view_c_layout(m_sendBufferSizesHost);

      int numItems;
      MPI_Reduce_scatter_block(&(m_sendBufferSizesHost(0)), &numItems, 1, MPI_INT, MPI_SUM, m_comm);

      int numProcs = m_sendBufferSizesHost.extent(0);
      std::vector<int> recvCounts(numProcs, 0);
      std::vector<MPI_Request> sendReqs(numProcs, MPI_REQUEST_NULL);

      for (int proc=0; proc < numProcs; ++proc)
      {
        if (m_sendBufferSizesHost(proc) > 0)
        {
          MPI_Isend(&(m_sendBufferSizesHost(proc)), 1, MPI_INT, proc, tag, m_comm, &(sendReqs[proc]));
        }
      }

      //TODO: this doesn't need to be a separate step.  We could do this
      //      as part of the receives, and then run a kernel to permute
      //      the data (or modify the API of DeviceBuffers to allow the data to be permuted)
      int recvCount = 0;
      while (recvCount != numItems)
      {
        MPI_Status status;
        int itemCount = 0;
        MPI_Recv(&itemCount, 1, MPI_INT, MPI_ANY_SOURCE, tag, m_comm, &status);
        int sender = status.MPI_SOURCE;

        recvCounts[sender] = itemCount;
        recvCount += itemCount;
      }

      MPI_Waitall(numProcs, sendReqs.data(), MPI_STATUSES_IGNORE);

      return recvCounts;
    }

    int postSends(std::vector<MPI_Request>& sendRequests, stk::MPITag tag)
    {
      size_t sendBufOffset = 0;
      int numSends = 0;
      const T* sendBufferData = m_useDeviceMPI ? m_sendBuffers.buffers.data() : m_sendBufferHost.data();
      for (int rank=0; rank < int(sendRequests.size()); ++rank)
      {
        if (m_sendBufferSizesHost(rank) > 0)
        {
          MPI_Isend(sendBufferData + sendBufOffset, m_sendBufferSizesHost(rank)*sizeof(T), MPI_CHAR, rank, tag, m_comm, &sendRequests[rank]);
          numSends++;
        }

        sendBufOffset += m_sendBufferSizesHost(rank);
      }

      return numSends;
    }


    DeviceBuffers m_sendBuffers;
    DeviceBuffers m_recvBuffers;

    BufferSizesHostView m_sendBufferSizesHost;
    BufferHostView m_sendBufferHost;

    BufferHostView m_recvBufferHost;

    ExecutionSpace m_execSpace;
    MPI_Comm m_comm;
    bool m_useDeviceMPI;
};

}
}

#endif
