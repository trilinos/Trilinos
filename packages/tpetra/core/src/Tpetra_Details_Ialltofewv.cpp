// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_Ialltofewv.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <mpi.h>
#include <Kokkos_Core.hpp>

#ifndef NDEBUG
#include <iostream>
#include <sstream>
#endif

namespace {

struct ProfilingRegion {
  ProfilingRegion()                             = delete;
  ProfilingRegion(const ProfilingRegion &other) = delete;
  ProfilingRegion(ProfilingRegion &&other)      = delete;

  ProfilingRegion(const std::string &name) {
    Kokkos::Profiling::pushRegion(name);
  }
  ~ProfilingRegion() {
    Kokkos::Profiling::popRegion();
  }
};

struct MemcpyArg {
  void *dst;
  void *src;
  size_t count;
};

template <typename T>
KOKKOS_INLINE_FUNCTION bool is_compatible(const MemcpyArg &arg) {
  return (0 == (uintptr_t(arg.dst) % sizeof(T))) && (0 == (uintptr_t(arg.src) & sizeof(T))) && (0 == (arg.count % sizeof(T)));
}

template <typename T, typename Member>
KOKKOS_INLINE_FUNCTION void team_memcpy_as(const Member &member, void *dst, void *const src, size_t count) {
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, count),
      [&](size_t i) {
        reinterpret_cast<T *>(dst)[i] = reinterpret_cast<T const *>(src)[i];
      });
}

template <typename Member>
KOKKOS_INLINE_FUNCTION void team_memcpy(const Member &member, MemcpyArg &arg) {
  if (is_compatible<uint64_t>(arg)) {
    team_memcpy_as<uint64_t>(member, arg.dst, arg.src, arg.count / sizeof(uint64_t));
  } else if (is_compatible<uint32_t>(arg)) {
    team_memcpy_as<uint32_t>(member, arg.dst, arg.src, arg.count / sizeof(uint32_t));
  } else {
    team_memcpy_as<uint8_t>(member, arg.dst, arg.src, arg.count);
  }
}

}  // namespace

namespace Tpetra::Details {

struct Ialltofewv::Cache::impl {
  impl()
    : rootBufDev("rootBufDev", 0)
    , rootBufHost("rootBufHost", 0)
    , aggBufDev("aggBufDev", 0)
    , aggBufHost("rootBufHost", 0)
    , argsDev("argsDev", 0)
    , argsHost("argsHost", 0)
    , rootBufGets_(0)
    , rootBufHits_(0)
    , aggBufGets_(0)
    , aggBufHits_(0)
    , argsGets_(0)
    , argsHits_(0)
    , rootBufDevSize_(0)
    , aggBufDevSize_(0)
    , argsDevSize_(0)
    , argsHostSize_(0)
    , rootBufHostSize_(0)
    , aggBufHostSize_(0) {}

  // cached views
  Kokkos::View<uint8_t *, typename Kokkos::DefaultExecutionSpace::memory_space> rootBufDev;
  Kokkos::View<uint8_t *, typename Kokkos::DefaultHostExecutionSpace::memory_space> rootBufHost;
  Kokkos::View<char *, typename Kokkos::DefaultExecutionSpace::memory_space> aggBufDev;
  Kokkos::View<char *, typename Kokkos::DefaultHostExecutionSpace::memory_space> aggBufHost;
  Kokkos::View<MemcpyArg *, typename Kokkos::DefaultExecutionSpace::memory_space> argsDev;
  Kokkos::View<MemcpyArg *, typename Kokkos::DefaultHostExecutionSpace::memory_space> argsHost;

  size_t rootBufGets_;
  size_t rootBufHits_;
  size_t aggBufGets_;
  size_t aggBufHits_;
  size_t argsGets_;
  size_t argsHits_;

  size_t rootBufDevSize_, aggBufDevSize_;
  size_t argsDevSize_, argsHostSize_;
  size_t rootBufHostSize_, aggBufHostSize_;

  template <typename ExecSpace>
  auto get_rootBuf(size_t size) {
    ++rootBufGets_;
    if constexpr (std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
      if (rootBufDev.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, rootBufDev, size);
        rootBufDevSize_ = size;
      } else {
        ++rootBufHits_;
      }
      return Kokkos::subview(rootBufDev, Kokkos::pair{size_t(0), size});
    } else {
      if (rootBufHost.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, rootBufHost, size);
        rootBufHostSize_ = size;
      } else {
        ++rootBufHits_;
      }
      return Kokkos::subview(rootBufHost, Kokkos::pair{size_t(0), size});
    }
  }

  template <typename ExecSpace>
  auto get_aggBuf(size_t size) {
    ++aggBufGets_;
    if constexpr (std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
      if (aggBufDev.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, aggBufDev, size);
        aggBufHostSize_ = size;
      } else {
        ++aggBufHits_;
      }
      return Kokkos::subview(aggBufDev, Kokkos::pair{size_t(0), size});
    } else {
      if (aggBufHost.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, aggBufHost, size);
        aggBufHostSize_ = size;
      } else {
        ++aggBufHits_;
      }
      return Kokkos::subview(aggBufHost, Kokkos::pair{size_t(0), size});
    }
  }

  template <typename ExecSpace>
  auto get_args(size_t size) {
    ++argsGets_;
    if constexpr (std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
      if (argsDev.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, argsDev, size);
        argsHostSize_ = size;
      } else {
        ++argsHits_;
      }
      return Kokkos::subview(argsDev, Kokkos::pair{size_t(0), size});
    } else {
      if (argsHost.extent(0) < size) {
        Kokkos::resize(Kokkos::WithoutInitializing, argsHost, size);
        argsHostSize_ = size;
      } else {
        ++argsHits_;
      }
      return Kokkos::subview(argsHost, Kokkos::pair{size_t(0), size});
    }
  }
};

Ialltofewv::Cache::Cache()  = default;
Ialltofewv::Cache::~Cache() = default;

namespace {
template <typename RecvExecSpace>
int wait_impl(Ialltofewv::Req &req, Ialltofewv::Cache &cache) {
  auto finalize = [&]() -> int {
    req.completed = true;
    return MPI_SUCCESS;
  };

  if (0 == req.nroots) {
    return finalize();
  }

  ProfilingRegion pr("alltofewv::wait");

  // lazy-init view cache
  if (!cache.pimpl) {
    cache.pimpl = std::make_shared<Ialltofewv::Cache::impl>();
  }

  const int rank = [&]() -> int {
    int _rank;
    MPI_Comm_rank(req.comm, &_rank);
    return _rank;
  }();

  const int size = [&]() -> int {
    int _size;
    MPI_Comm_size(req.comm, &_size);
    return _size;
  }();

  const size_t sendSize = [&]() -> size_t {
    int _size;
    MPI_Type_size(req.sendtype, &_size);
    return _size;
  }();

  const size_t recvSize = [&]() -> size_t {
    int _size;
    MPI_Type_size(req.recvtype, &_size);
    return _size;
  }();

  // is this rank a root? linear search - nroots expected to be small
  const bool isRoot = std::find(req.roots, req.roots + req.nroots, rank) != req.roots + req.nroots;

  const int AGG_TAG  = req.tag + 0;
  const int ROOT_TAG = req.tag + 1;

  // Balance the number of incoming messages at each phase:
  // Aggregation = size / naggs * nroots
  // Root =        naggs
  // so
  // size / naggs * nroots = naggs
  // size * nroots = naggs^2
  // naggs = sqrt(size * nroots)
  const int naggs = std::sqrt(size_t(size) * size_t(req.nroots)) + /*rounding*/ 0.5;

  // how many srcs go to each aggregator
  const int srcsPerAgg = (size + naggs - 1) / naggs;

  // the aggregator I send to
  const int myAgg = rank / srcsPerAgg * srcsPerAgg;

  // ensure aggregators know how much data each rank is sending to the root
  // [si * nroots + r1] -> how much si'th rank in group wants to send to root ri
  std::vector<int> groupSendCounts(size_t(req.nroots) * size_t(srcsPerAgg));
  std::vector<MPI_Request> reqs;
  if (rank == myAgg) {
    reqs.reserve(srcsPerAgg);
    // recv counts from each member of my group
    for (int si = 0; si < srcsPerAgg && si + rank < size; ++si) {
      MPI_Request rreq;

#ifndef NDEBUG
      if (size_t(si) * req.nroots + req.nroots > groupSendCounts.size()) {
        std::stringstream ss;
        ss << __FILE__ << ":" << __LINE__
           << " [" << rank << "] tpetra internal Ialltofewv error: OOB access in recv buffer\n";
        std::cerr << ss.str();
      }
#endif
      MPI_Irecv(&groupSendCounts[size_t(si) * size_t(req.nroots)], req.nroots, MPI_INT, si + rank,
                req.tag, req.comm, &rreq);
      reqs.push_back(rreq);
    }
  }
  // send sendcounts to aggregator
  MPI_Send(req.sendcounts, req.nroots, MPI_INT, myAgg, req.tag, req.comm);

  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
  reqs.resize(0);

  // at this point, in each aggregator, groupSendCounts holds the send counts
  // from each member of the aggregator. The first nroots entries are from the first rank,
  // the second nroots entries are from the second rank, etc

  // a temporary buffer to aggregate data. Data for a root is contiguous.
  auto aggBuf = cache.pimpl->get_aggBuf<RecvExecSpace>(0);
  std::vector<size_t> rootCount(req.nroots, 0);  // [ri] the count of data held for root ri
  if (rank == myAgg) {
    size_t aggBytes = 0;
    for (int si = 0; si < srcsPerAgg && si + rank < size; ++si) {
      for (int ri = 0; ri < req.nroots; ++ri) {
        int count = groupSendCounts[si * req.nroots + ri];
        rootCount[ri] += count;
        aggBytes += count * sendSize;
      }
    }

    aggBuf = cache.pimpl->get_aggBuf<RecvExecSpace>(aggBytes);
  }
  // now, on the aggregator ranks,
  // * aggBuf is resized to accomodate all incoming data
  // * rootCount holds how much data i hold for each root

  // Send the actual data to the aggregator
  if (rank == myAgg) {
    reqs.reserve(srcsPerAgg + req.nroots);
    // receive from all ranks in my group
    size_t displ = 0;
    // senders will send in root order, so we will recv in that order as well
    // this puts all data for a root contiguous in the aggregation buffer
    for (int ri = 0; ri < req.nroots; ++ri) {
      for (int si = 0; si < srcsPerAgg && si + rank < size; ++si) {
        // receive data for the ri'th root from the si'th sender
        const int count = groupSendCounts[si * req.nroots + ri];
        if (count) {
#ifndef NDEBUG
          if (displ + count * sendSize > aggBuf.size()) {
            std::stringstream ss;
            ss << __FILE__ << ":" << __LINE__
               << " [" << rank << "] tpetra internal Ialltofewv error: OOB access in send buffer\n";
            std::cerr << ss.str();
          }
#endif
          MPI_Request rreq;
          // &aggBuf(displ)
          MPI_Irecv(aggBuf.data() + displ, count, req.sendtype, si + rank, req.tag, req.comm, &rreq);
          reqs.push_back(rreq);
          displ += size_t(count) * sendSize;
        }
      }
    }
  } else {
    reqs.reserve(req.nroots);  // prepare for one send per root
  }

  // send data to aggregator
  for (int ri = 0; ri < req.nroots; ++ri) {
    const size_t displ = size_t(req.sdispls[ri]) * sendSize;
    const int count    = req.sendcounts[ri];
    if (count) {
      MPI_Request sreq;
      MPI_Isend(&reinterpret_cast<const char *>(req.sendbuf)[displ], req.sendcounts[ri],
                req.sendtype, myAgg, AGG_TAG, req.comm, &sreq);
      reqs.push_back(sreq);
    }
  }

  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
  reqs.resize(0);

  // if I am a root, receive data from each aggregator
  // The aggregator will send contiguous data, which we may need to spread out according to rdispls
  auto rootBuf = cache.pimpl->get_rootBuf<RecvExecSpace>(0);
  if (isRoot) {
    reqs.reserve(naggs);  // receive from each aggregator

    const size_t totalRecvd = recvSize * [&]() -> size_t {
      size_t acc = 0;
      for (int i = 0; i < size; ++i) {
        acc += req.recvcounts[i];
      }
      return acc;
    }();
    rootBuf = cache.pimpl->get_rootBuf<RecvExecSpace>(totalRecvd);

    // Receive data from each aggregator.
    // Aggregators send data in order of the ranks they're aggregating,
    // which is also the order the root needs in its recv buffer.
    size_t displ = 0;
    for (int aggSrc = 0; aggSrc < size; aggSrc += srcsPerAgg) {
      // tally up the total data to recv from the sending aggregator
      int count = 0;
      for (int origSrc = aggSrc;
           origSrc < aggSrc + srcsPerAgg && origSrc < size; ++origSrc) {
        count += req.recvcounts[origSrc];
      }

      if (count) {
        MPI_Request rreq;
        MPI_Irecv(rootBuf.data() + displ, count, req.recvtype, aggSrc, ROOT_TAG, req.comm, &rreq);
        reqs.push_back(rreq);
        displ += size_t(count) * recvSize;
      }
    }
  }

  // if I am an aggregator, forward data to the roots
  // To each root, send my data in order of the ranks that sent to me
  // which is the order the recvers expect
  if (rank == myAgg) {
    size_t displ = 0;
    for (int ri = 0; ri < req.nroots; ++ri) {
      const size_t count = rootCount[ri];
      if (count) {
        // &aggBuf[displ]
        MPI_Send(aggBuf.data() + displ, count, req.sendtype, req.roots[ri], ROOT_TAG, req.comm);
        displ += count * sendSize;
      }
    }
  }

  MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

  // at root, copy data from contiguous buffer into recv buffer
  if (isRoot) {
    // set up src and dst for each block
    auto args   = cache.pimpl->get_args<RecvExecSpace>(size);
    auto args_h = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, args);

    size_t srcOff = 0;
    for (int sRank = 0; sRank < size; ++sRank) {
      const size_t dstOff = req.rdispls[sRank] * recvSize;

      void *dst          = &reinterpret_cast<char *>(req.recvbuf)[dstOff];
      void *const src    = rootBuf.data() + srcOff;  // &rootBuf(srcOff);
      const size_t count = req.recvcounts[sRank] * recvSize;
      args_h(sRank)      = MemcpyArg{dst, src, count};

#ifndef NDEBUG
      if (srcOff + count > rootBuf.extent(0)) {
        std::stringstream ss;
        ss << __FILE__ << ":" << __LINE__ << " Tpetra internal Ialltofewv error: src access OOB in memcpy\n";
        std::cerr << ss.str();
      }
#endif
      srcOff += count;
    }

    // Actually copy the data
    Kokkos::deep_copy(args, args_h);
    using Policy = Kokkos::TeamPolicy<RecvExecSpace>;
    Policy policy(size, Kokkos::AUTO);
    Kokkos::parallel_for(
        "Tpetra::Details::Ialltofewv: apply rdispls to contiguous root buffer", policy,
        KOKKOS_LAMBDA(typename Policy::member_type member) {
          team_memcpy(member, args(member.league_rank()));
        });
    Kokkos::fence("Tpetra::Details::Ialltofewv: after apply rdispls to contiguous root buffer");
  }

  return finalize();
}
}  //  namespace

int Ialltofewv::wait(Req &req) {
  if (req.devAccess) {
    return wait_impl<Kokkos::DefaultExecutionSpace>(req, cache_);
  } else {
    return wait_impl<Kokkos::DefaultHostExecutionSpace>(req, cache_);
  }
}

int Ialltofewv::get_status(const Req &req, int *flag, MPI_Status * /*status*/) const {
  *flag = req.completed;
  return MPI_SUCCESS;
}

}  // namespace Tpetra::Details
