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
    ProfilingRegion() = delete;
    ProfilingRegion(const ProfilingRegion &other) = delete;
    ProfilingRegion(ProfilingRegion &&other) = delete;
    
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
    return (0 == (uintptr_t(arg.dst) % sizeof(T)))
        && (0 == (uintptr_t(arg.src) & sizeof(T)))
        && (0 == (arg.count % sizeof(T)));
}

template <typename T, typename Member>
KOKKOS_INLINE_FUNCTION void team_memcpy_as(const Member &member, void *dst, void *const src, size_t count) {
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(member, count),
        [&] (size_t i) {
            reinterpret_cast<T *>(dst)[i] = reinterpret_cast<T const *>(src)[i];
        }
    );
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

} // namespace

namespace Tpetra::Details {
  
struct Ialltofewv::Cache::impl {

  impl() : 
    rootBufDev("rootBufDev"), rootBufHost("rootBufHost"),
    aggBufDev("aggBufDev"), aggBufHost("rootBufHost"),
    argsDev("argsDev"), argsHost("argsHost"),
    rootBufGets_(0), rootBufHits_(0),
    aggBufGets_(0), aggBufHits_(0),
    argsGets_(0), argsHits_(0),
    rootBufDevSize_(0), aggBufDevSize_(0), 
    argsDevSize_(0), argsHostSize_(0),
    rootBufHostSize_(0), aggBufHostSize_(0)
{}

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
    if constexpr(std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
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
    if constexpr(std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
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
    if constexpr(std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>) {
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

Ialltofewv::Cache::Cache() : pimpl(std::make_shared<Cache::impl>()) {}
#ifdef NDEBUG
Ialltofewv::Cache::~Cache() = default;
#else
Ialltofewv::Cache::~Cache() {
  if (pimpl->rootBufGets_) {
    std::cerr << "rootBuf:" 
    << " " << pimpl->rootBufDevSize_ + pimpl->rootBufHostSize_
    << " " << pimpl->rootBufHits_ << "/" << pimpl->rootBufGets_ << "\n";
  }
  if (pimpl->aggBufGets_) {
    std::cerr << "aggBuf:" 
    << " " << pimpl->aggBufDevSize_ + pimpl->aggBufHostSize_
    << " " << pimpl->aggBufHits_ << "/" << pimpl->aggBufGets_ << "\n";
  }
  if (pimpl->argsGets_) {
    std::cerr << "args:" 
    << " " << pimpl->argsDevSize_ + pimpl->argsHostSize_
    << " " << pimpl->argsHits_ << "/" << pimpl->argsGets_ << "\n";
  }
}
#endif

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

  const int AGG_TAG = req.tag + 0;
  const int ROOT_TAG = req.tag + 1;

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
  
    // is this rank a root? linear search - nroots expected to be small
    const bool isRoot = std::find(req.roots, req.roots + req.nroots, rank) !=  req.roots + req.nroots;

    // ensure aggregators know how much data each rank is sending to the root
    // [si * nroots + r1] -> how much rank si wants to send to root ri
    std::vector<int> groupSendCounts(size_t(req.nroots) * size_t(srcsPerAgg));
    std::vector<MPI_Request> reqs;
    if (rank == myAgg) {
      reqs.reserve(srcsPerAgg);
      // recv counts from each member of my group
      for (int si = 0; si < srcsPerAgg && si + rank < size; ++si) {
        MPI_Request rreq;

#ifndef NDEBUG
        // if (size_t(si) * req.nroots + req.nroots > groupSendCounts.size()) {
        //   std::stringstream ss;
        //   ss << __FILE__ << ":" << __LINE__ 
        //      << " [" << rank << "]"
        //      << " OOB\n";
        //   std::cerr << ss.str();
        // }
#endif
#ifndef NDEBUG
    // {
    //   std::stringstream ss;
    //   ss << __FILE__ << ":" << __LINE__ 
    //   << " [" << rank << "] ph0 Irecv\n";
    //   std::cerr << ss.str();
    // }
#endif
        MPI_Irecv(&groupSendCounts[size_t(si) * size_t(req.nroots)], req.nroots, MPI_INT, si + rank, 
                  req.tag, req.comm, &rreq);
        reqs.push_back(rreq);
      }
    }
    // send sendcounts to aggregator
#ifndef NDEBUG
    // {
    //   std::stringstream ss;
    //   ss << __FILE__ << ":" << __LINE__ 
    //   << " [" << rank << "] ph0 Send\n";
    //   std::cerr << ss.str();
    // }
#endif
    MPI_Send(req.sendcounts, req.nroots, MPI_INT, myAgg, req.tag, req.comm);
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.resize(0);

#ifndef NDEBUG
    // if (rank == myAgg) {
    //   std::stringstream ss;
    //   ss << __FILE__ << ":" << __LINE__ 
    //      << " [" << rank << "] groupSendCounts=";
    //   for (auto e : groupSendCounts) {
    //     ss << e << " ";
    //   }
    //   ss << "\n";
    //   std::cerr << ss.str();
    // }
#endif

    // at this point, in each aggregator, groupSendCounts holds the send counts
    // from each member of the aggregator. The first nroots entries are from the first rank,
    // the second nroots entries are from the second rank, etc

    // a temporary buffer to aggregate data. Data for a root is contiguous.
#if 0
    Kokkos::View<char *, typename RecvExecSpace::memory_space>aggBuf("aggBuf");
#else
    auto aggBuf = cache.pimpl->get_aggBuf<RecvExecSpace>(0);
#endif
    std::vector<size_t> rootCount(req.nroots, 0); // [ri] the count of data held for root ri
    if (rank == myAgg) {
      size_t aggBytes = 0;
      for (int si = 0; si < srcsPerAgg && si + rank < size; ++si) {
        for (int ri = 0; ri < req.nroots; ++ri) {
          int count = groupSendCounts[si * req.nroots + ri];
          rootCount[ri] += count;
          aggBytes += count * sendSize;
        }
      }

  #ifndef NDEBUG
      // {
      //   std::stringstream ss;
      //   ss << __FILE__ << ":" << __LINE__ 
      //      << " [" << rank << "] rootCount=";
      //   for (auto e : rootCount) {
      //     ss << e << " ";
      //   }
      //   ss << "\n";
      //   std::cerr << ss.str();
      // }
  #endif

#ifndef NDEBUG
      // {
      //   std::stringstream ss;
      //   ss << __FILE__ << ":" << __LINE__ 
      //      << " [" << rank << "] aggBuf.resize(" << aggBytes << ")\n";
      //   std::cerr << ss.str();
      // }
  #endif
#if 0
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), aggBuf, aggBytes);
#else
      aggBuf = cache.pimpl->get_aggBuf<RecvExecSpace>(aggBytes);
#endif
    }
    // now, on the aggregator ranks,
    // * aggBuf is resized to accomodate all incoming data
    // * rootCount holds how much data i hold for each root
    


    // Send the actual data to the aggregator
    reqs.reserve(srcsPerAgg);
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
            // {
            //   std::stringstream ss;
            //   ss << __FILE__ << ":" << __LINE__ 
            //   << " [" << rank << "] ph1 recv(@" << displ << ", " << count << ", ..., " << si+rank << "\n";
            //   std::cerr << ss.str();
            // }
#endif
#ifndef NDEBUG
            if (displ + count * sendSize > aggBuf.size()) {
              std::stringstream ss;
              ss << __FILE__ << ":" << __LINE__ 
              << " [" << rank << "] OOB\n";
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
      reqs.reserve(req.nroots); // prepare for one send per root
    }
    
    // send data to aggregator
    for (int ri = 0; ri < req.nroots; ++ri) {
      const size_t displ = size_t(req.sdispls[ri]) * sendSize;
      const int count = req.sendcounts[ri];
      if (count) {
#ifndef NDEBUG
        // {
        //   std::stringstream ss;
        //   ss << __FILE__ << ":" << __LINE__ 
        //   << " [" << rank << "] ph1 Isend(@" << displ << ", " << count << ", ..., " << myAgg << "\n";
        //   std::cerr << ss.str();
        // }
#endif
        MPI_Request sreq;
        MPI_Isend(&reinterpret_cast<const char *>(req.sendbuf)[displ], req.sendcounts[ri],
        req.sendtype, myAgg, AGG_TAG, req.comm, &sreq);
        reqs.push_back(sreq);
      }
    }
  
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.resize(0);


    // if I am a root, recieve data from each aggregator
    // The aggregator will send contiguous data, which we may need to spread out according to rdispls
#if 0
    Kokkos::View<uint8_t *, typename RecvExecSpace::memory_space>rootBuf("rootBuf");
#else
    auto rootBuf = cache.pimpl->get_rootBuf<RecvExecSpace>(0);
#endif
    if (isRoot) {
      reqs.reserve(naggs); // receive from each aggregator

      const size_t totalRecvd = recvSize * [&]() -> size_t {
        size_t acc = 0;
        for (int i = 0; i < size; ++i) {
          acc += req.recvcounts[i];
        }
        return acc;
      }();
#if 0
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing), rootBuf, totalRecvd);
#else
      rootBuf = cache.pimpl->get_rootBuf<RecvExecSpace>(totalRecvd);
#endif

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
#ifndef NDEBUG
          // {
          //   std::stringstream ss;
          //   ss << __FILE__ << ":" << __LINE__ 
          //   << " [" << rank << "] ph2 Irecv(@" << displ << ", " << count << ", ..., " << aggSrc << "\n";
          //   std::cerr << ss.str();
          // }
#endif
          MPI_Request rreq;
          // &rootBuf(displ)
          MPI_Irecv(rootBuf.data() + displ, count, req.recvtype, aggSrc, ROOT_TAG,  req.comm, &rreq);
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
#ifndef NDEBUG
          // {
          //   std::stringstream ss;
          //   ss << __FILE__ << ":" << __LINE__ 
          //   << " [" << rank << "] ph2 send(@" << displ << ", " << count << ", ..., " << req.roots[ri] << "\n";
          //   std::cerr << ss.str();
          // }
#endif
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
#if 0
      Kokkos::View<MemcpyArg*> args(Kokkos::view_alloc("args", Kokkos::WithoutInitializing), size);
#else
      auto args = cache.pimpl->get_args<RecvExecSpace>(size);
#endif
      auto args_h = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, args);

      size_t srcOff = 0;
      for (int sRank = 0; sRank < size; ++sRank) {
        const size_t dstOff = req.rdispls[sRank] * recvSize;
        
        void *dst = &reinterpret_cast<char *>(req.recvbuf)[dstOff];
        void *const src = rootBuf.data() + srcOff; // &rootBuf(srcOff);
        const size_t count = req.recvcounts[sRank] * recvSize;
        args_h(sRank) = MemcpyArg{dst, src, count};

#ifndef NDEBUG
        if (srcOff + count > rootBuf.extent(0)) {
          std::stringstream ss;
          ss << __FILE__ << ":" << __LINE__ << " OOB!\n";
          std::cerr << ss.str();
        }
#endif
        srcOff += count;
      }

      // Actually copy the data
      Kokkos::deep_copy(args, args_h);
      using Policy = Kokkos::TeamPolicy<RecvExecSpace>;
      Policy policy(size, Kokkos::AUTO);
      Kokkos::parallel_for("fixup rdispl", policy, 
        KOKKOS_LAMBDA(typename Policy::member_type member){
          team_memcpy(member, args(member.league_rank()));
        }
      );
      Kokkos::fence("after fixup rdispl");

    }

    return finalize();
}
} //  namespace


int Ialltofewv::wait(Req &req) {
  if (req.devAccess) {
    return wait_impl<Kokkos::DefaultExecutionSpace>(req, cache_);
  } else {
    return wait_impl<Kokkos::DefaultHostExecutionSpace>(req, cache_);
  }
}

int Ialltofewv::get_status(const Req &req, int *flag, MPI_Status */*status*/) const {
  *flag = req.completed;
  return MPI_SUCCESS;
}


} // namespace Tpetra::Details