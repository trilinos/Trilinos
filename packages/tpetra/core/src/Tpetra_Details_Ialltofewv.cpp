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

template <typename ArgsView, typename Member>
struct ApplyRdispls {
  ArgsView args;

  KOKKOS_INLINE_FUNCTION void operator()(const Member &member) const {
    team_memcpy(member, args(member.league_rank()));
  }
};

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

Ialltofewv::Cache::Cache() = default;

// DistributorActor can be copied, but we don't want to share a single
// cache since each actor can have its own Ialltofewv in flight.
// This creates an empty cache in the copy.s
Ialltofewv::Cache::Cache(const Cache &) {}

// As with copy construction, assignment leaves this cache indepdendent.
Ialltofewv::Cache &Ialltofewv::Cache::operator=(const Cache &other) {
  if (this != &other) pimpl.reset();
  return *this;
}

Ialltofewv::Cache::~Cache() = default;

// Intermediate buffers and MPI requests for one asynchronous
// Ialltofewv. Inherited by actual imlpementation, which is specialized on the
// execution space.
struct Ialltofewv::State {
  virtual ~State()                = default;
  virtual int start()             = 0;
  virtual int progress()          = 0;
  virtual bool isComplete() const = 0;
  virtual MPI_Comm comm() const   = 0;
  // Copies completion state back to the request that Ialltofewv callers maintain.
  // Needed when progress on one Ialltofewv completes a different operation.
  virtual void updateReq(Req &req) const = 0;
};

namespace {
enum class Stage { counts,
                   payload,
                   forwarding,
                   complete };

int test_all(std::vector<MPI_Request> &reqs, int &flag) {
  if (reqs.empty()) {
    flag = 1;
    return MPI_SUCCESS;
  }
  return MPI_Testall(static_cast<int>(reqs.size()), reqs.data(), &flag, MPI_STATUSES_IGNORE);
}

template <typename RecvExecSpace>
class StateImpl : public Ialltofewv::State {
 public:
  using byte_view_type = Kokkos::View<char *, typename RecvExecSpace::memory_space>;
  using root_view_type = Kokkos::View<uint8_t *, typename RecvExecSpace::memory_space>;
  using args_view_type = Kokkos::View<MemcpyArg *, typename RecvExecSpace::memory_space>;

  StateImpl(Ialltofewv::Req &req, std::shared_ptr<Ialltofewv::Cache::impl> cache)
    : req_(req)
    , cache_(std::move(cache))
    , rank_(getRank())
    , size_(getSize())
    , sendSize_(getTypeSize(req_.sendtype))
    , recvSize_(getTypeSize(req_.recvtype))
    // Is this rank a root? Linear search because nroots is expected to be small.
    , isRoot_(std::find(req_.roots, req_.roots + req_.nroots, rank_) != req_.roots + req_.nroots)
    // Balance the number of incoming messages at each phase:
    // Aggregation = size / naggs * nroots
    // Root =        naggs
    // so size / naggs * nroots = naggs, and naggs = sqrt(size * nroots).
    , naggs_(std::max(1, int(std::sqrt(size_t(size_) * size_t(req_.nroots)) + 0.5)))
    // How many sources send to each aggregator.
    , srcsPerAgg_((size_ + naggs_ - 1) / naggs_)
    // The aggregator this rank sends to.
    , myAgg_(rank_ / srcsPerAgg_ * srcsPerAgg_)
    // [si * nroots + ri] is how much source si in this group sends to root ri.
    , groupSendCounts_(rank_ == myAgg_ ? size_t(req_.nroots) * size_t(srcsPerAgg_) : 0)
    , rootCount_(req_.nroots, 0) {}

  int start() override {
    if (req_.nroots == 0) {
      stage_         = Stage::complete;
      req_.completed = true;
      return MPI_SUCCESS;
    }

    int err = postRootReceives();
    if (err != MPI_SUCCESS) return setError(err);

    // Ensure aggregators know how much data each rank is sending to each root.
    if (rank_ == myAgg_) {
      countReqs_.reserve(srcsPerAgg_ + 1);
      // Receive counts from each member of this aggregation group.
      for (int si = 0; si < srcsPerAgg_ && si + rank_ < size_; ++si) {
        MPI_Request mpiReq;
        err = MPI_Irecv(&groupSendCounts_[size_t(si) * size_t(req_.nroots)], req_.nroots,
                        MPI_INT, si + rank_, req_.aggTag, req_.comm, &mpiReq);
        if (err != MPI_SUCCESS) return setError(err);
        countReqs_.push_back(mpiReq);
      }
    }

    // Send this rank's counts to its aggregator.
    MPI_Request mpiReq;
    err = MPI_Isend(req_.sendcounts, req_.nroots, MPI_INT, myAgg_, req_.aggTag,
                    req_.comm, &mpiReq);
    if (err != MPI_SUCCESS) return setError(err);
    countReqs_.push_back(mpiReq);
    return MPI_SUCCESS;
  }

  int progress() override {
    if (error_ != MPI_SUCCESS || stage_ == Stage::complete) return error_;

    int flag = 0;
    int err  = MPI_SUCCESS;
    if (stage_ == Stage::counts) {
      err = test_all(countReqs_, flag);
      if (err != MPI_SUCCESS) return setError(err);
      if (!flag) return MPI_SUCCESS;
      err = postPayload();
      if (err != MPI_SUCCESS) return setError(err);
    }

    if (stage_ == Stage::payload) {
      err = test_all(payloadReqs_, flag);
      if (err != MPI_SUCCESS) return setError(err);
      if (!flag) return MPI_SUCCESS;
      err = postForwards();
      if (err != MPI_SUCCESS) return setError(err);
    }

    if (stage_ == Stage::forwarding) {
      int rootFlag    = 0;
      int forwardFlag = 0;
      err             = test_all(rootReqs_, rootFlag);
      if (err != MPI_SUCCESS) return setError(err);
      err = test_all(forwardReqs_, forwardFlag);
      if (err != MPI_SUCCESS) return setError(err);
      if (!rootFlag || !forwardFlag) return MPI_SUCCESS;
      finish();
    }
    return MPI_SUCCESS;
  }

  bool isComplete() const override {
    return stage_ == Stage::complete;
  }

  MPI_Comm comm() const override {
    return req_.comm;
  }

  void updateReq(Ialltofewv::Req &req) const override {
    req.completed = isComplete();
  }

 private:
  int getRank() const {
    int rank;
    MPI_Comm_rank(req_.comm, &rank);
    return rank;
  }

  int getSize() const {
    int size;
    MPI_Comm_size(req_.comm, &size);
    return size;
  }

  static size_t getTypeSize(MPI_Datatype type) {
    int size;
    MPI_Type_size(type, &size);
    return size_t(size);
  }

  int setError(int err) {
    error_ = err;
    return err;
  }

  int postRootReceives() {
    if (!isRoot_) return MPI_SUCCESS;

    // An aggregator sends contiguous data, which the root later scatters
    // according to rdispls.
    size_t totalRecvd = 0;
    for (int i = 0; i < size_; ++i) totalRecvd += size_t(req_.recvcounts[i]) * recvSize_;
    rootBuf_ = cache_->get_rootBuf<RecvExecSpace>(totalRecvd);

    // Aggregators send data in rank order, which is also the order used in the
    // root's contiguous receive buffer.
    size_t displ = 0;
    for (int aggSrc = 0; aggSrc < size_; aggSrc += srcsPerAgg_) {
      // Tally the data received from this aggregator.
      int count = 0;
      for (int origSrc = aggSrc; origSrc < aggSrc + srcsPerAgg_ && origSrc < size_; ++origSrc) {
        count += req_.recvcounts[origSrc];
      }
      if (count) {
        MPI_Request mpiReq;
        const int err = MPI_Irecv(rootBuf_.data() + displ, count, req_.recvtype, aggSrc,
                                  req_.rootTag, req_.comm, &mpiReq);
        if (err != MPI_SUCCESS) return err;
        rootReqs_.push_back(mpiReq);
        displ += size_t(count) * recvSize_;
      }
    }
    return MPI_SUCCESS;
  }

  int postPayload() {
    stage_ = Stage::payload;
    if (rank_ == myAgg_) {
      // groupSendCounts now holds counts from every member of this aggregation
      // group. Size the aggregation buffer and record how much goes to each root.
      size_t aggBytes = 0;
      for (int si = 0; si < srcsPerAgg_ && si + rank_ < size_; ++si) {
        for (int ri = 0; ri < req_.nroots; ++ri) {
          const int count = groupSendCounts_[size_t(si) * size_t(req_.nroots) + size_t(ri)];
          rootCount_[ri] += count;
          aggBytes += size_t(count) * sendSize_;
        }
      }
      aggBuf_ = cache_->get_aggBuf<RecvExecSpace>(aggBytes);

      // Senders transmit in root order. Receive in the same order so all data
      // for each root is contiguous in the aggregation buffer.
      size_t displ = 0;
      for (int ri = 0; ri < req_.nroots; ++ri) {
        for (int si = 0; si < srcsPerAgg_ && si + rank_ < size_; ++si) {
          // Receive data for root ri from source si in this group.
          const int count = groupSendCounts_[size_t(si) * size_t(req_.nroots) + size_t(ri)];
          if (count) {
            MPI_Request mpiReq;
            const int err = MPI_Irecv(aggBuf_.data() + displ, count, req_.sendtype, si + rank_,
                                      req_.aggTag, req_.comm, &mpiReq);
            if (err != MPI_SUCCESS) return err;
            payloadReqs_.push_back(mpiReq);
            displ += size_t(count) * sendSize_;
          }
        }
      }
    }

    // Send this rank's data to its aggregator.
    for (int ri = 0; ri < req_.nroots; ++ri) {
      if (req_.sendcounts[ri]) {
        MPI_Request mpiReq;
        const size_t displ = size_t(req_.sdispls[ri]) * sendSize_;
        const int err      = MPI_Isend(reinterpret_cast<const char *>(req_.sendbuf) + displ,
                                       req_.sendcounts[ri], req_.sendtype, myAgg_, req_.aggTag,
                                       req_.comm, &mpiReq);
        if (err != MPI_SUCCESS) return err;
        payloadReqs_.push_back(mpiReq);
      }
    }
    return MPI_SUCCESS;
  }

  int postForwards() {
    stage_ = Stage::forwarding;
    if (rank_ != myAgg_) return MPI_SUCCESS;

    // Forward each root's contiguous data in source-rank order, which is the
    // order expected by the root's contiguous receive buffer.
    size_t displ = 0;
    for (int ri = 0; ri < req_.nroots; ++ri) {
      const size_t count = rootCount_[ri];
      if (count) {
        MPI_Request mpiReq;
        const int err = MPI_Isend(aggBuf_.data() + displ, int(count), req_.sendtype,
                                  req_.roots[ri], req_.rootTag, req_.comm, &mpiReq);
        if (err != MPI_SUCCESS) return err;
        forwardReqs_.push_back(mpiReq);
        displ += count * sendSize_;
      }
    }
    return MPI_SUCCESS;
  }

  void finish() {
    if (isRoot_) {
      // Spread the contiguous root buffer into recvbuf according to rdispls.
      // First set up the source and destination for each block.
      args_         = cache_->get_args<RecvExecSpace>(size_);
      auto args_h   = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, args_);
      size_t srcOff = 0;
      for (int srcRank = 0; srcRank < size_; ++srcRank) {
        const size_t dstOff = size_t(req_.rdispls[srcRank]) * recvSize_;
        const size_t count  = size_t(req_.recvcounts[srcRank]) * recvSize_;
        args_h(srcRank)     = MemcpyArg{reinterpret_cast<char *>(req_.recvbuf) + dstOff,
                                    rootBuf_.data() + srcOff, count};
        srcOff += count;
      }

      // Actually copy the data.
      Kokkos::deep_copy(args_, args_h);
      using Policy = Kokkos::TeamPolicy<RecvExecSpace>;
      Policy policy(size_, Kokkos::AUTO);
      auto args = args_;
      Kokkos::parallel_for(
          "Tpetra::Details::Ialltofewv: apply rdispls to contiguous root buffer", policy,
          ApplyRdispls<decltype(args), typename Policy::member_type>{args});
      Kokkos::fence("Tpetra::Details::Ialltofewv: after apply rdispls to contiguous root buffer");
    }
    stage_         = Stage::complete;
    req_.completed = true;
  }

  Ialltofewv::Req &req_;
  std::shared_ptr<Ialltofewv::Cache::impl> cache_;
  int rank_;
  int size_;
  size_t sendSize_;
  size_t recvSize_;
  bool isRoot_;
  int naggs_;
  int srcsPerAgg_;
  int myAgg_;
  Stage stage_ = Stage::counts;
  int error_   = MPI_SUCCESS;
  std::vector<int> groupSendCounts_;
  std::vector<size_t> rootCount_;
  byte_view_type aggBuf_;
  root_view_type rootBuf_;
  args_view_type args_;
  std::vector<MPI_Request> countReqs_;
  std::vector<MPI_Request> payloadReqs_;
  std::vector<MPI_Request> rootReqs_;
  std::vector<MPI_Request> forwardReqs_;
};

std::vector<std::weak_ptr<Ialltofewv::State>> &activeStates() {
  // Weak ownership lets each Req control its state lifetime while still making
  // every local operation available for cooperative progress.
  static std::vector<std::weak_ptr<Ialltofewv::State>> states;
  return states;
}

int progressAll(MPI_Comm comm) {
  int result   = MPI_SUCCESS;
  auto &states = activeStates();
  for (auto it = states.begin(); it != states.end();) {
    if (auto state = it->lock()) {
      int comparison       = MPI_UNEQUAL;
      const int compareErr = MPI_Comm_compare(comm, state->comm(), &comparison);
      if (result == MPI_SUCCESS && compareErr != MPI_SUCCESS) result = compareErr;
      if (comparison == MPI_IDENT) {
        const int err = state->progress();
        if (result == MPI_SUCCESS && err != MPI_SUCCESS) result = err;
      }
      if (state->isComplete()) {
        it = states.erase(it);
      } else {
        ++it;
      }
    } else {
      it = states.erase(it);
    }
  }
  return result;
}
}  // namespace

int Ialltofewv::start(Req &req) {
  if (!cache_.pimpl) cache_.pimpl = std::make_shared<Cache::impl>();

  if (req.devAccess) {
    req.state = std::make_shared<StateImpl<Kokkos::DefaultExecutionSpace>>(req, cache_.pimpl);
  } else {
    req.state = std::make_shared<StateImpl<Kokkos::DefaultHostExecutionSpace>>(req, cache_.pimpl);
  }

  const int err = req.state->start();
  req.state->updateReq(req);
  if (!req.completed) activeStates().push_back(req.state);
  return err;
}

int Ialltofewv::wait(Req &req) {
  ProfilingRegion pr("alltofewv::wait");
  while (!req.state->isComplete()) {
    const int err = progressAll(req.comm);
    if (err != MPI_SUCCESS) return err;
  }
  req.state->updateReq(req);
  return MPI_SUCCESS;
}

int Ialltofewv::get_status(const Req &req, int *flag, MPI_Status * /*status*/) const {
  const int err = progressAll(req.comm);
  *flag         = req.state->isComplete();
  return err;
}

}  // namespace Tpetra::Details
