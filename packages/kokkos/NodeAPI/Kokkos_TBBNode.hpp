#ifndef KOKKOS_TBBNODE_HPP_
#define KOKKOS_TBBNODE_HPP_

#include "Kokkos_StandardNodeMemoryModel.hpp"
#include "Kokkos_NodeHelpers.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>

#include <Teuchos_ParameterList.hpp>

#include <stdlib.h>

namespace Kokkos {

template <class WDPin>
struct BlockedRangeWDP {
  mutable WDPin &wd;
  inline BlockedRangeWDP(WDPin &in) : wd(in) {}
  inline void operator()(tbb::blocked_range<int> &rng) const { 
    for (int i=rng.begin(); i != rng.end(); ++i) wd.execute(i);
  }
};

template <class WDPin>
struct BlockedRangeWDPReducer {
  WDPin &wd;
  typename WDPin::ReductionType result;
  BlockedRangeWDPReducer(WDPin &in) : wd(in), result(WDPin::identity()) {}
  BlockedRangeWDPReducer(BlockedRangeWDPReducer &in, tbb::split) : wd(in.wd) {result = wd.identity();}
  void operator()(tbb::blocked_range<int> &rng)
  { 
    typename WDPin::ReductionType tmpi;
    int end = rng.end();
    for (int i=rng.begin(); i != end; ++i) {
      tmpi = wd.generate(i);
      result = wd.reduce( result, tmpi );
    }
  }
  inline void join( const BlockedRangeWDPReducer<WDPin> &other ) {
    result = wd.reduce( result, other.result );
  }
};

class TBBNode : public StandardNodeMemoryModel {
  public:

    TBBNode(Teuchos::ParameterList &pl) : alreadyInit_(false), tsi_(tbb::task_scheduler_init::deferred) {
      int numThreads = pl.get<int>("Num Threads",-1);
      if (numThreads >= 0) {
        init(numThreads);
      }
    }

    void init(int numThreads) {
      if (alreadyInit_) {
        tsi_.terminate();
      }
      // 
      if (numThreads >= 1) {
        tsi_.initialize(numThreads);
      }
      else {
        tsi_.initialize(tbb::task_scheduler_init::automatic);
      }
    }

    ~TBBNode() {}

    template <class WDP>
    static void parallel_for(int begin, int end, WDP wd) {
      BlockedRangeWDP<WDP> tbb_wd(wd);
      tbb::parallel_for(tbb::blocked_range<int>(begin,end), tbb_wd, tbb::auto_partitioner()); 
    }

    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd) {
      BlockedRangeWDPReducer<WDP> tbb_wd(wd);
      tbb::parallel_reduce(tbb::blocked_range<int>(begin,end), tbb_wd, tbb::auto_partitioner());
      return tbb_wd.result;
    }

  private:
    bool alreadyInit_;
    tbb::task_scheduler_init tsi_;

};

template <> class ArrayOfViewsHelper<TBBNode> : public ArrayOfViewsHelperTrivialImpl<TBBNode> {};

}

#endif
