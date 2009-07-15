#ifndef KOKKOS_TBBNODE_HPP_
#define KOKKOS_TBBNODE_HPP_

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <stdlib.h>

#include <Kokkos_StandardNodeMemoryModel.hpp>

#include <iostream> // debug

namespace Kokkos {

template <class WDPin>
struct BlockedRangeWDP {
  mutable WDPin wd;
  BlockedRangeWDP(WDPin &in) : wd(in) {}
  inline void operator()(tbb::blocked_range<int> &rng) const { 
    for (int i=rng.begin(); i != rng.end(); ++i) wd.execute(i);
  }
};

template <class WDPin>
struct BlockedRangeWDPReducer {
  WDPin wd;
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

class TBBNode : public StandardMemoryModel {
  public:

    TBBNode(int numThreads=0) : alreadyInit_(false) {
      init(numThreads);
    }

    void init(int numThreads) {
      if (alreadyInit_) {
        tsi_.terminate();
      }
      if (numThreads >= 1) {
        tsi_.initialize(numThreads);
      }
      else {
        tsi_.initialize(tbb::task_scheduler_init::automatic);
      }
    }

    ~TBBNode() {}

    template <class WDP>
    void parallel_for(int begin, int end, WDP wd) {
      BlockedRangeWDP<WDP> tbb_wd(wd);
      tbb::parallel_for(tbb::blocked_range<int>(begin,end), tbb_wd, tbb::auto_partitioner()); 
    }

    template <class WDP>
    typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd) {
      BlockedRangeWDPReducer<WDP> tbb_wd(wd);
      tbb::parallel_reduce(tbb::blocked_range<int>(begin,end), tbb_wd, tbb::auto_partitioner());
      return tbb_wd.result;
    }

  private:
    bool alreadyInit_;
    static tbb::task_scheduler_init tsi_;

};

}

#endif
