#ifndef KOKKOS_TBBNODE_HPP_
#define KOKKOS_TBBNODE_HPP_

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <stdlib.h>

#include <Kokkos_StandardNodeMemoryModel.hpp>

#include <iostream> // debug

template <class WDPin>
struct BlockedRangeWDP {
  mutable WDPin wd;
  BlockedRangeWDP(WDPin &in) : wd(in) {}
  inline void operator()(tbb::blocked_range<int> &rng) const
  { wd(rng.begin(),rng.end()); }
};

template <class WDPin>
struct BlockedRangeWDPReducer {
  WDPin wd;
  BlockedRangeWDPReducer(WDPin &in) : wd(in) {}
  BlockedRangeWDPReducer(BlockedRangeWDPReducer &in, tbb::split) : wd(in.wd) {wd.result = wd.identity();}
  void operator()(tbb::blocked_range<int> &rng)
  { 
    typename WDPin::ReductionType tmpres, tmpi;
    tmpres = wd.result;
    int end = rng.end();
    for (int i=rng.begin(); i != end; ++i) {
      tmpi = wd.generate(i);
      tmpres = wd.reduce( tmpres, tmpi );
    }
    wd.result = tmpres;
  }
  inline void join( const BlockedRangeWDPReducer<WDPin> &other ) {
    wd.result = wd.reduce( wd.result, other.wd.result );
  }
};

class TBBNode : public StandardMemoryModel {
  public:

    TBBNode(int numThreads=0) {
      if (numThreads >= 1) {
        tsi_.initialize(numThreads);
      }
      else {
        tsi_.initialize(tbb::task_scheduler_init::automatic);
      }
    }

    ~TBBNode() {}

    template <class WDP>
    void execute1D(int length, WDP wd) {
      BlockedRangeWDP<WDP> tbb_wd(wd);
      tbb::parallel_for(tbb::blocked_range<int>(0,length), tbb_wd, tbb::auto_partitioner()); 
    }

    template <class WDP>
    void reduce1D(int length, WDP &wd) {
      BlockedRangeWDPReducer<WDP> tbb_wd(wd);
      tbb::parallel_reduce(tbb::blocked_range<int>(0,length), tbb_wd, tbb::auto_partitioner());
      wd.result = tbb_wd.wd.result;  // have to put result from final tbb_wd into orginal wd
    }

  private:
    static tbb::task_scheduler_init tsi_;

};

#endif
