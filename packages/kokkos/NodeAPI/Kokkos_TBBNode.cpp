#include "Kokkos_TBBNode.hpp"
#include <Teuchos_ParameterList.hpp>
#include <iostream>

// tbb::task_scheduler_init Kokkos::TBBNode::tsi_(tbb::task_scheduler_init::deferred);

namespace Kokkos {

  TBBNode::TBBNode(Teuchos::ParameterList &pl) : alreadyInit_(false), tsi_(tbb::task_scheduler_init::deferred) {
    int numThreads = pl.get<int>("Num Threads",-1);
    int verbose = pl.get<int>("Verbose",0);
    if (numThreads >= 0) {
      if (verbose) {
        std::cout << "TBBNode initializing with numThreads == " << numThreads << std::endl;
      }
      init(numThreads);
    }
  }

  void TBBNode::init(int numThreads) {
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

  TBBNode::~TBBNode() {}

}
