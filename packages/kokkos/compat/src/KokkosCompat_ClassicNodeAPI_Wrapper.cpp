#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Kokkos {
  namespace Compat {
    KokkosThreadsWrapperNode::KokkosThreadsWrapperNode(Teuchos::ParameterList &pl) {
      ParameterList params = getDefaultParameters();
      params.setParameters(pl);
      const int curNumThreads = params.get<int>("Num Threads");
      int verboseInt = params.get<int>("Verbose");
      bool verbose = (verboseInt != 0);
      if (verbose) {
        std::cout << "ThreadsWrapperNode initializing with \"Num Threads\" = "
                  << curNumThreads << std::endl;
      }
      init (curNumThreads);
    }
    KokkosThreadsWrapperNode::KokkosThreadsWrapperNode() {
      init(1);
    };

    KokkosThreadsWrapperNode::~KokkosThreadsWrapperNode() {
      Kokkos::Threads::finalize();
    }
    ParameterList KokkosThreadsWrapperNode::getDefaultParameters() {
      ParameterList params;
      params.set("Verbose",     0);
      params.set("Num Threads", 1);
      return params;
    }

    void KokkosThreadsWrapperNode::init(int num_threads) {
      Kokkos::Threads::initialize(num_threads,1);
    }
  }
}
