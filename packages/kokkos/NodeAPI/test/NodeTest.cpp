#include "NodeTest.hpp"

int NodeTest::N          = 16384;
int NodeTest::numIters   = 10;
int NodeTest::numThreads = -1;
int NodeTest::verbose    = 0;
int NodeTest::cuda_dev   = 0;

namespace {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&NodeTest::N,"Vector length for tests, required >= 2.");
    clp.setOption("num-iters",&NodeTest::numIters,"Number of iterations in TimeTest.");
    clp.setOption("num-threads",&NodeTest::numThreads, "Number of threads. -1 for automatic.");
    clp.setOption("verbose",&NodeTest::verbose, "Node verbosity. Zero is quiet, non-zero is not.");
#ifdef HAVE_KOKKOS_CUDA
    clp.setOption("cuda-device",&NodeTest::cuda_dev,"CUDA device used for testing.");
#endif
  }

}
