#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "TestOps.hpp"

#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif

#ifdef HAVE_KOKKOS_CUDA
#include "Kokkos_CUDANode.hpp"
#endif

namespace KokkosTest {
  SerialNode serialnode;
#ifdef HAVE_KOKKOS_TBB
  TBBNode tbbnode;
#endif
#ifdef HAVE_KOKKOS_TBB
  CUDANode cudanode;
#endif
}

namespace {

  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Kokkos::SerialNode;
#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
#endif
#ifdef HAVE_KOKKOS_CUDA
  using Kokkos::CUDANode;
#endif

  int N = 100;
  template <class NODE>
  NODE & getNode() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  template <>
  SerialNode & getNode<SerialNode() {
    return KokkosTest::serialnode;
  }

#ifdef HAVE_KOKKOS_TBB
  template <>
  TBBNode & getNode<TBBNode() {
    return KokkosTest::tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_TBB
  template <>
  CUDANode & getNode<CUDANode() {
    return KokkosTest::cudanode;
  }
#endif

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
    Kokkos::SerialNode KokkosTest::serialnode;
#ifdef HAVE_KOKKOS_TBB
    {
      int tbb_nT = 0;
      clp.setOption("tbb-num-threads",&tbb_nT,"Number of TBB threads: 0 for automatic.");
      Kokkos::TBBNode KokkosTest::tbbnode(tbb_nT);
    }
#endif
#ifdef HAVE_KOKKOS_TBB
    {
      int cuda_dev = 0;
      int cuda_nT  = 64;
      int cuda_nB  = 64;
      int cuda_verb = 0;
      clp.setOption("cuda-device",&cuda_dev,"CUDA device used for testing.");
      clp.setOption("cuda-num-threads",&cuda_nT,"Number of CUDA threads.");
      clp.setOption("cuda-num-blocks",&cuda_nB,"Number of CUDA blocks.");
      clp.setOption("cuda-verbose",&cuda_verb,"CUDA verbosity level.");
      Kokkos::CUDANode KokkosTest::cudanode(cuda_dev,cuda_nB,cuda_nT,cuda_verb);
    }
#endif
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( NodeAPI, SumTest, SCALAR, NODE )
  {
    Time tAlloc("Alloc Time"), tInit("Init Op"), tSum("Sum Op"), tFree("Free Time");
    typedef Node::template buffer<SCALAR>::buffer_t x;
    NODE &node = getNode<NODE>();
    SCALAR result;
    {
      TimeMonitor localTimer(tAlloc);
      x = node.template allocBuffer<SCALAR>(N);
    }
    // set x[i] = i, i=0:N-1
    {
      InitOp<SCALAR,NODE> wdp;
      wdp.x = x;
      wdp.n = N;
      node.parallel_for(N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      SumOp<SCALAR,NODE> wdp;
      wdp.x = x;
      result = node.parallel_reduce(N,wdp);
    }
    SCALAR expectedResult = (SCALAR)(N) * (SCALAR)(N-1) / (SCALAR)(2);
    TEST_EQUALITY(result, expectedResult);
    {
      TimeMonitor localTimer(tFree);
      node.template freeBuffer<T>(x);
    }
  }

  // 
  // INSTANTIATIONS
  //

#define SERIAL_INSTANT(SCALAR) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, SumTest, SCALAR, SerialNode )

#ifdef HAVE_KOKKOS_TBB
#define TBB_INSTANT(SCALAR) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, SumTest, SCALAR, TBBNode )
#else
#define TBB_INSTANT(SCALAR) 
#endif

#ifdef HAVE_KOKKOS_CUDA
#define CUDA_INSTANT(SCALAR) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, SumTest, SCALAR, CUDANode )
#else
#define CUDA_INSTANT(SCALAR) 
#endif

#define UNIT_TEST_GROUP_SCALAR(SCALAR) \
  SERIAL_INSTANT(SCALAR) \
  TBB_INSTANT(SCALAR) \
  CUDA_INSTANT(SCALAR)

  UNIT_TEST_GROUP_SCALAR(int)
  UNIT_TEST_GROUP_SCALAR(float)
}
