#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Tuple.hpp>

#include <functional>
#include <algorithm>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "TestOps.hpp"

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOS_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

#define TOCBUF(arr) Teuchos::arcp_reinterpret_cast<const char>(arr)
#define TONCBUF(arr) Teuchos::arcp_reinterpret_cast<char>(arr)

namespace {

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::Time;
  using Teuchos::null;
  using Teuchos::TimeMonitor;
  using Kokkos::SerialNode;
  using Kokkos::ReadyBufferHelper;
  using Kokkos::ArrayOfViewsHelper;
  using Teuchos::tuple;

  int N = 100;
  int NumIters = 1000;
  template <class NODE>
  RCP<NODE> getNode() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  RCP<SerialNode> serialnode;
  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    ParameterList pl;
    if (serialnode == null) {
      serialnode = rcp(new SerialNode(pl));
    }
    return serialnode;
  }

#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
  int tbb_nT = 0;
  RCP<TBBNode> tbbnode;
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      ParameterList pl;
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
  using Kokkos::TPINode;
  int tpi_nT = 1;
  RCP<TPINode> tpinode;
  template <>
  RCP<TPINode> getNode<TPINode>() {
    if (tpinode == null) {
      ParameterList pl;
      pl.set("Num Threads",1);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }
#endif

#ifdef HAVE_KOKKOS_THRUST
  using Kokkos::ThrustGPUNode;
  int cuda_dev = 0;
  int cuda_verb = 0;
  RCP<ThrustGPUNode> thrustnode;
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      ParameterList pl; 
      pl.set("Device Number",cuda_dev);
      pl.set("Verbosity",cuda_verb);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
    clp.setOption("num-iters",&NumIters,"Number of iterations in TimeTest.");
    if (N < 2) N = 2;
#ifdef HAVE_KOKKOS_TBB
    {
      clp.setOption("tbb-num-threads",&tbb_nT,"Number of TBB threads: 0 for automatic.");
    }
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
    {
      clp.setOption("tpi-num-threads",&tpi_nT,"Number of TPI threads.");
    }
#endif
#ifdef HAVE_KOKKOS_CUDA
    {
      clp.setOption("cuda-device",&cuda_dev,"CUDA device used for testing.");
      clp.setOption("cuda-verbose",&cuda_verb,"CUDA verbosity level.");
    }
#endif
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST( AAAAA, Init_Nodes )
  {
#ifdef HAVE_KOKKOS_TBB
    out << "Initializing TBB node to " << tbb_nT << " threads." << std::endl;
    getNode<TBBNode>()->init(tbb_nT);
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
    out << "Initializing TPI node to " << tpi_nT << " threads." << std::endl;
    getNode<TPINode>()->init(tpi_nT);
#endif
#ifdef HAVE_KOKKOS_CUDA
    out << "Initializing CUDA device " << cuda_dev << std::endl;
    getNode<ThrustGPUNode>();
#endif
    TEST_EQUALITY(0,0);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NodeAPI, MemoryInitTest, NODE )
  {
    out << "Testing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    ArrayRCP<int> x, y;
    RCP<NODE> node = getNode<NODE>();
    ReadyBufferHelper<NODE> rbh(node);
    const unsigned int N = 4096;
    // const unsigned int N = 1000000;
    const unsigned int Xoff  = 1024;
    const unsigned int Xsize = N + 2*Xoff;  // 6144
    const unsigned int Yoffs[5] = {0, N, 2*N, 3*N, 4*N};
    // allocate two compute buffers: 
    x = node->template allocBuffer<int>(Xsize);
    y = node->template allocBuffer<int>(Yoffs[4]);
    // set x = {?  fours  ?} using copyTo
    {
      Array<int> fours(N);
      std::fill(fours.begin(), fours.end(), 4);
      node->template copyToBuffer<int>(N,fours(),x+Xoff);
    }
    // set y = {? ? ? fours} using copyBuffers
    node->template copyBuffers<int>(N,x+Xoff,y+Yoffs[3]);
    // set y = {? ? threes fours} using local viewBuffer(WriteOnly)
    {
      ArrayRCP<int> view = node->viewBufferNonConst(Kokkos::WriteOnly, N, y+Yoffs[2]);
      TEST_EQUALITY_CONST( view.size(), N );
      std::fill(view.begin(), view.end(), 3);
    }
    // set y = {ones ones threes fours} using init
    {
      InitOp<int> iop;
      rbh.begin();
      iop.x = rbh.template addNonConstBuffer<int>(y);
      rbh.end();
      const int beg = Yoffs[0],
                end = Yoffs[2];
      node->parallel_for(beg,end,iop);
    }
    // set y = {ones twos threes fours} using viewBuffer(ReadWrite) and local arithmetic
    // over-size the view, to ensure that we don't affect any more changes than necessary
    {
      ArrayRCP<int> view = node->viewBufferNonConst(Kokkos::ReadWrite, 2*N, y+Yoffs[1]);
      // add 1 to all elements of view
      std::transform(view.begin(), view.begin()+N, view.begin(), std::bind2nd(std::plus<int>(), 1));
    }
    // check the result
    { 
      Array<int> expected(N);
      ArrayRCP<const int> yview = node->template viewBuffer<int>(4*N,y);
      // y should be {ones  twos  threes  fours}
      // check ones
      std::fill(expected.begin(), expected.end(), 1);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[0],N) ); 
      // check twos
      std::fill(expected.begin(), expected.end(), 2);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[1],N) ); 
      // check threes
      std::fill(expected.begin(), expected.end(), 3);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[2],N) ); 
      // check fours
      std::fill(expected.begin(), expected.end(), 4);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[3],N) ); 
    }
    // free the allocations
    x = Teuchos::null;
    y = Teuchos::null;
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NodeAPI, ArrayOfViewsHelper, NODE )
  {
    RCP<NODE> node = getNode<NODE>();
    const int N = 5;
    ArrayRCP< ArrayRCP<int> > bufs, views;
    bufs = arcp< ArrayRCP<int> >(N);
    for (size_t i=0; i < N; ++i) {
      bufs[i] = node->template allocBuffer<int>(i+1);
    }
    // get views, set data, delete them
    views = ArrayOfViewsHelper<NODE>::template getArrayOfNonConstViews<int>(node, Kokkos::WriteOnly, bufs);
    TEST_EQUALITY_CONST( (unsigned int)views.size(), N );
    for (size_t i=0; i < N; ++i) {
      TEST_EQUALITY( (unsigned int)views[i].size(), i+1 );
      std::fill( views[i].begin(), views[i].end(), i+1 );
    }
    views = Teuchos::null; 
    // get views, verify data, delete them
    views = ArrayOfViewsHelper<NODE>::template getArrayOfNonConstViews<int>(node, Kokkos::ReadWrite, bufs);
    TEST_EQUALITY_CONST( (unsigned int)views.size(), N );
    for (size_t i=0; i < N; ++i) {
      TEST_EQUALITY( (unsigned int)views[i].size(), i+1 );
      Array<int> exp(i+1,i+1);      
      TEST_COMPARE_ARRAYS( (views[i])(), exp() ); 
    }
    views = Teuchos::null; 
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( NodeAPI, SumTest, SCALAR, NODE )
  {
    out << "Testing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    Time tAlloc("Alloc Time"), tInit("Init Op"), tSum("Sum Op"), tFree("Free Time");
    Teuchos::ArrayRCP<SCALAR> x;
    RCP<NODE> node = getNode<NODE>();
    ReadyBufferHelper<NODE> rbh(node);
    SCALAR result;
    {
      TimeMonitor localTimer(tAlloc);
      x = node->template allocBuffer<SCALAR>(N);
    }
    // set x[i] = 1, i=0:N-1
    {
      TimeMonitor localTimer(tInit);
      InitOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addNonConstBuffer<SCALAR>(x);
      rbh.end();
      node->parallel_for(0,N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      TimeMonitor localTimer(tSum);
      SumOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addConstBuffer<SCALAR>(x);
      rbh.end();
      result = node->parallel_reduce(0,N,wdp);
    }
    SCALAR expectedResult = (SCALAR)(N);
    TEST_EQUALITY(result, expectedResult);
    // compute sum x[i], i=1:N-2
    {
      TimeMonitor localTimer(tSum);
      SumOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addConstBuffer<SCALAR>(x);
      rbh.end();
      result = node->parallel_reduce(1,N-1,wdp);
    }
    expectedResult = (SCALAR)(N-2);
    TEST_EQUALITY(result, expectedResult);
    {
      TimeMonitor localTimer(tFree);
      x = Teuchos::null;
    }
    out << "allocBuffer Time: " << tAlloc.totalElapsedTime() << std::endl;
    out << "InitOp Time: " << tInit.totalElapsedTime() << std::endl;
    out << "SumOp Time: " << tSum.totalElapsedTime() << std::endl;
    out << "freeBuffer Time: " << tFree.totalElapsedTime() << std::endl;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NodeAPI, TimeTest, NODE )
  {
    out << "Testing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    Time tNoop("Null Op");
    RCP<NODE> node = getNode<NODE>();
    NullOp noop;
    int red = 0;
    {
      TimeMonitor localTimer(tNoop);
      for (int i=0; i<NumIters; ++i) {
        red = node->parallel_reduce(0,1,noop);
      }
    }
    TEST_EQUALITY_CONST(red,0);
    out << "NullOp Time: " << tNoop.totalElapsedTime() << std::endl;
    out << "    average: " << (int)(tNoop.totalElapsedTime() / (double)(NumIters) * 1000000000.0) << " ns" << std::endl;
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

#ifdef HAVE_KOKKOS_THREADPOOL
#define TPI_INSTANT(SCALAR) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, SumTest, SCALAR, TPINode )
#else
#define TPI_INSTANT(SCALAR) 
#endif

#ifdef HAVE_KOKKOS_CUDA
#define CUDA_INSTANT(SCALAR) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, SumTest, SCALAR, ThrustGPUNode )
#else
#define CUDA_INSTANT(SCALAR) 
#endif

#define UNIT_TEST_GROUP_SCALAR(SCALAR) \
  SERIAL_INSTANT(SCALAR) \
  TBB_INSTANT(SCALAR) \
  TPI_INSTANT(SCALAR) \
  CUDA_INSTANT(SCALAR)

  UNIT_TEST_GROUP_SCALAR(int)
  UNIT_TEST_GROUP_SCALAR(float)
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, TimeTest, SerialNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, MemoryInitTest, SerialNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, ArrayOfViewsHelper, SerialNode )
#ifdef HAVE_KOKKOS_TBB
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, TimeTest, TBBNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, MemoryInitTest, TBBNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, ArrayOfViewsHelper, TBBNode )
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, TimeTest, TPINode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, MemoryInitTest, TPINode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, ArrayOfViewsHelper, TPINode )
#endif
#ifdef HAVE_KOKKOS_CUDA
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, TimeTest, ThrustGPUNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, MemoryInitTest, ThrustGPUNode )
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, ArrayOfViewsHelper, ThrustGPUNode )
#endif

}
