#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "TestOps.hpp"

namespace {

  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Kokkos::DefaultNode;
  using Kokkos::ReadyBufferHelper;
  using Teuchos::tuple;
  using Teuchos::ArrayRCP;
  using Teuchos::Tuple;

  int N = 1000;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( NodeAPI, DefaultNodeTest )
  {
    Time tAlloc("Alloc Time"), tInit("Init Op"), tSum("Sum Op"), tFree("Free Time");
    typedef DefaultNode::DefaultNodeType NODE;
    typedef ArrayRCP<const char>  cbuf;
    typedef ArrayRCP<      char> ncbuf;
    Teuchos::ArrayRCP<int> x;
    NODE &node = DefaultNode::getDefaultNode();
    ReadyBufferHelper<NODE> rbh(node);
    out << "Default Node Type: " << Teuchos::typeName(node) << std::endl;
    int result;
    {
      TimeMonitor localTimer(tAlloc);
      x = node.allocBuffer<int>(N);
    }
    // set x[i] = 1, i=0:N-1
    {
      TimeMonitor localTimer(tInit);
      InitOp<int> wdp;
      rbh.begin();
      wdp.x = rbh.addNonConstBuffer(x);
      rbh.end();
      node.parallel_for(0,N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      TimeMonitor localTimer(tSum);
      SumOp<int> wdp;
      rbh.begin();
      wdp.x = rbh.addConstBuffer<int>(x);
      rbh.end();
      result = node.parallel_reduce(0,N,wdp);
    }
    int expectedResult = (int)(N);
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

}
