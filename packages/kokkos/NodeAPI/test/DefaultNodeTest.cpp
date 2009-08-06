#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_ConfigDefs.hpp"
#include "TestOps.hpp"

namespace {

  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Kokkos::DefaultNode;

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
    Teuchos::ArrayRCP<int> x;
    NODE &node = DefaultNode::getDefaultNode();
    out << "Default Node Type: " << Teuchos::typeName(node) << std::endl;
    int result;
    {
      TimeMonitor localTimer(tAlloc);
      x = node.allocBuffer<int>(N);
    }
    // set x[i] = 1, i=0:N-1
    {
      TimeMonitor localTimer(tInit);
      InitOp<int,NODE> wdp;
      wdp.x = x;
      node.parallel_for(0,N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      TimeMonitor localTimer(tSum);
      SumOp<int,NODE> wdp;
      wdp.x = x;
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
