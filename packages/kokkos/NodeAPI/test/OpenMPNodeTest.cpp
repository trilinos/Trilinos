#include "NodeTest.hpp"

#include "Kokkos_OpenMPNode.hpp"
#include <Teuchos_ScalarTraits.hpp>
#include <omp.h>

namespace {

  using Kokkos::OpenMPNode;
  RCP<OpenMPNode> ompNode_;

  template <>
  RCP<OpenMPNode> getNode<OpenMPNode>() {
    return ompNode_;
  }

  template <>
  void initNode<OpenMPNode>() {
    Teuchos::ParameterList plist;

    int useNumThreads = 1;
    if (NodeTest::numThreads != -1) {
      useNumThreads = NodeTest::numThreads;
    }
    plist.set<int>("Num Threads",useNumThreads);
    plist.set<int>("Verbose",NodeTest::verbose);
    ompNode_ = rcp(new OpenMPNode(plist));
  }

  template <>
  std::pair<double,double> nativeTimings<float,OpenMPNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,OpenMPNode init"), sTime("float,OpenMPNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        float *bptr = buff.getRawPtr();
#pragma omp parallel for default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          bptr[i] = 1.0f;
        }
      }
    }
    float sum = 0.0f;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        const float *bptr = buff.getRawPtr();
        sum = 0.0f;
#pragma omp parallel for reduction (+:sum) default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          sum += bptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,OpenMPNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,OpenMPNode init"), sTime("int,OpenMPNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        int *bptr = buff.getRawPtr();
#pragma omp parallel for default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          bptr[i] = 1.0f;
        }
      }
    }
    int sum = 0;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        const int *bptr = buff.getRawPtr();
        sum = 0;
#pragma omp parallel for reduction (+:sum) default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          sum += bptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(OpenMPNode)

}
