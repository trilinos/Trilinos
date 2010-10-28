#include "Kokkos_SerialNode.hpp"
#include "NodeTest.hpp"

namespace {

  using Kokkos::SerialNode;
  RCP<SerialNode> serialNode_;

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    return serialNode_;
  }

  template <>
  void initNode<SerialNode>() {
    Teuchos::ParameterList plist;
    serialNode_ = rcp(new SerialNode(plist));
  }

  template <>
  std::pair<double,double> nativeTimings<float,SerialNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,SerialNode init"), sTime("float,SerialNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      float *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        for (int i=0; i < N; ++i) {
          ptr[i] = 1;
        }
      }
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      const float *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        sum = ptr[0];
        for (int i=1; i < N; ++i) {
          sum += ptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,SerialNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,SerialNode init"), sTime("int,SerialNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      int *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        for (int i=0; i < N; ++i) {
          ptr[i] = 1;
        }
      }
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      const int *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        sum = ptr[0];
        for (int i=1; i < N; ++i) {
          sum += ptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(SerialNode)

}
