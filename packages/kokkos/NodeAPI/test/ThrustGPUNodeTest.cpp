#include "NodeTest.hpp"

#ifdef HAVE_KOKKOS_THRUST

#include "Kokkos_ThrustGPUNode.hpp"
#include "Kokkos_CUDA_util_inline_runtime.h"
#include <thrust/device_vector.h>

void thrust_float_init(thrust::device_vector<float> &buff);
float thrust_float_sum(const thrust::device_vector<float> &buff);
void thrust_int_init(thrust::device_vector<int> &buff);
int thrust_int_sum(const thrust::device_vector<int> &buff);

namespace {

  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustNode_;

  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    return thrustNode_;
  }

  template <>
  void initNode<ThrustGPUNode>() {
    Teuchos::ParameterList plist;
    plist.set<int>("Device Number",NodeTest::cuda_dev);
    plist.set<int>("Verbose",NodeTest::verbose);
    thrustNode_ = rcp(new ThrustGPUNode(plist));
  }

  template <>
  std::pair<double,double> nativeTimings<float,ThrustGPUNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,ThrustGPUNode init"), sTime("float,ThrustGPUNode sum");
    thrust::device_vector<float> buff(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        thrust_float_init(buff);
      }
      Kokkos::cutilSafeThreadSync();
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        sum = thrust_float_sum(buff);
      }
      Kokkos::cutilSafeThreadSync();
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,ThrustGPUNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,ThrustGPUNode init"), sTime("int,ThrustGPUNode sum");
    thrust::device_vector<int> buff(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        thrust_int_init(buff);
      }
      Kokkos::cutilSafeThreadSync();
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        sum = thrust_int_sum(buff);
      }
      Kokkos::cutilSafeThreadSync();
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(ThrustGPUNode);
}

#endif
