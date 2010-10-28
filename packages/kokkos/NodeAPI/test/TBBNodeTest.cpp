#include "NodeTest.hpp"

#ifdef HAVE_KOKKOS_TBB

#include "Kokkos_TBBNode.hpp"
#include <Teuchos_ScalarTraits.hpp>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace {

  using Kokkos::TBBNode;
  RCP<TBBNode> tbbNode_;

  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    return tbbNode_;
  }

  template <>
  void initNode<TBBNode>() {
    Teuchos::ParameterList plist;
    if (NodeTest::numThreads != -1) {
      plist.set<int>("Num Threads",NodeTest::numThreads);
    }
    plist.set<int>("Verbose",NodeTest::verbose);
    tbbNode_ = rcp(new TBBNode(plist));
  }

  template <class T>
  class TBBNodeTestInit {
  private:
    T * const a_;
  public:
    TBBNodeTestInit(T* a) : a_(a) {}
    void operator()(const tbb::blocked_range<int> &r) const {
      T * const my_a = a_;
      for (int i=r.begin(); i != r.end(); ++i) {
        my_a[i] = Teuchos::ScalarTraits<T>::one();
      }
    }
  };

  template <class T>
  class TBBNodeTestSum {
  private:
    const T * const a_;
  public:
    T sum;
    TBBNodeTestSum(const T* a) : a_(a), sum(Teuchos::ScalarTraits<T>::zero()) {}
    TBBNodeTestSum(TBBNodeTestSum<T> &other, tbb::split) : a_(other.a_), sum(Teuchos::ScalarTraits<T>::zero()) {}
    void join(const TBBNodeTestSum<T> &other) { sum += other.sum; }
    void operator()(const tbb::blocked_range<int> &r) {
      const T* const my_a = a_;
      for (int i=r.begin(); i != r.end(); ++i) {
        sum += my_a[i];
      }
    }
  };

  template <>
  std::pair<double,double> nativeTimings<float,TBBNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,TBBNode init"), sTime("float,TBBNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        tbb::parallel_for( tbb::blocked_range<int>(0,N), TBBNodeTestInit<float>(buff.getRawPtr()), tbb::auto_partitioner() );
      }
    }
    float sum = 0.0f;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TBBNodeTestSum<float> op(buff.getRawPtr());
        tbb::parallel_reduce( tbb::blocked_range<int>(0,N), op, tbb::auto_partitioner() );
        sum = op.sum;
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,TBBNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,TBBNode init"), sTime("int,TBBNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        tbb::parallel_for( tbb::blocked_range<int>(0,N), TBBNodeTestInit<int>(buff.getRawPtr()), tbb::auto_partitioner() );
      }
    }
    int sum = 0;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TBBNodeTestSum<int> op(buff.getRawPtr());
        tbb::parallel_reduce( tbb::blocked_range<int>(0,N), op, tbb::auto_partitioner() );
        sum = op.sum;
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(TBBNode)

}

#endif
