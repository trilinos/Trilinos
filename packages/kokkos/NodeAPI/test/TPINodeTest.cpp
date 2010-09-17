#include "NodeTest.hpp"

#ifdef HAVE_KOKKOS_THREADPOOL

#include "Kokkos_TPINode.hpp"

namespace {

  using Kokkos::TPINode;
  RCP<TPINode> tpiNode_;

  template <>
  RCP<TPINode> getNode<TPINode>() {
    return tpiNode_;
  }

  template <>
  void initNode<TPINode>() {
    Teuchos::ParameterList plist;
    if (NodeTest::numThreads != -1) {
      plist.set<int>("Num Threads",NodeTest::numThreads);
    }
    plist.set<int>("Verbose",NodeTest::verbose);
    tpiNode_ = rcp(new TPINode(plist));
  }

  template <class T>
  struct TPIInit {
    T * x;
    unsigned int N;
    static void work( TPI_Work * work ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
      unsigned int begin;
      const unsigned int max = ( w->N + work->count - 1 ) / work->count ;
      const unsigned int end = w->N - max * ( work->count - ( work->rank + 1 ) );
      if ( work->rank ) {
        begin  = end - max ;
      }
      else {
        begin  = 0 ;
      }
      T * const my_x = w->x;
      for (int i=begin; i != end; ++i) {my_x[i] = Teuchos::ScalarTraits<T>::one();}
    }
  };

  template <class T>
  struct TPISum {
    const T * x;
    unsigned int N;
    // reduction
    static void work( TPI_Work * work ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
      T * const dst = (T *) work->reduce;
      unsigned int begin;
      const unsigned int max = ( w->N + work->count - 1 ) / work->count ;
      const unsigned int end = w->N - max * ( work->count - ( work->rank + 1 ) );
      if ( work->rank ) {
        begin  = end - max ;
      }
      else {
        begin  = 0 ;
      }
      T * const my_x = w->x;
      for (int i=begin; i != end; ++i) {*dst += my_x[i];}
    }
    // initialization
    static void init( TPI_Work * work ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
      T * const dst = (T *) work->reduce ;
      (*dst) = Teuchos::ScalarTraits<T>::zero();
    }
    // combination
    static void join( TPI_Work * work , const void * arg_src ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
            T * const dst = (T *) work->reduce ;
      const T * const src = (const T *) arg_src ;
      (*dst) += (*src);
    }
  };

  template <>
  std::pair<double,double> nativeTimings<float,TPINode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,TPINode init"), sTime("float,TPINode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        TPIInit<float> data;
        data.x = buff.getRawPtr();
        data.N = N;
        TPI_Run_threads( &TPIInit<float>::work, &data, 0 );
      }
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TPISum<float> data;
        data.x = buff.getRawPtr();
        data.N = N;
        sum = 0.0f;
        TPI_Run_threads_reduce( &TPISum<float>::work, &data, &TPISum<float>::join, &TPISum<float>::init, sizeof(sum), &sum );
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,TPINode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,TPINode init"), sTime("int,TPINode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        TPIInit<int> data;
        data.x = buff.getRawPtr();
        data.N = N;
        TPI_Run_threads( &TPIInit<int>::work, &data, 0 );
      }
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TPISum<int> data;
        data.x = buff.getRawPtr();
        data.N = N;
        sum = 0;
        TPI_Run_threads_reduce( &TPISum<int>::work, &data, &TPISum<int>::join, &TPISum<int>::init, sizeof(sum), &sum );
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(TPINode);

}

#endif
