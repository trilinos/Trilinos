/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

//#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemv_Decl.hpp"
#include "KokkosBatched_Gemv_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched::Experimental;

namespace Test {

  template<typename T>
  struct ParamTag { 
    typedef T trans;
  };
 
  template<typename DeviceType,
           typename ViewType,
           typename ScalarType,
           typename ParamTagType, 
           typename AlgoTagType>
  struct Functor_TestBatchedSerialGemv {
    ViewType _a, _b, _c;
    
    ScalarType _alpha, _beta;
    
    KOKKOS_INLINE_FUNCTION
    Functor_TestBatchedSerialGemv(const ScalarType alpha, 
            const ViewType &a,
            const ViewType &b,
            const ScalarType beta,
            const ViewType &c)
      : _a(a), _b(b), _c(c), _alpha(alpha), _beta(beta) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const ParamTagType &, const int k) const {
      auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
      auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), 0);
      auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), 0);
      
      SerialGemv<typename ParamTagType::trans,
        AlgoTagType>::
        invoke(_alpha, aa, bb, _beta, cc);
    }
    
    inline
    void run() {
      Kokkos::RangePolicy<DeviceType,ParamTagType> policy(0, _c.extent(0));
      Kokkos::parallel_for(policy, *this);            
    }
  };
    
  template<typename DeviceType,
           typename ViewType,
           typename ScalarType,
           typename ParamTagType, 
           typename AlgoTagType>
  void impl_test_batched_gemv(const int N, const int BlkSize) {
    typedef typename ViewType::value_type value_type;
    typedef Kokkos::Details::ArithTraits<value_type> ats;

    /// randomized input testing views
    ScalarType alpha = 1.5, beta = 3.0;

    ViewType
      a0("a0", N, BlkSize, BlkSize), a1("a1", N, BlkSize, BlkSize),
      b0("b0", N, BlkSize, 1),       b1("b1", N, BlkSize, 1),
      c0("c0", N, BlkSize, 1),       c1("c1", N, BlkSize, 1);

    Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
    Kokkos::fill_random(a0, random, value_type(1.0));
    Kokkos::fill_random(b0, random, value_type(1.0));
    Kokkos::fill_random(c0, random, value_type(1.0));

    Kokkos::fence();

    Kokkos::deep_copy(a1, a0);
    Kokkos::deep_copy(b1, b0);
    Kokkos::deep_copy(c1, c0);

    /// test body
    Functor_TestBatchedSerialGemv<DeviceType,ViewType,ScalarType,
      ParamTagType,Algo::Gemv::Unblocked>(alpha, a0, b0, beta, c0).run();
    Functor_TestBatchedSerialGemv<DeviceType,ViewType,ScalarType,
      ParamTagType,AlgoTagType>(alpha, a1, b1, beta, c1).run();

    Kokkos::fence();

    /// for comparison send it to host
    typename ViewType::HostMirror c0_host = Kokkos::create_mirror_view(c0);
    typename ViewType::HostMirror c1_host = Kokkos::create_mirror_view(c1);

    Kokkos::deep_copy(c0_host, c0);
    Kokkos::deep_copy(c1_host, c1);

    /// check c0 = c1 ; this eps is about 10^-14
    typedef typename ats::mag_type mag_type;
    mag_type sum(1), diff(0);
    const mag_type eps = 1.0e3 * ats::epsilon();

    for (int k=0;k<N;++k) 
      for (int i=0;i<BlkSize;++i) 
        for (int j=0;j<1;++j) {
          sum  += ats::abs(c0_host(k,i,j));
          diff += ats::abs(c0_host(k,i,j)-c1_host(k,i,j));
        }
    EXPECT_NEAR_KK( diff/sum, 0, eps);
  }
}

template<typename DeviceType, 
         typename ValueType, 
         typename ScalarType,
         typename ParamTagType,
         typename AlgoTagType>
int test_batched_gemv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) 
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceType> ViewType;
    Test::impl_test_batched_gemv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(     0, 10);
    for (int i=0;i<10;++i) {                                                                                        
      //printf("Testing: LayoutLeft,  Blksize %d\n", i); 
      Test::impl_test_batched_gemv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(1024,  i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) 
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutRight,DeviceType> ViewType;
    Test::impl_test_batched_gemv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(     0, 10);
    for (int i=0;i<10;++i) {                                                                                        
      //printf("Testing: LayoutRight, Blksize %d\n", i); 
      Test::impl_test_batched_gemv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(1024,  i);
    }
  }
#endif
  
  return 0;
}

