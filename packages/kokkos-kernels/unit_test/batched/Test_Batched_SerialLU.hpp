/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

//#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched::Experimental;

namespace Test {

  template<typename DeviceType,
           typename ViewType,
           typename AlgoTagType>
  struct Functor {
    ViewType _a;

    KOKKOS_INLINE_FUNCTION
    Functor(const ViewType &a) 
      : _a(a) {} 

    KOKKOS_INLINE_FUNCTION
    void operator()(const int k) const {
      auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());

      const int iend = aa.dimension_0();
      for (int i=0;i<iend;++i)
        aa(i,i) = typename ViewType::value_type(10);

      SerialLU<AlgoTagType>::invoke(aa);
    }

    inline
    void run() {
      Kokkos::RangePolicy<DeviceType> policy(0, _a.dimension_0());
      Kokkos::parallel_for(policy, *this);
    }
  };

  template<typename DeviceType,
           typename ViewType,
           typename AlgoTagType>
  void impl_test_batched_lu(const int N, const int BlkSize) {
    typedef typename ViewType::value_type value_type;
    typedef Kokkos::Details::ArithTraits<value_type> ats;

    /// randomized input testing views
    ViewType
      a0("a0", N, BlkSize,BlkSize), a1("a1", N, BlkSize, BlkSize);

    Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
    Kokkos::fill_random(a0, random, value_type(1.0));

    Kokkos::deep_copy(a1, a0);

    Functor<DeviceType,ViewType,Algo::LU::Unblocked>(a0).run();
    Functor<DeviceType,ViewType,AlgoTagType>(a1).run();

    /// for comparison send it to host
    typename ViewType::HostMirror a0_host = Kokkos::create_mirror_view(a0);
    typename ViewType::HostMirror a1_host = Kokkos::create_mirror_view(a1);

    Kokkos::deep_copy(a0_host, a0);
    Kokkos::deep_copy(a1_host, a1);

    /// check b0 = b1 ; this eps is about 10^-14
    typedef typename ats::mag_type mag_type;
    mag_type sum(1), diff(0);
    const mag_type eps = 1.0e3 * ats::epsilon();

    for (int k=0;k<N;++k)
      for (int i=0;i<BlkSize;++i)
        for (int j=0;j<BlkSize;++j) {
          sum  += ats::abs(a0_host(k,i,j));
          diff += ats::abs(a0_host(k,i,j)-a1_host(k,i,j));
        }
    EXPECT_NEAR_KK( diff/sum, 0, eps);
  }
}


template<typename DeviceType,
         typename ValueType,
         typename AlgoTagType>
int test_batched_lu() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceType> ViewType;
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(     0, 10);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(    10, 15);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(  1024,  9);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(132231,  3);
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutRight,DeviceType> ViewType;
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(     0, 10);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(    10, 15);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(  1024,  9);
    Test::impl_test_batched_lu<DeviceType,ViewType,AlgoTagType>(132231,  3);
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_scalar_lu_float ) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestExecSpace,float,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_scalar_lu_double ) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestExecSpace,double,algo_tag_type>();
}
#endif


#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_scalar_lu_dcomplex ) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestExecSpace,Kokkos::complex<double>,algo_tag_type>();
}
#endif
