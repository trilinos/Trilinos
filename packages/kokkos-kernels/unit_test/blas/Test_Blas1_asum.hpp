#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_nrm1.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class Device>
  void impl_test_asum(int N) {

    typedef typename ViewTypeA::value_type ScalarA;
    typedef Kokkos::Details::ArithTraits<ScalarA> AT;

    typedef Kokkos::View<ScalarA*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeA::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeA;


    BaseTypeA b_a("A",N);

    ViewTypeA a = Kokkos::subview(b_a,Kokkos::ALL(),0);

    typename BaseTypeA::HostMirror h_b_a = Kokkos::create_mirror_view(b_a);

    typename ViewTypeA::HostMirror h_a = Kokkos::subview(h_b_a,Kokkos::ALL(),0);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_a,rand_pool,ScalarA(10));
    Kokkos::fence();
    Kokkos::deep_copy(h_b_a,b_a);

    typename ViewTypeA::const_type c_a = a;
    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    typename AT::mag_type expected_result = 0;
    for(int i=0;i<N;i++)
      expected_result += AT::abs(h_a(i));

    typename AT::mag_type nonconst_result = KokkosBlas::asum(a);
    EXPECT_NEAR_KK( nonconst_result, expected_result, eps*expected_result);

    typename AT::mag_type const_result = KokkosBlas::asum(c_a);
    EXPECT_NEAR_KK( const_result, expected_result, eps*expected_result);

  }


}

template<class ScalarA, class Device>
int test_asum() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_asum<view_type_a_ll, Device>(0);
  Test::impl_test_asum<view_type_a_ll, Device>(13);
  Test::impl_test_asum<view_type_a_ll, Device>(1024);
  //Test::impl_test_asum<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_asum<view_type_a_lr, Device>(0);
  Test::impl_test_asum<view_type_a_lr, Device>(13);
  Test::impl_test_asum<view_type_a_lr, Device>(1024);
  //Test::impl_test_asum<view_type_a_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_asum<view_type_a_ls, Device>(0);
  Test::impl_test_asum<view_type_a_ls, Device>(13);
  Test::impl_test_asum<view_type_a_ls, Device>(1024);
  //Test::impl_test_asum<view_type_a_ls, Device>(132231);
#endif

  return 1;
}


#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, asum_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::asum_float");
    test_asum<float,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, asum_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::asum_double");
    test_asum<double,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, asum_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::asum_complex_double");
    test_asum<Kokkos::complex<double>,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, asum_int ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::asum_int");
    test_asum<int,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif


