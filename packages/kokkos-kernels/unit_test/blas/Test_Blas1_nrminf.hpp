#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_nrminf.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class Device>
  void impl_test_nrminf(int N) {

    typedef typename ViewTypeA::non_const_value_type ScalarA;
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

    typename AT::mag_type expected_result = Kokkos::Details::ArithTraits<typename AT::mag_type>::min();
    for(int i=0;i<N;i++)
      if(AT::abs(h_a(i)) > expected_result) expected_result = AT::abs(h_a(i));

    if(N == 0) expected_result = typename AT::mag_type(0);

    typename AT::mag_type nonconst_result = KokkosBlas::nrminf(a);
    EXPECT_NEAR_KK( nonconst_result, expected_result, eps*expected_result);

    typename AT::mag_type const_result = KokkosBlas::nrminf(c_a);
    EXPECT_NEAR_KK( const_result, expected_result, eps*expected_result);

  }

  template<class ViewTypeA, class Device>
  void impl_test_nrminf_mv(int N, int K) {

    typedef typename ViewTypeA::non_const_value_type ScalarA;
    typedef Kokkos::Details::ArithTraits<ScalarA> AT;

    typedef multivector_layout_adapter<ViewTypeA> vfA_type;

    typename vfA_type::BaseType b_a("A",N,K);

    ViewTypeA a = vfA_type::view(b_a);

    typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;

    typename h_vfA_type::BaseType h_b_a = Kokkos::create_mirror_view(b_a);

    typename ViewTypeA::HostMirror h_a = h_vfA_type::view(h_b_a);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_a,rand_pool,ScalarA(10));

    Kokkos::fence();

    Kokkos::deep_copy(h_b_a,b_a);

    typename ViewTypeA::const_type c_a = a;

    typename AT::mag_type* expected_result = new typename AT::mag_type[K];
    for(int j=0;j<K;j++) {
      expected_result[j] = Kokkos::Details::ArithTraits<typename AT::mag_type>::min();
      for(int i=0;i<N;i++) {
        if(AT::abs(h_a(i,j)) > expected_result[j]) expected_result[j] = AT::abs(h_a(i,j));
      }
      if(N == 0) expected_result[j] = typename AT::mag_type(0);
    }

    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    Kokkos::View<typename AT::mag_type*,Kokkos::HostSpace> r("Dot::Result",K);

    KokkosBlas::nrminf(r,a);
    for(int k=0;k<K;k++) {
      typename AT::mag_type nonconst_result = r(k);
      typename AT::mag_type exp_result = expected_result[k];
      EXPECT_NEAR_KK( nonconst_result, exp_result, eps*exp_result);
    }

   /* KokkosBlas::nrminf(r,c_a);
    for(int k=0;k<K;k++) {
      typename AT::mag_type const_result = r(k);
      typename AT::mag_type exp_result = expected_result[k];
      EXPECT_NEAR_KK( const_result, exp_result, eps*exp_result);
    }
*/
    delete [] expected_result;
  }
}

template<class ScalarA, class Device>
int test_nrminf() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrminf<view_type_a_ll, Device>(0);
  Test::impl_test_nrminf<view_type_a_ll, Device>(13);
  Test::impl_test_nrminf<view_type_a_ll, Device>(1024);
  //Test::impl_test_nrminf<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrminf<view_type_a_lr, Device>(0);
  Test::impl_test_nrminf<view_type_a_lr, Device>(13);
  Test::impl_test_nrminf<view_type_a_lr, Device>(1024);
  //Test::impl_test_nrminf<view_type_a_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrminf<view_type_a_ls, Device>(0);
  Test::impl_test_nrminf<view_type_a_ls, Device>(13);
  Test::impl_test_nrminf<view_type_a_ls, Device>(1024);
  //Test::impl_test_nrminf<view_type_a_ls, Device>(132231);
#endif

  return 1;
}

template<class ScalarA, class Device>
int test_nrminf_mv() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrminf_mv<view_type_a_ll, Device>(0,5);
  Test::impl_test_nrminf_mv<view_type_a_ll, Device>(13,5);
  Test::impl_test_nrminf_mv<view_type_a_ll, Device>(1024,5);
  //Test::impl_test_nrminf_mv<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrminf_mv<view_type_a_lr, Device>(0,5);
  Test::impl_test_nrminf_mv<view_type_a_lr, Device>(13,5);
  Test::impl_test_nrminf_mv<view_type_a_lr, Device>(1024,5);
  //Test::impl_test_nrminf_mv<view_type_a_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrminf_mv<view_type_a_ls, Device>(0,5);
  Test::impl_test_nrminf_mv<view_type_a_ls, Device>(13,5);
  Test::impl_test_nrminf_mv<view_type_a_ls, Device>(1024,5);
  //Test::impl_test_nrminf_mv<view_type_a_ls, Device>(132231,5);
#endif

  return 1;}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrminf_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_float");
    test_nrminf<float,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, nrminf_mv_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_mvfloat");
    test_nrminf_mv<float,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrminf_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_double");
    test_nrminf<double,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, nrminf_mv_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_mv_double");
    test_nrminf_mv<double,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrminf_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_complex_double");
    test_nrminf<Kokkos::complex<double>,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, nrminf_mv_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_mv_complex_double");
    test_nrminf_mv<Kokkos::complex<double>,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrminf_int ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_int");
    test_nrminf<int,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, nrminf_mv_int ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrminf_mv_int");
    test_nrminf_mv<int,TestExecSpace> ();
  Kokkos::Profiling::popRegion();
}
#endif


