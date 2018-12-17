#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_nrm2.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class Device>
  void impl_test_nrm2(int N) {

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

    Kokkos::fill_random(b_a,rand_pool,ScalarA(1));

    Kokkos::fence();

    Kokkos::deep_copy(h_b_a,b_a);

    typename ViewTypeA::const_type c_a = a;
    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    typename AT::mag_type expected_result = 0;
    for(int i=0;i<N;i++)
    { expected_result += AT::abs(h_a(i))*AT::abs(h_a(i)); }
    expected_result = Kokkos::Details::ArithTraits<typename AT::mag_type>::sqrt(expected_result);

    typename AT::mag_type nonconst_result = KokkosBlas::nrm2(a);
    EXPECT_NEAR_KK( nonconst_result, expected_result, eps*expected_result);

    typename AT::mag_type const_result = KokkosBlas::nrm2(c_a);
    EXPECT_NEAR_KK( const_result, expected_result, eps*expected_result);

  }

  template<class ViewTypeA, class Device>
  void impl_test_nrm2_mv(int N, int K) {

    typedef typename ViewTypeA::value_type ScalarA;
    typedef Kokkos::Details::ArithTraits<ScalarA> AT;

    typedef multivector_layout_adapter<ViewTypeA> vfA_type;

    typename vfA_type::BaseType b_a("A",N,K);

    ViewTypeA a = vfA_type::view(b_a);

    typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;

    typename h_vfA_type::BaseType h_b_a = Kokkos::create_mirror_view(b_a);

    typename ViewTypeA::HostMirror h_a = h_vfA_type::view(h_b_a);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_a,rand_pool,ScalarA(1));

    Kokkos::fence();

    Kokkos::deep_copy(h_b_a,b_a);

    typename ViewTypeA::const_type c_a = a;

    typename AT::mag_type* expected_result = new typename AT::mag_type[K];
    for(int j=0;j<K;j++) {
      expected_result[j] = typename AT::mag_type();
      for(int i=0;i<N;i++)
      { expected_result[j] += AT::abs(h_a(i,j))*AT::abs(h_a(i,j)); }
      expected_result[j] = Kokkos::Details::ArithTraits<typename AT::mag_type>::sqrt(expected_result[j]);
    }

    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    Kokkos::View<typename AT::mag_type*,Kokkos::HostSpace> r("Dot::Result",K);

    KokkosBlas::nrm2(r,a);
    for(int k=0;k<K;k++) {
      typename AT::mag_type nonconst_result = r(k);
      EXPECT_NEAR_KK( nonconst_result, expected_result[k], eps*expected_result[k]);
    }

    KokkosBlas::nrm2(r,c_a);
    for(int k=0;k<K;k++) {
      typename AT::mag_type const_result = r(k);
      EXPECT_NEAR_KK( const_result, expected_result[k], eps*expected_result[k]);
    }

    delete [] expected_result;
  }
}

template<class ScalarA, class Device>
int test_nrm2() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2<view_type_a_ll, Device>(0);
  Test::impl_test_nrm2<view_type_a_ll, Device>(13);
  Test::impl_test_nrm2<view_type_a_ll, Device>(1024);
  Test::impl_test_nrm2<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2<view_type_a_lr, Device>(0);
  Test::impl_test_nrm2<view_type_a_lr, Device>(13);
  Test::impl_test_nrm2<view_type_a_lr, Device>(1024);
  Test::impl_test_nrm2<view_type_a_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm2<view_type_a_ls, Device>(0);
  Test::impl_test_nrm2<view_type_a_ls, Device>(13);
  Test::impl_test_nrm2<view_type_a_ls, Device>(1024);
  Test::impl_test_nrm2<view_type_a_ls, Device>(132231);
#endif

  return 1;
}

template<class ScalarA, class Device>
int test_nrm2_mv() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2_mv<view_type_a_ll, Device>(0,5);
  Test::impl_test_nrm2_mv<view_type_a_ll, Device>(13,5);
  Test::impl_test_nrm2_mv<view_type_a_ll, Device>(1024,5);
  Test::impl_test_nrm2_mv<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2_mv<view_type_a_lr, Device>(0,5);
  Test::impl_test_nrm2_mv<view_type_a_lr, Device>(13,5);
  Test::impl_test_nrm2_mv<view_type_a_lr, Device>(1024,5);
  Test::impl_test_nrm2_mv<view_type_a_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm2_mv<view_type_a_ls, Device>(0,5);
  Test::impl_test_nrm2_mv<view_type_a_ls, Device>(13,5);
  Test::impl_test_nrm2_mv<view_type_a_ls, Device>(1024,5);
  Test::impl_test_nrm2_mv<view_type_a_ls, Device>(132231,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrm2_float ) {
    test_nrm2<float,TestExecSpace> ();
}
TEST_F( TestCategory, nrm2_mv_float ) {
    test_nrm2_mv<float,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrm2_double ) {
    test_nrm2<double,TestExecSpace> ();
}
TEST_F( TestCategory, nrm2_mv_double ) {
    test_nrm2_mv<double,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrm2_complex_double ) {
    test_nrm2<Kokkos::complex<double>,TestExecSpace> ();
}
TEST_F( TestCategory, nrm2_mv_complex_double ) {
    test_nrm2_mv<Kokkos::complex<double>,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, nrm2_int ) {
    test_nrm2<int,TestExecSpace> ();
}
TEST_F( TestCategory, nrm2_mv_int ) {
    test_nrm2_mv<int,TestExecSpace> ();
}
#endif


