#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas2_gemv.hpp>
#include<KokkosBlas1_dot.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device>
  void impl_test_gemv(const char* mode, int M, int N) {

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeX::value_type ScalarX;
    typedef typename ViewTypeY::value_type ScalarY;

    typedef multivector_layout_adapter<ViewTypeA> vfA_type;
    typedef Kokkos::View<ScalarX*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeX::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeX;
    typedef Kokkos::View<ScalarY*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeY::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeY;


    ScalarA alpha = 3;
    ScalarX beta = 5;
    double eps = (std::is_same<ScalarY,float>::value ? 2*1e-5 : 1e-7);

    int ldx;
    int ldy;
    if(mode[0]=='N') {
      ldx = N;
      ldy = M;
    } else {
      ldx = M;
      ldy = N;
    }
    typename vfA_type::BaseType b_A("A", M, N);
    BaseTypeX b_x("X", ldx);
    BaseTypeY b_y("Y", ldy);
    BaseTypeY b_org_y("Org_Y", ldy);
    

    ViewTypeA A = vfA_type::view(b_A);
    ViewTypeX x = Kokkos::subview(b_x,Kokkos::ALL(),0);
    ViewTypeY y = Kokkos::subview(b_y,Kokkos::ALL(),0);
    typename ViewTypeX::const_type c_x = x;
    typename ViewTypeA::const_type c_A = A;

    typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;

    typename h_vfA_type::BaseType h_b_A = Kokkos::create_mirror_view(b_A);
    typename BaseTypeX::HostMirror h_b_x = Kokkos::create_mirror_view(b_x);
    typename BaseTypeY::HostMirror h_b_y = Kokkos::create_mirror_view(b_y);

    typename ViewTypeA::HostMirror h_A = h_vfA_type::view(h_b_A);
    typename ViewTypeX::HostMirror h_x = Kokkos::subview(h_b_x,Kokkos::ALL(),0);
    typename ViewTypeY::HostMirror h_y = Kokkos::subview(h_b_y,Kokkos::ALL(),0);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_x,rand_pool,ScalarX(10));
    Kokkos::fill_random(b_y,rand_pool,ScalarY(10));
    Kokkos::fill_random(b_A,rand_pool,ScalarA(10));

    Kokkos::fence();

    Kokkos::deep_copy(b_org_y,b_y);

    Kokkos::deep_copy(h_b_x,b_x);
    Kokkos::deep_copy(h_b_y,b_y);
    Kokkos::deep_copy(h_b_A,b_A);

    typedef Kokkos::Details::ArithTraits<typename ViewTypeA::non_const_value_type> KAT;
    ScalarY expected_result = KAT:: zero();
    if(mode[0] == 'N') {
      for(int i = 0; i < M; i++) {
        ScalarY y_i = KAT::zero ();
        for(int j = 0; j < N; j++) {
           y_i += h_A(i,j) * h_x(j);
        }
        expected_result += (beta * h_y(i) + alpha * y_i) * (beta * h_y(i) + alpha * y_i) ;
      }
    } else if(mode[0] == 'T') {
      for(int j = 0; j < N; j++) {
        ScalarY y_j = KAT::zero ();
        for(int i = 0; i < M; i++) {
           y_j += h_A(i,j) * h_x(i);
        }
        expected_result += (beta * h_y(j) + alpha * y_j) * (beta * h_y(j) + alpha * y_j) ;
      }
    } else if(mode[0] == 'C') {
      for(int j = 0; j < N; j++) {
        ScalarY y_j = KAT::zero ();
        for(int i = 0; i < M; i++) {
           y_j += KAT::conj (h_A(i,j)) * h_x(i);
        }
        expected_result += (beta * h_y(j) + alpha * y_j) * (beta * h_y(j) + alpha * y_j) ;
      }
    }

    KokkosBlas::gemv(mode, alpha, A, x, beta, y);
    ScalarY nonconst_nonconst_result = KokkosBlas::dot(y, y);
    EXPECT_NEAR_KK( nonconst_nonconst_result, expected_result, eps*expected_result);
 
    Kokkos::deep_copy(b_y, b_org_y);
    KokkosBlas::gemv(mode, alpha,A ,c_x, beta, y);
    ScalarY const_nonconst_result = KokkosBlas::dot(y, y);
    EXPECT_NEAR_KK( const_nonconst_result, expected_result, eps*expected_result);

    Kokkos::deep_copy(b_y, b_org_y);
    KokkosBlas::gemv(mode, alpha, c_A, c_x, beta, y);
    ScalarY const_const_result = KokkosBlas::dot(y, y);
    EXPECT_NEAR_KK( const_const_result, expected_result, eps*expected_result);
  }
}



template<class ScalarA, class ScalarX, class ScalarY, class Device>
int test_gemv(const char* mode) {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutLeft, Device> view_type_c_ll;
  #if 0
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,100,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,150);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,150,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,200);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,200,10);
  #endif
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,0,1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,13,1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,1024,1024);
  //Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,132231,1024);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode,0,1024);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode,13,1024);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode,1024,1024);
  //Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode,132231,1024);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,0,1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,13,1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,1024,1024);
  //Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,132231,1024);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(mode,1024,1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(mode,1024,1024);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemv_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_float");
    test_gemv<float,float,float,TestExecSpace> ("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_float");
    test_gemv<float,float,float,TestExecSpace> ("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemv_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_double");
    test_gemv<double,double,double,TestExecSpace> ("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_double");
    test_gemv<double,double,double,TestExecSpace> ("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemv_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_complex_double");
    test_gemv<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_complex_double");
    test_gemv<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("T");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_conj_complex_double");
    test_gemv<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemv_int ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_int");
    test_gemv<int,int,int,TestExecSpace> ("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_int");
    test_gemv<int,int,int,TestExecSpace> ("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F( TestCategory, gemv_double_int ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_double_int");
    test_gemv<double,int,float,TestExecSpace> ("N");
  Kokkos::Profiling::popRegion();

  //Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemvt_double_int");
  //  test_gemv<double,int,float,TestExecSpace> ("T");
  //Kokkos::Profiling::popRegion();
}
#endif
