#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_team_mult.hpp>
#include<KokkosBlas1_dot.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
  void impl_test_team_mult(int N) {
	  
    typedef Kokkos::TeamPolicy<Device>        team_policy ;
    typedef typename team_policy::member_type team_member ;

    //Launch M teams of the maximum number of threads per team
    int M = 4;
    const team_policy policy( M, Kokkos::AUTO );
    const int team_data_siz = (N%M == 0)?(N/M):(N/M+1);

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;
    typedef typename ViewTypeC::value_type ScalarC;

    typedef Kokkos::View<ScalarA*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeA::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeA;
    typedef Kokkos::View<ScalarB*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeB::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeB;
    typedef Kokkos::View<ScalarC*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeC::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeC;


    ScalarA a = 3;
    ScalarB b = 5;
    double eps = std::is_same<ScalarC,float>::value?2*1e-5:1e-7;

    BaseTypeA b_x("X",N);
    BaseTypeB b_y("Y",N);
    BaseTypeC b_z("Y",N);
    BaseTypeC b_org_z("Org_Z",N);
    

    ViewTypeA x = Kokkos::subview(b_x,Kokkos::ALL(),0);
    ViewTypeB y = Kokkos::subview(b_y,Kokkos::ALL(),0);
    ViewTypeC z = Kokkos::subview(b_z,Kokkos::ALL(),0);
    typename ViewTypeA::const_type c_x = x;
    typename ViewTypeB::const_type c_y = y;

    typename BaseTypeA::HostMirror h_b_x = Kokkos::create_mirror_view(b_x);
    typename BaseTypeB::HostMirror h_b_y = Kokkos::create_mirror_view(b_y);
    typename BaseTypeC::HostMirror h_b_z = Kokkos::create_mirror_view(b_z);

    typename ViewTypeA::HostMirror h_x = Kokkos::subview(h_b_x,Kokkos::ALL(),0);
    typename ViewTypeB::HostMirror h_y = Kokkos::subview(h_b_y,Kokkos::ALL(),0);
    typename ViewTypeC::HostMirror h_z = Kokkos::subview(h_b_z,Kokkos::ALL(),0);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_x,rand_pool,ScalarA(10));
    Kokkos::fill_random(b_y,rand_pool,ScalarB(10));
    Kokkos::fill_random(b_z,rand_pool,ScalarC(10));

    Kokkos::deep_copy(b_org_z,b_z);

    Kokkos::deep_copy(h_b_x,b_x);
    Kokkos::deep_copy(h_b_y,b_y);
    Kokkos::deep_copy(h_b_z,b_z);

    ScalarA expected_result = 0;
    for(int i=0;i<N;i++)
      expected_result += ScalarC(b*h_z(i) + a*h_x(i)*h_y(i)) * ScalarC(b*h_z(i) + a*h_x(i)*h_y(i));

    //KokkosBlas::mult(b,z,a,x,y);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), a, Kokkos::subview(x,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(y,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    ScalarC nonconst_nonconst_result = KokkosBlas::dot(z,z);
    EXPECT_NEAR_KK( nonconst_nonconst_result, expected_result, eps*expected_result);
 
    Kokkos::deep_copy(b_z,b_org_z);
    //KokkosBlas::mult(b,z,a,x,c_y);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), a, Kokkos::subview(x,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(c_y,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    ScalarC const_nonconst_result = KokkosBlas::dot(z,z);
    EXPECT_NEAR_KK( const_nonconst_result, expected_result, eps*expected_result);

    Kokkos::deep_copy(b_z,b_org_z);
    //KokkosBlas::mult(b,z,a,c_x,c_y);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), a, Kokkos::subview(c_x,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(c_y,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    ScalarC const_const_result = KokkosBlas::dot(z,z);
    EXPECT_NEAR_KK( const_const_result, expected_result, eps*expected_result);
  }

  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
  void impl_test_team_mult_mv(int N, int K) {

    typedef Kokkos::TeamPolicy<Device>        team_policy ;
    typedef typename team_policy::member_type team_member ;

    //Launch K teams of the maximum number of threads per team
    const team_policy policy( K, Kokkos::AUTO );

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;
    typedef typename ViewTypeC::value_type ScalarC;

    typedef Kokkos::View<ScalarA*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeA::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeA;
    typedef multivector_layout_adapter<ViewTypeB> vfB_type;
    typedef multivector_layout_adapter<ViewTypeC> vfC_type;

    BaseTypeA b_x("X",N);
    typename vfB_type::BaseType b_y("Y",N,K);
    typename vfC_type::BaseType b_z("Z",N,K);
    typename vfC_type::BaseType b_org_z("Z",N,K);

    ViewTypeA x = Kokkos::subview(b_x,Kokkos::ALL(),0);
    ViewTypeB y = vfB_type::view(b_y);
    ViewTypeC z = vfC_type::view(b_z);

    typedef multivector_layout_adapter<typename ViewTypeB::HostMirror> h_vfB_type;
    typedef multivector_layout_adapter<typename ViewTypeC::HostMirror> h_vfC_type;

    typename BaseTypeA::HostMirror h_b_x = Kokkos::create_mirror_view(b_x);
    typename h_vfB_type::BaseType h_b_y = Kokkos::create_mirror_view(b_y);
    typename h_vfC_type::BaseType h_b_z = Kokkos::create_mirror_view(b_z);

    typename ViewTypeA::HostMirror h_x = Kokkos::subview(h_b_x,Kokkos::ALL(),0);
    typename ViewTypeB::HostMirror h_y = h_vfB_type::view(h_b_y);
    typename ViewTypeC::HostMirror h_z = h_vfC_type::view(h_b_z);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_x,rand_pool,ScalarA(10));
    Kokkos::fill_random(b_y,rand_pool,ScalarB(10));
    Kokkos::fill_random(b_z,rand_pool,ScalarC(10));

    Kokkos::deep_copy(b_org_z,b_z);

    Kokkos::deep_copy(h_b_x,b_x);
    Kokkos::deep_copy(h_b_y,b_y);
    Kokkos::deep_copy(h_b_z,b_z);

    ScalarA a = 3;
    ScalarB b = 5;
    typename ViewTypeA::const_type c_x = x;
    typename ViewTypeB::const_type c_y = y;

    ScalarC* expected_result = new ScalarC[K];
    for(int j=0;j<K;j++) {
      expected_result[j] = ScalarC();
      for(int i=0;i<N;i++)
        expected_result[j] += ScalarC(b*h_z(i,j) + a*h_x(i)*h_y(i,j)) * ScalarC(b*h_z(i,j) + a*h_x(i)*h_y(i,j));
    }

    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    Kokkos::View<ScalarC*,Kokkos::HostSpace> r("Dot::Result",K);

    //KokkosBlas::mult(b,z,a,x,y);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z,Kokkos::ALL(),teamId), a, x, Kokkos::subview(y,Kokkos::ALL(),teamId));
    } );
    KokkosBlas::dot(r,z,z);
    for(int k=0;k<K;k++) {
      ScalarA nonconst_nonconst_result = r(k);
      EXPECT_NEAR_KK( nonconst_nonconst_result, expected_result[k], eps*expected_result[k]);
    }

    Kokkos::deep_copy(b_z,b_org_z);
    //KokkosBlas::mult(b,z,a,x,c_y);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z,Kokkos::ALL(),teamId), a, x, Kokkos::subview(c_y,Kokkos::ALL(),teamId));
    } );
    KokkosBlas::dot(r,z,z);
    for(int k=0;k<K;k++) {
      ScalarA const_non_const_result = r(k);
      EXPECT_NEAR_KK( const_non_const_result, expected_result[k], eps*expected_result[k]);
    }

    delete [] expected_result;
  }
}



template<class ScalarA, class ScalarB, class ScalarC, class Device>
int test_team_mult() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(124);
  //Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0);
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13);
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(124);
  //Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0);
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13);
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(124);
  //Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(124);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(124);
#endif

  return 1;
}

template<class ScalarA, class ScalarB, class ScalarC, class Device>
int test_team_mult_mv() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0,5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13,5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(124,5);
  //Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0,5);
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13,5);
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(124,5);
  //Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0,5);
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13,5);
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(124,5);
  //Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(124,5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(124,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_mult_float ) {
    test_team_mult<float,float,float,TestExecSpace> ();
}
TEST_F( TestCategory, team_mult_mv_float ) {
    test_team_mult_mv<float,float,float,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_mult_double ) {
    test_team_mult<double,double,double,TestExecSpace> ();
}
TEST_F( TestCategory, team_mult_mv_double ) {
    test_team_mult_mv<double,double,double,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_mult_complex_double ) {
    test_team_mult<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ();
}
TEST_F( TestCategory, team_mult_mv_complex_double ) {
    test_team_mult_mv<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_mult_int ) {
    test_team_mult<int,int,int,TestExecSpace> ();
}
TEST_F( TestCategory, team_mult_mv_int ) {
    test_team_mult_mv<int,int,int,TestExecSpace> ();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F( TestCategory, team_mult_double_int ) {
    test_team_mult<double,int,float,TestExecSpace> ();
}
TEST_F( TestCategory, team_mult_double_mv_int ) {
    test_team_mult_mv<double,int,float,TestExecSpace> ();
}
#endif
