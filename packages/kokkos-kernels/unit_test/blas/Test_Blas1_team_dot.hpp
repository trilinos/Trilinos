#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas1_team_dot.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {
  template<class ViewTypeA, class ViewTypeB, class Device>
  void impl_test_team_dot(int N) {

    typedef Kokkos::TeamPolicy<Device>        team_policy ;
    typedef typename team_policy::member_type team_member ;

    //Launch M teams of the maximum number of threads per team
    int M = 4;
    const team_policy policy( M, Kokkos::AUTO );
    const int team_data_siz = (N%M == 0)?(N/M):(N/M+1);

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;

    typedef Kokkos::View<ScalarA*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeA::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeA;
    typedef Kokkos::View<ScalarB*[2],
       typename std::conditional<
                std::is_same<typename ViewTypeB::array_layout,Kokkos::LayoutStride>::value,
                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,Device> BaseTypeB;


    BaseTypeA b_a("A",N);
    BaseTypeB b_b("B",N);

    ViewTypeA a = Kokkos::subview(b_a,Kokkos::ALL(),0);
    ViewTypeB b = Kokkos::subview(b_b,Kokkos::ALL(),0);

    typename BaseTypeA::HostMirror h_b_a = Kokkos::create_mirror_view(b_a);
    typename BaseTypeB::HostMirror h_b_b = Kokkos::create_mirror_view(b_b);

    typename ViewTypeA::HostMirror h_a = Kokkos::subview(h_b_a,Kokkos::ALL(),0);
    typename ViewTypeB::HostMirror h_b = Kokkos::subview(h_b_b,Kokkos::ALL(),0);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_a,rand_pool,ScalarA(10));
    Kokkos::fill_random(b_b,rand_pool,ScalarB(10));

    Kokkos::fence();

    Kokkos::deep_copy(h_b_a,b_a);
    Kokkos::deep_copy(h_b_b,b_b);

    ScalarA expected_result = 0;
    for(int i=0;i<N;i++)
      expected_result += h_a(i)*h_b(i);

    Kokkos::View<ScalarB*,Kokkos::HostSpace> r("PartialDots",M);
    Kokkos::View<ScalarB*,Device> d_r("PartialDots",M);

    //ScalarA nonconst_nonconst_result = KokkosBlas::dot(a,b);
    ScalarA nonconst_nonconst_result = 0;

    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(a,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(b,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<M;k++)
      nonconst_nonconst_result += r(k);

    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;
    EXPECT_NEAR_KK( nonconst_nonconst_result, expected_result, eps*expected_result);
    
	typename ViewTypeA::const_type c_a = a;
    typename ViewTypeB::const_type c_b = b;

    //ScalarA const_const_result = KokkosBlas::dot(c_a,c_b);
    ScalarA const_const_result = 0;

    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(c_a,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(c_b,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<M;k++)
      const_const_result += r(k);

    EXPECT_NEAR_KK( const_const_result, expected_result, eps*expected_result);

    //ScalarA nonconst_const_result = KokkosBlas::dot(a,c_b);
    ScalarA nonconst_const_result = 0;
    
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(a,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(c_b,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    Kokkos::deep_copy(r,d_r);    
    for(int k=0;k<M;k++)
      nonconst_const_result += r(k);

    EXPECT_NEAR_KK( nonconst_const_result, expected_result, eps*expected_result);

    //ScalarA const_nonconst_result = KokkosBlas::dot(c_a,b);
    ScalarA const_nonconst_result = 0;
    
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(c_a,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)), Kokkos::subview(b,Kokkos::make_pair(teamId*team_data_siz,(teamId < M-1)?(teamId+1)*team_data_siz:N)));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<M;k++)
      const_nonconst_result += r(k);

    EXPECT_NEAR_KK( const_nonconst_result, expected_result, eps*expected_result);
  }

  template<class ViewTypeA, class ViewTypeB, class Device>
  void impl_test_team_dot_mv(int N, int K) {

    typedef Kokkos::TeamPolicy<Device>        team_policy ;
    typedef typename team_policy::member_type team_member ;

    //Launch K teams of the maximum number of threads per team
    const team_policy policy( K, Kokkos::AUTO );

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;

    typedef multivector_layout_adapter<ViewTypeA> vfA_type;
    typedef multivector_layout_adapter<ViewTypeB> vfB_type;

    typename vfA_type::BaseType b_a("A",N,K);
    typename vfB_type::BaseType b_b("B",N,K);

    ViewTypeA a = vfA_type::view(b_a);
    ViewTypeB b = vfB_type::view(b_b);

    typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;
    typedef multivector_layout_adapter<typename ViewTypeB::HostMirror> h_vfB_type;

    typename h_vfA_type::BaseType h_b_a = Kokkos::create_mirror_view(b_a);
    typename h_vfB_type::BaseType h_b_b = Kokkos::create_mirror_view(b_b);

    typename ViewTypeA::HostMirror h_a = h_vfA_type::view(h_b_a);
    typename ViewTypeB::HostMirror h_b = h_vfB_type::view(h_b_b);

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    Kokkos::fill_random(b_a,rand_pool,ScalarA(10));
    Kokkos::fill_random(b_b,rand_pool,ScalarB(10));

    Kokkos::fence();

    Kokkos::deep_copy(h_b_a,b_a);
    Kokkos::deep_copy(h_b_b,b_b);

    typename ViewTypeA::const_type c_a = a;
    typename ViewTypeB::const_type c_b = b;

    ScalarA* expected_result = new ScalarA[K];
    for(int j=0;j<K;j++) {
      expected_result[j] = ScalarA();
      for(int i=0;i<N;i++)
        expected_result[j] += h_a(i,j)*h_b(i,j);
    }

    double eps = std::is_same<ScalarA,float>::value?2*1e-5:1e-7;

    Kokkos::View<ScalarB*,Kokkos::HostSpace> r("Dot::Result",K);
    Kokkos::View<ScalarB*,Device> d_r("Dot::Result",K);

    //KokkosBlas::dot(r,a,b);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(a,Kokkos::ALL(),teamId), Kokkos::subview(b,Kokkos::ALL(),teamId));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<K;k++) {
      ScalarA nonconst_nonconst_result = r(k);
      EXPECT_NEAR_KK( nonconst_nonconst_result, expected_result[k], eps*expected_result[k]);
    }

    //KokkosBlas::dot(r,c_a,c_b);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(c_a,Kokkos::ALL(),teamId), Kokkos::subview(c_b,Kokkos::ALL(),teamId));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<K;k++) {
      ScalarA const_const_result = r(k);
      EXPECT_NEAR_KK( const_const_result, expected_result[k], eps*expected_result[k]);
    }

    //KokkosBlas::dot(r,a,c_b);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(a,Kokkos::ALL(),teamId), Kokkos::subview(c_b,Kokkos::ALL(),teamId));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<K;k++) {
      ScalarA non_const_const_result = r(k);
      EXPECT_NEAR_KK( non_const_const_result, expected_result[k], eps*expected_result[k]);
    }

    //KokkosBlas::dot(r,c_a,b);
    Kokkos::parallel_for( "KokkosBlas::Test::TeamDot", policy, KOKKOS_LAMBDA ( const team_member &teamMember ) {
       const int teamId = teamMember.league_rank();
       d_r(teamId) = KokkosBlas::Experimental::dot(teamMember, Kokkos::subview(c_a,Kokkos::ALL(),teamId), Kokkos::subview(b,Kokkos::ALL(),teamId));
    } );
    Kokkos::deep_copy(r,d_r);
    for(int k=0;k<K;k++) {
      ScalarA const_non_const_result = r(k);
      EXPECT_NEAR_KK( const_non_const_result, expected_result[k], eps*expected_result[k]);
    }

    delete [] expected_result;
  }
}

template<class ScalarA, class ScalarB, class Device>
int test_team_dot() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_team_dot<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_team_dot<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_team_dot<view_type_a_ll, view_type_b_ll, Device>(124);
  //Test::impl_test_team_dot<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_team_dot<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_team_dot<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_team_dot<view_type_a_lr, view_type_b_lr, Device>(124);
  //Test::impl_test_team_dot<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_team_dot<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_team_dot<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_team_dot<view_type_a_ls, view_type_b_ls, Device>(124);
  //Test::impl_test_team_dot<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_dot<view_type_a_ls, view_type_b_ll, Device>(124);
  Test::impl_test_team_dot<view_type_a_ll, view_type_b_ls, Device>(124);
#endif

  return 1;
}

template<class ScalarA, class ScalarB, class Device>
int test_team_dot_mv() {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_team_dot_mv<view_type_a_ll, view_type_b_ll, Device>(0,5);
  Test::impl_test_team_dot_mv<view_type_a_ll, view_type_b_ll, Device>(13,5);
  Test::impl_test_team_dot_mv<view_type_a_ll, view_type_b_ll, Device>(124,5);
  //Test::impl_test_team_dot_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_team_dot_mv<view_type_a_lr, view_type_b_lr, Device>(0,5);
  Test::impl_test_team_dot_mv<view_type_a_lr, view_type_b_lr, Device>(13,5);
  Test::impl_test_team_dot_mv<view_type_a_lr, view_type_b_lr, Device>(124,5);
  //Test::impl_test_team_dot_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_team_dot_mv<view_type_a_ls, view_type_b_ls, Device>(0,5);
  Test::impl_test_team_dot_mv<view_type_a_ls, view_type_b_ls, Device>(13,5);
  Test::impl_test_team_dot_mv<view_type_a_ls, view_type_b_ls, Device>(124,5);
  //Test::impl_test_team_dot_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_dot_mv<view_type_a_ls, view_type_b_ll, Device>(124,5);
  Test::impl_test_team_dot_mv<view_type_a_ll, view_type_b_ls, Device>(124,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_dot_float ) {
    test_team_dot<float,float,TestExecSpace> ();
}
TEST_F( TestCategory, team_dot_mv_float ) {
    test_team_dot_mv<float,float,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_dot_double ) {
    test_team_dot<double,double,TestExecSpace> ();
}
TEST_F( TestCategory, team_dot_mv_double ) {
    test_team_dot_mv<double,double,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_dot_complex_double ) {
    test_team_dot<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ();
}
TEST_F( TestCategory, team_dot_mv_complex_double ) {
    test_team_dot_mv<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, team_dot_int ) {
    test_team_dot<int,int,TestExecSpace> ();
}
TEST_F( TestCategory, team_dot_mv_int ) {
    test_team_dot_mv<int,int,TestExecSpace> ();
}
#endif

/*#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F( TestCategory, team_dot_double_int ) {
    test_team_dot<double,int,TestExecSpace> ();
}
TEST_F( TestCategory, team_dot_mv_double_int ) {
    test_team_dot_mv<double,int,TestExecSpace> ();
}
#endif*/

