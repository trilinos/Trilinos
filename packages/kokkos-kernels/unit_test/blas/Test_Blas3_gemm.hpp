#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas3_gemm.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {

  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
  struct VanillaGEMM {
    bool A_t, B_t, A_c, B_c;
    int N,K;
    ViewTypeA A;
    ViewTypeB B;
    ViewTypeC C;

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;
    typedef typename ViewTypeC::value_type ScalarC;
    typedef Kokkos::Details::ArithTraits<ScalarC> APT;
    typedef typename APT::mag_type mag_type;
    ScalarA alpha;
    ScalarC beta;

    KOKKOS_INLINE_FUNCTION
    void operator() (const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) && !defined(__CUDA_ARCH__)
      int i = team.league_rank();
#else
      const int i = team.league_rank();
#endif
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,N), [&] (const int& j) {
        ScalarC C_ij = 0.0;

        // GNU 5.3, 5.4 and 6.1 (and maybe more) crash with another nested lambda here

#if defined(KOKKOS_COMPILER_GNU) && !defined(KOKKOS_COMPILER_NVCC)
        for(int k=0; k<K; k++) {
          ScalarA A_ik = A_t?(A_c?APT::conj(A(k,i)):A(k,i)):A(i,k);
          ScalarB B_kj = B_t?(B_c?APT::conj(B(j,k)):B(j,k)):B(k,j);
          C_ij += A_ik*B_kj;
        }
#else
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,K), [&] (const int& k, ScalarC& lsum) {
           ScalarA A_ik = A_t?(A_c?APT::conj(A(k,i)):A(k,i)):A(i,k);
           ScalarB B_kj = B_t?(B_c?APT::conj(B(j,k)):B(j,k)):B(k,j);
           lsum += A_ik*B_kj;
        },C_ij);
#endif

        C(i,j) = beta*C(i,j) + alpha*C_ij;
      });
    }
  };

  template<class ViewTypeC, class ExecutionSpace>
  struct DiffGEMM {
    int N;
    ViewTypeC C,C2;

    typedef typename ViewTypeC::value_type ScalarC;
    typedef Kokkos::Details::ArithTraits<ScalarC> APT;
    typedef typename APT::mag_type mag_type;

    KOKKOS_INLINE_FUNCTION
    void operator() (const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team, mag_type& diff) const {
      const int i = team.league_rank();
      mag_type diff_row = 0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,N), [&] (const int& j,mag_type& diff_ij) {
        //printf("A (%i %i) (%i %i) (%i %i)\n",C.extent(0),C.extent(1),C2.extent(0),C2.extent(1),i,j);
        diff_ij += APT::abs(C(i,j)-C2(i,j));
        //printf("B (%i %i) (%i %i) (%i %i)\n",C.extent(0),C.extent(1),C2.extent(0),C2.extent(1),i,j);
      },diff_row);
      Kokkos::single(Kokkos::PerTeam(team), [&] () {
        diff += diff_row;
      });
    }
  };

  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
  void impl_test_gemm(const char* TA, const char* TB, int M, int N, int K,
      typename ViewTypeA::value_type alpha,
      typename ViewTypeC::value_type beta) {


    bool A_t = (TA[0]!='N') && (TA[0]!='n');
    bool B_t = (TB[0]!='N') && (TB[0]!='n');
    bool A_c = (TA[0]=='C') || (TA[0]=='c');
    bool B_c = (TB[0]=='C') || (TB[0]=='c');
    typedef typename ViewTypeA::device_type::execution_space execution_space;
    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;
    typedef typename ViewTypeC::value_type ScalarC;
    typedef Kokkos::Details::ArithTraits<ScalarC> APT;
    typedef typename APT::mag_type mag_type;

    double machine_eps = APT::epsilon();

    ViewTypeA A("A",A_t?K:M,A_t?M:K);
    ViewTypeB B("B",B_t?N:K,B_t?K:N);
    ViewTypeC C("C",M,N);
    ViewTypeC C2("C",M,N);

    uint64_t seed = Kokkos::Impl::clock_tic();
    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

    Kokkos::fill_random(A,rand_pool,ScalarA(10));
    Kokkos::fill_random(B,rand_pool,ScalarB(10));
    Kokkos::fill_random(C,rand_pool,ScalarC(10));
    
    Kokkos::deep_copy(C2,C);

    Kokkos::fence();
 
    struct VanillaGEMM<ViewTypeA,ViewTypeB,ViewTypeC,execution_space> vgemm;
    vgemm.A_t = A_t; vgemm.B_t = B_t;
    vgemm.A_c = A_c; vgemm.B_c = B_c;
    vgemm.N = N;     vgemm.K = K;
    vgemm.A = A;     vgemm.B = B;
    vgemm.C = C2;
    vgemm.alpha = alpha;
    vgemm.beta = beta;

    Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), vgemm);

    KokkosBlas::gemm(TA,TB,alpha,A,B,beta,C);

    Kokkos::fence();

    mag_type diff_C = 0;
    struct DiffGEMM<ViewTypeC,execution_space> diffgemm;
    diffgemm.N = N;
    diffgemm.C = C;
    diffgemm.C2 = C2;

    Kokkos::parallel_reduce("KokkosBlas::Test::DiffGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), diffgemm, diff_C);

    if( N!=0 && M!=0 && K!=0 ) {
      double diff_C_average = diff_C/(N*M);
      // Expected Result: Random Walk in the least significant bit (i.e. ~ sqrt(K)*eps
      // eps scales with the total sum and has a factor in it for the accuracy of the operations ->
      // eps = K * 75 * machine_eps * 7
      double diff_C_expected = 1.0*sqrt(K)*K*75*machine_eps*7;

      //printf("Result: %e %e\n",diff_C_average,diff_C_expected);
      EXPECT_TRUE( (diff_C_average < 1.05*diff_C_expected ) );
    }
  }
}



template<class ScalarA, class ScalarB, class ScalarC, class Device>
int test_gemm(const char* mode, ScalarA alpha, ScalarB beta) {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_gemm<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(&mode[0],&mode[1],0,0,0,alpha,beta);
  Test::impl_test_gemm<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(&mode[0],&mode[1],13,15,17,alpha,beta);
  Test::impl_test_gemm<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(&mode[0],&mode[1],179,15,211,alpha,beta);
  Test::impl_test_gemm<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(&mode[0],&mode[1],12,3071,517,alpha,beta);
  Test::impl_test_gemm<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(&mode[0],&mode[1],1024,1024,2048,alpha,beta);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_gemm<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(&mode[0],&mode[1],0,0,0,alpha,beta);
  Test::impl_test_gemm<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(&mode[0],&mode[1],13,15,17,alpha,beta);
  Test::impl_test_gemm<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(&mode[0],&mode[1],179,15,211,alpha,beta);
  Test::impl_test_gemm<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(&mode[0],&mode[1],12,3071,517,alpha,beta);
  Test::impl_test_gemm<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(&mode[0],&mode[1],1024,1024,2048,alpha,beta);
#endif
/*
#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,0,1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,13,1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,1024,1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode,132231,1024);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(mode,1024,1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(mode,1024,1024);
#endif
*/
  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemm_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_float");
    float alpha = 5.0f;
    float beta = 3.0f;
    test_gemm<float,float,float,TestExecSpace> ("NN",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("TN",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("NT",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("TT",alpha,beta);

    alpha = 4.5f;
    beta = 0.0f;
    test_gemm<float,float,float,TestExecSpace> ("NN",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("TN",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("NT",alpha,beta);
    test_gemm<float,float,float,TestExecSpace> ("TT",alpha,beta);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemm_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_double");
    double alpha = 5.0;
    double beta = 3.0;
    test_gemm<double,double,double,TestExecSpace> ("NN",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("TN",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("NT",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("TT",alpha,beta);

    alpha = 4.5;
    beta = 0.0;
    test_gemm<double,double,double,TestExecSpace> ("NN",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("TN",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("NT",alpha,beta);
    test_gemm<double,double,double,TestExecSpace> ("TT",alpha,beta);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemm_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_complex_double");
    Kokkos::complex<double> alpha = 5.0;
    Kokkos::complex<double> beta = 3.0;
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("NN",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("CN",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("NC",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("CC",alpha,beta);

    alpha = Kokkos::complex<double>(4.5,0.0);
    beta = 0.0;
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("NN",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("CN",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("NC",alpha,beta);
    test_gemm<Kokkos::complex<double>,Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("CC",alpha,beta);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemm_complex_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_complex_float");
    Kokkos::complex<float> alpha = 5.0f;
    Kokkos::complex<float> beta = 3.0f;
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("NN",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("CN",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("NC",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("CC",alpha,beta);

    alpha = Kokkos::complex<float>(4.5f,0.0f);
    beta = 0.0;
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("NN",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("CN",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("NC",alpha,beta);
    test_gemm<Kokkos::complex<float>,Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("CC",alpha,beta);
  Kokkos::Profiling::popRegion();
}
#endif

/*
#if defined(KOKKOSKERNELS_INST_INT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gemm_int ) {
    test_gemm<int,int,int,TestExecSpace> ("N");
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F( TestCategory, gemm_double_int ) {
    test_gemm<double,int,float,TestExecSpace> ("N");
}
#endif
*/
