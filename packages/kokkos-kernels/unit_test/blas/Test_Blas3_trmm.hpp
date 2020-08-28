#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas3_trmm.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {

  template<class ViewTypeA, class ExecutionSpace>
  struct UnitDiagTRMM {
    ViewTypeA A_;
    using ScalarA = typename ViewTypeA::value_type;

    UnitDiagTRMM (const ViewTypeA& A) : A_(A) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int& i) const {
      A_(i,i) = ScalarA(1);
    }
  };
  template<class ViewTypeA, class ExecutionSpace>
  struct NonUnitDiagTRMM {
    ViewTypeA A_;
    using ScalarA = typename ViewTypeA::value_type;

    NonUnitDiagTRMM (const ViewTypeA& A) : A_(A) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int& i) const {
      A_(i,i) = A_(i,i)+10;
    }
  };
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

  template<class Scalar, class ViewTypeA, class ViewTypeB, class Device>
  void impl_test_trmm(const char* side, const char* uplo, const char* trans, const char* diag, 
                      int M, int N, Scalar alpha) {

    using execution_space = typename ViewTypeA::device_type::execution_space;
    using ScalarA         = typename ViewTypeA::value_type;
    using APT             = Kokkos::Details::ArithTraits<ScalarA>;
    using mag_type        = typename APT::mag_type;

    double machine_eps = APT::epsilon();
    const mag_type eps = 1.0e8 * machine_eps; //~1e-13 for double
    bool A_l = (side[0]=='L') || (side[0]=='l');
    int K = A_l?M:N;
    ViewTypeA A  ("A", K,K);
    ViewTypeB B  ("B", M,N);
    ViewTypeB B_expected ("B_expected", M,N);
    uint64_t seed = Kokkos::Impl::clock_tic();
    ScalarA beta       = ScalarA(0);

    //printf("KokkosBlas::trmm test for alpha %g, %c %c %c %c, M %d, N %d, eps %g, ViewType: %s\n", Kokkos::Details::ArithTraits<Scalar>::real(alpha),side[0],uplo[0],trans[0],diag[0],M,N,eps,typeid(ViewTypeA).name());

    typename ViewTypeA::HostMirror host_A  = Kokkos::create_mirror_view(A);
    typename ViewTypeB::HostMirror host_B_actual  = Kokkos::create_mirror_view(B);
    typename ViewTypeB::HostMirror host_B_expected = Kokkos::create_mirror_view(B_expected);

    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

    if((diag[0]=='U')||(diag[0]=='u')) {
      // Initialize A with deterministic random numbers
      Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
      using functor_type = UnitDiagTRMM<ViewTypeA,execution_space>;
      functor_type udtrmm(A);
      // Initialize As diag with 1s
      Kokkos::parallel_for("KokkosBlas::Test::UnitDiagTRMM", Kokkos::RangePolicy<execution_space>(0,K), udtrmm);
    } else {//(diag[0]=='N')||(diag[0]=='n')
      // Initialize A with random numbers
      Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
      using functor_type = NonUnitDiagTRMM<ViewTypeA,execution_space>;
      functor_type nudtrmm(A);
      // Initialize As diag with A(i,i)+10
      Kokkos::parallel_for("KokkosBlas::Test::NonUnitDiagTRMM", Kokkos::RangePolicy<execution_space>(0,K), nudtrmm);
    }
    Kokkos::fill_random(B, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
    Kokkos::fence();
    
    Kokkos::deep_copy(host_A,  A);
    // Make host_A a lower triangle
    if ((uplo[0]=='L')||(uplo[0]=='l')) {
      for (int i = 0; i < K-1; i++)
        for (int j = i+1; j < K; j++)
          host_A(i,j) = ScalarA(0);
    }
    else {
      // Make host_A a upper triangle
      for (int i = 1; i < K; i++)
        for (int j = 0; j < i; j++)
          host_A(i,j) = ScalarA(0); 
    }
    Kokkos::deep_copy(A, host_A);

    if (A_l){
      // B_expected = alpha * op(A) * B + beta * C = 1 * op(A) * B + 0 * C
      struct VanillaGEMM<ViewTypeB,ViewTypeA,ViewTypeB,execution_space> vgemm;
      vgemm.A_t = (trans[0]!='N') && (trans[0]!='n'); vgemm.B_t = false;
      vgemm.A_c = (trans[0]=='C') || (trans[0]=='c'); vgemm.B_c = false;
      vgemm.N = N;    vgemm.K = K;
      vgemm.A = A;    vgemm.B = B;
      vgemm.C = B_expected; // out
      vgemm.alpha = alpha;
      vgemm.beta = beta;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), vgemm);
    }
    else {
      // B_expected = alpha * B * op(A) + beta * C = 1 * B * op(A) + 0 * C
      struct VanillaGEMM<ViewTypeB,ViewTypeA,ViewTypeB,execution_space> vgemm;
      vgemm.A_t = false; vgemm.B_t = (trans[0]!='N') && (trans[0]!='n');
      vgemm.A_c = false; vgemm.B_c = (trans[0]=='C') || (trans[0]=='c');
      vgemm.N = N;     vgemm.K = K;
      vgemm.A = B;    vgemm.B = A;
      vgemm.C = B_expected; // out
      vgemm.alpha = alpha;
      vgemm.beta = beta;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), vgemm);
    }
    Kokkos::fence();
    Kokkos::deep_copy(host_B_expected, B_expected);

    KokkosBlas::trmm(side, uplo, trans, diag, alpha, A, B);
    Kokkos::fence();
    Kokkos::deep_copy(host_B_actual, B);

    bool test_flag = true;
    for (int i=0; i<M; i++) {
      for (int j=0; j<N; j++) {
        if ( APT::abs(host_B_actual(i,j) - host_B_expected(i,j)) > eps ) {
          test_flag = false;
          //printf("   Error: eps ( %g ), abs_result( %.15lf ) != abs_solution( %.15lf ) (abs result-solution %g) at (i %d, j %d)\n", eps, APT::abs(host_B_actual(i,j)), APT::abs(B_expected(i,j)), APT::abs(host_B_actual(i,j) - B_expected(i,j)), i, j);
          break;
        }
      }
      if (!test_flag) break;
    }
    ASSERT_EQ( test_flag, true );
  }
}

template<class ScalarA, class ScalarB, class Device>
int test_trmm(const char* mode, ScalarA alpha) {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll = Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device>;
  using view_type_b_ll = Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device>;
  Test::impl_test_trmm<ScalarA, view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],0,0,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],101,19,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],19,101,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],12,731,alpha);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_lr = Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device>;
  using view_type_b_lr = Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device>;
  Test::impl_test_trmm<ScalarA, view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],0,0,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],101,19,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],19,101,alpha);
  Test::impl_test_trmm<ScalarA, view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],12,731,alpha);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trmm_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_float");
    float alpha = 1.0f;
    test_trmm<float,float,TestExecSpace> ("LLNN",alpha);
    test_trmm<float,float,TestExecSpace> ("LLNU",alpha);
    test_trmm<float,float,TestExecSpace> ("LLTN",alpha);
    test_trmm<float,float,TestExecSpace> ("LLTU",alpha);
    test_trmm<float,float,TestExecSpace> ("LUNN",alpha);
    test_trmm<float,float,TestExecSpace> ("LUNU",alpha);
    test_trmm<float,float,TestExecSpace> ("LUTN",alpha);
    test_trmm<float,float,TestExecSpace> ("LUTU",alpha);

    test_trmm<float,float,TestExecSpace> ("RLNN",alpha);
    test_trmm<float,float,TestExecSpace> ("RLNU",alpha);
    test_trmm<float,float,TestExecSpace> ("RLTN",alpha);
    test_trmm<float,float,TestExecSpace> ("RLTU",alpha);
    test_trmm<float,float,TestExecSpace> ("RUNN",alpha);
    test_trmm<float,float,TestExecSpace> ("RUNU",alpha);
    test_trmm<float,float,TestExecSpace> ("RUTN",alpha);
    test_trmm<float,float,TestExecSpace> ("RUTU",alpha);

    alpha = 4.5f;
    test_trmm<float,float,TestExecSpace> ("LLNN",alpha);
    test_trmm<float,float,TestExecSpace> ("LLNU",alpha);
    test_trmm<float,float,TestExecSpace> ("LLTN",alpha);
    test_trmm<float,float,TestExecSpace> ("LLTU",alpha);
    test_trmm<float,float,TestExecSpace> ("LUNN",alpha);
    test_trmm<float,float,TestExecSpace> ("LUNU",alpha);
    test_trmm<float,float,TestExecSpace> ("LUTN",alpha);
    test_trmm<float,float,TestExecSpace> ("LUTU",alpha);

    test_trmm<float,float,TestExecSpace> ("RLNN",alpha);
    test_trmm<float,float,TestExecSpace> ("RLNU",alpha);
    test_trmm<float,float,TestExecSpace> ("RLTN",alpha);
    test_trmm<float,float,TestExecSpace> ("RLTU",alpha);
    test_trmm<float,float,TestExecSpace> ("RUNN",alpha);
    test_trmm<float,float,TestExecSpace> ("RUNU",alpha);
    test_trmm<float,float,TestExecSpace> ("RUTN",alpha);
    test_trmm<float,float,TestExecSpace> ("RUTU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trmm_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_double");
    double alpha = 1.0;
    test_trmm<double,double,TestExecSpace> ("LLNN",alpha);
    test_trmm<double,double,TestExecSpace> ("LLNU",alpha);
    test_trmm<double,double,TestExecSpace> ("LLTN",alpha);
    test_trmm<double,double,TestExecSpace> ("LLTU",alpha);
    test_trmm<double,double,TestExecSpace> ("LUNN",alpha);
    test_trmm<double,double,TestExecSpace> ("LUNU",alpha);
    test_trmm<double,double,TestExecSpace> ("LUTN",alpha);
    test_trmm<double,double,TestExecSpace> ("LUTU",alpha);

    test_trmm<double,double,TestExecSpace> ("RLNN",alpha);
    test_trmm<double,double,TestExecSpace> ("RLNU",alpha);
    test_trmm<double,double,TestExecSpace> ("RLTN",alpha);
    test_trmm<double,double,TestExecSpace> ("RLTU",alpha);
    test_trmm<double,double,TestExecSpace> ("RUNN",alpha);
    test_trmm<double,double,TestExecSpace> ("RUNU",alpha);
    test_trmm<double,double,TestExecSpace> ("RUTN",alpha);
    test_trmm<double,double,TestExecSpace> ("RUTU",alpha);

    alpha = 4.5;
    test_trmm<double,double,TestExecSpace> ("LLNN",alpha);
    test_trmm<double,double,TestExecSpace> ("LLNU",alpha);
    test_trmm<double,double,TestExecSpace> ("LLTN",alpha);
    test_trmm<double,double,TestExecSpace> ("LLTU",alpha);
    test_trmm<double,double,TestExecSpace> ("LUNN",alpha);
    test_trmm<double,double,TestExecSpace> ("LUNU",alpha);
    test_trmm<double,double,TestExecSpace> ("LUTN",alpha);
    test_trmm<double,double,TestExecSpace> ("LUTU",alpha);

    test_trmm<double,double,TestExecSpace> ("RLNN",alpha);
    test_trmm<double,double,TestExecSpace> ("RLNU",alpha);
    test_trmm<double,double,TestExecSpace> ("RLTN",alpha);
    test_trmm<double,double,TestExecSpace> ("RLTU",alpha);
    test_trmm<double,double,TestExecSpace> ("RUNN",alpha);
    test_trmm<double,double,TestExecSpace> ("RUNU",alpha);
    test_trmm<double,double,TestExecSpace> ("RUTN",alpha);
    test_trmm<double,double,TestExecSpace> ("RUTU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
///////////////// alpha 1.0 /////////////////
TEST_F( TestCategory, trmm_complex_double_LLNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNU",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCN",1.0);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCU",1.0);
  Kokkos::Profiling::popRegion();
}
///////////////// alpha 4.5 /////////////////
TEST_F( TestCategory, trmm_complex_double_LLNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LLCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LLCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_LUCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_LUCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RLCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RLCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUNN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUNU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUCN");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCN",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_double_RUCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_double_RUCU");
    test_trmm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCU",Kokkos::complex<double>(4.5,0.0));
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
///////////////// alpha 1.0 /////////////////
TEST_F( TestCategory, trmm_complex_float_LLNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUNN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUNU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNU",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUCN_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCN",1.0f);
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUCU_one ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCU",1.0f);
  Kokkos::Profiling::popRegion();
}
///////////////// alpha 4.5 /////////////////
TEST_F( TestCategory, trmm_complex_float_LLNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LLCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LLCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_LUCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_LUCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RLCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RLCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUNN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUNN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUNU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUNU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUCN_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUCN");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCN",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
TEST_F( TestCategory, trmm_complex_float_RUCU_fourfive ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trmm_complex_float_RUCU");
    test_trmm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCU",Kokkos::complex<float>(4.5f,0.0f));
  Kokkos::Profiling::popRegion();
}
#endif
