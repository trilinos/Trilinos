#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas3_trsm.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {

  template<class ViewTypeA, class ExecutionSpace>
  struct UnitDiagTRSM {
    ViewTypeA A_;
    using ScalarA = typename ViewTypeA::value_type;

    UnitDiagTRSM (const ViewTypeA& A) : A_(A) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int& i) const {
      A_(i,i) = ScalarA(1);
    }
  };
  template<class ViewTypeA, class ExecutionSpace>
  struct NonUnitDiagTRSM {
    ViewTypeA A_;
    using ScalarA = typename ViewTypeA::value_type;

    NonUnitDiagTRSM (const ViewTypeA& A) : A_(A) {}

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
  //
  //
  //

  template<class ViewTypeA, class ViewTypeB, class Device>
  void impl_test_trsm(const char* side, const char* uplo, const char* trans, const char* diag, 
                      int M, int N, typename ViewTypeA::value_type alpha) {

    using execution_space = typename ViewTypeA::device_type::execution_space;
    using ScalarA         = typename ViewTypeA::value_type;
    using APT             = Kokkos::Details::ArithTraits<ScalarA>;
    using mag_type        = typename APT::mag_type;
    
    double machine_eps = APT::epsilon();
    bool A_l = (side[0]=='L') || (side[0]=='l');
    int K = A_l?M:N;

    //printf("KokkosBlas::trsm test for alpha %lf, %c %c %c %c, M %d, N %d, eps %.12lf, ViewType: %s\n", double(APT::abs(alpha)),side[0],uplo[0],trans[0],diag[0],M,N,1.0e8*machine_eps,typeid(ViewTypeA).name());

    ViewTypeA A  ("A", K,K);
    ViewTypeB B  ("B", M,N);
    ViewTypeB X0 ("X0",M,N);

    typename ViewTypeA::HostMirror h_A  = Kokkos::create_mirror_view(A);
    typename ViewTypeB::HostMirror h_B  = Kokkos::create_mirror_view(B);
    typename ViewTypeB::HostMirror h_X0 = Kokkos::create_mirror_view(X0);

    uint64_t seed = Kokkos::Impl::clock_tic();
    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

    if((diag[0]=='U')||(diag[0]=='u')) {
      Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max()*0.1);
      using functor_type = UnitDiagTRSM<ViewTypeA,execution_space>;
      functor_type udtrsm(A);
      Kokkos::parallel_for("KokkosBlas::Test::UnitDiagTRSM", Kokkos::RangePolicy<execution_space>(0,K), udtrsm);
    } else {//(diag[0]=='N')||(diag[0]=='n')
      Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
      using functor_type = NonUnitDiagTRSM<ViewTypeA,execution_space>;
      functor_type nudtrsm(A);
      Kokkos::parallel_for("KokkosBlas::Test::NonUnitDiagTRSM", Kokkos::RangePolicy<execution_space>(0,K), nudtrsm);
    }
    Kokkos::fill_random(X0, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());

    Kokkos::deep_copy(h_A,  A);
    Kokkos::deep_copy(h_X0, X0);

    ScalarA alpha_trmm = ScalarA(1)/alpha;
    ScalarA beta       = ScalarA(0);

    Kokkos::fence();
 
    if ((uplo[0]=='L')||(uplo[0]=='l')) {
      for (int i = 0; i < K-1; i++)
        for (int j = i+1; j < K; j++)
          h_A(i,j) = ScalarA(0);
    }
    else {
      for (int i = 1; i < K; i++)
        for (int j = 0; j < i; j++)
          h_A(i,j) = ScalarA(0); 
    }

    Kokkos::deep_copy(A, h_A);

    if (A_l){
      struct VanillaGEMM<ViewTypeB,ViewTypeA,ViewTypeB,execution_space> vgemm;
      vgemm.A_t = (trans[0]!='N') && (trans[0]!='n'); vgemm.B_t = false;
      vgemm.A_c = (trans[0]=='C') || (trans[0]=='c'); vgemm.B_c = false;
      vgemm.N = N;    vgemm.K = K;
      vgemm.A = A;    vgemm.B = X0;
      vgemm.C = B;
      vgemm.alpha = alpha_trmm;
      vgemm.beta = beta;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), vgemm);
    }
    else {
      struct VanillaGEMM<ViewTypeB,ViewTypeA,ViewTypeB,execution_space> vgemm;
      vgemm.A_t = false; vgemm.B_t = (trans[0]!='N') && (trans[0]!='n');
      vgemm.A_c = false; vgemm.B_c = (trans[0]=='C') || (trans[0]=='c');
      vgemm.N = N;     vgemm.K = K;
      vgemm.A = X0;    vgemm.B = A;
      vgemm.C = B;
      vgemm.alpha = alpha_trmm;
      vgemm.beta = beta;
      Kokkos::parallel_for("KokkosBlas::Test::VanillaGEMM", Kokkos::TeamPolicy<execution_space>(M,Kokkos::AUTO,16), vgemm);
    }
    Kokkos::fence();

    KokkosBlas::trsm(side, uplo, trans, diag, alpha, A, B);

    Kokkos::fence();

    Kokkos::deep_copy(h_B, B);

    // Checking vs ref on CPU, this eps is about 10^-6
    const mag_type eps = 1.0e8 * machine_eps;
    bool test_flag = true;
    for (int i=0; i<M; i++) {
      for (int j=0; j<N; j++) {
        if ( APT::abs(h_B(i,j) - h_X0(i,j)) > eps ) {
          test_flag = false;
          //printf("   Error: abs_result( %.15lf ) != abs_solution( %.15lf ) (abs result-solution %.15lf) at (i %ld, j %ld)\n", APT::abs(h_B(i,j)), APT::abs(h_X0(i,j)), APT::abs(h_B(i,j) - h_X0(i,j)), i, j);
          break;
        }
      }
      if (!test_flag) break;
    }
    ASSERT_EQ( test_flag, true );
  }
}

template<class ScalarA, class ScalarB, class Device>
int test_trsm(const char* mode, ScalarA alpha) {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll = Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device>;
  using view_type_b_ll = Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device>;
  Test::impl_test_trsm<view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],0,0,alpha);
  Test::impl_test_trsm<view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],101,19,alpha);
  Test::impl_test_trsm<view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],19,101,alpha);
  Test::impl_test_trsm<view_type_a_ll, view_type_b_ll, Device>(&mode[0],&mode[1],&mode[2],&mode[3],343,201,alpha);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_lr = Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device>;
  using view_type_b_lr = Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device>;
  Test::impl_test_trsm<view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],0,0,alpha);
  Test::impl_test_trsm<view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],101,19,alpha);
  Test::impl_test_trsm<view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],19,101,alpha);
  Test::impl_test_trsm<view_type_a_lr, view_type_b_lr, Device>(&mode[0],&mode[1],&mode[2],&mode[3],343,201,alpha);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trsm_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trsm_float");
    float alpha = 1.0f;
    test_trsm<float,float,TestExecSpace> ("LLNN",alpha);
    test_trsm<float,float,TestExecSpace> ("LLNU",alpha);
    test_trsm<float,float,TestExecSpace> ("LLTN",alpha);
    test_trsm<float,float,TestExecSpace> ("LLTU",alpha);
    test_trsm<float,float,TestExecSpace> ("LUNN",alpha);
    test_trsm<float,float,TestExecSpace> ("LUNU",alpha);
    test_trsm<float,float,TestExecSpace> ("LUTN",alpha);
    test_trsm<float,float,TestExecSpace> ("LUTU",alpha);

    test_trsm<float,float,TestExecSpace> ("RLNN",alpha);
    test_trsm<float,float,TestExecSpace> ("RLNU",alpha);
    test_trsm<float,float,TestExecSpace> ("RLTN",alpha);
    test_trsm<float,float,TestExecSpace> ("RLTU",alpha);
    test_trsm<float,float,TestExecSpace> ("RUNN",alpha);
    test_trsm<float,float,TestExecSpace> ("RUNU",alpha);
    test_trsm<float,float,TestExecSpace> ("RUTN",alpha);
    test_trsm<float,float,TestExecSpace> ("RUTU",alpha);

    alpha = 4.5f;
    test_trsm<float,float,TestExecSpace> ("LLNN",alpha);
    test_trsm<float,float,TestExecSpace> ("LLNU",alpha);
    test_trsm<float,float,TestExecSpace> ("LLTN",alpha);
    test_trsm<float,float,TestExecSpace> ("LLTU",alpha);
    test_trsm<float,float,TestExecSpace> ("LUNN",alpha);
    test_trsm<float,float,TestExecSpace> ("LUNU",alpha);
    test_trsm<float,float,TestExecSpace> ("LUTN",alpha);
    test_trsm<float,float,TestExecSpace> ("LUTU",alpha);

    test_trsm<float,float,TestExecSpace> ("RLNN",alpha);
    test_trsm<float,float,TestExecSpace> ("RLNU",alpha);
    test_trsm<float,float,TestExecSpace> ("RLTN",alpha);
    test_trsm<float,float,TestExecSpace> ("RLTU",alpha);
    test_trsm<float,float,TestExecSpace> ("RUNN",alpha);
    test_trsm<float,float,TestExecSpace> ("RUNU",alpha);
    test_trsm<float,float,TestExecSpace> ("RUTN",alpha);
    test_trsm<float,float,TestExecSpace> ("RUTU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trsm_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trsm_double");
    double alpha = 1.0;
    test_trsm<double,double,TestExecSpace> ("LLNN",alpha);
    test_trsm<double,double,TestExecSpace> ("LLNU",alpha);
    test_trsm<double,double,TestExecSpace> ("LLTN",alpha);
    test_trsm<double,double,TestExecSpace> ("LLTU",alpha);
    test_trsm<double,double,TestExecSpace> ("LUNN",alpha);
    test_trsm<double,double,TestExecSpace> ("LUNU",alpha);
    test_trsm<double,double,TestExecSpace> ("LUTN",alpha);
    test_trsm<double,double,TestExecSpace> ("LUTU",alpha);

    test_trsm<double,double,TestExecSpace> ("RLNN",alpha);
    test_trsm<double,double,TestExecSpace> ("RLNU",alpha);
    test_trsm<double,double,TestExecSpace> ("RLTN",alpha);
    test_trsm<double,double,TestExecSpace> ("RLTU",alpha);
    test_trsm<double,double,TestExecSpace> ("RUNN",alpha);
    test_trsm<double,double,TestExecSpace> ("RUNU",alpha);
    test_trsm<double,double,TestExecSpace> ("RUTN",alpha);
    test_trsm<double,double,TestExecSpace> ("RUTU",alpha);

    alpha = 4.5;
    test_trsm<double,double,TestExecSpace> ("LLNN",alpha);
    test_trsm<double,double,TestExecSpace> ("LLNU",alpha);
    test_trsm<double,double,TestExecSpace> ("LLTN",alpha);
    test_trsm<double,double,TestExecSpace> ("LLTU",alpha);
    test_trsm<double,double,TestExecSpace> ("LUNN",alpha);
    test_trsm<double,double,TestExecSpace> ("LUNU",alpha);
    test_trsm<double,double,TestExecSpace> ("LUTN",alpha);
    test_trsm<double,double,TestExecSpace> ("LUTU",alpha);

    test_trsm<double,double,TestExecSpace> ("RLNN",alpha);
    test_trsm<double,double,TestExecSpace> ("RLNU",alpha);
    test_trsm<double,double,TestExecSpace> ("RLTN",alpha);
    test_trsm<double,double,TestExecSpace> ("RLTU",alpha);
    test_trsm<double,double,TestExecSpace> ("RUNN",alpha);
    test_trsm<double,double,TestExecSpace> ("RUNU",alpha);
    test_trsm<double,double,TestExecSpace> ("RUTN",alpha);
    test_trsm<double,double,TestExecSpace> ("RUTU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trsm_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trsm_complex_double");
    Kokkos::complex<double> alpha = 1.0;
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCU",alpha);
    
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCU",alpha);

    alpha = Kokkos::complex<double>(4.5,0.0);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LLCU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("LUCU",alpha);
    
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RLCU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUNU",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCN",alpha);
    test_trsm<Kokkos::complex<double>,Kokkos::complex<double>,TestExecSpace> ("RUCU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, trsm_complex_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::trsm_complex_float");
    Kokkos::complex<float> alpha = 1.0f;
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCU",alpha);

    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCU",alpha);

    alpha = Kokkos::complex<float>(4.5f,0.0f);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LLCU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("LUCU",alpha);

    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RLCU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUNU",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCN",alpha);
    test_trsm<Kokkos::complex<float>,Kokkos::complex<float>,TestExecSpace> ("RUCU",alpha);
  Kokkos::Profiling::popRegion();
}
#endif
