#include<gtest/gtest.h>
#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>

#include<KokkosBlas_gesv.hpp>
#include<KokkosBlas2_gemv.hpp>
#include<KokkosBlas3_gemm.hpp>
#include<KokkosKernels_TestUtils.hpp>

namespace Test {

template<class ViewTypeA, class ViewTypeB, class Device>
void impl_test_gesv(const char* mode, const char* padding, int N) {
    typedef typename Device::execution_space execution_space;
    typedef typename ViewTypeA::value_type ScalarA;
    typedef Kokkos::Details::ArithTraits<ScalarA> ats;

    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

    int ldda, lddb;

    if(padding[0]=='Y') {//rounded up to multiple of 32
      ldda = ((N+32-1)/32)*32;
      lddb = ldda;
    }
    else {
      ldda = N;
      lddb = N;
    }

    // Create device views
    ViewTypeA A ( "A",  ldda, N );
    ViewTypeB X0( "X0", N );
    ViewTypeB B ( "B",  lddb );

    // Create host mirrors of device views.
    typename ViewTypeB::HostMirror h_X0 = Kokkos::create_mirror_view(X0);
    typename ViewTypeB::HostMirror h_B  = Kokkos::create_mirror(B);
    
    // Initialize data.
    Kokkos::fill_random(A, rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());
    Kokkos::fill_random(X0,rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());

    // Generate RHS B = A*X0.
    ScalarA alpha = 1.0;
    ScalarA beta  = 0.0;

    KokkosBlas::gemv("N",alpha,A,X0,beta,B); Kokkos::fence();

    // Deep copy device view to host view.
    Kokkos::deep_copy( h_X0, X0 );

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
    if( std::is_same< typename Device::execution_space, Kokkos::Cuda >::value ) {
      // Allocate IPIV view on host
      typedef Kokkos::View<magma_int_t*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ViewTypeP;
      magma_int_t *ipiv_raw = nullptr;
      int Nt = 0;
      if(mode[0]=='Y') {
        Nt = N;
        magma_imalloc_cpu( &ipiv_raw, Nt );
      }
      ViewTypeP ipiv(ipiv_raw, Nt);
	  
      // Solve.
      KokkosBlas::gesv(A,B,ipiv);
      Kokkos::fence();
      
      // Get the solution vector.
      Kokkos::deep_copy( h_B, B );
      
      // Checking vs ref on CPU, this eps is about 10^-9
      typedef typename ats::mag_type mag_type;
      const mag_type eps = 1.0e7 * ats::epsilon();
      bool test_flag = true;
      for (int i=0; i<N; i++) {
        if ( ats::abs(h_B(i) - h_X0(i)) > eps ) {
          test_flag = false;
          //printf( "    Error %d, pivot %c, padding %c: result( %.15lf ) != solution( %.15lf ) at (%ld)\n", N, mode[0], padding[0], ats::abs(h_B(i)), ats::abs(h_X0(i)), i );		
          break;
        }
      }	
      ASSERT_EQ( test_flag, true );

      if(mode[0]=='Y') {
        magma_free_cpu( ipiv_raw );
      }
    }
#else
    // Allocate IPIV view on host
    typedef Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> ViewTypeP;
    int Nt = 0;
    if(mode[0]=='Y') Nt = N;
    ViewTypeP ipiv("IPIV", Nt);
	
    // Solve.
    KokkosBlas::gesv(A,B,ipiv);
    Kokkos::fence();

    // Get the solution vector.
    Kokkos::deep_copy( h_B, B );

    // Checking vs ref on CPU, this eps is about 10^-9
    typedef typename ats::mag_type mag_type;
    const mag_type eps = 1.0e7 * ats::epsilon();
    bool test_flag = true;
    for (int i=0; i<N; i++) {
      if ( ats::abs(h_B(i) - h_X0(i)) > eps ) {
        test_flag = false;
        //printf( "    Error %d, pivot %c, padding %c: result( %.15lf ) != solution( %.15lf ) at (%ld)\n", N, mode[0], padding[0], ats::abs(h_B(i)), ats::abs(h_X0(i)), i );		
        break;
      }
    }	
    ASSERT_EQ( test_flag, true );
#endif

  }

template<class ViewTypeA, class ViewTypeB, class Device>
void impl_test_gesv_mrhs(const char* mode, const char* padding, int N, int nrhs) {
    typedef typename Device::execution_space execution_space;
    typedef typename ViewTypeA::value_type ScalarA;
    typedef Kokkos::Details::ArithTraits<ScalarA> ats;

    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

    int ldda, lddb;
	
    if(padding[0]=='Y') {//rounded up to multiple of 32
      ldda = ((N+32-1)/32)*32;
      lddb = ldda;
    }
    else {
      ldda = N;
      lddb = N;
    }

    // Create device views
    ViewTypeA A ( "A",  ldda, N );
    ViewTypeB X0( "X0", N, nrhs );
    ViewTypeB B ( "B",  lddb, nrhs );

    // Create host mirrors of device views.
    typename ViewTypeB::HostMirror h_X0 = Kokkos::create_mirror_view( X0 );
    typename ViewTypeB::HostMirror h_B  = Kokkos::create_mirror( B );
    
    // Initialize data.
    Kokkos::fill_random(A, rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());
    Kokkos::fill_random(X0,rand_pool,Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,ScalarA >::max());

    // Generate RHS B = A*X0.
    ScalarA alpha = 1.0;
    ScalarA beta  = 0.0;

    KokkosBlas::gemm("N","N",alpha,A,X0,beta,B); Kokkos::fence();

    // Deep copy device view to host view.
    Kokkos::deep_copy( h_X0, X0 );

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
    if( std::is_same< typename Device::execution_space, Kokkos::Cuda >::value ) {
      // Allocate IPIV view on host
      typedef Kokkos::View<magma_int_t*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ViewTypeP;
      magma_int_t *ipiv_raw = nullptr;
      int Nt = 0;
      if(mode[0]=='Y') {
        Nt = N;
        magma_imalloc_cpu( &ipiv_raw, Nt );
      }
      ViewTypeP ipiv(ipiv_raw, Nt);
      
      // Solve.	
      KokkosBlas::gesv(A,B,ipiv);
      Kokkos::fence();
      
      // Get the solution vector.
      Kokkos::deep_copy( h_B, B );
      
      // Checking vs ref on CPU, this eps is about 10^-9
      typedef typename ats::mag_type mag_type;
      const mag_type eps = 1.0e7 * ats::epsilon();
      bool test_flag = true;
      for (int j=0; j<nrhs; j++) {
        for (int i=0; i<N; i++) {
          if ( ats::abs(h_B(i,j) - h_X0(i,j)) > eps ) {
            test_flag = false;
            //printf( "    Error %d, pivot %c, padding %c: result( %.15lf ) != solution( %.15lf ) at (%ld) at rhs %d\n", N, mode[0], padding[0], ats::abs(h_B(i,j)), ats::abs(h_X0(i,j)), i, j );		
            break;
          }
        }
        if (test_flag == false) break;
      }
      ASSERT_EQ( test_flag, true );

      if(mode[0]=='Y') {
        magma_free_cpu( ipiv_raw );
      }
    }
#else
    // Allocate IPIV view on host
    typedef Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> ViewTypeP;
    int Nt = 0;
    if(mode[0]=='Y') Nt = N;
    ViewTypeP ipiv("IPIV", Nt);

    // Solve.	
    KokkosBlas::gesv(A,B,ipiv);
    Kokkos::fence();

    // Get the solution vector.
    Kokkos::deep_copy( h_B, B );

    // Checking vs ref on CPU, this eps is about 10^-9
    typedef typename ats::mag_type mag_type;
    const mag_type eps = 1.0e7 * ats::epsilon();
    bool test_flag = true;
    for (int j=0; j<nrhs; j++) {
      for (int i=0; i<N; i++) {
        if ( ats::abs(h_B(i,j) - h_X0(i,j)) > eps ) {
          test_flag = false;
          //printf( "    Error %d, pivot %c, padding %c: result( %.15lf ) != solution( %.15lf ) at (%ld) at rhs %d\n", N, mode[0], padding[0], ats::abs(h_B(i,j)), ats::abs(h_X0(i,j)), i, j );		
          break;
        }
      }
      if (test_flag == false) break;
    }
    ASSERT_EQ( test_flag, true );
#endif

  }

}//namespace Test

template<class Scalar, class Device>
int test_gesv(const char* mode) {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<Scalar*,  Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 2);   //no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 13);  //no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 179); //no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 64);  //no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 1024);//no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "Y", 13);  //padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "Y", 179); //padding
#endif

/*
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*,  Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 2);   //no padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 13);  //no padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 179); //no padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 64);  //no padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 1024);//no padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "Y", 13);  //padding
  Test::impl_test_gesv<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "Y", 179); //padding
#endif
*/

  return 1;
}

template<class Scalar, class Device>
int test_gesv_mrhs(const char* mode) {

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 2,   5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 13,  5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 179, 5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 64,  5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "N", 1024,5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "Y", 13,  5);//padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device>(&mode[0], "Y", 179, 5);//padding
#endif

/*
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 2,   5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 13,  5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 179, 5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 64,  5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "N", 1024,5);//no padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "Y", 13,  5);//padding
  Test::impl_test_gesv_mrhs<view_type_a_lr, view_type_b_lr, Device>(&mode[0], "Y", 179, 5);//padding
#endif
*/

  return 1;
}

#if defined( KOKKOSKERNELS_ENABLE_TPL_MAGMA ) || defined (KOKKOSKERNELS_ENABLE_TPL_BLAS)

#if defined(KOKKOSKERNELS_INST_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gesv_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_float");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv<float,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv<float,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F( TestCategory, gesv_mrhs_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_mrhs_float");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv_mrhs<float,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv_mrhs<float,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gesv_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_double");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv<double,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv<double,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F( TestCategory, gesv_mrhs_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_mrhs_double");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv_mrhs<double,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv_mrhs<double,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gesv_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_complex_double");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv<Kokkos::complex<double>,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv<Kokkos::complex<double>,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F( TestCategory, gesv_mrhs_complex_double ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_mrhs_complex_double");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv_mrhs<Kokkos::complex<double>,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv_mrhs<Kokkos::complex<double>,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F( TestCategory, gesv_complex_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_complex_float");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv<Kokkos::complex<float>,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv<Kokkos::complex<float>,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F( TestCategory, gesv_mrhs_complex_float ) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gesv_mrhs_complex_float");
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  if( std::is_same< typename TestExecSpace::execution_space, Kokkos::Cuda >::value )
    test_gesv_mrhs<Kokkos::complex<float>,TestExecSpace> ("N");//No pivoting
#endif
    test_gesv_mrhs<Kokkos::complex<float>,TestExecSpace> ("Y");//Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#endif//KOKKOSKERNELS_ENABLE_TPL_MAGMA || KOKKOSKERNELS_ENABLE_TPL_BLAS