/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched::Experimental;

namespace Test {

  template<typename VectorTagType>
  void impl_test_batched_vector_arithmatic() {
    /// random data initialization
    typedef Vector<VectorTagType> vector_type;
    typedef typename vector_type::value_type value_type;    
    const int vector_length = vector_type::vector_length;
    
    typedef Kokkos::Details::ArithTraits<value_type> ats;
    typedef typename ats::mag_type mag_type;

    vector_type a, b, c;
    value_type alpha;
    mag_type beta;
    const value_type zero(0);

    Random<value_type> a_random;
    Random<mag_type> b_random;
    for (int iter=0;iter<100;++iter) {
      for (int k=0;k<vector_length;++k) {
        a[k] = a_random.value();
        b[k] = a_random.value();
        c[k] = zero;
      }
      alpha = a_random.value();
      beta  = b_random.value();

      const mag_type eps = 1.0e3 * ats::epsilon();

      {
        /// test : vec + vec
        c = a + b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]+b[k], eps*c[k]);      
      
        /// test : value + vec
        c = alpha + b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha+b[k], eps*c[k]);      
      
        /// test : vec + value
        c = b + alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k] + alpha, eps*c[k]);      

        /// test : vec + mag
        c = a + beta;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k] + beta, eps*c[k]);      

        /// test : mag + vec
        c = beta + a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], beta + a[k], eps*c[k]);      
      }
      {
        /// test : vec - vec
        c = a - b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]-b[k], eps*c[k]);      
      
        /// test : value - vec
        c = alpha - b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha-b[k], eps*c[k]);      
      
        /// test : vec + value
        c = b - alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]-alpha, eps*c[k]);      

        /// test : vec - mag
        c = a - beta;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k] - beta, eps*c[k]);      

        /// test : mag - vec
        c = beta - a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], beta - a[k], eps*c[k]);      
      }
      {
        /// test : vec * vec
        c = a * b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]*b[k], eps*c[k]);      
      
        /// test : value * vec
        c = alpha * b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha*b[k], eps*c[k]);      
      
        /// test : vec + value
        c = b * alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]*alpha, eps*c[k]);      

        /// test : vec * mag
        c = a * beta;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k] * beta, eps*c[k]);      

        /// test : mag * vec
        c = beta * a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], beta * a[k], eps*c[k]);      
      }
      {
        /// test : vec / vec
        c = a / b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]/b[k], eps*c[k]);      
        
        /// test : value / vec
        c = alpha / b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha/b[k], eps*c[k]);      
        
        /// test : vec / value
        c = b / alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]/alpha, eps*c[k]);      
        
        /// test : mag / vec
        c = beta / a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], beta/a[k], eps*c[k]);      
        
        /// test : vec / value
        c = a / beta;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]/beta, eps*c[k]);      
      }
      {
        /// test : vec  -vec
        c = -a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], -a[k], eps*c[k]);      
      }
    }    
  }
}

template<typename DeviceType,typename VectorTagType>
int test_batched_vector_arithmatic() {


  static_assert(Kokkos::Impl::SpaceAccessibility<DeviceType,Kokkos::HostSpace >::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_arithmatic<VectorTagType>();
  
  return 0;
}


///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_simd_float8 ) {
  typedef VectorTag<SIMD<float,TestExecSpace>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_simd_double4 ) {
  typedef VectorTag<SIMD<double,TestExecSpace>, 4> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_simd_dcomplex2 ) {
  typedef VectorTag<SIMD<Kokkos::complex<double>,TestExecSpace>, 2> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

///
/// AVX
///

#if defined(__AVX__) || defined(__AVX2__)
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_avx_float8 ) {
  typedef VectorTag<AVX<float,TestExecSpace>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_double4 ) {
  typedef VectorTag<AVX<double,TestExecSpace>, 4> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_dcomplex2 ) {
  typedef VectorTag<AVX<Kokkos::complex<double>,TestExecSpace>, 2> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif
#endif

///
/// AVX 512
///

#if defined(__AVX512F__)
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_avx_float16 ) {
  typedef VectorTag<AVX<float,TestExecSpace>, 16> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_double8 ) {
  typedef VectorTag<AVX<double,TestExecSpace>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_dcomplex4 ) {
  typedef VectorTag<AVX<Kokkos::complex<double>,TestExecSpace>, 4> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif
#endif
