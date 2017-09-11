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
    
    constexpr int vector_length = vector_type::vector_length;
    
    typedef Kokkos::Details::ArithTraits<value_type> ats;

    vector_type a, b, c;
    value_type alpha;
    const value_type zero(0);

    Random<value_type> random;
    for (int iter=0;iter<100;++iter) {
      for (int k=0;k<vector_length;++k) {
        a[k] = random.value();
        b[k] = random.value();
        c[k] = zero;
      }
      alpha = random.value();

      typedef typename ats::mag_type mag_type;
      const mag_type eps = 1.0e1 * ats::epsilon();

      {
        /// test : vec + vec
        c = a + b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]+b[k], eps);      
      
        /// test : scalar + vec
        c = alpha + b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha+b[k], eps);      
      
        /// test : vec + scalar
        c = b + alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k] + alpha, eps);      
      }
      {
        /// test : vec - vec
        c = a - b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]-b[k], eps);      
      
        /// test : scalar - vec
        c = alpha - b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha-b[k], eps);      
      
        /// test : vec + scalar
        c = b - alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]-alpha, eps);      
      }
      {
        /// test : vec * vec
        c = a * b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]*b[k], eps);      
      
        /// test : scalar * vec
        c = alpha * b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha*b[k], eps);      
      
        /// test : vec + scalar
        c = b * alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]*alpha, eps);      
      }
      {
        /// test : vec / vec
        c = a / b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], a[k]/b[k], eps);      
      
        /// test : scalar / vec
        c = alpha / b;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], alpha/b[k], eps);      
      
        /// test : vec + scalar
        c = b / alpha;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], b[k]/alpha, eps);      
      }
      {
        /// test : vec / vec
        c = -a;
        for (int k=0;k<vector_length;++k) 
          EXPECT_NEAR_KK( c[k], -a[k], eps);      
      }
    }    
  }
}

template<typename DeviceType,typename VectorTagType>
int test_batched_vector_arithmatic() {

  static_assert(std::is_same<typename DeviceType::memory_space,Kokkos::HostSpace>::value,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_arithmatic<VectorTagType>();
  
  return 0;
}


///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_simd_float8 ) {
  typedef VectorTag<SIMD<float>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_simd_double4 ) {
  typedef VectorTag<SIMD<double>, 4> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_arithmatic_simd_dcomplex2 ) {
//   typedef VectorTag<SIMD<Kokkos::complex<double> >, 2> vector_tag_type;
//   test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
// }
#endif

///
/// AVX
///

#if defined(__AVX__) || defined(__AVX2__)
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_avx_float8 ) {
  typedef VectorTag<AVX<float>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_double4 ) {
  typedef VectorTag<AVX<double>, 4> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_arithmatic_avx_dcomplex2 ) {
//   typedef VectorTag<AVX<Kokkos::complex<double> >, 2> vector_tag_type;
//   test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
// }
#endif
#endif

///
/// AVX 512
///

#if defined(__AVX512F__)
#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_arithmatic_avx_float16 ) {
  typedef VectorTag<AVX<float>, 16> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_arithmatic_avx_double8 ) {
  typedef VectorTag<AVX<double>, 8> vector_tag_type;
  test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_arithmatic_avx_dcomplex4 ) {
//   typedef VectorTag<AVX<Kokkos::complex<double> >, 4> vector_tag_type;
//   test_batched_vector_arithmatic<TestExecSpace,vector_tag_type>();
// }
#endif
#endif
