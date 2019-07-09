/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

  template<typename VectorTagType,int VectorLength>
  void impl_test_batched_vector_relation() {
    /// random data initialization
    typedef Vector<VectorTagType,VectorLength> vector_type;
    
    typedef typename vector_type::value_type value_type;    
    const int vector_length = vector_type::vector_length;
    
    //typedef Kokkos::Details::ArithTraits<value_type> ats;
    //typedef typename ats::mag_type mag_type;

    vector_type a, b;

    Random<value_type> random;
    for (int iter=0;iter<100;++iter) {
      for (int k=0;k<vector_length;++k) {
        a[k] = random.value();
        b[k] = random.value();
      }

      {

#undef CHECK
#define CHECK(op)                                       \
        {                                               \
          const auto comparison = a op b;               \
          for (int i=0;i<vector_length;++i)             \
            EXPECT_EQ( comparison[i], a[i] op b[i]);    \
        }
        
        CHECK(<);
        CHECK(>);
        CHECK(<=);
        CHECK(>=);
        CHECK(==);
        CHECK(!=);

#undef CHECK
#define CHECK(op)                                               \
        {                                                       \
          const auto comparison = a op value_type(0);           \
          for (int i=0;i<vector_length;++i)                     \
            EXPECT_EQ( comparison[i], a[i] op value_type(0));   \
        }
        
        CHECK(<);
        CHECK(>);
        CHECK(<=);
        CHECK(>=);
        CHECK(==);
        CHECK(!=);

#undef CHECK
#define CHECK(op)                                               \
        {                                                       \
          const auto comparison = value_type(0) op b;           \
          for (int i=0;i<vector_length;++i)                     \
            EXPECT_EQ( comparison[i], value_type(0) op b[i]);   \
        }
        
        CHECK(<);
        CHECK(>);
        CHECK(<=);
        CHECK(>=);
        CHECK(==);
        CHECK(!=);

#undef CHECK

      } // end test body
    } // end for
  } // impl
} // namespace

template<typename DeviceType,typename VectorTagType,int VectorLength>
int test_batched_vector_relation() {
  static_assert(Kokkos::Impl::SpaceAccessibility<DeviceType,Kokkos::HostSpace >::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_relation<VectorTagType,VectorLength>();
  
  return 0;
}


///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_relation_simd_float3 ) {
  test_batched_vector_relation<TestExecSpace,SIMD<float>,3>();
}
TEST_F( TestCategory, batched_vector_relation_simd_float8 ) {
  test_batched_vector_relation<TestExecSpace,SIMD<float>,8>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_relation_simd_double3 ) {
  test_batched_vector_relation<TestExecSpace,SIMD<double>,3>();
}
TEST_F( TestCategory, batched_vector_relation_simd_double4 ) {
  test_batched_vector_relation<TestExecSpace,SIMD<double>,4>();
}
#endif

/// comparison of complex variables is not defined

// #if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// TEST_F( TestCategory, batched_vector_relation_simd_scomplex4 ) {
//   test_batched_vector_relation<TestExecSpace,SIMD<Kokkos::complex<float> >,4>();
// }
// #endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_relation_simd_dcomplex2 ) {
//   test_batched_vector_relation<TestExecSpace,SIMD<Kokkos::complex<double> >,2>();
// }
// #endif
