/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

  template<typename VectorViewType>
  void impl_init_vector_view(const VectorViewType & a) {
    int cnt = 0;
    for (int i0=0,i0end=a.extent(0);i0<i0end;++i0) 
      for (int i1=0,i1end=a.extent(1);i1<i1end;++i1) 
        for (int i2=0,i2end=a.extent(2);i2<i2end;++i2) 
          for (int i3=0,i3end=a.extent(3);i3<i3end;++i3) 
            for (int i4=0,i4end=a.extent(4);i4<i4end;++i4) 
              for (int i5=0,i5end=a.extent(5);i5<i5end;++i5) 
                for (int i6=0,i6end=a.extent(6);i6<i6end;++i6) 
                  for (int i7=0,i7end=a.extent(7);i7<i7end;++i7) 
                    a.access(i0,i1,i2,i3,i4,i5,i6,i7) = cnt++;
  }
#define TEST_LOOP                                                       \
  for (int i0=0,i0end=b.extent(0);i0<i0end;++i0)                        \
    for (int i1=0,i1end=b.extent(1);i1<i1end;++i1)                      \
      for (int i2=0,i2end=b.extent(2);i2<i2end;++i2)                    \
        for (int i3=0,i3end=b.extent(3);i3<i3end;++i3)                  \
          for (int i4=0,i4end=b.extent(4);i4<i4end;++i4)                \
            for (int i5=0,i5end=b.extent(5);i5<i5end;++i5)              \
              for (int i6=0,i6end=b.extent(6);i6<i6end;++i6)            \
                for (int i7=0,i7end=b.extent(7);i7<i7end;++i7)          
  
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<0> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0/vl,i1,i2,i3,i4,i5,i6,i7)[i0%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<1> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1/vl,i2,i3,i4,i5,i6,i7)[i1%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<2> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2/vl,i3,i4,i5,i6,i7)[i2%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<3> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2,i3/vl,i4,i5,i6,i7)[i3%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<4> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2,i3,i4/vl,i5,i6,i7)[i4%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<5> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2,i3,i4,i5/vl,i6,i7)[i5%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<6> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2,i3,i4,i5,i6/vl,i7)[i6%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  template<typename VectorViewType>
  void impl_verify_vector_view(const VectorViewType & a, const SimdViewAccess<VectorViewType, PackDim<7> > & b) {
    typedef typename VectorViewType::value_type vector_type;
    constexpr int vl = vector_type::vector_length;
    typedef Kokkos::Details::ArithTraits<typename vector_type::value_type> ats;
    const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
    TEST_LOOP
      EXPECT_NEAR_KK( a.access(i0,i1,i2,i3,i4,i5,i6,i7/vl)[i7%vl], b(i0,i1,i2,i3,i4,i5,i6,i7), eps );
  }
  
  template<typename DeviceType,typename VectorTagType,int VectorLength>
  void impl_test_batched_vector_view() {
    /// random data initialization
    typedef Vector<VectorTagType,VectorLength> vector_type;
    
    //typedef typename vector_type::value_type value_type;    
    //const int vector_length = vector_type::vector_length;
    const int test_view_size = 4;
    { /// rank 1 array
      Kokkos::View<vector_type*,DeviceType> a("a", test_view_size);
      impl_init_vector_view(a);
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*,DeviceType>, PackDim<0> >(a));
    }
    { /// rank 2 array
      Kokkos::View<vector_type**,DeviceType> a("a", test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type**,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type**,DeviceType>, PackDim<1> >(a));
    }
    { /// rank 3 array
      Kokkos::View<vector_type***,DeviceType> a("a", test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type***,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type***,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type***,DeviceType>, PackDim<2> >(a));
    }
    { /// rank 4 array
      Kokkos::View<vector_type****,DeviceType> a("a", test_view_size, test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type****,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type****,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type****,DeviceType>, PackDim<2> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type****,DeviceType>, PackDim<3> >(a));
    }
    { /// rank 5 array
      Kokkos::View<vector_type*****,DeviceType> a("a", test_view_size, test_view_size, test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*****,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*****,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*****,DeviceType>, PackDim<2> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*****,DeviceType>, PackDim<3> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*****,DeviceType>, PackDim<4> >(a));
    }
    { /// rank 6 array
      Kokkos::View<vector_type******,DeviceType> a("a", test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<2> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<3> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<4> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type******,DeviceType>, PackDim<5> >(a));
    }
    { /// rank 7 array
      Kokkos::View<vector_type*******,DeviceType> a("a", test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<2> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<3> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<4> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<5> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type*******,DeviceType>, PackDim<6> >(a));
    }
    { /// rank 8 array
      Kokkos::View<vector_type********,DeviceType> a("a", test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size, test_view_size);
      impl_init_vector_view(a);
      
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<0> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<1> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<2> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<3> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<4> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<5> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<6> >(a));
      impl_verify_vector_view(a, SimdViewAccess<Kokkos::View<vector_type********,DeviceType>, PackDim<7> >(a));
    }
  }
}

template<typename DeviceType,typename VectorTagType,int VectorLength>
int test_batched_vector_view() {
  static_assert(Kokkos::Impl::SpaceAccessibility<DeviceType,Kokkos::HostSpace >::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_view<DeviceType,VectorTagType,VectorLength>();
  
  return 0;
}


///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F( TestCategory, batched_vector_view_simd_float8 ) {
  test_batched_vector_view<TestExecSpace,SIMD<float>,8>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F( TestCategory, batched_vector_view_simd_double4 ) {
  test_batched_vector_view<TestExecSpace,SIMD<double>,4>();
}
TEST_F( TestCategory, batched_vector_view_simd_double8 ) {
  test_batched_vector_view<TestExecSpace,SIMD<double>,8>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F( TestCategory, batched_vector_view_simd_scomplex4 ) {
  test_batched_vector_view<TestExecSpace,SIMD<Kokkos::complex<float> >,4>();
}
TEST_F( TestCategory, batched_vector_view_simd_scomplex8 ) {
  test_batched_vector_view<TestExecSpace,SIMD<Kokkos::complex<float> >,8>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F( TestCategory, batched_vector_view_simd_dcomplex2 ) {
  test_batched_vector_view<TestExecSpace,SIMD<Kokkos::complex<double> >,2>();
}
TEST_F( TestCategory, batched_vector_view_simd_dcomplex4 ) {
  test_batched_vector_view<TestExecSpace,SIMD<Kokkos::complex<double> >,4>();
}
#endif
