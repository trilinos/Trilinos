#ifndef STK_SIMD_INDEX_H
#define STK_SIMD_INDEX_H

#include <stk_simd/Simd.hpp>

namespace stk { namespace simd {

struct Index {
  KOKKOS_INLINE_FUNCTION explicit Index(int i) : index(i) {}
  friend KOKKOS_INLINE_FUNCTION int int_index(const Index&);
 private:
  int index;
};

KOKKOS_INLINE_FUNCTION int int_index(const Index& i) {
  return i.index;
}

KOKKOS_INLINE_FUNCTION int int_index(const int& i) {
  return i;
}

template <typename T>
struct IndexTraits {
  typedef simd::Double double_type;
  typedef simd::Float float_type;
};

template <>
struct IndexTraits<int> {
  typedef double double_type;
  typedef float float_type;
};

#ifdef KOKKOS_HAVE_CUDA
typedef int DeviceIndex;

template <typename T>
struct DeviceTraits {
  typedef typename stk::Traits<T>::base_type simd_type;
};
#else
typedef Index DeviceIndex;

template <typename T>
struct DeviceTraits {
  typedef typename stk::Traits<T>::simd_type simd_type;
};
#endif

typedef typename IndexTraits<DeviceIndex>::double_type DeviceDouble;
typedef typename IndexTraits<DeviceIndex>::float_type  DeviceFloat;
}}



#endif
