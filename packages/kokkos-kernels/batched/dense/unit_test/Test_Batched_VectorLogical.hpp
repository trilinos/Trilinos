//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// Note: Luc Berger-Vergiat 04/14/21
//       This test does not run on cuda or HIP
//       backends so the whole file is guarded
//       to ensure it is not included in these
//       backends unit-test

#if !defined(TEST_CUDA_BATCHED_DENSE_CPP) && !defined(TEST_HIP_BATCHED_DENSE_CPP) && \
    !defined(TEST_SYCL_BATCHED_DENSE_CPP) && !defined(TEST_OPENMPTARGET_BATCHED_DENSE_CPP)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

template <typename ValueType, int VectorLength>
void impl_test_batched_vector_logical() {
  /// random data initialization
  typedef Vector<SIMD<int>, VectorLength> vector_int_type;
  typedef ValueType value_type;
  const int vector_length = VectorLength;

  typedef Kokkos::ArithTraits<value_type> ats;
  typedef typename ats::mag_type mag_type;

  vector_int_type a, b;

  Random<mag_type> random;
  for (int iter = 0; iter < 100; ++iter) {
    for (int k = 0; k < vector_length; ++k) {
      a[k] = (random.value() > 0 ? 1 : -1);
      b[k] = (random.value() < 0 ? 1 : -1);
    }

    {
#undef CHECK
#define CHECK(op)                                                                   \
  {                                                                                 \
    const auto comparison = a op b;                                                 \
    for (int i = 0; i < vector_length; ++i) EXPECT_EQ(comparison[i], a[i] op b[i]); \
  }

      CHECK(||);
      CHECK(&&);

#undef CHECK
#define CHECK(op)                                                                \
  {                                                                              \
    const auto comparison = a op 0;                                              \
    for (int i = 0; i < vector_length; ++i) EXPECT_EQ(comparison[i], a[i] op 0); \
  }

      CHECK(||);
      CHECK(&&);

#undef CHECK
#define CHECK(op)                                                                \
  {                                                                              \
    const auto comparison = 0 op b;                                              \
    for (int i = 0; i < vector_length; ++i) EXPECT_EQ(comparison[i], 0 op b[i]); \
  }

      CHECK(||);
      CHECK(&&);

#undef CHECK

    }  // end test body
  }    // end for
}  // impl
}  // namespace Test

template <typename DeviceType, typename ValueType, int VectorLength>
int test_batched_vector_logical() {
  static_assert(Kokkos::SpaceAccessibility<DeviceType, Kokkos::HostSpace>::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_logical<ValueType, VectorLength>();

  return 0;
}

///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_vector_logical_simd_float3) { test_batched_vector_logical<TestDevice, float, 3>(); }
TEST_F(TestCategory, batched_vector_logical_simd_float8) { test_batched_vector_logical<TestDevice, float, 8>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_vector_logical_simd_double3) { test_batched_vector_logical<TestDevice, double, 3>(); }
TEST_F(TestCategory, batched_vector_logical_simd_double4) { test_batched_vector_logical<TestDevice, double, 4>(); }
#endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// TEST_F( TestCategory, batched_vector_logical_simd_scomplex3 ) {
//   test_batched_vector_logical<TestDevice,Kokkos::complex<float>,3>();
// }
// TEST_F( TestCategory, batched_vector_logical_simd_scomplex4 ) {
//   test_batched_vector_logical<TestDevice,Kokkos::complex<float>,4>();
// }
// #endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_logical_simd_dcomplex3 ) {
//   test_batched_vector_logical<TestDevice,Kokkos::complex<double>,3>();
// }
// TEST_F( TestCategory, batched_vector_logical_simd_dcomplex2 ) {
//   test_batched_vector_logical<TestDevice,Kokkos::complex<double>,2>();
// }
// #endif

#endif  // check to not include this in a device test
