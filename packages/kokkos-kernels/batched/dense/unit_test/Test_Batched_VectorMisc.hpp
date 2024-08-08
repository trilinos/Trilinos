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

template <typename VectorTagType, int VectorLength>
void impl_test_batched_vector_misc() {
  /// random data initialization
  typedef Vector<VectorTagType, VectorLength> vector_type;

  typedef typename vector_type::value_type value_type;
  const int vector_length = vector_type::vector_length;

  typedef Kokkos::ArithTraits<value_type> ats;
  typedef typename ats::mag_type mag_type;

  vector_type a, b, c;
  const mag_type eps = 1.0e3 * ats::epsilon();

  Random<value_type> random;
  for (int iter = 0; iter < 100; ++iter) {
    for (int k = 0; k < vector_length; ++k) {
      a[k] = random.value();
      b[k] = random.value();
    }

    {
      c = conditional_assign(a < b, a, b);
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? a[i] : b[i];
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }

      c = 0;
      conditional_assign(c, a < b, a, b);
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? a[i] : b[i];
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }
    }
    {
      c = conditional_assign(a < b, a, value_type(0));
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? a[i] : 0;
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }

      c = 0;
      conditional_assign(c, a < b, a, value_type(0));
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? a[i] : 0;
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }
    }
    {
      c = conditional_assign(a < b, value_type(0), b);
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? 0 : b[i];
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }

      c = 0;
      conditional_assign(c, a < b, value_type(0), b);
      for (int i = 0; i < vector_length; ++i) {
        const auto cc = a[i] < b[i] ? 0 : b[i];
        EXPECT_NEAR_KK(c[i], cc, eps * c[i]);
      }
    }

    {
      typedef Vector<SIMD<bool>, VectorLength> vector_bool_type;
      vector_bool_type cond_all_true, cond_all_false, cond_alternate;

      for (int i = 0; i < vector_length; ++i) {
        cond_all_true[i]  = true;
        cond_all_false[i] = false;
        cond_alternate[i] = i % 2;
      }
      bool all_true, any_true;

      all_true = true;
      any_true = false;
      for (int i = 0; i < vector_length; ++i) {
        all_true &= cond_all_true[i];
        any_true |= cond_all_true[i];
      }
      EXPECT_EQ(all_true, true);
      EXPECT_EQ(any_true, true);

      all_true = true;
      any_true = false;
      for (int i = 0; i < vector_length; ++i) {
        all_true &= cond_all_false[i];
        any_true |= cond_all_false[i];
      }
      EXPECT_EQ(all_true, false);
      EXPECT_EQ(any_true, false);

      all_true = true;
      any_true = false;
      for (int i = 0; i < vector_length; ++i) {
        all_true &= cond_alternate[i];
        any_true |= cond_alternate[i];
      }
      EXPECT_EQ(all_true, false);
      EXPECT_EQ(any_true, true);
    }
    {
      value_type min_a = a[0], max_a = a[0], sum_a = 0, prod_a = 1;
      for (int i = 0; i < vector_length; ++i) {
        min_a = min(min_a, a[i]);
        max_a = max(max_a, a[i]);
        sum_a += a[i];
        prod_a *= a[i];
      }
      EXPECT_NEAR_KK(min_a, min(a), eps * min_a);
      EXPECT_NEAR_KK(max_a, max(a), eps * max_a);
      EXPECT_NEAR_KK(sum_a, sum(a), eps * sum_a);
      EXPECT_NEAR_KK(prod_a, prod(a), eps * prod_a);
    }  // end test body
  }    // end for
}  // impl
}  // namespace Test

template <typename DeviceType, typename VectorTagType, int VectorLength>
int test_batched_vector_misc() {
  static_assert(Kokkos::SpaceAccessibility<DeviceType, Kokkos::HostSpace>::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_misc<VectorTagType, VectorLength>();

  return 0;
}

///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_vector_misc_simd_float3) { test_batched_vector_misc<TestDevice, SIMD<float>, 3>(); }
TEST_F(TestCategory, batched_vector_misc_simd_float8) { test_batched_vector_misc<TestDevice, SIMD<float>, 8>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_vector_misc_simd_double3) { test_batched_vector_misc<TestDevice, SIMD<double>, 3>(); }
TEST_F(TestCategory, batched_vector_misc_simd_double4) { test_batched_vector_misc<TestDevice, SIMD<double>, 4>(); }
#endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// TEST_F( TestCategory, batched_vector_misc_simd_scomplex3 ) {
//   test_batched_vector_misc<TestDevice,SIMD<Kokkos::complex<float> >,3>();
// }
// TEST_F( TestCategory, batched_vector_misc_simd_scomplex4 ) {
//   test_batched_vector_misc<TestDevice,SIMD<Kokkos::complex<float> >,4>();
// }
// #endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_misc_simd_dcomplex3 ) {
//   test_batched_vector_misc<TestDevice,SIMD<Kokkos::complex<double> >,3>();
// }
// TEST_F( TestCategory, batched_vector_misc_simd_dcomplex2 ) {
//   test_batched_vector_misc<TestDevice,SIMD<Kokkos::complex<double> >,2>();
// }
// #endif

#endif  // check to not include this in a device test
