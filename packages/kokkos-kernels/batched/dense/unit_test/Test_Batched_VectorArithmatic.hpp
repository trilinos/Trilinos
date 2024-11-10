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
void impl_test_complex_real_imag_value() {
  typedef Vector<VectorTagType, VectorLength> vector_type;

  typedef typename vector_type::value_type value_type;
  const int vector_length = vector_type::vector_length;

  vector_type a;

  for (int k = 0; k < vector_length; ++k) {
    a[k].real() = k * 3 + 1;
    a[k].imag() = k * 5 + 4;
  }

  const auto a_real = Kokkos::ArithTraits<vector_type>::real(a);
  const auto a_imag = Kokkos::ArithTraits<vector_type>::imag(a);

  typedef Kokkos::ArithTraits<value_type> ats;
  const typename ats::mag_type eps = 1.0e3 * ats::epsilon();
  for (int k = 0; k < vector_length; ++k) {
    EXPECT_NEAR(a[k].real(), a_real[k], eps);
    EXPECT_NEAR(a[k].imag(), a_imag[k], eps);
  }
}

template <typename VectorTagType, int VectorLength>
void impl_test_batched_vector_arithmatic() {
  /// random data initialization
  typedef Vector<VectorTagType, VectorLength> vector_type;

  typedef typename vector_type::value_type value_type;
  const int vector_length = vector_type::vector_length;

  typedef Kokkos::ArithTraits<value_type> ats;
  typedef typename ats::mag_type mag_type;

  vector_type a, b, c;
  value_type alpha;
  mag_type beta;
  const value_type zero(0);

  Random<value_type> a_random;
  Random<mag_type> b_random;
  for (int iter = 0; iter < 100; ++iter) {
    for (int k = 0; k < vector_length; ++k) {
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
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] + b[k]), eps * ats::abs(c[k]));

      /// test : value + vec
      c = alpha + b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(alpha + b[k]), eps * ats::abs(c[k]));

      /// test : vec + value
      c = b + alpha;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(b[k] + alpha), eps * ats::abs(c[k]));

      /// test : vec + mag
      c = a + beta;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] + beta), eps * ats::abs(c[k]));

      /// test : mag + vec
      c = beta + a;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(beta + a[k]), eps * ats::abs(c[k]));
    }
    {
      /// test : vec - vec
      c = a - b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] - b[k]), eps * ats::abs(c[k]));

      /// test : value - vec
      c = alpha - b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(alpha - b[k]), eps * ats::abs(c[k]));

      /// test : vec + value
      c = b - alpha;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(b[k] - alpha), eps * ats::abs(c[k]));

      /// test : vec - mag
      c = a - beta;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] - beta), eps * ats::abs(c[k]));

      /// test : mag - vec
      c = beta - a;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(beta - a[k]), eps * ats::abs(c[k]));
    }
    {
      /// test : vec * vec
      c = a * b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] * b[k]), eps * ats::abs(c[k]));

      /// test : value * vec
      c = alpha * b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(alpha * b[k]), eps * ats::abs(c[k]));

      /// test : vec + value
      c = b * alpha;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(b[k] * alpha), eps * ats::abs(c[k]));

      /// test : vec * mag
      c = a * beta;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] * beta), eps * ats::abs(c[k]));

      /// test : mag * vec
      c = beta * a;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(beta * a[k]), eps * ats::abs(c[k]));
    }
    {
      /// test : vec / vec
      c = a / b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] / b[k]), eps * ats::abs(c[k]));

      /// test : value / vec
      c = alpha / b;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(alpha / b[k]), eps * ats::abs(c[k]));

      /// test : vec / value
      c = b / alpha;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(b[k] / alpha), eps * ats::abs(c[k]));

      /// test : mag / vec
      c = beta / a;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(beta / a[k]), eps * ats::abs(c[k]));

      /// test : vec / value
      c = a / beta;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] / beta), eps * ats::abs(c[k]));
    }
    {
      /// test : vec  -vec
      c = -a;
      for (int k = 0; k < vector_length; ++k) EXPECT_NEAR(ats::abs(c[k]), ats::abs(-a[k]), eps * ats::abs(c[k]));
    }
#if defined(__DO_NOT_TEST__)
    {
      /// test : add radial
      const mag_type tiny = 1.0;

      c = vector_type(0);
      c += -vector_type(tiny) * vector_type(a < 0);
      c += vector_type(tiny) * vector_type(a >= 0);

      for (int k = 0; k < vector_length; ++k)
        EXPECT_NEAR(ats::abs(c[k]), ats::abs(a[k] < 0 ? -tiny : tiny), eps * ats::abs(c[k]));
    }
#endif
  }
}
}  // namespace Test

template <typename DeviceType, typename VectorTagType, int VectorLength>
int test_batched_vector_arithmatic() {
  static_assert(Kokkos::SpaceAccessibility<DeviceType, Kokkos::HostSpace>::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_arithmatic<VectorTagType, VectorLength>();

  return 0;
}
template <typename DeviceType, typename VectorTagType, int VectorLength>
int test_batched_complex_real_imag_value() {
  static_assert(Kokkos::SpaceAccessibility<DeviceType, Kokkos::HostSpace>::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_complex_real_imag_value<VectorTagType, VectorLength>();

  return 0;
}

///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_vector_arithmatic_simd_float3) {
  test_batched_vector_arithmatic<TestDevice, SIMD<float>, 3>();
}
TEST_F(TestCategory, batched_vector_arithmatic_simd_float4) {
  test_batched_vector_arithmatic<TestDevice, SIMD<float>, 4>();
}
// avx
TEST_F(TestCategory, batched_vector_arithmatic_simd_float8) {
  test_batched_vector_arithmatic<TestDevice, SIMD<float>, 8>();
}
// avx 512
TEST_F(TestCategory, batched_vector_arithmatic_simd_float16) {
  test_batched_vector_arithmatic<TestDevice, SIMD<float>, 16>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_vector_arithmatic_simd_double3) {
  test_batched_vector_arithmatic<TestDevice, SIMD<double>, 3>();
}
// avx
TEST_F(TestCategory, batched_vector_arithmatic_simd_double4) {
  test_batched_vector_arithmatic<TestDevice, SIMD<double>, 4>();
}
// avx 512
TEST_F(TestCategory, batched_vector_arithmatic_simd_double8) {
  test_batched_vector_arithmatic<TestDevice, SIMD<double>, 8>();
}
#endif

#define __DO_NOT_TEST__
#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, batched_vector_arithmatic_simd_scomplex3) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<float> >, 3>();
}
// avx
TEST_F(TestCategory, batched_vector_arithmatic_simd_scomplex4) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<float> >, 4>();
}
// avx 512
TEST_F(TestCategory, batched_vector_arithmatic_simd_scomplex8) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<float> >, 8>();
}

TEST_F(TestCategory, batched_vector_scomplex_real_imag_value3) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<float> >, 3>();
}
// avx
TEST_F(TestCategory, batched_vector_scomplex_real_imag_value2) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<float> >, 2>();
}
// avx 512
TEST_F(TestCategory, batched_vector_scomplex_real_imag_value4) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<float> >, 4>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_vector_arithmatic_simd_dcomplex3) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<double> >, 3>();
}
// avx
TEST_F(TestCategory, batched_vector_arithmatic_simd_dcomplex2) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<double> >, 2>();
}
// avx 512
TEST_F(TestCategory, batched_vector_arithmatic_simd_dcomplex4) {
  test_batched_vector_arithmatic<TestDevice, SIMD<Kokkos::complex<double> >, 4>();
}

TEST_F(TestCategory, batched_vector_dcomplex_real_imag_value3) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<double> >, 3>();
}
// avx
TEST_F(TestCategory, batched_vector_dcomplex_real_imag_value2) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<double> >, 2>();
}
// avx 512
TEST_F(TestCategory, batched_vector_dcomplex_real_imag_value4) {
  test_batched_complex_real_imag_value<TestDevice, SIMD<Kokkos::complex<double> >, 4>();
}
#endif
#undef __DO_NOT_TEST__

#endif  // check to not include this in a device test
