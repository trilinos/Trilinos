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
void impl_test_batched_vector_math() {
  /// random data initialization
  typedef Vector<VectorTagType, VectorLength> vector_type;

  typedef typename vector_type::value_type value_type;
  const int vector_length = vector_type::vector_length;

  typedef Kokkos::ArithTraits<value_type> ats;
  typedef typename ats::mag_type mag_type;

  vector_type a, b, aref, bref;
  const value_type one(1.0);
  const value_type half(0.5);
  const value_type maxmin(1.0e-7);
  const mag_type eps = 1.0e3 * ats::epsilon();

  Random<value_type> random;
  for (int iter = 0; iter < 100; ++iter) {
    for (int k = 0; k < vector_length; ++k) {
      const auto aval = (random.value() + half);
      const auto bval = (random.value() + half);

      aref[k] = max(min(aval, one), maxmin);
      bref[k] = max(min(bval, one), maxmin);
    }

    {
#undef CHECK
#define CHECK(op)                                                                              \
  {                                                                                            \
    a = op(aref);                                                                              \
    for (int i = 0; i < vector_length; ++i) EXPECT_NEAR_KK(a[i], ats::op(aref[i]), eps* a[i]); \
  }

      CHECK(sqrt);
      CHECK(cbrt);
      CHECK(log);
      CHECK(exp);
      CHECK(sin);
      CHECK(cos);
      CHECK(tan);
      CHECK(sinh);
      CHECK(cosh);
      CHECK(tanh);
      CHECK(asin);
      CHECK(acos);
      CHECK(atan);

#undef CHECK
#define CHECK                                                                                            \
  {                                                                                                      \
    a = pow(aref, bref);                                                                                 \
    for (int i = 0; i < vector_length; ++i) EXPECT_NEAR_KK(a[i], ats::pow(aref[i], bref[i]), eps* a[i]); \
  }                                                                                                      \
  CHECK;

#undef CHECK
#define CHECK(op)                                                                                    \
  {                                                                                                  \
    mag_type beta = mag_type(3.2);                                                                   \
    a             = op(aref, beta);                                                                  \
    for (int i = 0; i < vector_length; ++i) EXPECT_NEAR_KK(a[i], ats::op(aref[i], beta), eps* a[i]); \
  }

      CHECK(pow);

#undef CHECK
#define CHECK(op)                                                                                     \
  {                                                                                                   \
    value_type alpha = random.value() + 2.0;                                                          \
    a                = op(alpha, bref);                                                               \
    for (int i = 0; i < vector_length; ++i) EXPECT_NEAR_KK(a[i], ats::op(alpha, bref[i]), eps* a[i]); \
  }

      CHECK(pow);
#undef CHECK
    }  // end test body
  }    // end for
}  // impl
}  // namespace Test

template <typename DeviceType, typename VectorTagType, int VectorLength>
int test_batched_vector_math() {
  static_assert(Kokkos::SpaceAccessibility<DeviceType, Kokkos::HostSpace>::accessible,
                "vector datatype is only tested on host space");
  Test::impl_test_batched_vector_math<VectorTagType, VectorLength>();

  return 0;
}

// template<typename ValueType>
// int test_complex_pow() {
//   typedef Kokkos::ArithTraits<Kokkos::complex<ValueType> > ats;
//   typedef typename ats::mag_type mag_type;

//   const mag_type eps = 1.0e3 * ats::epsilon();

//   mag_type a2 = 4.5;
//   Kokkos::complex<mag_type> a0(3.2, -1.4), a1(1.2, 2.3);
//   std::complex<mag_type> b0(3.2, -1.4), b1(1.2, 2.3);

//   Test::EXPECT_NEAR_KK( ats::pow(a0,a1), std::pow(b0,b1), eps );
//   Test::EXPECT_NEAR_KK( ats::pow(a0,a2), std::pow(b0,a2), eps );

//   return 0;
// }

///
/// SIMD
///

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_vector_math_simd_float3) { test_batched_vector_math<TestDevice, SIMD<float>, 3>(); }
TEST_F(TestCategory, batched_vector_math_simd_float8) { test_batched_vector_math<TestDevice, SIMD<float>, 8>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_vector_math_simd_double3) { test_batched_vector_math<TestDevice, SIMD<double>, 3>(); }
TEST_F(TestCategory, batched_vector_math_simd_double4) { test_batched_vector_math<TestDevice, SIMD<double>, 4>(); }
#endif

// using namespace Test;

// #if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// TEST_F( TestCategory, batched_vector_math_simd_scomplex3 ) {
//   test_complex_pow<float>();
//   test_batched_vector_math<TestDevice,SIMD<Kokkos::complex<float> >,3>();
// }
// TEST_F( TestCategory, batched_vector_math_simd_scomplex4 ) {
//   test_batched_vector_math<TestDevice,SIMD<Kokkos::complex<float> >,4>();
// }
// #endif

// #if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// TEST_F( TestCategory, batched_vector_math_simd_dcomplex3 ) {
//   test_complex_pow<double>();
//   test_batched_vector_math<TestDevice,SIMD<Kokkos::complex<double> >,3>();
// }
// TEST_F( TestCategory, batched_vector_math_simd_dcomplex2 ) {
//   test_batched_vector_math<TestDevice,SIMD<Kokkos::complex<double> >,2>();
// }
// #endif

#endif  // check to not include this in a device test
