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
#include <KokkosBlas_util.hpp>
#include <KokkosKernels_TestUtils.hpp>  // for ETI test guards
// Note: include serial gemv before util so it knows if CompactMKL is available
#include <KokkosBlas2_gemv.hpp>
#include <Test_Blas2_gemv_util.hpp>

namespace Test {

template <class AType, class XType, class YType, class ScalarType, class AlgoTag>
struct SerialGEMVOp : public GemvOpBase<AType, XType, YType, ScalarType> {
  using params = GemvOpBase<AType, XType, YType, ScalarType>;

  SerialGEMVOp(char trans_, ScalarType alpha_, AType A_, XType x_, ScalarType beta_, YType y_)
      : params(trans_, alpha_, A_, x_, beta_, y_) {}

  template <typename TeamMember>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& member) const {
    KokkosBlas::Experimental::Gemv<KokkosBlas::Mode::Serial, AlgoTag>::invoke(
        member, params::trans, params::alpha, params::A, params::x, params::beta, params::y);
  }
};

struct SerialGemvFactory {
  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device, class ScalarType>
  using functor_type = SerialGEMVOp<ViewTypeA, ViewTypeX, ViewTypeY, ScalarType, AlgoTag>;

  using algorithms = std::tuple<KokkosBlas::Algo::Gemv::Unblocked, KokkosBlas::Algo::Gemv::Blocked>;
};

#ifdef __KOKKOSBLAS_ENABLE_INTEL_MKL_COMPACT__
struct SerialMKLGemvFactory {
  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device, class ScalarType>
  using functor_type = SerialGEMVOp<ViewTypeA, ViewTypeX, ViewTypeY, ScalarType, AlgoTag>;

  using algorithms = std::tuple<KokkosBlas::Algo::Gemv::CompactMKL>;
};
#endif

}  // namespace Test

#define TEST_SERIAL_CASE4(N, A, X, Y, SC) TEST_CASE4(serial, SerialGemvFactory, N, A, X, Y, SC)
#define TEST_SERIAL_CASE2(N, S, SC) TEST_CASE2(serial, SerialGemvFactory, N, S, SC)
#define TEST_SERIAL_CASE(N, S) TEST_CASE(serial, SerialGemvFactory, N, S)

#ifdef KOKKOSKERNELS_TEST_FLOAT
TEST_SERIAL_CASE(float, float)
// MKL vector types
#ifdef __KOKKOSBLAS_ENABLE_INTEL_MKL_COMPACT__
using simd_float_sse    = ::Test::simd_vector<float, 4>;
using simd_float_avx    = ::Test::simd_vector<float, 8>;
using simd_float_avx512 = ::Test::simd_vector<float, 16>;
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_float_sse, simd_float_sse, float)
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_float_avx, simd_float_avx, float)
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_float_avx512, simd_float_avx512, float)
#endif
#endif

#ifdef KOKKOSKERNELS_TEST_DOUBLE
TEST_SERIAL_CASE(double, double)
// MKL vector types
#ifdef __KOKKOSBLAS_ENABLE_INTEL_MKL_COMPACT__
using simd_double_sse    = ::Test::simd_vector<double, 2>;
using simd_double_avx    = ::Test::simd_vector<double, 4>;
using simd_double_avx512 = ::Test::simd_vector<double, 8>;
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_double_sse, simd_double_sse, double)
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_double_avx, simd_double_avx, double)
TEST_CASE2(serial, SerialMKLGemvFactory, mkl_double_avx512, simd_double_avx512, double)
#endif
#endif

#ifdef KOKKOSKERNELS_TEST_COMPLEX_DOUBLE
TEST_SERIAL_CASE(complex_double, Kokkos::complex<double>)
#endif

#ifdef KOKKOSKERNELS_TEST_COMPLEX_FLOAT
TEST_SERIAL_CASE(complex_float, Kokkos::complex<float>)
#endif

#ifdef KOKKOSKERNELS_TEST_INT
TEST_SERIAL_CASE(int, int)
#endif

#ifdef KOKKOSKERNELS_TEST_ALL_TYPES
// test mixed scalar types (void -> default alpha/beta)
TEST_SERIAL_CASE4(mixed, double, int, float, void)

// test arbitrary double alpha/beta with complex<double> values
TEST_SERIAL_CASE2(alphabeta, Kokkos::complex<double>, double)
#endif

#undef TEST_SERIAL_CASE4
#undef TEST_SERIAL_CASE2
#undef TEST_SERIAL_CASE
