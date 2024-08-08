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
// Note: Luc Berger-Vergiat 04/14/21
//       This tests uses KOKKOS_LAMBDA so we need
//       to make sure that these are enabled in
//       the CUDA backend before including this test.
#if !defined(TEST_CUDA_BLAS_CPP) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)

#include <KokkosBlas_util.hpp>
#include <KokkosKernels_TestUtils.hpp>  // for test/inst guards
// Note: include serial gemv before util so it knows if CompactMKL is available
#include <Test_Blas2_gemv_util.hpp>
#include <KokkosBlas2_gemv.hpp>

namespace Test {

template <class AType, class XType, class YType, class ScalarType, class AlgoTag>
struct TeamVectorGEMVOp : public GemvOpBase<AType, XType, YType, ScalarType> {
  using params = GemvOpBase<AType, XType, YType, ScalarType>;

  TeamVectorGEMVOp(char trans_, ScalarType alpha_, AType A_, XType x_, ScalarType beta_, YType y_)
      : params(trans_, alpha_, A_, x_, beta_, y_) {}

  template <typename TeamMember>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& member) const {
    KokkosBlas::Experimental::Gemv<KokkosBlas::Mode::TeamVector, AlgoTag>::invoke(
        member, params::trans, params::alpha, params::A, params::x, params::beta, params::y);
  }
};

struct TeamVectorGemvFactory {
  template <class AlgoTag, class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device, class ScalarType>
  using functor_type = TeamVectorGEMVOp<ViewTypeA, ViewTypeX, ViewTypeY, ScalarType, AlgoTag>;

  // no Blocked implementation
  using algorithms = std::tuple<KokkosBlas::Algo::Gemv::Unblocked>;
};

}  // namespace Test

#define TEST_TEAMVECTOR_CASE4(N, A, X, Y, SC) TEST_CASE4(teamvector, TeamVectorGemvFactory, N, A, X, Y, SC)
#define TEST_TEAMVECTOR_CASE2(N, S, SC) TEST_CASE2(teamvector, TeamVectorGemvFactory, N, S, SC)
#define TEST_TEAMVECTOR_CASE(N, S) TEST_CASE(teamvector, TeamVectorGemvFactory, N, S)

#ifdef KOKKOSKERNELS_TEST_FLOAT
TEST_TEAMVECTOR_CASE(float, float)
#endif

#ifdef KOKKOSKERNELS_TEST_DOUBLE
TEST_TEAMVECTOR_CASE(double, double)
#endif

#ifdef KOKKOSKERNELS_TEST_COMPLEX_DOUBLE
TEST_TEAMVECTOR_CASE(complex_double, Kokkos::complex<double>)
#endif

#ifdef KOKKOSKERNELS_TEST_COMPLEX_FLOAT
TEST_TEAMVECTOR_CASE(complex_float, Kokkos::complex<float>)
#endif

#ifdef KOKKOSKERNELS_TEST_INT
TEST_TEAMVECTOR_CASE(int, int)
#endif

#ifdef KOKKOSKERNELS_TEST_ALL_TYPES
// test mixed scalar types (void -> default alpha/beta)
TEST_TEAMVECTOR_CASE4(mixed, double, int, float, void)

// test arbitrary double alpha/beta with complex<double> values
TEST_TEAMVECTOR_CASE2(alphabeta, Kokkos::complex<double>, double)
#endif

#undef TEST_TEAMVECTOR_CASE4
#undef TEST_TEAMVECTOR_CASE2
#undef TEST_TEAMVECTOR_CASE

#endif  // Check for lambda availability on CUDA backend
