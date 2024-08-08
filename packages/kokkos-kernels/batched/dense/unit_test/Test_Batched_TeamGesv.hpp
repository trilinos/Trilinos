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
/// \author Kim Liegeois (knliege@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Gesv.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBlas2_serial_gemv_impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamGesv {

template <typename DeviceType, typename MatrixType, typename VectorType, typename AlgoTagType>
struct Functor_TestBatchedTeamGesv {
  using execution_space = typename DeviceType::execution_space;
  const MatrixType _A;
  const VectorType _X;
  const VectorType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamGesv(const MatrixType &A, const VectorType &X, const VectorType &B) : _A(A), _X(X), _B(B) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int matrix_id = static_cast<int>(member.league_rank());
    auto A              = Kokkos::subview(_A, matrix_id, Kokkos::ALL, Kokkos::ALL);
    auto x              = Kokkos::subview(_X, matrix_id, Kokkos::ALL);
    auto b              = Kokkos::subview(_B, matrix_id, Kokkos::ALL);

    member.team_barrier();
    KokkosBatched::TeamGesv<MemberType, AlgoTagType>::invoke(member, A, x, b);
    member.team_barrier();
  }

  inline void run() {
    typedef typename VectorType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamGesv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::TeamPolicy<execution_space> policy(_X.extent(0), Kokkos::AUTO(), Kokkos::AUTO());

    using MatrixViewType = Kokkos::View<typename MatrixType::non_const_value_type **, typename MatrixType::array_layout,
                                        typename MatrixType::execution_space>;

    const int n    = _A.extent(1);
    size_t bytes_0 = MatrixViewType::shmem_size(n, n + 4);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_0));

    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename MatrixType, typename VectorType, typename AlgoTagType>
void impl_test_batched_gesv(const int N, const int BlkSize) {
  typedef typename MatrixType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  using MagnitudeType = typename Kokkos::ArithTraits<value_type>::mag_type;
  using NormViewType  = Kokkos::View<MagnitudeType *, Kokkos::LayoutLeft, DeviceType>;

  NormViewType sqr_norm_j("sqr_norm_j", N);
  auto sqr_norm_j_host = Kokkos::create_mirror_view(sqr_norm_j);

  MatrixType A("A", N, BlkSize, BlkSize), A2("A", N, BlkSize, BlkSize);
  VectorType B("b", N, BlkSize), B2("b", N, BlkSize), X("x", N, BlkSize);

  create_tridiagonal_batched_matrices(A, B);
  Kokkos::deep_copy(A2, A);
  Kokkos::deep_copy(B2, B);

  auto A_host = Kokkos::create_mirror_view(A2);
  auto B_host = Kokkos::create_mirror_view(B2);
  auto X_host = Kokkos::create_mirror_view(X);

  Kokkos::deep_copy(A_host, A2);
  Kokkos::deep_copy(B_host, B2);

  Kokkos::fence();

  Functor_TestBatchedTeamGesv<DeviceType, MatrixType, VectorType, AlgoTagType>(A, X, B).run();

  Kokkos::fence();

  Kokkos::deep_copy(X_host, X);

  for (int l = 0; l < N; ++l)
    KokkosBlas::SerialGemv<Trans::NoTranspose, KokkosBlas::Algo::Gemv::Unblocked>::invoke(
        -1, Kokkos::subview(A_host, l, Kokkos::ALL, Kokkos::ALL), Kokkos::subview(X_host, l, Kokkos::ALL), 1,
        Kokkos::subview(B_host, l, Kokkos::ALL));

  KokkosBatched::SerialDot<Trans::NoTranspose>::invoke(B_host, B_host, sqr_norm_j_host);

  const MagnitudeType eps = 1.0e3 * ats::epsilon();

  for (int l = 0; l < N; ++l) EXPECT_NEAR_KK(sqr_norm_j_host(l), 0, eps);
}
}  // namespace TeamGesv
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_team_gesv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> MatrixType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> VectorType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamGesv::impl_test_batched_gesv<DeviceType, MatrixType, VectorType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> MatrixType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> VectorType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamGesv::impl_test_batched_gesv<DeviceType, MatrixType, VectorType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
