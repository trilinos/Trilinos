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
#include "KokkosBatched_CG.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosBatched_CrsMatrix.hpp"
#include "Test_Batched_SparseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamVectorCG {

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType,
          typename KrylovHandleType>
struct Functor_TestBatchedTeamVectorCG {
  using execution_space = typename DeviceType::execution_space;
  const ValuesViewType _D;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;
  const int _N_team;
  KrylovHandleType handle;

  Functor_TestBatchedTeamVectorCG(const ValuesViewType &D, const IntView &r, const IntView &c, const VectorViewType &X,
                                  const VectorViewType &B, const int N_team)
      : _D(D), _r(r), _c(c), _X(X), _B(B), _N_team(N_team), handle(KrylovHandleType(_D.extent(0), _N_team)) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int first_matrix = static_cast<int>(member.league_rank()) * _N_team;
    const int N            = _D.extent(0);
    const int last_matrix =
        (static_cast<int>(member.league_rank() + 1) * _N_team < N ? static_cast<int>(member.league_rank() + 1) * _N_team
                                                                  : N);

    auto d = Kokkos::subview(_D, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto x = Kokkos::subview(_X, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto b = Kokkos::subview(_B, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    using Operator = KokkosBatched::CrsMatrix<ValuesViewType, IntView>;

    Operator A(d, _r, _c);

    KokkosBatched::TeamVectorCG<MemberType>::template invoke<Operator, VectorViewType>(member, A, b, x, handle);
  }

  inline void run() {
    typedef typename ValuesViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVectorCG");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::TeamPolicy<execution_space> policy(_D.extent(0) / _N_team, Kokkos::AUTO(), Kokkos::AUTO());

    size_t bytes_0 = ValuesViewType::shmem_size(_N_team, _X.extent(1));
    size_t bytes_1 = ValuesViewType::shmem_size(_N_team, 1);
    policy.set_scratch_size(0, Kokkos::PerTeam(4 * bytes_0 + 5 * bytes_1));

    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType>
void impl_test_batched_CG(const int N, const int BlkSize, const int N_team) {
  typedef typename ValuesViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  const int nnz = (BlkSize - 2) * 3 + 2 * 2;

  VectorViewType X("x0", N, BlkSize);
  VectorViewType R("r0", N, BlkSize);
  VectorViewType B("b", N, BlkSize);
  ValuesViewType D("D", N, nnz);
  IntView r("r", BlkSize + 1);
  IntView c("c", nnz);

  using ScalarType = typename ValuesViewType::non_const_value_type;
  using Layout     = typename ValuesViewType::array_layout;
  using EXSP       = typename ValuesViewType::execution_space;

  using MagnitudeType = typename Kokkos::ArithTraits<ScalarType>::mag_type;
  using NormViewType  = Kokkos::View<MagnitudeType *, Layout, EXSP>;

  using Norm2DViewType   = Kokkos::View<MagnitudeType **, Layout, EXSP>;
  using Scalar3DViewType = Kokkos::View<ScalarType ***, Layout, EXSP>;
  using IntViewType      = Kokkos::View<int *, Layout, EXSP>;

  using KrylovHandleType = KrylovHandle<Norm2DViewType, IntViewType, Scalar3DViewType>;

  NormViewType sqr_norm_0("sqr_norm_0", N);
  NormViewType sqr_norm_j("sqr_norm_j", N);

  create_tridiagonal_batched_matrices(nnz, BlkSize, N, r, c, D, X, B);

  // Compute initial norm

  Kokkos::deep_copy(R, B);

  auto sqr_norm_0_host = Kokkos::create_mirror_view(sqr_norm_0);
  auto sqr_norm_j_host = Kokkos::create_mirror_view(sqr_norm_j);
  auto R_host          = Kokkos::create_mirror_view(R);
  auto X_host          = Kokkos::create_mirror_view(X);
  auto D_host          = Kokkos::create_mirror_view(D);
  auto r_host          = Kokkos::create_mirror_view(r);
  auto c_host          = Kokkos::create_mirror_view(c);

  Kokkos::deep_copy(R, B);
  Kokkos::deep_copy(R_host, R);
  Kokkos::deep_copy(X_host, X);

  Kokkos::deep_copy(c_host, c);
  Kokkos::deep_copy(r_host, r);
  Kokkos::deep_copy(D_host, D);

  KokkosBatched::SerialSpmv<Trans::NoTranspose>::template invoke<
      typename ValuesViewType::HostMirror, typename IntView::HostMirror, typename VectorViewType::HostMirror,
      typename VectorViewType::HostMirror, 1>(-1, D_host, r_host, c_host, X_host, 1, R_host);
  KokkosBatched::SerialDot<Trans::NoTranspose>::invoke(R_host, R_host, sqr_norm_0_host);
  Functor_TestBatchedTeamVectorCG<DeviceType, ValuesViewType, IntView, VectorViewType, KrylovHandleType>(D, r, c, X, B,
                                                                                                         N_team)
      .run();

  Kokkos::fence();

  Kokkos::deep_copy(R, B);
  Kokkos::deep_copy(R_host, R);
  Kokkos::deep_copy(X_host, X);

  KokkosBatched::SerialSpmv<Trans::NoTranspose>::template invoke<
      typename ValuesViewType::HostMirror, typename IntView::HostMirror, typename VectorViewType::HostMirror,
      typename VectorViewType::HostMirror, 1>(-1, D_host, r_host, c_host, X_host, 1, R_host);
  KokkosBatched::SerialDot<Trans::NoTranspose>::invoke(R_host, R_host, sqr_norm_j_host);

  const MagnitudeType eps = 1.0e3 * ats::epsilon();

  for (int l = 0; l < N; ++l) EXPECT_NEAR_KK(sqr_norm_j_host(l) / sqr_norm_0_host(l), 0, eps);
}
}  // namespace TeamVectorCG
}  // namespace Test

template <typename DeviceType, typename ValueType>
int test_batched_teamvector_CG() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutLeft, DeviceType> IntView;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> VectorViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorCG::impl_test_batched_CG<DeviceType, ViewType, IntView, VectorViewType>(1024, i, 2);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutRight, DeviceType> IntView;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> VectorViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorCG::impl_test_batched_CG<DeviceType, ViewType, IntView, VectorViewType>(1024, i, 2);
    }
  }
#endif

  return 0;
}
