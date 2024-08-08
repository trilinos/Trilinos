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

// #include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_Spmv_TeamVector_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

#include "Test_Batched_SparseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamVectorSpmv {

template <typename T>
struct ParamTag {
  typedef T trans;
};

template <typename DeviceType, typename ParamTagType, typename ValuesViewType, typename IntView, typename xViewType,
          typename yViewType, typename alphaViewType, typename betaViewType, int dobeta>
struct Functor_TestBatchedTeamVectorSpmv {
  using execution_space = typename DeviceType::execution_space;
  const alphaViewType _alpha;
  const ValuesViewType _D;
  const IntView _r;
  const IntView _c;
  const xViewType _X;
  const betaViewType _beta;
  const yViewType _Y;
  const int _N_team;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorSpmv(const alphaViewType &alpha, const ValuesViewType &D, const IntView &r,
                                    const IntView &c, const xViewType &X, const betaViewType &beta, const yViewType &Y,
                                    const int N_team)
      : _alpha(alpha), _D(D), _r(r), _c(c), _X(X), _beta(beta), _Y(Y), _N_team(N_team) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ParamTagType &, const MemberType &member) const {
    const int first_matrix = static_cast<int>(member.league_rank()) * _N_team;
    const int N            = _D.extent(0);
    const int last_matrix =
        (static_cast<int>(member.league_rank() + 1) * _N_team < N ? static_cast<int>(member.league_rank() + 1) * _N_team
                                                                  : N);

    auto alpha = Kokkos::subview(_alpha, Kokkos::make_pair(first_matrix, last_matrix));
    auto d     = Kokkos::subview(_D, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto x     = Kokkos::subview(_X, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto beta  = Kokkos::subview(_beta, Kokkos::make_pair(first_matrix, last_matrix));
    auto y     = Kokkos::subview(_Y, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    if (last_matrix != N)
      KokkosBatched::TeamVectorSpmv<MemberType, typename ParamTagType::trans, 2>::template invoke<
          ValuesViewType, IntView, xViewType, yViewType, alphaViewType, betaViewType, dobeta>(member, alpha, d, _r, _c,
                                                                                              x, beta, y);
    else
      KokkosBatched::TeamVectorSpmv<MemberType, typename ParamTagType::trans>::template invoke<
          ValuesViewType, IntView, xViewType, yViewType, alphaViewType, betaViewType, dobeta>(member, alpha, d, _r, _c,
                                                                                              x, beta, y);
  }

  inline void run() {
    typedef typename ValuesViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVectorSpmv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::TeamPolicy<execution_space, ParamTagType> policy(ceil(static_cast<double>(_D.extent(0)) / _N_team),
                                                             Kokkos::AUTO(), Kokkos::AUTO());
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ParamTagType, typename ValuesViewType, typename IntView, typename xViewType,
          typename yViewType, typename alphaViewType, typename betaViewType, int dobeta>
void impl_test_batched_spmv(const int N, const int BlkSize, const int N_team) {
  typedef typename ValuesViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  const int nnz = (BlkSize - 2) * 3 + 2 * 2;

  xViewType X0("x0", N, BlkSize), X1("x1", N, BlkSize);
  yViewType Y0("y0", N, BlkSize), Y1("y1", N, BlkSize);
  ValuesViewType D("D", N, nnz);
  IntView r("r", BlkSize + 1);
  IntView c("c", nnz);

  alphaViewType alpha("alpha", N);
  betaViewType beta("beta", N);

  Kokkos::deep_copy(alpha, value_type(1.0));
  Kokkos::deep_copy(beta, value_type(1.0));

  create_tridiagonal_batched_matrices(nnz, BlkSize, N, r, c, D, X0, Y0);

  Kokkos::deep_copy(X1, X0);
  Kokkos::deep_copy(Y1, Y0);

  /// test body
  auto alpha_host = Kokkos::create_mirror_view(alpha);
  auto beta_host  = Kokkos::create_mirror_view(beta);
  auto X0_host    = Kokkos::create_mirror_view(X0);
  auto Y0_host    = Kokkos::create_mirror_view(Y0);

  Kokkos::deep_copy(alpha_host, alpha);
  Kokkos::deep_copy(beta_host, beta);
  Kokkos::deep_copy(X0_host, X0);
  Kokkos::deep_copy(Y0_host, Y0);

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) {
      if (dobeta == 0)
        Y0_host(l, i) = value_type(0.0);
      else
        Y0_host(l, i) *= beta_host(l);
      if (i != 0 && i != (BlkSize - 1))
        Y0_host(l, i) += alpha_host(l) * (2 * X0_host(l, i) - X0_host(l, i - 1) - X0_host(l, i + 1));
      else if (i == 0)
        Y0_host(l, i) += alpha_host(l) * (2 * X0_host(l, i) - X0_host(l, i + 1));
      else
        Y0_host(l, i) += alpha_host(l) * (2 * X0_host(l, i) - X0_host(l, i - 1));
    }

  Functor_TestBatchedTeamVectorSpmv<DeviceType, ParamTagType, ValuesViewType, IntView, xViewType, yViewType,
                                    alphaViewType, betaViewType, dobeta>(alpha, D, r, c, X1, beta, Y1, N_team)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  auto Y1_host = Kokkos::create_mirror_view(Y1);

  Kokkos::deep_copy(Y1_host, Y1);

  /// check c0 = c1 ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) {
      sum += ats::abs(Y0_host(l, i));
      diff += ats::abs(Y0_host(l, i) - Y1_host(l, i));
    }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}
}  // namespace TeamVectorSpmv
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType>
int test_batched_teamvector_spmv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutLeft, DeviceType> IntView;
    typedef Kokkos::View<ScalarType *, Kokkos::LayoutLeft, DeviceType> alphaViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorSpmv::impl_test_batched_spmv<DeviceType, ParamTagType, ViewType, IntView, ViewType, ViewType,
                                                   alphaViewType, alphaViewType, 0>(1025, i, 2);
    }
    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorSpmv::impl_test_batched_spmv<DeviceType, ParamTagType, ViewType, IntView, ViewType, ViewType,
                                                   alphaViewType, alphaViewType, 1>(1025, i, 2);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutRight, DeviceType> IntView;
    typedef Kokkos::View<ScalarType *, Kokkos::LayoutRight, DeviceType> alphaViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorSpmv::impl_test_batched_spmv<DeviceType, ParamTagType, ViewType, IntView, ViewType, ViewType,
                                                   alphaViewType, alphaViewType, 0>(1025, i, 2);
    }
    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorSpmv::impl_test_batched_spmv<DeviceType, ParamTagType, ViewType, IntView, ViewType, ViewType,
                                                   alphaViewType, alphaViewType, 1>(1025, i, 2);
    }
  }
#endif

  return 0;
}
