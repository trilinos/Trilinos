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
#include "KokkosBatched_GMRES.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosBatched_CrsMatrix.hpp"
#include "Test_Batched_SparseUtils.hpp"
#include "KokkosBatched_JacobiPrec.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamVectorGMRES {

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType,
          typename KrylovHandleType>
struct Functor_TestBatchedTeamVectorGMRES {
  using execution_space = typename DeviceType::execution_space;
  const ValuesViewType _D;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;
  const VectorViewType _Diag;
  const int _N_team;
  KrylovHandleType _handle;

  Functor_TestBatchedTeamVectorGMRES(const ValuesViewType &D, const IntView &r, const IntView &c,
                                     const VectorViewType &X, const VectorViewType &B, const VectorViewType &diag,
                                     const int N_team, KrylovHandleType &handle)
      : _D(D), _r(r), _c(c), _X(X), _B(B), _Diag(diag), _N_team(N_team), _handle(handle) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int first_matrix = static_cast<int>(member.league_rank()) * _N_team;
    const int N            = _D.extent(0);
    const int last_matrix =
        (static_cast<int>(member.league_rank() + 1) * _N_team < N ? static_cast<int>(member.league_rank() + 1) * _N_team
                                                                  : N);

    auto d    = Kokkos::subview(_D, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto diag = Kokkos::subview(_Diag, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto x    = Kokkos::subview(_X, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto b    = Kokkos::subview(_B, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    using Operator     = KokkosBatched::CrsMatrix<ValuesViewType, IntView>;
    using PrecOperator = KokkosBatched::JacobiPrec<ValuesViewType>;

    Operator A(d, _r, _c);
    PrecOperator P(diag);
    P.setComputedInverse();

    KokkosBatched::TeamVectorGMRES<MemberType>::template invoke<Operator, VectorViewType>(member, A, b, x, P, _handle);
  }

  inline void run() {
    typedef typename ValuesViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVectorGMRES");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::TeamPolicy<execution_space> policy(_D.extent(0) / _N_team, Kokkos::AUTO(), Kokkos::AUTO());

    const int N                 = _D.extent(0);
    const int n                 = _X.extent(1);
    const int maximum_iteration = _handle.get_max_iteration();

    _handle.set_ortho_strategy(0);
    _handle.set_compute_last_residual(false);
    _handle.set_tolerance(1e-8);

    _handle.Arnoldi_view =
        typename KrylovHandleType::ArnoldiViewType("", N, maximum_iteration, n + maximum_iteration + 3);

    using ScalarType = typename ValuesViewType::non_const_value_type;
    using Layout     = typename ValuesViewType::array_layout;
    using EXSP       = typename ValuesViewType::execution_space;

    using ViewType2D = Kokkos::View<ScalarType **, Layout, EXSP>;

    size_t bytes_1D   = ViewType2D::shmem_size(_N_team, 1);
    size_t bytes_2D_1 = ViewType2D::shmem_size(_N_team, _X.extent(1));
    size_t bytes_2D_2 = ViewType2D::shmem_size(_N_team, maximum_iteration + 1);

    size_t bytes_row_ptr = IntView::shmem_size(_r.extent(0));
    size_t bytes_col_idc = IntView::shmem_size(_c.extent(0));

    size_t bytes_int  = bytes_row_ptr + bytes_col_idc;
    size_t bytes_diag = bytes_2D_1;
    size_t bytes_tmp  = 2 * bytes_2D_1 + 2 * bytes_1D + bytes_2D_2;
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_tmp + bytes_diag + bytes_int));

    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType>
void impl_test_batched_GMRES(const int N, const int BlkSize, const int N_team) {
  typedef typename ValuesViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  const int nnz = (BlkSize - 2) * 3 + 2 * 2;

  VectorViewType X("x0", N, BlkSize);
  VectorViewType R("r0", N, BlkSize);
  VectorViewType B("b", N, BlkSize);
  ValuesViewType D("D", N, nnz);
  ValuesViewType Diag("Diag", N, BlkSize);
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

  {
    auto diag_values_host = Kokkos::create_mirror_view(Diag);
    auto values_host      = Kokkos::create_mirror_view(D);
    auto row_ptr_host     = Kokkos::create_mirror_view(r);
    auto colIndices_host  = Kokkos::create_mirror_view(c);

    Kokkos::deep_copy(values_host, D);
    Kokkos::deep_copy(row_ptr_host, r);
    Kokkos::deep_copy(colIndices_host, c);

    int current_index;
    for (int i = 0; i < BlkSize; ++i) {
      for (current_index = row_ptr_host(i); current_index < row_ptr_host(i + 1); ++current_index) {
        if (colIndices_host(current_index) == i) break;
      }
      for (int j = 0; j < N; ++j) diag_values_host(j, i) = values_host(j, current_index);
    }

    Kokkos::deep_copy(Diag, diag_values_host);
  }

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

  const int n_iterations = 10;
  KrylovHandleType handle(N, N_team, n_iterations);

  KokkosBatched::SerialSpmv<Trans::NoTranspose>::template invoke<
      typename ValuesViewType::HostMirror, typename IntView::HostMirror, typename VectorViewType::HostMirror,
      typename VectorViewType::HostMirror, 1>(-1, D_host, r_host, c_host, X_host, 1, R_host);
  KokkosBatched::SerialDot<Trans::NoTranspose>::invoke(R_host, R_host, sqr_norm_0_host);
  Functor_TestBatchedTeamVectorGMRES<DeviceType, ValuesViewType, IntView, VectorViewType, KrylovHandleType>(
      D, r, c, X, B, Diag, N_team, handle)
      .run();

  Kokkos::fence();

  Kokkos::deep_copy(R, B);
  Kokkos::deep_copy(R_host, R);
  Kokkos::deep_copy(X_host, X);

  KokkosBatched::SerialSpmv<Trans::NoTranspose>::template invoke<
      typename ValuesViewType::HostMirror, typename IntView::HostMirror, typename VectorViewType::HostMirror,
      typename VectorViewType::HostMirror, 1>(-1, D_host, r_host, c_host, X_host, 1, R_host);
  KokkosBatched::SerialDot<Trans::NoTranspose>::invoke(R_host, R_host, sqr_norm_j_host);

  const MagnitudeType eps = 1.0e5 * ats::epsilon();

  for (int l = 0; l < N; ++l) EXPECT_NEAR_KK(std::sqrt(sqr_norm_j_host(l)) / std::sqrt(sqr_norm_0_host(l)), 0, eps);
}
}  // namespace TeamVectorGMRES
}  // namespace Test

template <typename DeviceType, typename ValueType>
int test_batched_teamvector_GMRES() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutLeft, DeviceType> IntView;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> VectorViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorGMRES::impl_test_batched_GMRES<DeviceType, ViewType, IntView, VectorViewType>(1024, i, 2);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> ViewType;
    typedef Kokkos::View<int *, Kokkos::LayoutRight, DeviceType> IntView;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> VectorViewType;

    for (int i = 3; i < 10; ++i) {
      Test::TeamVectorGMRES::impl_test_batched_GMRES<DeviceType, ViewType, IntView, VectorViewType>(1024, i, 2);
    }
  }
#endif

  return 0;
}
