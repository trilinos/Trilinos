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

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_ApplyPivot_Decl.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_UTV_Decl.hpp"
#include "KokkosBatched_SolveUTV_Decl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

template <typename DeviceType, typename MatrixViewType, typename VectorViewType, typename PivViewType,
          typename WorkViewType, typename AlgoTagType>
struct Functor_TestBatchedTeamVectorSolveUTV2 {
  using execution_space = typename DeviceType::execution_space;
  MatrixViewType _r, _a, _acopy, _u, _v;
  PivViewType _p;
  VectorViewType _x, _b;
  WorkViewType _w;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorSolveUTV2(const MatrixViewType &r, const MatrixViewType &a, const MatrixViewType &acopy,
                                         const MatrixViewType &u, const MatrixViewType &v, const PivViewType &p,
                                         const VectorViewType &x, const VectorViewType &b, const WorkViewType &w)
      : _r(r), _a(a), _acopy(acopy), _u(u), _v(v), _p(p), _x(x), _b(b), _w(w) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    typedef typename MatrixViewType::non_const_value_type value_type;
    const value_type one(1), zero(0), add_this(10);

    const int k = member.league_rank();
    auto rr     = Kokkos::subview(_r, k, Kokkos::ALL(), Kokkos::ALL());
    auto ac     = Kokkos::subview(_acopy, k, Kokkos::ALL(), Kokkos::ALL());

    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto uu = Kokkos::subview(_u, k, Kokkos::ALL(), Kokkos::ALL());
    auto vv = Kokkos::subview(_v, k, Kokkos::ALL(), Kokkos::ALL());

    auto pp = Kokkos::subview(_p, k, Kokkos::ALL());

    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto xx = Kokkos::subview(_x, k, Kokkos::ALL(), Kokkos::ALL());

    auto ww = Kokkos::subview(_w, k, Kokkos::ALL());

    // make diagonal dominant and set xx = 1,2,3,4,5
    const int m = aa.extent(0), r = rr.extent(1);
    if (m <= r) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) {
        aa(i, i) += add_this;
        for (int j = 0; j < 2; ++j) xx(i, j) = (i + 1);
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m * m), [=](const int &ij) {
        const int i = ij / m, j = ij % m;
        value_type tmp(0);
        for (int l = 0; l < r; ++l) tmp += rr(i, l) * rr(j, l);
        aa(i, j) = tmp;
      });
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) {
        for (int j = 0; j < 2; ++j) xx(i, j) = (i + 1);
      });
    }
    member.team_barrier();  // finish writing aa, xx

    /// copy for verification
    TeamVectorCopy<MemberType, Trans::NoTranspose>::invoke(member, aa, ac);

    /// bb = AA*xx
    KokkosBatched::TeamVectorGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
        member, one, aa, xx, zero, bb);
    member.team_barrier();

    /// Solving Ax = b using UTV transformation
    /// A P^T P x = b
    /// UTV P x = b;

    /// UTV = A P^T
    int matrix_rank(0);
    TeamVectorUTV<MemberType, AlgoTagType>::invoke(member, aa, pp, uu, vv, ww, matrix_rank);
    member.team_barrier();

    TeamVectorSolveUTV<MemberType, AlgoTagType>::invoke(member, matrix_rank, uu, aa, vv, pp, xx, bb, ww);
  }

  inline void run() {
    typedef typename MatrixViewType::non_const_value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVectorSolveUTV");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    const int league_size = _a.extent(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);

    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename MatrixViewType, typename VectorViewType, typename PivViewType,
          typename WorkViewType, typename AlgoTagType>
void impl_test_batched_solve_utv2(const int N, const int BlkSize) {
  typedef typename MatrixViewType::non_const_value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;
  // const value_type one(1);
  /// randomized input testing views
  MatrixViewType r("r", N, BlkSize, 3);
  MatrixViewType a("a", N, BlkSize, BlkSize);
  MatrixViewType acopy("copy", N, BlkSize, BlkSize);
  MatrixViewType u("u", N, BlkSize, BlkSize);
  MatrixViewType v("v", N, BlkSize, BlkSize);
  PivViewType p("p", N, BlkSize);
  VectorViewType x("x", N, BlkSize, 2);
  VectorViewType b("b", N, BlkSize, 2);
  WorkViewType w("w", N, 3 * BlkSize * 2);

  Kokkos::fence();

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  if (BlkSize <= 3)
    Kokkos::fill_random(a, random, value_type(1.0));
  else
    Kokkos::fill_random(r, random, value_type(1.0));

  Kokkos::fence();

  Functor_TestBatchedTeamVectorSolveUTV2<DeviceType, MatrixViewType, VectorViewType, PivViewType, WorkViewType,
                                         AlgoTagType>(r, a, acopy, u, v, p, x, b, w)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  auto a_host = Kokkos::create_mirror_view(acopy);
  auto x_host = Kokkos::create_mirror_view(x);
  auto b_host = Kokkos::create_mirror_view(b);
  auto w_host = Kokkos::create_mirror_view(w);

  Kokkos::deep_copy(a_host, acopy);
  Kokkos::deep_copy(x_host, x);
  Kokkos::deep_copy(b_host, b);

  /// this is least square; we cannot expect high
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 1e3 * ats::epsilon();

  for (int k = 0; k < N; ++k) {
    mag_type residual(0), norm(0);

    for (int l = 0; l < 2; ++l) {
      for (int i = 0; i < BlkSize; ++i) {
        value_type tmp(0);
        for (int j = 0; j < BlkSize; ++j) {
          tmp += a_host(k, i, j) * x_host(k, j, l);
        }
        w_host(k, i) = tmp - b_host(k, i, l);
      }
      for (int i = 0; i < BlkSize; ++i) {
        value_type tmp(0);
        for (int j = 0; j < BlkSize; ++j) {
          tmp += a_host(k, j, i) * w_host(k, j);
        }
        residual += ats::abs(tmp) * ats::abs(tmp);
        norm += ats::abs(b_host(k, i, l)) * ats::abs(b_host(k, i, l));
      }
    }
    // printf("norm %e, residual %e, rel res %e\n", norm, residual,
    // residual/norm);
    EXPECT_NEAR_KK(residual / norm, mag_type(0), eps);
  }
}
}  // namespace Test

template <typename DeviceType, typename ValueType, typename IntType, typename AlgoTagType>
int test_batched_solve_utv2() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> MatrixViewType;
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> VectorViewType;
    typedef Kokkos::View<IntType **, Kokkos::LayoutLeft, DeviceType> PivViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WorkViewType;
    Test::impl_test_batched_solve_utv2<DeviceType, MatrixViewType, VectorViewType, PivViewType, WorkViewType,
                                       AlgoTagType>(0, 10);
    for (int i = 1; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::impl_test_batched_solve_utv2<DeviceType, MatrixViewType, VectorViewType, PivViewType, WorkViewType,
                                         AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> MatrixViewType;
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> VectorViewType;
    typedef Kokkos::View<IntType **, Kokkos::LayoutRight, DeviceType> PivViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WorkViewType;
    Test::impl_test_batched_solve_utv2<DeviceType, MatrixViewType, VectorViewType, PivViewType, WorkViewType,
                                       AlgoTagType>(0, 10);
    for (int i = 1; i < 10; ++i) {
      // printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::impl_test_batched_solve_utv2<DeviceType, MatrixViewType, VectorViewType, PivViewType, WorkViewType,
                                         AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
