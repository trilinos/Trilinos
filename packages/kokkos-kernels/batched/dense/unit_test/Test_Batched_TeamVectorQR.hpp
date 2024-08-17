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

#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBlas2_team_gemv_spec.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_QR_Decl.hpp"
#include "KokkosBatched_ApplyQ_Decl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

template <typename DeviceType, typename MatrixViewType, typename VectorViewType, typename WorkViewType,
          typename AlgoTagType>
struct Functor_TestBatchedTeamVectorQR {
  using execution_space = typename DeviceType::execution_space;
  MatrixViewType _a;
  VectorViewType _x, _b, _t;
  WorkViewType _w;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorQR(const MatrixViewType &a, const VectorViewType &x, const VectorViewType &b,
                                  const VectorViewType &t, const WorkViewType &w)
      : _a(a), _x(x), _b(b), _t(t), _w(w) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    typedef typename MatrixViewType::non_const_value_type value_type;
    const value_type one(1), zero(0), add_this(10);

    const int k = member.league_rank();
    auto aa     = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb     = Kokkos::subview(_b, k, Kokkos::ALL());
    auto xx     = Kokkos::subview(_x, k, Kokkos::ALL());
    auto tt     = Kokkos::subview(_t, k, Kokkos::ALL());
    auto ww     = Kokkos::subview(_w, k, Kokkos::ALL());

    // make diagonal dominant
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, aa.extent(0)), [&](const int &i) { aa(i, i) += add_this; });

    /// xx = 1
    KokkosBlas::TeamVectorSet<MemberType>::invoke(member, one, xx);
    member.team_barrier();

    /// bb = AA*xx
    KokkosBlas::TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(member, one, aa, xx, zero,
                                                                                              bb);
    member.team_barrier();

    /// AA = QR
    TeamVectorQR<MemberType, AlgoTagType>::invoke(member, aa, tt, ww);
    member.team_barrier();

    /// xx = bb;
    TeamVectorCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, bb, xx);
    member.team_barrier();

    /// xx = Q^{T}xx;
    TeamVectorApplyQ<MemberType, Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(member, aa, tt, xx, ww);
    member.team_barrier();

    /// xx = R^{-1} xx
    TeamVectorTrsv<MemberType, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Trsv::Unblocked>::invoke(
        member, one, aa, xx);
  }

  inline void run() {
    typedef typename MatrixViewType::non_const_value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVectorQR");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    const int league_size = _a.extent(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);

    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename MatrixViewType, typename VectorViewType, typename WorkViewType,
          typename AlgoTagType>
void impl_test_batched_qr(const int N, const int BlkSize) {
  typedef typename MatrixViewType::non_const_value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;
  const value_type one(1);
  /// randomized input testing views
  MatrixViewType a("a", N, BlkSize, BlkSize);
  VectorViewType x("x", N, BlkSize);
  VectorViewType b("b", N, BlkSize);
  VectorViewType t("t", N, BlkSize);
  WorkViewType w("w", N, BlkSize);

  Kokkos::fence();

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a, random, value_type(1.0));

  Kokkos::fence();

  Functor_TestBatchedTeamVectorQR<DeviceType, MatrixViewType, VectorViewType, WorkViewType, AlgoTagType>(a, x, b, t, w)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  typename VectorViewType::HostMirror x_host = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(x_host, x);

  /// check x = 1; this eps is about 1e-14
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 1e3 * ats::epsilon();

  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < BlkSize; ++i) {
      const mag_type sum  = ats::abs(x_host(k, i));
      const mag_type diff = ats::abs(x_host(k, i) - one);
      EXPECT_NEAR_KK(diff / sum, mag_type(0), eps);
      // printf("k = %d, i = %d, sum %e diff %e \n", k, i, sum, diff );
    }
  }
}
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_qr() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> MatrixViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> VectorViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WorkViewType;
    Test::impl_test_batched_qr<DeviceType, MatrixViewType, VectorViewType, WorkViewType, AlgoTagType>(0, 10);
    for (int i = 1; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::impl_test_batched_qr<DeviceType, MatrixViewType, VectorViewType, WorkViewType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> MatrixViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> VectorViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WorkViewType;
    Test::impl_test_batched_qr<DeviceType, MatrixViewType, VectorViewType, WorkViewType, AlgoTagType>(0, 10);
    for (int i = 1; i < 10; ++i) {
      // printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::impl_test_batched_qr<DeviceType, MatrixViewType, VectorViewType, WorkViewType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
