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
/// \author Vinh Dang (vqdang@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

// #include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Team_Impl.hpp"
#include "KokkosBatched_InverseLU_Decl.hpp"
// #include "KokkosBatched_InverseLU_Team_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamInverseLU {

template <typename TA, typename TB>
struct ParamTag {
  typedef TA transA;
  typedef TB transB;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_BatchedTeamGemm {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a, _b, _c;

  ScalarType _alpha, _beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedTeamGemm(const ScalarType alpha, const ViewType &a, const ViewType &b, const ScalarType beta,
                          const ViewType &c)
      : _a(a), _b(b), _c(c), _alpha(alpha), _beta(beta) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ParamTagType &, const MemberType &member) const {
    const int k = member.league_rank();

    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), Kokkos::ALL());

    if (member.team_rank() == 0) {
      for (int i = 0; i < static_cast<int>(aa.extent(0)); ++i) aa(i, i) += 10.0;
    }
    member.team_barrier();

    KokkosBatched::TeamGemm<MemberType, typename ParamTagType::transA, typename ParamTagType::transB,
                            AlgoTagType>::invoke(member, _alpha, aa, bb, _beta, cc);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    const int league_size = _c.extent(0);
    Kokkos::TeamPolicy<execution_space, ParamTagType> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for((name + "::GemmFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename AlgoTagType>
struct Functor_BatchedTeamLU {
  ViewType _a;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedTeamLU(const ViewType &a) : _a(a) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int k = member.league_rank();
    auto aa     = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());

    if (member.team_rank() == 0) {
      for (int i = 0; i < static_cast<int>(aa.extent(0)); ++i) aa(i, i) += 10.0;
    }
    member.team_barrier();

    TeamLU<MemberType, AlgoTagType>::invoke(member, aa);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    const int league_size = _a.extent(0);
    Kokkos::TeamPolicy<DeviceType> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for((name + "::LUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename AViewType, typename WViewType, typename AlgoTagType>
struct Functor_TestBatchedTeamInverseLU {
  AViewType _a;
  WViewType _w;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamInverseLU(const AViewType &a, const WViewType &w) : _a(a), _w(w) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int k = member.league_rank();
    auto aa     = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto ww     = Kokkos::subview(_w, k, Kokkos::ALL());

    KokkosBatched::TeamInverseLU<MemberType, AlgoTagType>::invoke(member, aa, ww);
  }

  inline void run() {
    typedef typename AViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    const int league_size = _a.extent(0);
    Kokkos::TeamPolicy<DeviceType> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for((name + "::InverseLUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename AViewType, typename WViewType, typename AlgoTagType>
void impl_test_batched_inverselu(const int N, const int BlkSize) {
  typedef typename AViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  /// randomized input testing views
  AViewType a0("a0", N, BlkSize, BlkSize);
  AViewType a1("a1", N, BlkSize, BlkSize);
  WViewType w("w", N, BlkSize * BlkSize);
  AViewType c0("c0", N, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a1, a0);
  Kokkos::deep_copy(w, value_type(0.0));

  Functor_BatchedTeamLU<DeviceType, AViewType, AlgoTagType>(a1).run();

  Functor_TestBatchedTeamInverseLU<DeviceType, AViewType, WViewType, AlgoTagType>(a1, w).run();

  value_type alpha = 1.0, beta = 0.0;
  typedef ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;

  Functor_BatchedTeamGemm<DeviceType, AViewType, value_type, param_tag_type, AlgoTagType>(alpha, a0, a1, beta, c0)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  typename AViewType::HostMirror c0_host = Kokkos::create_mirror_view(c0);

  Kokkos::deep_copy(c0_host, c0);

  /// check identity matrix ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum_diag(0), sum_all(0), sum_diag_ref(N * BlkSize);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < BlkSize; ++j) {
        sum_all += ats::abs(c0_host(k, i, j));
        if (i == j) sum_diag += ats::abs(c0_host(k, i, j));
      }
  EXPECT_NEAR_KK(sum_all - sum_diag, 0, eps);
  EXPECT_NEAR_KK(sum_diag - sum_diag_ref, 0, eps);
}
}  // namespace TeamInverseLU
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_team_inverselu() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> AViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WViewType;
    Test::TeamInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::TeamInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> AViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WViewType;
    Test::TeamInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::TeamInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
