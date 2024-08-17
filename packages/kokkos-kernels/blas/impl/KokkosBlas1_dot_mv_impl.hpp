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
#ifndef KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_
#define KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <KokkosBlas_util.hpp>
#include <sstream>

namespace KokkosBlas {
namespace Impl {

template <class ExecSpace, class RV, class XV, class YV, class size_type>
struct Dot_MV_Functor {
  using Scalar   = typename RV::non_const_value_type;
  using IPT      = Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>;
  using dot_type = typename IPT::dot_type;
  using KAT      = Kokkos::ArithTraits<dot_type>;

  using TeamMem = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  RV r;
  XV x;
  YV y;

  size_type teamsPerDot;  // number of teams collectively performing a dot product

  Dot_MV_Functor(const RV& r_, const XV& x_, const YV& y_, int teamsPerDot_)
      : r(r_), x(x_), y(y_), teamsPerDot(teamsPerDot_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamMem& t) const {
    size_type globalRank = t.league_rank();
    size_type localRank  = globalRank % teamsPerDot;
    size_type i          = globalRank / teamsPerDot;
    size_type xcol       = x.extent(1) == 1 ? 0 : i;
    size_type ycol       = y.extent(1) == 1 ? 0 : i;

    dot_type localResult = KAT::zero();
    size_type begin      = localRank * (x.extent(0) / teamsPerDot);
    size_type end        = (localRank + 1) * (x.extent(0) / teamsPerDot);
    if (localRank == teamsPerDot - 1) end = x.extent(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(t, begin, end),
        [&](size_type k, dot_type& update) {
          Kokkos::Details::updateDot(update, x.access(k, xcol), y.access(k, ycol));
        },
        localResult);

    Kokkos::single(Kokkos::PerTeam(t), [&]() { Kokkos::atomic_add(&r(i), Scalar(localResult)); });
  }
};

// Main version: the result view is accessible from execution space, so it can
// be computed in-place
template <class execution_space, class RV, class XV, class YV, class size_type>
void MV_Dot_Invoke(
    const execution_space& space, const RV& r, const XV& x, const YV& y,
    typename std::enable_if<Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible>::type* =
        nullptr) {
  size_type numDots = std::max(x.extent(1), y.extent(1));
  if (x.extent(0) != y.extent(0)) {
    std::ostringstream oss;
    oss << "KokkosBlas::dot (rank-2): x and y have different lengths (" << x.extent(0) << " and " << y.extent(0) << ")";
    throw std::runtime_error(oss.str());
  }
  if ((x.extent(1) != size_t(1) && x.extent(1) != size_t(numDots)) ||
      (y.extent(1) != size_t(1) && y.extent(1) != size_t(numDots))) {
    std::ostringstream oss;
    oss << "KokkosBlas::dot (rank-2): x and y have incompatible numbers of "
           "columns ("
        << x.extent(1) << " and " << y.extent(1) << ")";
    throw std::runtime_error(oss.str());
  }
  if (r.extent(0) != size_t(numDots)) {
    std::ostringstream oss;
    oss << "KokkosBlas::dot (rank-2): result vector has wrong length (" << r.extent(0) << ", but " << numDots
        << " dot products will be computed)";
    throw std::runtime_error(oss.str());
  }
  // Zero out the result vector
  Kokkos::deep_copy(space, r, Kokkos::ArithTraits<typename RV::non_const_value_type>::zero());
  size_type teamsPerDot;
  KokkosBlas::Impl::multipleReductionWorkDistribution<execution_space, size_type>(x.extent(0), numDots, teamsPerDot);
  size_type numTeams = numDots * teamsPerDot;
  Kokkos::TeamPolicy<execution_space> pol(space, numTeams, Kokkos::AUTO);
  Kokkos::parallel_for("Dot_MV", pol, Dot_MV_Functor<execution_space, RV, XV, YV, size_type>(r, x, y, teamsPerDot));
}

// Version for when a temporary result view is needed (implemented in terms of
// the other version)
template <class execution_space, class RV, class XV, class YV, class size_type>
void MV_Dot_Invoke(
    const execution_space& space, const RV& r, const XV& x, const YV& y,
    typename std::enable_if<
        !Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible>::type* = nullptr) {
  Kokkos::View<typename RV::non_const_value_type*, typename XV::memory_space> tempResult(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Dot_MV temp result"), r.extent(0));
  MV_Dot_Invoke<execution_space, decltype(tempResult), XV, YV, size_type>(space, tempResult, x, y);
  Kokkos::deep_copy(space, r, tempResult);
  space.fence();
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_
