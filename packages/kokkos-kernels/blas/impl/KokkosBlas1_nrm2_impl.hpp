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
#ifndef KOKKOSBLAS1_NRM2_IMPL_HPP_
#define KOKKOSBLAS1_NRM2_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_nrm2_spec.hpp>
#include <KokkosBlas_util.hpp>

namespace KokkosBlas {
namespace Impl {

//
// nrm2_squared
//

/// \brief 2-norm (squared) functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template <class RV, class XV, class SizeType = typename XV::size_type>
struct V_Nrm2_Functor {
  typedef SizeType size_type;
  typedef typename XV::non_const_value_type xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::ArithTraits<typename IPT::mag_type> AT;
  typedef typename IPT::mag_type value_type;

  typename XV::const_type m_x;
  bool m_take_sqrt;

  V_Nrm2_Functor(const XV& x, bool take_sqrt) : m_x(x), m_take_sqrt(take_sqrt) {
    static_assert(Kokkos::is_view<RV>::value,
                  "KokkosBlas::Impl::V_Nrm2_Functor: "
                  "R is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XV>::value,
                  "KokkosBlas::Impl::V_Nrm2_Functor: "
                  "X is not a Kokkos::View.");
    static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                  "KokkosBlas::Impl::V_Nrm2_Functor: R is const.  "
                  "It must be nonconst, because it is an output argument "
                  "(we have to be able to write to its entries).");
    static_assert(RV::rank == 0 && XV::rank == 1,
                  "KokkosBlas::Impl::V_Nrm2_Functor: "
                  "RV must have rank 0 and XV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i, value_type& sum) const {
    const typename IPT::mag_type tmp = IPT::norm(m_x(i));
    sum += tmp * tmp;
  }

  KOKKOS_INLINE_FUNCTION void init(value_type& update) const { update = AT::zero(); }

  KOKKOS_INLINE_FUNCTION void join(value_type& update, const value_type& source) const { update += source; }

  KOKKOS_INLINE_FUNCTION void final(value_type& update) const {
    if (m_take_sqrt) update = Kokkos::ArithTraits<typename RV::non_const_value_type>::sqrt(update);
  }
};

/// \brief Column-wise 2-norm functor for multivectors; works for
///   any layout, but best performance with LayoutLeft.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template <class ExecSpace, class RV, class XV, class size_type>
struct Nrm2_MV_Functor {
  typedef typename RV::non_const_value_type rvalue_type;
  typedef typename XV::non_const_value_type xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::ArithTraits<typename IPT::mag_type> AT;
  typedef typename IPT::mag_type value_type;

  using TeamMem = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  RV r;
  XV x;

  size_type teamsPerVec;  // number of teams collectively performing a dot product

  Nrm2_MV_Functor(const RV& r_, const XV& x_, int teamsPerVec_) : r(r_), x(x_), teamsPerVec(teamsPerVec_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamMem& t) const {
    size_type globalRank = t.league_rank();
    size_type localRank  = globalRank % teamsPerVec;
    size_type i          = globalRank / teamsPerVec;

    value_type localResult = AT::zero();
    size_type begin        = localRank * (x.extent(0) / teamsPerVec);
    size_type end          = (localRank + 1) * (x.extent(0) / teamsPerVec);
    if (localRank == teamsPerVec - 1) end = x.extent(0);

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(t, begin, end),
        [&](size_type k, value_type& update) {
          value_type tmp = IPT::norm(x(k, i));
          update += tmp * tmp;
        },
        localResult);

    Kokkos::single(Kokkos::PerTeam(t), [&]() { Kokkos::atomic_add(&r(i), rvalue_type(localResult)); });
  }
};

/// \brief Compute the 2-norm (or its square) of the single vector (1-D
///   View) X, and store the result in the 0-D View r.
template <class execution_space, class RV, class XV, class SizeType>
void V_Nrm2_Invoke(const execution_space& space, const RV& r, const XV& X, const bool& take_sqrt) {
  const SizeType numRows = static_cast<SizeType>(X.extent(0));
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  using functor_type = V_Nrm2_Functor<RV, XV, SizeType>;
  functor_type op(X, take_sqrt);
  Kokkos::parallel_reduce("KokkosBlas::Nrm2::S0", policy, op, r);
}

/// \brief Compute the 2-norms (or their square) of the columns of the
///   multivector (2-D View) X, and store result(s) in the 1-D View r.
// Main version: the result view is accessible from execution space, so it can
// be computed in-place
template <class execution_space, class RV, class XV, class size_type>
void MV_Nrm2_Invoke(
    const execution_space& space, const RV& r, const XV& x, bool take_sqrt,
    typename std::enable_if<Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible>::type* =
        nullptr) {
  if (r.extent(0) != x.extent(1)) {
    std::ostringstream oss;
    oss << "KokkosBlas::nrm2 (rank-2): result vector has wrong length (" << r.extent(0) << ", but x has " << x.extent(1)
        << " columns)";
    throw std::runtime_error(oss.str());
  }
  // Zero out the result vector
  Kokkos::deep_copy(space, r, Kokkos::ArithTraits<typename RV::non_const_value_type>::zero());
  size_type teamsPerVec;
  KokkosBlas::Impl::multipleReductionWorkDistribution<execution_space, size_type>(x.extent(0), x.extent(1),
                                                                                  teamsPerVec);
  size_type numTeams = x.extent(1) * teamsPerVec;
  Kokkos::TeamPolicy<execution_space> pol(space, numTeams, Kokkos::AUTO);
  Kokkos::parallel_for("KokkosBlas1::Nrm2::S1", pol,
                       Nrm2_MV_Functor<execution_space, RV, XV, size_type>(r, x, teamsPerVec));
  if (take_sqrt) {
    Kokkos::parallel_for("KokkosBlas1::Nrm2::Sqrt", Kokkos::RangePolicy<execution_space>(space, 0, r.extent(0)),
                         TakeSqrtFunctor<RV>(r));
  }
}

// Version for when a temporary result view is needed (implemented in terms of
// the other version)
template <class execution_space, class RV, class XV, class size_type>
void MV_Nrm2_Invoke(
    const execution_space& space, const RV& r, const XV& x, bool take_sqrt,
    typename std::enable_if<
        !Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible>::type* = nullptr) {
  Kokkos::View<typename RV::non_const_value_type*, typename XV::memory_space> tempResult(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Nrm2 temp result"), r.extent(0));
  MV_Nrm2_Invoke<execution_space, decltype(tempResult), XV, size_type>(space, tempResult, x, take_sqrt);
  Kokkos::deep_copy(space, r, tempResult);
  space.fence();
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_NRM2_IMPL_HPP_
