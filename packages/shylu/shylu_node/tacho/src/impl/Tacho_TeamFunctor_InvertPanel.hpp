// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_TEAMFUNCTOR_INVERT_PANEL_HPP__
#define __TACHO_TEAMFUNCTOR_INVERT_PANEL_HPP__

/// \file Tacho_TeamFunctor_InvertPanel.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_InvertPanel {
public:
  typedef Kokkos::pair<ordinal_type, ordinal_type> range_type;

  typedef SupernodeInfoType supernode_info_type;
  typedef typename supernode_info_type::supernode_type supernode_type;

  typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
  typedef typename supernode_info_type::size_type_array size_type_array;

  typedef typename supernode_info_type::value_type value_type;
  typedef typename supernode_info_type::value_type_array value_type_array;
  typedef typename supernode_info_type::value_type_matrix value_type_matrix;

private:
  supernode_info_type _info;
  ordinal_type_array _prepare_mode;
  ordinal_type _scratch_level;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_InvertPanel() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_InvertPanel(const supernode_info_type &info, const ordinal_type_array &prepare_mode,
                          const ordinal_type scratch_level)
      : _info(info), _prepare_mode(prepare_mode), _scratch_level(scratch_level) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void copyAndSetIdentity(MemberType &member, const value_type use_this_one, value_type_matrix A,
                                                 value_type_matrix P) const {
    const value_type zero(0);
    const ordinal_type m = A.extent(0);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m * m), [&](const ordinal_type &k) {
      const ordinal_type i = k % m;
      const ordinal_type j = k / m;
      A(i, j) = i <= j ? P(i, j) : zero;
      P(i, j) = i == j ? use_this_one : zero;
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void invert(MemberType &member, const value_type use_this_one, const value_type_matrix &A,
                                     const value_type_matrix &P) const {
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), use_this_one, A,
                                                                            P);
  }

  template <int Var> struct VariantTag {};

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const VariantTag<0> &, const MemberType &member) const {
    // dummy
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const VariantTag<1> &, const MemberType &member) const {
    const ordinal_type sid = member.league_rank();
    const ordinal_type mode = _prepare_mode(sid);
    if (mode == 1) {
      const auto s = _info.supernodes(sid);
      const ordinal_type m = s.m;
      UnmanagedViewType<value_type_matrix> A(member.team_scratch(_scratch_level), m, m);
      UnmanagedViewType<value_type_matrix> P(s.buf, m, m);
      const value_type one(1);
      copyAndSetIdentity(member, one, A, P);
      invert(member, one, A, P);
    } else {
      // skip
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const VariantTag<2> &, const MemberType &member) const {
    const ordinal_type sid = member.league_rank();
    const ordinal_type mode = _prepare_mode(sid);
    if (mode == 1) {
      const auto s = _info.supernodes(sid);
      const ordinal_type m = s.m, n = s.n;
      UnmanagedViewType<value_type_matrix> A(member.team_scratch(_scratch_level), m, m);
      UnmanagedViewType<value_type_matrix> P(s.buf, m, n);
      const value_type minus_one(-1);
      copyAndSetIdentity(member, minus_one, A, P);
      invert(member, minus_one, A, P);
    } else {
      // skip
    }
  }
};

} // namespace Tacho

#endif
