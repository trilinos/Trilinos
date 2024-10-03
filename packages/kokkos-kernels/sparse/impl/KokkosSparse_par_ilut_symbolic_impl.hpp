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

#ifndef KOKKOSSPARSE_IMPL_PAR_ILUT_SYMBOLIC_HPP_
#define KOKKOSSPARSE_IMPL_PAR_ILUT_SYMBOLIC_HPP_

/// \file KokkosSparse_par_ilut_symbolic_impl.hpp
/// \brief Implementation of the symbolic phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_par_ilut_handle.hpp>
#include <KokkosSparse_par_ilut_numeric_impl.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosKernels_Error.hpp>

// #define SYMBOLIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

template <class IlutHandle, class ARowMapType, class AEntriesType, class LRowMapType, class URowMapType>
void ilut_symbolic(IlutHandle& thandle, const ARowMapType& A_row_map_d, const AEntriesType& A_entries_d,
                   LRowMapType& L_row_map_d, URowMapType& U_row_map_d) {
  using execution_space = typename ARowMapType::execution_space;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;
  using size_type       = typename IlutHandle::size_type;
  using Ilut            = IlutWrap<IlutHandle>;

  const size_type a_nrows = A_row_map_d.extent(0);
  const size_type nrows   = a_nrows > 0 ? (a_nrows - 1) : 0;
  thandle.set_nrows(nrows);

  const auto policy = thandle.get_default_team_policy();

  // Sizing for the initial L/U approximation
  Kokkos::parallel_for(
      "symbolic sizing", policy, KOKKOS_LAMBDA(const member_type& team) {
        const auto row_idx = team.league_rank();

        const auto row_nnz_begin = A_row_map_d(row_idx);
        const auto row_nnz_end   = A_row_map_d(row_idx + 1);

        size_type nnzsL_temp = 0, nnzsU_temp = 0;
        // Multi-reductions are not supported at the TeamThread level
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, row_nnz_begin, row_nnz_end),
            [&](const size_type nnz, size_type& nnzsL_inner) {
              const auto col_idx = A_entries_d(nnz);
              nnzsL_inner += col_idx < row_idx;
            },
            nnzsL_temp);

        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, row_nnz_begin, row_nnz_end),
            [&](const size_type nnz, size_type& nnzsU_inner) {
              const auto col_idx = A_entries_d(nnz);
              nnzsU_inner += col_idx > row_idx;
            },
            nnzsU_temp);

        team.team_barrier();

        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          L_row_map_d(row_idx) = nnzsL_temp + 1;
          U_row_map_d(row_idx) = nnzsU_temp + 1;
        });
      });

  const size_type nnzsL = Ilut::prefix_sum(L_row_map_d);
  const size_type nnzsU = Ilut::prefix_sum(U_row_map_d);

  // Set symbolic info on handle
  thandle.set_nnzL(nnzsL);
  thandle.set_nnzU(nnzsU);
  thandle.set_symbolic_complete();

}  // end ilut_symbolic

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#endif
