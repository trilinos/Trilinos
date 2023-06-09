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

#ifndef KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_

/// \file KokkosSparse_spiluk_numeric_impl.hpp
/// \brief Implementation(s) of the numeric phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_spiluk_handle.hpp>

//#define NUMERIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

// struct UnsortedTag {};

template <class ARowMapType, class AEntriesType, class AValuesType,
          class LRowMapType, class LEntriesType, class LValuesType,
          class URowMapType, class UEntriesType, class UValuesType,
          class LevelViewType, class WorkViewType, class nnz_lno_t>
struct ILUKLvlSchedRPNumericFunctor {
  using lno_t    = typename AEntriesType::non_const_value_type;
  using scalar_t = typename AValuesType::non_const_value_type;
  ARowMapType A_row_map;
  AEntriesType A_entries;
  AValuesType A_values;
  LRowMapType L_row_map;
  LEntriesType L_entries;
  LValuesType L_values;
  URowMapType U_row_map;
  UEntriesType U_entries;
  UValuesType U_values;
  LevelViewType level_idx;
  WorkViewType iw;
  nnz_lno_t lev_start;

  ILUKLvlSchedRPNumericFunctor(
      const ARowMapType &A_row_map_, const AEntriesType &A_entries_,
      const AValuesType &A_values_, const LRowMapType &L_row_map_,
      const LEntriesType &L_entries_, LValuesType &L_values_,
      const URowMapType &U_row_map_, const UEntriesType &U_entries_,
      UValuesType &U_values_, const LevelViewType &level_idx_,
      WorkViewType &iw_, const nnz_lno_t &lev_start_)
      : A_row_map(A_row_map_),
        A_entries(A_entries_),
        A_values(A_values_),
        L_row_map(L_row_map_),
        L_entries(L_entries_),
        L_values(L_values_),
        U_row_map(U_row_map_),
        U_entries(U_entries_),
        U_values(U_values_),
        level_idx(level_idx_),
        iw(iw_),
        lev_start(lev_start_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const lno_t i) const {
    auto rowid = level_idx(i);
    auto tid   = i - lev_start;
    auto k1    = L_row_map(rowid);
    auto k2    = L_row_map(rowid + 1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2 - 1; ++k) {
#else
    for (auto k = k1; k < k2; ++k) {
#endif
      auto col     = L_entries(k);
      L_values(k)  = 0.0;
      iw(tid, col) = k;
    }
#ifdef KEEP_DIAG
    L_values(k2 - 1) = scalar_t(1.0);
#endif

    k1 = U_row_map(rowid);
    k2 = U_row_map(rowid + 1);
    for (auto k = k1; k < k2; ++k) {
      auto col     = U_entries(k);
      U_values(k)  = 0.0;
      iw(tid, col) = k;
    }

    // Unpack the ith row of A
    k1 = A_row_map(rowid);
    k2 = A_row_map(rowid + 1);
    for (auto k = k1; k < k2; ++k) {
      auto col  = A_entries(k);
      auto ipos = iw(tid, col);
      if (col < rowid)
        L_values(ipos) = A_values(k);
      else
        U_values(ipos) = A_values(k);
    }

    // Eliminate prev rows
    k1 = L_row_map(rowid);
    k2 = L_row_map(rowid + 1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2 - 1; ++k) {
#else
    for (auto k = k1; k < k2; ++k) {
#endif
      auto prev_row = L_entries(k);
#ifdef KEEP_DIAG
      auto fact = L_values(k) / U_values(U_row_map(prev_row));
#else
      auto fact = L_values(k) * U_values(U_row_map(prev_row));
#endif
      L_values(k) = fact;
      for (auto kk = U_row_map(prev_row) + 1; kk < U_row_map(prev_row + 1);
           ++kk) {
        auto col  = U_entries(kk);
        auto ipos = iw(tid, col);
        if (ipos == -1) continue;
        auto lxu = -U_values(kk) * fact;
        if (col < rowid)
          L_values(ipos) += lxu;
        else
          U_values(ipos) += lxu;
      }  // end for kk
    }    // end for k

#ifdef KEEP_DIAG
    if (U_values(iw(tid, rowid)) == 0.0) {
      U_values(iw(tid, rowid)) = 1e6;
    }
#else
    if (U_values(iw(tid, rowid)) == 0.0) {
      U_values(iw(tid, rowid)) = 1e6;
    } else {
      U_values(iw(tid, rowid)) = 1.0 / U_values(iw(tid, rowid));
    }
#endif

    // Reset
    k1 = L_row_map(rowid);
    k2 = L_row_map(rowid + 1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2 - 1; ++k)
#else
    for (auto k = k1; k < k2; ++k)
#endif
      iw(tid, L_entries(k)) = -1;

    k1 = U_row_map(rowid);
    k2 = U_row_map(rowid + 1);
    for (auto k = k1; k < k2; ++k) iw(tid, U_entries(k)) = -1;
  }
};

template <class ARowMapType, class AEntriesType, class AValuesType,
          class LRowMapType, class LEntriesType, class LValuesType,
          class URowMapType, class UEntriesType, class UValuesType,
          class LevelViewType, class WorkViewType, class nnz_lno_t>
struct ILUKLvlSchedTP1NumericFunctor {
  using execution_space = typename ARowMapType::execution_space;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;
  using size_type       = typename ARowMapType::non_const_value_type;
  using lno_t           = typename AEntriesType::non_const_value_type;
  using scalar_t        = typename AValuesType::non_const_value_type;

  ARowMapType A_row_map;
  AEntriesType A_entries;
  AValuesType A_values;
  LRowMapType L_row_map;
  LEntriesType L_entries;
  LValuesType L_values;
  URowMapType U_row_map;
  UEntriesType U_entries;
  UValuesType U_values;
  LevelViewType level_idx;
  WorkViewType iw;
  nnz_lno_t lev_start;

  ILUKLvlSchedTP1NumericFunctor(
      const ARowMapType &A_row_map_, const AEntriesType &A_entries_,
      const AValuesType &A_values_, const LRowMapType &L_row_map_,
      const LEntriesType &L_entries_, LValuesType &L_values_,
      const URowMapType &U_row_map_, const UEntriesType &U_entries_,
      UValuesType &U_values_, const LevelViewType &level_idx_,
      WorkViewType &iw_, const nnz_lno_t &lev_start_)
      : A_row_map(A_row_map_),
        A_entries(A_entries_),
        A_values(A_values_),
        L_row_map(L_row_map_),
        L_entries(L_entries_),
        L_values(L_values_),
        U_row_map(U_row_map_),
        U_entries(U_entries_),
        U_values(U_values_),
        level_idx(level_idx_),
        iw(iw_),
        lev_start(lev_start_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &team) const {
    nnz_lno_t my_team = static_cast<nnz_lno_t>(team.league_rank());
    nnz_lno_t rowid =
        static_cast<nnz_lno_t>(level_idx(my_team + lev_start));  // map to rowid

    size_type k1 = static_cast<size_type>(L_row_map(rowid));
    size_type k2 = static_cast<size_type>(L_row_map(rowid + 1));
#ifdef KEEP_DIAG
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2 - 1),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(L_entries(k));
                           L_values(k)   = 0.0;
                           iw(my_team, col) = static_cast<nnz_lno_t>(k);
                         });
#else
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(L_entries(k));
                           L_values(k)   = 0.0;
                           iw(my_team, col) = static_cast<nnz_lno_t>(k);
                         });
#endif

#ifdef KEEP_DIAG
    // if (my_thread == 0) L_values(k2 - 1) = scalar_t(1.0);
    Kokkos::single(Kokkos::PerTeam(team),
                   [&]() { L_values(k2 - 1) = scalar_t(1.0); });
#endif

    team.team_barrier();

    k1 = static_cast<size_type>(U_row_map(rowid));
    k2 = static_cast<size_type>(U_row_map(rowid + 1));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(U_entries(k));
                           U_values(k)   = 0.0;
                           iw(my_team, col) = static_cast<nnz_lno_t>(k);
                         });

    team.team_barrier();

    // Unpack the ith row of A
    k1 = static_cast<size_type>(A_row_map(rowid));
    k2 = static_cast<size_type>(A_row_map(rowid + 1));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(A_entries(k));
                           nnz_lno_t ipos = iw(my_team, col);
                           if (col < rowid)
                             L_values(ipos) = A_values(k);
                           else
                             U_values(ipos) = A_values(k);
                         });

    team.team_barrier();

    // Eliminate prev rows
    k1 = static_cast<size_type>(L_row_map(rowid));
    k2 = static_cast<size_type>(L_row_map(rowid + 1));
#ifdef KEEP_DIAG
    for (size_type k = k1; k < k2 - 1; k++)
#else
    for (size_type k = k1; k < k2; k++)
#endif
    {
      nnz_lno_t prev_row = L_entries(k);

      scalar_t fact = scalar_t(0.0);
      Kokkos::single(
          Kokkos::PerTeam(team),
          [&](scalar_t &tmp_fact) {
#ifdef KEEP_DIAG
            tmp_fact = L_values(k) / U_values(U_row_map(prev_row));
#else
            tmp_fact = L_values(k) * U_values(U_row_map(prev_row));
#endif
            L_values(k) = tmp_fact;
          },
          fact);

      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, U_row_map(prev_row) + 1,
                                  U_row_map(prev_row + 1)),
          [&](const size_type kk) {
            nnz_lno_t col  = static_cast<nnz_lno_t>(U_entries(kk));
            nnz_lno_t ipos = iw(my_team, col);
            auto lxu       = -U_values(kk) * fact;
            if (ipos != -1) {
              if (col < rowid)
                L_values(ipos) += lxu;
              else
                U_values(ipos) += lxu;
            }
          });  // end for kk

      team.team_barrier();
    }  // end for k

    // if (my_thread == 0) {
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      nnz_lno_t ipos = iw(my_team, rowid);
#ifdef KEEP_DIAG
      if (U_values(ipos) == 0.0) {
        U_values(ipos) = 1e6;
      }
#else
      if (U_values(ipos) == 0.0) {
        U_values(ipos) = 1e6;
      } else {
        U_values(ipos) = 1.0 / U_values(ipos);
      }
#endif
    });
    //}

    team.team_barrier();

    // Reset
    k1 = static_cast<size_type>(L_row_map(rowid));
    k2 = static_cast<size_type>(L_row_map(rowid + 1));
#ifdef KEEP_DIAG
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2 - 1),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(L_entries(k));
                           iw(my_team, col) = -1;
                         });
#else
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(L_entries(k));
                           iw(my_team, col) = -1;
                         });
#endif

    k1 = static_cast<size_type>(U_row_map(rowid));
    k2 = static_cast<size_type>(U_row_map(rowid + 1));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2),
                         [&](const size_type k) {
                           nnz_lno_t col = static_cast<nnz_lno_t>(U_entries(k));
                           iw(my_team, col) = -1;
                         });
  }
};

template <class IlukHandle, class ARowMapType, class AEntriesType,
          class AValuesType, class LRowMapType, class LEntriesType,
          class LValuesType, class URowMapType, class UEntriesType,
          class UValuesType>
void iluk_numeric(IlukHandle &thandle, const ARowMapType &A_row_map,
                  const AEntriesType &A_entries, const AValuesType &A_values,
                  const LRowMapType &L_row_map, const LEntriesType &L_entries,
                  LValuesType &L_values, const URowMapType &U_row_map,
                  const UEntriesType &U_entries, UValuesType &U_values) {
  using execution_space         = typename IlukHandle::execution_space;
  using size_type               = typename IlukHandle::size_type;
  using nnz_lno_t               = typename IlukHandle::nnz_lno_t;
  using HandleDeviceEntriesType = typename IlukHandle::nnz_lno_view_t;
  using WorkViewType            = typename IlukHandle::work_view_t;
  using LevelHostViewType       = typename IlukHandle::nnz_lno_view_host_t;

  size_type nlevels = thandle.get_num_levels();
  int team_size     = thandle.get_team_size();

  LevelHostViewType level_ptr_h     = thandle.get_host_level_ptr();
  HandleDeviceEntriesType level_idx = thandle.get_level_idx();

  LevelHostViewType level_nchunks_h, level_nrowsperchunk_h;
  WorkViewType iw;

  //{
  if (thandle.get_algorithm() ==
      KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1) {
    level_nchunks_h       = thandle.get_level_nchunks();
    level_nrowsperchunk_h = thandle.get_level_nrowsperchunk();
  }
  iw = thandle.get_iw();

  // Main loop must be performed sequential. Question: Try out Cuda's graph
  // stuff to reduce kernel launch overhead
  for (size_type lvl = 0; lvl < nlevels; ++lvl) {
    nnz_lno_t lev_start = level_ptr_h(lvl);
    nnz_lno_t lev_end   = level_ptr_h(lvl + 1);

    if ((lev_end - lev_start) != 0) {
      if (thandle.get_algorithm() ==
          KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP) {
        Kokkos::parallel_for(
            "parfor_fixed_lvl",
            Kokkos::RangePolicy<execution_space>(lev_start, lev_end),
            ILUKLvlSchedRPNumericFunctor<
                ARowMapType, AEntriesType, AValuesType, LRowMapType,
                LEntriesType, LValuesType, URowMapType, UEntriesType,
                UValuesType, HandleDeviceEntriesType, WorkViewType, nnz_lno_t>(
                A_row_map, A_entries, A_values, L_row_map, L_entries, L_values,
                U_row_map, U_entries, U_values, level_idx, iw, lev_start));
      } else if (thandle.get_algorithm() ==
                 KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1) {
        using policy_type = Kokkos::TeamPolicy<execution_space>;

        nnz_lno_t lvl_rowid_start = 0;
        nnz_lno_t lvl_nrows_chunk;
        for (int chunkid = 0; chunkid < level_nchunks_h(lvl); chunkid++) {
          if ((lvl_rowid_start + level_nrowsperchunk_h(lvl)) >
              (lev_end - lev_start))
            lvl_nrows_chunk = (lev_end - lev_start) - lvl_rowid_start;
          else
            lvl_nrows_chunk = level_nrowsperchunk_h(lvl);

          ILUKLvlSchedTP1NumericFunctor<
              ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType,
              LValuesType, URowMapType, UEntriesType, UValuesType,
              HandleDeviceEntriesType, WorkViewType, nnz_lno_t>
              tstf(A_row_map, A_entries, A_values, L_row_map, L_entries,
                   L_values, U_row_map, U_entries, U_values, level_idx, iw,
                   lev_start + lvl_rowid_start);

          if (team_size == -1)
            Kokkos::parallel_for(
                "parfor_tp1", policy_type(lvl_nrows_chunk, Kokkos::AUTO), tstf);
          else
            Kokkos::parallel_for("parfor_tp1",
                                 policy_type(lvl_nrows_chunk, team_size), tstf);
          Kokkos::fence();
          lvl_rowid_start += lvl_nrows_chunk;
        }
      }
    }  // end if
  }    // end for lvl
  //}

// Output check
#ifdef NUMERIC_OUTPUT_INFO
  std::cout << "  iluk_numeric result: " << std::endl;

  std::cout << "  nnzL: " << thandle.get_nnzL() << std::endl;
  std::cout << "  L_row_map = ";
  for (size_type i = 0; i < thandle.get_nrows() + 1; ++i) {
    std::cout << L_row_map(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "  L_entries = ";
  for (size_type i = 0; i < thandle.get_nnzL(); ++i) {
    std::cout << L_entries(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "  L_values = ";
  for (size_type i = 0; i < thandle.get_nnzL(); ++i) {
    std::cout << L_values(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "  nnzU: " << thandle.get_nnzU() << std::endl;
  std::cout << "  U_row_map = ";
  for (size_type i = 0; i < thandle.get_nrows() + 1; ++i) {
    std::cout << U_row_map(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "  U_entries = ";
  for (size_type i = 0; i < thandle.get_nnzU(); ++i) {
    std::cout << U_entries(i) << " ";
  }
  std::cout << std::endl;

  std::cout << "  U_values = ";
  for (size_type i = 0; i < thandle.get_nnzU(); ++i) {
    std::cout << U_values(i) << " ";
  }
  std::cout << std::endl;
#endif

}  // end iluk_numeric

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#endif
