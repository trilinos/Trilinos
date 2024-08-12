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

#ifndef KOKKOSSPARSE_MDF_IMPL_HPP_
#define KOKKOSSPARSE_MDF_IMPL_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include "KokkosKernels_Sorting.hpp"
#include "KokkosSparse_findRelOffset.hpp"
#include <type_traits>
#include "Kokkos_ArithTraits.hpp"

namespace KokkosSparse {
namespace Impl {

template <typename crs_matrix_type>
struct MDF_types {
  using scalar_type     = typename crs_matrix_type::value_type;
  using KAS             = typename Kokkos::ArithTraits<scalar_type>;
  using scalar_mag_type = typename KAS::mag_type;
  using values_mag_type = Kokkos::View<scalar_mag_type*, Kokkos::LayoutRight, typename crs_matrix_type::device_type,
                                       typename crs_matrix_type::memory_traits>;
};

template <class crs_matrix_type>
struct MDF_count_lower {
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using size_type    = typename crs_matrix_type::ordinal_type;
  using value_type   = typename crs_matrix_type::size_type;
  using KAV          = typename Kokkos::ArithTraits<value_type>;

  crs_matrix_type A;
  col_ind_type permutation;
  col_ind_type permutation_inv;

  using execution_space = typename crs_matrix_type::execution_space;
  using team_policy_t   = Kokkos::TeamPolicy<execution_space>;
  using team_member_t   = typename team_policy_t::member_type;

  MDF_count_lower(crs_matrix_type A_, col_ind_type permutation_, col_ind_type permutation_inv_)
      : A(A_), permutation(permutation_), permutation_inv(permutation_inv_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t team, value_type& update) const {
    const auto rowIdx  = team.league_rank();
    const auto rowView = A.graph.rowConst(rowIdx);

    value_type local_contrib = KAV::zero();
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, rowView.length),
        [&](const size_type entryIdx, value_type& partial) {
          if (rowView(entryIdx) <= rowIdx) partial += 1;
        },
        Kokkos::Sum<value_type, execution_space>(local_contrib));

    Kokkos::single(Kokkos::PerTeam(team), [&] {
      permutation(rowIdx)     = rowIdx;
      permutation_inv(rowIdx) = rowIdx;
      update += local_contrib;
    });
  }
};  // MDF_count_lower

template <class crs_matrix_type, bool is_initial_fill>
struct MDF_discarded_fill_norm {
  using device_type = typename crs_matrix_type::device_type;

  using static_crs_graph_type = typename crs_matrix_type::StaticCrsGraphType;
  using col_ind_type          = typename static_crs_graph_type::entries_type::non_const_type;
  using values_type           = typename crs_matrix_type::values_type::non_const_type;
  using values_mag_type       = typename MDF_types<crs_matrix_type>::values_mag_type;
  using size_type             = typename crs_matrix_type::size_type;
  using ordinal_type          = typename crs_matrix_type::ordinal_type;
  using scalar_type           = typename crs_matrix_type::value_type;
  using KAS                   = typename Kokkos::ArithTraits<scalar_type>;
  using scalar_mag_type       = typename KAS::mag_type;
  using KAM                   = typename Kokkos::ArithTraits<scalar_mag_type>;
  using permutation_set_type  = Kokkos::UnorderedMap<ordinal_type, void, device_type>;

  crs_matrix_type A, At;
  ordinal_type factorization_step;
  col_ind_type permutation;
  permutation_set_type permutation_set;
  col_ind_type update_list;

  values_mag_type discarded_fill;
  col_ind_type deficiency;
  int verbosity;

  MDF_discarded_fill_norm(crs_matrix_type A_, crs_matrix_type At_, ordinal_type factorization_step_,
                          col_ind_type permutation_, permutation_set_type permutation_set_,
                          values_mag_type discarded_fill_, col_ind_type deficiency_, int verbosity_,
                          col_ind_type update_list_ = col_ind_type{})
      : A(A_),
        At(At_),
        factorization_step(factorization_step_),
        permutation(permutation_),
        permutation_set(permutation_set_),
        update_list(update_list_),
        discarded_fill(discarded_fill_),
        deficiency(deficiency_),
        verbosity(verbosity_){};

  using execution_space = typename crs_matrix_type::execution_space;
  using team_policy_t   = Kokkos::TeamPolicy<execution_space>;
  using team_member_t   = typename team_policy_t::member_type;

  struct DiscNormReducer {
    using reducer = DiscNormReducer;
    struct value_type {
      scalar_mag_type discarded_norm;
      ordinal_type numFillEntries;
      scalar_type diag_val;
    };
    using result_view_type = Kokkos::View<value_type, execution_space>;

   private:
    result_view_type value;

   public:
    KOKKOS_INLINE_FUNCTION
    DiscNormReducer(value_type& value_) : value(&value_) {}

    KOKKOS_INLINE_FUNCTION
    static void join(value_type& dest, const value_type& src) {
      dest.discarded_norm += src.discarded_norm;
      dest.numFillEntries += src.numFillEntries;
      if (dest.diag_val == KAS::zero()) dest.diag_val = src.diag_val;
    }

    KOKKOS_INLINE_FUNCTION
    static void init(value_type& val) {
      val.discarded_norm = Kokkos::reduction_identity<scalar_mag_type>::sum();
      val.numFillEntries = Kokkos::reduction_identity<ordinal_type>::sum();
      val.diag_val       = KAS::zero();
    }

    KOKKOS_INLINE_FUNCTION
    static value_type init() {
      value_type out;
      init(out);
      return out;
    }

    KOKKOS_INLINE_FUNCTION
    value_type& reference() const { return *value.data(); }

    KOKKOS_INLINE_FUNCTION
    result_view_type view() const { return value; }
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(team_member_t team) const {
    const ordinal_type rowIdx =
        is_initial_fill ? permutation(team.league_rank()) : permutation(update_list(team.league_rank()));
    const auto colView = At.rowConst(rowIdx);
    const auto rowView = A.rowConst(rowIdx);

    using reduction_val_t         = typename DiscNormReducer::value_type;
    reduction_val_t reduction_val = DiscNormReducer::init();
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, colView.length),
        [&](const size_type alpha, reduction_val_t& running_disc_norm) {
          const ordinal_type fillRowIdx = colView.colidx(alpha);

          // Record diagonal term
          if (fillRowIdx == rowIdx) {
            Kokkos::single(Kokkos::PerThread(team), [&] { running_disc_norm.diag_val = colView.value(alpha); });
            return;
          }

          // Check if row already eliminated
          if constexpr (!is_initial_fill) {
            if (permutation_set.exists(fillRowIdx)) return;
          }

          const auto fillRowView              = A.rowConst(fillRowIdx);
          reduction_val_t local_reduction_val = DiscNormReducer::init();
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(team, rowView.length),
              [&](const ordinal_type beta, reduction_val_t& vect_running_disc_norm) {
                const ordinal_type fillColIdx = rowView.colidx(beta);

                if (fillColIdx == rowIdx) return;

                if constexpr (!is_initial_fill) {
                  if (permutation_set.exists(fillColIdx)) return;
                }

                bool entryIsDiscarded = true;
                for (ordinal_type gamma = 0; gamma < fillRowView.length; ++gamma) {
                  if (fillRowView.colidx(gamma) == fillColIdx) {
                    entryIsDiscarded = false;
                  }
                }
                if (entryIsDiscarded) {
                  vect_running_disc_norm.numFillEntries += 1;
                  vect_running_disc_norm.discarded_norm += KAS::abs(colView.value(alpha) * rowView.value(beta)) *
                                                           KAS::abs(colView.value(alpha) * rowView.value(beta));
                }
              },
              DiscNormReducer(local_reduction_val));

          Kokkos::single(Kokkos::PerThread(team), [&] {
            running_disc_norm.discarded_norm += local_reduction_val.discarded_norm;
            running_disc_norm.numFillEntries += local_reduction_val.numFillEntries;
          });
        },
        DiscNormReducer(reduction_val));

    Kokkos::single(Kokkos::PerTeam(team), [&] {
      const scalar_mag_type& discard_norm = reduction_val.discarded_norm;
      const ordinal_type& numFillEntries  = reduction_val.numFillEntries;
      const scalar_type& diag_val         = reduction_val.diag_val;

      // TODO add a check on `diag_val == zero`
      discarded_fill(rowIdx) = discard_norm / KAS::abs(diag_val * diag_val);
      deficiency(rowIdx)     = numFillEntries;
    });
  }
};  // MDF_discarded_fill_norm

template <class crs_matrix_type>
struct MDF_select_row {
  using values_type     = typename crs_matrix_type::values_type::non_const_type;
  using col_ind_type    = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using row_map_type    = typename crs_matrix_type::StaticCrsGraphType::row_map_type;
  using size_type       = typename crs_matrix_type::size_type;
  using ordinal_type    = typename crs_matrix_type::ordinal_type;
  using scalar_type     = typename crs_matrix_type::value_type;
  using values_mag_type = typename MDF_types<crs_matrix_type>::values_mag_type;

  // type used to perform the reduction
  // do not confuse it with scalar_type!
  using value_type = typename crs_matrix_type::ordinal_type;

  value_type factorization_step;
  values_mag_type discarded_fill;
  col_ind_type deficiency;
  row_map_type row_map;
  col_ind_type permutation;

  MDF_select_row(value_type factorization_step_, values_mag_type discarded_fill_, col_ind_type deficiency_,
                 row_map_type row_map_, col_ind_type permutation_)
      : factorization_step(factorization_step_),
        discarded_fill(discarded_fill_),
        deficiency(deficiency_),
        row_map(row_map_),
        permutation(permutation_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type src, ordinal_type& dst) const {
    const ordinal_type src_perm   = permutation(src);
    const ordinal_type dst_perm   = permutation(dst);
    const ordinal_type degree_src = row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst = row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if (discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src == degree_dst) && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const {
    const ordinal_type src_perm   = permutation(src);
    const ordinal_type dst_perm   = permutation(dst);
    const ordinal_type degree_src = row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst = row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if (discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) && (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src == degree_dst) && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& dst) const { dst = factorization_step; }

};  // MDF_select_row

template <class view_type, class ordinal_type>
KOKKOS_INLINE_FUNCTION bool sorted_view_contains(const view_type& values, const ordinal_type size,
                                                 typename view_type::const_value_type search_val) {
  return KokkosSparse::findRelOffset(values, size, search_val, size, true) != size;
}

template <class crs_matrix_type>
struct MDF_factorize_row {
  using device_type     = typename crs_matrix_type::device_type;
  using execution_space = typename crs_matrix_type::execution_space;
  using team_policy_t   = Kokkos::TeamPolicy<execution_space>;
  using team_member_t   = typename team_policy_t::member_type;

  using row_map_type         = typename crs_matrix_type::StaticCrsGraphType::row_map_type::non_const_type;
  using col_ind_type         = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type          = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type         = typename crs_matrix_type::ordinal_type;
  using size_type            = typename crs_matrix_type::size_type;
  using value_type           = typename crs_matrix_type::value_type;
  using values_mag_type      = typename MDF_types<crs_matrix_type>::values_mag_type;
  using value_mag_type       = typename values_mag_type::value_type;
  using permutation_set_type = Kokkos::UnorderedMap<ordinal_type, void, device_type>;

  crs_matrix_type A, At;

  row_map_type row_mapL;
  col_ind_type entriesL;
  values_type valuesL;

  row_map_type row_mapU;
  col_ind_type entriesU;
  values_type valuesU;

  col_ind_type permutation, permutation_inv;
  permutation_set_type permutation_set;
  values_mag_type discarded_fill;
  col_ind_type factored;
  ordinal_type selected_row_idx, factorization_step;

  col_ind_type update_list;

  int verbosity;

  MDF_factorize_row(crs_matrix_type A_, crs_matrix_type At_, row_map_type row_mapL_, col_ind_type entriesL_,
                    values_type valuesL_, row_map_type row_mapU_, col_ind_type entriesU_, values_type valuesU_,
                    col_ind_type permutation_, col_ind_type permutation_inv_, permutation_set_type permutation_set_,
                    values_mag_type discarded_fill_, col_ind_type factored_, ordinal_type selected_row_idx_,
                    ordinal_type factorization_step_, col_ind_type& update_list_, int verbosity_)
      : A(A_),
        At(At_),
        row_mapL(row_mapL_),
        entriesL(entriesL_),
        valuesL(valuesL_),
        row_mapU(row_mapU_),
        entriesU(entriesU_),
        valuesU(valuesU_),
        permutation(permutation_),
        permutation_inv(permutation_inv_),
        permutation_set(permutation_set_),
        discarded_fill(discarded_fill_),
        factored(factored_),
        selected_row_idx(selected_row_idx_),
        factorization_step(factorization_step_),
        update_list(update_list_),
        verbosity(verbosity_){};

  // Phase 2, do facrotization
  KOKKOS_INLINE_FUNCTION
  void operator()(team_member_t team) const {
    const auto alpha                = team.league_rank();
    const ordinal_type selected_row = permutation(factorization_step);
    const auto colView              = At.rowConst(selected_row);

    const auto rowInd = colView.colidx(alpha);
    if (rowInd == selected_row) return;

    if (permutation_set.exists(rowInd)) return;

    // Only one of the values will match selected so can just sum all contribs
    const auto rowView = A.rowConst(selected_row);
    value_type diag    = Kokkos::ArithTraits<value_type>::zero();
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, rowView.length),
        [&](const size_type ind, value_type& running_diag) {
          if (rowView.colidx(ind) == selected_row) running_diag = rowView.value(ind);
        },
        Kokkos::Sum<value_type, execution_space>(diag));

    // Extract alpha and beta vectors
    // Then insert alpha*beta/diag_val if the corresponding
    // entry in A is non-zero.
    auto fillRowView = A.row(rowInd);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowView.length), [&](const ordinal_type beta) {
      const auto colInd = rowView.colidx(beta);

      if (colInd == selected_row) return;

      if (permutation_set.exists(colInd)) return;

      const auto subVal = colView.value(alpha) * rowView.value(beta) / diag;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, fillRowView.length), [&](const ordinal_type gamma) {
        if (colInd == fillRowView.colidx(gamma)) {
          Kokkos::atomic_sub(&fillRowView.value(gamma), subVal);
        }
      });

      auto fillColView = At.row(colInd);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, fillColView.length), [&](const ordinal_type delt) {
        if (rowInd == fillColView.colidx(delt)) {
          Kokkos::atomic_sub(&fillColView.value(delt), subVal);
        }
      });
    });
  }
};

template <class crs_matrix_type>
struct MDF_compute_list_length {
  using device_type     = typename crs_matrix_type::device_type;
  using execution_space = typename crs_matrix_type::execution_space;
  using team_policy_t   = Kokkos::TeamPolicy<execution_space>;
  using team_member_t   = typename team_policy_t::member_type;

  using row_map_type    = typename crs_matrix_type::StaticCrsGraphType::row_map_type::non_const_type;
  using col_ind_type    = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type     = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type    = typename crs_matrix_type::ordinal_type;
  using size_type       = typename crs_matrix_type::size_type;
  using value_type      = typename crs_matrix_type::value_type;
  using values_mag_type = typename MDF_types<crs_matrix_type>::values_mag_type;
  using value_mag_type  = typename values_mag_type::value_type;

  using permutation_set_type = Kokkos::UnorderedMap<ordinal_type, void, device_type>;

  crs_matrix_type A, At;

  row_map_type row_mapL;
  col_ind_type entriesL;
  values_type valuesL;

  row_map_type row_mapU;
  col_ind_type entriesU;
  values_type valuesU;

  col_ind_type permutation, permutation_inv;
  permutation_set_type permutation_set;
  values_mag_type discarded_fill;
  col_ind_type factored;
  ordinal_type selected_row_idx, factorization_step;

  col_ind_type update_list;

  int verbosity;

  MDF_compute_list_length(crs_matrix_type A_, crs_matrix_type At_, row_map_type row_mapL_, col_ind_type entriesL_,
                          values_type valuesL_, row_map_type row_mapU_, col_ind_type entriesU_, values_type valuesU_,
                          col_ind_type permutation_, col_ind_type permutation_inv_,
                          permutation_set_type permutation_set_, values_mag_type discarded_fill_,
                          col_ind_type factored_, ordinal_type selected_row_idx_, ordinal_type factorization_step_,
                          col_ind_type& update_list_, int verbosity_)
      : A(A_),
        At(At_),
        row_mapL(row_mapL_),
        entriesL(entriesL_),
        valuesL(valuesL_),
        row_mapU(row_mapU_),
        entriesU(entriesU_),
        valuesU(valuesU_),
        permutation(permutation_),
        permutation_inv(permutation_inv_),
        permutation_set(permutation_set_),
        discarded_fill(discarded_fill_),
        factored(factored_),
        selected_row_idx(selected_row_idx_),
        factorization_step(factorization_step_),
        update_list(update_list_),
        verbosity(verbosity_){};

  // Phase 1, update list length
  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t team, ordinal_type& update_list_len, ordinal_type& selected_row_len) const {
    ordinal_type selected_row = 0;

    size_type U_entryIdx = row_mapU(factorization_step);
    size_type L_entryIdx = row_mapL(factorization_step);

    Kokkos::single(Kokkos::PerTeam(team), [&] {
      selected_row                 = permutation(selected_row_idx);
      discarded_fill(selected_row) = Kokkos::ArithTraits<value_mag_type>::max();

      // Swap entries in permutation vectors
      permutation(selected_row_idx)                    = permutation(factorization_step);
      permutation(factorization_step)                  = selected_row;
      permutation_inv(permutation(factorization_step)) = factorization_step;
      permutation_inv(permutation(selected_row_idx))   = selected_row_idx;

      // Diagonal value of L
      entriesL(L_entryIdx) = selected_row;
      valuesL(L_entryIdx)  = Kokkos::ArithTraits<value_type>::one();

      // Insert into permutation set for later
      const auto res = permutation_set.insert(selected_row);
      (void)res;  // avoid unused error
      assert(res.success());
    });
    ++L_entryIdx;

    // Only one thread has the selected row
    team.team_reduce(Kokkos::Max<ordinal_type, execution_space>(selected_row));
    const auto rowView = A.rowConst(selected_row);
    const auto colView = At.rowConst(selected_row);

    // Insert the upper part of the selected row in U
    // including the diagonal term.
    ordinal_type updateIdx = 0;
    value_type diag        = Kokkos::ArithTraits<value_type>::zero();
    {
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, rowView.length),
                            [&](const size_type alpha, ordinal_type& running_update, bool is_final) {
                              const auto colInd = rowView.colidx(alpha);
                              if ((colInd != selected_row) && (factored(colInd) != 1)) {
                                if (is_final) {
                                  update_list(running_update) = colInd;
                                  ++updateIdx;
                                }
                                ++running_update;
                              }
                            }
                            // ,updateIdx
      );

      // Until https://github.com/kokkos/kokkos/issues/6259 is resolved, do
      // reduction outside of parallel_scan
      team.team_reduce(Kokkos::Sum<ordinal_type, execution_space>(updateIdx));

      // Sort update list
      KokkosKernels::TeamBitonicSort(&update_list(0), updateIdx, team);
    }
    {
      size_type numEntrU = 0;
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, rowView.length),
                            [&](const size_type alpha, size_type& running_nEntr, bool is_final) {
                              const auto colInd = rowView.colidx(alpha);
                              if (permutation_inv(colInd) >= factorization_step) {
                                if (is_final) {
                                  ++numEntrU;
                                  entriesU(U_entryIdx + running_nEntr) = colInd;
                                  valuesU(U_entryIdx + running_nEntr)  = rowView.value(alpha);
                                  if (colInd == selected_row) diag = rowView.value(alpha);
                                }
                                ++running_nEntr;
                              }
                            }
                            // , numEntrU
      );

      // Until https://github.com/kokkos/kokkos/issues/6259 is resolved, do
      // reduction outside of parallel_scan
      team.team_reduce(Kokkos::Sum<size_type, execution_space>(numEntrU));

      U_entryIdx += numEntrU;
    }

    // Only one thread found diagonal so just sum over all
    team.team_reduce(Kokkos::Sum<value_type, execution_space>(diag));

    // Insert the lower part of the selected column of A
    // divided by its the diagonal value to obtain a unit
    // diagonal value in L.
    {
      size_type numEntrL = 0;
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, colView.length),
                            [&](const size_type alpha, size_type& running_nEntr, bool is_final) {
                              const auto rowInd = colView.colidx(alpha);
                              if (permutation_inv(rowInd) > factorization_step) {
                                if (is_final) {
                                  ++numEntrL;
                                  entriesL(L_entryIdx + running_nEntr) = rowInd;
                                  valuesL(L_entryIdx + running_nEntr)  = colView.value(alpha) / diag;
                                }
                                ++running_nEntr;
                              }
                            }
                            // , numEntrL
      );

      // Until https://github.com/kokkos/kokkos/issues/6259 is resolved, do
      // reduction outside of parallel_scan
      team.team_reduce(Kokkos::Sum<size_type, execution_space>(numEntrL));

      L_entryIdx += numEntrL;
    }
    {
      ordinal_type numUpdateL = 0;
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, colView.length),
                            [&](const size_type alpha, ordinal_type& running_update, bool is_final) {
                              const auto rowInd = colView.colidx(alpha);
                              if ((rowInd != selected_row) && (factored(rowInd) != 1)) {
                                // updateIdx currently holds the rows that were updated. don't add
                                // duplicates
                                const size_type& update_rows = updateIdx;

                                const bool already_updated = sorted_view_contains(update_list, update_rows, rowInd);

                                if (!already_updated) {
                                  // Cannot make use of vector ranges until
                                  // https://github.com/kokkos/kokkos/issues/6259 is resolved
                                  // Kokkos::single(Kokkos::PerThread(team),[&]{
                                  if (is_final) {
                                    update_list(updateIdx + running_update) = rowInd;
                                    ++numUpdateL;
                                  }
                                  ++running_update;
                                  // });
                                }
                              }
                            }
                            // , numUpdateL
      );

      // Until https://github.com/kokkos/kokkos/issues/6259 is resolved, do
      // reduction outside of parallel_scan
      team.team_reduce(Kokkos::Sum<ordinal_type, execution_space>(numUpdateL));

      updateIdx += numUpdateL;
    }

    Kokkos::single(Kokkos::PerTeam(team), [&] {
      row_mapU(factorization_step + 1) = U_entryIdx;
      row_mapL(factorization_step + 1) = L_entryIdx;

      update_list_len  = updateIdx;
      selected_row_len = rowView.length;

      factored(selected_row) = 1;
    });
  }
};

template <class col_ind_type>
struct MDF_reindex_matrix {
  col_ind_type permutation_inv;
  col_ind_type entries;

  MDF_reindex_matrix(col_ind_type permutation_inv_, col_ind_type entries_)
      : permutation_inv(permutation_inv_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int entryIdx) const { entries(entryIdx) = permutation_inv(entries(entryIdx)); }
};

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSSPARSE_MDF_IMPL_HPP_
