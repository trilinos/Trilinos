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

namespace KokkosSparse {
namespace Impl {

template <class crs_matrix_type>
struct MDF_count_lower {
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using size_type  = typename crs_matrix_type::ordinal_type;
  using value_type = typename crs_matrix_type::size_type;

  crs_matrix_type A;
  col_ind_type permutation;
  col_ind_type permutation_inv;

  MDF_count_lower(crs_matrix_type A_, col_ind_type permutation_,
                  col_ind_type permutation_inv_)
      : A(A_), permutation(permutation_), permutation_inv(permutation_inv_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type rowIdx, value_type& update) const {
    permutation(rowIdx)     = rowIdx;
    permutation_inv(rowIdx) = rowIdx;
    for (value_type entryIdx = A.graph.row_map(rowIdx);
         entryIdx < A.graph.row_map(rowIdx + 1); ++entryIdx) {
      if (A.graph.entries(entryIdx) <= rowIdx) {
        update += 1;
      }
    }
  }

};  // MDF_count_lower

template <class crs_matrix_type>
struct MDF_discarded_fill_norm {
  using static_crs_graph_type = typename crs_matrix_type::StaticCrsGraphType;
  using col_ind_type =
      typename static_crs_graph_type::entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using scalar_type  = typename crs_matrix_type::value_type;
  using KAS          = typename Kokkos::ArithTraits<scalar_type>;

  const scalar_type zero = KAS::zero();

  crs_matrix_type A, At;
  ordinal_type factorization_step;
  col_ind_type permutation;

  values_type discarded_fill;
  col_ind_type deficiency;
  int verbosity;

  MDF_discarded_fill_norm(crs_matrix_type A_, crs_matrix_type At_,
                          ordinal_type factorization_step_,
                          col_ind_type permutation_,
                          values_type discarded_fill_, col_ind_type deficiency_,
                          int verbosity_)
      : A(A_),
        At(At_),
        factorization_step(factorization_step_),
        permutation(permutation_),
        discarded_fill(discarded_fill_),
        deficiency(deficiency_),
        verbosity(verbosity_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type i) const {
    ordinal_type rowIdx      = permutation(i);
    scalar_type discard_norm = zero, diag_val = zero;
    bool entryIsDiscarded       = true;
    ordinal_type numFillEntries = 0;
    for (size_type alphaIdx = At.graph.row_map(rowIdx);
         alphaIdx < At.graph.row_map(rowIdx + 1); ++alphaIdx) {
      ordinal_type fillRowIdx = At.graph.entries(alphaIdx);
      bool row_not_eliminated = true;
      for (ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
        if (fillRowIdx == permutation(stepIdx)) {
          row_not_eliminated = false;
        }
      }

      if (fillRowIdx != rowIdx && row_not_eliminated) {
        for (size_type betaIdx = A.graph.row_map(rowIdx);
             betaIdx < A.graph.row_map(rowIdx + 1); ++betaIdx) {
          ordinal_type fillColIdx = A.graph.entries(betaIdx);
          bool col_not_eliminated = true;
          for (ordinal_type stepIdx = 0; stepIdx < factorization_step;
               ++stepIdx) {
            if (fillColIdx == permutation(stepIdx)) {
              col_not_eliminated = false;
            }
          }

          if (fillColIdx != rowIdx && col_not_eliminated) {
            entryIsDiscarded = true;
            for (size_type entryIdx = A.graph.row_map(fillRowIdx);
                 entryIdx < A.graph.row_map(fillRowIdx + 1); ++entryIdx) {
              if (A.graph.entries(entryIdx) == fillColIdx) {
                entryIsDiscarded = false;
              }
            }
            if (entryIsDiscarded) {
              numFillEntries += 1;
              discard_norm +=
                  KAS::abs(At.values(alphaIdx) * A.values(betaIdx)) *
                  KAS::abs(At.values(alphaIdx) * A.values(betaIdx));
              if (verbosity > 1) {
                KOKKOS_IMPL_DO_NOT_USE_PRINTF(
                    "Adding value A[%d,%d]=%f to discard norm of row %d\n",
                    int(At.graph.entries(alphaIdx)),
                    int(A.graph.entries(betaIdx)),
                    KAS::abs(At.values(alphaIdx) * A.values(betaIdx)) *
                        KAS::abs(At.values(alphaIdx) * A.values(betaIdx)),
                    int(rowIdx));
              }
            }
          }
        }
      } else if (fillRowIdx == rowIdx) {
        diag_val = At.values(alphaIdx);
        if (verbosity > 1) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "Row %d diagonal value dected, values(%d)=%f\n", int(rowIdx),
              int(alphaIdx), At.values(alphaIdx));
        }
      }
    }

    // TODO add a check on `diag_val == zero`
    discard_norm           = discard_norm / (diag_val * diag_val);
    discarded_fill(rowIdx) = discard_norm;
    deficiency(rowIdx)     = numFillEntries;
    if (verbosity > 0) {
      const ordinal_type degree = ordinal_type(A.graph.row_map(rowIdx + 1) -
                                               A.graph.row_map(rowIdx) - 1);
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "Row %d has discarded fill of %f, deficiency of %d and degree %d\n",
          static_cast<int>(rowIdx),
          static_cast<double>(KAS::sqrt(discard_norm)),
          static_cast<int>(deficiency(rowIdx)), static_cast<int>(degree));
    }
  }

};  // MDF_discarded_fill_norm

template <class crs_matrix_type>
struct MDF_selective_discarded_fill_norm {
  using static_crs_graph_type = typename crs_matrix_type::StaticCrsGraphType;
  using col_ind_type =
      typename static_crs_graph_type::entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using scalar_type  = typename crs_matrix_type::value_type;
  using KAS          = typename Kokkos::ArithTraits<scalar_type>;

  const scalar_type zero = KAS::zero();

  crs_matrix_type A, At;
  ordinal_type factorization_step;
  col_ind_type permutation;
  col_ind_type update_list;

  values_type discarded_fill;
  col_ind_type deficiency;
  int verbosity;

  MDF_selective_discarded_fill_norm(crs_matrix_type A_, crs_matrix_type At_,
                                    ordinal_type factorization_step_,
                                    col_ind_type permutation_,
                                    col_ind_type update_list_,
                                    values_type discarded_fill_,
                                    col_ind_type deficiency_, int verbosity_)
      : A(A_),
        At(At_),
        factorization_step(factorization_step_),
        permutation(permutation_),
        update_list(update_list_),
        discarded_fill(discarded_fill_),
        deficiency(deficiency_),
        verbosity(verbosity_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type i) const {
    ordinal_type rowIdx      = permutation(update_list(i));
    scalar_type discard_norm = zero, diag_val = zero;
    bool entryIsDiscarded       = true;
    ordinal_type numFillEntries = 0;
    for (size_type alphaIdx = At.graph.row_map(rowIdx);
         alphaIdx < At.graph.row_map(rowIdx + 1); ++alphaIdx) {
      ordinal_type fillRowIdx = At.graph.entries(alphaIdx);
      bool row_not_eliminated = true;
      for (ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
        if (fillRowIdx == permutation(stepIdx)) {
          row_not_eliminated = false;
        }
      }

      if (fillRowIdx != rowIdx && row_not_eliminated) {
        for (size_type betaIdx = A.graph.row_map(rowIdx);
             betaIdx < A.graph.row_map(rowIdx + 1); ++betaIdx) {
          ordinal_type fillColIdx = A.graph.entries(betaIdx);
          bool col_not_eliminated = true;
          for (ordinal_type stepIdx = 0; stepIdx < factorization_step;
               ++stepIdx) {
            if (fillColIdx == permutation(stepIdx)) {
              col_not_eliminated = false;
            }
          }

          if (fillColIdx != rowIdx && col_not_eliminated) {
            entryIsDiscarded = true;
            for (size_type entryIdx = A.graph.row_map(fillRowIdx);
                 entryIdx < A.graph.row_map(fillRowIdx + 1); ++entryIdx) {
              if (A.graph.entries(entryIdx) == fillColIdx) {
                entryIsDiscarded = false;
              }
            }
            if (entryIsDiscarded) {
              numFillEntries += 1;
              discard_norm +=
                  KAS::abs(At.values(alphaIdx) * A.values(betaIdx)) *
                  KAS::abs(At.values(alphaIdx) * A.values(betaIdx));
              if (verbosity > 1) {
                KOKKOS_IMPL_DO_NOT_USE_PRINTF(
                    "Adding value A[%d,%d]=%f to discard norm of row %d\n",
                    static_cast<int>(At.graph.entries(alphaIdx)),
                    static_cast<int>(A.graph.entries(betaIdx)),
                    static_cast<double>(
                        KAS::abs(At.values(alphaIdx) * A.values(betaIdx)) *
                        KAS::abs(At.values(alphaIdx) * A.values(betaIdx))),
                    static_cast<int>(rowIdx));
              }
            }
          }
        }
      } else if (fillRowIdx == rowIdx) {
        diag_val = At.values(alphaIdx);
        if (verbosity > 1) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "Row %d diagonal value dected, values(%d)=%f\n",
              static_cast<int>(rowIdx), static_cast<int>(alphaIdx),
              static_cast<double>(At.values(alphaIdx)));
        }
      }
    }

    // TODO add a check on `diag_val == zero`
    discard_norm           = discard_norm / (diag_val * diag_val);
    discarded_fill(rowIdx) = discard_norm;
    deficiency(rowIdx)     = numFillEntries;
    if (verbosity > 0) {
      const ordinal_type degree = ordinal_type(A.graph.row_map(rowIdx + 1) -
                                               A.graph.row_map(rowIdx) - 1);
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "Row %d has discarded fill of %f, deficiency of %d and degree %d\n",
          static_cast<int>(rowIdx),
          static_cast<double>(KAS::sqrt(discard_norm)),
          static_cast<int>(deficiency(rowIdx)), static_cast<int>(degree));
    }
  }

};  // MDF_selective_discarded_fill_norm

template <class crs_matrix_type>
struct MDF_select_row {
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using row_map_type =
      typename crs_matrix_type::StaticCrsGraphType::row_map_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using scalar_type  = typename crs_matrix_type::value_type;

  // type used to perform the reduction
  // do not confuse it with scalar_type!
  using value_type = typename crs_matrix_type::ordinal_type;

  value_type factorization_step;
  values_type discarded_fill;
  col_ind_type deficiency;
  row_map_type row_map;
  col_ind_type permutation;

  MDF_select_row(value_type factorization_step_, values_type discarded_fill_,
                 col_ind_type deficiency_, row_map_type row_map_,
                 col_ind_type permutation_)
      : factorization_step(factorization_step_),
        discarded_fill(discarded_fill_),
        deficiency(deficiency_),
        row_map(row_map_),
        permutation(permutation_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type src, ordinal_type& dst) const {
    const ordinal_type src_perm = permutation(src);
    const ordinal_type dst_perm = permutation(dst);
    const ordinal_type degree_src =
        row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst =
        row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if (discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src == degree_dst) && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const {
    const ordinal_type src_perm = permutation(src);
    const ordinal_type dst_perm = permutation(dst);
    const ordinal_type degree_src =
        row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst =
        row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if (discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if ((discarded_fill(src_perm) == discarded_fill(dst_perm)) &&
        (deficiency(src_perm) == deficiency(dst_perm)) &&
        (degree_src == degree_dst) && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& dst) const {
    dst = Kokkos::ArithTraits<ordinal_type>::zero();
  }

};  // MDF_select_row

template <class crs_matrix_type>
struct MDF_factorize_row {
  using row_map_type = typename crs_matrix_type::StaticCrsGraphType::
      row_map_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using size_type    = typename crs_matrix_type::size_type;
  using value_type   = typename crs_matrix_type::value_type;

  crs_matrix_type A, At;

  row_map_type row_mapL;
  col_ind_type entriesL;
  values_type valuesL;

  row_map_type row_mapU;
  col_ind_type entriesU;
  values_type valuesU;

  col_ind_type permutation, permutation_inv;
  values_type discarded_fill;
  col_ind_type factored;
  ordinal_type selected_row_idx, factorization_step;

  int verbosity;

  MDF_factorize_row(crs_matrix_type A_, crs_matrix_type At_,
                    row_map_type row_mapL_, col_ind_type entriesL_,
                    values_type valuesL_, row_map_type row_mapU_,
                    col_ind_type entriesU_, values_type valuesU_,
                    col_ind_type permutation_, col_ind_type permutation_inv_,
                    values_type discarded_fill_, col_ind_type factored_,
                    ordinal_type selected_row_idx_,
                    ordinal_type factorization_step_, int verbosity_)
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
        discarded_fill(discarded_fill_),
        factored(factored_),
        selected_row_idx(selected_row_idx_),
        factorization_step(factorization_step_),
        verbosity(verbosity_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type /* idx */) const {
    const ordinal_type selected_row = permutation(selected_row_idx);
    discarded_fill(selected_row)    = Kokkos::ArithTraits<value_type>::max();

    // Swap entries in permutation vectors
    permutation(selected_row_idx)   = permutation(factorization_step);
    permutation(factorization_step) = selected_row;
    permutation_inv(permutation(factorization_step)) = factorization_step;
    permutation_inv(permutation(selected_row_idx))   = selected_row_idx;

    if (verbosity > 0) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("Permutation vector: { ");
      for (ordinal_type rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%d ",
                                      static_cast<int>(permutation(rowIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}\n");
    }

    // Insert the upper part of the selected row in U
    // including the diagonal term.
    value_type diag      = Kokkos::ArithTraits<value_type>::zero();
    size_type U_entryIdx = row_mapU(factorization_step);
    for (size_type entryIdx = A.graph.row_map(selected_row);
         entryIdx < A.graph.row_map(selected_row + 1); ++entryIdx) {
      if (permutation_inv(A.graph.entries(entryIdx)) >= factorization_step) {
        entriesU(U_entryIdx) = A.graph.entries(entryIdx);
        valuesU(U_entryIdx)  = A.values(entryIdx);
        ++U_entryIdx;
        if (A.graph.entries(entryIdx) == selected_row) {
          diag = A.values(entryIdx);
        }
      }
    }
    row_mapU(factorization_step + 1) = U_entryIdx;

    if (verbosity > 0) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("Diagonal values of row %d is %f\n",
                                    static_cast<int>(selected_row),
                                    static_cast<double>(diag));
    }

    if (verbosity > 2) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("U, row_map={ ");
      for (ordinal_type rowIdx = 0; rowIdx < factorization_step + 1; ++rowIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%d ",
                                      static_cast<int>(row_mapU(rowIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}, entries={ ");
      for (size_type entryIdx = row_mapU(0);
           entryIdx < row_mapU(factorization_step + 1); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%d ",
                                      static_cast<int>(entriesU(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}, values={ ");
      for (size_type entryIdx = row_mapU(0);
           entryIdx < row_mapU(factorization_step + 1); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f ",
                                      static_cast<double>(valuesU(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}\n");
    }

    // Insert the lower part of the selected column of A
    // divided by its the diagonal value to obtain a unit
    // diagonal value in L.
    size_type L_entryIdx = row_mapL(factorization_step);
    entriesL(L_entryIdx) = selected_row;
    valuesL(L_entryIdx)  = Kokkos::ArithTraits<value_type>::one();
    ++L_entryIdx;
    for (size_type entryIdx = At.graph.row_map(selected_row);
         entryIdx < At.graph.row_map(selected_row + 1); ++entryIdx) {
      if (permutation_inv(At.graph.entries(entryIdx)) > factorization_step) {
        entriesL(L_entryIdx) = At.graph.entries(entryIdx);
        valuesL(L_entryIdx)  = At.values(entryIdx) / diag;
        ++L_entryIdx;
      }
    }
    row_mapL(factorization_step + 1) = L_entryIdx;

    if (verbosity > 2) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "L(%d), [row_map(%d), row_map(%d)[ = [%d, %d[, entries={ ",
          static_cast<int>(factorization_step),
          static_cast<int>(factorization_step),
          static_cast<int>(factorization_step + 1),
          static_cast<int>(row_mapL(factorization_step)),
          static_cast<int>(row_mapL(factorization_step + 1)));
      for (size_type entryIdx = row_mapL(factorization_step);
           entryIdx < row_mapL(factorization_step + 1); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%d ",
                                      static_cast<int>(entriesL(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}, values={ ");
      for (size_type entryIdx = row_mapL(factorization_step);
           entryIdx < row_mapL(factorization_step + 1); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f ",
                                      static_cast<double>(valuesL(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}\n");
    }

    // If this was the last row no need to update A and At!
    if (factorization_step == A.numRows() - 1) {
      return;
    }

    // Finally we want to update A and At with the values
    // that where not discarded during factorization.
    // Note: this is almost the same operation as computing
    // the norm of the discarded fill...

    // First step: find the diagonal entry in selected_row
    value_type diag_val = Kokkos::ArithTraits<value_type>::zero();
    for (size_type entryIdx = A.graph.row_map(selected_row);
         entryIdx < A.graph.row_map(selected_row + 1); ++entryIdx) {
      ordinal_type colIdx = A.graph.entries(entryIdx);
      if (selected_row == colIdx) {
        diag_val = A.values(entryIdx);
      }
    }

    // Extract alpha and beta vectors
    // Then insert alpha*beta/diag_val if the corresponding
    // entry in A is non-zero.
    for (size_type alphaIdx = At.graph.row_map(selected_row);
         alphaIdx < At.graph.row_map(selected_row + 1); ++alphaIdx) {
      ordinal_type fillRowIdx = At.graph.entries(alphaIdx);
      bool row_not_eliminated = true;
      for (ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
        if (fillRowIdx == permutation(stepIdx)) {
          row_not_eliminated = false;
        }
      }

      if ((fillRowIdx != selected_row) && row_not_eliminated) {
        for (size_type betaIdx = A.graph.row_map(selected_row);
             betaIdx < A.graph.row_map(selected_row + 1); ++betaIdx) {
          ordinal_type fillColIdx = A.graph.entries(betaIdx);
          bool col_not_eliminated = true;
          for (ordinal_type stepIdx = 0; stepIdx < factorization_step;
               ++stepIdx) {
            if (fillColIdx == permutation(stepIdx)) {
              col_not_eliminated = false;
            }
          }

          if ((fillColIdx != selected_row) && col_not_eliminated) {
            for (size_type entryIdx = A.graph.row_map(fillRowIdx);
                 entryIdx < A.graph.row_map(fillRowIdx + 1); ++entryIdx) {
              if (A.graph.entries(entryIdx) == fillColIdx) {
                A.values(entryIdx) -=
                    At.values(alphaIdx) * A.values(betaIdx) / diag_val;

                if (verbosity > 1) {
                  KOKKOS_IMPL_DO_NOT_USE_PRINTF(
                      "A[%d, %d] -= %f\n", static_cast<int>(fillRowIdx),
                      static_cast<int>(fillColIdx),
                      static_cast<double>(At.values(alphaIdx) *
                                          A.values(betaIdx) / diag_val));
                }
              }
            }

            for (size_type entryIdx = At.graph.row_map(fillColIdx);
                 entryIdx < At.graph.row_map(fillColIdx + 1); ++entryIdx) {
              if (At.graph.entries(entryIdx) == fillRowIdx) {
                At.values(entryIdx) -=
                    At.values(alphaIdx) * A.values(betaIdx) / diag_val;
              }
            }
          }
        }
      }
    }

    factored(selected_row) = 1;

    if (verbosity > 0) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("New values in A: { ");
      for (size_type entryIdx = 0; entryIdx < A.nnz(); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f ",
                                      static_cast<double>(A.values(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}\n");
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("New values in At: { ");
      for (size_type entryIdx = 0; entryIdx < At.nnz(); ++entryIdx) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("%f ",
                                      static_cast<double>(At.values(entryIdx)));
      }
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("}\n");
    }
  }  // operator()

};  // MDF_factorize_row

template <class crs_matrix_type>
struct MDF_compute_list_length {
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::
      entries_type::non_const_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using size_type    = typename crs_matrix_type::size_type;

  ordinal_type selected_row_idx;
  crs_matrix_type A;
  crs_matrix_type At;
  col_ind_type permutation;
  col_ind_type factored;
  col_ind_type update_list_length;
  col_ind_type update_list;

  MDF_compute_list_length(const ordinal_type rowIdx_, const crs_matrix_type& A_,
                          const crs_matrix_type& At_,
                          const col_ind_type& permutation_,
                          const col_ind_type factored_,
                          col_ind_type& update_list_length_,
                          col_ind_type& update_list_)
      : selected_row_idx(rowIdx_),
        A(A_),
        At(At_),
        permutation(permutation_),
        factored(factored_),
        update_list_length(update_list_length_),
        update_list(update_list_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type /*idx*/) const {
    const ordinal_type selected_row = permutation(selected_row_idx);

    size_type updateIdx = 0;
    for (size_type entryIdx = A.graph.row_map(selected_row);
         entryIdx < A.graph.row_map(selected_row + 1); ++entryIdx) {
      if ((A.graph.entries(entryIdx) != selected_row) &&
          (factored(A.graph.entries(entryIdx)) != 1)) {
        update_list(updateIdx) = A.graph.entries(entryIdx);
        ++updateIdx;
      }
    }
    size_type update_rows = updateIdx;
    for (size_type entryIdx = At.graph.row_map(selected_row);
         entryIdx < At.graph.row_map(selected_row + 1); ++entryIdx) {
      if ((At.graph.entries(entryIdx) != selected_row) &&
          (factored(A.graph.entries(entryIdx)) != 1)) {
        bool already_updated = false;
        for (size_type checkIdx = 0; checkIdx < update_rows; ++checkIdx) {
          if (At.graph.entries(entryIdx) == update_list(checkIdx)) {
            already_updated = true;
            break;
          }
        }
        if (already_updated == false) {
          update_list(updateIdx) = At.graph.entries(entryIdx);
          ++updateIdx;
        }
      }
    }
    update_list_length(0) = updateIdx;
  }
};

template <class col_ind_type>
struct MDF_reindex_matrix {
  col_ind_type permutation_inv;
  col_ind_type entries;

  MDF_reindex_matrix(col_ind_type permutation_inv_, col_ind_type entries_)
      : permutation_inv(permutation_inv_), entries(entries_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int entryIdx) const {
    entries(entryIdx) = permutation_inv(entries(entryIdx));
  }
};

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSSPARSE_MDF_IMPL_HPP_
