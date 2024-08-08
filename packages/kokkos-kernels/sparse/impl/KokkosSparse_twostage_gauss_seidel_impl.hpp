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

#ifndef _KOKKOS_TWOSTAGE_GS_IMP_HPP
#define _KOKKOS_TWOSTAGE_GS_IMP_HPP

#include <KokkosKernels_config.h>
#include "Kokkos_Core.hpp"

// Blas Kernels
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_mult.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosKernels_SimpleUtils.hpp"

// needed for classical GS
#include "KokkosSparse_sptrsv.hpp"
#include "KokkosSparse_Utils.hpp"

#include "KokkosSparse_gauss_seidel_handle.hpp"

// #define KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS

namespace KokkosSparse {
namespace Impl {

template <typename HandleType, typename input_row_map_view_t, typename input_entries_view_t,
          typename input_values_view_t>
class TwostageGaussSeidel {
 public:
  using TwoStageGaussSeidelHandleType = typename HandleType::TwoStageGaussSeidelHandleType;
  using execution_space               = typename HandleType::HandleExecSpace;
  using memory_space                  = typename TwoStageGaussSeidelHandleType::memory_space;

  using const_scalar_t  = typename TwoStageGaussSeidelHandleType::const_scalar_t;
  using scalar_t        = typename TwoStageGaussSeidelHandleType::scalar_t;
  using const_ordinal_t = typename TwoStageGaussSeidelHandleType::const_ordinal_t;
  using ordinal_t       = typename TwoStageGaussSeidelHandleType::ordinal_t;
  using size_type       = typename TwoStageGaussSeidelHandleType::size_type;

  using const_row_map_view_t = typename TwoStageGaussSeidelHandleType::const_row_map_view_t;
  using row_map_view_t       = typename TwoStageGaussSeidelHandleType::row_map_view_t;
  using entries_view_t       = typename TwoStageGaussSeidelHandleType::entries_view_t;
  using values_view_t        = typename TwoStageGaussSeidelHandleType::values_view_t;

  using crsmat_t = typename TwoStageGaussSeidelHandleType::crsmat_t;
  using graph_t  = typename TwoStageGaussSeidelHandleType::graph_t;

  using range_type = Kokkos::pair<int, int>;

  // to wrap input (rowmap, colind, values) into crsmat
  using input_device_t =
      Kokkos::Device<typename input_row_map_view_t::execution_space, typename input_row_map_view_t::memory_space>;
  using input_memory_t  = typename input_values_view_t::memory_traits;
  using input_scalar_t  = typename input_values_view_t::value_type;
  using input_ordinal_t = typename input_entries_view_t::value_type;
  using input_size_t    = typename input_row_map_view_t::value_type;
  using input_crsmat_t =
      KokkosSparse::CrsMatrix<input_scalar_t, input_ordinal_t, input_device_t, input_memory_t, input_size_t>;
  using input_graph_t          = typename input_crsmat_t::StaticCrsGraphType;
  using single_vector_view_t   = Kokkos::View<scalar_t *, Kokkos::LayoutLeft, input_device_t, input_memory_t>;
  using internal_vector_view_t = typename TwoStageGaussSeidelHandleType::vector_view_t;

  using ST    = Kokkos::ArithTraits<scalar_t>;
  using mag_t = typename ST::mag_type;

 private:
  HandleType *handle;

  // Get the specialized TwostageGaussSeidel handle from the main handle
  TwoStageGaussSeidelHandleType *get_gs_handle() {
    auto gsHandle = dynamic_cast<TwoStageGaussSeidelHandleType *>(this->handle->get_gs_handle());
    if (!gsHandle) {
      throw std::runtime_error(
          "TwostageGaussSeidel: GS handle has not been created, or is set up "
          "for Cluster GS.");
    }
    return gsHandle;
  }

  bool diagos_given;
  const_ordinal_t num_rows, num_cols;

  input_row_map_view_t rowmap_view;
  input_entries_view_t column_view;
  input_values_view_t values_view;
  input_values_view_t d_invert_view;

  // --------------------------------------------------------- //
 public:
  // tag for counting nnz
  struct Tag_countNnzL {};
  struct Tag_countNnzU {};
  // tag for inserting entries
  struct Tag_entriesLU {};
  // tag for inserting values
  struct Tag_valuesLU {};
  // tag for computing residual norm
  struct Tag_normR {};

  template <typename output_row_map_view_t, typename output_entries_view_t, typename output_values_view_t>
  struct TwostageGaussSeidel_functor {
   public:
    // input
    bool two_stage;
    bool compact_form;
    bool diagos_given;
    const_ordinal_t num_rows;
    input_row_map_view_t rowmap_view;
    input_entries_view_t column_view;
    input_values_view_t values_view;
    input_values_view_t d_invert_view;
    // output (lower)
    output_row_map_view_t row_map;
    output_entries_view_t entries;
    output_values_view_t values;
    output_values_view_t diags;
    // output (upper)
    output_row_map_view_t row_map2;
    output_entries_view_t entries2;
    output_values_view_t values2;
    // output (complement of U+D)
    output_row_map_view_t row_map_a;
    output_entries_view_t entries_a;
    output_values_view_t values_a;
    output_values_view_t diags_a;
    // output (complement of L+D)
    output_row_map_view_t row_map_a2;
    output_entries_view_t entries_a2;
    output_values_view_t values_a2;
    // for computing residual norm with
    bool forward_sweep;
    internal_vector_view_t localX;
    internal_vector_view_t localB;

    // ------------------------------------------------------- //
    // Constructors
    // for counting nnz
    TwostageGaussSeidel_functor(bool two_stage_, bool compact_form_, const_ordinal_t num_rows_,
                                input_row_map_view_t rowmap_view_, input_entries_view_t column_view_,
                                output_row_map_view_t row_map_, output_row_map_view_t row_map_a_)
        : two_stage(two_stage_),
          compact_form(compact_form_),
          diagos_given(false),
          num_rows(num_rows_),
          // input
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(),
          d_invert_view(),
          // output
          row_map(row_map_),
          entries(),
          values(),
          diags(),
          row_map2(),
          entries2(),
          values2(),
          // output (complement)
          row_map_a(row_map_a_) {}

    // for storing booth L&U entries
    TwostageGaussSeidel_functor(bool two_stage_, bool compact_form_, const_ordinal_t num_rows_,
                                input_row_map_view_t rowmap_view_, input_entries_view_t column_view_,
                                output_row_map_view_t row_map_, output_entries_view_t entries_,
                                output_row_map_view_t row_map2_, output_entries_view_t entries2_,
                                output_row_map_view_t row_map_a_, output_entries_view_t entries_a_,
                                output_row_map_view_t row_map_a2_, output_entries_view_t entries_a2_)
        : two_stage(two_stage_),
          compact_form(compact_form_),
          diagos_given(false),
          num_rows(num_rows_),
          // input Crs
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(),
          // input diag
          d_invert_view(),
          // output CrsL
          row_map(row_map_),
          entries(entries_),
          values(),
          // output diag
          diags(),
          // output CrsU
          row_map2(row_map2_),
          entries2(entries2_),
          values2(),
          // output complement of U+D
          row_map_a(row_map_a_),
          entries_a(entries_a_),
          values_a(),
          // output complement of L+D
          row_map_a2(row_map_a2_),
          entries_a2(entries_a2_) {}

    // for storing both L&U values (with D extracted)
    TwostageGaussSeidel_functor(bool two_stage_, bool compact_form_, bool diagos_given_, const_ordinal_t num_rows_,
                                input_row_map_view_t rowmap_view_, input_entries_view_t column_view_,
                                input_values_view_t values_view_, input_values_view_t d_invert_view_,
                                output_row_map_view_t row_map_, output_values_view_t values_,
                                output_values_view_t diags_, output_row_map_view_t row_map2_,
                                output_values_view_t values2_, output_row_map_view_t row_map_a_,
                                output_values_view_t values_a_, output_values_view_t diags_a_,
                                output_row_map_view_t row_map_a2_, output_values_view_t values_a2_)
        : two_stage(two_stage_),
          compact_form(compact_form_),
          diagos_given(diagos_given_),
          num_rows(num_rows_),
          // input Crs
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(values_view_),
          // input diag
          d_invert_view(d_invert_view_),
          // output CrsL
          row_map(row_map_),
          entries(),
          values(values_),
          // output diag
          diags(diags_),
          // output CrsU
          row_map2(row_map2_),
          entries2(),
          values2(values2_),
          // output complement of U
          row_map_a(row_map_a_),
          entries_a(),
          values_a(values_a_),
          diags_a(diags_a_),
          // output complement of L
          row_map_a2(row_map_a2_),
          entries_a2(),
          values_a2(values_a2_) {}

    // for computing residual norm
    TwostageGaussSeidel_functor(bool forward_sweep_, const_ordinal_t num_rows_, input_row_map_view_t rowmap_view_,
                                input_entries_view_t column_view_, input_values_view_t values_view_,
                                output_values_view_t diags_, internal_vector_view_t localX_,
                                internal_vector_view_t localB_)
        : two_stage(false),
          compact_form(false),
          diagos_given(false),
          num_rows(num_rows_),
          // input Crs
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(values_view_),
          d_invert_view(),
          row_map(),
          entries(),
          values(),
          diags(diags_),
          // input vectors
          forward_sweep(forward_sweep_),
          localX(localX_),
          localB(localB_) {}

    // ------------------------------------------------------- //
    // functor for counting nnzL (with parallel_reduce)
    KOKKOS_INLINE_FUNCTION
    void operator()(const Tag_countNnzL &, const ordinal_t i, ordinal_t &nnz) const {
      ordinal_t nnz_i = 0;
      for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
        if (column_view(k) < i) {
          nnz_i++;
        } else if (!two_stage && column_view(k) == i) {
          nnz_i++;
        }
      }
      row_map(i + 1) = nnz_i;
      if (i == 0) {
        row_map(0) = 0;
      }
      if (compact_form) {
        // complement of L+D
        row_map_a(i + 1) = (rowmap_view(i + 1) - rowmap_view(i)) - nnz_i;
        if (two_stage) {
          // two-stage iterates with L (no D)
          row_map_a(i + 1)--;
        }
        if (i == 0) {
          row_map_a(0) = 0;
        }
      }
      nnz += nnz_i;
    }

    // ------------------------------------------------------- //
    // functor for counting nnzU (with parallel_reduce)
    KOKKOS_INLINE_FUNCTION
    void operator()(const Tag_countNnzU &, const ordinal_t i, ordinal_t &nnz) const {
      ordinal_t nnz_i = 0;
      for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
        if (column_view(k) > i && column_view(k) < num_rows) {
          nnz_i++;
        } else if (!two_stage && column_view(k) == i) {
          nnz_i++;
        }
      }
      row_map(i + 1) = nnz_i;
      if (i == 0) {
        row_map(0) = 0;
      }
      if (compact_form) {
        // complement of U+D
        row_map_a(i + 1) = (rowmap_view(i + 1) - rowmap_view(i)) - nnz_i;
        if (two_stage) {
          // two-stage iterates with U (no D)
          row_map_a(i + 1)--;
        }
        if (i == 0) {
          row_map_a(0) = 0;
        }
      }
      nnz += nnz_i;
    }

    // ------------------------------------------------------- //
    // functor for storing entriesL and entriesU (with parallel_for)
    KOKKOS_INLINE_FUNCTION
    void operator()(const Tag_entriesLU &, const ordinal_t i) const {
      ordinal_t nnzL  = row_map(i);
      ordinal_t nnzU  = row_map2(i);
      ordinal_t nnzLa = 0;
      ordinal_t nnzUa = 0;
      if (compact_form) {
        nnzLa = row_map_a(i);
        nnzUa = row_map_a2(i);
      }
      if (!two_stage) {
        // NOTE: Kokkos' sptrsv assumes diagonal of U to be at the start
        entries2(nnzU) = i;
        nnzU++;
      }
      for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
        if (column_view(k) < i) {
          // L
          entries(nnzL) = column_view(k);
          nnzL++;
          if (compact_form) {
            // complement of U+D
            entries_a(nnzLa) = column_view(k);
            nnzLa++;
          }
        } else if (column_view(k) > i) {
          if (column_view(k) < num_rows) {
            // U
            entries2(nnzU) = column_view(k);
            nnzU++;
            if (compact_form) {
              // complement of L+D
              entries_a2(nnzUa) = column_view(k);
              nnzUa++;
            }
          } else if (compact_form) {
            // complement of U+D
            entries_a(nnzLa) = column_view(k);
            nnzLa++;
            // complement of L+D
            entries_a2(nnzUa) = column_view(k);
            nnzUa++;
          }
        }
      }
      if (!two_stage) {
        // NOTE: Kokkos' sptrsv assumes diagonal of L to be at the end
        entries(nnzL) = i;
        nnzL++;
      }
    }

    // functor for storing both valuesL & valuesU (with parallel_for)
    KOKKOS_INLINE_FUNCTION
    void operator()(const Tag_valuesLU &, const ordinal_t i) const {
      const_scalar_t one = Kokkos::ArithTraits<scalar_t>::one();
      ordinal_t nnzL     = row_map(i);
      ordinal_t nnzU     = row_map2(i);
      ordinal_t nnzLa    = 0;
      ordinal_t nnzUa    = 0;
      if (compact_form) {
        nnzLa = row_map_a(i);
        nnzUa = row_map_a2(i);
      }
      if (!two_stage) {
        // Kokkos' sptrsv assumes diagonal U to come at the start, so increment
        // nnzU
        nnzU++;
      }
      for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
        if (column_view(k) < i) {
          // save L (without diag)
          values(nnzL) = values_view(k);
          nnzL++;
          if (compact_form) {
            // complement of U+D
            values_a(nnzLa) = values_view(k);
            nnzLa++;
          }
        } else if (column_view(k) == i) {
          // save D
          if (diagos_given) {
            // as inverse
            diags(i) = d_invert_view(i);
          } else {
            // as original
            diags(i) = values_view(k);
          }
          if (compact_form) {
            diags_a(i) = values_view(k);
          }
        } else {
          if (column_view(k) < num_rows) {
            // save U (without diag)
            values2(nnzU) = values_view(k);
            nnzU++;
            if (compact_form) {
              // complement of L+D
              values_a2(nnzUa) = values_view(k);
              nnzUa++;
            }
          } else if (compact_form) {
            // complement of U+D
            values_a(nnzLa) = values_view(k);
            nnzLa++;
            // complement of L+D
            values_a2(nnzUa) = values_view(k);
            nnzUa++;
          }
        }
      }
      if (!two_stage) {
        // if using sptrsv, add diagonals in L and U
        // > Kokkos' sptrsv assumes diagonal of L and U to come at end and start
        nnzU = row_map2(i);
        if (diagos_given) {
          values2(nnzU) = one / diags(i);
          values(nnzL)  = one / diags(i);
        } else {
          values2(nnzU) = diags(i);
          values(nnzL)  = diags(i);
        }
      }
      if (two_stage) {
        if (!diagos_given) {
          // when diag is provided, it is already provided as inverse
          diags(i) = one / diags(i);
        }
        // compute inv(D)*L (apply row-scaling to valueL)
        for (size_type k = row_map(i); k < row_map(i + 1); k++) {
          values(k) *= diags(i);
        }
        // compute inv(D)*U (apply row-scaling to valueU)
        for (size_type k = row_map2(i); k < row_map2(i + 1); k++) {
          values2(k) *= diags(i);
        }
      }
    }

    // ------------------------------------------------------- //
    // functor for computing residual norm (with parallel_reduce)
    KOKKOS_INLINE_FUNCTION
    void operator()(const Tag_normR &, const ordinal_t i, mag_t &normR) const {
      scalar_t normRi = localB(i, 0);
      if (forward_sweep) {
        // compute R(i) = B(i) - (L+D)(i,:)*X
        for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
          if (column_view(k) <= i) {
            normRi -= values_view(k) * localX(column_view(k), 0);
          }
        }
      } else {
        // compute R(i) = B(i) - (D+U)(i,:)*X
        for (size_type k = rowmap_view(i); k < rowmap_view(i + 1); k++) {
          if (column_view(k) >= i && column_view(k) < num_rows) {
            normRi -= values_view(k) * localX(column_view(k), 0);
          }
        }
      }
      normR += ST::abs(normRi * normRi);
    }
  };
  // --------------------------------------------------------- //

 public:
  /**
   * \brief constructor
   */
  // for symbolic (wihout values)
  TwostageGaussSeidel(HandleType *handle_, const_ordinal_t num_rows_, const_ordinal_t num_cols_,
                      input_row_map_view_t rowmap_view_, input_entries_view_t column_view_)
      : handle(handle_),
        diagos_given(false),
        num_rows(num_rows_),
        num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(),
        d_invert_view() {}

  // for numeric/solve (with values)
  TwostageGaussSeidel(HandleType *handle_, const_ordinal_t num_rows_, const_ordinal_t num_cols_,
                      input_row_map_view_t rowmap_view_, input_entries_view_t column_view_,
                      input_values_view_t values_view_)
      : handle(handle_),
        diagos_given(false),
        num_rows(num_rows_),
        num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(values_view_),
        d_invert_view() {}

  // for numeric/solve (with values and diagonal)
  TwostageGaussSeidel(HandleType *handle_, const_ordinal_t num_rows_, const_ordinal_t num_cols_,
                      input_row_map_view_t rowmap_view_, input_entries_view_t column_view_,
                      input_values_view_t values_view_, input_values_view_t d_invert_view_)
      : handle(handle_),
        diagos_given(true),
        num_rows(num_rows_),
        num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(values_view_),
        d_invert_view(d_invert_view_) {}

  /**
   * Symbolic setup
   */
  void initialize_symbolic() {
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    double tic;
    Kokkos::Timer timer;
    Kokkos::fence();
    tic = timer.seconds();
#endif
    auto *gsHandle        = get_gs_handle();
    bool two_stage        = gsHandle->isTwoStage();
    bool compact_form     = gsHandle->isCompactForm();
    GSDirection direction = gsHandle->getSweepDirection();
    using GS_Functor_t    = TwostageGaussSeidel_functor<row_map_view_t, entries_view_t, values_view_t>;
    // count nnz in local L & U matrices (rowmap_viewL/rowmap_viewU stores
    // offsets for each row)
    ordinal_t nnzA = column_view.extent(0);
    ordinal_t nnzL = 0;  // lower-part of diagonal block
    ordinal_t nnzU = 0;  // upper-part of diagonal block
    row_map_view_t rowmap_viewL("row_mapL",
                                num_rows + 1);  // lower-part of diagonal block
    row_map_view_t rowmap_viewU("row_mapU",
                                num_rows + 1);  // upper-part of diagonal block
    row_map_view_t rowmap_viewLa("row_mapLa",
                                 num_rows + 1);  // complement of U+D
    row_map_view_t rowmap_viewUa("row_mapUa",
                                 num_rows + 1);  // complement of L+D
    if (direction == GS_FORWARD || direction == GS_SYMMETRIC) {
      using range_policy = Kokkos::RangePolicy<Tag_countNnzL, execution_space>;
      Kokkos::parallel_reduce(
          "nnzL", range_policy(0, num_rows),
          GS_Functor_t(two_stage, compact_form, num_rows, rowmap_view, column_view, rowmap_viewL, rowmap_viewUa), nnzL);
    }
    if (direction == GS_BACKWARD || direction == GS_SYMMETRIC) {
      using range_policy = Kokkos::RangePolicy<Tag_countNnzU, execution_space>;
      Kokkos::parallel_reduce(
          "nnzU", range_policy(0, num_rows),
          GS_Functor_t(two_stage, compact_form, num_rows, rowmap_view, column_view, rowmap_viewU, rowmap_viewLa), nnzU);
    }
    ordinal_t nnzLa = 0;  // complement of U+D
    ordinal_t nnzUa = 0;  // complement of L+D
    if (compact_form) {
      nnzLa = nnzA - nnzU;
      nnzUa = nnzA - nnzL;
      if (two_stage) {
        // two-stage iterates with L or U (no D)
        nnzLa -= num_rows;
        nnzUa -= num_rows;
      }
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << std::endl << "TWO-STAGE GS::SYMBOLIC::COUNT-NNZ TIME : " << tic << std::endl;
    timer.reset();
#endif
    // shift ptr so that it now contains offsets (combine it with the previous
    // functor calls?)
    if (direction == GS_FORWARD || direction == GS_SYMMETRIC) {
      KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>(1 + num_rows, rowmap_viewL);
      if (compact_form) {
        KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>(1 + num_rows, rowmap_viewLa);
      }
    }
    if (direction == GS_BACKWARD || direction == GS_SYMMETRIC) {
      KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>(1 + num_rows, rowmap_viewU);
      if (compact_form) {
        KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>(1 + num_rows, rowmap_viewUa);
      }
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << "TWO-STAGE GS::SYMBOLIC::COMP-PTR TIME  : " << tic << std::endl;
    timer.reset();
#endif
    // allocate memory to store local D
    values_view_t viewD(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diags"), num_rows);

    // allocate memory to store local L
    entries_view_t column_viewL(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesL"), nnzL);
    values_view_t values_viewL(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesL"), nnzL);

    // allocate memory to store local U
    entries_view_t column_viewU(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesU"), nnzU);
    values_view_t values_viewU(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesU"), nnzU);

    // allocate memory to store complement of U+D
    entries_view_t column_viewLa(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesLa"), nnzLa);
    values_view_t values_viewLa(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesLa"), nnzLa);

    // allocate memory to store complement of L+D
    entries_view_t column_viewUa(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesUa"), nnzUa);
    values_view_t values_viewUa(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesUa"), nnzUa);
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << "TWO-STAGE GS::SYMBOLIC::ALLOCATE TIME  : " << tic << std::endl;
    timer.reset();
#endif

    {
      // extract local L & U structures (for computing (L+D)^{-1} or (D+U)^{-1})
      using range_policy = Kokkos::RangePolicy<Tag_entriesLU, execution_space>;
      Kokkos::parallel_for("entriesLU", range_policy(0, num_rows),
                           GS_Functor_t(two_stage, compact_form, num_rows, rowmap_view, column_view, rowmap_viewL,
                                        column_viewL, rowmap_viewU, column_viewU,
                                        //
                                        rowmap_viewLa, column_viewLa, rowmap_viewUa, column_viewUa));
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << "TWO-STAGE GS::SYMBOLIC::INSERT TIME    : " << tic << std::endl;
    timer.reset();
#endif

    // construct CrsMat with them
    graph_t graphL(column_viewL, rowmap_viewL);
    graph_t graphU(column_viewU, rowmap_viewU);
    crsmat_t crsmatL("L", num_rows, values_viewL, graphL);
    crsmat_t crsmatU("U", num_rows, values_viewU, graphU);

    // store them in handle
    gsHandle->setL(crsmatL);
    gsHandle->setU(crsmatU);
    gsHandle->setD(viewD);

    if (compact_form) {
      // construct complements
      graph_t graphLa(column_viewLa, rowmap_viewLa);
      graph_t graphUa(column_viewUa, rowmap_viewUa);
      crsmat_t crsmatLa("La", num_rows, values_viewLa, graphLa);
      crsmat_t crsmatUa("Ua", num_rows, values_viewUa, graphUa);

      // store them in handle
      gsHandle->setLa(crsmatLa);
      gsHandle->setUa(crsmatUa);

      values_view_t viewDa(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diags"), num_rows);
      gsHandle->setDa(viewDa);
    }

    if (!(gsHandle->isTwoStage())) {
      // create SpTRSV handles for classical GS
      using namespace KokkosSparse::Experimental;
      auto sptrsv_algo = handle->get_gs_sptrsvL_handle()->get_sptrsv_handle()->get_algorithm();
      if (sptrsv_algo != SPTRSVAlgorithm::SPTRSV_CUSPARSE) {  // symbolic with CuSparse needs
                                                              // values
        sptrsv_symbolic(handle->get_gs_sptrsvL_handle(), rowmap_viewL, crsmatL.graph.entries);
        sptrsv_symbolic(handle->get_gs_sptrsvU_handle(), rowmap_viewU, crsmatU.graph.entries);
      }
    }
  }

  /**
   * Numerical setup
   */
  void initialize_numeric() {
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    double tic;
    Kokkos::Timer timer;
    Kokkos::fence();
    timer.reset();
#endif
    using GS_Functor_t = TwostageGaussSeidel_functor<const_row_map_view_t, entries_view_t, values_view_t>;

    auto *gsHandle    = get_gs_handle();
    bool two_stage    = gsHandle->isTwoStage();
    bool compact_form = gsHandle->isCompactForm();

    // load local D from handle
    auto viewD  = gsHandle->getD();
    auto viewDa = gsHandle->getDa();

    // load local L from handle
    auto crsmatL      = gsHandle->getL();
    auto values_viewL = crsmatL.values;
    auto rowmap_viewL = crsmatL.graph.row_map;
    auto column_viewL = crsmatL.graph.entries;

    // load local U from handle
    auto crsmatU      = gsHandle->getU();
    auto values_viewU = crsmatU.values;
    auto rowmap_viewU = crsmatU.graph.row_map;
    auto column_viewU = crsmatU.graph.entries;

    // load complement of U+D from handle
    auto crsmatLa      = gsHandle->getLa();
    auto values_viewLa = crsmatLa.values;
    auto rowmap_viewLa = crsmatLa.graph.row_map;

    // load complement of L+D from handle
    auto crsmatUa      = gsHandle->getUa();
    auto values_viewUa = crsmatUa.values;
    auto rowmap_viewUa = crsmatUa.graph.row_map;

    // extract local L, D & U matrices
    using range_policy = Kokkos::RangePolicy<Tag_valuesLU, execution_space>;
    Kokkos::parallel_for(
        "valueLU", range_policy(0, num_rows),
        GS_Functor_t(two_stage, compact_form, diagos_given, num_rows, rowmap_view, column_view, values_view,
                     d_invert_view, rowmap_viewL, values_viewL, viewD, rowmap_viewU, values_viewU, rowmap_viewLa,
                     values_viewLa, viewDa, rowmap_viewUa, values_viewUa));
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << std::endl << "TWO-STAGE GS::NUMERIC::INSERT LU TIME : " << tic << std::endl;
    timer.reset();
#endif

    if (!(gsHandle->isTwoStage())) {
      using namespace KokkosSparse::Experimental;
      auto sptrsv_algo = handle->get_gs_sptrsvL_handle()->get_sptrsv_handle()->get_algorithm();
      if (sptrsv_algo == SPTRSVAlgorithm::SPTRSV_CUSPARSE) {  // symbolic with CuSparse needs
                                                              // values
        // CuSparse needs matrix sorted by column indexes for each row
        // TODO: may need to move this to symbolic/numeric of sptrsv
        KokkosSparse::sort_crs_matrix<execution_space, const_row_map_view_t, entries_view_t, values_view_t>(
            rowmap_viewL, column_viewL, values_viewL);
        KokkosSparse::sort_crs_matrix<execution_space, const_row_map_view_t, entries_view_t, values_view_t>(
            rowmap_viewU, column_viewU, values_viewU);

        // now do symbolic
        sptrsv_symbolic(handle->get_gs_sptrsvL_handle(), rowmap_viewL, crsmatL.graph.entries, values_viewL);
        sptrsv_symbolic(handle->get_gs_sptrsvU_handle(), rowmap_viewU, crsmatU.graph.entries, values_viewU);
      }
    }
  }

  /**
   * Apply solve
   */
  template <typename x_value_array_type, typename y_value_array_type>
  void apply(x_value_array_type localX,  // in/out
             y_value_array_type localB,  // in
             bool init_zero_x_vector = false, int numIter = 1, scalar_t omega = ST::one(), bool apply_forward = true,
             bool apply_backward = true, bool /*update_y_vector*/ = true) {
    const_scalar_t one  = Kokkos::ArithTraits<scalar_t>::one();
    const_scalar_t zero = Kokkos::ArithTraits<scalar_t>::zero();
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    double tic;
    Kokkos::Timer timer;
    Kokkos::fence();
    tic = timer.seconds();
#endif

    //
    auto *gsHandle    = get_gs_handle();
    bool two_stage    = gsHandle->isTwoStage();
    bool compact_form = gsHandle->isCompactForm();
    scalar_t gamma    = gsHandle->getInnerDampFactor();

    GSDirection direction = gsHandle->getSweepDirection();
    if (apply_forward && apply_backward) {
      direction = GS_SYMMETRIC;
    } else if (apply_forward) {
      direction = GS_FORWARD;
    } else if (apply_backward) {
      direction = GS_BACKWARD;
    } else {
      return;
    }

    // load auxiliary matrices from handle
    auto localD   = gsHandle->getD();
    auto crsmatL  = gsHandle->getL();  // lower-part of diagonal block
    auto crsmatU  = gsHandle->getU();  // upper-part of diagonal block
    auto localDa  = gsHandle->getDa();
    auto crsmatLa = gsHandle->getLa();  // complement of L+D (used only for compact form)
    auto crsmatUa = gsHandle->getUa();  // complement of U+D (used only for compact form)

    // wratp A into crsmat
    input_crsmat_t crsmatA("A", num_rows, num_cols, values_view.extent(0), values_view, rowmap_view, column_view);
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    Kokkos::fence();
    tic = timer.seconds();
    std::cout << std::endl
              << "TWO-STAGE GS::APPLY with " << numIter << " outer GS sweeps with omega = " << omega << ", and "
              << gsHandle->getNumInnerSweeps() << " inner JR sweeps, with gamma = " << gamma << " (numRows=" << num_rows
              << ")" << std::endl;
    std::cout << std::endl << "TWO-STAGE GS::APPLY::CREATE CRS_A TIME : " << tic << std::endl;
    timer.reset();
#endif

    // load auxiliary vectors
    int nrhs = localX.extent(1);
    gsHandle->initVectors(num_rows, nrhs);
    auto localR = gsHandle->getVectorR();
    auto localT = gsHandle->getVectorT();
    auto localZ = gsHandle->getVectorZ();

    // outer Gauss-Seidel iteration
    int NumOuterSweeps = gsHandle->getNumOuterSweeps();
    int NumInnerSweeps = gsHandle->getNumInnerSweeps();
    int NumSweeps      = (NumOuterSweeps > numIter ? NumOuterSweeps : numIter);
    if (direction == GS_SYMMETRIC) {
      NumSweeps *= 2;
    }
    if (init_zero_x_vector) {
      KokkosKernels::Impl::zero_vector<x_value_array_type, execution_space>(nrhs, localX);
    }
    for (int sweep = 0; sweep < NumSweeps; ++sweep) {
      bool forward_sweep = (direction == GS_FORWARD || (direction == GS_SYMMETRIC && sweep % 2 == 0));
      // compute residual vector
      KokkosBlas::scal(localR, one, localB);
      if (sweep > 0 || !init_zero_x_vector) {
        if (compact_form) {
          if (forward_sweep) {
            // R = B - U*x
            KokkosSparse::spmv("N", scalar_t(-one), crsmatUa, localX, one, localR);
          } else {
            // R = B - L*x
            KokkosSparse::spmv("N", scalar_t(-one), crsmatLa, localX, one, localR);
          }
          if (omega != one) {
            // R = B - (U + (1-1/omega)D)*x
            scalar_t omega2 = (one / omega - one);
            auto localY     = Kokkos::subview(localX, range_type(0, num_rows), Kokkos::ALL());
            KokkosBlas::mult(zero, localZ, one, localDa, localY);
            KokkosBlas::axpy(omega2, localZ, localR);
          }
        } else {  // not compact_form
          // R = B - A*x
          KokkosSparse::spmv("N", scalar_t(-one), crsmatA, localX, one, localR);
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
          {
            auto localRj = Kokkos::subview(localR, Kokkos::ALL(), range_type(0, 1));
            single_vector_view_t Rj(localRj.data(), num_rows);
            std::cout << "norm(GS)-" << sweep << " " << KokkosBlas::nrm2(Rj) << " ("
                      << (forward_sweep ? "forward" : "backward") << ")" << std::endl;
          }
#endif
        }
      }
      if (!two_stage) {
        // ===== sparse-triangular solve =====
        // TODO: omega is not supported here
        //       (L + D is extracted in initialize_numeric,
        //        but (omega*L + D)^{-1} needs to be applied with omega passed
        //        into apply)
        //       hence, omega = one
        if (omega != one) {
          throw std::invalid_argument(
              " *** TwostageGaussSeidel::apply with omega != one is not "
              "supported with sptrsv ***\n");
        }
        if (forward_sweep) {
          // Z = (omega * L + D)^{-1} * R
          // NOTE: need to go over RHSs
          using namespace KokkosSparse::Experimental;
          for (int j = 0; j < nrhs; j++) {
            auto localRj = Kokkos::subview(localR, Kokkos::ALL(), range_type(j, j + 1));
            auto localZj = Kokkos::subview(localZ, Kokkos::ALL(), range_type(j, j + 1));
            single_vector_view_t Rj(localRj.data(), num_rows);
            single_vector_view_t Zj(localZj.data(), num_rows);
            sptrsv_solve(handle->get_gs_sptrsvL_handle(), crsmatL.graph.row_map, crsmatL.graph.entries, crsmatL.values,
                         Rj, Zj);
          }
        } else {
          using namespace KokkosSparse::Experimental;
          // Z = (omega * U + D)^{-1} * R
          // NOTE: need to go over RHSs
          for (int j = 0; j < nrhs; j++) {
            auto localRj = Kokkos::subview(localR, Kokkos::ALL(), range_type(j, j + 1));
            auto localZj = Kokkos::subview(localZ, Kokkos::ALL(), range_type(j, j + 1));
            single_vector_view_t Rj(localRj.data(), num_rows);
            single_vector_view_t Zj(localZj.data(), num_rows);
            sptrsv_solve(handle->get_gs_sptrsvU_handle(), crsmatU.graph.row_map, crsmatU.graph.entries, crsmatU.values,
                         Rj, Zj);
          }
        }

        // update solution (no omega)
        auto localY = Kokkos::subview(localX, range_type(0, num_rows), Kokkos::ALL());
        if (compact_form) {
          // Y = omega * Z
          KokkosBlas::scal(localY, one, localZ);
        } else {
          // Y = Y + omega * Z
          KokkosBlas::axpy(one, localZ, localY);
        }
      } else {
        // ====== inner Jacobi-Richardson =====
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        // compute initial residual norm
        // > compute RHS for the inner loop, R = B - A*x
        internal_vector_view_t tempR("tempR", num_rows, 1);
        KokkosBlas::scal(tempR, one, localB);
        KokkosSparse::spmv("N", scalar_t(-one), crsmatA, localX, one, tempR);
        // > initial vector for the inner loop is zero
        Kokkos::deep_copy(localZ, zero);
        using Norm_Functor_t = TwostageGaussSeidel_functor<row_map_view_t, entries_view_t, values_view_t>;
        using range_policy   = Kokkos::RangePolicy<Tag_normR, execution_space>;
        {
          mag_t normR = zero;
          Kokkos::parallel_reduce(
              "normR", range_policy(0, num_rows),
              Norm_Functor_t(forward_sweep, num_rows, rowmap_view, column_view, values_view, localD, localZ, tempR),
              normR);
          std::cout << "> norm(JR)-" << 0 << " " << sqrt(normR) << std::endl;
        }
#endif
        // compute starting vector: Z = D^{-1}*R (Z is correction, i.e., output
        // of JR)
        if (NumInnerSweeps == 0) {
          // this is Jacobi-Richardson X_{k+1} := X_{k} + D^{-1}(b-A*X_{k})
          // copy to localZ (output of JR iteration)

          // row-scale: (D^{-1}*L)*Y = D^{-1}*B
          // compute Z := D^{-1}*R
          KokkosBlas::mult(zero, localZ, one, localD, localR);
          // apply inner damping factor, if not one
          if (gamma != one) {
            // Z = gamma * Z
            KokkosBlas::scal(localZ, gamma, localZ);
          }
        } else {
          // copy to localT (workspace used to save D^{-1}*R for JR iteration)
          KokkosBlas::mult(zero, localT, one, localD, localR);
          // initialize Jacobi-Richardson (using R as workspace for JR
          // iteration)
          KokkosBlas::scal(localR, one, localT);

          // apply inner damping factor, if not one
          if (gamma != one) {
            // R = gamma * R
            KokkosBlas::scal(localR, gamma, localR);
          }
        }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        {
          // compute residual norm of the starting vector (D^{-1}R)
          mag_t normR = zero;
          Kokkos::parallel_reduce(
              "normR", range_policy(0, num_rows),
              Norm_Functor_t(forward_sweep, num_rows, rowmap_view, column_view, values_view, localD, localT, tempR),
              normR);
          std::cout << "> norm(JR)-" << 1 << " " << sqrt(normR) << std::endl;
        }
#endif
        // inner Jacobi-Richardson:
        for (int ii = 0; ii < NumInnerSweeps; ii++) {
          // T = D^{-1}*R, and L = D^{-1}*L and U = D^{-1}*U
          // copy T into Z
          KokkosBlas::scal(localZ, one, localT);
          if (forward_sweep) {
            // Z = Z - L*R
            KokkosSparse::spmv("N", scalar_t(-omega), crsmatL, localR, one, localZ);
          } else {
            // Z = R - U*T
            KokkosSparse::spmv("N", scalar_t(-omega), crsmatU, localR, one, localZ);
          }
          // apply inner damping factor, if not one
          if (gamma != one) {
            // Z = gamma * Z
            KokkosBlas::scal(localZ, gamma, localZ);
            // Z = Z + (one - one/gamma) * R
            scalar_t gamma2 = one - gamma;
            KokkosBlas::axpy(gamma2, localR, localZ);
          }
          if (ii + 1 < NumInnerSweeps) {
            // reinitialize (R to be Z)
            KokkosBlas::scal(localR, one, localZ);
          }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
          {
            // compute residual norm(r - (L+D)*y)
            mag_t normR = zero;
            Kokkos::parallel_reduce(
                "normR", range_policy(0, num_rows),
                Norm_Functor_t(forward_sweep, num_rows, rowmap_view, column_view, values_view, localD, localZ, tempR),
                normR);
            std::cout << "> norm(JR)-" << 2 + ii << " " << sqrt(normR) << std::endl;
          }
#endif
        }  // end of inner Jacobi Richardson

        // update solution
        auto localY = Kokkos::subview(localX, range_type(0, num_rows), Kokkos::ALL());
        if (compact_form) {
          // Y := omega * z
          KokkosBlas::scal(localY, omega, localZ);
        } else {
          // Y := X + omega * Z
          KokkosBlas::axpy(omega, localZ, localY);
        }
      }  // end of inner GS sweep
    }    // end of outer GS sweep
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
    {
      // R = B - A*x
      KokkosBlas::scal(localR, one, localB);
      KokkosSparse::spmv("N", scalar_t(-one), crsmatA, localX, one, localR);
      auto localRj = Kokkos::subview(localR, Kokkos::ALL(), range_type(0, 1));
      single_vector_view_t Rj(localRj.data(), num_rows);
      std::cout << "norm(GS)-" << NumSweeps << " " << KokkosBlas::nrm2(Rj) << std::endl;
    }
#endif
  }
};
}  // namespace Impl
}  // namespace KokkosSparse
#endif
