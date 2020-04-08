/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
#include "KokkosKernels_SparseUtils.hpp"

#include "KokkosSparse_gauss_seidel_handle.hpp"

#define KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV

namespace KokkosSparse{
  namespace Impl{

    template <typename HandleType, typename input_row_map_view_t, typename input_entries_view_t, typename input_values_view_t>
    class TwostageGaussSeidel{

    public:
      using TwoStageGaussSeidelHandleType = typename HandleType::TwoStageGaussSeidelHandleType;
      using execution_space = typename HandleType::HandleExecSpace;
      using memory_space    = typename TwoStageGaussSeidelHandleType::memory_space;

      using  const_scalar_t = typename TwoStageGaussSeidelHandleType::const_scalar_t;
      using        scalar_t = typename TwoStageGaussSeidelHandleType::scalar_t;
      using const_ordinal_t = typename TwoStageGaussSeidelHandleType::const_ordinal_t;
      using       ordinal_t = typename TwoStageGaussSeidelHandleType::ordinal_t;
      using       size_type = typename TwoStageGaussSeidelHandleType::size_type;

      using  const_row_map_view_t = typename TwoStageGaussSeidelHandleType::const_row_map_view_t;
      using  row_map_view_t       = typename TwoStageGaussSeidelHandleType::row_map_view_t;
      using  entries_view_t       = typename TwoStageGaussSeidelHandleType::entries_view_t;
      using   values_view_t       = typename TwoStageGaussSeidelHandleType::values_view_t;

      using       crsmat_t      = typename TwoStageGaussSeidelHandleType::crsmat_t;
      using        graph_t      = typename TwoStageGaussSeidelHandleType::graph_t;

      using range_type = Kokkos::pair<int, int>;

      // to wrap input (rowmap, colind, values) into crsmat
      using input_device_t  = Kokkos::Device<typename input_row_map_view_t::execution_space, typename input_row_map_view_t::memory_space>;
      using input_memory_t  = typename input_values_view_t::memory_traits;
      using input_scalar_t  = typename input_values_view_t::value_type;
      using input_ordinal_t = typename input_entries_view_t::value_type;
      using input_size_t    = typename input_row_map_view_t::value_type;
      using input_crsmat_t  = KokkosSparse::CrsMatrix <input_scalar_t,
                                                       input_ordinal_t,
                                                       input_device_t,
                                                       input_memory_t,
                                                       input_size_t>;
      using input_graph_t  = typename input_crsmat_t::StaticCrsGraphType;
      using single_vector_view_t = Kokkos::View<scalar_t*, Kokkos::LayoutLeft, input_device_t, input_memory_t>;

    private:
      HandleType *handle;

      //Get the specialized TwostageGaussSeidel handle from the main handle
      TwoStageGaussSeidelHandleType* get_gs_handle()
      {
        auto gsHandle = dynamic_cast<TwoStageGaussSeidelHandleType*>(this->handle->get_gs_handle());
        if(!gsHandle)
        {
          throw std::runtime_error("TwostageGaussSeidel: GS handle has not been created, or is set up for Cluster GS.");
        }
        return gsHandle;
      }

      bool diagos_given;
      const_ordinal_t num_rows, num_cols;

      input_row_map_view_t rowmap_view;
      input_entries_view_t column_view;
      input_values_view_t  values_view;
      input_values_view_t  d_invert_view;

    // --------------------------------------------------------- //
    public:
      // tag for counting nnz
      struct Tag_countNnzL{};
      struct Tag_countNnzU{};
      // tag for inserting entries
      struct Tag_entriesL{};
      struct Tag_entriesU{};
      struct Tag_entriesLU{};
      // tag for inserting values
      struct Tag_valuesL{};
      struct Tag_valuesU{};
      struct Tag_valuesLU{};

      template <typename output_row_map_view_t, 
                typename output_entries_view_t, 
                typename output_values_view_t>
      struct TwostageGaussSeidel_functor {

        public:
        // input
        bool two_stage;
        bool diagos_given;
        const_ordinal_t num_rows;
        input_row_map_view_t rowmap_view;
        input_entries_view_t column_view;
        input_values_view_t  values_view;
        input_values_view_t  d_invert_view;
        // output
        output_row_map_view_t  row_map;
        output_entries_view_t  entries;
        output_values_view_t   values;
        output_values_view_t   diags;
        // output
        output_row_map_view_t  row_map2;
        output_entries_view_t  entries2;
        output_values_view_t   values2;

        // ------------------------------------------------------- //
        // Constructors
        // for counting nnz
        TwostageGaussSeidel_functor (
                  bool two_stage_,
                  const_ordinal_t num_rows_,
                  input_row_map_view_t  rowmap_view_,
                  input_entries_view_t  column_view_,
                  output_row_map_view_t row_map_) :
          two_stage(two_stage_),
          diagos_given(false),
          num_rows(num_rows_),
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(),
          d_invert_view(),
          row_map(row_map_)
        {}

        // for storing entries
        TwostageGaussSeidel_functor (
                  bool two_stage_,
                  const_ordinal_t num_rows_,
                  input_row_map_view_t  rowmap_view_,
                  input_entries_view_t  column_view_,
                  output_row_map_view_t row_map_,
                  output_entries_view_t entries_) :
          two_stage(two_stage_),
          diagos_given(false),
          num_rows(num_rows_),
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(),
          d_invert_view(),
          row_map(row_map_),
          entries(entries_)
        {}

        // for storing values
        TwostageGaussSeidel_functor (
                  bool two_stage_,
                  const_ordinal_t num_rows_,
                  input_row_map_view_t  rowmap_view_,
                  input_entries_view_t  column_view_,
                  input_values_view_t   values_view_,
                  output_row_map_view_t row_map_,
                  output_values_view_t  values_,
                  output_values_view_t  diags_) :
          two_stage(two_stage_),
          diagos_given(false),
          num_rows(num_rows_),
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(values_view_),
          d_invert_view(),
          row_map(row_map_),
          entries(),
          values(values_),
          diags(diags_)
        {}

        // for storing booth L&U entries
        TwostageGaussSeidel_functor (
                  bool two_stage_,
                  const_ordinal_t num_rows_,
                  input_row_map_view_t  rowmap_view_,
                  input_entries_view_t  column_view_,
                  output_row_map_view_t row_map_,
                  output_entries_view_t entries_,
                  output_row_map_view_t row_map2_,
                  output_entries_view_t entries2_) :
          two_stage(two_stage_),
          diagos_given(false),
          num_rows(num_rows_),
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(),
          d_invert_view(),
          row_map(row_map_),
          entries(entries_),
          values(),
          diags(),
          row_map2(row_map2_),
          entries2(entries2_)
        {}

        // for storing both L&U values (with D extracted)
        TwostageGaussSeidel_functor (
                  bool two_stage_,
                  bool diagos_given_,
                  const_ordinal_t num_rows_,
                  input_row_map_view_t  rowmap_view_,
                  input_entries_view_t  column_view_,
                  input_values_view_t   values_view_,
                  input_values_view_t   d_invert_view_,
                  output_row_map_view_t row_map_,
                  output_entries_view_t entries_,
                  output_values_view_t  values_,
                  output_values_view_t  diags_,
                  output_row_map_view_t row_map2_,
                  output_entries_view_t entries2_,
                  output_values_view_t  values2_) :
          two_stage(two_stage_),
          diagos_given(diagos_given_),
          num_rows(num_rows_),
          rowmap_view(rowmap_view_),
          column_view(column_view_),
          values_view(values_view_),
          d_invert_view(d_invert_view_),
          row_map(row_map_),
          entries(entries_),
          values(values_),
          diags(diags_),
          row_map2(row_map2_),
          entries2(entries2_),
          values2(values2_)
        {}


        // ------------------------------------------------------- //
        // functor for counting nnzL (with parallel_reduce)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_countNnzL&, const ordinal_t i, ordinal_t &nnz) const
        {
          ordinal_t nnz_i = 0;
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) < i) {
              nnz_i ++;
            } else if(!two_stage && column_view (k) == i) {
              nnz_i ++;
            }
          }
          row_map (i+1) = nnz_i;
          if (i == 0) {
            row_map (0) = 0;
          }
          nnz +=  nnz_i;
        }

        // functor for storing entriesL (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_entriesL&, const ordinal_t i) const
        {
          ordinal_t nnz = row_map (i);
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) < i) {
              entries (nnz) = column_view (k);
              nnz ++;
            } else if(!two_stage && column_view (k) == i) {
              entries (nnz) = column_view (k);
              nnz ++;
            }
          }
        }

        // functor for storing valuesL (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_valuesL&, const ordinal_t i) const
        {
          const_scalar_t one = Kokkos::Details::ArithTraits<scalar_t>::one ();
          ordinal_t nnz = row_map (i);
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) < i) {
              values (nnz) = values_view (k);
              nnz ++;
            } else if (column_view (k) == i) {
              if (two_stage) {
                if (diagos_given) {
                  diags (i) = d_invert_view (i);
                } else {
                  diags (i) = one / values_view (k);
                }
              } else {
                values (nnz) = values_view (k);
                nnz ++;
              }
            }
          }
          #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
          if (two_stage) {
            for (size_type k = row_map (i); k < nnz; k++) {
              values (k) *= diags (i);
            }
          }
          #endif
        }


        // ------------------------------------------------------- //
        // functor for counting nnzU (with parallel_reduce)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_countNnzU&, const ordinal_t i, ordinal_t &nnz) const
        {
          ordinal_t nnz_i = 0;
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) > i && column_view (k) < num_rows) {
              nnz_i ++;
            } else if(!two_stage && column_view (k) == i) {
              nnz_i ++;
            }
          }
          row_map (i+1) = nnz_i;
          if (i == 0) {
            row_map (0) = 0;
          }
          nnz +=  nnz_i;
        }

        // functor for storing entriesU (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_entriesU&, const ordinal_t i) const
        {
          ordinal_t nnz = row_map (i);
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) > i && column_view (k) < num_rows) {
              entries (nnz) = column_view (k);
              nnz ++;
            } else if(!two_stage && column_view (k) == i) {
              entries (nnz) = column_view (k);
              nnz ++;
            }
          }
        }

        // functor for storing valuesU (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_valuesU&, const ordinal_t i) const
        {
          const_scalar_t one = Kokkos::Details::ArithTraits<scalar_t>::one ();
          ordinal_t nnz = row_map (i);
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) == i) {
              if (two_stage) {
                if (diagos_given) {
                  diags (i) = d_invert_view (i);
                } else {
                  diags (i) = one / values_view (k);
                }
              } else {
                values (nnz) = values_view (k);
                nnz ++;
              }
            } else if (column_view (k) > i && column_view (k) < num_rows) {
              values (nnz) = values_view (k);
              nnz ++;
            }
          }
          #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
          if (two_stage) {
            for (size_type k = row_map (i); k < nnz; k++) {
              values (k) *= diags (i);
            }
          }
          #endif
        }

        // ------------------------------------------------------- //
        // functor for storing entriesL and entriesU (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_entriesLU&, const ordinal_t i) const
        {
          ordinal_t nnzL = row_map (i);
          ordinal_t nnzU = row_map2 (i);
          if (!two_stage) {
            // NOTE: Kokkos' sptrsv assumes diagonal of U to be at the start
            entries2 (nnzU) = i;
            nnzU ++;
          }
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) < i) {
              entries (nnzL) = column_view (k);
              nnzL ++;
            } else if (column_view (k) > i && column_view (k) < num_rows) {
              entries2 (nnzU) = column_view (k);
              nnzU ++;
            }
          }
          if (!two_stage) {
            // NOTE: Kokkos' sptrsv assumes diagonal of L to be at the end
            entries (nnzL) = i;
            nnzL ++;
          }
        }

        // functor for storing both valuesL & valuesU (with parallel_for)
        KOKKOS_INLINE_FUNCTION
        void operator()(const Tag_valuesLU&, const ordinal_t i) const
        {
          const_scalar_t one = Kokkos::Details::ArithTraits<scalar_t>::one ();
          ordinal_t nnzL = row_map (i);
          ordinal_t nnzU = row_map2 (i);
          if (!two_stage) {
            // Kokkos' sptrsv assumes diagonal U to come at the start, so increment nnzU
            nnzU ++;
          }
          for (size_type k = rowmap_view (i); k < rowmap_view (i+1); k++) {
            if (column_view (k) < i) {
              // save L (without diag)
              values (nnzL) = values_view (k);
              nnzL ++;
            } else if (column_view (k) == i) {
              // save D
              if (diagos_given) {
                // as inverse
                diags (i) = d_invert_view (i);
              } else {
                // as original
                diags (i) = values_view (k);
              }
            } else if (column_view (k) < num_rows) {
              // save U (without diag)
              values2 (nnzU) = values_view (k);
              nnzU ++;
            }
          }
          if (!two_stage) {
            // if using sptrsv, add diagonals in L and U
            // > Kokkos' sptrsv assumes diagonal of L and U to come at end and start
            nnzU = row_map2 (i);
            if (diagos_given) {
              values2 (nnzU) = one / diags (i);
              values (nnzL) = one / diags (i);
            } else {
              values2 (nnzU) = diags (i);
              values (nnzL) = diags (i);
            }
          }
          #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
          if (two_stage) {
            if (!diagos_given) {
              // when diag is provided, it is already provided as inverse
              diags (i) = one / diags (i);
            }
            // compute inv(D)*L
            for (size_type k = row_map (i); k < row_map (i+1); k++) {
              values (k) *= diags (i);
            }
            for (size_type k = row_map2 (i); k < row_map2 (i+1); k++) {
              values2 (k) *= diags (i);
            }
          }
          #endif
        }
      };
    // --------------------------------------------------------- //


    public:
      /**
       * \brief constructor
       */
      // for symbolic (wihout values)
      TwostageGaussSeidel(HandleType *handle_,
                  const_ordinal_t num_rows_,
                  const_ordinal_t num_cols_,
                  input_row_map_view_t rowmap_view_,
                  input_entries_view_t column_view_) :
        handle(handle_),
        diagos_given(false),
        num_rows(num_rows_), num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(),
        d_invert_view() {}

      // for numeric/solve (with values)
      TwostageGaussSeidel (HandleType *handle_,
                           const_ordinal_t num_rows_,
                           const_ordinal_t num_cols_,
                           input_row_map_view_t rowmap_view_,
                           input_entries_view_t column_view_,
                           input_values_view_t values_view_) :
        handle(handle_),
        diagos_given(false),
        num_rows(num_rows_), num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(values_view_),
        d_invert_view() {}

      // for numeric/solve (with values and diagonal)
      TwostageGaussSeidel (HandleType *handle_,
                           const_ordinal_t num_rows_,
                           const_ordinal_t num_cols_,
                           input_row_map_view_t rowmap_view_,
                           input_entries_view_t column_view_,
                           input_values_view_t values_view_,
                           input_values_view_t d_invert_view_) :
        handle(handle_),
        diagos_given(true),
        num_rows(num_rows_), num_cols(num_cols_),
        rowmap_view(rowmap_view_),
        column_view(column_view_),
        values_view(values_view_),
        d_invert_view(d_invert_view_) {}


      /**
       * Symbolic setup
       */
      void initialize_symbolic ()
      {
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        double tic;
        Kokkos::Impl::Timer timer;
        Kokkos::fence();
        tic = timer.seconds ();
#endif
        auto *gsHandle = get_gs_handle();
        bool two_stage = gsHandle->isTwoStage ();
        GSDirection direction = gsHandle->getSweepDirection ();
        using GS_Functor_t = TwostageGaussSeidel_functor<row_map_view_t, entries_view_t, values_view_t>;
        // count nnz in local L & U matrices (rowmap_viewL/rowmap_viewU stores offsets for each row)
        ordinal_t nnzL = 0;
        ordinal_t nnzU = 0;
        row_map_view_t  rowmap_viewL ("row_mapL", num_rows+1);
        row_map_view_t  rowmap_viewU ("row_mapU", num_rows+1);
        if (direction == GS_FORWARD || direction == GS_SYMMETRIC) {
          using range_policy = Kokkos::RangePolicy <Tag_countNnzL, execution_space>;
          Kokkos::parallel_reduce ("nnzL", range_policy (0, num_rows),
                                   GS_Functor_t (two_stage, num_rows, rowmap_view, column_view,
                                                                      rowmap_viewL),
                                   nnzL);
        }
        if (direction == GS_BACKWARD || direction == GS_SYMMETRIC) {
          using range_policy = Kokkos::RangePolicy <Tag_countNnzU, execution_space>;
          Kokkos::parallel_reduce ("nnzU", range_policy (0, num_rows),
                                   GS_Functor_t (two_stage, num_rows, rowmap_view, column_view,
                                                                      rowmap_viewU),
                                   nnzU);
        }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << std::endl << "TWO-STAGE GS::SYMBOLIC::COUNT-NNZ TIME : " << tic << std::endl;
        timer.reset();
#endif
        // shift ptr so that it now contains offsets (combine it with the previous functor calls?)
        if (direction == GS_FORWARD || direction == GS_SYMMETRIC) {
          KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<row_map_view_t, execution_space> 
            (1+num_rows, rowmap_viewL);
        }
        if (direction == GS_BACKWARD || direction == GS_SYMMETRIC) {
          KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<row_map_view_t, execution_space> 
            (1+num_rows, rowmap_viewU);
        }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << "TWO-STAGE GS::SYMBOLIC::COMP-PTR TIME  : " << tic << std::endl;
        timer.reset();
#endif
        // allocate memory to store local D
        values_view_t viewD (Kokkos::ViewAllocateWithoutInitializing("diags"), num_rows);

        // allocate memory to store local L
        entries_view_t  column_viewL (Kokkos::ViewAllocateWithoutInitializing("entriesL"), nnzL);
        values_view_t   values_viewL (Kokkos::ViewAllocateWithoutInitializing("valuesL"),  nnzL);

        // allocate memory to store local U
        entries_view_t  column_viewU (Kokkos::ViewAllocateWithoutInitializing("entriesU"), nnzU);
        values_view_t   values_viewU (Kokkos::ViewAllocateWithoutInitializing("valuesU"),  nnzU);
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << "TWO-STAGE GS::SYMBOLIC::ALLOCATE TIME  : " << tic << std::endl;
        timer.reset();
#endif

        {
          // extract local L & U structures (for computing (L+D)^{-1} or (D+U)^{-1})
          using range_policy = Kokkos::RangePolicy <Tag_entriesLU, execution_space>;
          Kokkos::parallel_for ("entryLU", range_policy (0, num_rows),
                                GS_Functor_t (two_stage, num_rows, rowmap_view, column_view,
                                                                   rowmap_viewL, column_viewL,
                                                                   rowmap_viewU, column_viewU));
        }
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << "TWO-STAGE GS::SYMBOLIC::INSERT TIME    : " << tic << std::endl;
        timer.reset();
#endif

        // construct CrsMat with them
        graph_t graphL (column_viewL, rowmap_viewL);
        graph_t graphU (column_viewU, rowmap_viewU);
        crsmat_t crsmatL ("L", num_rows, values_viewL, graphL);
        crsmat_t crsmatU ("U", num_rows, values_viewU, graphU);

        // store them in handle
        gsHandle->setL (crsmatL);
        gsHandle->setU (crsmatU);
        gsHandle->setD (viewD);
        if (!(gsHandle->isTwoStage ())) {
          // create SpTRSV handles for classical GS
          using namespace KokkosSparse::Experimental;
          auto sptrsv_algo = handle->get_gs_sptrsvL_handle()->get_sptrsv_handle()->get_algorithm();
          if (sptrsv_algo != SPTRSVAlgorithm::SPTRSV_CUSPARSE) { // symbolic with CuSparse needs values
            sptrsv_symbolic (handle->get_gs_sptrsvL_handle(), rowmap_viewL, crsmatL.graph.entries);
            sptrsv_symbolic (handle->get_gs_sptrsvU_handle(), rowmap_viewU, crsmatU.graph.entries);
          }
        }
      }


      /**
       * Numerical setup
       */
      void initialize_numeric ()
      {
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        double tic;
        Kokkos::Impl::Timer timer;
        Kokkos::fence();
        timer.reset();
#endif
        using GS_Functor_t = TwostageGaussSeidel_functor<const_row_map_view_t, entries_view_t, values_view_t>;

        auto *gsHandle = get_gs_handle();
        bool two_stage = gsHandle->isTwoStage ();

        // load local D from handle
        auto viewD = gsHandle->getD ();

        // load local L from handle
        auto crsmatL = gsHandle->getL ();
        auto values_viewL = crsmatL.values;
        auto rowmap_viewL = crsmatL.graph.row_map;
        auto column_viewL = crsmatL.graph.entries;

        // load local U from handle
        auto crsmatU = gsHandle->getU ();
        auto values_viewU = crsmatU.values;
        auto rowmap_viewU = crsmatU.graph.row_map;
        auto column_viewU = crsmatU.graph.entries;

        // extract local L, D & U matrices
        using range_policy = Kokkos::RangePolicy <Tag_valuesLU, execution_space>;
        Kokkos::parallel_for ("valueLU", range_policy (0, num_rows),
                              GS_Functor_t (two_stage, diagos_given, num_rows,
                                            rowmap_view, column_view, values_view, d_invert_view,
                                            rowmap_viewL, column_viewL, values_viewL, viewD,
                                            rowmap_viewU, column_viewU, values_viewU));
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << std::endl << "TWO-STAGE GS::NUMERIC::INSERT LU TIME : " << tic << std::endl;
        timer.reset();
#endif

        if (!(gsHandle->isTwoStage ())) {
          using namespace KokkosSparse::Experimental;
          auto sptrsv_algo = handle->get_gs_sptrsvL_handle()->get_sptrsv_handle()->get_algorithm();
          if (sptrsv_algo == SPTRSVAlgorithm::SPTRSV_CUSPARSE) { // symbolic with CuSparse needs values
            // CuSparse needs matrix sorted by column indexes for each row
            // TODO: may need to move this to symbolic/numeric of sptrsv
            KokkosKernels::Impl::sort_crs_matrix <execution_space, const_row_map_view_t, entries_view_t, values_view_t>
              (rowmap_viewL, column_viewL, values_viewL);
            KokkosKernels::Impl::sort_crs_matrix <execution_space, const_row_map_view_t, entries_view_t, values_view_t>
              (rowmap_viewU, column_viewU, values_viewU);

            // now do symbolic
            sptrsv_symbolic (handle->get_gs_sptrsvL_handle(), rowmap_viewL, crsmatL.graph.entries, values_viewL);
            sptrsv_symbolic (handle->get_gs_sptrsvU_handle(), rowmap_viewU, crsmatU.graph.entries, values_viewU);
          }
        }
      }


      /**
       * Apply solve
       */
      template <typename x_value_array_type, typename y_value_array_type>
      void apply (x_value_array_type localX, // in/out
                  y_value_array_type localB, // in
                  bool init_zero_x_vector = false,
                  int numIter = 1,
                  scalar_t omega = Kokkos::Details::ArithTraits<scalar_t>::one(),
                  bool apply_forward = true,
                  bool apply_backward = true,
                  bool update_y_vector = true)
      {
        const_scalar_t one = Kokkos::Details::ArithTraits<scalar_t>::one ();
        const_scalar_t zero = Kokkos::Details::ArithTraits<scalar_t>::zero ();
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        double tic;
        Kokkos::Impl::Timer timer;
        Kokkos::fence();
        tic = timer.seconds ();
#endif

        //
        auto *gsHandle = get_gs_handle();
        bool two_stage = gsHandle->isTwoStage ();
        GSDirection direction = gsHandle->getSweepDirection ();
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
        auto localD = gsHandle->getD ();
        auto crsmatL = gsHandle->getL ();
        auto crsmatU = gsHandle->getU ();

        // wratp A into crsmat
        input_graph_t graphA (column_view, rowmap_view);
        input_crsmat_t crsmatA ("A", num_rows, values_view, graphA);
#ifdef KOKKOSSPARSE_IMPL_TIME_TWOSTAGE_GS
        Kokkos::fence();
        tic = timer.seconds ();
        std::cout << std::endl << "TWO-STAGE GS::APPLY::CREATE CRS_A TIME : " << tic << std::endl;
        timer.reset();
#endif

        // load auxiliary vectors
        int nrows = num_rows;
        int nrhs = localX.extent (1);
        gsHandle->initVectors (nrows, nrhs);
        auto localR = gsHandle->getVectorR ();
        auto localT = gsHandle->getVectorT ();
        auto localZ = gsHandle->getVectorZ ();

        // outer Gauss-Seidel iteration
        int NumSweeps = numIter;
        int NumInnerSweeps = gsHandle->getNumInnerSweeps ();
        if (direction == GS_SYMMETRIC) {
          NumSweeps *= 2;
        }
        if (init_zero_x_vector) {
          KokkosKernels::Impl::zero_vector<x_value_array_type, execution_space>(nrhs, localX);
        }
        for (int sweep = 0; sweep < NumSweeps; ++sweep) {
          // R = B - A*x
          KokkosBlas::scal (localR, one, localB);
          if (sweep > 0 || !init_zero_x_vector) {
            KokkosSparse::
            spmv ("N", scalar_t(-one), crsmatA,
                                       localX,
                                 one,  localR);
          }
          if (!two_stage) { // ===== sparse-triangular solve =====
            if (direction == GS_FORWARD ||
               (direction == GS_SYMMETRIC && sweep%2 == 0)) {
              // Z = (L+D)^{-1} * R
              // NOTE: need to go over RHSs
              using namespace KokkosSparse::Experimental;
              for (int j = 0; j < nrhs; j++) {
                auto localRj = Kokkos::subview (localR, Kokkos::ALL (), range_type (j, j+1));
                auto localZj = Kokkos::subview (localZ, Kokkos::ALL (), range_type (j, j+1));
                single_vector_view_t Rj (localRj.data (), nrows);
                single_vector_view_t Zj (localZj.data (), nrows);
                sptrsv_solve (handle->get_gs_sptrsvL_handle(), crsmatL.graph.row_map, crsmatL.graph.entries, crsmatL.values, Rj, Zj);
              }
            } else {
              using namespace KokkosSparse::Experimental;
              // Z = (U+D)^{-1} * R
              // NOTE: need to go over RHSs
              for (int j = 0; j < nrhs; j++) {
                auto localRj = Kokkos::subview (localR, Kokkos::ALL (), range_type (j, j+1));
                auto localZj = Kokkos::subview (localZ, Kokkos::ALL (), range_type (j, j+1));
                single_vector_view_t Rj (localRj.data (), nrows);
                single_vector_view_t Zj (localZj.data (), nrows);
                sptrsv_solve (handle->get_gs_sptrsvU_handle(), crsmatU.graph.row_map, crsmatU.graph.entries, crsmatU.values, Rj, Zj);
              }
            }
          } else { // ====== inner Jacobi-Richardson =====
            // compute starting vector: Z = D^{-1}*R (Z is correction, i.e., output of JR)
            #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
            if (NumInnerSweeps == 0) {
              // this is Jacobi-Richardson X_{k+1} := X_{k} + D^{-1}(b-A*X_{k})
              // copy to localZ (output of JR iteration)
              KokkosBlas::mult (zero, localZ,
                                one,  localD, localR);
            } else {
              // copy to localT (workspace used to save D^{-1}*R for JR iteration)
              KokkosBlas::mult (zero, localT,
                                one,  localD, localR);
              // initialize Jacobi-Richardson (using R as workspace for JR iteration)
              KokkosBlas::scal (localR, one, localT);
            }
            #else
            KokkosBlas::mult (zero, localT,
                              one,  localD, localR);
            #endif
            // inner Jacobi-Richardson:
            for (int ii = 0; ii < NumInnerSweeps; ii++) {
              #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
              // T = D^{-1}*R, and L = D^{-1}*L and U = D^{-1}*U
              // copy T into Z
              KokkosBlas::scal (localZ, one, localT);
              #else
              // Z = R
              KokkosBlas::scal (localZ, one, localR);
              #endif
              if (direction == GS_FORWARD ||
                 (direction == GS_SYMMETRIC && sweep%2 == 0)) {
                // Z = Z - L*R
                KokkosSparse::
                spmv("N", scalar_t(-one), crsmatL,
                                          localR,
                                    one, localZ);
              }
              else {
                // Z = R - U*T
                KokkosSparse::
                spmv("N", scalar_t(-one), crsmatU,
                                          localR,
                                    one, localZ);
              }
              #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
              if (ii+1 < NumInnerSweeps) {
                // reinitialize (R to be Z)
                KokkosBlas::scal (localR, one, localZ);
              }
              #else
              // T = D^{-1}*Z
              KokkosBlas::mult (zero, localT,
                                one,  localD, localZ);
              #endif
            } // end of inner Jacobi Richardson
          }
          // Y = X + T
          auto localY = Kokkos::subview (localX, range_type(0, nrows), Kokkos::ALL ());
          #if defined(KOKKOSSPARSE_IMPL_TWOSTAGE_GS_MERGE_SPMV)
          KokkosBlas::axpy (one, localZ, localY);
          #else
          KokkosBlas::axpy (one, localT, localY);
          #endif
        } // end of outer GS sweep
      }
    };
  }
}
#endif
