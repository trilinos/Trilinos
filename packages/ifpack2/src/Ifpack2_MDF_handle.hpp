/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_RILUK_decl.hpp
/// \brief Declaration of MDF interface

#ifndef IFPACK2_MDF_HANDLE_HPP
#define IFPACK2_MDF_HANDLE_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <iostream>

#include "Tpetra_Core.hpp"
// #include "Tpetra_MatrixIO.hpp"
// #include "MatrixMarket_Tpetra.hpp"
// #include "TpetraExt_MatrixMatrix.hpp"

#include "KokkosKernels_Sorting.hpp"

// #include "Ifpack2_UnitTestHelpers.hpp"
// #include "Ifpack2_MDF.hpp"

namespace Ifpack2 {
namespace MDFImpl {


template <class crs_matrix_type>
struct MDF_count_lower {

  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using size_type    = typename crs_matrix_type::ordinal_type;
  using value_type   = typename crs_matrix_type::size_type;

  crs_matrix_type A;
  col_ind_type permutation;
  col_ind_type permutation_inv;

  MDF_count_lower(crs_matrix_type A_, col_ind_type permutation_, col_ind_type permutation_inv_) :
    A(A_), permutation(permutation_), permutation_inv(permutation_inv_) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type rowIdx, value_type& update) const {
    permutation(rowIdx) = rowIdx;
    permutation_inv(rowIdx) = rowIdx;
    for(value_type entryIdx = A.graph.row_map(rowIdx); entryIdx < A.graph.row_map(rowIdx + 1); ++entryIdx) {
      if(A.graph.entries(entryIdx) <= rowIdx) {
        update += 1;
      }
    }
  }

}; // MDF_count_lower

template <class crs_matrix_type>
struct MDF_discarded_fill_norm {

  using static_crs_graph_type = typename crs_matrix_type::StaticCrsGraphType;
  using col_ind_type          = typename static_crs_graph_type::entries_type::non_const_type;
  using values_type           = typename crs_matrix_type::values_type::non_const_type;
  using size_type             = typename crs_matrix_type::size_type;
  using ordinal_type          = typename crs_matrix_type::ordinal_type;
  using scalar_type           = typename crs_matrix_type::value_type;
  using KAS                   = typename Kokkos::ArithTraits<scalar_type>;

  const scalar_type zero = KAS::zero();

  crs_matrix_type A, At;
  ordinal_type factorization_step;
  col_ind_type permutation;

  values_type  discarded_fill;
  col_ind_type deficiency;
  int verbosity;

  MDF_discarded_fill_norm(crs_matrix_type A_, crs_matrix_type At_, ordinal_type factorization_step_,
                          col_ind_type permutation_,
                          values_type  discarded_fill_, col_ind_type deficiency_, int verbosity_)
    : A(A_), At(At_), factorization_step(factorization_step_), permutation(permutation_),
      discarded_fill(discarded_fill_), deficiency(deficiency_), verbosity(verbosity_) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type i) const {
    ordinal_type rowIdx = permutation(i);
    scalar_type discard_norm = zero, diag_val = zero;
    bool entryIsDiscarded = true;
    ordinal_type numFillEntries = 0;
    for(size_type alphaIdx = At.graph.row_map(rowIdx); alphaIdx < At.graph.row_map(rowIdx + 1); ++alphaIdx) {
      ordinal_type fillRowIdx = At.graph.entries(alphaIdx);
      bool row_not_eliminated = true;
      for(ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
        if(fillRowIdx == permutation(stepIdx)) {
          row_not_eliminated = false;
        }
      }

      if(fillRowIdx != rowIdx && row_not_eliminated) {
        for(size_type betaIdx = A.graph.row_map(rowIdx); betaIdx < A.graph.row_map(rowIdx + 1); ++betaIdx) {
          ordinal_type fillColIdx = A.graph.entries(betaIdx);
          bool col_not_eliminated = true;
          for(ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
            if(fillColIdx == permutation(stepIdx)) {
              col_not_eliminated = false;
            }
          }

          if(fillColIdx != rowIdx && col_not_eliminated) {
            entryIsDiscarded = true;
            for(size_type entryIdx = A.graph.row_map(fillRowIdx); entryIdx < A.graph.row_map(fillRowIdx + 1); ++entryIdx) {
              if(A.graph.entries(entryIdx) == fillColIdx) {entryIsDiscarded = false;}
            }
            if(entryIsDiscarded) {
              numFillEntries += 1;
              discard_norm += KAS::abs(At.values(alphaIdx)*A.values(betaIdx))*KAS::abs(At.values(alphaIdx)*A.values(betaIdx));
              if(verbosity > 1) {
                printf("Adding value A[%d,%d]=%f to discard norm of row %d\n",
                       int(At.graph.entries(alphaIdx)), int(A.graph.entries(betaIdx)),
                       KAS::abs(At.values(alphaIdx)*A.values(betaIdx))*KAS::abs(At.values(alphaIdx)*A.values(betaIdx)),
                       int(rowIdx));
              }
            }
          }
        }
      } else if(fillRowIdx == rowIdx) {
        diag_val = At.values(alphaIdx);
        if(verbosity > 1) {
          printf("Row %d diagonal value dected, values(%d)=%f\n", int(rowIdx), int(alphaIdx),
                 At.values(alphaIdx));
        }
      }
    }

    // TODO add a check on `diag_val == zero`
    discard_norm              = discard_norm / (diag_val*diag_val);
    discarded_fill(rowIdx)    = discard_norm;
    deficiency(rowIdx)        = numFillEntries;
    if(verbosity > 0) {
      const ordinal_type degree = ordinal_type(A.graph.row_map(rowIdx + 1) - A.graph.row_map(rowIdx) - 1);
      printf("Row %d has discarded fill of %f, deficiency of %d and degree %d\n", rowIdx,
             KAS::sqrt(discard_norm), deficiency(rowIdx), degree);
    }
  }

}; // MDF_discarded_fill_norm

template <class crs_matrix_type>
struct MDF_select_row{

  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using row_map_type = typename crs_matrix_type::StaticCrsGraphType::row_map_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using scalar_type  = typename crs_matrix_type::value_type;

  // type used to perform the reduction
  // do not confuse it with scalar_type!
  using value_type = typename crs_matrix_type::ordinal_type;

  value_type   factorization_step;
  values_type  discarded_fill;
  col_ind_type deficiency;
  row_map_type row_map;
  col_ind_type permutation;

  MDF_select_row(value_type factorization_step_, values_type  discarded_fill_,
                 col_ind_type deficiency_, row_map_type row_map_, col_ind_type permutation_)
    : factorization_step(factorization_step_), discarded_fill(discarded_fill_),
      deficiency(deficiency_), row_map(row_map_), permutation(permutation_) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type src, ordinal_type& dst) const {
    const ordinal_type src_perm = permutation(src);
    const ordinal_type dst_perm = permutation(dst);
    const ordinal_type degree_src = row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst = row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if(discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) == deficiency(dst_perm))
       && (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) == deficiency(dst_perm))
       && (degree_src == degree_dst)
       && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dst,
            const volatile value_type& src) const {
    const ordinal_type src_perm = permutation(src);
    const ordinal_type dst_perm = permutation(dst);
    const ordinal_type degree_src = row_map(src_perm + 1) - row_map(src_perm) - 1;
    const ordinal_type degree_dst = row_map(dst_perm + 1) - row_map(dst_perm) - 1;

    if(discarded_fill(src_perm) < discarded_fill(dst_perm)) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) < deficiency(dst_perm))) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) == deficiency(dst_perm))
       && (degree_src < degree_dst)) {
      dst = src;
      return;
    }

    if((discarded_fill(src_perm) == discarded_fill(dst_perm))
       && (deficiency(src_perm) == deficiency(dst_perm))
       && (degree_src == degree_dst)
       && (src_perm < dst_perm)) {
      dst = src;
      return;
    }

    return;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& dst) const {
    dst = Kokkos::ArithTraits<ordinal_type>::zero();
  }

}; // MDF_select_row

template<class crs_matrix_type>
struct MDF_factorize_row{

  using row_map_type = typename crs_matrix_type::StaticCrsGraphType::row_map_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using size_type    = typename crs_matrix_type::size_type;
  using value_type   = typename crs_matrix_type::value_type;

  crs_matrix_type A, At;

  row_map_type row_mapL;
  col_ind_type entriesL;
  values_type  valuesL;

  row_map_type row_mapU;
  col_ind_type entriesU;
  values_type  valuesU;

  col_ind_type permutation, permutation_inv;
  ordinal_type selected_row_idx, factorization_step;

  int verbosity;

  MDF_factorize_row(crs_matrix_type A_, crs_matrix_type At_,
                    row_map_type row_mapL_, col_ind_type entriesL_, values_type valuesL_,
                    row_map_type row_mapU_, col_ind_type entriesU_, values_type valuesU_,
                    col_ind_type permutation_, col_ind_type permutation_inv_,
                    ordinal_type selected_row_idx_, ordinal_type factorization_step_,
                    int verbosity_)
    : A(A_), At(At_),
      row_mapL(row_mapL_), entriesL(entriesL_), valuesL(valuesL_),
    row_mapU(row_mapU_), entriesU(entriesU_), valuesU(valuesU_),
    permutation(permutation_), permutation_inv(permutation_inv_),
    selected_row_idx(selected_row_idx_), factorization_step(factorization_step_),
    verbosity(verbosity_) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type /* idx */) const {
    const ordinal_type selected_row = permutation(selected_row_idx);

    // Swap entries in permutation vectors
    permutation(selected_row_idx) = permutation(factorization_step);
    permutation(factorization_step) = selected_row;
    permutation_inv(permutation(factorization_step)) = factorization_step;
    permutation_inv(permutation(selected_row_idx)) = selected_row_idx;

    if(verbosity > 0) {
      printf("Permutation vector: { ");
      for(ordinal_type rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        printf("%d ", permutation(rowIdx));
      }
      printf("}\n");
    }

    // Insert the upper part of the selected row in U
    // including the diagonal term.
    value_type   diag;
    size_type U_entryIdx = row_mapU(factorization_step);
    for(size_type entryIdx = A.graph.row_map(selected_row);
        entryIdx < A.graph.row_map(selected_row + 1);
        ++entryIdx) {
      if(permutation_inv(A.graph.entries(entryIdx)) >= factorization_step) {
        if (A.values(entryIdx) != Kokkos::ArithTraits<value_type>::zero())
        {
          entriesU(U_entryIdx) = A.graph.entries(entryIdx);
          valuesU(U_entryIdx) = A.values(entryIdx);
          ++U_entryIdx;
        }
        if(A.graph.entries(entryIdx) == selected_row) {
          diag = A.values(entryIdx);
        }
      }
    }
    row_mapU(factorization_step + 1) = U_entryIdx;

    if(verbosity > 0) {
      printf("Diagonal values of row %d is %f\n", selected_row, diag);
    }

    if(verbosity > 2) {
      printf("U, row_map={ ");
      for(ordinal_type rowIdx = 0; rowIdx < factorization_step + 1; ++rowIdx) {
        printf("%d ", int(row_mapU(rowIdx)));
      }
      printf("}, entries={ ");
      for(size_type entryIdx = row_mapU(0); entryIdx < row_mapU(factorization_step + 1); ++entryIdx) {
        printf("%d ", int(entriesU(entryIdx)));
      }
      printf("}, values={ ");
      for(size_type entryIdx = row_mapU(0); entryIdx < row_mapU(factorization_step + 1); ++entryIdx) {
        printf("%f ", valuesU(entryIdx));
      }
      printf("}\n");
    }

    // Insert the lower part of the selected column of A
    // divided by its the diagonal value to obtain a unit
    // diagonal value in L.
    size_type L_entryIdx = row_mapL(factorization_step);
    entriesL(L_entryIdx) = selected_row;
    valuesL(L_entryIdx)  = Kokkos::ArithTraits<value_type>::one();
    ++L_entryIdx;
    for(size_type entryIdx = At.graph.row_map(selected_row);
        entryIdx < At.graph.row_map(selected_row + 1);
        ++entryIdx) {
      if(permutation_inv(At.graph.entries(entryIdx)) > factorization_step && At.values(entryIdx) != Kokkos::ArithTraits<value_type>::zero()) {
        entriesL(L_entryIdx) = At.graph.entries(entryIdx);
        valuesL(L_entryIdx) = At.values(entryIdx) / diag;
        ++L_entryIdx;
      }
    }
    row_mapL(factorization_step + 1) = L_entryIdx;

    if(verbosity > 2) {
      printf("L(%d), [row_map(%d), row_map(%d)[ = [%d, %d[, entries={ ",
             int(factorization_step), int(factorization_step), int(factorization_step+1),
             int(row_mapL(factorization_step)), int(row_mapL(factorization_step+1)));
      for(size_type entryIdx = row_mapL(factorization_step); entryIdx < row_mapL(factorization_step + 1); ++entryIdx) {
        printf("%d ", int(entriesL(entryIdx)));
      }
      printf("}, values={ ");
      for(size_type entryIdx = row_mapL(factorization_step); entryIdx < row_mapL(factorization_step + 1); ++entryIdx) {
        printf("%f ", valuesL(entryIdx));
      }
      printf("}\n");
    }

    // If this was the last row no need to update A and At!
    if(factorization_step == A.numRows()-1) {return;}

    // Finally we want to update A and At with the values
    // that where not discarded during factorization.
    // Note: this is almost the same operation as computing
    // the norm of the discarded fill...

    // First step: find the diagonal entry in selected_row
    value_type diag_val;
    for(size_type entryIdx = A.graph.row_map(selected_row); entryIdx < A.graph.row_map(selected_row + 1); ++entryIdx) {
      ordinal_type colIdx = A.graph.entries(entryIdx);
      if(selected_row == colIdx) {
        diag_val = A.values(entryIdx);
      }
    }

    // Extract alpha and beta vectors
    // Then insert alpha*beta/diag_val if the corresponding
    // entry in A is non-zero.
    for(size_type alphaIdx = At.graph.row_map(selected_row); alphaIdx < At.graph.row_map(selected_row + 1); ++alphaIdx) {
      ordinal_type fillRowIdx = At.graph.entries(alphaIdx);
      bool row_not_eliminated = true;
      for(ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
        if(fillRowIdx == permutation(stepIdx)) {
          row_not_eliminated = false;
        }
      }

      if((fillRowIdx != selected_row) && row_not_eliminated) {
        for(size_type betaIdx = A.graph.row_map(selected_row); betaIdx < A.graph.row_map(selected_row + 1); ++betaIdx) {
          ordinal_type fillColIdx = A.graph.entries(betaIdx);
          bool col_not_eliminated = true;
          for(ordinal_type stepIdx = 0; stepIdx < factorization_step; ++stepIdx) {
            if(fillColIdx == permutation(stepIdx)) {
              col_not_eliminated = false;
            }
          }

          if((fillColIdx != selected_row) && col_not_eliminated) {
            for(size_type entryIdx = A.graph.row_map(fillRowIdx); entryIdx < A.graph.row_map(fillRowIdx + 1); ++entryIdx) {
              if(A.graph.entries(entryIdx) == fillColIdx) {
                A.values(entryIdx) -= At.values(alphaIdx)*A.values(betaIdx) / diag_val;

                if(verbosity > 1) {
                  printf("A[%d, %d] -= %f\n", int(fillRowIdx), int(fillColIdx),
                         At.values(alphaIdx)*A.values(betaIdx) / diag_val);
                }
              }
            }

            for(size_type entryIdx = At.graph.row_map(fillColIdx); entryIdx < At.graph.row_map(fillColIdx + 1); ++entryIdx) {
              if(At.graph.entries(entryIdx) == fillRowIdx) {
                At.values(entryIdx) -= At.values(alphaIdx)*A.values(betaIdx) / diag_val;
              }
            }
          }
        }
      }
    }

    if(verbosity > 0) {
      printf("New values in A: { ");
      for(size_type entryIdx = 0; entryIdx < A.nnz(); ++entryIdx) {
        printf("%f ", A.values(entryIdx));
      }
      printf("}\n");
      printf("New values in At: { ");
      for(size_type entryIdx = 0; entryIdx < At.nnz(); ++entryIdx) {
        printf("%f ", At.values(entryIdx));
      }
      printf("}\n");
    }
  } // operator()

}; // MDF_factorize_row

template<class col_ind_type>
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

template<class matrix_type>
struct MDF_handle {
  using crs_matrix_type = matrix_type;
  using execution_space = typename matrix_type::execution_space;
  using row_map_type    = typename crs_matrix_type::StaticCrsGraphType::row_map_type::non_const_type;
  using col_ind_type    = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type     = typename crs_matrix_type::values_type::non_const_type;
  using size_type       = typename crs_matrix_type::size_type;
  using ordinal_type    = typename crs_matrix_type::ordinal_type;

  ordinal_type numRows;

  // Views needed to construct L and U
  // at the end of the numerical phase.
  row_map_type row_mapL, row_mapU;
  col_ind_type entriesL, entriesU;
  values_type  valuesL, valuesU;

  // Row permutation that defines
  // the MDF ordering or order of
  // elimination during the factorization.
  col_ind_type permutation, permutation_inv;

  int verbosity;


  MDF_handle(const crs_matrix_type A) : numRows(A.numRows()),
                                        permutation(col_ind_type("row permutation", A.numRows())),
                                        permutation_inv(col_ind_type("inverse row permutation", A.numRows())),
                                        verbosity(0) {};

  void set_verbosity(const int verbosity_level) {verbosity = verbosity_level;}

  void allocate_data(const size_type nnzL,
                     const size_type nnzU) {
    
    if(verbosity>0)
      printf("Allocating L (%d) and U (%d)\n",nnzL,nnzU);

    // Allocate L
    row_mapL = row_map_type("MDF_handle::rowMapL", numRows + 1);
    entriesL = col_ind_type("MDF_handle::entriesL", nnzL);
    valuesL  = values_type("MDF_handle::valuesL",   nnzL);

    // Allocate U
    row_mapU = row_map_type("MDF_handle::rowMapU", numRows + 1);
    entriesU = col_ind_type("MDF_handle::entriesU", nnzU);
    valuesU  = values_type("MDF_handle::valuesU",   nnzU);
  }

  col_ind_type get_permutation() {return permutation;}
  col_ind_type get_permutation_inv() {return permutation_inv;}

  void sort_factors() {
    KokkosKernels::sort_crs_matrix<execution_space, row_map_type, col_ind_type, values_type>
      (row_mapL, entriesL, valuesL);
    KokkosKernels::sort_crs_matrix<execution_space, row_map_type, col_ind_type, values_type>
      (row_mapU, entriesU, valuesU);
  }

  crs_matrix_type getL() {
    crs_matrix_type Lt("L", numRows, numRows, entriesL.extent(0),
                       valuesL, row_mapL, entriesL);

    crs_matrix_type L = KokkosKernels::Impl::transpose_matrix<crs_matrix_type>(Lt);
    KokkosKernels::sort_crs_matrix<crs_matrix_type>(L);
    return L;
  }

  crs_matrix_type getU() {
    crs_matrix_type U("U", numRows, numRows, entriesU.extent(0),
                      valuesU, row_mapU, entriesU);
    KokkosKernels::sort_crs_matrix<crs_matrix_type>(U);
    return U;
  }

};

template<class crs_matrix_type, class MDF_handle>
void mdf_symbolic_phase(const crs_matrix_type& A, MDF_handle& handle) {
  using row_map_type = typename crs_matrix_type::StaticCrsGraphType::row_map_type::non_const_type;
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using size_type    = typename crs_matrix_type::size_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;

  using execution_space = typename crs_matrix_type::execution_space;
  using range_policy_type = Kokkos::RangePolicy<ordinal_type, execution_space>;

  // Symbolic phase:
  // compute transpose of A for easy access to columns of A
  // allocate temporaries
  // allocate L and U
  size_type nnzL = 0, nnzU = 0;
  range_policy_type setupPolicy(0, A.numRows());
  MDF_count_lower<crs_matrix_type> compute_nnzL(A, handle.permutation, handle.permutation_inv);
  Kokkos::parallel_reduce(range_policy_type(0, A.numRows()), compute_nnzL, nnzL);
  nnzU = A.nnz() - nnzL + A.numRows();
  handle.allocate_data(nnzL, nnzU);

  if(handle.verbosity > 0) {
    printf("MDF symbolic:  nnzL = %d, nnzU = %d\n",
           static_cast<int>(nnzL),
           static_cast<int>(nnzU));
  }

  return;
} // mdf_symbolic_phase

template<class crs_matrix_type, class MDF_handle>
void mdf_numeric_phase(const crs_matrix_type& A, MDF_handle& handle) {
  using col_ind_type = typename crs_matrix_type::StaticCrsGraphType::entries_type::non_const_type;
  using values_type  = typename crs_matrix_type::values_type::non_const_type;
  using ordinal_type = typename crs_matrix_type::ordinal_type;
  using value_type   = typename crs_matrix_type::value_type;

  using execution_space = typename crs_matrix_type::execution_space;
  using range_policy_type = Kokkos::RangePolicy<ordinal_type, execution_space>;

  // Numerical phase:
  // loop over rows
  //   compute discarded fill of each row
  //   selected pivot based on MDF
  //   factorize pivot row of A
  crs_matrix_type Atmp = crs_matrix_type("A fill", A);
  crs_matrix_type At = KokkosKernels::Impl::transpose_matrix<crs_matrix_type>(A);
  KokkosKernels::sort_crs_matrix<crs_matrix_type>(At);
  values_type  discarded_fill("discarded fill", A.numRows());
  col_ind_type deficiency("deficiency", A.numRows());

  const int verbosity_level = handle.verbosity;
  for(ordinal_type factorization_step = 0; factorization_step < A.numRows(); ++factorization_step) {
    if(verbosity_level > 0) {
      printf("\n\nFactorization step %d\n\n", static_cast<int>(factorization_step));
    }

    range_policy_type stepPolicy(factorization_step, Atmp.numRows());
    Kokkos::deep_copy(discarded_fill, Kokkos::ArithTraits<value_type>::max());
    Kokkos::deep_copy(deficiency, Kokkos::ArithTraits<ordinal_type>::max());
    MDF_discarded_fill_norm<crs_matrix_type> MDF_df_norm(Atmp, At, factorization_step,
                                                         handle.permutation,
                                                         discarded_fill, deficiency,
                                                         verbosity_level);
    Kokkos::parallel_for(stepPolicy, MDF_df_norm);

    ordinal_type selected_row_idx = 0;
    MDF_select_row<crs_matrix_type> MDF_row_selector(factorization_step, discarded_fill,
                                                     deficiency, Atmp.graph.row_map,
                                                     handle.permutation);
    Kokkos::parallel_reduce(stepPolicy, MDF_row_selector, selected_row_idx);

    MDF_factorize_row<crs_matrix_type> factorize_row(Atmp, At,
                                                     handle.row_mapL, handle.entriesL, handle.valuesL,
                                                     handle.row_mapU, handle.entriesU, handle.valuesU,
                                                     handle.permutation, handle.permutation_inv,
                                                     selected_row_idx, factorization_step,
                                                     verbosity_level);
    Kokkos::parallel_for(range_policy_type(0, 1), factorize_row);

    if(verbosity_level > 0) {
      printf("\n");
    }
  }

  if(verbosity_level > 0) printf("Reindexing U\n");
  MDF_reindex_matrix<col_ind_type> reindex_U(handle.permutation_inv, handle.entriesU);
  Kokkos::parallel_for(range_policy_type(0, handle.entriesU.extent(0)),
                       reindex_U);

  if(verbosity_level > 0) printf("Reindexing L\n");
  MDF_reindex_matrix<col_ind_type> reindex_L(handle.permutation_inv, handle.entriesL);
  Kokkos::parallel_for(range_policy_type(0, handle.entriesL.extent(0)),
                       reindex_L);

  return;
} // mdf_numeric_phase

} // namespace MDFImpl
} // namespace Ifpack2

#endif /* IFPACK2_MDF_HANDLE_HPP */
