// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKCOMPUTERES_AND_SOLVE_DECL_HPP
#define IFPACK2_BLOCKCOMPUTERES_AND_SOLVE_DECL_HPP

#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_BlockCrsMatrix_decl.hpp>
#include "Ifpack2_BlockHelper.hpp"
#include "Ifpack2_BlockHelper_ETI.hpp"

namespace Ifpack2::BlockHelperDetails {

template <typename ExecSpace, typename DiagOffsets, typename Rowptrs,
          typename Entries>
DiagOffsets findDiagOffsets(const Rowptrs& rowptrs, const Entries& entries,
                            size_t nrows, int blocksize) {
  DiagOffsets diag_offsets(do_not_initialize_tag("btdm.diag_offsets"), nrows);
  int err1 = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace>(0, nrows),
      KOKKOS_LAMBDA(size_t row, int& err2) {
        auto rowBegin = rowptrs(row);
        auto rowEnd   = rowptrs(row + 1);
        for (size_t j = rowBegin; j < rowEnd; j++) {
          if (size_t(entries(j)) == row) {
            diag_offsets(row) = j * blocksize * blocksize;
            return;
          }
        }
        err2++;
      },
      err1);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      err1, "Ifpack2 BTD: at least one row has no diagonal entry");
  return diag_offsets;
}

template <typename MatrixType,
          typename ImplTagType = typename BlockTriDiContainerDetails::ImplTag<typename MatrixType::scalar_type>::type>
struct ComputeResidualAndSolve;

template <typename MatrixType>
struct ComputeResidualAndSolve<MatrixType, BlockTriDiContainerDetails::ImplSimdTag> {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using size_type          = typename impl_type::size_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
  using magnitude_type     = typename impl_type::magnitude_type;
  /// views
  using local_ordinal_type_1d_view =
      typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type =
      typename impl_type::tpetra_block_access_view_type;  // block crs (layout
                                                          // right)
  using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra =
      typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector
                                                            // (layout left)
  using btdm_scalar_type_3d_view = typename impl_type::btdm_scalar_type_3d_view;
  using btdm_scalar_type_4d_view = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view              = typename impl_type::i64_3d_view;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  ComputeResidualAndSolve(const AmD<MatrixType>& amd_,
                          const btdm_scalar_type_3d_view& d_inv_,
                          const impl_scalar_type_1d_view& W_,
                          const local_ordinal_type& blocksize_requested_,
                          const impl_scalar_type& damping_factor_)
    : amd(amd_)
    , blocksize_requested(blocksize_requested_)
    , d_inv(d_inv_)
    , W(W_)
    , damping_factor(damping_factor_) {}

  void run_y_zero(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& y_);
  void run_single_pass(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& x_remote_,
      const impl_scalar_type_2d_view_tpetra& y_);
  void run_pass1_of_2(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& y_);
  void run_pass2_of_2(
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& x_remote_,
      const impl_scalar_type_2d_view_tpetra& y_);

 private:
  // AmD information
  const AmD<MatrixType> amd;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // diagonal block inverses
  const btdm_scalar_type_3d_view d_inv;

  // squared update norms
  const impl_scalar_type_1d_view W;

  impl_scalar_type damping_factor;
};

template <typename MatrixType>
struct ComputeResidualAndSolve<MatrixType, BlockTriDiContainerDetails::ImplNotAvailTag> {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using size_type          = typename impl_type::size_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
  using magnitude_type     = typename impl_type::magnitude_type;
  /// views
  using local_ordinal_type_1d_view =
      typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type =
      typename impl_type::tpetra_block_access_view_type;  // block crs (layout
                                                          // right)
  using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra =
      typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector
                                                            // (layout left)
  using btdm_scalar_type_3d_view = typename impl_type::btdm_scalar_type_3d_view;
  using btdm_scalar_type_4d_view = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view              = typename impl_type::i64_3d_view;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  ComputeResidualAndSolve(const AmD<MatrixType>& amd_,
                          const btdm_scalar_type_3d_view& d_inv_,
                          const impl_scalar_type_1d_view& W_,
                          const local_ordinal_type& blocksize_requested_,
                          const impl_scalar_type& damping_factor_)
    : amd(amd_)
    , blocksize_requested(blocksize_requested_)
    , d_inv(d_inv_)
    , W(W_)
    , damping_factor(damping_factor_) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: BlockTriDiContainer and related classes are not available for this scalar_type");
  }

  void run_y_zero(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& y_) {}
  void run_single_pass(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& x_remote_,
      const impl_scalar_type_2d_view_tpetra& y_) {}
  void run_pass1_of_2(
      const Const<impl_scalar_type_2d_view_tpetra>& b_,
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& y_) {}
  void run_pass2_of_2(
      const impl_scalar_type_2d_view_tpetra& x_,
      const impl_scalar_type_2d_view_tpetra& x_remote_,
      const impl_scalar_type_2d_view_tpetra& y_) {}

 private:
  // AmD information
  const AmD<MatrixType> amd;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // diagonal block inverses
  const btdm_scalar_type_3d_view d_inv;

  // squared update norms
  const impl_scalar_type_1d_view W;

  impl_scalar_type damping_factor;
};

}  // namespace Ifpack2::BlockHelperDetails

#endif
