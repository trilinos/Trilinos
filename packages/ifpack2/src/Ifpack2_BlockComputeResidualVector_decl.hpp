// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKCOMPUTERES_VECTOR_DECL_HPP
#define IFPACK2_BLOCKCOMPUTERES_VECTOR_DECL_HPP

#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_BlockCrsMatrix_decl.hpp>
#include "Ifpack2_BlockHelper.hpp"

namespace Ifpack2::BlockHelperDetails {

template <typename MatrixType>
struct ComputeResidualVector {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type  = typename impl_type::local_ordinal_type;
  using size_type           = typename impl_type::size_type;
  using impl_scalar_type    = typename impl_type::impl_scalar_type;
  using magnitude_type      = typename impl_type::magnitude_type;
  using btdm_scalar_type    = typename impl_type::btdm_scalar_type;
  using btdm_magnitude_type = typename impl_type::btdm_magnitude_type;
  /// views
  using local_ordinal_type_1d_view      = typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view               = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type   = typename impl_type::tpetra_block_access_view_type;  // block crs (layout right)
  using impl_scalar_type_1d_view        = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra = typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector (layout left)
  using vector_type_3d_view             = typename impl_type::vector_type_3d_view;
  using btdm_scalar_type_4d_view        = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view                     = typename impl_type::i64_3d_view;
  static constexpr int vector_length    = impl_type::vector_length;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  // AmD information
  const ConstUnmanaged<size_type_1d_view> rowptr, rowptr_remote;
  const ConstUnmanaged<local_ordinal_type_1d_view> colindsub, colindsub_remote;
  const ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

  // block crs graph information
  // for cuda (kokkos crs graph uses a different size_type from size_t)
  const ConstUnmanaged<Kokkos::View<size_t *, node_device_type>> A_block_rowptr;
  const ConstUnmanaged<Kokkos::View<size_t *, node_device_type>> A_point_rowptr;
  const ConstUnmanaged<Kokkos::View<local_ordinal_type *, node_device_type>> A_colind;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // part interface
  const ConstUnmanaged<local_ordinal_type_1d_view> part2packrowidx0;
  const ConstUnmanaged<local_ordinal_type_1d_view> part2rowidx0;
  const ConstUnmanaged<local_ordinal_type_1d_view> rowidx2part;
  const ConstUnmanaged<local_ordinal_type_1d_view> partptr;
  const ConstUnmanaged<local_ordinal_type_1d_view> lclrow;
  const ConstUnmanaged<local_ordinal_type_1d_view> dm2cm;

  // block offsets
  const ConstUnmanaged<i64_3d_view> A_x_offsets;
  const ConstUnmanaged<i64_3d_view> A_x_offsets_remote;

  const bool is_dm2cm_active;
  const bool hasBlockCrsMatrix;

  template <typename LocalCrsGraphType>
  ComputeResidualVector(const AmD<MatrixType> &amd,
                        const LocalCrsGraphType &block_graph,
                        const LocalCrsGraphType &point_graph,
                        const local_ordinal_type &blocksize_requested_,
                        const PartInterface<MatrixType> &interf,
                        const local_ordinal_type_1d_view &dm2cm_,
                        bool hasBlockCrsMatrix_)
    : rowptr(amd.rowptr)
    , rowptr_remote(amd.rowptr_remote)
    , colindsub(amd.A_colindsub)
    , colindsub_remote(amd.A_colindsub_remote)
    , tpetra_values(amd.tpetra_values)
    , A_block_rowptr(block_graph.row_map)
    , A_point_rowptr(point_graph.row_map)
    , A_colind(block_graph.entries)
    , blocksize_requested(blocksize_requested_)
    , part2packrowidx0(interf.part2packrowidx0)
    , part2rowidx0(interf.part2rowidx0)
    , rowidx2part(interf.rowidx2part)
    , partptr(interf.partptr)
    , lclrow(interf.lclrow)
    , dm2cm(dm2cm_)
    , A_x_offsets(amd.A_x_offsets)
    , A_x_offsets_remote(amd.A_x_offsets_remote)
    , is_dm2cm_active(dm2cm_.span() > 0)
    , hasBlockCrsMatrix(hasBlockCrsMatrix_) {}

  // Precompute offsets of each A and x entry to speed up residual.
  // (Applies for hasBlockCrsMatrix == true and OverlapTag/AsyncTag)
  // Reading A, x take up to 4, 6 levels of indirection respectively,
  // but precomputing the offsets reduces it to 2 for both.
  //
  // This function allocates and populates these members of AmD:
  // A_x_offsets, A_x_offsets_remote
  static void precompute_A_x_offsets(
      AmD<MatrixType> &amd,
      const PartInterface<MatrixType> &interf,
      const Teuchos::RCP<const typename ImplType<MatrixType>::tpetra_crs_graph_type> &g,
      const typename ImplType<MatrixType>::local_ordinal_type_1d_view &dm2cm,
      int blocksize,
      bool ownedRemoteSeparate);

  // y = b - Rx; seq method
  void run(const impl_scalar_type_2d_view_tpetra &y_,
           const Const<impl_scalar_type_2d_view_tpetra> &b_,
           const impl_scalar_type_2d_view_tpetra &x_);

  // y = b - R (x , x_remote)
  void run(const vector_type_3d_view &y_packed_,
           const Const<impl_scalar_type_2d_view_tpetra> &b_,
           const impl_scalar_type_2d_view_tpetra &x_,
           const impl_scalar_type_2d_view_tpetra &x_remote_);

  // y = b - R (y , y_remote)
  void run(const vector_type_3d_view &y_packed_,
           const Const<impl_scalar_type_2d_view_tpetra> &b_,
           const impl_scalar_type_2d_view_tpetra &x_,
           const impl_scalar_type_2d_view_tpetra &x_remote_,
           const bool compute_owned);
};

}  // namespace Ifpack2::BlockHelperDetails

#endif
