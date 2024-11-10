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
#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosSparse_gauss_seidel_spec.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {

namespace Experimental {

///
/// @brief Gauss-Seidel preconditioner setup (first phase, based on sparsity
/// pattern only)
///
/// @tparam ExecutionSpace This kernels execution space type.
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @param space The execution space instance this kernel will be run on.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @pre   <tt>handle->create_gs_handle(...)</tt> has been called previously
///
template <typename ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
void gauss_seidel_symbolic(const ExecutionSpace &space, KernelHandle *handle,
                           typename KernelHandle::const_nnz_lno_t num_rows,
                           typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                           lno_nnz_view_t_ entries, bool is_graph_symmetric = true) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: Size type of the matrix "
                "should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: lno type of the matrix "
                "should be same as kernelHandle lno_t.");

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  // const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  // Internal_alno_row_view_t_ const_a_r  = row_map;
  // Internal_alno_nnz_view_t_ const_a_l  = entries;
  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_SYMBOLIC<ExecutionSpace, const_handle_type, Internal_alno_row_view_t_,
                        Internal_alno_nnz_view_t_>::gauss_seidel_symbolic(space, &tmp_handle, num_rows, num_cols,
                                                                          const_a_r, const_a_l, is_graph_symmetric);
}

///
/// @brief Gauss-Seidel preconditioner setup (first phase, based on sparsity
/// pattern only)
///
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @pre   <tt>handle->create_gs_handle(...)</tt> has been called previously
///
template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
void gauss_seidel_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                           typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                           lno_nnz_view_t_ entries, bool is_graph_symmetric = true) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  gauss_seidel_symbolic(my_exec_space, handle, num_rows, num_cols, row_map, entries, is_graph_symmetric);
}

///
/// @brief Block Gauss-Seidel preconditioner setup (first phase, based on
/// sparsity pattern only)
///
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param block_size The number of degrees of freedom per block
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @pre   <tt>handle->create_gs_handle(...)</tt> has been called previously
template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
void block_gauss_seidel_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                 typename KernelHandle::const_nnz_lno_t num_cols,
                                 typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map,
                                 lno_nnz_view_t_ entries, bool is_graph_symmetric = true) {
  auto gsHandle = handle->get_point_gs_handle();
  if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
    throw std::runtime_error(
        "Block versions of Gauss-Seidel are incompatible with algorithm "
        "GS_CLUSTER");
  }
  gsHandle->set_block_size(block_size);

  gauss_seidel_symbolic(handle, num_rows, num_cols, row_map, entries, is_graph_symmetric);
}

///
/// @brief Gauss-Seidel preconditioner setup (second phase, based on matrix's
/// numeric values)
///
/// @tparam ExecutionSpace This kernels execution space type.
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @param space The execution space instance this kernel will be run on.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
///
template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
          typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle,
                          typename KernelHandle::const_nnz_lno_t num_rows,
                          typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                          lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric = true) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: Size type of the matrix "
                "should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: lno type of the matrix "
                "should be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::gauss_seidel_symbolic: scalar type of the matrix should "
      "be same as kernelHandle scalar_t.");

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  // const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_v(values.data(), values.extent(0));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_NUMERIC<ExecutionSpace, const_handle_type, format, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                       Internal_ascalar_nnz_view_t_>::gauss_seidel_numeric(space, &tmp_handle, num_rows, num_cols,
                                                                           const_a_r, const_a_l, const_a_v,
                                                                           is_graph_symmetric);
}

///
/// @brief Gauss-Seidel preconditioner setup (second phase, based on matrix's
/// numeric values). This version accepts the matrix's inverse diagonal from the
/// user.
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type. The user-provided
/// inverse diagonal must share this type.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @remark If the inverse diagonal is not already available, it's best to call
/// the version of <tt>gauss_seidel_numeric</tt> that
///         doesn't take it as an argument. The inverse diagonal will be
///         computed internally.
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                          typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                          lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric = true) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  gauss_seidel_numeric<decltype(my_exec_space), format>(my_exec_space, handle, num_rows, num_cols, row_map, entries,
                                                        values, is_graph_symmetric);
}

///
/// @brief Gauss-Seidel preconditioner setup (second phase, based on matrix's
/// numeric values). This version accepts the matrix's inverse diagonal from the
/// user.
///
/// @tparam ExecutionSpace This kernels execution space type.
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type. The user-provided
/// inverse diagonal must share this type.
/// @param space The execution space instance this kernel will be run on.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param given_inverse_diagonal The inverse (reciprocal) of diagonal
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @remark If the inverse diagonal is not already available, it's best to call
/// the version of <tt>gauss_seidel_numeric</tt> that
///         doesn't take it as an argument. The inverse diagonal will be
///         computed internally.
template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
          typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle,
                          typename KernelHandle::const_nnz_lno_t num_rows,
                          typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                          lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal,
                          bool is_graph_symmetric = true) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: Size type of the matrix "
                "should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::gauss_seidel_symbolic: lno type of the matrix "
                "should be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::gauss_seidel_symbolic: scalar type of the matrix should "
      "be same as kernelHandle scalar_t.");

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  // const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_v(values.data(), values.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_d(given_inverse_diagonal.data(), given_inverse_diagonal.extent(0));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_NUMERIC<ExecutionSpace, const_handle_type, format, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                       Internal_ascalar_nnz_view_t_>::gauss_seidel_numeric(space, &tmp_handle, num_rows, num_cols,
                                                                           const_a_r, const_a_l, const_a_v, const_a_d,
                                                                           is_graph_symmetric);
}

///
/// @brief Gauss-Seidel preconditioner setup (second phase, based on matrix's
/// numeric values). This version accepts the matrix's inverse diagonal from the
/// user.
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type. The user-provided
/// inverse diagonal must share this type.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param given_inverse_diagonal The inverse (reciprocal) of diagonal
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
/// @remark If the inverse diagonal is not already available, it's best to call
/// the version of <tt>gauss_seidel_numeric</tt> that
///         doesn't take it as an argument. The inverse diagonal will be
///         computed internally.
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                          typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                          lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal,
                          bool is_graph_symmetric = true) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  gauss_seidel_numeric<decltype(my_exec_space), format>(my_exec_space, handle, num_rows, num_cols, row_map, entries,
                                                        values, given_inverse_diagonal, is_graph_symmetric);
}

///
/// @brief Block Gauss-Seidel preconditioner setup (second phase, based on
/// matrix's numeric values)
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @param handle handle A KokkosKernelsHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param block_size The number of degrees of freedom per block
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param is_graph_symmetric Whether the upper-left <tt>num_rows x
/// num_rows</tt> submatrix of A is structurally symmetric
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::BSR, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void block_gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                typename KernelHandle::const_nnz_lno_t num_cols,
                                typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map,
                                lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric = true) {
  auto gsHandle = handle->get_point_gs_handle();
  if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
    throw std::runtime_error(
        "Block versions of Gauss-Seidel are incompatible with algorithm "
        "GS_CLUSTER");
  }
  gsHandle->set_block_size(block_size);

  gauss_seidel_numeric<format>(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
}

///
/// @brief Apply symmetric (forward + backward) Gauss-Seidel preconditioner to
/// system AX=Y
///
/// @tparam ExecutionSpace This kernels execution space type.
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param space The execution space instance this kernel will be run
/// on. NOTE: Currently only used for GS_DEFAULT.
/// @param handle handle A KokkosKernelsHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run (forward and backward counts as 1)
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
          typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
          typename x_scalar_view_t, typename y_scalar_view_t>
void symmetric_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle,
                                  typename KernelHandle::const_nnz_lno_t num_rows,
                                  typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                  lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec,
                                  y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector,
                                  typename KernelHandle::nnz_scalar_t omega, int numIter) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: Size type of the "
                "matrix should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: lno type of the "
                "matrix should be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the matrix "
      "should be same as kernelHandle scalar_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename y_scalar_view_t::const_value_type>::value,
      "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of "
      "the y-vector should be same as kernelHandle scalar_t.");

  static_assert(std::is_same<typename KernelHandle::nnz_scalar_t, typename x_scalar_view_t::value_type>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the x-vector "
                "should be same as kernelHandle non-const scalar_t.");

  static_assert(!std::is_same<typename lno_row_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: row_map must have "
                "a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename lno_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: entries must have "
                "a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename scalar_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::symmetric_gauss_seidel_apply: values must have "
                "a contiguous layout (Left or Right, not Stride)");

  // Check compatibility of #vectors
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::symmetric_gauss_seidel_apply: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  typedef Kokkos::View<typename y_scalar_view_t::const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
                       typename y_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_yscalar_nnz_view_t_;

  typedef Kokkos::View<typename x_scalar_view_t::non_const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
                       typename x_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_xscalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_v(values.data(), values.extent(0));

  Internal_xscalar_nnz_view_t_ nonconst_x_v(x_lhs_output_vec.data(), x_lhs_output_vec.extent(0),
                                            x_lhs_output_vec.extent(1));
  Internal_yscalar_nnz_view_t_ const_y_v(y_rhs_input_vec.data(), y_rhs_input_vec.extent(0), y_rhs_input_vec.extent(1));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_APPLY<ExecutionSpace, const_handle_type, format, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                     Internal_ascalar_nnz_view_t_, Internal_xscalar_nnz_view_t_,
                     Internal_yscalar_nnz_view_t_>::gauss_seidel_apply(space, &tmp_handle, num_rows, num_cols,
                                                                       const_a_r, const_a_l, const_a_v, nonconst_x_v,
                                                                       const_y_v, init_zero_x_vector, update_y_vector,
                                                                       omega, numIter, true, true);
}

///
/// @brief Apply symmetric (forward + backward) Gauss-Seidel preconditioner to
/// system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle handle A KokkosKernelsHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run (forward and backward counts as 1)
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void symmetric_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                  typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                  lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec,
                                  y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector,
                                  typename KernelHandle::nnz_scalar_t omega, int numIter) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  symmetric_gauss_seidel_apply<decltype(my_exec_space), format>(my_exec_space, handle, num_rows, num_cols, row_map,
                                                                entries, values, x_lhs_output_vec, y_rhs_input_vec,
                                                                init_zero_x_vector, update_y_vector, omega, numIter);
}

///
/// @brief Apply symmetric (forward + backward) Block Gauss-Seidel
/// preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle handle A KokkosKernelsHandle instance.
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param block_size The number of degrees of freedom per block
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run (forward and backward counts as 1)
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::BSR, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void symmetric_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                        typename KernelHandle::const_nnz_lno_t num_cols,
                                        typename KernelHandle::const_nnz_lno_t block_size,

                                        lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                        x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                        bool init_zero_x_vector, bool update_y_vector,
                                        typename KernelHandle::nnz_scalar_t omega, int numIter) {
  // Check compatibility of dimensions at run time.
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::symmetric_block_gauss_seidel_apply: Dimensions of X "
          "and Y do not match: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  auto gsHandle = handle->get_point_gs_handle();
  if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
    throw std::runtime_error(
        "Block versions of Gauss-Seidel are incompatible with algorithm "
        "GS_CLUSTER");
  }

  gsHandle->set_block_size(block_size);
  symmetric_gauss_seidel_apply<format>(handle, num_rows, num_cols, row_map, entries, values, x_lhs_output_vec,
                                       y_rhs_input_vec, init_zero_x_vector, update_y_vector, omega, numIter);
}

///
/// @brief Apply forward Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param space The execution space instance this kernel will be run
/// on. NOTE: Currently only used for GS_DEFAULT.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
          class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
          typename x_scalar_view_t, typename y_scalar_view_t>
void forward_sweep_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle,
                                      typename KernelHandle::const_nnz_lno_t num_rows,
                                      typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                      lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                      x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                      bool init_zero_x_vector, bool update_y_vector,
                                      typename KernelHandle::nnz_scalar_t omega, int numIter) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: Size type of "
                "the matrix should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: lno type of "
                "the matrix should be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::forward_sweep_gauss_seidel_apply: scalar type of the "
      "matrix should be same as kernelHandle scalar_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename y_scalar_view_t::const_value_type>::value,
      "KokkosSparse::forward_sweep_gauss_seidel_apply: scalar type "
      "of the y-vector should be same as kernelHandle scalar_t.");

  static_assert(std::is_same<typename KernelHandle::nnz_scalar_t, typename x_scalar_view_t::value_type>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: scalar type of the "
                "x-vector should be same as kernelHandle non-const scalar_t.");

  static_assert(!std::is_same<typename lno_row_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: row_map must "
                "have a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename lno_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: entries must "
                "have a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename scalar_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::forward_sweep_gauss_seidel_apply: values must "
                "have a contiguous layout (Left or Right, not Stride)");

  // Check compatibility of dimensions at run time.
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::forward_sweep_gauss_seidel_apply: Dimensions of X and "
          "Y do not match: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  // const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  typedef Kokkos::View<typename y_scalar_view_t::const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
                       typename y_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_yscalar_nnz_view_t_;

  typedef Kokkos::View<typename x_scalar_view_t::non_const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
                       typename x_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_xscalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_v(values.data(), values.extent(0));

  Internal_xscalar_nnz_view_t_ nonconst_x_v(x_lhs_output_vec.data(), x_lhs_output_vec.extent(0),
                                            x_lhs_output_vec.extent(1));
  Internal_yscalar_nnz_view_t_ const_y_v(y_rhs_input_vec.data(), y_rhs_input_vec.extent(0), y_rhs_input_vec.extent(1));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_APPLY<ExecutionSpace, const_handle_type, format, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                     Internal_ascalar_nnz_view_t_, Internal_xscalar_nnz_view_t_,
                     Internal_yscalar_nnz_view_t_>::gauss_seidel_apply(space, &tmp_handle, num_rows, num_cols,
                                                                       const_a_r, const_a_l, const_a_v, nonconst_x_v,
                                                                       const_y_v, init_zero_x_vector, update_y_vector,
                                                                       omega, numIter, true, false);
}

///
/// @brief Apply forward Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, class KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void forward_sweep_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                      typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                      lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                      x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                      bool init_zero_x_vector, bool update_y_vector,
                                      typename KernelHandle::nnz_scalar_t omega, int numIter) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  forward_sweep_gauss_seidel_apply<decltype(my_exec_space), format>(
      my_exec_space, handle, num_rows, num_cols, row_map, entries, values, x_lhs_output_vec, y_rhs_input_vec,
      init_zero_x_vector, update_y_vector, omega, numIter);
}

///
/// @brief Apply forward Block Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param block_size The number of degrees of freedom per block
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::BSR, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void forward_sweep_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                            typename KernelHandle::const_nnz_lno_t num_cols,
                                            typename KernelHandle::const_nnz_lno_t block_size,

                                            lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                            x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                            bool init_zero_x_vector, bool update_y_vector,
                                            typename KernelHandle::nnz_scalar_t omega, int numIter) {
  // Check compatibility of dimensions at run time.
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::forward_sweep_block_gauss_seidel_apply: Dimensions of "
          "X and Y do not match: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  auto gsHandle = handle->get_point_gs_handle();
  if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
    throw std::runtime_error(
        "Block versions of Gauss-Seidel are incompatible with algorithm "
        "GS_CLUSTER");
  }
  gsHandle->set_block_size(block_size);
  forward_sweep_gauss_seidel_apply<format>(handle, num_rows, num_cols, row_map, entries, values, x_lhs_output_vec,
                                           y_rhs_input_vec, init_zero_x_vector, update_y_vector, omega, numIter);
}

///
/// @brief Apply backward Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam ExecutionSpace This kernels execution space type.
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param space The execution space instance this kernel will be run
/// on. NOTE: Currently only used for GS_DEFAULT.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
          class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
          typename x_scalar_view_t, typename y_scalar_view_t>
void backward_sweep_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle,
                                       typename KernelHandle::const_nnz_lno_t num_rows,
                                       typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                       lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                       x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                       bool init_zero_x_vector, bool update_y_vector,
                                       typename KernelHandle::nnz_scalar_t omega, int numIter) {
  static_assert(std::is_same<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: Size type of "
                "the matrix should be same as kernelHandle sizetype.");

  static_assert(std::is_same<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: lno type of "
                "the matrix should be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::backward_sweep_gauss_seidel_apply: scalar type of the "
      "matrix should be same as kernelHandle scalar_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename y_scalar_view_t::const_value_type>::value,
      "KokkosSparse::backward_sweep_gauss_seidel_apply: scalar type "
      "of the y-vector should be same as kernelHandle scalar_t.");

  static_assert(std::is_same<typename KernelHandle::nnz_scalar_t, typename x_scalar_view_t::value_type>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: scalar type of the "
                "x-vector should be same as kernelHandle non-const scalar_t.");

  static_assert(!std::is_same<typename lno_row_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: row_map must "
                "have a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename lno_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: entries must "
                "have a contiguous layout (Left or Right, not Stride)");
  static_assert(!std::is_same<typename scalar_nnz_view_t_::array_layout, Kokkos::LayoutStride>::value,
                "KokkosSparse::backward_sweep_gauss_seidel_apply: values must "
                "have a contiguous layout (Left or Right, not Stride)");

  // Check compatibility of dimensions at run time.
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::backward_sweep_gauss_seidel_apply: Dimensions of X "
          "and Y do not match: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  // const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  typedef Kokkos::View<typename y_scalar_view_t::const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
                       typename y_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_yscalar_nnz_view_t_;

  typedef Kokkos::View<typename x_scalar_view_t::non_const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
                       typename x_scalar_view_t::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_xscalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_map.data(), row_map.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entries.data(), entries.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_v(values.data(), values.extent(0));

  Internal_xscalar_nnz_view_t_ nonconst_x_v(x_lhs_output_vec.data(), x_lhs_output_vec.extent(0),
                                            x_lhs_output_vec.extent(1));
  Internal_yscalar_nnz_view_t_ const_y_v(y_rhs_input_vec.data(), y_rhs_input_vec.extent(0), y_rhs_input_vec.extent(1));

  using namespace KokkosSparse::Impl;

  GAUSS_SEIDEL_APPLY<ExecutionSpace, const_handle_type, format, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                     Internal_ascalar_nnz_view_t_, Internal_xscalar_nnz_view_t_,
                     Internal_yscalar_nnz_view_t_>::gauss_seidel_apply(space, &tmp_handle, num_rows, num_cols,
                                                                       const_a_r, const_a_l, const_a_v, nonconst_x_v,
                                                                       const_y_v, init_zero_x_vector, update_y_vector,
                                                                       omega, numIter, false, true);
}

///
/// @brief Apply backward Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, class KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void backward_sweep_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                       typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                       lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                       x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                       bool init_zero_x_vector, bool update_y_vector,
                                       typename KernelHandle::nnz_scalar_t omega, int numIter) {
  auto my_exec_space = handle->get_gs_handle()->get_execution_space();
  backward_sweep_gauss_seidel_apply<decltype(my_exec_space), format>(
      my_exec_space, handle, num_rows, num_cols, row_map, entries, values, x_lhs_output_vec, y_rhs_input_vec,
      init_zero_x_vector, update_y_vector, omega, numIter);
}

///
/// @brief Apply backward Block Gauss-Seidel preconditioner to system AX=Y
///
/// @tparam format The matrix storage format, CRS or BSR
/// @tparam KernelHandle A specialization of
/// KokkosKernels::Experimental::KokkosKernelsHandle
/// @tparam lno_row_view_t_ The matrix's rowmap type
/// @tparam lno_nnz_view_t_ The matrix's entries type
/// @tparam scalar_nnz_view_t_ The matrix's values type
/// @tparam x_scalar_view_t The type of the X (left-hand side, unknown) vector.
/// May be rank-1 or rank-2 View.
/// @tparam y_scalar_view_t The type of the Y (right-hand side) vector. May be
/// rank-1 or rank-2 View.
/// @param handle KernelHandle instance
/// @param num_rows Number of rows in the matrix
/// @param num_cols Number of columns in the matrix
/// @param block_size The number of degrees of freedom per block
/// @param row_map The matrix's rowmap
/// @param entries The matrix's entries
/// @param values The matrix's values
/// @param x_lhs_output_vec The X (left-hand side, unknown) vector
/// @param y_rhs_input_vec The Y (right-hand side) vector
/// @param init_zero_x_vector Whether to zero out X before applying
/// @param update_y_vector Whether Y has changed since the last call to apply
/// @param omega The damping factor for successive over-relaxation
/// @param numIter How many iterations to run
/// @pre   <tt>x_lhs_output_vec.extent(0) == num_cols</tt>
/// @pre   <tt>y_rhs_input_vec.extent(0) == num_rows</tt>
/// @pre   <tt>x_lhs_output_vec.extent(1) == y_rhs_input_vec.extent(1)</tt>
///
template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::BSR, typename KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
          typename y_scalar_view_t>
void backward_sweep_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                             typename KernelHandle::const_nnz_lno_t num_cols,
                                             typename KernelHandle::const_nnz_lno_t block_size,

                                             lno_row_view_t_ row_map, lno_nnz_view_t_ entries,
                                             scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec,
                                             y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector,
                                             bool update_y_vector, typename KernelHandle::nnz_scalar_t omega,
                                             int numIter) {
  // Check compatibility of dimensions at run time.
  if (x_lhs_output_vec.extent(1) != y_rhs_input_vec.extent(1)) {
    std::ostringstream os;
    os << "KokkosSparse::backward_sweep_block_gauss_seidel_apply: Dimensions "
          "of X and Y do not match: "
       << "X has " << x_lhs_output_vec.extent(1) << "columns, Y has " << y_rhs_input_vec.extent(1) << " columns.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  auto gsHandle = handle->get_point_gs_handle();
  if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
    throw std::runtime_error(
        "Block versions of Gauss-Seidel are incompatible with algorithm "
        "GS_CLUSTER");
  }
  gsHandle->set_block_size(block_size);
  backward_sweep_gauss_seidel_apply<format>(handle, num_rows, num_cols, row_map, entries, values, x_lhs_output_vec,
                                            y_rhs_input_vec, init_zero_x_vector, update_y_vector, omega, numIter);
}
}  // namespace Experimental
}  // namespace KokkosSparse
#endif
