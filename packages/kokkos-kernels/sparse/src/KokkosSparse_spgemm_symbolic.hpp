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
#ifndef _KOKKOS_SPGEMM_SYMBOLIC_HPP
#define _KOKKOS_SPGEMM_SYMBOLIC_HPP

#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_spgemm_symbolic_spec.hpp"
#include "KokkosSparse_Utils.hpp"

namespace KokkosSparse {

namespace Experimental {

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename blno_row_view_t_,
          typename blno_nnz_view_t_, typename clno_row_view_t_>
void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                     typename KernelHandle::const_nnz_lno_t n, typename KernelHandle::const_nnz_lno_t k,
                     alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, bool transposeA, blno_row_view_t_ row_mapB,
                     blno_nnz_view_t_ entriesB, bool transposeB, clno_row_view_t_ row_mapC,
                     bool computeRowptrs = false) {
  static_assert(
      std::is_same<typename clno_row_view_t_::value_type, typename clno_row_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Output matrix rowmap must be non-const.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename alno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of left handside matrix should "
      "be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename blno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of right handside matrix "
      "should be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename clno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of output matrix should be "
      "same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename alno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: lno type of left handside matrix should "
      "be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename blno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: lno type of right handside matrix should "
      "be same as kernelHandle lno_t.");

  if (transposeA || transposeB) {
    throw std::runtime_error(
        "SpGEMM is not implemented for Transposes yet. "
        "If you need this case please let kokkos-kernels developers know.\n");
  }

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KernelHandle::HandleExecSpace c_exec_t;
  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;
  typedef typename Kokkos::Device<c_exec_t, c_temp_t> UniformDevice_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t,
                                                                    c_persist_t>
      const_handle_type;
  const_handle_type tmp_handle(*handle);

  typedef Kokkos::View<typename alno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_row_view_t_>::array_layout,
                       UniformDevice_t,  //      typename
                                         //      alno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //      typename
                                         //      alno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename blno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout,
                       UniformDevice_t,  //      typename
                                         //      blno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_row_view_t_;

  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //      typename
                                         //      blno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_nnz_view_t_;

  // static assert clno_row_view_t_ cannot be const type.
  typedef Kokkos::View<typename clno_row_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout,
                       UniformDevice_t,  //      typename
                                         //      clno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_clno_row_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_mapA.data(), row_mapA.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entriesA.data(), entriesA.extent(0));
  Internal_blno_row_view_t_ const_b_r(row_mapB.data(), row_mapB.extent(0));
  Internal_blno_nnz_view_t_ const_b_l(entriesB.data(), entriesB.extent(0));
  Internal_clno_row_view_t_ c_r(row_mapC.data(), row_mapC.extent(0));

  // Verify that graphs A and B are sorted.
  // This test is designed to be as efficient as possible, but still skip
  // it in a release build.
  //
  // Temporary fix for Trilinos issue #11655: Only perform this check if a TPL
  // is to be called. The KokkosKernels (non-TPL) implementation does not
  // actually require sorted indices yet. And Tpetra uses size_type = size_t, so
  // it will (currently) not be calling a TPL path.
#ifndef NDEBUG
  if constexpr (KokkosSparse::Impl::spgemm_symbolic_tpl_spec_avail<
                    const_handle_type, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_blno_row_view_t_,
                    Internal_blno_nnz_view_t_, Internal_clno_row_view_t_>::value) {
    if (!KokkosSparse::Impl::isCrsGraphSorted(const_a_r, const_a_l))
      throw std::runtime_error(
          "KokkosSparse::spgemm_symbolic: entries of A are not sorted within "
          "rows. May use KokkosSparse::sort_crs_matrix to sort it.");
    if (!KokkosSparse::Impl::isCrsGraphSorted(const_b_r, const_b_l))
      throw std::runtime_error(
          "KokkosSparse::spgemm_symbolic: entries of B are not sorted within "
          "rows. May use KokkosSparse::sort_crs_matrix to sort it.");
  }
#endif

  auto spgemmHandle = tmp_handle.get_spgemm_handle();

  if (!spgemmHandle) {
    throw std::invalid_argument(
        "KokkosSparse::spgemm_symbolic: the given KernelHandle does not have "
        "an SpGEMM handle associated with it.");
  }

  if (!spgemmHandle->checkMatrixIdentitiesSymbolic(const_a_r, const_a_l, const_b_r, const_b_l)) {
    throw std::invalid_argument(
        "KokkosSparse::spgemm_symbolic: once used, an spgemm handle cannot be "
        "reused for a product with a different sparsity pattern.\n"
        "The rowptrs and entries of A and B must be identical to those "
        "passed to the first spgemm_symbolic call.");
  }

  auto algo = spgemmHandle->get_algorithm_type();

  if (algo == SPGEMM_DEBUG || algo == SPGEMM_SERIAL) {
    // Never call a TPL if serial/debug is requested (this is needed for
    // testing)
    KokkosSparse::Impl::SPGEMM_SYMBOLIC<const_handle_type,  // KernelHandle,
                                        Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_blno_row_view_t_,
                                        Internal_blno_nnz_view_t_, Internal_clno_row_view_t_,
                                        false>::spgemm_symbolic(&tmp_handle,  // handle,
                                                                m, n, k, const_a_r, const_a_l, transposeA, const_b_r,
                                                                const_b_l, transposeB, c_r, computeRowptrs);
  } else {
    KokkosSparse::Impl::SPGEMM_SYMBOLIC<const_handle_type,  // KernelHandle,
                                        Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_blno_row_view_t_,
                                        Internal_blno_nnz_view_t_,
                                        Internal_clno_row_view_t_>::spgemm_symbolic(&tmp_handle,  // handle,
                                                                                    m, n, k, const_a_r, const_a_l,
                                                                                    transposeA, const_b_r, const_b_l,
                                                                                    transposeB, c_r, computeRowptrs);
  }
}

}  // namespace Experimental
}  // namespace KokkosSparse
#endif
