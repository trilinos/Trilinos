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
#ifndef _KOKKOS_SPGEMM_NUMERIC_HPP
#define _KOKKOS_SPGEMM_NUMERIC_HPP

#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_spgemm_numeric_spec.hpp"
#include "KokkosSparse_bspgemm_numeric_spec.hpp"

namespace KokkosSparse {

namespace Experimental {

//
// NOTE: block_dim = 1 for CRS-formated views
//       block_dim >= 1 for BSR-formatted views (bs=1 BSR is CRS)
//
// NOTE: Block CRS format is not yet supported !
//
template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename ascalar_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_,
          typename clno_nnz_view_t_, typename cscalar_nnz_view_t_>
void spgemm_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                    typename KernelHandle::const_nnz_lno_t n, typename KernelHandle::const_nnz_lno_t k,
                    alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, ascalar_nnz_view_t_ valuesA,

                    bool transposeA, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB, bscalar_nnz_view_t_ valuesB,
                    bool transposeB, clno_row_view_t_ row_mapC, clno_nnz_view_t_ &entriesC,
                    cscalar_nnz_view_t_ &valuesC,

                    typename KernelHandle::const_nnz_lno_t block_dim = 1) {
  static_assert(
      std::is_same<typename clno_nnz_view_t_::value_type, typename clno_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_numeric: Output matrix entriesView must be "
      "non-const.");

  static_assert(
      std::is_same<typename cscalar_nnz_view_t_::value_type, typename cscalar_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_numeric: Output matrix scalar view must be "
      "non-const.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename alno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: Size type of left handside matrix should "
      "be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename blno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: Size type of right handside matrix should "
      "be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename clno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: Size type of output matrix should be same "
      "as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename alno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: lno type of left handside matrix should "
      "be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename blno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: lno type of right handside matrix should "
      "be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename clno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: lno type of output matrix should be same "
      "as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename ascalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: scalar type of left handside matrix "
      "should be same as kernelHandle scalar.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename bscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: scalar type of right handside matrix "
      "should be same as kernelHandle scalar.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename cscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_numeric: scalar type of output matrix should be "
      "same as kernelHandle scalar.");

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
                       UniformDevice_t,  // typename
                                         // alno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout,
                       UniformDevice_t,  // typename
                                         // alno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename ascalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<ascalar_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //       typename ascalar_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  typedef Kokkos::View<typename blno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout,
                       UniformDevice_t,  //       typename
                                         //       blno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_row_view_t_;

  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //       typename
                                         //       blno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_nnz_view_t_;

  typedef Kokkos::View<typename bscalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<bscalar_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //       typename bscalar_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_bscalar_nnz_view_t_;

  // static assert clno_row_view_t_ can be const type (row map is fixed after
  // symbolic phase).
  typedef Kokkos::View<typename clno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout,
                       UniformDevice_t,  //       typename
                                         //       clno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_clno_row_view_t_;

  typedef Kokkos::View<typename clno_nnz_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //       typename
                                         //       clno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_clno_nnz_view_t_;

  typedef Kokkos::View<typename cscalar_nnz_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<cscalar_nnz_view_t_>::array_layout,
                       UniformDevice_t,  //       typename cscalar_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_cscalar_nnz_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_mapA.data(), row_mapA.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entriesA.data(), entriesA.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_s(valuesA.data(), valuesA.extent(0));
  Internal_blno_row_view_t_ const_b_r(row_mapB.data(), row_mapB.extent(0));
  Internal_blno_nnz_view_t_ const_b_l(entriesB.data(), entriesB.extent(0));
  Internal_bscalar_nnz_view_t_ const_b_s(valuesB.data(), valuesB.extent(0));
  Internal_clno_row_view_t_ const_c_r(row_mapC.data(), row_mapC.extent(0));
  Internal_clno_nnz_view_t_ nonconst_c_l(entriesC.data(), entriesC.extent(0));
  Internal_cscalar_nnz_view_t_ nonconst_c_s(valuesC.data(), valuesC.extent(0));

  if (block_dim > 1) {
    KokkosSparse::Impl::BSPGEMM_NUMERIC<
        const_handle_type, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_,
        Internal_blno_row_view_t_, Internal_blno_nnz_view_t_, Internal_bscalar_nnz_view_t_, Internal_clno_row_view_t_,
        Internal_clno_nnz_view_t_, Internal_cscalar_nnz_view_t_>::bspgemm_numeric(&tmp_handle, m, n, k, block_dim,
                                                                                  const_a_r, const_a_l, const_a_s,
                                                                                  transposeA, const_b_r, const_b_l,
                                                                                  const_b_s, transposeB, const_c_r,
                                                                                  nonconst_c_l, nonconst_c_s);
    return;
  }

  auto spgemmHandle = tmp_handle.get_spgemm_handle();

  if (!spgemmHandle) {
    throw std::invalid_argument(
        "KokkosSparse::spgemm_numeric: the given KernelHandle does not have "
        "an SpGEMM handle associated with it.");
  }

  if (!spgemmHandle->checkMatrixIdentitiesNumeric(const_a_r, const_a_l, const_b_r, const_b_l)) {
    throw std::invalid_argument(
        "KokkosSparse::spgemm_numeric: once used, an spgemm handle cannot be "
        "reused for a product with a different sparsity pattern.\n"
        "The rowptrs and entries of A and B must be identical to those "
        "passed to the first spgemm_symbolic and spgemm_numeric calls.");
  }

  auto algo = spgemmHandle->get_algorithm_type();

  if (algo == SPGEMM_DEBUG || algo == SPGEMM_SERIAL) {
    // Never call a TPL if serial/debug is requested (this is needed for
    // testing)
    KokkosSparse::Impl::SPGEMM_NUMERIC<
        const_handle_type,  // KernelHandle,
        Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_, Internal_blno_row_view_t_,
        Internal_blno_nnz_view_t_, Internal_bscalar_nnz_view_t_, Internal_clno_row_view_t_, Internal_clno_nnz_view_t_,
        Internal_cscalar_nnz_view_t_,
        false>::spgemm_numeric(&tmp_handle,  // handle,
                               m, n, k, const_a_r, const_a_l, const_a_s, transposeA, const_b_r, const_b_l, const_b_s,
                               transposeB, const_c_r, nonconst_c_l, nonconst_c_s);
  } else {
    KokkosSparse::Impl::SPGEMM_NUMERIC<
        const_handle_type,  // KernelHandle,
        Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_, Internal_blno_row_view_t_,
        Internal_blno_nnz_view_t_, Internal_bscalar_nnz_view_t_, Internal_clno_row_view_t_, Internal_clno_nnz_view_t_,
        Internal_cscalar_nnz_view_t_>::spgemm_numeric(&tmp_handle,  // handle,
                                                      m, n, k, const_a_r, const_a_l, const_a_s, transposeA, const_b_r,
                                                      const_b_l, const_b_s, transposeB, const_c_r, nonconst_c_l,
                                                      nonconst_c_s);
  }
}

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
