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
#ifndef _KOKKOS_SPGEMM_JACOBI_HPP
#define _KOKKOS_SPGEMM_JACOBI_HPP

#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_spgemm_jacobi_spec.hpp"

namespace KokkosSparse {

namespace Experimental {

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename ascalar_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_,
          typename clno_nnz_view_t_, typename cscalar_nnz_view_t_, typename dinv_view_t_>
void spgemm_jacobi(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                   typename KernelHandle::const_nnz_lno_t n, typename KernelHandle::const_nnz_lno_t k,

                   alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, ascalar_nnz_view_t_ valuesA,

                   bool transposeA, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB, bscalar_nnz_view_t_ valuesB,

                   bool transposeB, clno_row_view_t_ row_mapC, clno_nnz_view_t_ &entriesC, cscalar_nnz_view_t_ &valuesC,

                   typename cscalar_nnz_view_t_::const_value_type omega, dinv_view_t_ dinv) {
  static_assert(
      std::is_same<typename clno_row_view_t_::value_type, typename clno_row_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Output matrix rowmap must be non-const.");

  static_assert(
      std::is_same<typename clno_nnz_view_t_::value_type, typename clno_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Output matrix entriesView must be "
      "non-const.");

  static_assert(
      std::is_same<typename cscalar_nnz_view_t_::value_type, typename cscalar_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Output matrix scalar view must be "
      "non-const.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename alno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Size type of left handside matrix should "
      "be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename blno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Size type of right handside matrix should "
      "be same as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_size_type, typename clno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: Size type of output matrix should be same "
      "as kernelHandle sizetype.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename alno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: lno type of left handside matrix should be "
      "same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename blno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: lno type of right handside matrix should "
      "be same as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_lno_t, typename clno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: lno type of output matrix should be same "
      "as kernelHandle lno_t.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename ascalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: scalar type of left handside matrix should "
      "be same as kernelHandle scalar.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename bscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: scalar type of right handside matrix "
      "should be same as kernelHandle scalar.");

  static_assert(
      std::is_same<typename KernelHandle::const_nnz_scalar_t, typename cscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_jacobi: scalar type of output matrix should be "
      "same as kernelHandle scalar.");

  if (transposeA || transposeB) {
    throw std::runtime_error(
        "spgemm-jacobi does not support transposed multiply. "
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
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_row_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_row_view_t_;

  typedef Kokkos::View<typename alno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_alno_nnz_view_t_;

  typedef Kokkos::View<typename ascalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<ascalar_nnz_view_t_>::array_layout,
                       UniformDevice_t, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_ascalar_nnz_view_t_;

  typedef Kokkos::View<typename blno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_row_view_t_;

  typedef Kokkos::View<typename blno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_blno_nnz_view_t_;

  typedef Kokkos::View<typename bscalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<bscalar_nnz_view_t_>::array_layout,
                       UniformDevice_t, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_bscalar_nnz_view_t_;

  typedef Kokkos::View<typename clno_row_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_clno_row_view_t_;

  typedef Kokkos::View<typename clno_nnz_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<clno_nnz_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_clno_nnz_view_t_;

  typedef Kokkos::View<typename cscalar_nnz_view_t_::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<cscalar_nnz_view_t_>::array_layout,
                       UniformDevice_t, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_cscalar_nnz_view_t_;

  typedef Kokkos::View<typename dinv_view_t_::const_value_type **,
                       typename KokkosKernels::Impl::GetUnifiedLayout<dinv_view_t_>::array_layout, UniformDevice_t,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_dinv_view_t_;

  Internal_alno_row_view_t_ const_a_r(row_mapA.data(), row_mapA.extent(0));
  Internal_alno_nnz_view_t_ const_a_l(entriesA.data(), entriesA.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_s(valuesA.data(), valuesA.extent(0));
  Internal_blno_row_view_t_ const_b_r(row_mapB.data(), row_mapB.extent(0));
  Internal_blno_nnz_view_t_ const_b_l(entriesB.data(), entriesB.extent(0));
  Internal_bscalar_nnz_view_t_ const_b_s(valuesB.data(), valuesB.extent(0));
  Internal_clno_row_view_t_ nonconst_c_r(row_mapC.data(), row_mapC.extent(0));
  Internal_clno_nnz_view_t_ nonconst_c_l(entriesC.data(), entriesC.extent(0));
  Internal_cscalar_nnz_view_t_ nonconst_c_s(valuesC.data(), valuesC.extent(0));
  Internal_dinv_view_t_ const_d_s(dinv.data(), dinv.extent(0), dinv.extent(1));

  KokkosSparse::Impl::SPGEMM_JACOBI<const_handle_type, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_,
                                    Internal_ascalar_nnz_view_t_, Internal_blno_row_view_t_, Internal_blno_nnz_view_t_,
                                    Internal_bscalar_nnz_view_t_, Internal_clno_row_view_t_, Internal_clno_nnz_view_t_,
                                    Internal_cscalar_nnz_view_t_,
                                    Internal_dinv_view_t_>::spgemm_jacobi(&tmp_handle, m, n, k, const_a_r, const_a_l,
                                                                          const_a_s, transposeA, const_b_r, const_b_l,
                                                                          const_b_s, transposeB, nonconst_c_r,
                                                                          nonconst_c_l, nonconst_c_s, omega, const_d_s);
}

}  // namespace Experimental
}  // namespace KokkosSparse
#endif
