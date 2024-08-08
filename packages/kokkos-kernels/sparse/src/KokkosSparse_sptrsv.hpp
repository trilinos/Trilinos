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

/// \file KokkosSparse_sptrsv.hpp
/// \brief Parallel sparse triangular solve
///
/// This file provides KokkosSparse::sptrsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPTRSV_HPP_
#define KOKKOSSPARSE_SPTRSV_HPP_

#include <type_traits>

// #include "KokkosSparse_sptrsv_handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_sptrsv_symbolic_spec.hpp"
#include "KokkosSparse_sptrsv_solve_spec.hpp"

#include "KokkosSparse_sptrsv_cuSPARSE_impl.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_SPTRSV_SAME_TYPE(A, B) \
  std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

/**
 * @brief sptrsv symbolic phase for linear system Ax=b
 *
 * @tparam ExecutionSpace This kernels execution space type
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @param space The execution space instance this kernel will run on
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 */
template <typename ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
void sptrsv_symbolic(const ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap,
                     lno_nnz_view_t_ entries) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
                "sptrsv_symbolic: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
                "sptrsv_symbolic: A entry type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");

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
                       typename lno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      RowMap_Internal;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      Entries_Internal;

#ifdef KK_TRISOLVE_TIMERS
  Kokkos::Timer timer_sptrsv;
#endif
  RowMap_Internal rowmap_i   = rowmap;
  Entries_Internal entries_i = entries;

  KokkosSparse::Impl::SPTRSV_SYMBOLIC<ExecutionSpace, const_handle_type, RowMap_Internal,
                                      Entries_Internal>::sptrsv_symbolic(space, &tmp_handle, rowmap_i, entries_i);

#ifdef KK_TRISOLVE_TIMERS
  std::cout << "     > sptrsv_symbolic time = " << timer_sptrsv.seconds() << std::endl;
#endif
}  // sptrsv_symbolic

/**
 * @brief sptrsv symbolic phase for linear system Ax=b
 *
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 */
template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
void sptrsv_symbolic(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries) {
  using ExecutionSpace = typename KernelHandle::HandleExecSpace;
  auto my_exec_space   = ExecutionSpace();
  sptrsv_symbolic(my_exec_space, handle, rowmap, entries);
}

/**
 * @brief sptrsv symbolic phase for linear system Ax=b
 *
 * @tparam ExecutionSpace This kernels execution space type
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @param space The execution space instance this kernel will run on
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 * @param values The CRS matrix's (A) values
 */
template <typename ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_,
          typename scalar_nnz_view_t_>
void sptrsv_symbolic(ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries,
                     scalar_nnz_view_t_ values) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_type;

  static_assert(std::is_same_v<ExecutionSpace, typename KernelHandle::HandleExecSpace>,
                "sptrsv_symbolic: ExecutionSpace and HandleExecSpace need to match!");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
                "sptrsv_symbolic: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
                "sptrsv_symbolic: A entry type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
                "sptrsv_symbolic: A scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

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

#ifdef KK_TRISOLVE_TIMERS
  Kokkos::Timer timer_sptrsv;
#endif
  auto sptrsv_handle = handle->get_sptrsv_handle();
  if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
      using RowMap_Internal =
          Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                       typename lno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

      using Entries_Internal =
          Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

      using Values_Internal =
          Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

      RowMap_Internal rowmap_i   = rowmap;
      Entries_Internal entries_i = entries;
      Values_Internal values_i   = values;

      typedef typename KernelHandle::SPTRSVHandleType sptrsvHandleType;
      sptrsvHandleType *sh = handle->get_sptrsv_handle();
      auto nrows           = sh->get_nrows();

      KokkosSparse::Impl::sptrsvcuSPARSE_symbolic<ExecutionSpace, sptrsvHandleType, RowMap_Internal, Entries_Internal,
                                                  Values_Internal>(space, sh, nrows, rowmap_i, entries_i, values_i,
                                                                   false);
    } else {
      (void)values;
      KokkosSparse::Experimental::sptrsv_symbolic(space, handle, rowmap, entries);
    }

#else  // We better go to the native implementation
    (void)values;
    KokkosSparse::Experimental::sptrsv_symbolic(space, handle, rowmap, entries);
#endif
  } else {
    (void)values;
    KokkosSparse::Experimental::sptrsv_symbolic(space, handle, rowmap, entries);
  }
#ifdef KK_TRISOLVE_TIMERS
  std::cout << "     + sptrsv_symbolic time = " << timer_sptrsv.seconds() << std::endl;
#endif
}  // sptrsv_symbolic

/**
 * @brief sptrsv symbolic phase for linear system Ax=b
 *
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 * @param values The CRS matrix's (A) values
 */
template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
void sptrsv_symbolic(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values) {
  using ExecutionSpace = typename KernelHandle::HandleExecSpace;
  auto my_exec_space   = ExecutionSpace();

  sptrsv_symbolic(my_exec_space, handle, rowmap, entries, values);
}

/**
 * @brief sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam ExecutionSpace This kernels execution space
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @tparam scalar_nnz_view_t_ The CRS matrix's (A) values type
 * @tparam BType The b vector type
 * @tparam XType The x vector type
 * @param space The execution space instance this kernel will be run on
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 * @param values The CRS matrix's (A) values
 * @param b The b vector
 * @param x The x vector
 */
template <typename ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_,
          typename scalar_nnz_view_t_, class BType, class XType>
void sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries,
                  scalar_nnz_view_t_ values, BType b, XType x) {
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_lno_t ordinal_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_type;

  static_assert(std::is_same_v<ExecutionSpace, typename KernelHandle::HandleExecSpace>,
                "sptrsv solve: ExecutionSpace and HandleExecSpace need to match");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
                "sptrsv_solve: A size_type must match KernelHandle size_type "
                "(const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
                "sptrsv_solve: A entry type must match KernelHandle entry type (aka "
                "nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
                "sptrsv_solve: A scalar type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<BType>::value, "sptrsv: b is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XType>::value, "sptrsv: x is not a Kokkos::View.");
  static_assert((int)BType::rank == (int)XType::rank, "sptrsv: The ranks of b and x do not match.");
  static_assert(BType::rank == 1, "sptrsv: b and x must both either have rank 1.");
  static_assert(std::is_same<typename XType::value_type, typename XType::non_const_value_type>::value,
                "sptrsv: The output x must be nonconst.");
  static_assert(std::is_same<typename BType::device_type, typename XType::device_type>::value,
                "sptrsv: Views BType and XType have different device_types.");
  static_assert(std::is_same<typename BType::device_type::execution_space,
                             typename KernelHandle::SPTRSVHandleType::execution_space>::value,
                "sptrsv: KernelHandle and Views have different execution spaces.");
  static_assert(std::is_same<typename lno_row_view_t_::device_type, typename lno_nnz_view_t_::device_type>::value,
                "sptrsv: rowmap and entries have different device types.");
  static_assert(std::is_same<typename lno_row_view_t_::device_type, typename scalar_nnz_view_t_::device_type>::value,
                "sptrsv: rowmap and values have different device types.");

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
                       typename lno_row_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      RowMap_Internal;

  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                       typename lno_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      Entries_Internal;

  typedef Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                       typename scalar_nnz_view_t_::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      Values_Internal;

  typedef Kokkos::View<typename BType::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<BType>::array_layout, typename BType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      BType_Internal;

  typedef Kokkos::View<typename XType::non_const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<XType>::array_layout, typename XType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XType_Internal;

  RowMap_Internal rowmap_i   = rowmap;
  Entries_Internal entries_i = entries;
  Values_Internal values_i   = values;

  BType_Internal b_i = b;
  XType_Internal x_i = x;

  auto sptrsv_handle = handle->get_sptrsv_handle();
  if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
      typedef typename KernelHandle::SPTRSVHandleType sptrsvHandleType;
      sptrsvHandleType *sh = handle->get_sptrsv_handle();
      auto nrows           = sh->get_nrows();

      KokkosSparse::Impl::sptrsvcuSPARSE_solve<ExecutionSpace, sptrsvHandleType, RowMap_Internal, Entries_Internal,
                                               Values_Internal, BType_Internal, XType_Internal>(
          space, sh, nrows, rowmap_i, entries_i, values_i, b_i, x_i, false);
    } else {
      KokkosSparse::Impl::SPTRSV_SOLVE<ExecutionSpace, const_handle_type, RowMap_Internal, Entries_Internal,
                                       Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve(space,
                                                                                                      &tmp_handle,
                                                                                                      rowmap_i,
                                                                                                      entries_i,
                                                                                                      values_i, b_i,
                                                                                                      x_i);
    }
#else
    KokkosSparse::Impl::SPTRSV_SOLVE<ExecutionSpace, const_handle_type, RowMap_Internal, Entries_Internal,
                                     Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve(space, &tmp_handle,
                                                                                                    rowmap_i, entries_i,
                                                                                                    values_i, b_i, x_i);
#endif
  } else {
    KokkosSparse::Impl::SPTRSV_SOLVE<ExecutionSpace, const_handle_type, RowMap_Internal, Entries_Internal,
                                     Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve(space, &tmp_handle,
                                                                                                    rowmap_i, entries_i,
                                                                                                    values_i, b_i, x_i);
  }

}  // sptrsv_solve

/**
 * @brief sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam lno_row_view_t_ The CRS matrix's (A) rowmap type
 * @tparam lno_nnz_view_t_ The CRS matrix's (A) entries type
 * @tparam scalar_nnz_view_t_ The CRS matrix's (A) values type
 * @tparam BType The b vector type
 * @tparam XType The x vector type
 * @param handle KernelHandle instance
 * @param rowmap The CRS matrix's (A) rowmap
 * @param entries The CRS matrix's (A) entries
 * @param values The CRS matrix's (A) values
 * @param b The b vector
 * @param x The x vector
 */
template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
          class BType, class XType>
void sptrsv_solve(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                  BType b, XType x) {
  using ExecutionSpace = typename KernelHandle::HandleExecSpace;
  auto my_exec_space   = ExecutionSpace();
  sptrsv_solve(my_exec_space, handle, rowmap, entries, values, b, x);
}

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) || defined(DOXY)
/**
 * @brief Supernodal sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam ExecutionSpace This kernels execution space
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam XType The x and b vector type
 * @param space The execution space instance this kernel will run on
 * @param handle KernelHandle instance
 * @param x The x vector
 * @param b The b vector
 */
template <typename ExecutionSpace, typename KernelHandle, class XType>
void sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, XType x, XType b) {
  auto crsmat  = handle->get_sptrsv_handle()->get_crsmat();
  auto values  = crsmat.values;
  auto graph   = crsmat.graph;
  auto row_map = graph.row_map;
  auto entries = graph.entries;

  if (!(handle->get_sptrsv_handle()->is_numeric_complete())) {
    std::cout << std::endl
              << " ** needs to call sptrsv_compute before calling sptrsv_solve **" << std::endl
              << std::endl;
    return;
  }

  if (handle->is_sptrsv_lower_tri()) {
    // apply forward pivoting
    Kokkos::deep_copy(space, x, b);

    // the fifth argument (i.e., first x) is not used
    sptrsv_solve(space, handle, row_map, entries, values, x, x);
  } else {
    // the fifth argument (i.e., first x) is not used
    sptrsv_solve(space, handle, row_map, entries, values, b, b);

    // apply backward pivoting
    Kokkos::deep_copy(space, x, b);
  }
}

/**
 * @brief Supernodal sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam XType The x and b vector type
 * @param handle KernelHandle instance
 * @param x The x vector
 * @param b The b vector
 */
template <typename KernelHandle, class XType>
void sptrsv_solve(KernelHandle *handle, XType x, XType b) {
  using ExecutionSpace = typename KernelHandle::HandleExecSpace;
  auto my_exec_space   = ExecutionSpace();
  sptrsv_solve(my_exec_space, handle, x, b);
}

/**
 * @brief Supernodal sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam ExecutionSpace This kernels execution space
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam XType The x and b vector type
 * @param space The execution space instance this kernel will run on
 * @param handleL KernelHandle instance for lower triangular matrix
 * @param handleU KernelHandle instance for upper triangular matrix
 * @param x The x vector
 * @param b The b vector
 */
template <typename ExecutionSpace, typename KernelHandle, class XType>
void sptrsv_solve(ExecutionSpace &space, KernelHandle *handleL, KernelHandle *handleU, XType x, XType b) {
  // Lower-triangular solve
  sptrsv_solve(space, handleL, x, b);

  // copy the solution to rhs
  Kokkos::deep_copy(space, b, x);

  // uper-triangular solve
  sptrsv_solve(space, handleU, x, b);
}

/**
 * @brief Supernodal sptrsv solve phase of x for linear system Ax=b
 *
 * @tparam KernelHandle A specialization of
 * KokkosKernels::Experimental::KokkosKernelsHandle
 * @tparam XType The x and b vector type
 * @param handleL KernelHandle instance for lower triangular matrix
 * @param handleU KernelHandle instance for upper triangular matrix
 * @param x The x vector
 * @param b The b vector
 */
template <typename KernelHandle, class XType>
void sptrsv_solve(KernelHandle *handleL, KernelHandle *handleU, XType x, XType b) {
  using ExecutionSpace = typename KernelHandle::HandleExecSpace;
  auto my_exec_space   = ExecutionSpace();
  sptrsv_solve(my_exec_space, handleL, handleU, x, b);
}
#endif

template <class ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_,
          typename scalar_nnz_view_t_, class BType, class XType>
void sptrsv_solve_streams(const std::vector<ExecutionSpace> &execspace_v, const std::vector<KernelHandle *> &handle_v,
                          const std::vector<lno_row_view_t_> &rowmap_v, const std::vector<lno_nnz_view_t_> &entries_v,
                          const std::vector<scalar_nnz_view_t_> &values_v, const std::vector<BType> &b_v,
                          std::vector<XType> &x_v) {
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;
  using scalar_type  = typename KernelHandle::nnz_scalar_t;

  static_assert(Kokkos::is_execution_space<ExecutionSpace>::value, "ExecutionSpace is not valid");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename lno_row_view_t_::memory_space>::accessible,
                "sptrsv_solve_streams: ExecutionSpace cannot access data in "
                "lno_row_view_t_");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename lno_nnz_view_t_::memory_space>::accessible,
                "sptrsv_solve_streams: ExecutionSpace cannot access data in "
                "lno_nnz_view_t_");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename scalar_nnz_view_t_::memory_space>::accessible,
                "sptrsv_solve_streams: ExecutionSpace cannot access data in "
                "scalar_nnz_view_t_");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename BType::memory_space>::accessible,
                "sptrsv_solve_streams: ExecutionSpace cannot access data in BType");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XType::memory_space>::accessible,
                "sptrsv_solve_streams: ExecutionSpace cannot access data in XType");

  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
                "sptrsv_solve_streams: A size_type must match KernelHandle "
                "size_type (const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
                "sptrsv_solve_streams: A entry type must match KernelHandle entry type "
                "(aka nnz_lno_t, and const doesn't matter)");
  static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
                "sptrsv_solve_streams: A scalar type must match KernelHandle "
                "entry type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(Kokkos::is_view<BType>::value, "sptrsv_solve_streams: b is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XType>::value, "sptrsv_solve_streams: x is not a Kokkos::View.");
  static_assert((int)BType::rank == (int)XType::rank, "sptrsv_solve_streams: The ranks of b and x do not match.");
  static_assert(BType::rank == 1, "sptrsv_solve_streams: b and x must both either have rank 1.");
  static_assert(std::is_same<typename XType::value_type, typename XType::non_const_value_type>::value,
                "sptrsv_solve_streams: The output x must be nonconst.");
  static_assert(std::is_same<typename BType::device_type, typename XType::device_type>::value,
                "sptrsv_solve_streams: Views BType and XType have different "
                "device_types.");
  static_assert(std::is_same<ExecutionSpace, typename KernelHandle::SPTRSVHandleType::execution_space>::value,
                "sptrsv_solve_streams: KernelHandle's execution space is different from "
                "ExecutionSpace.");
  static_assert(std::is_same<typename BType::device_type::execution_space,
                             typename KernelHandle::SPTRSVHandleType::execution_space>::value,
                "sptrsv_solve_streams: KernelHandle and Views have different execution "
                "spaces.");
  static_assert(std::is_same<typename lno_row_view_t_::device_type, typename lno_nnz_view_t_::device_type>::value,
                "sptrsv_solve_streams: rowmap and entries have different device types.");
  static_assert(std::is_same<typename lno_row_view_t_::device_type, typename scalar_nnz_view_t_::device_type>::value,
                "sptrsv_solve_streams: rowmap and values have different device types.");

  // Check sizes of vectors
  if (execspace_v.size() != handle_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. handle_v.size() " << handle_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != rowmap_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. rowmap_v.size() " << rowmap_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != entries_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. entries_v.size() " << entries_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != values_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. values_v.size() " << values_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != b_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. b_v.size() " << b_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (execspace_v.size() != x_v.size()) {
    std::ostringstream os;
    os << "KokkosSparse::Experimental::sptrsv_solve_streams: vector sizes "
          "must match -- execspace_v.size() "
       << execspace_v.size() << " vs. x_v.size() " << x_v.size();
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using c_size_t    = typename KernelHandle::const_size_type;
  using c_lno_t     = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t  = typename KernelHandle::const_nnz_scalar_t;
  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  using const_handle_type = typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t,
                                                                                      c_exec_t, c_temp_t, c_persist_t>;

  using RowMap_Internal = Kokkos::View<typename lno_row_view_t_::const_value_type *,
                                       typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
                                       typename lno_row_view_t_::device_type,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using Entries_Internal = Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
                                        typename lno_nnz_view_t_::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using Values_Internal = Kokkos::View<typename scalar_nnz_view_t_::const_value_type *,
                                       typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
                                       typename scalar_nnz_view_t_::device_type,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using BType_Internal =
      Kokkos::View<typename BType::const_value_type *,
                   typename KokkosKernels::Impl::GetUnifiedLayout<BType>::array_layout, typename BType::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using XType_Internal = Kokkos::View<typename XType::non_const_value_type *,
                                      typename KokkosKernels::Impl::GetUnifiedLayout<XType>::array_layout,
                                      typename XType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  std::vector<const_handle_type> handle_i_v(execspace_v.size());
  std::vector<RowMap_Internal> rowmap_i_v(execspace_v.size());
  std::vector<Entries_Internal> entries_i_v(execspace_v.size());
  std::vector<Values_Internal> values_i_v(execspace_v.size());
  std::vector<BType_Internal> b_i_v(execspace_v.size());
  std::vector<XType_Internal> x_i_v(execspace_v.size());

  for (int i = 0; i < static_cast<int>(execspace_v.size()); i++) {
    handle_i_v[i]  = const_handle_type(*(handle_v[i]));
    rowmap_i_v[i]  = rowmap_v[i];
    entries_i_v[i] = entries_v[i];
    values_i_v[i]  = values_v[i];
    b_i_v[i]       = b_v[i];
    x_i_v[i]       = x_v[i];
  }

  if (handle_v[0]->get_sptrsv_handle()->get_algorithm() ==
      KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    // NOTE: assume all streams use the same SPTRSV_CUSPARSE algo.
    KokkosSparse::Impl::sptrsvcuSPARSE_solve_streams<ExecutionSpace, const_handle_type, RowMap_Internal,
                                                     Entries_Internal, Values_Internal, BType_Internal, XType_Internal>(
        execspace_v, handle_i_v, rowmap_i_v, entries_i_v, values_i_v, b_i_v, x_i_v, false);
#else
    KokkosSparse::Impl::SPTRSV_SOLVE<ExecutionSpace, const_handle_type, RowMap_Internal, Entries_Internal,
                                     Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve_streams(execspace_v,
                                                                                                            handle_i_v,
                                                                                                            rowmap_i_v,
                                                                                                            entries_i_v,
                                                                                                            values_i_v,
                                                                                                            b_i_v,
                                                                                                            x_i_v);
#endif
  } else {
    KokkosSparse::Impl::SPTRSV_SOLVE<ExecutionSpace, const_handle_type, RowMap_Internal, Entries_Internal,
                                     Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve_streams(execspace_v,
                                                                                                            handle_i_v,
                                                                                                            rowmap_i_v,
                                                                                                            entries_i_v,
                                                                                                            values_i_v,
                                                                                                            b_i_v,
                                                                                                            x_i_v);
  }

}  // sptrsv_solve_streams

}  // namespace Experimental
}  // namespace KokkosSparse

#undef KOKKOSKERNELS_SPTRSV_SAME_TYPE

#endif  // KOKKOSSPARSE_SPTRSV_HPP_
