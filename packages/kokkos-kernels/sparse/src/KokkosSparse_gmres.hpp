/*
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
*/

/// \file KokkosSparse_gmres.hpp
/// \brief GMRES Ax = b solver
///
/// This file provides KokkosSparse::gmres.  This function performs a
/// local (no MPI) solve of Ax = b for sparse A. It is expected that A is in
/// compressed row sparse ("Crs") format.
///
/// This algorithm is described in the paper:
/// GMRES - A Generalized Minimal Residual Algorithm for Solving Nonsymmetric
/// Linear Systems - Saad, Schultz
///
/// For more info, see example/gmres/README.md

#ifndef KOKKOSSPARSE_GMRES_HPP_
#define KOKKOSSPARSE_GMRES_HPP_

#include <type_traits>

#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"
#include "KokkosSparse_gmres_spec.hpp"
#include "KokkosSparse_Preconditioner.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_GMRES_SAME_TYPE(A, B) \
  std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

/// @brief
/// @tparam KernelHandle
/// @tparam AMatrix
/// @tparam BType
/// @tparam XType
/// @param handle
/// @param A
/// @param B
/// @param X
/// @param precond
template <typename KernelHandle, typename AMatrix, typename BType, typename XType>
void gmres(KernelHandle* handle, AMatrix& A, BType& B, XType& X, Preconditioner<AMatrix>* precond = nullptr) {
  using scalar_type  = typename KernelHandle::nnz_scalar_t;
  using size_type    = typename KernelHandle::size_type;
  using ordinal_type = typename KernelHandle::nnz_lno_t;

  static_assert(KOKKOSKERNELS_GMRES_SAME_TYPE(typename BType::value_type, scalar_type),
                "gmres: B scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_GMRES_SAME_TYPE(typename XType::value_type, scalar_type),
                "gmres: X scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_GMRES_SAME_TYPE(typename AMatrix::value_type, scalar_type),
                "gmres: A scalar type must match KernelHandle entry "
                "type (aka nnz_scalar_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_GMRES_SAME_TYPE(typename AMatrix::ordinal_type, ordinal_type),
                "gmres: A ordinal type must match KernelHandle entry "
                "type (aka nnz_lno_t, and const doesn't matter)");

  static_assert(KOKKOSKERNELS_GMRES_SAME_TYPE(typename AMatrix::size_type, size_type),
                "gmres: A size type must match KernelHandle entry "
                "type (aka size_type, and const doesn't matter)");

  static_assert(
      KokkosSparse::is_crs_matrix<AMatrix>::value || KokkosSparse::Experimental::is_bsr_matrix<AMatrix>::value,
      "gmres: A is not a CRS or BSR matrix.");
  static_assert(Kokkos::is_view<BType>::value, "gmres: B is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XType>::value, "gmres: X is not a Kokkos::View.");

  static_assert(BType::rank == 1, "gmres: B must have rank 1");
  static_assert(XType::rank == 1, "gmres: X must have rank 1");

  static_assert(std::is_same<typename XType::value_type, typename XType::non_const_value_type>::value,
                "gmres: The output X must be nonconst.");

  static_assert(std::is_same<typename XType::device_type, typename BType::device_type>::value,
                "gmres: X and B have different device types.");

  static_assert(std::is_same<typename AMatrix::device_type, typename BType::device_type>::value,
                "gmres: A and B have different device types.");

  using c_size_t   = typename KernelHandle::const_size_type;
  using c_lno_t    = typename KernelHandle::const_nnz_lno_t;
  using c_scalar_t = typename KernelHandle::const_nnz_scalar_t;

  using c_exec_t    = typename KernelHandle::HandleExecSpace;
  using c_temp_t    = typename KernelHandle::HandleTempMemorySpace;
  using c_persist_t = typename KernelHandle::HandlePersistentMemorySpace;

  if ((X.extent(0) != B.extent(0)) || (static_cast<size_t>(A.numPointCols()) != static_cast<size_t>(X.extent(0))) ||
      (static_cast<size_t>(A.numPointRows()) != static_cast<size_t>(B.extent(0)))) {
    std::ostringstream os;
    os << "KokkosSparse::gmres: Dimensions do not match: "
       << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << X.extent(0) << ", b: " << B.extent(0);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using const_handle_type = typename KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t,
                                                                                      c_exec_t, c_temp_t, c_persist_t>;

  const_handle_type tmp_handle(*handle);

  using AMatrix_Bsr_Internal =
      KokkosSparse::Experimental::BsrMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                                            typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                            typename AMatrix::const_size_type>;

  using AMatrix_Internal = std::conditional_t<
      KokkosSparse::is_crs_matrix<AMatrix>::value,
      KokkosSparse::CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                              typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              typename AMatrix::const_size_type>,
      AMatrix_Bsr_Internal>;

  using B_Internal =
      Kokkos::View<typename BType::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<BType>::array_layout, typename BType::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using X_Internal =
      Kokkos::View<typename XType::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayout<XType>::array_layout, typename XType::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

  using Precond_Internal = Preconditioner<AMatrix_Internal>;

  AMatrix_Internal A_i(A);
  B_Internal b_i = B;
  X_Internal x_i = X;

  Precond_Internal* precond_i = reinterpret_cast<Precond_Internal*>(precond);

  KokkosSparse::Impl::GMRES<const_handle_type, typename AMatrix_Internal::value_type,
                            typename AMatrix_Internal::ordinal_type, typename AMatrix_Internal::device_type,
                            typename AMatrix_Internal::memory_traits, typename AMatrix_Internal::size_type, B_Internal,
                            X_Internal>::gmres(&tmp_handle, A_i, b_i, x_i, precond_i);

}  // gmres

}  // namespace Experimental
}  // namespace KokkosSparse

#undef KOKKOSKERNELS_GMRES_SAME_TYPE

#endif  // KOKKOSSPARSE_GMRES_HPP_
