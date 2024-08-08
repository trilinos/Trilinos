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

/// \file KokkosLapack_gesv.hpp
/// \brief Local dense linear solve
///
/// This file provides KokkosLapack::gesv. This function performs a
/// local (no MPI) dense linear solve on a system of linear equations
/// A * X = B where A is a general N-by-N matrix and X and B are N-by-NRHS
/// matrices.

#ifndef KOKKOSLAPACK_GESV_HPP_
#define KOKKOSLAPACK_GESV_HPP_

#include <type_traits>

#include "KokkosLapack_gesv_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief Solve the dense linear equation system A*X = B.
///
/// \tparam ExecutionSpace the space where the kernel will run.
/// \tparam AMatrix Input matrix/Output LU, as a 2-D Kokkos::View.
/// \tparam BXMV Input (right-hand side)/Output (solution) (multi)vector, as a
/// 1-D or 2-D Kokkos::View.
/// \tparam IPIVV Output pivot indices, as a 1-D Kokkos::View
///
/// \param space [in] execution space instance used to specified how to execute
///   the gesv kernels.
/// \param A [in,out] On entry, the N-by-N matrix to be solved. On exit, the
/// factors L and U from
///   the factorization A = P*L*U; the unit diagonal elements of L are not
///   stored.
/// \param B [in,out] On entry, the right hand side (multi)vector B. On exit,
/// the solution (multi)vector X.
/// \param IPIV [out] On exit, the pivot indices (for partial pivoting).
/// If the View extents are zero and its data pointer is NULL, pivoting is not
/// used.
///
template <class ExecutionSpace, class AMatrix, class BXMV, class IPIVV>
void gesv(const ExecutionSpace& space, const AMatrix& A, const BXMV& B, const IPIVV& IPIV) {
  // NOTE: Currently, KokkosLapack::gesv only supports LAPACK, MAGMA and
  // rocSOLVER TPLs.
  //       MAGMA/rocSOLVER TPL should be enabled to call the MAGMA/rocSOLVER GPU
  //       interface for device views LAPACK TPL should be enabled to call the
  //       LAPACK interface for host views

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename BXMV::memory_space>::accessible);
#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  if constexpr (!std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
    static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename IPIVV::memory_space>::accessible);
  }
#else
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename IPIVV::memory_space>::accessible);
#endif
  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::gesv: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<BXMV>::value, "KokkosLapack::gesv: B must be a Kokkos::View.");
  static_assert(Kokkos::is_view<IPIVV>::value, "KokkosLapack::gesv: IPIV must be a Kokkos::View.");
  static_assert(static_cast<int>(AMatrix::rank) == 2, "KokkosLapack::gesv: A must have rank 2.");
  static_assert(static_cast<int>(BXMV::rank) == 1 || static_cast<int>(BXMV::rank) == 2,
                "KokkosLapack::gesv: B must have either rank 1 or rank 2.");
  static_assert(static_cast<int>(IPIVV::rank) == 1, "KokkosLapack::gesv: IPIV must have rank 1.");

  int64_t IPIV0 = IPIV.extent(0);
  int64_t A0    = A.extent(0);
  int64_t A1    = A.extent(1);
  int64_t B0    = B.extent(0);

  // Check validity of pivot argument
  bool valid_pivot = (IPIV0 == A1) || ((IPIV0 == 0) && (IPIV.data() == nullptr));
  if (!(valid_pivot)) {
    std::ostringstream os;
    os << "KokkosLapack::gesv: IPIV: " << IPIV0 << ". "
       << "Valid options include zero-extent 1-D view (no pivoting), or 1-D "
          "View with size of "
       << A0 << " (partial pivoting).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check for no pivoting case. Only MAGMA supports no pivoting interface
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA   // have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // and have LAPACK TPL
  if ((!std::is_same<typename AMatrix::device_type::memory_space, Kokkos::CudaSpace>::value) && (IPIV0 == 0) &&
      (IPIV.data() == nullptr)) {
    std::ostringstream os;
    os << "KokkosLapack::gesv: IPIV: " << IPIV0 << ". "
       << "LAPACK TPL does not support no pivoting.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
#endif
#else                                   // not have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // but have LAPACK TPL
  if ((IPIV0 == 0) && (IPIV.data() == nullptr)) {
    std::ostringstream os;
    os << "KokkosLapack::gesv: IPIV: " << IPIV0 << ". "
       << "LAPACK TPL does not support no pivoting.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
#endif
#endif

  // Check compatibility of dimensions at run time.
  if ((A0 < A1) || (A0 != B0)) {
    std::ostringstream os;
    os << "KokkosLapack::gesv: Dimensions of A, and B do not match: "
       << " A: " << A.extent(0) << " x " << A.extent(1) << " B: " << B.extent(0) << " x " << B.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  typedef Kokkos::View<typename AMatrix::non_const_value_type**, typename AMatrix::array_layout,
                       typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      AMatrix_Internal;
  typedef Kokkos::View<typename BXMV::non_const_value_type**, typename BXMV::array_layout, typename BXMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      BXMV_Internal;
  typedef Kokkos::View<typename IPIVV::non_const_value_type*, typename IPIVV::array_layout, typename IPIVV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      IPIVV_Internal;
  AMatrix_Internal A_i = A;
  // BXMV_Internal B_i = B;
  IPIVV_Internal IPIV_i = IPIV;

  if (BXMV::rank == 1) {
    auto B_i = BXMV_Internal(B.data(), B.extent(0), 1);
    KokkosLapack::Impl::GESV<ExecutionSpace, AMatrix_Internal, BXMV_Internal, IPIVV_Internal>::gesv(space, A_i, B_i,
                                                                                                    IPIV_i);
  } else {  // BXMV::rank == 2
    auto B_i = BXMV_Internal(B.data(), B.extent(0), B.extent(1));
    KokkosLapack::Impl::GESV<ExecutionSpace, AMatrix_Internal, BXMV_Internal, IPIVV_Internal>::gesv(space, A_i, B_i,
                                                                                                    IPIV_i);
  }
}

/// \brief Solve the dense linear equation system A*X = B.
///
/// \tparam AMatrix Input matrix/Output LU, as a 2-D Kokkos::View.
/// \tparam BXMV Input (right-hand side)/Output (solution) (multi)vector, as a
/// 1-D or 2-D Kokkos::View.
/// \tparam IPIVV Output pivot indices, as a 1-D Kokkos::View
///
/// \param A [in,out] On entry, the N-by-N matrix to be solved. On exit, the
/// factors L and U from
///   the factorization A = P*L*U; the unit diagonal elements of L are not
///   stored.
/// \param B [in,out] On entry, the right hand side (multi)vector B. On exit,
/// the solution (multi)vector X.
/// \param IPIV [out] On exit, the pivot indices (for partial pivoting).
/// If the View extents are zero and its data pointer is NULL, pivoting is not
/// used.
///
template <class AMatrix, class BXMV, class IPIVV>
void gesv(const AMatrix& A, const BXMV& B, const IPIVV& IPIV) {
  typename AMatrix::execution_space space{};
  gesv(space, A, B, IPIV);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GESV_HPP_
