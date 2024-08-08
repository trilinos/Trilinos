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

/// \file KokkosLapack_svd.hpp
/// \brief Singular Value Decomposition (SVD)
///
/// This file provides KokkosLapack::svd. This function performs a
/// local (no MPI) singular value decomposition of the input matrix A
/// and returns the singular values and vectors dedending on input flags.

#ifndef KOKKOSLAPACK_SVD_HPP_
#define KOKKOSLAPACK_SVD_HPP_

#include <type_traits>

#include "KokkosLapack_svd_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

// clang-format off
/// \brief Compute the Singular Value Decomposition of A = U*S*Vt
///
/// \tparam ExecutionSpace the space where the kernel will run.
/// \tparam AMatrix (mxn) matrix as a rank-2 Kokkos::View.
/// \tparam SVector min(m,n) vector as a rank-1 Kokkos::View
/// \tparam UMatrix (mxm) matrix as a rank-2 Kokkos::View
/// \tparam VMatrix (nxn) matrix as a rank-2 Kokkos::View
///
/// \param space [in] execution space instance used to specified how to execute
///   the svd kernels.
/// \param jobu [in] flag to control the computation of the left singular
/// vectors when set to: 'A' all vectors are computed, 'S' the first min(m,n)
/// singular vectors are computed, 'O' the first min(m,n) singular vectors are
/// overwritten into A, 'N' no singular vectors are computed.
/// \param jobvt [in] flag to control the computation of the right singular
/// vectors when set to: 'A' all vectors are computed, 'S' the first min(m,n)
/// singular vectors are computed, 'O' the first min(m,n) singular vectors are
/// overwritten into A, 'N' no singular vectors are computed.
/// \param A [in] An m-by-n matrix to be decomposed using its singular values.
/// \param S [out] Vector of the min(m, n) singular values of A.
/// \param U [out] the first min(m, n) columns of U are the left singular
/// vectors of A.
/// \param Vt [out] the first min(m, n) columns of Vt are the right singular
/// vectors of A.
///
// clang-format on
template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
void svd(const ExecutionSpace& space, const char jobu[], const char jobvt[], const AMatrix& A, const SVector& S,
         const UMatrix& U, const VMatrix& Vt) {
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename SVector::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename UMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename VMatrix::memory_space>::accessible);
  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::svd: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<SVector>::value, "KokkosLapack::svd: S must be a Kokkos::View.");
  static_assert(Kokkos::is_view<UMatrix>::value, "KokkosLapack::svd: U must be a Kokkos::View.");
  static_assert(Kokkos::is_view<VMatrix>::value, "KokkosLapack::svd: Vt must be a Kokkos::View.");
  static_assert(AMatrix::rank() == 2, "KokkosLapack::svd: A must have rank 2.");
  static_assert(SVector::rank() == 1, "KokkosLapack::svd: S must have rank 1.");
  static_assert(UMatrix::rank() == 2, "KokkosLapack::svd: U must have rank 2.");
  static_assert(VMatrix::rank() == 2, "KokkosLapack::svd: Vt must have rank 2.");

  int64_t m     = A.extent(0);
  int64_t n     = A.extent(1);
  int64_t rankA = Kokkos::min(m, n);

  // No work to do since the matrix is empty...
  // Also do not send a matrix with size zero
  // to Lapack TPLs or they will complain!
  if ((m == 0) || (n == 0)) {
    return;
  }

  // Check the jobu and jobvt control flags
  // The only valid options there are 'A', 'S', 'O' and 'N'
  const bool is_jobu_invalid = !((jobu[0] == 'A') || (jobu[0] == 'a') || (jobu[0] == 'S') || (jobu[0] == 's') ||
                                 (jobu[0] == 'O') || (jobu[0] == 'o') || (jobu[0] == 'N') || (jobu[0] == 'n'));

  const bool is_jobvt_invalid = !((jobvt[0] == 'A') || (jobvt[0] == 'a') || (jobvt[0] == 'S') || (jobvt[0] == 's') ||
                                  (jobvt[0] == 'O') || (jobvt[0] == 'o') || (jobvt[0] == 'N') || (jobvt[0] == 'n'));

  if (is_jobu_invalid && is_jobvt_invalid) {
    std::ostringstream oss;
    oss << "KokkosLapack::svd: both jobu and jobvt are invalid!\n"
        << "Possible values are A, S, O or N, submitted values are " << jobu[0] << " and " << jobvt[0] << "\n";
    KokkosKernels::Impl::throw_runtime_exception(oss.str());
  }
  if (is_jobu_invalid) {
    std::ostringstream oss;
    oss << "KokkosLapack::svd: jobu is invalid!\n"
        << "Possible values are A, S, O or N, submitted value is " << jobu[0] << "\n";
    KokkosKernels::Impl::throw_runtime_exception(oss.str());
  }
  if (is_jobvt_invalid) {
    std::ostringstream oss;
    oss << "KokkosLapack::svd: jobvt is invalid!\n"
        << "Possible values are A, S, O or N, submitted value is " << jobvt[0] << "\n";
    KokkosKernels::Impl::throw_runtime_exception(oss.str());
  }

  if (((jobu[0] == 'O') || (jobu[0] == 'o')) && ((jobvt[0] == 'O') || (jobvt[0] == 'o'))) {
    std::ostringstream oss;
    oss << "KokkosLapack::svd: jobu and jobvt cannot be O at the same time!\n";
    KokkosKernels::Impl::throw_runtime_exception(oss.str());
  }

  // Check validity of output views sizes
  // Note that of jobu/jobvt are set to O or N
  // then the associated matrix does not need storage
  bool is_extent_invalid = false;
  std::ostringstream os;
  if (S.extent_int(0) != rankA) {
    is_extent_invalid = true;
    os << "KokkosLapack::svd: S has extent " << S.extent(0) << ", instead of " << rankA << ".\n";
  }
  if ((jobu[0] == 'A') || (jobu[0] == 'a') || (jobu[0] == 'S') || (jobu[0] == 's')) {
    if (U.extent_int(0) != m || U.extent_int(1) != m) {
      is_extent_invalid = true;
      os << "KokkosLapack::svd: U has extents (" << U.extent(0) << ", " << U.extent(1) << ") instead of (" << m << ", "
         << m << ").\n";
    }
  }
  if ((jobvt[0] == 'A') || (jobvt[0] == 'a') || (jobvt[0] == 'S') || (jobvt[0] == 's')) {
    if (Vt.extent_int(0) != n || Vt.extent_int(1) != n) {
      is_extent_invalid = true;
      os << "KokkosLapack::svd: V has extents (" << Vt.extent(0) << ", " << Vt.extent(1) << ") instead of (" << n
         << ", " << n << ").\n";
    }
  }
  if (is_extent_invalid) {
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
  if (std::is_same_v<ExecutionSpace, Kokkos::Cuda> && (A.extent(0) < A.extent(1))) {
    throw std::runtime_error(
        "CUSOLVER does not support SVD for matrices with more columns "
        "than rows, you can transpose you matrix first then compute "
        "SVD of that transpose: At=VSUt, and swap the output U and Vt"
        " and transpose them to recover the desired SVD.");
  }
#endif

  using AMatrix_Internal = Kokkos::View<typename AMatrix::non_const_value_type**, typename AMatrix::array_layout,
                                        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using SVector_Internal = Kokkos::View<typename SVector::non_const_value_type*, typename SVector::array_layout,
                                        typename SVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using UMatrix_Internal = Kokkos::View<typename UMatrix::non_const_value_type**, typename UMatrix::array_layout,
                                        typename UMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using VMatrix_Internal = Kokkos::View<typename VMatrix::non_const_value_type**, typename VMatrix::array_layout,
                                        typename VMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  AMatrix_Internal A_i  = A;
  SVector_Internal S_i  = S;
  UMatrix_Internal U_i  = U;
  VMatrix_Internal Vt_i = Vt;

  KokkosLapack::Impl::SVD<ExecutionSpace, AMatrix_Internal, SVector_Internal, UMatrix_Internal, VMatrix_Internal>::svd(
      space, jobu, jobvt, A_i, S_i, U_i, Vt_i);
}

// clang-format off
/// \brief Compute the Singular Value Decomposition of A = U*S*Vt
///
/// \tparam AMatrix (mxn) matrix as a rank-2 Kokkos::View.
/// \tparam SVector min(m,n) vector as a rank-1 Kokkos::View
/// \tparam UMatrix (mxm) matrix as a rank-2 Kokkos::View
/// \tparam VMatrix (nxn) matrix as a rank-2 Kokkos::View
///
/// \param jobu [in] flag to control the computation of the left singular
/// vectors when set to: 'A' all vectors are computed, 'S' the first min(m,n)
/// singular vectors are computed, 'O' the first min(m,n) singular vectors are
/// overwritten into A, 'N' no singular vectors are computed.
/// \param jobvt [in] flag to control the computation of the right singular
/// vectors when set to: 'A' all vectors are computed, 'S' the first min(m,n)
/// singular vectors are computed, 'O' the first min(m,n) singular vectors are
/// overwritten into A, 'N' no singular vectors are computed.
/// \param A [in] An m-by-n matrix to be decomposed using its singular values.
/// \param S [out] Vector of the min(m, n) singular values of A.
/// \param U [out] the first min(m, n) columns of U are the left singular
/// vectors of A.
/// \param Vt [out] the first min(m, n) columns of Vt are the right singular
/// vectors of A.
///
// clang-format on
template <class AMatrix, class SVector, class UMatrix, class VMatrix>
void svd(const char jobu[], const char jobvt[], const AMatrix& A, const SVector& S, const UMatrix& U,
         const VMatrix& Vt) {
  typename AMatrix::execution_space space{};
  svd(space, jobu, jobvt, A, S, U, Vt);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_SVD_HPP_
