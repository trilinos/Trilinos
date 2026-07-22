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

/// \file KokkosLapack_geqrf.hpp
/// \brief QR factorization
///
/// This file provides KokkosLapack::geqrf. This function performs a
/// local (no MPI) QR factorization of a M-by-N matrix A.

#ifndef KOKKOSLAPACK_GEQRF_HPP_
#define KOKKOSLAPACK_GEQRF_HPP_

#include <type_traits>

#include "KokkosLapack_geqrf_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief Computes a QR factorization of a matrix A
///
/// \tparam ExecutionSpace The space where the kernel will run.
/// \tparam AMatrix        Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray       Type of array Tau, as a 1-D Kokkos::View.
/// \tparam InfoArray      Type of array Info, as a 1-D Kokkos::View.
///
/// \param space [in] Execution space instance used to specified how to execute
///                   the geqrf kernels.
/// \param A [in,out] On entry, the M-by-N matrix to be factorized.
///                   On exit, the elements on and above the diagonal contain
///                   the min(M,N)-by-N upper trapezoidal matrix R (R is upper
///                   triangular if M >= N); the elements below the diagonal,
///                   with the array Tau, represent the unitary matrix Q as a
///                   product of min(M,N) elementary reflectors. The matrix Q
///                   is represented as a product of elementary reflectors
///                     Q = H(1) H(2) . . . H(k), where k = min(M,N).
///                   Each H(i) has the form
///                     H(i) = I - Tau(i) * v * v**H,
///                   where v is a vector with v(1:i-1) = 0 and v(i) = 1;
///                   v(i+1:M) is stored on exit in A(i+1:M,i).
/// \param Tau [out]  One-dimensional array of size min(M,N) that contains the
///                   scalar factors of the elementary reflectors.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
void geqrf(const ExecutionSpace& space, const AMatrix& A, const TauArray& Tau, const InfoArray& Info) {
  // NOTE: Currently, KokkosLapack::geqrf only supports LAPACK, MAGMA and
  // rocSOLVER TPLs.
  //       MAGMA/rocSOLVER TPL should be enabled to call the MAGMA/rocSOLVER GPU
  //       interface for device views LAPACK TPL should be enabled to call the
  //       LAPACK interface for host views

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename TauArray::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename InfoArray::memory_space>::accessible);

  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::geqrf: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<TauArray>::value, "KokkosLapack::geqrf: Tau must be Kokkos::View.");
  static_assert(Kokkos::is_view<InfoArray>::value, "KokkosLapack::geqrf: Info must be Kokkos::View.");

  static_assert(static_cast<int>(AMatrix::rank) == 2, "KokkosLapack::geqrf: A must have rank 2.");
  static_assert(static_cast<int>(TauArray::rank) == 1, "KokkosLapack::geqrf: Tau must have rank 1.");
  static_assert(static_cast<int>(InfoArray::rank) == 1, "KokkosLapack::geqrf: Info must have rank 1.");

  static_assert(std::is_same_v<typename InfoArray::non_const_value_type, int>,
                "KokkosLapack::geqrf: Info must be an array of integers.");

  const int64_t m     = A.extent(0);
  const int64_t n     = A.extent(1);
  const int64_t tau0  = Tau.extent(0);
  const int64_t info0 = Info.extent(0);

  // Check validity of dimensions
  if (tau0 != std::min(m, n)) {
    std::ostringstream os;
    os << "KokkosLapack::geqrf: length of Tau must be equal to min(m,n): "
       << " A: " << m << " x " << n << ", Tau length = " << tau0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (info0 < 1) {
    std::ostringstream os;
    os << "KokkosLapack::geqrf: length of Info must be at least 1, Info length = " << info0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using AMatrix_Internal   = Kokkos::View<typename AMatrix::non_const_value_type**, typename AMatrix::array_layout,
                                        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using TauArray_Internal  = Kokkos::View<typename TauArray::non_const_value_type*, typename TauArray::array_layout,
                                         typename TauArray::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InfoArray_Internal = Kokkos::View<typename InfoArray::non_const_value_type*, typename InfoArray::array_layout,
                                          typename InfoArray::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  AMatrix_Internal A_i      = A;
  TauArray_Internal Tau_i   = Tau;
  InfoArray_Internal Info_i = Info;

  KokkosLapack::Impl::GEQRF<ExecutionSpace, AMatrix_Internal, TauArray_Internal, InfoArray_Internal>::geqrf(
      space, A_i, Tau_i, Info_i);
}

/// \brief Computes a QR factorization of a matrix A
///
/// \tparam AMatrix   Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray  Type of array Tau, as a 1-D Kokkos::View.
/// \tparam InfoArray Type of array Info, as a 1-D Kokkos::View.
///
/// \param A [in,out] On entry, the M-by-N matrix to be factorized.
///                   On exit, the elements on and above the diagonal contain
///                   the min(M,N)-by-N upper trapezoidal matrix R (R is upper
///                   triangular if M >= N); the elements below the diagonal,
///                   with the array Tau, represent the unitary matrix Q as a
///                   product of min(M,N) elementary reflectors. The matrix Q
///                   is represented as a product of elementary reflectors
///                     Q = H(1) H(2) . . . H(k), where k = min(M,N).
///                   Each H(i) has the form
///                     H(i) = I - Tau(i) * v * v**H,
///                   where v is a vector with v(1:i-1) = 0 and v(i) = 1;
///                   v(i+1:M) is stored on exit in A(i+1:M,i).
/// \param Tau [out]  One-dimensional array of size min(M,N) that contains the
///                   scalar factors of the elementary reflectors.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class AMatrix, class TauArray, class InfoArray>
void geqrf(const AMatrix& A, const TauArray& Tau, const InfoArray& Info) {
  typename AMatrix::execution_space space{};
  geqrf(space, A, Tau, Info);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GEQRF_HPP_
