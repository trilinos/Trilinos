// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosLapack_gegqr.hpp
/// \brief compute Q from the QR factorization
///
/// This file provides KokkosLapack::gegqr. This function performs a
/// local (no MPI) computation of Q using the output from geqrf.

#ifndef KOKKOSLAPACK_GEGQR_HPP_
#define KOKKOSLAPACK_GEGQR_HPP_

#include <type_traits>

#include "KokkosLapack_gegqr_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief Computes the Q factor from a QR decomposition
///
/// \tparam ExecutionSpace The space where the kernel will run.
/// \tparam AMatrix        Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray       Type of array Tau, as a 1-D Kokkos::View.
/// \tparam InfoArray      Type of array Info, as a 1-D Kokkos::View.
///
/// \param space [in] Execution space instance used to specified how to execute
///                   the gegqr kernels.
/// \param k [in]     The number of reflectors to use to compute Q.
/// \param A [in-out] The i-th column must contain the vector which defines the
///                   elementary reflector H(i), for i = 1,2,...,k, as returned by
///                   GEQRF in the first k columns of its array argument A. On return
///                   A contains Q.
/// \param Tau [in]   One-dimensional array of size k. TAU(i) must contain the scalar
///                   factor of the elementary reflector H(i), as returned by GEQRF.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
void gegqr(const ExecutionSpace& space, const int k, const AMatrix& A, const TauArray& Tau, const InfoArray& Info) {
  // NOTE: Currently, KokkosLapack::gegqr only supports LAPACK, cuSOLVER and
  // rocSOLVER TPLs.

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename TauArray::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename InfoArray::memory_space>::accessible);

  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::gegqr: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<TauArray>::value, "KokkosLapack::gegqr: Tau must be Kokkos::View.");
  static_assert(Kokkos::is_view<InfoArray>::value, "KokkosLapack::gegqr: Info must be Kokkos::View.");

  static_assert(static_cast<int>(AMatrix::rank) == 2, "KokkosLapack::gegqr: A must have rank 2.");
  static_assert(static_cast<int>(TauArray::rank) == 1, "KokkosLapack::gegqr: Tau must have rank 1.");
  static_assert(static_cast<int>(InfoArray::rank) == 1, "KokkosLapack::gegqr: Info must have rank 1.");

  static_assert(std::is_same_v<typename InfoArray::non_const_value_type, int>,
                "KokkosLapack::gegqr: Info must be an array of integers.");

  static_assert(std::is_same_v<typename AMatrix::value_type, typename AMatrix::non_const_value_type>,
                "KokkosLapack::gegqr: AMatrix must store non const values.");

  const int64_t m     = A.extent(0);
  const int64_t n     = A.extent(1);
  const int64_t tau0  = Tau.extent(0);
  const int64_t info0 = Info.extent(0);

  if (m < n) {
    std::ostringstream os;
    os << "KokkosLapack::gegqr: m must be larger or equal to n, m=" << m << ", n=" << n;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (n < k || k < 0) {
    std::ostringstream os;
    os << "KokkosLapack::gegqr: k= " << k << ", must be positive and less or equal to n=" << n;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check validity of dimensions
  if (tau0 != std::min(m, n)) {
    std::ostringstream os;
    os << "KokkosLapack::gegqr: length of Tau must be equal to min(m,n): "
       << " A: " << m << " x " << n << ", Tau length = " << tau0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (info0 < 1) {
    std::ostringstream os;
    os << "KokkosLapack::gegqr: length of Info must be at least 1, Info length = " << info0;
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

  KokkosLapack::Impl::GEGQR<ExecutionSpace, AMatrix_Internal, TauArray_Internal, InfoArray_Internal>::gegqr(
      space, k, A_i, Tau_i, Info_i);
}

/// \brief Computes the Q factor from a QR decomposition
///
/// \tparam AMatrix        Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray       Type of array Tau, as a 1-D Kokkos::View.
/// \tparam InfoArray      Type of array Info, as a 1-D Kokkos::View.
///
/// \param space [in] Execution space instance used to specified how to execute
///                   the gegqr kernels.
/// \param A [in-out] The i-th column must contain the vector which defines the
///                   elementary reflector H(i), for i = 1,2,...,k, as returned by
///                   GEQRF in the first k columns of its array argument A. On return
///                   A contains Q.
/// \param Tau [in]   One-dimensional array of size k. TAU(i) must contain the scalar
///                   factor of the elementary reflector H(i), as returned by GEQRF.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class AMatrix, class TauArray, class InfoArray>
void gegqr(const int k, const AMatrix& A, const TauArray& Tau, const InfoArray& Info) {
  typename AMatrix::execution_space space{};
  gegqr(space, k, A, Tau, Info);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GEGQR_HPP_
