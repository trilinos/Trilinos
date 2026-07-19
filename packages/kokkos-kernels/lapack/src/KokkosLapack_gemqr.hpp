// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosLapack_gemqr.hpp
/// \brief QR multiply by Q factor
///
/// This file provides KokkosLapack::gemqr. This function performs a
/// local (no MPI) multiplication of Q, computed by geqrf, and a matrix.

#ifndef KOKKOSLAPACK_GEMQR_HPP_
#define KOKKOSLAPACK_GEMQR_HPP_

#include <type_traits>

#include "KokkosLapack_gemqr_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief Multiplies matrix C with the Q factor, from a QR decomposition
///
/// \tparam ExecutionSpace The space where the kernel will run.
/// \tparam AMatrix        Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray       Type of array Tau, as a 1-D Kokkos::View.
/// \tparam CMatrix        Type of matrix C, as a 2-D Kokkos::View.
/// \tparam InfoArray      Type of array Info, as a 1-D Kokkos::View.
///
/// \param space [in] Execution space instance used to specified how to execute
///                   the gemqr kernels.
/// \param side [in]  The side of C to be used to multiply by Q
/// \param trans [in] Operation applied to Q for the multiplcation: none, transpose
///                   or hermitian
/// \param A [in]     The i-th column must contain the vector which defines the
///                   elementary reflector H(i), for i = 1,2,...,k, as returned by
///                   GEQRF in the first k columns of its array argument A.
/// \param Tau [in]   One-dimensional array of size k. TAU(i) must contain the scalar
///                   factor of the elementary reflector H(i), as returned by GEQRF.
/// \param C [in,out] On entry, the M-by-N matrix C.
///                   On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class ExecutionSpace, class AMatrix, class TauArray, class CMatrix, class InfoArray>
void gemqr(const ExecutionSpace& space, const char side[], const char trans[], const AMatrix& A, const TauArray& Tau,
           const CMatrix& C, const InfoArray& Info) {
  // NOTE: Currently, KokkosLapack::gemqr only supports LAPACK, cuSOLVER and
  // rocSOLVER TPLs.

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename TauArray::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename CMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename InfoArray::memory_space>::accessible);

  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::gemqr: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<TauArray>::value, "KokkosLapack::gemqr: Tau must be Kokkos::View.");
  static_assert(Kokkos::is_view<CMatrix>::value, "KokkosLapack::gemqr: C must be a Kokkos::View.");
  static_assert(Kokkos::is_view<InfoArray>::value, "KokkosLapack::gemqr: Info must be Kokkos::View.");

  static_assert(static_cast<int>(AMatrix::rank) == 2, "KokkosLapack::gemqr: A must have rank 2.");
  static_assert(static_cast<int>(TauArray::rank) == 1, "KokkosLapack::gemqr: Tau must have rank 1.");
  static_assert(static_cast<int>(CMatrix::rank) == 2, "KokkosLapack::gemqr: C must have rank 2.");
  static_assert(static_cast<int>(InfoArray::rank) == 1, "KokkosLapack::gemqr: Info must have rank 1.");

  static_assert(std::is_same_v<typename InfoArray::non_const_value_type, int>,
                "KokkosLapack::gemqr: Info must be an array of integers.");

  if (side == nullptr || (side[0] != 'L' && side[0] != 'l' && side[0] != 'R' && side[0] != 'r')) {
    std::ostringstream os;
    os << "KokkosLapack::gemrf: side must be \"L\", \"l\", \"R\" or \"r\"";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (trans == nullptr || (trans[0] != 'N' && trans[0] != 'n' && trans[0] != 'T' && trans[0] != 't' &&
                           trans[0] != 'C' && trans[0] != 'c')) {
    std::ostringstream os;
    os << "KokkosLapack::gemrf: trans must be \"N\", \"n\", \"T\", \"t\", \"C\" or \"c\"";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  const int64_t m     = A.extent(0);
  const int64_t n     = A.extent(1);
  const int64_t tau0  = Tau.extent(0);
  const int64_t info0 = Info.extent(0);

  // Check validity of dimensions
  if (tau0 != std::min(m, n)) {
    std::ostringstream os;
    os << "KokkosLapack::gemqr: length of Tau must be equal to min(m,n): "
       << " A: " << m << " x " << n << ", Tau length = " << tau0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (info0 < 1) {
    std::ostringstream os;
    os << "KokkosLapack::gemqr: length of Info must be at least 1, Info length = " << info0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if ((side[0] == 'L' || side[0] == 'l') && C.extent_int(0) != m) {
    std::ostringstream os;
    os << "KokkosLapack::gemqr: multiplying on the left but A.extent(0) != C.extent(0)";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if ((side[0] == 'R' || side[0] == 'r') && C.extent_int(0) != n) {
    std::ostringstream os;
    os << "KokkosLapack::gemqr: multiplying on the right but A.extent(1) != C.extent(0)";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using AMatrix_Internal   = Kokkos::View<typename AMatrix::non_const_value_type**, typename AMatrix::array_layout,
                                        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using TauArray_Internal  = Kokkos::View<typename TauArray::non_const_value_type*, typename TauArray::array_layout,
                                         typename TauArray::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using CMatrix_Internal   = Kokkos::View<typename CMatrix::non_const_value_type**, typename CMatrix::array_layout,
                                        typename CMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InfoArray_Internal = Kokkos::View<typename InfoArray::non_const_value_type*, typename InfoArray::array_layout,
                                          typename InfoArray::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  AMatrix_Internal A_i      = A;
  TauArray_Internal Tau_i   = Tau;
  CMatrix_Internal C_i      = C;
  InfoArray_Internal Info_i = Info;

  KokkosLapack::Impl::GEMQR<ExecutionSpace, AMatrix_Internal, TauArray_Internal, CMatrix_Internal,
                            InfoArray_Internal>::gemqr(space, side, trans, A_i, Tau_i, C_i, Info_i);
}

/// \brief Multiplies matrix C with the Q factor, from a QR decomposition
///
/// \tparam AMatrix        Type of matrix A, as a 2-D Kokkos::View.
/// \tparam TauArray       Type of array Tau, as a 1-D Kokkos::View.
/// \tparam CMatrix        Type of matrix C, as a 2-D Kokkos::View.
/// \tparam InfoArray      Type of array Info, as a 1-D Kokkos::View.
///
/// \param side [in]  The side of C to be used to multiply by Q
/// \param trans [in] Operation applied to Q for the multiplication: none, transpose
///                   or hermitian
/// \param A [in]     The i-th column must contain the vector which defines the
///                   elementary reflector H(i), for i = 1,2,...,k, as returned by
///                   GEQRF in the first k columns of its array argument A.
/// \param Tau [in]   One-dimensional array of size k. TAU(i) must contain the scalar
///                   factor of the elementary reflector H(i), as returned by GEQRF.
/// \param C [in,out] On entry, the M-by-N matrix C.
///                   On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
template <class AMatrix, class TauArray, class CMatrix, class InfoArray>
void gemqr(const char side[], const char trans[], const AMatrix& A, const TauArray& Tau, const CMatrix& C,
           const InfoArray& Info) {
  typename AMatrix::execution_space space{};
  gemqr(space, side, trans, A, Tau, C, Info);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GEMQR_HPP_
