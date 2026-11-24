// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosLapack_getrf.hpp
/// \brief LU factorization
///
/// This file provides KokkosLapack::getrf. This function performs a
/// local (no MPI) LU factorization of a M-by-N matrix A.

#ifndef KOKKOSLAPACK_GETRF_HPP_
#define KOKKOSLAPACK_GETRF_HPP_

#include <type_traits>

#include <KokkosLapack_getrf_spec.hpp>
#include <KokkosKernels_Error.hpp>
#include <KokkosKernels_helpers.hpp>

namespace KokkosLapack {

/// \brief Computes a LU factorization of a matrix A
///
/// \tparam ExecutionSpace The space where the kernel will run.
/// \tparam AMatrix   Type of matrix A, as a rank-2 Kokkos::View.
/// \tparam IpivView  Type of array Ipiv, as a rank-1 Kokkos::View.
/// \tparam InfoView  Type of array Info, as a rank-1 Kokkos::View.
///
/// \param A [in,out] On entry, the M-by-N matrix to be decomposed.
///                   On exit, the L and U factors with L diagonal assumed
///                   to be unit.
/// \param Ipiv [out] One-dimensional array of size min(M,N) that contains the
///                   pivot ordering selected for the decomposition.
/// \param Info [out] One-dimensional array of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class ExecutionSpace, class AMatrix, class IpivView, class InfoView>
void getrf(const ExecutionSpace& space, const AMatrix& A, const IpivView& Ipiv, const InfoView& Info) {
  // NOTE: Currently, KokkosLapack::getrf only supports LAPACK, cuSOLVER and
  // rocSOLVER TPLs.
  //       cuSOLVER/rocSOLVER TPL should be enabled to call the cuSOLVER/rocSOLVER GPU
  //       interface for device views LAPACK TPL should be enabled to call the
  //       LAPACK interface for host views

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename IpivView::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename InfoView::memory_space>::accessible);

  static_assert(Kokkos::is_view<AMatrix>::value, "KokkosLapack::getrf: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view<IpivView>::value, "KokkosLapack::getrf: Ipiv must be Kokkos::View.");
  static_assert(Kokkos::is_view<InfoView>::value, "KokkosLapack::getrf: Info must be Kokkos::View.");

  static_assert(static_cast<int>(AMatrix::rank) == 2, "KokkosLapack::getrf: A must have rank 2.");
  static_assert(static_cast<int>(IpivView::rank) == 1, "KokkosLapack::getrf: Ipiv must have rank 1.");
  static_assert(static_cast<int>(InfoView::rank) == 1, "KokkosLapack::getrf: Info must have rank 1.");

  static_assert(std::is_same_v<typename InfoView::non_const_value_type, int>,
                "KokkosLapack::getrf: Info must be an array of integers.");

  const int64_t m     = A.extent(0);
  const int64_t n     = A.extent(1);
  const int64_t ipiv0 = Ipiv.extent(0);
  const int64_t info0 = Info.extent(0);

  // Check validity of dimensions
  if (ipiv0 != std::min(m, n)) {
    std::ostringstream os;
    os << "KokkosLapack::getrf: length of Ipiv must be equal to min(m,n): "
       << " A: " << m << " x " << n << ", Ipiv length = " << ipiv0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (info0 < 1) {
    std::ostringstream os;
    os << "KokkosLapack::getrf: length of Info must be at least 1, Info length = " << info0;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check for possible quick return
  if (A.extent(0) == 0 || A.extent(1) == 0) {
    Kokkos::deep_copy(space, Info, 0);
    return;
  }

  using ALayout          = typename AMatrix::array_layout;
  using AMatrix_Internal = Kokkos::View<typename AMatrix::non_const_value_type**, ALayout,
                                        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using IpivView_Internal =
      Kokkos::View<typename IpivView::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<IpivView, ALayout>::array_layout,
                   typename IpivView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InfoView_Internal =
      Kokkos::View<typename InfoView::non_const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<InfoView, ALayout>::array_layout,
                   typename InfoView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  AMatrix_Internal A_i     = A;
  IpivView_Internal Ipiv_i = Ipiv;
  InfoView_Internal Info_i = Info;

  KokkosLapack::Impl::GETRF<ExecutionSpace, AMatrix_Internal, IpivView_Internal, InfoView_Internal>::getrf(
      space, A_i, Ipiv_i, Info_i);
}

/// \brief Computes a LU factorization of a matrix A
///
/// \tparam AMatrix   Type of matrix A, as a rank-2 Kokkos::View.
/// \tparam IpivView  Type of array Ipiv, as a rank-1 Kokkos::View.
/// \tparam InfoView  Type of array Info, as a rank-1 Kokkos::View.
///
/// \param A [in,out] On entry, the M-by-N matrix to be decomposed.
///                   On exit, the L and U factors with L diagaonal assumed
///                   to be unit.
/// \param Ipiv [out] rank-1 view of size min(M,N) that contains the
///                   pivot ordering selected for the decomposition.
/// \param Info [out] rank-1 view of integers and of size 1:
///                   Info[0] = 0: successful exit
///                   Info[0] < 0: if equal to '-i', the i-th argument had an
///                                illegal value
///
template <class AMatrix, class IpivView, class InfoView>
void getrf(const AMatrix& A, const IpivView& Ipiv, const InfoView& Info) {
  typename AMatrix::execution_space space{};
  getrf(space, A, Ipiv, Info);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GETRF_HPP_
