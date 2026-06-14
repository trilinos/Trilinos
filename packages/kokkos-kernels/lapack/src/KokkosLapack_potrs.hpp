// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_POTRS_HPP_
#define KOKKOSLAPACK_POTRS_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosLapack_potrs_spec.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief solves a system of linear equations with a Hermitian
///        positive definite matrix A using the Cholesky factorization
///
/// Solves a system of linear equations :math:`A X = B` where :math:`A` is a
/// symmetric (or Hermitian) positive definite matrix whose Cholesky factorization
/// has already been computed by :code:`KokkosLapack::potrf`.
///
/// Given the factorization produced by ``potrf``:
///    A = U^H U    if uplo = 'U'
///    A = L L^H    if uplo = 'L'
/// potrs solves for `X` by two triangular solves, overwriting `B`
/// with the solution `X`.
///
/// \tparam ExecutionSpace The space where the kernel will run.
/// \tparam AViewType [in] Type of matrix A, as a 2-D Kokkos::View (LayoutLeft!)
/// \tparam BViewType [in] Type of matrix B, as a 2-D Kokkos::View (LayoutLeft!)
///
/// \param space [in] Execution space instance used to specify how to execute
///                   the potrf kernels.
/// \param uplo  [in] 'U': upper triangle of A is stored; 'L': lower triangle.
/// \param A     [in] Square 2-D Kokkos::View of dimension N x N. Comes from potrf
/// \param B     [in,out] On entry, the right hand side matrix B. On exit, the solution matrix X.
///
template <class ExecutionSpace, class AViewType, class BViewType>
void potrs(const ExecutionSpace& space, const char uplo[], const AViewType& A, BViewType& B) {
  static_assert(Kokkos::is_execution_space_v<ExecutionSpace>,
                "KokkosLapack::potrs: ExecutionSpace must be a valid Kokkos execution space");
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosLapack::potrs: AViewType must be a Kokkos::View");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosLapack::potrs: BViewType must be a Kokkos::View");
  static_assert(static_cast<int>(AViewType::rank) == 2, "KokkosLapack::potrs: A must have rank 2.");
  static_assert(static_cast<int>(BViewType::rank) == 2, "KokkosLapack::potrs: B must have rank 2.");
  static_assert(!std::is_const_v<typename BViewType::value_type>,
                "KokkosLapack::potrs: B should not have const value type");
  static_assert(std::is_same_v<typename AViewType::array_layout, Kokkos::LayoutLeft>,
                "KokkosLapack::potrs: A must have Kokkos::LayoutLeft (column-major) layout, "
                "as required by LAPACK/cuSOLVER/rocSOLVER.");
  static_assert(std::is_same_v<typename BViewType::array_layout, Kokkos::LayoutLeft>,
                "KokkosLapack::potrs: A must have Kokkos::LayoutLeft (column-major) layout, "
                "as required by LAPACK/cuSOLVER/rocSOLVER.");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AViewType::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename BViewType::memory_space>::accessible);

  KK_REQUIRE_MSG(A.extent(0) == A.extent(1),
                 "KokkosLapack::potrs: A must be square, got " << A.extent(0) << " x " << A.extent(1));
  KK_REQUIRE_MSG(A.extent(1) == B.extent(0), "KokkosLapack::potrs: A dim 1 must be compatible with B dim 0, got A1 "
                                                 << A.extent(1) << " vs B0" << B.extent(0));

  // Convert views to unmanaged
  using AViewInternalType = Kokkos::View<typename AViewType::const_data_type, typename AViewType::array_layout,
                                         typename AViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using BViewInternalType = Kokkos::View<typename BViewType::data_type, typename BViewType::array_layout,
                                         typename BViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  AViewInternalType uA(A);
  BViewInternalType uB(B);

  Impl::Potrs<ExecutionSpace, AViewInternalType, BViewInternalType>::potrs(space, uplo, uA, uB);
}

// Overload without execution space (uses default)
template <class AViewType, class BViewType>
void potrs(const char uplo[], const AViewType& A, BViewType& B) {
  potrs(typename AViewType::execution_space{}, uplo, A, B);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_POTRS_HPP_
