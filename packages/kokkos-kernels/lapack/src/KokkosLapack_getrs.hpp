// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosLapack_getrs.hpp
/// \brief Solve a dense linear system using an existing LU factorization.

#ifndef KOKKOSLAPACK_GETRS_HPP_
#define KOKKOSLAPACK_GETRS_HPP_

#include <type_traits>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Error.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosLapack_getrs_spec.hpp>

namespace KokkosLapack {

/// \brief Solve op(A) * X = B using the factors produced by getrf.
///
/// \param space [in] Execution space instance used for the operation.
/// \param trans [in] "N" for A, "T" for transpose(A), or "C" for
///                   conjugate-transpose(A). Lowercase is also accepted.
/// \param A [in] Square rank-2 LayoutLeft view containing the L and U
///               factors from getrf; the diagonal of L is implicitly one.
///               Const element types are accepted.
/// \param Ipiv [in] Rank-1 view of signed integral pivot indices of length N,
///                  using the one-based row swaps returned by getrf.
///                  Const element types are accepted. No-pivot mode is not supported.
/// \param B [in,out] Rank-2 LayoutLeft view of size N-by-NRHS. On entry,
///                   the right-hand sides; on successful exit, the solutions.
///                   Its scalar type must match A and its elements must be writable.
/// \param Info [out] Writable rank-1 int view with at least one element.
///                   Info(0) is zero on success. Additional elements are untouched.
///
/// All views must be accessible from ExecutionSpace. A and Ipiv are not modified
/// and must not alias the outputs. The caller must provide a successful,
/// nonsingular getrf factorization; pivot values and singularity are not checked.
/// Invalid dimensions or trans throw before modifying outputs. If N or NRHS is
/// zero, Info(0) is set to zero using space and B is untouched. Callers must
/// synchronize space before reading results on the host when needed.
///
/// \note Host LAPACK supports float, double, and Kokkos::complex variants,
///       with int pivots in HostSpace on Serial, OpenMP, and Threads.
///       Unsupported combinations report a missing implementation.
template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
void getrs(const ExecutionSpace& space, const char trans[], const AMatrix& A, const IpivView& Ipiv, const BMatrix& B,
           const InfoView& Info) {
  static_assert(Kokkos::is_execution_space_v<ExecutionSpace>,
                "KokkosLapack::getrs: ExecutionSpace must be a Kokkos execution space.");
  static_assert(Kokkos::is_view_v<AMatrix>, "KokkosLapack::getrs: A must be a Kokkos::View.");
  static_assert(Kokkos::is_view_v<IpivView>, "KokkosLapack::getrs: Ipiv must be a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BMatrix>, "KokkosLapack::getrs: B must be a Kokkos::View.");
  static_assert(Kokkos::is_view_v<InfoView>, "KokkosLapack::getrs: Info must be a Kokkos::View.");
  static_assert(AMatrix::rank == 2, "KokkosLapack::getrs: A must have rank 2.");
  static_assert(BMatrix::rank == 2, "KokkosLapack::getrs: B must have rank 2.");
  static_assert(IpivView::rank == 1, "KokkosLapack::getrs: Ipiv must have rank 1.");
  static_assert(InfoView::rank == 1, "KokkosLapack::getrs: Info must have rank 1.");
  static_assert(std::is_same_v<typename AMatrix::array_layout, Kokkos::LayoutLeft>,
                "KokkosLapack::getrs: A must have LayoutLeft.");
  static_assert(std::is_same_v<typename BMatrix::array_layout, Kokkos::LayoutLeft>,
                "KokkosLapack::getrs: B must have LayoutLeft.");
  static_assert(!std::is_const_v<typename BMatrix::value_type>, "KokkosLapack::getrs: B must be writable.");

  static_assert(std::is_integral_v<typename IpivView::non_const_value_type>,
                "KokkosLapack::getrs: Ipiv must contain integers.");
  static_assert(std::is_same_v<typename InfoView::value_type, int>,
                "KokkosLapack::getrs: Info must contain writable int elements.");

  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename IpivView::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename BMatrix::memory_space>::accessible);
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename InfoView::memory_space>::accessible);

  if (trans == nullptr || !(trans[0] == 'N' || trans[0] == 'n' || trans[0] == 'T' || trans[0] == 't' ||
                            trans[0] == 'C' || trans[0] == 'c')) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: trans must be N, T, or C.");
  }

  if (A.extent(0) != A.extent(1)) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: A must be square.");
  }
  if (B.extent(0) != A.extent(0)) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: B.extent(0) must be equal to A.extent(0).");
  }
  if (Ipiv.extent(0) != A.extent(0)) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: Ipiv must have length A.extent(0).");
  }
  if (Info.extent(0) < 1) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: Info must have at least one element.");
  }

  // Check for possiblity of a quick return
  if (A.extent(0) == 0 || B.extent(1) == 0) {
    Kokkos::deep_copy(space, Kokkos::subview(Info, 0), 0);
    return;
  }

  // Perform some type unification on the views
  // to hit more ETI and TPL paths
  using ALayout         = typename AMatrix::array_layout;
  using AMatrixInternal = Kokkos::View<typename AMatrix::const_data_type, ALayout,
                                       Kokkos::Device<ExecutionSpace, typename AMatrix::memory_space>,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using IpivViewInternal =
      Kokkos::View<typename IpivView::const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<IpivView, ALayout>::array_layout,
                   Kokkos::Device<ExecutionSpace, typename IpivView::memory_space>,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using BMatrixInternal = Kokkos::View<typename BMatrix::non_const_data_type, typename BMatrix::array_layout,
                                       Kokkos::Device<ExecutionSpace, typename BMatrix::memory_space>,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InfoViewInternal =
      Kokkos::View<typename InfoView::non_const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<InfoView, ALayout>::array_layout,
                   Kokkos::Device<ExecutionSpace, typename InfoView::memory_space>,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  AMatrixInternal A_i(A);
  IpivViewInternal Ipiv_i(Ipiv);
  BMatrixInternal B_i(B);
  InfoViewInternal Info_i(Info);
  const char trans_i[] = {trans[0] == 'n' ? 'N' : trans[0] == 't' ? 'T' : trans[0] == 'c' ? 'C' : trans[0], '\0'};
  Impl::GETRS<ExecutionSpace, AMatrixInternal, IpivViewInternal, BMatrixInternal, InfoViewInternal>::getrs(
      space, trans_i, A_i, Ipiv_i, B_i, Info_i);
}

/// \brief GETRS overload using A's default execution space instance.
/// \see getrs(const ExecutionSpace&, const char[], const AMatrix&, const IpivView&, const BMatrix&, const InfoView&)
template <class AMatrix, class IpivView, class BMatrix, class InfoView>
void getrs(const char trans[], const AMatrix& A, const IpivView& Ipiv, const BMatrix& B, const InfoView& Info) {
  getrs(typename AMatrix::execution_space{}, trans, A, Ipiv, B, Info);
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GETRS_HPP_
