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

#ifndef KOKKOSBLAS2_SYR_HPP_
#define KOKKOSBLAS2_SYR_HPP_

#include "KokkosKernels_helpers.hpp"

#include <KokkosBlas2_syr_spec.hpp>

namespace KokkosBlas {

/// \brief Rank-1 update (just lower portion or just upper portion) of a
///        symmetric/Hermitian matrix: A = A + alpha * x * x^{T,H}.
///
///        Important note 1: this routine encapsulates the syr() and her()
///        routines specified in BLAS documentations. It has the purpose of
///        updating a symmetric (or Hermitian) matrix A in such a way that
///        it continues to be symmetric (or Hermitian). Therefore, in
///        Hermitian cases, the parameter alpha must be real.
///
///        Important note 2: however, this routine will honor all parameters
///        passed to it, even if A is not symmetric or not Hermitian, and
///        even if a complex alpha is supplied in Hermitian cases. Moreover,
///        this routine will always compute either the lower portion or the
///        upper portion (per user's request) of the final matrix A. So, in
///        order to obtain meaningful results, the user must make sure to
///        follow the conditions specified in the "important note 1" above.
///
///        Important note 3: if TPL is enabled, this routine will call the
///        third party library BLAS routines whenever the parameters passed
///        are consistent with the parameters expected by the corresponding
///        TPL routine. If not, then this routine will route the execution
///        to the kokkos-kernels implementation, thus honoring all
///        parameters passed, as stated in the "important note 2" above.
///
/// \tparam ExecutionSpace The type of execution space
/// \tparam XViewType      Input vector, as a 1-D Kokkos::View
/// \tparam AViewType      Input/Output matrix, as a 2-D Kokkos::View
///
/// \param space [in]     Execution space instance on which to run the kernel.
///                       This may contain information about which stream to
///                       run on.
/// \param trans [in]     "T" or "t" for transpose, "H" or "h" for Hermitian.
///                       Only the first character is taken into account.
/// \param uplo  [in]     "U" or "u" for upper portion, "L" or "l" for lower
///                       portion. Only the first character is taken into
///                       account.
/// \param alpha [in]     Input coefficient of x * x^{T,H}
/// \param x     [in]     Input vector, as a 1-D Kokkos::View
/// \param A     [in/out] Output matrix, as a nonconst 2-D Kokkos::View
template <class ExecutionSpace, class XViewType, class AViewType>
void syr(const ExecutionSpace& space, const char trans[], const char uplo[],
         const typename AViewType::const_value_type& alpha, const XViewType& x, const AViewType& A) {
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AViewType::memory_space>::accessible,
                "AViewType memory space must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XViewType::memory_space>::accessible,
                "XViewType memory space must be accessible from ExecutionSpace");

  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");

  static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
  static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if ((A.extent(0) != x.extent(0)) || (A.extent(1) != x.extent(0))) {
    std::ostringstream os;
    os << "KokkosBlas::syr: Dimensions of A, x: "
       << "A is " << A.extent(0) << " by " << A.extent(1) << ", x has size " << x.extent(0);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if ((trans[0] == 'T') || (trans[0] == 't') || (trans[0] == 'H') || (trans[0] == 'h')) {
    // Ok
  } else {
    std::ostringstream os;
    os << "KokkosBlas2::syr(): invalid trans[0] = '" << trans[0] << "'. It must be equal to 'T' or 't' or 'H' or 'h'";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if ((uplo[0] == 'U') || (uplo[0] == 'u') || (uplo[0] == 'L') || (uplo[0] == 'l')) {
    // Ok
  } else {
    std::ostringstream oss;
    oss << "KokkosBlas2::syr(): invalid uplo[0] = " << uplo[0] << "'. It must be equal to 'U' or 'u' or 'L' or 'l'";
    throw std::runtime_error(oss.str());
  }

  if ((A.extent(0) == 0) || (A.extent(1) == 0)) {
    return;
  }

  using ALayout = typename AViewType::array_layout;

  // Minimize the number of Impl::SYR instantiations, by standardizing
  // on particular View specializations for its template parameters.
  using XVT = Kokkos::View<typename XViewType::const_value_type*,
                           typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<XViewType, ALayout>::array_layout,
                           typename XViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  using AVT = Kokkos::View<typename AViewType::non_const_value_type**, ALayout, typename AViewType::device_type,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  Impl::SYR<ExecutionSpace, XVT, AVT>::syr(space, trans, uplo, alpha, x, A);
}

/// \brief Rank-1 update (just lower portion or just upper portion) of a
///        symmetric/Hermitian matrix: A = A + alpha * x * x^{T,H}.
///
///        Important note 1: this routine encapsulates the syr() and her()
///        routines specified in BLAS documentations. It has the purpose of
///        updating a symmetric (or Hermitian) matrix A in such a way that
///        it continues to be symmetric (or Hermitian). Therefore, in
///        Hermitian cases, the parameter alpha must be real.
///
///        Important note 2: however, this routine will honor all parameters
///        passed to it, even if A is not symmetric or not Hermitian, and
///        even if a complex alpha is supplied in Hermitian cases. Moreover,
///        this routine will always compute either the lower portion or the
///        upper portion (per user's request) of the final matrix A. So, in
///        order to obtain meaningful results, the user must make sure to
///        follow the conditions specified in the "important note 1" above.
///
///        Important note 3: if TPL is enabled, this routine will call the
///        third party library BLAS routines whenever the parameters passed
///        are consistent with the parameters expected by the corresponding
///        TPL routine. If not, then this routine will route the execution
///        to the kokkos-kernels implementation, thus honoring all
///        parameters passed, as stated in the "important note 2" above.
///
/// \tparam XViewType Input vector, as a 1-D Kokkos::View
/// \tparam AViewType Input/Output matrix, as a 2-D Kokkos::View
///
/// \param trans [in]     "T" or "t" for transpose, "H" or "h" for Hermitian.
///                       Only the first character is taken into account.
/// \param uplo  [in]     "U" or "u" for upper portion, "L" or "l" for lower
///                       portion. Only the first character is taken into
///                       account.
/// \param alpha [in]     Input coefficient of x * x^{T,H}
/// \param x     [in]     Input vector, as a 1-D Kokkos::View
/// \param A     [in/out] Output matrix, as a nonconst 2-D Kokkos::View
template <class XViewType, class AViewType>
void syr(const char trans[], const char uplo[], const typename AViewType::const_value_type& alpha, const XViewType& x,
         const AViewType& A) {
  const typename AViewType::execution_space space = typename AViewType::execution_space();
  syr<typename AViewType::execution_space, XViewType, AViewType>(space, trans, uplo, alpha, x, A);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS2_SYR_HPP_
