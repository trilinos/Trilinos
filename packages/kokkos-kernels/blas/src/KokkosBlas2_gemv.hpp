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
#ifndef KOKKOSBLAS2_GEMV_HPP_
#define KOKKOSBLAS2_GEMV_HPP_

/// \file KokkosBlas2_gemv.hpp
/// \brief BLAS 2 kernels specifically optimized for typical
///   Tpetra::MultiVector use cases.

#include <KokkosBlas2_gemv_spec.hpp>
#include <KokkosBlas2_serial_gemv.hpp>
#include <KokkosBlas2_team_gemv.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>
#include <sstream>
#include <type_traits>  // requires C++11, but so does Kokkos

namespace KokkosBlas {

/// \brief Dense matrix-vector multiply: y = beta*y + alpha*A*x.
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam XViewType Input vector, as a 1-D Kokkos::View
/// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
/// \tparam AlphaCoeffType Type of input coefficient alpha
/// \tparam BetaCoeffType Type of input coefficient beta
///
/// \param space [in] execution space instance on which to run the
///   kernel. This may contain information about which stream to
///   run on.
/// \param trans [in] "N" for non-transpose, "T" for transpose, "C"
///   for conjugate transpose.  All characters after the first are
///   ignored.  This works just like the BLAS routines.
/// \param alpha [in] Input coefficient of A*x
/// \param A [in] Input matrix, as a 2-D Kokkos::View
/// \param x [in] Input vector, as a 1-D Kokkos::View
/// \param beta [in] Input coefficient of y
/// \param y [in/out] Output vector, as a nonconst 1-D Kokkos::View
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType>
void gemv(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
          const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta, const YViewType& y) {
  static_assert(Kokkos::is_execution_space_v<ExecutionSpace>,
                "KokkosBlas::gemv: ExecutionSpace must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<AViewType>::value, "KokkosBlas::gemv: AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<XViewType>::value, "KokkosBlas::gemv: XViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "KokkosBlas::gemv: YViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2, "KokkosBlas::gemv: AViewType must have rank 2.");
  static_assert(static_cast<int>(XViewType::rank) == 1, "KokkosBlas::gemv: XViewType must have rank 1.");
  static_assert(static_cast<int>(YViewType::rank) == 1, "KokkosBlas::gemv: YViewType must have rank 1.");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AViewType::memory_space>::accessible,
                "KokkosBlas::gemv: AViewType must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XViewType::memory_space>::accessible,
                "KokkosBlas::gemv: XViewType must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename YViewType::memory_space>::accessible,
                "KokkosBlas::gemv: YViewType must be accessible from ExecutionSpace");

  // Check compatibility of dimensions at run time.
  if (trans[0] == 'N' || trans[0] == 'n') {
    if (A.extent(0) != y.extent(0) || A.extent(1) != x.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::gemv: Dimensions of A, x, and y do not match: "
         << "A: " << A.extent(0) << " x " << A.extent(1) << ", x: " << x.extent(0) << ", y: " << y.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else if (trans[0] == 'T' || trans[0] == 't' || trans[0] == 'C' || trans[0] == 'c' || trans[0] == 'H' ||
             trans[0] == 'h') {
    if (A.extent(1) != y.extent(0) || A.extent(0) != x.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::dot: Dimensions of A, x, and y do not match: "
         << "A: " << A.extent(0) << " x " << A.extent(1) << ", x: " << x.extent(0) << ", y: " << y.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    std::ostringstream os;
    os << "KokkosBlas::gemv: trans[0] = '" << trans[0]
       << "'.  Valid values "
          "include 'N' (No transpose), 'T' (Transpose), and 'C' (Conjugate "
          "transpose).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using ALayout = typename AViewType::array_layout;

  // Minimize the number of Impl::GEMV instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  typedef Kokkos::View<typename AViewType::const_value_type**, ALayout, typename AViewType::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      AVT;
  typedef Kokkos::View<typename XViewType::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<XViewType, ALayout>::array_layout,
                       typename XViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XVT;
  typedef Kokkos::View<typename YViewType::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<YViewType, ALayout>::array_layout,
                       typename YViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YVT;

  // Degenerate case is essentially same as scal - use fallback impl
  // to avoid potential (unlikely?) circular dependence issues by including
  // other KokkosBlas headers
  bool useFallback = A.extent(0) == 0 || A.extent(1) == 0;
  // If A is LayoutRight and we have the BLAS, cuBLAS or rocBLAS TPL, use
  // fallback because those only support LayoutLeft
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  useFallback = useFallback ||
                (tolower(*trans) == 'c' && std::is_same<typename AViewType::array_layout, Kokkos::LayoutRight>::value &&
                 std::is_same<typename AViewType::memory_space, Kokkos::CudaSpace>::value);
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
  useFallback = useFallback ||
                (tolower(*trans) == 'c' && std::is_same<typename AViewType::array_layout, Kokkos::LayoutRight>::value &&
                 std::is_same<typename AViewType::memory_space, Kokkos::HIPSpace>::value);
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
  useFallback = useFallback ||
                (tolower(*trans) == 'c' && std::is_same<typename AViewType::array_layout, Kokkos::LayoutRight>::value &&
                 std::is_same<typename AViewType::memory_space, Kokkos::HostSpace>::value);
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#ifdef KOKKOS_ENABLE_SYCL
  // oneMKL supports both row-major and column-major of A
  // but only supports oneapi::mkl::transpose::nontrans op
  useFallback =
      useFallback || !std::is_same_v<typename AViewType::memory_space, Kokkos::Experimental::SYCLDeviceUSMSpace>;
#endif
#endif

  if (useFallback) {
    const bool eti_spec_avail = KokkosBlas::Impl::gemv_eti_spec_avail<ExecutionSpace, AVT, XVT, YVT>::value;
    typedef Impl::GEMV<ExecutionSpace, AVT, XVT, YVT, false, eti_spec_avail> fallback_impl_type;
    fallback_impl_type::gemv(space, trans, alpha, A, x, beta, y);
  } else {
    typedef Impl::GEMV<ExecutionSpace, AVT, XVT, YVT> impl_type;
    impl_type::gemv(space, trans, alpha, A, x, beta, y);
  }
}

/// \brief Dense matrix-vector multiply: y = beta*y + alpha*A*x.
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam XViewType Input vector, as a 1-D Kokkos::View
/// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
/// \tparam AlphaCoeffType Type of input coefficient alpha
/// \tparam BetaCoeffType Type of input coefficient beta
///
/// \param trans [in] "N" for non-transpose, "T" for transpose, "C"
///   for conjugate transpose.  All characters after the first are
///   ignored.  This works just like the BLAS routines.
/// \param alpha [in] Input coefficient of A*x
/// \param A [in] Input matrix, as a 2-D Kokkos::View
/// \param x [in] Input vector, as a 1-D Kokkos::View
/// \param beta [in] Input coefficient of y
/// \param y [in/out] Output vector, as a nonconst 1-D Kokkos::View
template <class AViewType, class XViewType, class YViewType>
void gemv(const char trans[], typename AViewType::const_value_type& alpha, const AViewType& A, const XViewType& x,
          typename YViewType::const_value_type& beta, const YViewType& y) {
  gemv(typename AViewType::execution_space{}, trans, alpha, A, x, beta, y);
}

namespace Experimental {
///
/// Selective Interface
///
template <class ArgMode, class ArgAlgo>
struct Gemv {
  template <class MemberType, class MatrixType, class XVector, class YVector, class ScalarType>
  static void KOKKOS_INLINE_FUNCTION invoke(const MemberType& member, const char trans, const ScalarType& alpha,
                                            const MatrixType& A, const XVector& x, const ScalarType& beta,
                                            const YVector& y);
};

template <class ArgAlgo>
struct Gemv<Mode::Serial, ArgAlgo> {
  template <class MemberType, class MatrixType, class XVector, class YVector, class ScalarType>
  static void KOKKOS_INLINE_FUNCTION invoke(const MemberType& /*member*/, const char trans, const ScalarType& alpha,
                                            const MatrixType& A, const XVector& x, const ScalarType& beta,
                                            const YVector& y) {
    serial_gemv<ArgAlgo>(trans, alpha, A, x, beta, y);
  }
};

template <class ArgAlgo>
struct Gemv<Mode::Team, ArgAlgo> {
  template <class MemberType, class MatrixType, class XVector, class YVector, class ScalarType>
  static void KOKKOS_INLINE_FUNCTION invoke(const MemberType& member, const char trans, const ScalarType& alpha,
                                            const MatrixType& A, const XVector& x, const ScalarType& beta,
                                            const YVector& y) {
    team_gemv<ArgAlgo>(member, trans, alpha, A, x, beta, y);
  }
};

template <class ArgAlgo>
struct Gemv<Mode::TeamVector, ArgAlgo> {
  template <class MemberType, class MatrixType, class XVector, class YVector, class ScalarType>
  static void KOKKOS_INLINE_FUNCTION invoke(const MemberType& member, const char trans, const ScalarType& alpha,
                                            const MatrixType& A, const XVector& x, const ScalarType& beta,
                                            const YVector& y) {
    teamvector_gemv<ArgAlgo>(member, trans, alpha, A, x, beta, y);
  }
};

}  // namespace Experimental
}  // namespace KokkosBlas

#endif  // KOKKOS_BLAS2_MV_HPP_
