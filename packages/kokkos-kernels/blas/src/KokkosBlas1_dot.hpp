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

#ifndef KOKKOSBLAS1_DOT_HPP_
#define KOKKOSBLAS1_DOT_HPP_

#include <KokkosBlas1_dot_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Return the dot product of the two vectors x and y.
///
/// \tparam execution_space the Kokkos execution space where the kernel
///         will be executed.
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
/// \tparam YVector Type of the second vector y; a 1-D Kokkos::View.
///
/// \param space [in] an execution space instance that may specify
///                   in which stream/queue the kernel will be executed.
/// \param x [in] Input 1-D View.
/// \param y [in] Input 1-D View.
///
/// \return The dot product result; a single value.
template <class execution_space, class XVector, class YVector,
          typename std::enable_if<Kokkos::is_execution_space_v<execution_space>, int>::type = 0>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type dot(
    const execution_space& space, const XVector& x, const YVector& y) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::dot: execution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<XVector>::value, "KokkosBlas::dot: XVector must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
                "KokkosBlas::dot: XVector must be accessible from execution_space");
  static_assert(Kokkos::is_view<YVector>::value, "KokkosBlas::dot: YVector must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YVector::memory_space>::accessible,
                "KokkosBlas::dot: YVector must be accessible from execution_space");
  static_assert((int)XVector::rank == (int)YVector::rank, "KokkosBlas::dot: Vector ranks do not match.");
  static_assert(XVector::rank == 1,
                "KokkosBlas::dot: "
                "Both Vector inputs must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (x.extent(0) != y.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::dot: Dimensions do not match: "
       << ", x: " << x.extent(0) << " x 1"
       << ", y: " << y.extent(0) << " x 1";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using XVector_Internal = Kokkos::View<typename XVector::const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                                        typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using YVector_Internal = Kokkos::View<typename YVector::const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
                                        typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using dot_type = typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type;
  // result_type is usually just dot_type, except:
  //  if dot_type is float, result_type is double
  //  if dot_type is complex<float>, result_type is complex<double>
  // These special cases are to maintain accuracy.
  using result_type = typename KokkosBlas::Impl::DotAccumulatingScalar<dot_type>::type;
  using RVector_Internal =
      Kokkos::View<dot_type, default_layout, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using RVector_Result =
      Kokkos::View<result_type, default_layout, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  XVector_Internal X = x;
  YVector_Internal Y = y;

  bool useFallback = false;
  if (useFallback) {
    // Even though RVector is the template parameter, Dot::dot has an overload
    // that accepts RVector_Internal (with the special accumulator, if dot_type
    // is 32-bit precision). Impl::Dot needs to support both cases, and it's
    // easier to do this with overloading than by extending the ETI to deal with
    // two different scalar types.
    result_type result{};
    RVector_Result R = RVector_Result(&result);
    Impl::DotSpecialAccumulator<execution_space, RVector_Internal, XVector_Internal, YVector_Internal>::dot(space, R, X,
                                                                                                            Y);
    space.fence();
    // mfh 22 Jan 2020: We need the line below because
    // Kokkos::complex<T> lacks a constructor that takes a
    // Kokkos::complex<U> with U != T.
    return Kokkos::Details::CastPossiblyComplex<dot_type, result_type>::cast(result);
  } else {
    dot_type result{};
    RVector_Internal R = RVector_Internal(&result);
    Impl::Dot<execution_space, RVector_Internal, XVector_Internal, YVector_Internal>::dot(space, R, X, Y);
    space.fence();
    return Kokkos::Details::CastPossiblyComplex<dot_type, result_type>::cast(result);
  }
}

/// \brief Return the dot product of the two vectors x and y.
///
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVector.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
/// \tparam YVector Type of the second vector y; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
/// \param y [in] Input 1-D View.
///
/// \return The dot product result; a single value.
template <class XVector, class YVector>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type dot(
    const XVector& x, const YVector& y) {
  return dot(typename XVector::execution_space{}, x, y);
}

/// \brief Compute the column-wise dot products of two multivectors.
///
/// This function is non-blocking and thread-safe.
///
/// \tparam execution_space the Kokkos execution space where the kernel
///         will be executed.
/// \tparam RV 0-D resp. 1-D output View
/// \tparam XMV 1-D resp. 2-D input View
/// \tparam YMV 1-D resp. 2-D input View
///
/// \param space [in] an execution space instance that may specify
///                   in which stream/queue the kernel will be executed.
/// \param R [out] Output 1-D or 0-D View to which to write results.
/// \param X [in] Input 2-D or 1-D View.
/// \param Y [in] Input 2-D or 1-D View.
///
/// This function implements a few different use cases:
/// <ul>
/// <li> If X and Y are both 1-D, then this is a single dot product.
///   R must be 0-D (a View of a single value). </li>
/// <li> If X and Y are both 2-D, then this function computes their
///   dot products columnwise.  R must be 1-D. </li>
/// <li> If X is 2-D and Y is 1-D, then this function computes the dot
///   product of each column of X, with Y, in turn.  R must be
///   1-D. </li>
/// <li> If X is 1-D and Y is 2-D, then this function computes the dot
///   product X with each column of Y, in turn.  R must be 1-D. </li>
/// </ul>
///
/// \note To implementers: We use enable_if here so that the compiler
///   doesn't confuse this version of dot() with the three-argument
///   version of dot() in Kokkos_Blas1.hpp.
template <class execution_space, class RV, class XMV, class YMV>
void dot(const execution_space& space, const RV& R, const XMV& X, const YMV& Y,
         typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::dot: excution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::dot: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::dot: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::dot: XMV must be accessible from execution_space.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::dot: "
                "Y is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible,
                "KokkosBlas::dot: XMV must be accessible from execution_space.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::dot: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(RV::rank == 0 || RV::rank == 1, "KokkosBlas::dot: R must have rank 0 or 1.");
  static_assert(XMV::rank == 1 || XMV::rank == 2, "KokkosBlas::dot: X must have rank 1 or 2.");
  static_assert(YMV::rank == 1 || YMV::rank == 2, "KokkosBlas::dot: Y must have rank 1 or 2.");
  static_assert(
      (XMV::rank == 2 && YMV::rank == 2 && RV::rank == 1) || (XMV::rank == 1 && YMV::rank == 1 && RV::rank == 0) ||
          (XMV::rank == 2 && YMV::rank == 1 && RV::rank == 1) || (XMV::rank == 1 && YMV::rank == 2 && RV::rank == 1),
      "KokkosBlas::dot: Ranks of RV, XMV, and YMV don't match.  "
      "See this function's documentation for the allowed "
      "combinations of ranks.");

  // Check compatibility of dimensions at run time.

  // Regardless of ranks of X and Y, their numbers of rows must match.
  bool dimsMatch = true;
  if (X.extent(0) != Y.extent(0)) {
    dimsMatch = false;
  } else if (X.extent(1) != Y.extent(1) && X.extent(1) != 1 && Y.extent(1) != 1) {
    // Numbers of columns don't match, and neither X nor Y have one column.
    dimsMatch = false;
  }
  const auto maxNumCols = X.extent(1) > Y.extent(1) ? X.extent(1) : Y.extent(1);
  if (RV::rank == 1 && R.extent(0) != maxNumCols) {
    dimsMatch = false;
  }

  if (!dimsMatch) {
    std::ostringstream os;
    os << "KokkosBlas::dot: Dimensions of R, X, and Y do not match: ";
    if (RV::rank == 1) {
      os << "R: " << R.extent(0) << " x " << X.extent(1) << ", ";
    }
    os << "X: " << X.extent(0) << " x " << X.extent(1) << ", Y: " << Y.extent(0) << " x " << Y.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Create unmanaged versions of the input Views.
  using UnifiedXLayout  = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedRVLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<RV, UnifiedXLayout>::array_layout;

  typedef Kokkos::View<typename std::conditional<RV::rank == 0, typename RV::non_const_value_type,
                                                 typename RV::non_const_value_type*>::type,
                       UnifiedRVLayout, typename RV::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      RV_Internal;
  typedef Kokkos::View<typename std::conditional<XMV::rank == 1, typename XMV::const_value_type*,
                                                 typename XMV::const_value_type**>::type,
                       UnifiedXLayout, typename XMV::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      XMV_Internal;
  typedef Kokkos::View<typename std::conditional<YMV::rank == 1, typename YMV::const_value_type*,
                                                 typename YMV::const_value_type**>::type,
                       typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout, typename YMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      YMV_Internal;

  RV_Internal R_internal  = R;
  XMV_Internal X_internal = X;
  YMV_Internal Y_internal = Y;

  Impl::Dot<execution_space, RV_Internal, XMV_Internal, YMV_Internal>::dot(space, R_internal, X_internal, Y_internal);
}

/// \brief Compute the column-wise dot products of two multivectors.
///
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVM.
///
/// \tparam RV 0-D resp. 1-D output View
/// \tparam XMV 1-D resp. 2-D input View
/// \tparam YMV 1-D resp. 2-D input View
///
/// \param R [out] Output 1-D or 0-D View to which to write results.
/// \param X [in] Input 2-D or 1-D View.
/// \param Y [in] Input 2-D or 1-D View.
///
/// This function implements a few different use cases:
/// <ul>
/// <li> If X and Y are both 1-D, then this is a single dot product.
///   R must be 0-D (a View of a single value). </li>
/// <li> If X and Y are both 2-D, then this function computes their
///   dot products columnwise.  R must be 1-D. </li>
/// <li> If X is 2-D and Y is 1-D, then this function computes the dot
///   product of each column of X, with Y, in turn.  R must be
///   1-D. </li>
/// <li> If X is 1-D and Y is 2-D, then this function computes the dot
///   product X with each column of Y, in turn.  R must be 1-D. </li>
/// </ul>
///
/// \note To implementers: We use enable_if here so that the compiler
///   doesn't confuse this version of dot() with the three-argument
///   version of dot() in Kokkos_Blas1.hpp.
template <class RV, class XMV, class YMV>
void dot(const RV& R, const XMV& X, const YMV& Y, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  dot(typename XMV::execution_space{}, R, X, Y);
}
}  // namespace KokkosBlas

#endif
