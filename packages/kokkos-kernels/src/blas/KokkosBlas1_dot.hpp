/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS1_DOT_HPP_
#define KOKKOSBLAS1_DOT_HPP_

#include<KokkosBlas1_dot_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace KokkosBlas {

/// \brief Return the dot product of the two vectors x and y.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
/// \tparam YVector Type of the second vector y; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
/// \param y [in] Input 1-D View.
///
/// \return The dot product result; a single value.
template<class XVector,class YVector>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
dot (const XVector& x, const YVector& y)
{
  static_assert (Kokkos::Impl::is_view<XVector>::value,
                 "KokkosBlas::dot: XVector must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YVector>::value,
                 "KokkosBlas::dot: YVector must be a Kokkos::View.");
  static_assert ((int) XVector::rank == (int) YVector::rank,
                 "KokkosBlas::dot: Vector ranks do not match.");
  static_assert (XVector::rank == 1, "KokkosBlas::dot: "
                 "Both Vector inputs must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (x.extent(0) != y.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::dot: Dimensions do not match: "
       << ", x: " << x.extent(0) << " x 1"
       << ", y: " << y.extent(0) << " x 1";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }


  typedef Kokkos::View<typename XVector::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
    typename XVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XVector_Internal;
  typedef Kokkos::View<typename YVector::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
    typename YVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_Internal;

  using dot_type =
    typename Kokkos::Details::InnerProductSpaceTraits<
      typename XVector::non_const_value_type>::dot_type;
  // Some platforms, such as Mac Clang, seem to get poor accuracy with
  // float and complex<float>.  Work around some Trilinos test
  // failures by using a higher-precision type for intermediate dot
  // product sums.
  constexpr bool is_complex_float =
    std::is_same<dot_type, Kokkos::complex<float>>::value;
  constexpr bool is_real_float = std::is_same<dot_type, float>::value;
  using result_type = typename std::conditional<is_complex_float,
    Kokkos::complex<double>,
    typename std::conditional<is_real_float,
      double,
      dot_type
      >::type
    >::type;
  using RVector_Internal = Kokkos::View<result_type,
    Kokkos::LayoutLeft,
    Kokkos::HostSpace,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  result_type result {};
  RVector_Internal R = RVector_Internal(&result);
  XVector_Internal X = x;
  YVector_Internal Y = y;

  Impl::Dot<RVector_Internal,XVector_Internal,YVector_Internal>::dot (R,X,Y);
  Kokkos::fence();
  // mfh 22 Jan 2020: We need the line below because
  // Kokkos::complex<T> lacks a constructor that takes a
  // Kokkos::complex<U> with U != T.
  return Kokkos::Details::CastPossiblyComplex<dot_type, result_type>::cast(result);
}

/// \brief Compute the column-wise dot products of two multivectors.
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
template<class RV, class XMV, class YMV>
void
dot (const RV& R, const XMV& X, const YMV& Y,
     typename std::enable_if<Kokkos::Impl::is_view<RV>::value, int>::type = 0)
{
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::dot: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::dot: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::dot: "
                 "Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::dot: R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RV::rank == 0 || RV::rank == 1,
                 "KokkosBlas::dot: R must have rank 0 or 1.");
  static_assert (XMV::rank == 1 || XMV::rank == 2,
                 "KokkosBlas::dot: X must have rank 1 or 2.");
  static_assert (YMV::rank == 1 || YMV::rank == 2,
                 "KokkosBlas::dot: Y must have rank 1 or 2.");
  static_assert ((XMV::rank == 2 && YMV::rank == 2 && RV::rank == 1) ||
                 (XMV::rank == 1 && YMV::rank == 1 && RV::rank == 0) ||
                 (XMV::rank == 2 && YMV::rank == 1 && RV::rank == 1) ||
                 (XMV::rank == 1 && YMV::rank == 2 && RV::rank == 1),
                 "KokkosBlas::dot: Ranks of RV, XMV, and YMV don't match.  "
                 "See this function's documentation for the allowed "
                 "combinations of ranks.");

  // Check compatibility of dimensions at run time.

  // Regardless of ranks of X and Y, their numbers of rows must match.
  bool dimsMatch = true;
  if (X.extent(0) != Y.extent(0)) {
    dimsMatch = false;
  }
  else if (X.extent(1) != Y.extent(1) &&
           X.extent(1) != 1 &&
           Y.extent(1) != 1) {
    // Numbers of columns don't match, and neither X nor Y have one column.
    dimsMatch = false;
  }
  const auto maxNumCols = X.extent(1) > Y.extent(1) ?
    X.extent(1) : Y.extent(1);
  if (RV::rank == 1 && R.extent(0) != maxNumCols) {
    dimsMatch = false;
  }

  if (! dimsMatch) {
    std::ostringstream os;
    os << "KokkosBlas::dot: Dimensions of R, X, and Y do not match: ";
    if (RV::rank == 1) {
      os << "R: " << R.extent(0) << " x " << X.extent(1) << ", ";
    }
    os << "X: " << X.extent(0) << " x " << X.extent(1)
       << ", Y: " << Y.extent(0) << " x " << Y.extent(1);
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<RV>::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      YMV::rank == 1,
      typename YMV::const_value_type*,
      typename YMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;
  YMV_Internal Y_internal = Y;

  Impl::Dot<RV_Internal, XMV_Internal, YMV_Internal>::dot(R_internal, X_internal, Y_internal);
}
}

#endif
