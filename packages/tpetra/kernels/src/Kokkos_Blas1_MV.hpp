/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_BLAS1_MV_HPP_
#define KOKKOS_BLAS1_MV_HPP_

#include <Kokkos_Blas1_MV_impl_abs.hpp>
#include <Kokkos_Blas1_MV_impl_axpby.hpp>
#include <Kokkos_Blas1_MV_impl_dot.hpp>
#include <Kokkos_Blas1_MV_impl_fill.hpp>
#include <Kokkos_Blas1_MV_impl_nrm1.hpp>
#include <Kokkos_Blas1_MV_impl_nrm2.hpp>
#include <Kokkos_Blas1_MV_impl_nrmInf.hpp>
#include <Kokkos_Blas1_MV_impl_recip.hpp>
#include <Kokkos_Blas1_MV_impl_scal.hpp>
#include <Kokkos_Blas1_MV_impl_update.hpp>

#ifdef KOKKOS_HAVE_CXX11
#  include <type_traits>
#endif // KOKKOS_HAVE_CXX11

// Define this if you want some of the methods below to print out
// information about the types that they use internally.  This is
// useful for figuring out whether they are using the right
// specialization of the kernel.
//
// #define TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO 1

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
#  include <cxxabi.h>

namespace {
  // Print the type of t in demangled (human-readable) form.  This is
  // useful for debugging whether KokkosBlas::axpby() is calling one
  // of the full specializations, or the generic version of the
  // kernel.  See the following link for information about the
  // demangling function we call below:
  //
  // https://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
  template<class T>
  std::string demangledTypeName (const T& t) {
    int status = 0;
    return abi::__cxa_demangle (typeid (t).name (), 0, 0, &status);
  }
}
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO


namespace KokkosBlas {

/// \brief Compute the column-wise dot products of two multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam YMV 2-D input View
///
/// \param dots [out] Output 1-D View to which to write results.
/// \param x [in] Input 2-D View.
/// \param y [in] Input 2-D View.
template<class RV, class XMV, class YMV>
void
dot (const RV& dots, const XMV& x, const YMV& y)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::dot (MultiVector): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::dot (MultiVector): "
                 "The first input argument x is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::dot (MultiVector): "
                 "The second input argument y is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::dot (MultiVector): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RV::rank == 1, "KokkosBlas::dot (MultiVector): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::dot (MultiVector): "
                 "The first input argument x must have rank 2.");
  static_assert (YMV::rank == 2, "KokkosBlas::dot (MultiVector): "
                 "The second input argument y must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // RV, XMV, and YMV must be Kokkos::View specializations.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<YMV>::value>::type YMVIsNotView;

  // RV must be nonconst (else it can't be an output argument).
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;

  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Dot_MultiVectorRanksDontMatch;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Dot_XMultiVectorRankNot2;
  typedef typename
    Kokkos::Impl::StaticAssert<YMV::rank == 2 >::type Blas1_Dot_YMultiVectorRankNot2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (x.dimension_0 () != y.dimension_0 () ||
      x.dimension_1 () != y.dimension_1 () ||
      dots.dimension_0 () != x.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::dot (MultiVector): Dimensions do not match: "
       << "dots: " << dots.dimension_0 () << " x 1"
       << ", x: " << x.dimension_0 () << " x " << x.dimension_1 ()
       << ", y: " << y.dimension_0 () << " x " << y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;
  typedef Kokkos::View<typename YMV::const_value_type**,
    typename YMV::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YMV::specialize> YMV_Internal;

  RV_Internal dots_i = dots;
  XMV_Internal x_i = x;
  YMV_Internal y_i = y;

  return Impl::Dot_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize,
    typename YMV_Internal::value_type**,
    typename YMV_Internal::array_layout,
    typename YMV_Internal::device_type,
    typename YMV_Internal::memory_traits,
    typename YMV_Internal::specialize
    >::dot (dots_i, x_i, y_i);
}

/// \brief Compute the dot product of X(:,X_col) and Y(:,Y_col), and
///   store result in r(0).
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam YMV 2-D input View
///
/// \param dots [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
/// \param X_col [in] Column of X to use.
/// \param Y [in] Input 2-D View.
/// \param Y_col [in] Column of Y to use.
template<class RV, class XMV, class YMV>
void
dot (const RV& dots, const size_t dots_col, const XMV& X, const size_t X_col, const YMV& Y, const size_t Y_col)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::dot (MV, 5 arg): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::dot (MV, 5 arg): "
                 "The first input argument x is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::dot (MV, 5 arg): "
                 "The second input argument y is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::dot (MV, 5 arg): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RV::rank == 1, "KokkosBlas::dot (MV, 5 arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::dot (MV, 5 arg): "
                 "The first input argument x must have rank 2.");
  static_assert (YMV::rank == 2, "KokkosBlas::dot (MV, 5 arg): "
                 "The second input argument y must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // RV, XMV, and YMV must be Kokkos::View specializations.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<YMV>::value>::type YMVIsNotView;

  // RV must be nonconst (else it can't be an output argument).
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;

  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Dot_MultiVectorRanksDontMatch;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Dot_XMultiVectorRankNot2;
  typedef typename
    Kokkos::Impl::StaticAssert<YMV::rank == 2 >::type Blas1_Dot_YMultiVectorRankNot2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != Y.dimension_0 () ||
      dots.dimension_0 () < typename RV::size_type (1)) {
    std::ostringstream os;
    os << "KokkosBlas::dot (MV, 5 arg): Dimensions do not match: "
       << "dots: " << dots.dimension_0 () << " x 1"
       << ", x: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", y: " << Y.dimension_0 () << " x " << Y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;
  typedef Kokkos::View<typename YMV::const_value_type**,
    typename YMV::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YMV::specialize> YMV_Internal;

  RV_Internal dots_i = dots;
  XMV_Internal X_internal = X;
  YMV_Internal Y_internal = Y;

  return Impl::Dot_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize,
    typename YMV_Internal::value_type**,
    typename YMV_Internal::array_layout,
    typename YMV_Internal::device_type,
    typename YMV_Internal::memory_traits,
    typename YMV_Internal::specialize
      >::dot (dots_i, dots_col, X_internal, X_col, Y_internal, Y_col);
}

/// \brief Fill the multivector (2-D View) X with the given value.
///
/// \tparam XMV 2-D output View
///
/// \param X [out] Output 2-D View.
/// \param val [in] Value with which to fill the entries of X.
template<class XMV>
void
fill (const XMV& X, const typename XMV::non_const_value_type& val)
{
#ifdef KOKKOS_HAVE_CXX11
  // XMV must be a Kokkos::View specialization.
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::fill (MultiVector): "
                 "The first argument X is not a Kokkos::View.");
  // XMV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename XMV::value_type, typename XMV::non_const_value_type>::value,
                 "KokkosBlas::fill (MultiVector): X is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  // XMV must have rank 2.
  static_assert (XMV::rank == 2, "KokkosBlas::fill (MultiVector): "
                 "The first argument X must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // XMV must be a Kokkos::View specialization.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMVIsNotView;
  // XMV must be nonconst (else it can't be an output argument).
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename XMV::value_type,
      typename XMV::non_const_value_type>::value>::type XMV_is_const;
  // XMV must have rank 2.
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type XMV_not_rank_2;
#endif // KOKKOS_HAVE_CXX11

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename XMV::non_const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

#ifdef KOKKOS_HAVE_CXX11
  // XMV_Internal must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename XMV_Internal::value_type,
                   typename XMV_Internal::non_const_value_type>::value,
                 "KokkosBlas::fill (MultiVector): XMV_Internal is const.  "
                 "Please report this bug to the Tpetra developers.");
  // XMV_Internal must have rank 2.
  static_assert (XMV_Internal::rank == 2, "KokkosBlas::fill (MultiVector): "
                 "KokkosBlas::fill (MultiVector): "
                 "XMV_Internal does not have rank 2.  "
                 "Please report this bug to the Tpetra developers.");
#endif // KOKKOS_HAVE_CXX11

  XMV_Internal X_internal = X;

  Impl::Fill_MV<
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
    >::fill (X_internal, val);
}


/// \brief Compute the squares of 2-norms of the columns of the
///   multivector (2-D View) X.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
template<class RV, class XMV>
void
nrm2_squared (const RV& norms, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm2_squared (MV, 2-arg): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm2_squared (MV, 2-arg): "
                 "The first input argument X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm2_squared (MV, 2-arg): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrm2_squared (MV, 2-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrm2_squared (MV, 2-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Nrm2_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Nrm2_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (norms.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrm2_squared (MV, 2-arg): Dimensions do not match: "
       << "norms: " << norms.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_i = norms;
  XMV_Internal X_internal = X;

  Impl::Nrm2_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrm2_squared (norms_i, X_internal);
}

/// \brief Compute the 2-norm of X(:,X_col), and store the result in
///   norms(norms_col).
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
/// \param X_col [i] Column of X to process.
template<class RV, class XMV>
void
nrm2_squared (const RV& norms, const size_t norms_col, const XMV& X, const size_t X_col)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm2_squared "
                 "(MV, 3-arg): The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm2_squared "
                 "(MV, 3-arg): The first input argument X is not a "
                 "Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm2_squared (MV, 3-arg): The output argument "
                 "is const.  It must be nonconst, because it is an output "
                 "argument (we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrm2_squared (MV, 3-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrm2_squared (MV, 3-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Nrm2_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Nrm2_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_internal = norms;
  XMV_Internal X_internal = X;

  Impl::Nrm2_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrm2_squared (norms_internal, norms_col, X_internal, X_col);
}

/// \brief Compute the 1-norms of the columns of the multivector (2-D
///   View) X.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
template<class RV, class XMV>
void
nrm1 (const RV& norms, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm1 (MV, 2-arg): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm1 (MV, 2-arg): "
                 "The first input argument X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm1 (MV, 2-arg): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrm1 (MV, 2-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrm1 (MV, 2-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Nrm1_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Nrm1_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (norms.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrm1 (MV, 2-arg): Dimensions do not match: "
       << "norms: " << norms.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_internal = norms;
  XMV_Internal X_internal = X;

  Impl::Nrm1_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrm1 (norms_internal, X_internal);
}


/// \brief Compute the 1-norm of X(:,X_col), and store the result in
///   norms(norms_col).
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
/// \param X_col [i] Column of X to process.
template<class RV, class XMV>
void
nrm1 (const RV& norms, const size_t norms_col, const XMV& X, const size_t X_col)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm1 "
                 "(MV, 3-arg): The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm1 "
                 "(MV, 3-arg): The first input argument X is not a "
                 "Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm1 (MV, 3-arg): The output argument "
                 "is const.  It must be nonconst, because it is an output "
                 "argument (we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrm1 (MV, 3-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrm1 (MV, 3-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Nrm1_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Nrm1_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_internal = norms;
  XMV_Internal X_internal = X;

  Impl::Nrm1_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrm1 (norms_internal, norms_col, X_internal, X_col);
}

/// \brief Compute the inf-norms of the columns of the multivector
///   (2-D View) X.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
template<class RV, class XMV>
void
nrmInf (const RV& norms, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrmInf (MV, 2-arg): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrmInf (MV, 2-arg): "
                 "The first input argument X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrmInf (MV, 2-arg): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrmInf (MV, 2-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrmInf (MV, 2-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_NrmInf_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_NrmInf_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (norms.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrmInf (MV, 2-arg): Dimensions do not match: "
       << "norms: " << norms.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_internal = norms;
  XMV_Internal X_internal = X;

  Impl::NrmInf_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrmInf (norms_internal, X_internal);
}

/// \brief Compute the inf-norm of X(:,X_col), and store the result in
///   norms(norms_col).
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
///
/// \param norms [out] Output 1-D View to which to write results.
/// \param X [in] Input 2-D View.
/// \param X_col [i] Column of X to process.
template<class RV, class XMV>
void
nrmInf (const RV& norms, const size_t norms_col, const XMV& X, const size_t X_col)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrmInf "
                 "(MV, 3-arg): The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrmInf "
                 "(MV, 3-arg): The first input argument X is not a "
                 "Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrmInf (MV, 3-arg): The output argument "
                 "is const.  It must be nonconst, because it is an output "
                 "argument (we have to be able to write to its entries).");
  // RV must have rank 1, and XMV must have rank 2.
  static_assert (RV::rank == 1, "KokkosBlas::nrmInf (MV, 3-arg): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::nrmInf (MV, 3-arg): "
                 "The first input argument x must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMV_is_not_Kokkos_View;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;
  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_NrmInf_RV_rank_not_1;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_NrmInf_XMV_rank_not_2;
#endif // KOKKOS_HAVE_CXX11

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal norms_internal = norms;
  XMV_Internal X_internal = X;

  Impl::NrmInf_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize
      >::nrmInf (norms_internal, norms_col, X_internal, X_col);
}


namespace { // (anonymous)

// axpby() accepts both scalar coefficients a and b, and vector
// coefficients (apply one for each column of the input multivectors).
// This traits class helps axpby() select the correct specialization
// of AV and BV (the type of a resp. b) for invoking the
// implementation.
template<class T, bool isView>
struct GetInternalTypeForAxpby {
  typedef T type;
};

template<class T>
struct GetInternalTypeForAxpby<T, true> {
  typedef Kokkos::View<typename T::const_value_type*,
                       typename T::array_layout,
                       typename T::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       typename T::specialize> type;
};

  // typedef typename Kokkos::Impl::if_c<
  //   Kokkos::Impl::is_view<BV>::value, // if BV is a Kokkos::View,
  //   Kokkos::View<typename BV::const_value_type*,
  //     typename BV::array_layout,
  //     typename BV::device_type,
  //     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
  //     typename BV::specialize>, // use the unmanaged version of BV,
  //   BV >::type BV_Internal; // else just use BV.

} // namespace (anonymous)


template<class RMV, class AV, class XMV, class BV, class YMV>
void
axpby (const RMV& R, const AV& a, const XMV& X, const BV& b, const YMV& Y)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::axpby (MV): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::axpby (MV): "
                 "The first input argument x is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::axpby (MV): "
                 "The second input argument y is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                 "KokkosBlas::axpby (MV): R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RMV::rank == XMV::rank, "KokkosBlas::axpby (MV): "
                 "R and X must have the same rank.");
  static_assert (RMV::rank == YMV::rank, "KokkosBlas::axpby (MV): "
                 "R and Y must have the same rank.");
  static_assert (RMV::rank == 1 || RMV::rank == 2, "KokkosBlas::axpby (MV): "
                 "RMV, XMV, and YMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != Y.dimension_0 () ||
      X.dimension_1 () != Y.dimension_1 () ||
      X.dimension_0 () != R.dimension_0 () ||
      X.dimension_1 () != R.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::axpby (MV): Dimensions of R, X, and Y do not match: "
       << "R: " << R.dimension_0 () << " x " << R.dimension_1 ()
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", Y: " << Y.dimension_0 () << " x " << Y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV, XMV, and YMV
  // may be rank 1 or rank 2.  AV and BV may be either rank-1 Views,
  // or scalar values.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RMV::rank == 1,
      typename RMV::non_const_value_type*,
      typename RMV::non_const_value_type** >::type,
    typename RMV::array_layout,
    typename RMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RMV::specialize> RMV_Internal;
  typedef typename GetInternalTypeForAxpby<
    AV, Kokkos::Impl::is_view<AV>::value>::type AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::non_const_value_type*,
      typename XMV::non_const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;
  typedef typename GetInternalTypeForAxpby<
    BV, Kokkos::Impl::is_view<BV>::value>::type BV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      YMV::rank == 1,
      typename YMV::non_const_value_type*,
      typename YMV::non_const_value_type** >::type,
    typename YMV::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YMV::specialize> YMV_Internal;

  RMV_Internal R_internal = R;
  AV_Internal  a_internal = a;
  XMV_Internal X_internal = X;
  BV_Internal  b_internal = b;
  YMV_Internal Y_internal = Y;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::axpby:" << endl
       << "  RMV_Internal: " << demangledTypeName (R_internal) << endl
       << "  AV_Internal: " << demangledTypeName (a_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << "  BV_Internal: " << demangledTypeName (b_internal) << endl
       << "  YMV_Internal: " << demangledTypeName (Y_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  return Impl::Axpby<RMV_Internal, AV_Internal, XMV_Internal,
    BV_Internal, YMV_Internal>::axpby (R_internal, a_internal, X_internal,
                                       b_internal, Y_internal);
}


template<class RMV, class AV, class XMV>
void
scal (const RMV& R, const AV& a, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::scal (MV): "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::scal (MV): "
                 "X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                 "KokkosBlas::scal (MV): R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RMV::rank == XMV::rank, "KokkosBlas::scal (MV): "
                 "R and X must have the same rank.");
  static_assert (RMV::rank == 1 || RMV::rank == 2, "KokkosBlas::scal (MV): "
                 "RMV and XMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != R.dimension_0 () ||
      X.dimension_1 () != R.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::scal (MV): Dimensions of R and X do not match: "
       << "R: " << R.dimension_0 () << " x " << R.dimension_1 ()
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.  AV may be either a rank-1 View, or a scalar
  // value.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RMV::rank == 1,
      typename RMV::non_const_value_type*,
      typename RMV::non_const_value_type** >::type,
    typename RMV::array_layout,
    typename RMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RMV::specialize> RMV_Internal;
  typedef typename GetInternalTypeForAxpby<
    AV, Kokkos::Impl::is_view<AV>::value>::type AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::non_const_value_type*,
      typename XMV::non_const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RMV_Internal R_internal = R;
  AV_Internal  a_internal = a;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::scal:" << endl
       << "  RMV_Internal: " << demangledTypeName (R_internal) << endl
       << "  AV_Internal: " << demangledTypeName (a_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Scal<RMV_Internal, AV_Internal, XMV_Internal>::scal (R_internal, a_internal, X_internal);
}


/// \brief R(i,j) = abs(X(i,j))
///
/// Replace each entry in R with the absolute value (magnitude) of the
/// corresponding entry in X.
template<class RMV, class XMV>
void
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RMV and XMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::abs (MV): "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::abs (MV): "
                 "X is not a Kokkos::View.");
  // RMV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                 "KokkosBlas::abs (MV): R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RMV::rank == XMV::rank, "KokkosBlas::abs (MV): "
                 "R and X must have the same rank.");
  static_assert (RMV::rank == 1 || RMV::rank == 2, "KokkosBlas::abs (MV): "
                 "RMV and XMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != R.dimension_0 () ||
      X.dimension_1 () != R.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::abs (MV): Dimensions of R and X do not match: "
       << "R: " << R.dimension_0 () << " x " << R.dimension_1 ()
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RMV::rank == 1,
      typename RMV::non_const_value_type*,
      typename RMV::non_const_value_type** >::type,
    typename RMV::array_layout,
    typename RMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RMV::specialize> RMV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::non_const_value_type*,
      typename XMV::non_const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RMV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::abs:" << endl
       << "  RMV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Abs<RMV_Internal, XMV_Internal>::abs (R_internal, X_internal);
}


//! Compute Z := alpha*X + beta*Y + gamma*Z.
template<class XMV, class YMV, class ZMV>
void
update (const typename XMV::non_const_value_type& alpha, const XMV& X,
        const typename YMV::non_const_value_type& beta, const YMV& Y,
        const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
{
#ifdef KOKKOS_HAVE_CXX11
  // XMV, YMV, and ZMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::update (MV): "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::update (MV): "
                 "Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<ZMV>::value, "KokkosBlas::update (MV): "
                 "Z is not a Kokkos::View.");

  // ZMV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename ZMV::value_type, typename ZMV::non_const_value_type>::value,
                 "KokkosBlas::update (MV): Z is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (ZMV::rank == XMV::rank, "KokkosBlas::update (MV): "
                 "X and Z must have the same rank.");
  static_assert (ZMV::rank == YMV::rank, "KokkosBlas::update (MV): "
                 "Y and Z must have the same rank.");
  static_assert (ZMV::rank == 1 || ZMV::rank == 2, "KokkosBlas::update (MV): "
                 "XMV, YMV, and ZMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != Y.dimension_0 () ||
      X.dimension_1 () != Y.dimension_1 () ||
      X.dimension_0 () != Z.dimension_0 () ||
      X.dimension_1 () != Z.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::update (MV): Dimensions of X, Y, and Z do not match: "
       << "Z: " << Z.dimension_0 () << " x " << Z.dimension_1 ()
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", Y: " << Y.dimension_0 () << " x " << Y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  XMV, YMV, and ZMV
  // may be rank 1 or rank 2, but they must all have the same rank.

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      YMV::rank == 1,
      typename YMV::const_value_type*,
      typename YMV::const_value_type** >::type,
    typename YMV::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YMV::specialize> YMV_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      ZMV::rank == 1,
      typename ZMV::non_const_value_type*,
      typename ZMV::non_const_value_type** >::type,
    typename ZMV::array_layout,
    typename ZMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename ZMV::specialize> ZMV_Internal;

  XMV_Internal X_internal = X;
  YMV_Internal Y_internal = Y;
  ZMV_Internal Z_internal = Z;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::update:" << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << "  YMV_Internal: " << demangledTypeName (Y_internal) << endl
       << "  ZMV_Internal: " << demangledTypeName (Z_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  return Impl::Update<XMV_Internal, YMV_Internal,
    ZMV_Internal>::update (alpha, X_internal, beta, Y_internal,
                           gamma, Z_internal);
}


/// \brief R(i,j) = 1 / X(i,j)
///
/// Replace each entry in R with the reciprocal of the corresponding
/// entry in X.
template<class RMV, class XMV>
void
reciprocal (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RMV and XMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::reciprocal (MV): "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::reciprocal (MV): "
                 "X is not a Kokkos::View.");
  // RMV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                 "KokkosBlas::reciprocal (MV): R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RMV::rank == XMV::rank, "KokkosBlas::reciprocal (MV): "
                 "R and X must have the same rank.");
  static_assert (RMV::rank == 1 || RMV::rank == 2, "KokkosBlas::reciprocal (MV): "
                 "RMV and XMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != R.dimension_0 () ||
      X.dimension_1 () != R.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::reciprocal (MV): Dimensions of R and X do not match: "
       << "R: " << R.dimension_0 () << " x " << R.dimension_1 ()
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RMV::rank == 1,
      typename RMV::non_const_value_type*,
      typename RMV::non_const_value_type** >::type,
    typename RMV::array_layout,
    typename RMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RMV::specialize> RMV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::non_const_value_type*,
      typename XMV::non_const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RMV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::reciprocal:" << endl
       << "  RMV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Reciprocal<RMV_Internal, XMV_Internal>::reciprocal (R_internal, X_internal);
}


} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_HPP_
