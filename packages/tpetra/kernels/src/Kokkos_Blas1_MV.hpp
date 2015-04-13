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
#include <Kokkos_Blas1_MV_impl_mult.hpp>
#include <Kokkos_Blas1_MV_impl_nrm1.hpp>
#include <Kokkos_Blas1_MV_impl_nrm2.hpp>
#include <Kokkos_Blas1_MV_impl_nrm2w.hpp>
#include <Kokkos_Blas1_MV_impl_nrmInf.hpp>
#include <Kokkos_Blas1_MV_impl_recip.hpp>
#include <Kokkos_Blas1_MV_impl_scal.hpp>
#include <Kokkos_Blas1_MV_impl_sum.hpp>
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

namespace Impl {
// parallel_for functor for computing the square root, in place, of a
// one-dimensional View.  This is useful for following the MPI
// all-reduce that computes the square of the two-norms of the local
// columns of a Tpetra::MultiVector.
//
// mfh 14 Jul 2014: Carter says that, for now, the standard idiom for
// operating on a single scalar value on the device, is to run in a
// parallel_for with N = 1.
//
// FIXME (mfh 14 Jul 2014): If we could assume C++11, this functor
// would go away.
template<class ViewType>
class SquareRootFunctor {
public:
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;

  SquareRootFunctor (const ViewType& theView) : theView_ (theView) {}

  KOKKOS_INLINE_FUNCTION void operator() (const size_type i) const {
    typedef typename ViewType::value_type value_type;
    theView_(i) = Kokkos::Details::ArithTraits<value_type>::sqrt (theView_(i));
  }

private:
  ViewType theView_;
};
}

/// \brief Compute the column-wise dot products of two multivectors.
///
/// \tparam RV 0-D resp. 1-D output View
/// \tparam XMV 1-D resp. 2-D input View
/// \tparam YMV 1-D resp. 2-D input View
///
/// \param dots [out] Output 1-D View to which to write results.
/// \param x [in] Input 2-D View.
/// \param y [in] Input 2-D View.
template<class RV, class XMV, class YMV>
void
dot (const RV& R, const XMV& X, const YMV& Y)
{
#ifdef KOKKOS_HAVE_CXX11
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
  static_assert (XMV::rank == YMV::rank, "KokkosBlas::dot: "
                 "X and Y must have the same rank.");
  static_assert ((RV::rank == 0 && XMV::rank == 1) ||
                 (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::dot: Either RV has rank 0 and XMV and YMV have "
                 "rank 1, or RV has rank 1 and XMV and YMV have rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != Y.dimension_0 () ||
      X.dimension_1 () != Y.dimension_1 () ||
      (RV::rank == 1 && R.dimension_0 () != X.dimension_1 ())) {
    std::ostringstream os;
    os << "KokkosBlas::dot: Dimensions of R, X, and Y do not match: ";
    if (RV::rank == 1) {
      os << "R: " << R.dimension_0 () << " x " << X.dimension_1 () << ", ";
    }
    os << "X: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", Y: " << Y.dimension_0 () << " x " << Y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
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

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;
  YMV_Internal Y_internal = Y;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::dot:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << "  YMV_Internal: " << demangledTypeName (Y_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Dot_MV<RV_Internal, XMV_Internal, YMV_Internal>::dot (R_internal, X_internal, Y_internal);
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

/// \brief Compute the squares of the 2-norm(s) of the column(s) of
///   the multivector or single vector X.
///
/// \tparam RV 0-D or 1-D output View
/// \tparam XMV 1-D or 2-D input View
///
/// \param R [out] Output View to which to write results.
/// \param X [in] Input View.
template<class RV, class XMV>
void
nrm2_squared (const RV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm2_squared: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm2_squared: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm2_squared: The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((RV::rank == 0 && XMV::rank == 1) ||
                 (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::nrm2_squared: Either RV has rank 0 and XMV has "
                 "rank 1, or RV has rank 1 and XMV has rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (RV::rank == 1 && R.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrm2_squared: Dimensions do not match: "
       << "R: " << R.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::nrm2_squared:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Nrm2_MV<RV_Internal, XMV_Internal>::nrm2_squared (R_internal, X_internal);
}


/// \brief Compute the 1-norm(s) of the column(s) of the multivector
///   or single vector X.
///
/// \tparam RV 0-D or 1-D output View
/// \tparam XMV 1-D or 2-D input View
///
/// \param R [out] Output View to which to write results.
/// \param X [in] Input View.
template<class RV, class XMV>
void
nrm1 (const RV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm1: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm1: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm1: The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((RV::rank == 0 && XMV::rank == 1) ||
                 (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::nrm1: Either RV has rank 0 and XMV has "
                 "rank 1, or RV has rank 1 and XMV has rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (RV::rank == 1 && R.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrm1: Dimensions do not match: "
       << "R: " << R.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::nrm1:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Nrm1_MV<RV_Internal, XMV_Internal>::nrm1 (R_internal, X_internal);
}

/// \brief Compute the inf-norm(s) of the column(s) of the multivector
///   or single vector X.
///
/// \tparam RV 0-D or 1-D output View
/// \tparam XMV 1-D or 2-D input View
///
/// \param R [out] Output View to which to write results.
/// \param X [in] Input View.
template<class RV, class XMV>
void
nrmInf (const RV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrmInf: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrmInf: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrmInf: The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((RV::rank == 0 && XMV::rank == 1) ||
                 (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::nrmInf: Either RV has rank 0 and XMV has "
                 "rank 1, or RV has rank 1 and XMV has rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (RV::rank == 1 && R.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrmInf: Dimensions do not match: "
       << "R: " << R.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::nrmInf:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::NrmInf_MV<RV_Internal, XMV_Internal>::nrmInf (R_internal, X_internal);
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


template<class AV, class XMV, class BV, class YMV>
void
axpby (const AV& a, const XMV& X, const BV& b, const YMV& Y)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::axpby (MV): "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::axpby (MV): "
                 "Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::axpby (MV): Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (YMV::rank == XMV::rank, "KokkosBlas::axpby (MV): "
                 "X and Y must have the same rank.");
  static_assert (YMV::rank == 1 || YMV::rank == 2, "KokkosBlas::axpby (MV): "
                 "XMV and YMV must either have rank 1 or rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != Y.dimension_0 () ||
      X.dimension_1 () != Y.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::axpby (MV): Dimensions of X and Y do not match: "
       << "X: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", Y: " << Y.dimension_0 () << " x " << Y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  XMV and YMV may be
  // rank 1 or rank 2.  AV and BV may be either rank-1 Views, or
  // scalar values.
  typedef typename GetInternalTypeForAxpby<
    AV, Kokkos::Impl::is_view<AV>::value>::type AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
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

  AV_Internal  a_internal = a;
  XMV_Internal X_internal = X;
  BV_Internal  b_internal = b;
  YMV_Internal Y_internal = Y;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::axpby:" << endl
       << "  AV_Internal: " << demangledTypeName (a_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << "  BV_Internal: " << demangledTypeName (b_internal) << endl
       << "  YMV_Internal: " << demangledTypeName (Y_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  return Impl::Axpby<AV_Internal, XMV_Internal, BV_Internal,
    YMV_Internal>::axpby (a_internal, X_internal, b_internal, Y_internal);
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

//! Compute sum of each column of X, and write to corresponding entry of R.
template<class RV, class XMV>
void
sum (const RV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  // RMV and XMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::sum: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::sum: "
                 "X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::sum: R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((RV::rank == 0 && XMV::rank == 1) || (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::sum: Ranks of R and X do not match.");
#endif // KOKKOS_HAVE_CXX11

  // TODO Check compatibility of dimensions at run time.

  // Create unmanaged versions of the input Views.  XMV may be rank 1
  // or rank 2, and RV may be rank 0 or rank 1.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::sum:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Sum<RV_Internal, XMV_Internal>::sum (R_internal, X_internal);
}

/// \brief Compute the "weighted 2-norm" of each column of X, using
///   weights in the corresponding column of W, and write the result
///   to the corresponding entry of R.
///
/// For single vectors X and W, the "weighted 2-norm" is the 2-norm of
/// the entry-wise quotient X(i) / W(i).  X and W have the same type.
template<class RV, class XMV>
void
nrm2w_squared (const RV& R, const XMV& X, const XMV& W)
{
#ifdef KOKKOS_HAVE_CXX11
  // RMV and XMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::nrm2w_squared: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::nrm2w_squared: "
                 "X is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::nrm2w_squared: R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((RV::rank == 0 && XMV::rank == 1) || (RV::rank == 1 && XMV::rank == 2),
                 "KokkosBlas::nrm2w_squared: Ranks of R and X do not match.  "
                 "If R has rank 0, X and W must have rank 1.  "
                 "If R has rank 1, X and W must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (X.dimension_0 () != W.dimension_0 () ||
      X.dimension_1 () != W.dimension_1 () ||
      R.dimension_0 () != X.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::nrm2w_squared: Dimensions do not match: "
       << "R: " << R.dimension_0 () << " x 1"
       << ", X: " << X.dimension_0 () << " x " << X.dimension_1 ()
       << ", W: " << W.dimension_0 () << " x " << W.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  XMV may be rank 1
  // or rank 2, and RV may be rank 0 or rank 1.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;
  XMV_Internal W_internal = W;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::nrm2w_squared:" << endl
       << "  RV_Internal: " << demangledTypeName (R_internal) << endl
       << "  XMV_Internal: " << demangledTypeName (X_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Nrm2w<RV_Internal, XMV_Internal>::nrm2w_squared (R_internal, X_internal, W_internal);
}


template<class CMV, class AV, class BMV>
void
mult (typename CMV::const_value_type& c,
      const CMV& C,
      typename AV::const_value_type& ab,
      const AV& A,
      const BMV& B)
{
#ifdef KOKKOS_HAVE_CXX11
  // CMV and AMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<CMV>::value, "KokkosBlas::mult: "
                 "C is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::mult: "
                 "A is not a Kokkos::View.");
  // CMV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename CMV::value_type,
                   typename CMV::non_const_value_type>::value,
                 "KokkosBlas::mult: C is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((BMV::rank == 1 && CMV::rank == 1) ||
                 (BMV::rank == 2 && CMV::rank == 2),
                 "KokkosBlas::mult: C and B must be either both rank 1, "
                 "or both rank 2.");
  static_assert (AV::rank == 1, "KokkosBlas::mult: A must have rank 1.");
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (C.dimension_0 () != A.dimension_0 () ||
      C.dimension_0 () != B.dimension_0 () ||
      C.dimension_1 () != B.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::mult: Dimensions do not match: "
       << "C: " << C.dimension_0 () << " x " << C.dimension_1 ()
       << ", A: " << A.dimension_0 () << " x " << A.dimension_0 ()
       << ", B: " << B.dimension_0 () << " x " << B.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      CMV::rank == 1,
      typename CMV::non_const_value_type*,
      typename CMV::non_const_value_type** >::type,
    typename CMV::array_layout,
    typename CMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename CMV::specialize> CMV_Internal;
  typedef Kokkos::View<
    typename AV::const_value_type*,
    typename AV::array_layout,
    typename AV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename AV::specialize> AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      BMV::rank == 1,
      typename BMV::const_value_type*,
      typename BMV::const_value_type** >::type,
    typename BMV::array_layout,
    typename BMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename BMV::specialize> BMV_Internal;

  CMV_Internal C_internal = C;
  AV_Internal A_internal = A;
  BMV_Internal B_internal = B;

#ifdef TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO
  using std::cerr;
  using std::endl;
  cerr << "KokkosBlas::mult:" << endl
       << "  CMV_Internal: " << demangledTypeName (C_internal) << endl
       << "  AV_Internal: " << demangledTypeName (A_internal) << endl
       << "  BMV_Internal: " << demangledTypeName (B_internal) << endl
       << endl;
#endif // TPETRAKERNELS_PRINT_DEMANGLED_TYPE_INFO

  Impl::Mult<CMV_Internal, AV_Internal, BMV_Internal>::mult (c, C_internal, ab,
                                                             A_internal, B_internal);
}

} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_HPP_
