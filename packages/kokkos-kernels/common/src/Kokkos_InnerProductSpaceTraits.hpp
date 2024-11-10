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

#ifndef KOKKOS_INNERPRODUCTSPACETRAITS_HPP
#define KOKKOS_INNERPRODUCTSPACETRAITS_HPP

/// \file Kokkos_InnerProductSpaceTraits.hpp
/// \brief Declaration and definition of
/// Kokkos::Details::InnerProductSpaceTraits

#include "Kokkos_ArithTraits.hpp"

namespace Kokkos {
namespace Details {

/// \class InnerProductSpaceTraits
/// \brief Traits class for inner product space operations on type T.
/// \tparam T "Scalar" type of interest
///
/// \warning This interface is not ready for exposure to users yet.
///   Users beware!
///
/// \warning This class is an implementation detail of Kokkos and
///   Tpetra.  You only need to deal with this class if you are
///   implementing a new Scalar type for Tpetra, or if you are using
///   Kokkos to implement your own linear algebra library.
///
/// This class defines typedefs and functions necessary for
/// implementing inner product space operations on finite-dimensional
/// vector spaces over T.  A <t>finite-dimensional vector space over
/// T</i> of dimension N is \f$T^N\f$, supplied with vector space
/// operations.  An example would be the space of length N vectors of
/// \c double; a more interesting (and motivating) example would be
/// the space of length N vectors of one of Sacado's automatic
/// differentiation (AD) types.
///
/// \section Kokkos_IPST_math Mathematical definition
///
/// An <i>inner product space</i> is a vector space that has an inner
/// product.  For example, we can define a "Euclidean" inner product
/// on \f$T^N\f$ as follows:
/// \[
///   dot(x,y) := \sum_k dot(x_k, y_k).
/// \]
/// If T is a real-valued type, then \f$dot(x_k, y_k)\f$ is just the
/// ordinary product of real numbers \f$x_k * y_k\f$.  However, we do
/// not require this.  For example, if T is one of Sacado's AD types,
/// we may wish to define \f$dot(x_k, y_k)\f$ to return a single
/// scalar value, rather than an AD value.  This is one way to use
/// vectors of AD values in a Krylov solver, for example.  (This is
/// also why we put "Euclidean" in quotes above.)
///
/// Inner product spaces have a norm induced by their inner product.
/// For example, given the above inner product definition, we can
/// define a "Euclidean" norm ("two-norm") on \f$T^N\f$ as follows:
/// \[
///   norm(x) := \sum_k norm(dot(x_k, x_k)).
/// \]
/// The result has the type of the absolute value of the product of
/// two T values.  We can express this as
/// <tt>InnerProductSpaceTraits<InnerProductSpaceTraits<T>::dot_type>::mag_type</tt>.
/// For example, if <tt>InnerProductSpaceTraits<T>::dot_type</tt> is
/// complex, then its \c mag_type is real.  The "norm" of a single
/// complex value here is just its absolute value (i.e., its
/// magnitude).
///
/// One could also define other kinds of norms using this class.  For
/// example, a "one-norm" analogue would look like this:
/// \[
///   norm(x,1) := \sum_k norm(x_k).
/// \]
/// Note the use of the "norm" of a single element of the vector, that
/// is, a single value of type T.  That could be a different function
/// than the "norm" that takes \c dot_type, which the above example
/// uses.
///
/// \section Kokkos_IPST_dev Use of class methods in device kernels
///
/// All class methods must be suitable for parallel kernels, if the
/// type T itself is suitable for parallel kernels.  In particular,
/// specializations for types T that make sense to use on a CUDA
/// device must mark all class methods as device (and host) functions,
/// using the KOKKOS_FORCEINLINE_FUNCTION macro.  All class methods must be
/// callable both inside and outside a parallel kernel (for CUDA, this
/// means they must be marked as both device and host functions).
///
/// \section Kokkos_IPST_unsupp Unsupported types on CUDA devices
///
/// CUDA does not support long double or std::complex<T> in device
/// functions.  This class does have specializations for these types,
/// but the class methods therein are not marked as device functions.
///
/// \section Kokkos_IPST_new Adding a specialization for a new type T
///
/// You must first add a specialization of Kokkos::ArithTraits<T>.  Please
/// note that if CUDA does not support using T in device functions,
/// then you must <i>not</t> mark norm() or dot() as device functions
/// in your specialization.  (Simply omit the KOKKOS_FORCEINLINE_FUNCTION
/// macro.)  If CUDA <i>does</i> support using T in device functions,
/// you <i>must</i> mark norm() and dot() as device functions in order
/// to use them in device functions.
template <class T>
class InnerProductSpaceTraits {
 public:
  //! The type T itself.
  typedef T val_type;

  //! The type returned by norm(x) for a value x of type val_type.
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;

  //! The type returned by dot(x,y) for values x and y of type val_type.
  typedef val_type dot_type;

  //! The "norm" (absolute value or magnitude) of a value x of type val_type.
  static KOKKOS_FORCEINLINE_FUNCTION mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  /// \brief The "dot product" of two values x and y of type val_type.
  ///
  /// This default implementation should suffice unless val_type is
  /// complex.  In that case, see the partial specialization for
  /// Kokkos::complex below to see our convention for which input gets
  /// conjugated.
  static KOKKOS_FORCEINLINE_FUNCTION dot_type dot(const val_type& x, const val_type& y) { return x * y; }
};

/// \brief Partial specialization for long double.
///
/// \warning CUDA does not support long double in device functions.
template <>
struct InnerProductSpaceTraits<long double> {
  typedef long double val_type;
  typedef Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static dot_type dot(const val_type& x, const val_type& y) { return x * y; }
};

//! Partial specialization for Kokkos::complex<T>.
template <class T>
class InnerProductSpaceTraits<Kokkos::complex<T>> {
 public:
  typedef Kokkos::complex<T> val_type;
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static KOKKOS_FORCEINLINE_FUNCTION mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static KOKKOS_FORCEINLINE_FUNCTION dot_type dot(const val_type& x, const val_type& y) { return Kokkos::conj(x) * y; }
};

/// \brief Partial specialization for std::complex<T>.
///
/// \warning CUDA does not support std::complex<T> in device
///   functions.
template <class T>
struct InnerProductSpaceTraits<std::complex<T>> {
  typedef std::complex<T> val_type;
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static dot_type dot(const val_type& x, const val_type& y) { return std::conj(x) * y; }
};

#ifdef HAVE_KOKKOSKERNELS_QUADMATH

/// \brief Partial specialization for __float128.
///
/// CUDA does not support __float128 in device functions, so none of
/// the class methods in this specialization are marked as device
/// functions.
template <>
struct InnerProductSpaceTraits<__float128> {
  typedef __float128 val_type;
  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static dot_type dot(const val_type& x, const val_type& y) { return x * y; }
};

#endif  // HAVE_KOKKOSKERNELS_QUADMATH

// dd_real and qd_real are floating-point types provided by the QD
// library of David Bailey (LBNL):
//
// http://crd-legacy.lbl.gov/~dhbailey/mpdist/
//
// dd_real uses two doubles (128 bits), and qd_real uses four doubles
// (256 bits).
//
// Kokkos does <i>not</i> currently support these types in device
// functions.  It should be possible to use Kokkos' support for
// aggregate types to implement device function support for dd_real
// and qd_real, but we have not done this yet (as of 07 Jan 2014).
// Hence, the class methods of the Kokkos::ArithTraits specializations for
// dd_real and qd_real are not marked as device functions.
#ifdef HAVE_KOKKOS_QD
template <>
struct InnerProductSpaceTraits<dd_real> {
  typedef dd_real val_type;
  typedef Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static dot_type dot(const val_type& x, const val_type& y) { return x * y; }
};

template <>
struct InnerProductSpaceTraits<qd_real> {
  typedef qd_real val_type;
  typedef Kokkos::ArithTraits<val_type>::mag_type mag_type;
  typedef val_type dot_type;

  static mag_type norm(const val_type& x) { return Kokkos::ArithTraits<val_type>::abs(x); }
  static dot_type dot(const val_type& x, const val_type& y) { return x * y; }
};
#endif  // HAVE_KOKKOS_QD

template <class ResultType, class InputType1, class InputType2>
KOKKOS_INLINE_FUNCTION void updateDot(ResultType& sum, const InputType1& x, const InputType2& y) {
  // FIXME (mfh 22 Jan 2020) We should actually pick the type with the
  // greater precision.
  sum += InnerProductSpaceTraits<InputType1>::dot(x, y);
}

KOKKOS_INLINE_FUNCTION void updateDot(double& sum, const double x, const double y) { sum += x * y; }

KOKKOS_INLINE_FUNCTION void updateDot(double& sum, const float x, const float y) { sum += x * y; }

// This exists because complex<float> += complex<double> is not defined.
KOKKOS_INLINE_FUNCTION void updateDot(Kokkos::complex<double>& sum, const Kokkos::complex<float> x,
                                      const Kokkos::complex<float> y) {
  const auto tmp = Kokkos::conj(x) * y;
  sum += Kokkos::complex<double>(tmp.real(), tmp.imag());
}

// This exists in case people call the overload of KokkosBlas::dot
// that takes an output View, and the output View has element type
// Kokkos::complex<float>.
KOKKOS_INLINE_FUNCTION void updateDot(Kokkos::complex<float>& sum, const Kokkos::complex<float> x,
                                      const Kokkos::complex<float> y) {
  sum += Kokkos::conj(x) * y;
}

// This exists because Kokkos::complex<double> =
// Kokkos::complex<float> is not defined.
template <class Out, class In>
struct CastPossiblyComplex {
  static Out cast(const In& x) { return x; }
};

template <class OutReal, class InReal>
struct CastPossiblyComplex<Kokkos::complex<OutReal>, Kokkos::complex<InReal>> {
  static Kokkos::complex<OutReal> cast(const Kokkos::complex<InReal>& x) {
    return {static_cast<OutReal>(x.real()), static_cast<OutReal>(x.imag())};
  }
};

}  // namespace Details
}  // namespace Kokkos

#endif  // KOKKOS_INNERPRODUCTSPACETRAITS_HPP
