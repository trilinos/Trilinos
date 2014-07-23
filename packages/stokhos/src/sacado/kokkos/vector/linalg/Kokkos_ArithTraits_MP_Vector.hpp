// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_ARITHTRAITS_MP_VECTOR_HPP
#define KOKKOS_ARITHTRAITS_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_ArithTraits.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::ArithTraits for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Details {

template <typename S>
class ArithTraits< Sacado::MP::Vector<S> > {
public:
  typedef Sacado::MP::Vector<S> val_type;

  typedef typename val_type::value_type base_value_type;
  typedef typename val_type::ordinal_type ordinal_type;
  typedef ArithTraits<base_value_type> BAT;

  typedef typename BAT::mag_type mag_type;
  //typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = BAT::is_signed;
  static const bool is_integer = BAT::is_integer;
  static const bool is_exact = BAT::is_exact;
  static const bool is_complex = BAT::is_complex;

  static KOKKOS_FORCEINLINE_FUNCTION bool isInf (const val_type& x) {
    bool res = false;
    for (ordinal_type i=0; i<x.size(); ++i)
      res = res || BAT::isInf(x.fastAccessCoeff(i));
    return res;
  }
  static KOKKOS_FORCEINLINE_FUNCTION bool isNan (const val_type& x) {
   bool res = false;
    for (ordinal_type i=0; i<x.size(); ++i)
      res = res || BAT::isInf(x.fastAccessCoeff(i));
    return res;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) {
    //return std::abs(x);
    const ordinal_type sz = x.size();
    mag_type n = mag_type(0);
    for (ordinal_type i=0; i<sz; ++i)
      n += BAT::abs( x.fastAccessCoeff(i) );
    return n;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type zero () {
    return val_type(0.0);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type one () {
    return val_type(1.0);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type min () {
    return BAT::min();
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type max () {
    return BAT::max();
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type real (const val_type& x) {
    const ordinal_type sz = x.size();
    val_type y(sz, base_value_type(0.0));
    for (ordinal_type i=0; i<sz; ++i)
      y.fastAccessCoeff(i) = BAT::real(x.fastAccessCoeff(i));
    return y;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type imag (const val_type& x) {
    const ordinal_type sz = x.size();
    val_type y(sz, base_value_type(0.0));
    for (ordinal_type i=0; i<sz; ++i)
      y.fastAccessCoeff(i) = BAT::imag(x.fastAccessCoeff(i));
    return y;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj (const val_type& x) {
    const ordinal_type sz = x.size();
    val_type y(sz, base_value_type(0.0));
    for (ordinal_type i=0; i<sz; ++i)
      y.fastAccessCoeff(i) = BAT::conj(x.fastAccessCoeff(i));
    return y;
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type pow (const val_type& x,
                                              const val_type& y) {
    return std::pow(x, y);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt (const val_type& x) {
    return std::sqrt(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log (const val_type& x) {
    return std::log(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type log10 (const val_type& x) {
    return std::log10(x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type nan () {
    return BAT::nan();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon () {
    return BAT::epsilon();
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef typename BAT::halfPrecision base_half_precision;
  typedef typename BAT::doublePrecision base_double_precision;
  typedef typename Sacado::mpl::apply<S,ordinal_type,base_half_precision>::type half_storage;
  typedef typename Sacado::mpl::apply<S,ordinal_type,base_double_precision>::type double_storage;
  typedef Sacado::MP::Vector<half_storage> halfPrecision;
  typedef Sacado::MP::Vector<double_storage> doublePrecision;
  static const bool isComplex = is_complex;
  static const bool isOrdinal = is_integer;
  static const bool isComparable = BAT::isComparable;
  static const bool hasMachineParameters = BAT::hasMachineParameters;
  static bool isnaninf (const val_type& x) {
    return isNan (x) || isInf (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type magnitude (const val_type& x) {
    return abs (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate (const val_type& x) {
    return conj (x);
  }
  static std::string name () {
    return Sacado::StringName<val_type>::eval();
  }
  static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot (const val_type& x) {
    return sqrt (x);
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type eps () {
    return epsilon ();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type sfmin () {
    return BAT::sfmin();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int base () {
    return BAT::base();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type prec () {
    return BAT::prec();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int t () {
    return BAT::t();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rnd () {
    return BAT::rnd();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emin () {
    return BAT::emin();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmin () {
    return BAT::rmin();
  }
  static KOKKOS_FORCEINLINE_FUNCTION int emax () {
    return BAT::emax();
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type rmax () {
    return BAT::rmax();
  }
};

}
}

#endif /* #ifndef KOKKOS_ARITHTRAITS_MP_VECTOR_HPP */
