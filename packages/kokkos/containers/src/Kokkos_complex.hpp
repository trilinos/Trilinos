/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_COMPLEX_HPP
#define KOKKOS_COMPLEX_HPP

#include <Kokkos_Core.hpp>
#include <iostream>

namespace Kokkos {

/// \class complex
/// \tparam RealType
///
/// Subset of std::complex, that works as the result of a
/// Kokkos::parallel_reduce.
template<class RealType>
class complex {
private:
  RealType re_, im_;

public:
  typedef RealType value_type;

  KOKKOS_INLINE_FUNCTION complex () :
    re_ (0.0), im_ (0.0)
  {}

  KOKKOS_INLINE_FUNCTION complex (const complex<RealType>& src) :
    re_ (src.real ()), im_ (src.imag ())
  {}

  KOKKOS_INLINE_FUNCTION complex (const RealType& val) :
    re_ (val), im_ (0.0)
  {}

  template<class RealType1, class RealType2>
  KOKKOS_INLINE_FUNCTION complex (const RealType1& re, const RealType2& im) :
    re_ (re), im_ (im)
  {}

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const complex<RealType>& src) {
    re_ = src.real ();
    im_ = src.imag ();
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const RealType& val) {
    re_ = val;
    im_ = static_cast<RealType> (0.0);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const int& val) {
    re_ = static_cast<RealType> (val);
    im_ = static_cast<RealType> (0.0);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION RealType imag () const {
    return im_;
  }

  KOKKOS_INLINE_FUNCTION RealType real () const {
    return re_;
  }

  KOKKOS_INLINE_FUNCTION RealType imag () const volatile {
    return im_;
  }

  KOKKOS_INLINE_FUNCTION RealType real () const volatile {
    return re_;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator += (const complex<RealType>& src) {
    re_ += src.real ();
    im_ += src.imag ();
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile complex<RealType>& src) volatile {
    re_ += src.real ();
    im_ += src.imag ();
  }

  KOKKOS_INLINE_FUNCTION operator RealType () {
    return re_;
  }
};


template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator + (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () + y.real (), x.imag () + y.imag ());
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator - (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () - y.real (), x.imag () - y.imag ());
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator - (const complex<RealType>& x) {
  return complex<RealType> (-x.real (), -x.imag ());
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator * (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () * y.real () - x.imag () * y.imag (),
                            x.real () * y.imag () + x.imag () * y.real ());
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
RealType abs (const complex<RealType>& x) {
  return x.real () * x.real () + x.imag () * x.imag ();
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType> conj (const complex<RealType>& x) {
  return complex<RealType> (x.real (), -x.imag ());
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator / (const complex<RealType>& x, const complex<RealType>& y) {
  const RealType abs_y = abs (y);
  const complex<RealType> y_conj_scaled (y.real () / abs_y, -y.imag () / abs_y);
  return x * y_conj_scaled;
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator == (const complex<RealType>& x, const complex<RealType> y) {
  return x.real () == y.real () && x.imag () == y.imag ();
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const complex<RealType>& x, const complex<RealType> y) {
  return x.real () != y.real () || x.imag () != y.imag ();
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const complex<RealType>& x, const RealType& y) {
  return x.real () != y || x.imag () != static_cast<RealType> (0.0);
}

template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const RealType& x, const complex<RealType>& y) {
  return y != x;
}

template<class RealType>
std::ostream& operator << (std::ostream& os, const complex<RealType>& x) {
  os << x.real () << (x.imag () < 0.0 ? " - " : " + ") << x.imag () << "i";
  return os;
}

} // namespace Kokkos

#endif // KOKKOS_COMPLEX_HPP
