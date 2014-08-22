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

#ifndef KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP
#define KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::InnerProductSpaceTraits for Sacado::MP::Vector
// scalar type
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Details {

template <typename S>
class InnerProductSpaceTraits< Sacado::MP::Vector<S> > {
public:
  typedef Sacado::MP::Vector<S> val_type;

  typedef typename val_type::value_type base_value_type;
  typedef typename val_type::ordinal_type ordinal_type;
  typedef InnerProductSpaceTraits<base_value_type> BIT;
  typedef typename BIT::dot_type base_dot_type;

  typedef typename ArithTraits<val_type>::mag_type mag_type;
  typedef base_dot_type dot_type;

  static KOKKOS_FORCEINLINE_FUNCTION
  mag_type norm (const val_type& x) {
    //return ArithTraits<val_type>::abs (x);
    const ordinal_type sz = x.size();
    mag_type nrm = mag_type(0);
    for (ordinal_type i=0; i<sz; ++i) {
      const mag_type n = BIT::norm( x.fastAccessCoeff(i) );
      nrm += n*n;
    }
    return std::sqrt(nrm);
  }

  static KOKKOS_FORCEINLINE_FUNCTION
  dot_type dot (const val_type& x, const val_type& y) {
    const ordinal_type xsz = x.size();
    const ordinal_type ysz = y.size();
    const ordinal_type sz = xsz > ysz ? xsz : ysz;

    dot_type r = dot_type(0);
    if (x.hasFastAccess(sz) && y.hasFastAccess(sz))
      for (ordinal_type i=0; i<sz; ++i)
        r += BIT::dot( x.fastAccessCoeff(i), y.fastAccessCoeff(i) );
    else
      for (ordinal_type i=0; i<sz; ++i)
        r += BIT::dot( x.coeff(i), y.coeff(i) );

    return r;
  }

};

template <typename S>
class InnerProductSpaceTraits< const Sacado::MP::Vector<S> > {
public:
  typedef Sacado::MP::Vector<S> val_type;

  typedef typename val_type::value_type base_value_type;
  typedef typename val_type::ordinal_type ordinal_type;
  typedef InnerProductSpaceTraits<base_value_type> BIT;
  typedef typename BIT::dot_type base_dot_type;

  typedef typename ArithTraits<val_type>::mag_type mag_type;
  typedef base_dot_type dot_type;

  static KOKKOS_FORCEINLINE_FUNCTION
  mag_type norm (const val_type& x) {
    //return ArithTraits<val_type>::abs (x);
    const ordinal_type sz = x.size();
    mag_type nrm = mag_type(0);
    for (ordinal_type i=0; i<sz; ++i) {
      const mag_type n = BIT::norm( x.fastAccessCoeff(i) );
      nrm += n*n;
    }
    return std::sqrt(nrm);
  }

  static KOKKOS_FORCEINLINE_FUNCTION
  dot_type dot (const val_type& x, const val_type& y) {
    const ordinal_type xsz = x.size();
    const ordinal_type ysz = y.size();
    const ordinal_type sz = xsz > ysz ? xsz : ysz;

    dot_type r = dot_type(0);
    if (x.hasFastAccess(sz) && y.hasFastAccess(sz))
      for (ordinal_type i=0; i<sz; ++i)
        r += BIT::dot( x.fastAccessCoeff(i), y.fastAccessCoeff(i) );
    else
      for (ordinal_type i=0; i<sz; ++i)
        r += BIT::dot( x.coeff(i), y.coeff(i) );

    return r;
  }

};

}
}

#endif /* #ifndef KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP */
