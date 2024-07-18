// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP
#define KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_STOKHOS_ENSEMBLE_REDUCT)

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

  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
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

  typedef typename Kokkos::ArithTraits<val_type>::mag_type mag_type;
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

#endif

#endif /* #ifndef KOKKOS_INNER_PRODUCT_SPACE_TRAITS_MP_VECTOR_HPP */
