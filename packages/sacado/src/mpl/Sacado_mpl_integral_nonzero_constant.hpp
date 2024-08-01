// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_INTEGRAL_NONZERO_CONSTANT_HPP
#define SACADO_MPL_INTEGRAL_NONZERO_CONSTANT_HPP

#include "Sacado_ConfigDefs.h"

namespace Sacado {

  namespace mpl {

    template< typename T , T v , bool NonZero = ( v != T(0) ) >
    struct integral_nonzero_constant
    {
      enum { value = T(v) };
      typedef T value_type ;
      typedef integral_nonzero_constant<T,v> type ;
      SACADO_INLINE_FUNCTION integral_nonzero_constant( const T & ) {}
    };

    template< typename T , T zero >
    struct integral_nonzero_constant<T,zero,false>
    {
      const T value ;
      typedef T value_type ;
      typedef integral_nonzero_constant<T,0> type ;
      SACADO_INLINE_FUNCTION integral_nonzero_constant( const T & v ) :
        value(v) {}
    };

  }

}

#endif // SACADO_MPL_INTEGRAL_NONZERO_CONSTANT_HPP
