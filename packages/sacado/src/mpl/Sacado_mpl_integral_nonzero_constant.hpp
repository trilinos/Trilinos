// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
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
      KOKKOS_INLINE_FUNCTION integral_nonzero_constant( const T & ) {}
    };

    template< typename T , T zero >
    struct integral_nonzero_constant<T,zero,false>
    {
      const T value ;
      typedef T value_type ;
      typedef integral_nonzero_constant<T,0> type ;
      KOKKOS_INLINE_FUNCTION integral_nonzero_constant( const T & v ) :
        value(v) {}
    };

  }

}

#endif // SACADO_MPL_INTEGRAL_NONZERO_CONSTANT_HPP
