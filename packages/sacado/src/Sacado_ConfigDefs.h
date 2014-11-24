/*
@HEADER
*************************************************************************

                         Sacado Package
               Copyright (2006) Sandia Corporation

Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
the U.S. Government retains certain rights in this software.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
USA
Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
(etphipp@sandia.gov).

************************************************************************
@HEADER
*/

#ifndef SACADO_CONFIGDEFS_H
#define SACADO_CONFIGDEFS_H

#ifndef __cplusplus
#define __cplusplus
#endif

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and
 * need to be undef'd here to avoid warnings when this file is included from
 * another package.
 * KL 11/25/02
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#ifndef TRILINOS_NO_CONFIG_H
#include <Sacado_config.h>
#endif

/* Kokkos macros */

#ifdef HAVE_SACADO_KOKKOSCORE
#include "Kokkos_Macros.hpp"
#endif

/* Define them even if Kokkos isn't enabled */

#ifndef KOKKOS_FUNCTION
#define KOKKOS_FUNCTION /* */
#endif

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#ifndef KOKKOS_FORCEINLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION  inline
#endif

#include "Sacado_mpl_enable_if.hpp"
#include "Sacado_mpl_is_same.hpp"
#include "Sacado_mpl_is_convertible.hpp"

/* Define some macros useful for disabling template function overloads */
#define SACADO_ENABLE_IF_SAME(TYPE1, TYPE2, RETURN_TYPE)              \
  typename Sacado::mpl::enable_if<Sacado::mpl::is_same< TYPE1 , TYPE2 >, RETURN_TYPE >::type
#define SACADO_ENABLE_EXPR_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_SAME(typename Expr<S>::value_type, value_type, RETURN_TYPE)
#define SACADO_ENABLE_EXPR_CTOR_DEF SACADO_ENABLE_EXPR_FUNC(void*)
#define SACADO_ENABLE_EXPR_CTOR_DECL SACADO_ENABLE_EXPR_CTOR_DEF = 0

#define SACADO_ENABLE_IF_CONVERTIBLE(TYPE1, TYPE2, RETURN_TYPE)              \
  typename Sacado::mpl::enable_if<Sacado::mpl::is_convertible< TYPE1 , TYPE2 >, RETURN_TYPE >::type
#define SACADO_ENABLE_VALUE_FUNC(RETURN_TYPE) \
  SACADO_ENABLE_IF_CONVERTIBLE(S, value_type, RETURN_TYPE)
#define SACADO_ENABLE_VALUE_CTOR_DEF SACADO_ENABLE_VALUE_FUNC(void*)
#define SACADO_ENABLE_VALUE_CTOR_DECL SACADO_ENABLE_VALUE_CTOR_DEF = 0

#endif /* SACADO_CONFIGDEFS_H */
