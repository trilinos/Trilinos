// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

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

#if defined(HAVE_SACADO_KOKKOS)
#include "Kokkos_Macros.hpp"

#ifndef SACADO_FUNCTION
#define SACADO_FUNCTION KOKKOS_FUNCTION
#endif

#ifndef SACADO_DEFAULTED_FUNCTION
#define SACADO_DEFAULTED_FUNCTION KOKKOS_DEFAULTED_FUNCTION
#endif

#ifndef SACADO_INLINE_FUNCTION
#define SACADO_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#endif

#ifndef SACADO_FORCEINLINE_FUNCTION
#define SACADO_FORCEINLINE_FUNCTION  KOKKOS_FORCEINLINE_FUNCTION
#endif

#else
/* Define them even if Kokkos isn't enabled */

#ifndef SACADO_FUNCTION
#define SACADO_FUNCTION /* */
#endif

#ifndef SACADO_DEFAULTED_FUNCTION
#define SACADO_DEFAULTED_FUNCTION /* */
#endif

#ifndef SACADO_INLINE_FUNCTION
#define SACADO_INLINE_FUNCTION inline
#endif

#ifndef SACADO_FORCEINLINE_FUNCTION
#define SACADO_FORCEINLINE_FUNCTION  inline
#endif

#endif

/* Determine if the new fad design is supported.  Requies C++11,
   and if gcc, version 4.8 or greater.
*/
#if defined(__GNUC__) && !defined(__clang__)
#  if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8) )
#    define SACADO_ENABLE_NEW_DESIGN 1
#  endif
#else
#  define SACADO_ENABLE_NEW_DESIGN 1
#endif

#endif /* SACADO_CONFIGDEFS_H */
