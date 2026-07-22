// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CONFIGDEFS_H
#define STOKHOS_CONFIGDEFS_H

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
#include <Stokhos_config.h>
#endif

// Whether to use the Static-Fixed-Storage specialization of Sacado::MP::Vector,
// which does not use expression templates (default is to not)
#define STOKHOS_USE_MP_VECTOR_SFS_SPEC 0

#endif /* STOKHOS_CONFIGDEFS_H */
