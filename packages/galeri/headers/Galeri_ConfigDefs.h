
// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_CONFIGDEFS_H
#define GALERI_CONFIGDEFS_H

#ifndef __cplusplus
#define __cplusplus
#endif

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
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

#include <Galeri_config.h>

#ifdef HAVE_MPI
#ifdef HAVE_GALERI_EPETRA
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#endif
#endif

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <string>
#include <iostream>

/* Every line that begins with 'using' should eventually be dependent
   on some check within the configure script */

#include <cmath>

#include <iomanip>

using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

#endif /* GALERI_CONFIGDEFS_H */
