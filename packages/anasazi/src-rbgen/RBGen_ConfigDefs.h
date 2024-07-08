// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_CONFIGDEFS_H
#define RBGEN_CONFIGDEFS_H

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

#include <Anasazi_config.h>

#ifdef HAVE_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#endif

#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <cmath>

// Define DefultTracebackMode (HAVE_WARNING_MESSAGES and HAVE_FATAL_MESSAGES can be defined
// via the configure script command line)

#endif /* RBGEN_CONFIGDEFS_H */
