// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_CONFIGDEFS_HPP
#define BELOS_CONFIGDEFS_HPP

/*! \file BelosConfigDefs.hpp
    \brief Belos header file which uses auto-configuration information to include necessary C++ headers.
*/

#ifndef __cplusplus
#define __cplusplus
#endif

#ifndef TRILINOS_NO_CONFIG_H

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

#include <Belos_config.h>

#ifdef HAVE_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#endif

#include "Teuchos_ConfigDefs.hpp"

#else

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT) || defined(TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>

#endif 

#endif /*ndef HAVE_CONFIG_H*/


/* Define some macros */
#define BELOS_MAX(x,y) (( (x) > (y) ) ? (x)  : (y) )     /* max function  */
#define BELOS_MIN(x,y) (( (x) < (y) ) ? (x)  : (y) )     /* min function  */
#define BELOS_SGN(x)   (( (x) < 0.0 ) ? -1.0 : 1.0 )     /* sign function */

namespace Belos { std::string Belos_Version(); }

// This include file defines macros to avoid warnings under CUDA.  See github issue #1133.
#include "Teuchos_CompilerCodeTweakMacros.hpp"

#endif /* BELOS_CONFIGDEFS_HPP */
