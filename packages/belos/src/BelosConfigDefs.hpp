//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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

#endif /* BELOS_CONFIGDEFS_HPP */
