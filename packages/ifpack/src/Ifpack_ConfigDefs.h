/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef _IFPACK_CONFIGDEFS_H_
#define _IFPACK_CONFIGDEFS_H_

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
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

#include <Ifpack_config.h>

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

#ifdef HAVE_MPI

#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif

#endif

#include <cstdio>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

// prints out an error message if variable is not zero,
// and returns this value.
#define IFPACK_CHK_ERR(ifpack_err) \
{ if (ifpack_err < 0) { \
  std::cerr << "IFPACK ERROR " << ifpack_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return(ifpack_err);  } }

// prints out an error message if variable is not zero,
// and returns void
#define IFPACK_CHK_ERRV(ifpack_err) \
{ if (ifpack_err < 0) { \
  std::cerr << "IFPACK ERROR " << ifpack_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return;  } }
// prints out an error message if variable is not zero,
// and returns false
#define IFPACK_CHK_ERRB(ifpack_err) \
{ if (ifpack_err < 0) { \
  std::cerr << "IFPACK ERROR " << ifpack_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return false;  } }
// prints out an error message and returns
#define IFPACK_RETURN(ifpack_err) \
{ if (ifpack_err < 0) { \
  std::cerr << "IFPACK ERROR " << ifpack_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
                       } return(ifpack_err); }
// prints out an error message if its *global* variable is not zero,
// and returns this *global* value.
#define IFPACK_CHK_GLOBAL_ERR(ifpack_err) \
{ int local_err = ifpack_err; \
  int global_min_err = 0; \
  Comm().MinAll(&local_err, &global_min_err, 1); \
  if (global_min_err < 0) { \
    return global_min_err; \
  } }

#define IFPACK_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)  /* sign function */
#define IFPACK_ABS(x) (((x) > 0.0) ? (x) : (-x))  /* abs function */

#endif /*_IFPACK_CONFIGDEFS_H_*/
