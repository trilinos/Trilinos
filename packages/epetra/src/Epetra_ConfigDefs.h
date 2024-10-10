/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
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
*/

#ifndef EPETRA_CONFIGDEFS_H
#define EPETRA_CONFIGDEFS_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#ifndef __cplusplus
#define __cplusplus
#endif

#include <algorithm>
// windows defines these but the ones from algorithm
// are not macros, so this will let them be used
// from algorithm
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#define EPETRA_MAX(x,y) std::max(x,y) /* max function  */
#define EPETRA_MIN(x,y) std::min(x,y)/* min function  */
//#define EPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
//#define EPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define EPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)  /* sign function */

const double Epetra_MinDouble = 1.0E-100;
const double Epetra_MaxDouble = 1.0E+100;
const double Epetra_Overflow = 1.79E308; // Used to test if equilibration should be done.
const double Epetra_Underflow = 2.23E-308;

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

#include <Epetra_config.h>

#ifdef HAVE_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#endif

#include <cstdlib>

#include <cstdio>
using std::sprintf;
using std::sscanf;
using std::FILE;
using std::fopen;
using std::fclose;
using std::fgets;
using std::fprintf;

#include <cassert>

#include <cstring>
#include <string>

#include <iostream>

#include <sstream>

#include <cmath>
using std::rand;
using std::fabs;
using std::atoi;
using std::atof;
#ifndef __IBMCPP__
using std::abs;
#endif
using std::pow;
using std::sqrt;
using std::asin;
using std::sin;
using std::cos;
using std::ceil;
using std::floor;

#include <iomanip>

//using std::string;
using std::memcpy;
using std::strcpy;
using std::strcmp;
using std::strlen;
using std::strchr;
using std::strtok;
/*
For xlC 12.1, malloc and related functions are provided outside the std
namespace, so the below three using statements cause conflicting declarations.
*/
#ifndef __clang__
#if !defined __IBMCPP__ || ( __IBMCPP__ != 1210 )
using std::realloc;
using std::malloc;
using std::free;
#endif
#endif

//using std::istream;
//using std::ostream;
//using std::cerr;
//using std::cout;
//using std::endl;
//using std::flush;

using std::abort;
using std::exit;

/*-----------------------------------------------------------------------
  Must refine the following up to #else TRILINOS_NO_CONFIG_H is defined
  -----------------------------------------------------------------------*/

#ifdef EPETRA_SIMULATE_BOOL
#ifdef bool
#undef bool
#endif
#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif

#define bool int
#define true 1
#define false 0

#endif

#ifndef HAVE_FORMAT_IO
const bool Epetra_FormatStdout = false; // Set true if the ostream << operator should format output
#else
const bool Epetra_FormatStdout = true;
#endif

// Define DefultTracebackMode (HAVE_WARNING_MESSAGES and HAVE_FATAL_MESSAGES can be defined
// via the configure script command line)

#ifdef HAVE_WARNING_MESSAGES
const int DefaultTracebackMode = 2;
#elif defined HAVE_FATAL_MESSAGES
const int DefaultTracebackMode = 1;
#else
const int DefaultTracebackMode = 0;
#endif

#ifndef HAVE_FORTRAN_SUPPORT
#ifndef FORTRAN_DISABLED
#define FORTRAN_DISABLED
#endif
#endif

#else /*TRILINOS_NO_CONFIG_H is defined*/

#ifndef __cplusplus
#define __cplusplus
#endif

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string>

//using std::string;
//using std::istream;
//using std::ostream;
//using std::cerr;
//using std::cout;
//using std::endl;
//using std::flush;

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>

//using std::string;
//using std::istream;
//using std::ostream;
//using std::cerr;
//using std::cout;
//using std::endl;
//using std::flush;

#endif



#ifdef EPETRA_SIMULATE_BOOL
#ifdef bool
#undef bool
#endif
#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif

#define bool int
#define true 1
#define false 0

#endif

const bool Epetra_FormatStdout = true; // Set true if the ostream << operator should format output
const int DefaultTracebackMode = 1;

#endif /*TRILINOS_NO_CONFIG_H*/

// Delete any previous definition of EPETRA_NO_ERROR_REPORTS

#ifdef EPETRA_CHK_ERR
#undef EPETRA_CHK_ERR
#endif
#ifdef EPETRA_CHK_PTR
#undef EPETRA_CHK_PTR
#endif
#ifdef EPETRA_CHK_REF
#undef EPETRA_CHK_REF
#endif

// Great little macro obtained from Alan Williams (modified for dynamic switching on/off)

#define EPETRA_CHK_ERR(a) { { int epetra_err = a; \
                              if ((epetra_err < 0 && Epetra_Object::GetTracebackMode() > 0) || \
                                  (epetra_err > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
                      Epetra_Object::GetTracebackStream() << "Epetra ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; }\
                      if (epetra_err != 0) return(epetra_err);  }\
                   }

// Extension of same macro for pointer, returns zero if bad

#define EPETRA_CHK_PTR(a) { if (a == 0 && Epetra_Object::GetTracebackMode() > 0) { \
                      Epetra_Object::GetTracebackStream() << "Epetra returning zero pointer " << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; } \
                      return(a); \
                   }
// Extension of same macro for reference, returns a default reference

#define EPETRA_CHK_REF(a) { if (Epetra_Object::GetTracebackMode() > 0) {\
                      Epetra_Object::GetTracebackStream() << "Epetra returning default reference " << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; } \
                      return(a); \
                   }

#include "Epetra_DLLExportMacro.h"

#endif /* EPETRA_CONFIGDEFS_H */
