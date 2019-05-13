
//@HEADER
// ************************************************************************
//
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef EPETRA_CONFIGDEFS_H
#define EPETRA_CONFIGDEFS_H

#if !defined(TRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS)
#ifdef __GNUC__
#  warning "The package Claps is deprecated in April 2019; it will be removed from Trilinos in May 2019."
#endif
#endif

#ifndef __cplusplus
#define __cplusplus
#endif

#define EPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define EPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
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

#include <cassert>

#include <string>

#include <iostream>

/* Every line that begins with 'using' should eventually be dependent
   on some check within the configure script */


#ifndef TFLOP
#include <cmath>
using namespace std;
#else /* TFLOP defined */
#include <iomanip>
using std::string;
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
#endif

// Claps_Version() method
namespace Claps {
  string Claps_Version();
}

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

// RAB: 2002/1/25: Define EPETRA_ANSI_CPP as an argument to the compiler!
//#undef EPETRA_ANSI_CPP // Do not use ANSI/ISO C++ (curently just checked for I/O functions)
/*#ifdef EPETRA_ANSI_CPP
  typedef std::ios_base::fmtflags   Epetra_fmtflags;
  #else
  typedef long int                  Epetra_fmtflags;
  #endif*/

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
using namespace std;

#elif defined(TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
using std::string;
#include <iostream>
#include <iomanip>
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

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

// RAB: 2002/1/25: Define EPETRA_ANSI_CPP as an argument to the compiler!
//#undef EPETRA_ANSI_CPP // Do not use ANSI/ISO C++ (curently just checked for I/O functions)
/*#ifdef EPETRA_ANSI_CPP
  typedef std::ios_base::fmtflags   Epetra_fmtflags;
  #else
  typedef long int                  Epetra_fmtflags;
  #endif */

const bool Epetra_FormatStdout = true; // Set true if the ostream << operator should format output
const int DefaultTracebackMode = 1;

#endif /*ndef TRILINOS_NO_CONFIG_H*/

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
                      cerr << "Epetra ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; }\
                      if (epetra_err != 0) return(epetra_err);  }\
                   }

// Extension of same macro for pointer, returns zero if bad

#define EPETRA_CHK_PTR(a) { if (a == 0 && Epetra_Object::GetTracebackMode() > 0) { \
                      cerr << "Epetra returning zero pointer " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }
// Extension of same macro for reference, returns a default reference

#define EPETRA_CHK_REF(a) { if (Epetra_Object::GetTracebackMode() > 0) {\
                      cerr << "Epetra returning default reference " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }

#endif /* EPETRA_CONFIGDEFS_H */
