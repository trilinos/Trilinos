// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_CONFIGDEFS_HPP_
#define _TEUCHOS_CONFIGDEFS_HPP_

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

#ifdef HAVE_CONFIG_H
#include "Teuchos_config.h"
#endif

// configs from ScalarTraits

//  If we're using a Sun compiler, version earlier than 5.0,
//  then complex isn't available.
#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define NO_COMPLEX
#endif

//  If we're using the tflops Portland Group compiler, then complex isn't
//  available. (As of July 21, 2000. abw)
#if defined(__PGI) && defined(__i386)
#define NO_COMPLEX
#endif

#if defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500
#define TYPENAME
#else
#define TYPENAME typename
#endif

// end of ScalarTraits configs

#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif
  
#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>   
#endif

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>   
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>   
#endif

#ifdef HAVE_VECTOR
#include <vector>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#endif

#ifdef HAVE_STRING
#include <string>
#else 
#include <string.h>
#endif

#ifndef TFLOP
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>  
#endif
#else /* TFLOP defined */
#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>  
#endif
#ifdef HAVE_STRING
using std::string;
#endif
#ifdef HAVE_IOSTREAM
using std::istream;
using std::ostream;  
using std::cerr;
using std::cout;
using std::cerr;
using std::cout;
using std::endl;
using std::swap;  
#endif

#endif

#ifdef TEUCHOS_SIMULATE_BOOL

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

#define TEUCHOS_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define TEUCHOS_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define TEUCHOS_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)  /* sign function */


const double Teuchos_MinDouble = 1.0E-100;
const double Teuchos_MaxDouble = 1.0E+100;
const double Teuchos_Overflow = 1.79E308; // Used to test if equilibration should be done.
const double Teuchos_Underflow = 2.23E-308;

// Delete any previous definition of TEUCHOS_NO_ERROR_REPORTS

#ifdef TEUCHOS_CHK_ERR
#undef TEUCHOS_CHK_ERR
#endif
#ifdef TEUCHOS_CHK_PTR
#undef TEUCHOS_CHK_PTR
#endif
#ifdef TEUCHOS_CHK_REF
#undef TEUCHOS_CHK_REF
#endif

// Make error report silent by defining TEUCHOS_NO_ERROR_REPORTS

#define TEUCHOS_CHK_ERR(a) { int tpetra_err = a; if (tpetra_err != 0)  return(tpetra_err);}
#define TEUCHOS_CHK_PTR(a) { return(a);}
#define TEUCHOS_CHK_REF(a) { return(a);}
const int Teuchos_DefaultTracebackMode = 1; // Default value for traceback behavior

//} // namespace Teuchos

#endif // _TEUCHOS_CONFIGDEFS_HPP_
