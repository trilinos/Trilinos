#ifndef _TPETRA_CONFIGDEFS_HPP_
#define _TPETRA_CONFIGDEFS_HPP_

#ifndef __cplusplus
#define __cplusplus
#endif // ifndef __cplusplus

/* this section undefines all the things autotools defines for us that we wish it didn't. */

#ifdef PACKAGE
#undef PACKAGE
#endif // ifdef PACKAGE

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif // ifdef PACKAGE_NAME

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif // ifdef PACKAGE_BUGREPORT

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif // ifdef PACKAGE_STRING

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif // ifdef PACKAGE_TARNAME

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif // ifdef PACKAGE_VERSION

#ifdef VERSION
#undef VERSION
#endif // ifdef VERSION

// end of undoing autoconf's work section

#include <Tpetra_config.h>

#ifdef HAVE_MPI
#ifndef TPETRA_MPI
#define TPETRA_MPI
#endif // ifndef TPETRA_MPI
#endif // ifdef HAVE_MPI

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif // ifdef HAVE_CSTDLIB

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif // ifdef HAVE_CSTDIO

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>
#endif // ifdef HAVE_CASSERT

#ifdef HAVE_STRING
#include <string>
#else
#include <string.h>
#endif // ifdef HAVE_STRING

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif // ifdef HAVE_IOSTREAM

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif // ifdef HAVE_CMATH

/* 
this "using namespace" is a bad thing. (Not a Bad Thing, mind you, just a bad thing.) 
replace real soon now.
*/
using namespace std;

// other miscellaneous Tpetra configs:

#define TPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define TPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define TPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)    /* sign function */

const int Tpetra_DefaultTracebackMode = 1; // Default value for traceback behavior

#endif // _TPETRA_CONFIGDEFS_HPP_
