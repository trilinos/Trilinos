// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_CONFIGDEFS_HPP
#define TPETRA_CONFIGDEFS_HPP

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

#ifdef HAVE_MAP
#include <map>
#else
#include <map.h>
#endif // ifdef HAVE_MAP

#ifdef HAVE_VECTOR
#include <vector>
#else
#include <vector.h>
#endif // ifdef HAVE_VECTOR

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGO_H)
#include <algo.h>
#else
#include <algorithm.h>
#endif // ifdef HAVE_ALGORITHM

#ifdef HAVE_NUMERIC
#include <numeric>
#else
#include <algo.h>
#endif // ifdef HAVE_NUMERIC

/* 
this "using namespace" is a bad thing. (Not a Bad Thing, mind you, just a bad thing.) 
replace real soon now.
*/
using namespace std;

// other miscellaneous Tpetra configs:

#define TPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define TPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
// The TPETRA_SGN macro is a Bad Thing, since it uses floating-point literals.
//#define TPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)    /* sign function */

const int Tpetra_DefaultTracebackMode = 1; // Default value for traceback behavior

#endif // TPETRA_CONFIGDEFS_HPP
