// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_CONFIGDEFS_HPP_
#define _TEUCHOS_CONFIGDEFS_HPP_

#ifndef __cplusplus
#define __cplusplus
#endif

#ifdef HAVE_CONFIG_H

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

#include "Teuchos_config.h"

/******************************************************************************
 *   Choose header file flavor: either ANSI-style (no .h, e.g. <iostream>) or
 * old-style (with .h, e.g., <iostream.h>). 
 * KL 9/26/03
 *****************************************************************************/

#if HAVE_CSTDIO
#include <cstdio>
#elif HAVE_STDIO_H
#include <stdio.h>
#else
#error "Found neither cstdio nor stdio.h"
#endif

#if HAVE_CSTDARG
#include <cstdarg>
#elif HAVE_STDARG_H
#include <stdarg.h>
#else
#error "Found neither cstdarg nor stdarg.h"
#endif

#if HAVE_CSTDLIB
#include <cstdlib>
#elif HAVE_STDLIB_H
#include <stdlib.h>
#else
#error "Found neither cstdlib nor stdlib.h"
#endif

#if HAVE_STRING
#include <string>
#elif HAVE_STRING_H
#include <string.h>
#else
#error "Found neither string nor string.h"
#endif

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#else
#error "Found neither iostream nor iostream.h"
#endif

#if HAVE_IOSTREAM
#include <fstream>
#elif HAVE_IOSTREAM_H
#include <fstream.h>
#else
#error "Found neither fstream nor fstream.h"
#endif

#if HAVE_STDEXCEPT
#include <stdexcept>
#elif HAVE_STDEXCEPT_H
#include <stdexcept.h>
#else
#error "Found neither stdexcept nor stdexcept.h"
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>   
#endif

#ifdef HAVE_COMPLEX
#include <complex>
#elif defined(HAVE_COMPLEX_H)
#include <complex.h>
#endif

#ifdef HAVE_VECTOR
#include <vector>
#include <deque> // RAB: 2003/11/10: If you have <vector> you should have <deque>
#elif defined(HAVE_VECTOR_H)
#include <vector.h>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGORITHM_H)
#include <algorithm.h>
#endif

#ifdef HAVE_MAP
#include <map>
#elif defined(HAVE_MAP_H)
#include <map.h>
#endif

#ifdef HAVE_LIST
#include <list>
#elif defined(HAVE_LIST_H)
#include <list.h>
#endif

#ifdef HAVE_STRING
#include <string>
#elif defined(HAVE_STRING_H) 
#include <string.h>
#elif defined(HAVE_STRINGS_H)
#include <strings.h>
#endif

#ifdef HAVE_TYPEINFO
#include <typeinfo>
#endif

/******************************************************************************
 * Choose string stream type: preferably std::ostringstream, otherwise
 * ostringstream or (gasp!) ostrstream. 
 *
 * Ross is going to write ANSI-compatible stringstreams to replace ostrstream
 * on our non-compliant platforms.
 *
 * KL 09/26/03
 ******************************************************************************/

#if HAVE_SSTREAM
#include <sstream>
typedef std::ostringstream TeuchosOStringStream;
#define TEUCHOS_OSTRINGSTREAM_GET_C_STR(OSS) (OSS).str().c_str()
#elif HAVE_SSTREAM_H
#include <sstream.h>
typedef ostringstream TeuchosOStringStream;
#define TEUCHOS_OSTRINGSTREAM_GET_C_STR(OSS) (OSS).str().c_str()
#elif HAVE_STRSTREAM
#include <strstream>
typedef std::ostrstream TeuchosOStringStream;
#define TEUCHOS_OSTRINGSTREAM_GET_C_STR(OSS) (OSS).str()
#elif HAVE_STRSTREAM_H
#include <strstream.h>
typedef ostrstream TeuchosOStringStream;
#define TEUCHOS_OSTRINGSTREAM_GET_C_STR(OSS) (OSS).str()
#else
#error "Found neither sstream, sstream.h, strstream.h, nor strstream"
#endif

#ifndef TFLOP
#if HAVE_CMATH
#include <cmath>
#elif HAVE_MATH_H
#include <math.h>
#else
#error "Found neither cmath nor math.h"
#endif
using namespace std;
#else /*TFLOP defined */
#ifdef HAVE_STRING
using std::string;
#endif
#ifdef HAVE_IOSTREAM
using std::istream;
using std::ostream;  
using std::cerr;
using std::cout;
using std::endl;
#endif
#ifdef HAVE_COMPLEX
using std::complex;
#endif
#endif /* TFLOP */

// RAB: 20031002: Added this for all platforms in addition to TFLOPS?
#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>  
#endif

#else /* fallback for the amazingly unlikely event we have no HAVE_CONFIG_H! */

#ifndef __cplusplus
#define __cplusplus
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <typeinfo>

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
using namespace std;

#elif defined(TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
using std::string;
#include <iomanip>
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
#include <list>

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
using namespace std;

#endif

#endif /* end HAVE_CONFIG_H */



/* Define bool in case we are running on an ancient, quirky, or merely
 * stupid compiler having no built-in bool.  WARNING: THIS IS EXTREMELY
 * DANGEROUS, BECAUSE OTHER CODES TO WHICH WE LINK MAY HACK THEIR BOOL
 * DEFINITIONS IN DIFFERENT WAYS, E.G. VIA A TYPEDEF. If that happens, 
 * a whole bunch of link errors will result. Avoid this hack if at all 
 * possible.
 */ 

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

/* 
 *
 */
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

#define TEUCHOS_CHK_ERR(a) { if (a != 0)  return(a);}
#define TEUCHOS_CHK_PTR(a) { return(a);}
#define TEUCHOS_CHK_REF(a) { return(a);}
const int Teuchos_DefaultTracebackMode = 1; // Default value for traceback behavior

#endif // _TEUCHOS_CONFIGDEFS_HPP_
