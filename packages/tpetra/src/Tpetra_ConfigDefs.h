/*Paul
24-July-2002 Initial writeup. All system includes go here, and a few other things too.
06-August-2002 Changed to images (nothing changed).
09-Oct-2002 Changed name to Tpetra_Compiler_Directives.h
30-Oct-2002 Changed name to Tpetra_ConfigDefs.h (hopefully done renaming now).
04-Dec-2002 Moved configs out of Tpetra_ScalarTraits.h.
*/

#ifndef _TPETRA_CONFIGDEFS_H_
#define _TPETRA_CONFIGDEFS_H_

#ifndef __cplusplus
#define __cplusplus
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

#endif // end of SGI || SGI64 || SGI32 || CPLANT


#ifdef TPETRA_SIMULATE_BOOL

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

#define TPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define TPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define TPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)  /* sign function */


const double Tpetra_MinDouble = 1.0E-100;
const double Tpetra_MaxDouble = 1.0E+100;
const double Tpetra_Overflow = 1.79E308; // Used to test if equilibration should be done.
const double Tpetra_Underflow = 2.23E-308;

// RAB: 2002/1/25: Define TPETRA_ANSI_CPP as an argument to the compiler!
//#undef TPETRA_ANSI_CPP // Do not use ANSI/ISO C++ (curently just checked for I/O functions)
/*#ifdef TPETRA_ANSI_CPP
typedef std::ios_base::fmtflags Tpetra_fmtflags;
#else
typedef long int Tpetra_fmtflags;
#endif */

const bool Tpetra_FormatStdout = true; // Set true if the ostream << operator should format output

// Delete any previous definition of TPETRA_NO_ERROR_REPORTS

#ifdef TPETRA_CHK_ERR
#undef TPETRA_CHK_ERR
#endif
#ifdef TPETRA_CHK_PTR
#undef TPETRA_CHK_PTR
#endif
#ifdef TPETRA_CHK_REF
#undef TPETRA_CHK_REF
#endif

// Make error report silent by defining TPETRA_NO_ERROR_REPORTS

#define TPETRA_CHK_ERR(a) { int tpetra_err = a; if (tpetra_err != 0)  return(tpetra_err);}
#define TPETRA_CHK_PTR(a) { return(a);}
#define TPETRA_CHK_REF(a) { return(a);}
const int Tpetra_DefaultTracebackMode = 1; // Default value for traceback behavior

//} // namespace Tpetra

#endif // _TPETRA_CONFIGDEFS_H_
