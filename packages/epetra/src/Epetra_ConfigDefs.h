
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_CONFIGDEFS_H_
#define _EPETRA_CONFIGDEFS_H_

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

#ifdef HAVE_CONFIG_H
#include <Epetra_config.h>

#ifdef HAVE_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#endif

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>
#endif

#ifdef HAVE_STRING
#include <string>
#else
#include <string.h>
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

/* Every line that begins with 'using' should eventually be dependent
   on some check within the configure script */
   

#ifndef TFLOP
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif
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

/*-----------------------------------------------------------------------
Must refine the following up to #else HAVE_CONFIG_H is not defined
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

#else /*HAVE_CONFIG_H is not defined*/

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

#endif /*HAVE_CONFIG_H*/

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
                      if (epetra_err !=0) return(epetra_err);  }\
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

#endif /*_EPETRA_CONFIGDEFS_H_*/
