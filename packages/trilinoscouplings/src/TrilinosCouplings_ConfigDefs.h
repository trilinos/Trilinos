#ifndef TRILINOSCOUPLINGS_CONFIGDEFS_H
#define TRILINOSCOUPLINGS_CONFIGDEFS_H

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

#include <Trilinos_Couplings_config.h>

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
using std::endl;
#endif
#endif

// Define DefultTracebackMode (HAVE_WARNING_MESSAGES and HAVE_FATAL_MESSAGES can be defined
// via the configure script command line)

#endif /* ifndef TRILINOSCOUPLINGS_CONFIGDEFS_H */
