#ifndef _TPETRA_CONFIGDEFS_HPP_
#define _TPETRA_CONFIGDEFS_HPP_

#ifndef __cplusplus
#define __cplusplus
#endif

#ifdef HAVE_CONFIG_H

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

#include <Tpetra_config.h>

#ifdef HAVE_MPI
#ifndef TPETRA_MPI
#define TPETRA_MPI
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


#else /* HAVE_CONFIG_H is not defined */

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

#endif // end of HAVE_CONFIG conditional

#define TPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define TPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define TPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)    /* sign function */

const double Tpetra_MinDouble = 1.0E-100;
const double Tpetra_MaxDouble = 1.0E+100;
const double Tpetra_Overflow = 1.79E308; // Used to test if equilibration should be done.
const double Tpetra_Underflow = 2.23E-308;

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

#endif // end of TPETRA_SIMULATE_BOOL

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

#endif // _TPETRA_CONFIGDEFS_HPP_
