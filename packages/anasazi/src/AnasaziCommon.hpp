/*  This script was originated in Epetra (I believe) and is useful for defining a 
common header file when the software package may not have been autoconfigured.
-- HKT  07/14/03
*/

#ifndef ANASAZI_COMMON_HPP
#define ANASAZI_COMMON_HPP

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

#include <Anasazi_config.h>

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

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif
using namespace std;

#else /*HAVE_CONFIG_H is not defined*/

#include <iostream>
#include <string>

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT) || defined (TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#endif

#include <vector>
#include <map>
#include <deque>
#include <algorithm>

using namespace std;

#endif /*HAVE_CONFIG_H*/

#endif /*ANASAZI_COMMON_HPP*/
