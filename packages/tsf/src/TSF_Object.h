#ifndef _TSF_OBJECT_H_
#define _TSF_OBJECT_H_


#ifndef __cplusplus
#define __cplusplus
#endif

#undef NDEBUG             // make sure asserts are enabled

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(SOLARIS)
#include <stdlib.h>
#include <assert.h>
#include <iostream.h>
#include <strstream.h>
#include <string.h>
#include <math.h>
#include <iomanip.h>

#else

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

#endif

#ifdef TSF_MPI
#include <mpi.h>
#endif

#ifdef TSF_SIMULATE_BOOL
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

#endif /* _TSF_OBJECT_H_ */
