#ifndef _PETRA_PETRA_H_
#define _PETRA_PETRA_H_

#undef PETRA_LEVELSCHEDULING

#ifndef __cplusplus
#define __cplusplus
#endif

#undef NDEBUG             // make sure asserts are enabled

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(SOLARIS) || defined(TFLOP) || defined(CPLANT)
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;

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

#define maxfn(x,y) (( x > y ) ? x : y)     /* max function  */
#define minfn(x,y) (( x < y ) ? x : y)     /* min function  */
#define sgnfn(x) ((x < 0.0) ? -1.0 : 1.0)  /* sign function */

#define Petra_name2(a,b) a ## b

#ifdef CRAY
#define F77NAME(x) x
#else
#define F77NAME(x) Petra_name2(x,_)
#endif

#ifdef PETRA_MPI
#include <mpi.h>
#endif

#ifdef PETRA_SIMULATE_BOOL
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

#define Petra_error(msg,code) \
{ cerr << "PETRA: " << msg << "    CODE " << code << endl; exit(1); }


const double Petra_MinDouble = 1.0E-100;
const double Petra_MaxDouble = 1.0E+100;

//! Data Access Mode Selector
/*! \enum Petra_DataAccess
    If set to Copy, user data will be copied at construction.
    If set to View, user data will be encapsulated and used throughout
    the life of the object.
*/
enum Petra_DataAccess {Copy, /*!< User data will be copied at
                                  construction. */
                       View /*!< User data will be encapsulated and
                                 used throughout the life of the object. */
                       };

//! Combine Mode Selector
/*! \enum Petra_CombineMode 
    If set to Add, components on the receiving processor will be added
    together.    If set to Zero, off-processor components will be ignored.
    If set to Insert, off-processor components will replace existing
    components on the receiving processor.
    If set to Average, off-processor components will be averaged with
    existing components on the receiving processor.
*/

enum Petra_CombineMode {Add,    /*!< Components on the receiving processor
                                     will be added together. */
                        Zero,   /*!< Off-processor components will be
                                     ignored. */
                        Insert, /*!< Off-processor components will
                                     be inserted into locations on
                                     receiving processor. */
                        Replace, /*!< Off-processor components will
                                     replace existing components on the 
                                     receiving processor. */
                        Average /*!< Off-processor components will be
                                     averaged with existing components 
                                     on the receiving processor. */
                        };


#endif /* _PETRA_PETRA_H_ */
