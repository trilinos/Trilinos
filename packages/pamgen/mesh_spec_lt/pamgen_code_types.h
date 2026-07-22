// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* note this header file must use C, not C++, syntax */

#ifndef code_typesH
#define code_typesH

#include <float.h> 
#include <limits.h>

/* typedefs to allow easy modification of word size for a particular
 * application or computer */

#define REAL_MAX   DBL_MAX
#define REAL_MIN   DBL_MIN
#define REAL_EPSILON DBL_EPSILON
/* These are useful if any of your variables will end up in the single precision exodus file and set to some max. */
#define FLOAT_MAX   FLT_MAX
#define FLOAT_MIN   FLT_MIN
#define FLOAT_EPSILON FLT_EPSILON

typedef unsigned int        uint;

#ifdef CRAY
typedef float                Real;
#elif defined MAC_SEMANTIC
typedef float                Real;
#else
typedef double                Real;
#endif

/* namespace macros to use */

/* this typedef is used to declare strings that can never be changed, e.g. 
 * message strings */

typedef const char *const        ConstString;

#ifndef NULL
#define NULL 0
#endif

#endif
