/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 1996 Sandia Corporation, Albuquerque, NM.

#ifndef stk_percept_code_types_hpp
#define stk_percept_code_types_hpp

#include <sstream>
#include <stdexcept>

#define THROW(message) do {                     \
    std::ostringstream msg;                     \
    msg <<  message;                            \
    throw std::runtime_error(msg.str());        \
  } while (0)



// If one of these two macros is not defined, default the precision to double
#if !defined FOUR_BYTE_REAL && !defined EIGHT_BYTE_REAL
#define EIGHT_BYTE_REAL
#endif

#if !defined FOUR_BYTE_LONG && !defined EIGHT_BYTE_LONG
#define FOUR_BYTE_LONG
#endif

/* typedefs to allow easy modification of word size for a particular
 * application or computer */

#include <limits.h>

typedef void		Void;

typedef int		Int;
typedef unsigned int	UInt;
#define Int_MAX		INT_MAX
#define Int_MIN		INT_MIN
#define UInt_MAX	UINT_MAX

#ifdef FOUR_BYTE_LONG

typedef long		Long;
typedef unsigned long	ULong;
#define Long_MAX	LONG_MAX
#define Long_MIN	LONG_MIN
#define ULong_MAX	ULONG_MAX

#elif defined EIGHT_BYTE_LONG

typedef long long		Long;
typedef unsigned long long	ULong;
#define Long_MAX	LLONG_MAX
#define Long_MIN	LLONG_MIN
#define ULong_MAX	ULLONG_MAX

#endif

#include <float.h>

#ifdef FOUR_BYTE_REAL

typedef float		Real;
#define Real_MAX	FLT_MAX
#define Real_MIN	FLT_MIN
#define Real_DIG	FLT_DIG
#define Real_EPSILON	FLT_EPSILON

#define REAL_MAX	FLT_MAX
#define REAL_MIN	FLT_MIN
#define REAL_DIG	FLT_DIG
#define REAL_EPSILON	FLT_EPSILON

#elif defined EIGHT_BYTE_REAL

typedef double		Real;
#define Real_MAX	DBL_MAX
#define Real_MIN	DBL_MIN
#define Real_DIG	DBL_DIG
#define Real_EPSILON	DBL_EPSILON

#define REAL_MAX	DBL_MAX
#define REAL_MIN	DBL_MIN
#define REAL_DIG	DBL_DIG
#define REAL_EPSILON	DBL_EPSILON

#endif

#include <complex>

//typedef std::complex <Real>		Complex;

/* using C++ */
#ifdef __cplusplus

#else

#endif // __cplusplus

#endif // stk_percept_code_types_hpp
