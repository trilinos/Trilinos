#ifndef INC_llnlmath_h
#define INC_llnlmath_h

namespace CVODE {
/******************************************************************
 *                                                                *
 * File          : llnlmath.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 4 May 1998                                     *
 *----------------------------------------------------------------*
 * This is the header file for a C math library. The routines     *
 * listed here work with the type float as defined in llnltyps.h.  *
 * To do single precision floating point arithmetic, set the type *
 * float to be float. To do double precision arithmetic, set the   *
 * type float to be double. The default implementations for        *
 * RPowerR and RSqrt call standard math library functions which   *
 * do double precision arithmetic. If this is unacceptable when   *
 * float is float, then the user should re-implement these two     *
 * routines by calling single precision routines available on     *
 * his/her machine.                                               *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "llnltyps.h"


/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 *                                                                *
 * MAX(A, B) returns the maximum of A and B.                      *
 *                                                                *
 * ABS(A) returns the absolute value of A.                        *
 *                                                                *
 * SQR(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/
#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS(A)    ((A < 0) ? -(A) : (A))
#endif

#ifndef SQR
#define SQR(A)    ((A) * (A))
#endif

/******************************************************************
 *                                                                *
 * Function : UnitRoundoff                                        *
 * Usage    : float uround;                                        *
 *            uround = UnitRoundoff();                            *
 *----------------------------------------------------------------*
 * UnitRoundoff returns the unit roundoff u for float floating     *
 * point arithmetic, where u is defined to be the smallest        *
 * positive float such that 1.0 + u != 1.0.                        *
 *                                                                *
 ******************************************************************/

float UnitRoundoff(void);


/******************************************************************
 *                                                                *
 * Function : RPowerI                                             *
 * Usage    : int exponent;                                       *
 *            float base, ans;                                     *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerI returns the value base^exponent, where base is a float  *
 * and exponent is an int.                                        *
 *                                                                *
 ******************************************************************/

float RPowerI(float base, int exponent);


/******************************************************************
 *                                                                *
 * Function : RPowerR                                             *
 * Usage    : float base, exponent, ans;                           *
 *            ans = RPowerR(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerR returns the value base^exponent, where both base and   *
 * exponent are reals. If base < 0.0, then RPowerR returns 0.0.   *
 *                                                                *
 ******************************************************************/

float RPowerR(float base, float exponent);


/******************************************************************
 *                                                                *
 * Function : RSqrt                                               *
 * Usage    : float sqrt_x;                                        *
 *            sqrt_x = RSqrt(x);                                  *
 *----------------------------------------------------------------*
 * RSqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
 * returns 0.0.                                                   *
 *                                                                *
 ******************************************************************/

float RSqrt(float x);

} // namespace CVODE 

#endif // INC_llnlmath_h

#ifdef __cplusplus
}
#endif

