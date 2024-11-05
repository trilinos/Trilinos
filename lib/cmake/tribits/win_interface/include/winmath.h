// @HEADER
// *****************************************************************************
//            TriBITS: Tribal Build, Integrate, and Test System
//
// Copyright 2013-2016 NTESS and the TriBITS contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef WINMATH_H
#define WINMATH_H

/**********************************************************************

  acosh.c -

  $Author$
  $Date$
  created at: Fri Apr 12 00:34:17 JST 2002

  public domain rewrite of acosh(3), asinh(3) and atanh(3)

**********************************************************************/

#include <errno.h>
#include <float.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* DBL_MANT_DIG must be less than 4 times of bits of int */
#ifndef DBL_MANT_DIG
#define DBL_MANT_DIG 53         /* in this case, at least 12 digit precision */
#endif
#define BIG_CRITERIA_BIT (1<<DBL_MANT_DIG/2)
#if BIG_CRITERIA_BIT > 0
#define BIG_CRITERIA (1.0*BIG_CRITERIA_BIT)
#else
#define BIG_CRITERIA (1.0*(1<<DBL_MANT_DIG/4)*(1<<(DBL_MANT_DIG/2+1-DBL_MANT_DIG/4)))
#endif
#define SMALL_CRITERIA_BIT (1<<(DBL_MANT_DIG/3))
#if SMALL_CRITERIA_BIT > 0
#define SMALL_CRITERIA (1.0/SMALL_CRITERIA_BIT)
#else
#define SMALL_CRITERIA (1.0*(1<<DBL_MANT_DIG/4)*(1<<(DBL_MANT_DIG/3+1-DBL_MANT_DIG/4)))
#endif

inline double acosh(double x)
{
    if (x < 1)
        x = -1;                 /* NaN */
    else if (x == 1)
        return 0;
    else if (x > BIG_CRITERIA)
        x += x;
    else
        x += sqrt((x + 1) * (x - 1));
    return log(x);
}

inline double asinh(double x)
{
    int neg = x < 0;
    double z = fabs(x);

    if (z < SMALL_CRITERIA) return x;
    if (z < (1.0/(1<<DBL_MANT_DIG/5))) {
        double x2 = z * z;
        z *= 1 + x2 * (-1.0/6.0 + x2 * 3.0/40.0);
    }
    else if (z > BIG_CRITERIA) {
        z = log(z + z);
    }
    else {
        z = log(z + sqrt(z * z + 1));
    }
    if (neg) z = -z;
    return z;
}

inline double atanh(double x)
{
    int neg = x < 0;
    double z = fabs(x);

    if (z < SMALL_CRITERIA) return x;
    z = log(z > 1 ? -1 : (1 + z) / (1 - z)) / 2;
    if (neg) z = -z;
    return z;
}

inline double round(double val)
{
    return floor(val + 0.5);
}

inline void srand48(double seed)
{
  srand(seed);
}

inline double drand48()
{
  return (double(rand()) / RAND_MAX);
}

inline float tgammaf(float z){
  static unsigned int c = 9;
  static float lanczos_coefficients[] = {1.000000000000000174663f,
                                         5716.400188274341379136f,
                                         -14815.30426768413909044f,
                                         14291.49277657478554025f,
                                         -6348.160217641458813289f,
                                         1301.608286058321874105f,
                                         -108.1767053514369634679f,
                                         2.605696505611755827729f,
                                         -0.7423452510201416151527e-2f,
                                         0.5384136432509564062961e-7f,
                                         -0.4023533141268236372067e-8f
                                        };
  float return_val = 0.0f;

  z -= 1.0f;
  float temp = z + c + 0.5f;
  float A = lanczos_coefficients[0];
  int i = 0;
  for(i = 1; i < c+2; ++i){
    A += lanczos_coefficients[i]/(z+i);
  }
  return_val = sqrt(2 * M_PI)* pow(temp, z + 0.5f) * exp(-temp) * A;

  return return_val;
}

inline double tgamma(double z){
  static unsigned int c = 9;
  static double lanczos_coefficients[] = {1.000000000000000174663,
                                         5716.400188274341379136,
                                         -14815.30426768413909044,
                                         14291.49277657478554025,
                                         -6348.160217641458813289,
                                         1301.608286058321874105,
                                         -108.1767053514369634679,
                                         2.605696505611755827729,
                                         -0.7423452510201416151527e-2,
                                         0.5384136432509564062961e-7,
                                         -0.4023533141268236372067e-8
                                        };
  double return_val = 0.0;

  z -= 1.0;
  double temp = z + c + 0.5;
  double A = lanczos_coefficients[0];
  int i = 0;
  for(i = 1; i < c+2; ++i){
    A += lanczos_coefficients[i]/(z+i);
  }
  return_val = sqrt(2 * M_PI)* pow(temp, z + 0.5) * exp(-temp) * A;

  return return_val;
}

inline long double tgammal(long double z){
  static unsigned int c = 9;
  static long double lanczos_coefficients[] = {1.000000000000000174663,
                                         5716.400188274341379136,
                                         -14815.30426768413909044,
                                         14291.49277657478554025,
                                         -6348.160217641458813289,
                                         1301.608286058321874105,
                                         -108.1767053514369634679,
                                         2.605696505611755827729,
                                         -0.7423452510201416151527e-2,
                                         0.5384136432509564062961e-7,
                                         -0.4023533141268236372067e-8
                                        };
  long double return_val = 0.0;

  z -= 1.0;
  long double temp = z + c + 0.5;
  long double A = lanczos_coefficients[0];
  int i = 0;
  for(i = 1; i < c+2; ++i){
    A += lanczos_coefficients[i]/(z+i);
  }
  return_val = sqrt(2 * M_PI)* pow(temp, z + 0.5) * exp(-temp) * A;

  return return_val;
}

// This function was adapted from a public domain implementation of erf
// which is available at http://www.johndcook.com/cpp_erf.html. The only
// changes made were to change the type from double to float.
inline float erff(float x)
{
    // constants
    float a1 =  0.254829592f;
    float a2 = -0.284496736f;
    float a3 =  1.421413741f;
    float a4 = -1.453152027f;
    float a5 =  1.061405429f;
    float p  =  0.3275911f;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    float t = 1.0f/(1.0f + p*x);
    float y = 1.0f - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

#endif
