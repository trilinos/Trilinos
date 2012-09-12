// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_CMATH_HPP
#define SACADO_CMATH_HPP

#include <cmath>        // for most math functions
#include <algorithm>	// for std::min and std::max
#include "Sacado_ConfigDefs.h"

// Define some math functions that aren't usually in cmath
#if !( defined(HAS_C99_TR1_CMATH) || defined(_GLIBCXX_USE_C99_MATH_TR1) || defined(USER_DISABLE_SACADO_TR1_CMATH) )
namespace std {
  inline float acosh(float x) { 
    return std::log(x + std::sqrt(x*x - float(1.0))); }
  inline float asinh(float x) { 
    return std::log(x + std::sqrt(x*x + float(1.0))); }
  inline float atanh(float x) { 
    return float(0.5)*std::log((float(1.0)+x)/(float(1.0)-x)); }

  inline double acosh(double x) { 
    return std::log(x + std::sqrt(x*x - double(1.0))); }
  inline double asinh(double x) { 
    return std::log(x + std::sqrt(x*x + double(1.0))); }
  inline double atanh(double x) { 
    return double(0.5)*std::log((double(1.0)+x)/(double(1.0)-x)); }
}
#endif // HAS_C99_TR1_CMATH

#endif // SACADO_CMATH_HPP
