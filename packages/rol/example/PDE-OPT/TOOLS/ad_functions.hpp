// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  ad_functions.hpp
    \brief Provides a generic interface for using automatic differentiation
           with univariate and bivariate functions. The intention of these
           components is make implementing PDE residual terms as simple 
           as writing a function that evaluates c(u,z)

*/

#ifndef PDEOPT_AD_FUNCTIONS_HPP
#define PDEOPT_AD_FUNCTIONS_HPP

#include "Sacado.hpp"
#include <vector>

namespace AD {

namespace Univariate {

/* The univariate Evaluator operates on objects of which follow the pattern
 *
 * template<class T,class U>
 * class Function {
 * public:
 *   static U eval( const U&, const std::vector<T>& ); 
 * };
 *
 * where T is going to be a Real type (e.g. double) and the std::vector stores 
 * Real parameters used in evaluating an expression. U will be a Real or AD
 * type depending on what level Evaluator is calling it.
 *
 *
 * The template parameters are
 *
 * N       - order of derivative to evaluate
 * F       - Function to be evaluated
 * Real    - unchanging type of the parameters appearing the the function
 * ScalarT - changing type depending on N
 *
 */
template<int N, template<class,class> class F, class Real, class ScalarT=Real>
class Evaluator {
public:
  using Fad   = Sacado::Fad::SFad<ScalarT,1>;
  using PARAM = std::vector<Real>;

  static ScalarT eval( const ScalarT &x, const PARAM& param ) {
      Fad x_fad(1,0,x);
      Fad f_fad = Evaluator<N-1,F,Real,Fad>::eval(x_fad, param);
      return f_fad.dx(0);
  } 
}; // class Evaluator

// Specialize the evaluator for stopping case
template<template<class,class> class F, class Real, class ScalarT>
class Evaluator<0,F,Real,ScalarT> {
public:
  using PARAM = std::vector<Real>;

  static ScalarT eval( const ScalarT &x, const PARAM &param ) {
    return F<Real,ScalarT>::eval(x,param);
  }
}; // class Evaluator


// An example function implementation
//
// Sample usage
// std::vector<RealT> coeff({1.0,2.0,1.0});         // f(x) = x^2 + 2x + 1
// Evaluator<0,PowerSeries,RealT>::eval(x,coeff);   // evaluates f(x)
// Evaluator<1,PowerSeries,RealT>::eval(x,coeff);   // evaluates f'(x)

template<class Real,class ScalarT>
class PowerSeries {
public:
  using COEFF = std::vector<Real>;
  static ScalarT eval( const ScalarT &x, const COEFF &coeff ) {
    ScalarT sum(0);
    for( auto c: coeff ) {
      sum *=x; sum +=c;
    }
    return sum;
  }
};

} // namespace Univariate

namespace Bivariate {

/* The bivariate evaluator is similar to the above, except that it operates on 
 * objects whose eval function takes two arguments before the parameters
 *
 * The template parameters are
 *
 * Nx      - order of derivative w.r.t the first argument to evaluate
 * Ny      - order of derivative w.r.t the second argument to evaluate 
 * F       - Function to be evaluated
 * Real    - unchanging type of the parameters appearing the the function
 * ScalarT - changing type depending on N
 *
 */

// General case
template<int Nx, int Ny, template<class,class> class F, class Real, class ScalarT=Real> 
class Evaluator {
public:
  using Fad    = Sacado::Fad::SFad<ScalarT,2>;
  using PARAM  = std::vector<Real>;
  
  static ScalarT eval( const ScalarT &x, const ScalarT &y, const PARAM &param ) {
    Fad x_fad(2,0,x);
    Fad y_fad(2,1,y);
    Fad f_fad = Evaluator<Nx,Ny-1,F,Real,Fad>::eval(x_fad,y_fad,param);
    return f_fad.dx(1); 
  }
}; // Evaluator


// Edge case 1
template<int Nx, template<class,class> class F, class Real, class ScalarT>
class Evaluator<Nx, 0, F, Real, ScalarT> {
public:
  using Fad   = Sacado::Fad::SFad<ScalarT,2>;
  using PARAM = std::vector<Real>;

  static ScalarT eval( const ScalarT &x, const ScalarT &y, const PARAM &param ) {
    Fad x_fad(2,0,x);
    Fad y_fad(2,1,y);
    Fad f_fad = Evaluator<Nx-1, 0, F,Real,Fad>::eval(x_fad,y_fad,coeff);
    return f_fad.dx(0);
  }
}; 

// Edge case 2
template<int Ny, template<class,class> class F, class Real, class ScalarT>
class Evaluator<0, Ny, F, Real, ScalarT> {
public:
  using Fad   = Sacado::Fad::SFad<ScalarT,2>;
  using PARAM = std::vector<Real>;

  static ScalarT eval( const ScalarT &x, const ScalarT &y, const PARAM& param ) {
    Fad x_fad(2,0,x);
    Fad y_fad(2,1,y);
    Fad f_fad = Evaluator<0, Ny-1, F,Real,Fad>::eval(x_fad,y_fad,param);
    return f_fad.dx(1);
  }
};

// Vertex case
template<template<class,class> class F, class Real, class ScalarT>
class Evaluator<0, 0, F, Real, ScalarT> {
public:
  using PARAM = std::vector<Real>;

  static ScalarT eval( const ScalarT &x, const ScalarT &y, const PARAM& param ) {
    return F<Real,ScalarT>::eval(x,y,param);  
  }
};

// An example function implementation
// 
// Sample usage
// std::vector<RealT> coeff({1.0,0.0,2.0,0.0,1.0,-3.0}); // f(x,y) = x^2 + 2y^2 + y - 3
// Evaluator<1,1,Quadratic,RealT>::eval(x, y, coeff);    // evaluates f_{xy}(x,y)
//
template<class Real, class ScalarT>
class Quadratic {
public:
  static ScalarT eval( const ScalarT &x, const ScalarT &y, const std::vector<Real> &coeff ) {
    return coeff[0]*x*x + coeff[1]*x*y + coeff[2]*y*y + coeff[3]*x + coeff[4]*y + coeff[5];
  }
};

}

} // namespace AD


#endif // PDEOPT_AD_FUNCTIONS_HPP
