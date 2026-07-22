// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  ad_functions.hpp
    \brief Provides a generic interface for using automatic differentiation
           with univariate and bivariate functions. The intention of these
           components is make implementing PDE residual terms as simple 
           as writing a function that evaluates c(u,z)

*/

#ifndef PDEOPT_AD_FUNCTIONS_HPP
#define PDEOPT_AD_FUNCTIONS_HPP

#include "template_tools.hpp"

namespace AD {

// Prototype function that returns a scalar quantity
// Takes one or more templated containers with subscript access
template<class Param, template<class> class Array>  
struct ScalarFunction {

  using Real = ElementType<Param>;

  template<class X> 
  static ResultType<Real,X>
  eval( const Param &param,
        const Array<X> &x ) {
    return ResultType<Real,X>(0);
  }

  template<class X, class Y> 
  static ResultType<Real, X, Y> 
  eval( const Param<Real> &param,
        const Array<X> &x, 
        const Array<Y> &y ) {
    return ResultType<Real,X,Y>(0);
  }

}; // ScalarFunction


//struct<class Param, template<class> class Array, class Real, class ScalarX>


// Prototype function that modifies an array in-place
template<class Param, template<class> class Array>
struct ArrayFunction {

  using Real = ElementType<Param>;

  template<class ScalarX>
  static void eval( const Param &param,
                    Array<ResultType<Real,ScalarX>> &result,
                    const Array<ScalarX> &x ) { }

  template<class ScalarX, class ScalarY>
  static void eval( const Param &param,
                    Array<ResultType<Real,ScalarX>> &result,
                    const Array<ScalarX> &x,
                    const Array<ScalarY> &y ) { }

};




} // namespace AD


#endif // PDEOPT_AD_FUNCTIONS_HPP
