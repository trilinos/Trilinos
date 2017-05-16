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
