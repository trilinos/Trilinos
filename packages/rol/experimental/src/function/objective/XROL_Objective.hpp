
// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
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

#pragma once

#include "XROL_ObjectiveVectors.hpp"

namespace XROL {

template<class X> class Objective;  

// Forward declaration of check methods

template<class X> 
std::vector<std::vector<magnitude_t<X>>>
checkGradient( Objective<X>& obj, const X& x, const dual_t<X>& g, 
                 const X& d, std::ostream& os, Teuchos::ParameterList& parlist );
template<class X>
std::vector<std::vector<magnitude_t<X>>>
checkHessVec( Objective<X>& obj, const X& x, const dual_t<X>& hv, const X& v, 
                   std::ostream& os, Teuchos::ParameterList& parlist );

template<class X>
std::vector<magnitude_t<X>>
checkHessSym( Objective<X>& obj, const X& x, const dual_t<X>& hv, const X& v, 
              const X &w, std::ostream &os, Teuchos::ParameterList& parlist );
 
/** \file  XROL_Objective.hpp
    \brief Defines the interfaces for objective functions 
 */

// Visitor type object that allows access to functionality not found
// in the base class
// template<class X> struct Objective_ExtendedInterface;

// Do-nothing general type
// struct ObjectiveParameters{};


/** \class XROL::Objective
    \brief Basic abstract objective function with default implementations
*/

template<class X> 
class Objective {

//protected:

//  const decltype(auto) getParameters( void ) const {}

public: 

  virtual ~Objective(){}

//  virtual void setParameters( std::unique_ptr<ObjectiveParameters> param ) {
//    ignore(param);
//  }

//  virtual void access( Objective_ExtendedInterface<V,dual_t<V>>& objEI ) {
//    ignore(objEI);
//  }
  virtual void update( const X& x );

  virtual magnitude_t<X> value( const X& x, magnitude_t<X>& tol ) = 0;
   
  virtual void gradient( dual_t<X>& g, const X& x, magnitude_t<X>& tol );

  virtual magnitude_t<X> dirDeriv( const X& x, const X& d, magnitude_t<X>& tol );

  virtual void hessVec( dual_t<X>& hv, const X& v, const X& x, magnitude_t<X>& tol );

  virtual void invHessVec( X& hv, const dual_t<X>& v, const X& x, magnitude_t<X>& tol );

  virtual void precond( X& Pv, const dual_t<X>& v, const X& x, magnitude_t<X>& tol );


  // Friend functions for finite difference checks - give access to ObjectiveVectors
  template<class V>
  friend std::vector<std::vector<magnitude_t<V>>>
  checkGradient( Objective<V>& obj, const V& x, const dual_t<V>& g, 
                 const V& d, std::ostream& os, Teuchos::ParameterList& parlist );


  template<class V>
  friend std::vector<std::vector<magnitude_t<V>>>
  checkHessVec( Objective<V>& obj, const V& x, const dual_t<V>& hv, const V& v,
                std::ostream& os, Teuchos::ParameterList& parlist );

  
  template<class V>
  friend std::vector<magnitude_t<V>>
  checkHessSym( Objective<V>& obj, const V& x, const dual_t<V>& hv, const V& v,
                const V& w, std::ostream& os, Teuchos::ParameterList& parlist );


private:

  std::unique_ptr<ObjectiveVectors<X>> vec_;

};

} // namespace XROL

