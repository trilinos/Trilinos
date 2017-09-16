
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

#include "XROL.hpp"


namespace XROL {


template<class V> struct ObjectiveVisitor;

struct ObjectiveParameters;

template<class V> 
struct Objective {

private:

  std::unique_ptr<ObjectiveParameters> param_;

protected:

  const decltype(auto) getParameters( void ) const {
    return std::move(param_);
  }

public: 

  virtual void setParameters( const std::unique_ptr<ObjectiveParameters> &param ) {
    param_ = std::move(param);
  }

  Objective( const std::unique_ptr<ObjectiveParmeters> &param ) 
    : param_(std::move(param)) {}

  virtual ~Objective() {}

  virtual void accept( ObjectiveVisitor<V>& visitor ) = 0;

  virtual void update( /* args */ ) {}
 
  virtual auto value( const V& x ) = 0;
   
  virtual void gradient( V& g, const V& x ) { ignore(g,x); }

  virtual auto dirDeriv( const V& x, const V& d ) { ignore(x,d); }

  virtual void hessVec( V& hv, const V& v, const V& x ) { ignore(hv,v,x); } 

  virtual void invHessVec( V& hv, const V& v, const V& x ) { ignore(hv,v,x); }

  virtual void precond( V& Pv, const V& v, const V& x ) { ignore(Pv,v,x); }

  virtual auto checkGradient( const V& x, const V& d, std::ostream &os,
                              Teuchos::ParameterList &parlist ) { 
    ignore(x,d,os,parlist);  
  }

  virtual auto checkGradient( const V& x, const V &g, const V& d, std::ostream &os,
                              Teuchos::ParameterList &parlist ) {    
    ignore(x,g,d,os,parlist)
  }

  virtual auto checkHessVec( const V& x, const V&v, std::ostream &os,
                             Teuchos::ParmeterList &parlist ) {
    ignore(x,v,os,parlist);
  }


  virtual auto checkHessVec( const V& x, const V& hv, const V&v, std::ostream &os,
                             Teuchos::ParmeterList &parlist ) {
    ignore(x,hv,v,os,parlist);
  }
  
  virtual auto checkHessSym( const V& x, const V& v, const V& w, std::ostream & os,
                             Teuchos::parameterList& parlist ) {
    ignore(x,v,w,os,parlist);
  }             

  virtual auto checkHessSym( const V& x, const V& hv, const V& v, const V& w, 
                             std::ostream & os, Teuchos::parameterList& parlist ) {
    ignore(x,hv,v,w,os,parlist);
  }             

};


} // namespace XROL

