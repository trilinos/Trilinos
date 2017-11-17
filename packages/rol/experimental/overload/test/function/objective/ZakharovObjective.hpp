
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
#include "XROL_StdVector.hpp"
#include "XROL.hpp"

using dvector = std::vector<double>;

namespace Zakharov {

template<class V> 
class Objective : public XROL::Objective<V> {
public:

  using Real = XROL::element_t<V>;

  Objective( std::unique_ptr<V> k ) : 
    k_(std::move(k)) {}

  Real value( const V& x, Real& tol ) override {
    auto xdotx = XROL::dot(x,x);
    auto kdotx = XROL::dot(*k_,x);
    return xdotx + std::pow(kdotx,2)/4.0 + std::pow(kdotx,4)/16.0;
  }

  void gradient( V& g, const V& x, Real& tol ) override {
    auto kdotx = XROL::dot(x,*k_);
    auto coeff = 0.25*(2.0*kdotx+std::pow(kdotx,3.0));
    XROL::set(g,x);
    XROL::scale(g,2.0);
    XROL::axpy(g,coeff,*k_);
  }

  Real dirDeriv(const V& x, const V& d, Real& tol ) override {
    auto kdotd = XROL::dot(d,*k_);
    auto kdotx = XROL::dot(x,*k_);
    auto xdotd = XROL::dot(x,d);

    auto coeff = 0.25*(2.0*kdotx+std::pow(kdotx,3.0));

    return 2*xdotd + coeff*kdotd;
 
  }

  void hessVec( V& hv, const V& v, const V& x, Real& tol ) override {
    auto kdotx = XROL::dot(x,*k_);
    auto kdotv = XROL::dot(v,*k_);
    auto coeff = 0.25*(2.0+3.0*std::pow(kdotx,2.0))*kdotv;

    XROL::set(hv,v);
    XROL::scale(hv,2.0);
    XROL::axpy(hv,coeff,*k_);
  }

  void invHessVec( V& ihv, const V& v, const V& x, Real& tol ) override {

    auto kdotv = XROL::dot(*k_,v);
    auto kdotx = XROL::dot(*k_,x);
    auto kdotk = XROL::dot(*k_,*k_);
    auto coeff = -kdotv/(2.0*kdotk+16.0/(2.0+3.0*std::pow(kdotx,2.0)));

    XROL::set(ihv,v);
    XROL::scale(ihv,0.5);
    XROL::axpy(ihv,coeff,*k_);
  }

private:
  std::unique_ptr<V> k_;

};

//template<class V>
//auto make_objective( std::unique_ptr<V> k ) {
//  return std::move(std::make_unique<Objective<V>>(std::move(k)));  
//}



} // namespace Zakharov

