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

#ifndef ROL_PLUSFUNCTION_HPP
#define ROL_PLUSFUNCTION_HPP

#include "ROL_Types.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_Distribution.hpp"
#include "ROL_DistributionFactory.hpp"

namespace ROL {

template<class Real>
class PlusFunction : public PositiveFunction<Real> {
private:
  ROL::Ptr<Distribution<Real> > dist_;
  Real param_;

public: 
  PlusFunction(ROL::Ptr<Distribution<Real> > &dist, Real param = 1.) : dist_(dist) {
    param_ = ((param <= 0) ? 1.e-2 : param);
  }

  PlusFunction(ROL::ParameterList &parlist) {
    Real param(1.e-1), zero(0), one(1);
    ROL::ParameterList pfList;
    if (parlist.isSublist("Plus Function")) {
      param = parlist.sublist("Plus Function").get("Smoothing Parameter",1.);
      pfList = parlist.sublist("Plus Function");
    }
    else {
      param = parlist.get("Smoothing Parameter",1.);
      pfList = parlist;
    }
    param_ = ((param <= zero) ? one : param);
    dist_  = DistributionFactory<Real>(pfList);
  }

  Real evaluate(Real input, int deriv) {
    Real val = 0.0;
    switch(deriv) {
      case 0: val = param_*dist_->integrateCDF(input/param_);   break;
      case 1: val = dist_->evaluateCDF(input/param_);           break;
      case 2: val = dist_->evaluatePDF(input/param_)/param_;    break;
    }
    return val;
  }
 
  void test(Real x) {
    // FIRST DERIVATIVE
    Real vx = evaluate(x,0);
    Real vy = 0.0;
    Real dv = evaluate(x,1);
    Real t = 1.0;
    Real diff = 0.0;
    Real err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK PLUS FUNCTION: p'(x) with x = " 
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "p'(x)"
                            << std::setw(20) << "(p(x+t)-p(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = evaluate(x+t,0);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // SECOND DERIVATIVE
    vx = evaluate(x,1);
    vy = 0.0;
    dv = evaluate(x,2);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK PLUS FUNCTION: p''(x) with x = " 
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "p''(x)"
                            << std::setw(20) << "(p'(x+t)-p'(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = evaluate(x+t,1);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
  } 
};

}

#endif
