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

#ifndef ROL_PATHBASEDTARGETLEVEL_H
#define ROL_PATHBASEDTARGETLEVEL_H

/** \class ROL::PathBasedTargetLevel
    \brief Provides an implementation of path-based target leve line search.
*/

#include "ROL_LineSearch.hpp"

namespace ROL { 

template<class Real>
class PathBasedTargetLevel : public LineSearch<Real> {
private:
  Teuchos::RCP<Vector<Real> > xnew_; 

  Real min_value_;
  Real rec_value_;
  Real target_;
  Real delta_;
  Real sigma_;
  Real bound_;

public:

  virtual ~PathBasedTargetLevel() {}

  // Constructor
  PathBasedTargetLevel( Teuchos::ParameterList &parlist ) 
    : LineSearch<Real>(parlist), min_value_(ROL::ROL_OVERFLOW), rec_value_(ROL::ROL_OVERFLOW), 
      target_(0.0), sigma_(0.0) {
    delta_ = parlist.get("Path-Based Target Level: Target Relaxation Parameter",0.1);
    bound_ = parlist.get("Path-Based Target Level: Upper Bound on Path Length",1.0);
  }

  void initialize(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g, 
                  Objective<Real> &obj, BoundConstraint<Real> &con) {
    LineSearch<Real>::initialize(x,g,obj,con);
    xnew_ = x.clone();
  }

  // Run Iteration scaled line search
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);
    ls_neval = 0;
    ls_ngrad = 0;
    // Update target objective value
    if ( fval < min_value_ ) {
      min_value_ = fval;
    }
    target_ = rec_value_ - 0.5*delta_;
    if ( fval < target_ ) {
      rec_value_ = min_value_; 
      sigma_ = 0.0;
    }
    else {
      if ( sigma_ > bound_ ) {
        rec_value_ = min_value_;
        sigma_ = 0.0;
        delta_ *= 0.5;
      }
    }
    target_ = rec_value_ - delta_;
    // Get line-search parameter
    alpha = (fval - target_)/std::abs(gs);
    // Update iterate
    LineSearch<Real>::updateIterate(*xnew_,x,s,alpha,con);
    // Compute objective function value
    obj.update(*xnew_);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
    // Update sigma 
    sigma_ += alpha*std::sqrt(std::abs(gs));
  }
};

}

#endif
