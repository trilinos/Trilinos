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

#ifndef ROL_BACKTRACKING_U_H
#define ROL_BACKTRACKING_U_H

#include "ROL_LineSearch_U.hpp"

/** \class ROL::BackTracking_U
    \brief Implements a simple back tracking line search.
*/

namespace ROL {

template<typename Real>
class BackTracking_U : public LineSearch_U<Real> {
private:
  Real rho_;
  Ptr<Vector<Real>> xnew_; 

  using LineSearch_U<Real>::getInitialAlpha;
  using LineSearch_U<Real>::status;

public:

  BackTracking_U(ParameterList &parlist) : LineSearch_U<Real>(parlist) {
    const Real half(0.5);
    rho_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Backtracking Rate",half);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) override {
    LineSearch_U<Real>::initialize(x,g);
    xnew_ = x.clone();
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj ) override {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ls_neval = 0;
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj);
    // Update iterate
    xnew_->set(x); xnew_->axpy(alpha,s);
    // Get objective value at xnew
    Real fold = fval;
    obj.update(*xnew_,UpdateType::Trial);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
    // Perform backtracking
    while ( !status(LINESEARCH_U_BACKTRACKING,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew_,s,obj) ) {
      alpha *= rho_;
      // Update iterate
      xnew_->set(x); xnew_->axpy(alpha,s);
      // Get objective value at xnew
      obj.update(*xnew_,UpdateType::Trial);
      fval = obj.value(*xnew_,tol);
      ls_neval++;
    }
  }
}; // class ROL::BackTracking_U

} // namespace ROL

#endif
