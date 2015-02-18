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

#ifndef ROL_DOGLEG_H
#define ROL_DOGLEG_H

/** \class ROL::DogLeg
    \brief Provides interface for dog leg trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

namespace ROL { 

template<class Real>
class DogLeg : public TrustRegion<Real> {
private:

  Teuchos::RCP<CauchyPoint<Real> > cpt_;

  Teuchos::RCP<Vector<Real> > s_;
  Teuchos::RCP<Vector<Real> > Hp_;

  Real pRed_;

public:

  // Constructor
  DogLeg( Teuchos::ParameterList &parlist ) : TrustRegion<Real>(parlist), pRed_(0.0) {
    cpt_ = Teuchos::rcp(new CauchyPoint<Real>(parlist));
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    s_  = s.clone();
    Hp_ = g.clone();
  }

  void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
            const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) { 
    Real tol = std::sqrt(ROL_EPSILON);
    // Compute quasi-Newton step
    pObj.reducedInvHessVec(*s_,grad,x,grad,x,tol);
    s_->scale(-1.0);
    Real sNnorm = s_->norm();
    Real gsN    = s_->dot(grad.dual());
    bool negCurv = false;
    if ( gsN >= 0.0 ) {
      negCurv = true;
    }

    if ( negCurv ) {
      cpt_->run(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
      pRed_ = cpt_->getPredictedReduction();
      iflag = 2;
    }  
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) {        // Use the quasi-Newton step
        s.set(*s_); 
        snorm = sNnorm;
        pRed_ = -0.5*gsN;
        iflag = 0;
      }
      else {                      // quasi-Newton step is outside of trust region
        pObj.reducedHessVec(*Hp_,grad,x,grad,x,tol);
        Real alpha  = 0.0;
        Real beta   = 0.0;
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = grad.dot(*Hp_);
        Real gamma  = gnorm2/gBg;
        if ( gamma*gnorm >= del || gBg <= 0.0 ) {
            alpha = 0.0;
            beta  = del/gnorm;
            s.set(grad.dual()); 
            s.scale(-beta); 
            snorm = del;
            iflag = 2;
        }
        else {
          Real a = sNnorm*sNnorm + 2.0*gamma*gsN + gamma*gamma*gnorm2;
          Real b = -gamma*gsN - gamma*gamma*gnorm2;
          Real c = gamma*gamma*gnorm2 - del*del;
          alpha  = (-b + sqrt(b*b - a*c))/a;
          beta   = gamma*(1.0-alpha);
          s.set(grad.dual());
          s.scale(-beta);
          s.axpy(alpha,*s_);
          snorm = del;
          iflag = 1;
        }
        pRed_ = (alpha*(0.5*alpha-1)*gsN - 0.5*beta*beta*gBg + beta*(1-alpha)*gnorm2);
      }
    }
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }
};

}

#endif
