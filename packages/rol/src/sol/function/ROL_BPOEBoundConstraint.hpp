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

#ifndef ROL_BPOE_BOUND_CONSTRAINT_H
#define ROL_BPOE_BOUND_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_CVaRVector.hpp"

namespace ROL {

template <class Real>
class BPOEBoundConstraint : public BoundConstraint<Real> {
private:
  Teuchos::RCP<BoundConstraint<Real> > bc_;

public:
  BPOEBoundConstraint(void) : bc_(Teuchos::null) {}

  BPOEBoundConstraint(Teuchos::RCP<BoundConstraint<Real> > &bc) : BoundConstraint<Real>(), bc_(bc) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<const Vector<Real> > xv
        = (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->update(*xv,flag,iter);
    }
  }

  void project( Vector<Real> &x ) {
    Real xvar = Teuchos::dyn_cast<CVaRVector<Real> >(x).getVaR();
    xvar = ((xvar > 0.0) ? xvar : 0.0);
    (Teuchos::dyn_cast<CVaRVector<Real> >(x)).setVaR(xvar);
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > xvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(x)).getVector());
      bc_->project(*xvec);
      (Teuchos::dyn_cast<CVaRVector<Real> >(x)).setVector(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneUpperActive(*vvec,*xvec,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneUpperActive(*vvec,*gvec,*xvec,eps);
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    Real xvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x)).getVaR();
    if ( xvar <= eps ) {
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).setVaR(0.0);
    } 
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec
        = (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneLowerActive(*vvec,*xvec,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    Real gvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g)).getVaR();
    Real xvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x)).getVaR();
    if ( xvar <= eps && gvar > 0.0 ) {
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).setVaR(0.0);
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec
        = (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec
        = (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneLowerActive(*vvec,*gvec,*xvec,eps);
    }
  } 

  void setVectorToUpperBound( Vector<Real> &u ) {
    Real uvar = 0.1*ROL_OVERFLOW;
    (Teuchos::dyn_cast<CVaRVector<Real> >(u)).setVaR(uvar);
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > uvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(u)).getVector());
      bc_->setVectorToUpperBound(*uvec);
    }
  }

  void setVectorToLowerBound( Vector<Real> &l ) {
    (Teuchos::dyn_cast<CVaRVector<Real> >(l)).setVaR(0.0);
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > lvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(l)).getVector());
      bc_->setVectorToLowerBound(*lvec);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    Real xvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x)).getVaR();
    if ( xvar <= eps ) {
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).setVaR(0.0);
    } 
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneActive(*vvec,*xvec,eps);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    Real gvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g)).getVaR();
    Real xvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x)).getVaR();
    if ( xvar <= eps && gvar > 0.0 ) {
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).setVaR(0.0);
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneActive(*vvec,*gvec,*xvec,eps);
    }
  }

  bool isFeasible( const Vector<Real> &v ) { 
    bool flag = false;
    Real vvar = Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(v)).getVaR();
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<const Vector<Real> > vvec
        = (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      if ( bc_->isActivated() ) {
        flag = ((bc_->isFeasible(*vvec)) && (vvar >= 0.0));
      }
    }
    else {
      if ( bc_->isActivated() ) {
        flag = (vvar >= 0.0);
      }
    }
    return flag;
  }

}; // class BPOEBoundConstraint

} // namespace ROL

#endif
