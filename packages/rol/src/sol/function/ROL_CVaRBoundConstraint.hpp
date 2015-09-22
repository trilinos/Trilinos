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

#ifndef ROL_CVAR_BOUND_CONSTRAINT_H
#define ROL_CVAR_BOUND_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_CVaRVector.hpp"

namespace ROL {

template <class Real>
class CVaRBoundConstraint : public BoundConstraint<Real> {
private:
  Teuchos::RCP<BoundConstraint<Real> > bc_;

public:
  CVaRBoundConstraint(Teuchos::RCP<BoundConstraint<Real> > &bc) : BoundConstraint<Real>(), bc_(bc) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->update(*xv,flag,iter);
  }

  void project( Vector<Real> &x ) {
    Teuchos::RCP<Vector<Real> > xv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(x)).getVector());
    bc_->project(*xv);
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneUpperActive(*vv,*xv,eps);
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > gv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneUpperActive(*vv,*gv,*xv,eps);
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneLowerActive(*vv,*xv,eps);
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > gv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneLowerActive(*vv,*gv,*xv,eps);
} 

  void setVectorToUpperBound( Vector<Real> &u ) {
    Teuchos::RCP<Vector<Real> > uv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(u)).getVector());
    bc_->setVectorToUpperBound(*uv);
    Real uvar = 0.1*ROL_OVERFLOW;
    (Teuchos::dyn_cast<CVaRVector<Real> >(u)).setVaR(uvar);
  }

  void setVectorToLowerBound( Vector<Real> &l ) {
    Teuchos::RCP<Vector<Real> > lv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(l)).getVector());
    bc_->setVectorToLowerBound(*lv);
    Real lvar = -0.1*ROL_OVERFLOW;
    (Teuchos::dyn_cast<CVaRVector<Real> >(l)).setVaR(lvar);
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneActive(*vv,*xv,eps);
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    Teuchos::RCP<Vector<Real> > vv = Teuchos::rcp_const_cast<Vector<Real> >(
      (Teuchos::dyn_cast<CVaRVector<Real> >(v)).getVector());
    Teuchos::RCP<const Vector<Real> > gv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
    Teuchos::RCP<const Vector<Real> > xv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
    bc_->pruneActive(*vv,*gv,*xv,eps);
  }

  bool isFeasible( const Vector<Real> &v ) { 
    Teuchos::RCP<const Vector<Real> > vv =
        (Teuchos::dyn_cast<CVaRVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
    if ( bc_->isActivated() ) {
      return bc_->isFeasible(*vv);
    }
    else {
      return false;
    }
  }

}; // class CVaRBoundConstraint

} // namespace ROL

#endif
