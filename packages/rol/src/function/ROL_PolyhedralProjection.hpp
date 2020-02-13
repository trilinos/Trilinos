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


#ifndef ROL_POLYHEDRALPROJECTION_H
#define ROL_POLYHEDRALPROJECTION_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_KrylovFactory.hpp"

namespace ROL {

template<typename Real>
class PolyhedralProjection {
private:
  const Ptr<BoundConstraint<Real>> bnd_;
  const Ptr<Constraint<Real>>      con_;
  const Ptr<Vector<Real>>          mul_;

  int dim_;
  Ptr<Krylov<Real>> krylov_;
  Ptr<Vector<Real>> xnew_, xdual_, xprim_, lnew_, dlam_, res_;
  Real b_, mul1_, dlam1_;

public:
  virtual ~PolyhedralProjection() {}

  PolyhedralProjection(const Ptr<const Vector<Real>>    &x,
                       const Ptr<BoundConstraint<Real>> &bnd,
                       const Ptr<Constraint<Real>>      &con = nullPtr,
                       const Ptr<Vector<Real>>          &mul = nullPtr)
    : bnd_(bnd), con_(con), mul_(mul) {
    xnew_  = x->clone();
    xprim_ = x->clone();
    xdual_ = x->dual().clone();

    dim_ = 0;
    if (mul_ != nullPtr) {
      dim_ = mul->dimension();
      res_ = mul->dual().clone();
    }
    if (dim_ == 1) {
      mul1_  = static_cast<Real>(0);
      dlam1_ = static_cast<Real>(2);
      // con.value(x) = xprim_->dot(x) + b_
      Real tol(std::sqrt(ROL_EPSILON<Real>()));
      xprim_->zero();
      con_->value(*res_,*xprim_,tol);
      b_ = res_->dot(*res_->basis(0));
      mul_->setScalar(static_cast<Real>(1));
      con_->applyAdjointJacobian(*xdual_,*mul_,*x,tol);
      xprim_->set(xdual_->dual());
    }
    else if (dim_ > 1) {
      lnew_ = mul->clone();
      dlam_ = mul->clone();

      ParameterList list;
      list.sublist("General").sublist("Krylov").set("Type",               "CG");
      list.sublist("General").sublist("Krylov").set("Absolute Tolerance", 1e-6);
      list.sublist("General").sublist("Krylov").set("Relative Tolerance", 1e-4);
      list.sublist("General").sublist("Krylov").set("Iteration Limit",    dim_);
      list.sublist("General").set("Inexact Hessian-Times-A-Vector",      false);
      krylov_ = KrylovFactory<Real>(list);
    }
  }

  void project(Vector<Real> &x) {
    if (con_ == nullPtr) {
      bnd_->project(x);
    }
    else {
      if (dim_==1) {
        mul1_  = static_cast<Real>(0);
        dlam1_ = static_cast<Real>(2);
        //dlam1_ = static_cast<Real>(1)+std::abs(mul1_);
        project_1d(x,mul1_,dlam1_);
      }
      else {
        //std::cout << std::scientific << std::setprecision(16);
        //x.print(std::cout);
        project_nd(x,*mul_,*dlam_);
        //x.print(std::cout);
      }
    }
  }

  const Ptr<Constraint<Real>> getLinearConstraint(void) const {
    return con_;
  }

  const Ptr<Vector<Real>> getMultiplier(void) const {
    return mul_;
  }

private:

  Real residual_1d(const Vector<Real> &x) const {
    return xprim_->dot(x) + b_;
  }

  void update_primal_1d(Vector<Real> &y, const Vector<Real> &x, const Real lam) const {
    y.set(x);
    y.axpy(lam,*xprim_);
    bnd_->project(y);
  }

  void project_1d(Vector<Real> &x, Real &lam, Real &dlam) const {
    const Real atol(1e-4*std::sqrt(ROL_EPSILON<Real>())), rtol(1e-2);
    const Real zero(0), one(1), two(2), c1(0.1), c2(0.75), c3(0.25);
    Real lam1(0), lam2(0), lam3(0), r1(0), r2(0), s(0);
    update_primal_1d(*xnew_,x,lam);
    Real r = residual_1d(*xnew_);
    const Real ctol = std::min(atol,rtol*std::abs(r));
    int cnt = 0, maxit = 1000;
    // Bracketing phase
    if ( r < zero ) {
      lam1 = lam;
      lam += dlam;
      r1   = r;
      update_primal_1d(*xnew_,x,lam);
      r    = residual_1d(*xnew_);
      while ( r < zero && cnt < maxit ) {
        lam1  = lam;
        r1    = r;
        s     = std::max(r1/r-one,c1);
        dlam += dlam/s;
        lam  += dlam;
        update_primal_1d(*xnew_,x,lam);
        r     = residual_1d(*xnew_);
	cnt++;
      }
      lam2 = lam;
      r2   = r;
    }
    else {
      lam2 = lam;
      r2   = r;
      lam -= dlam;
      update_primal_1d(*xnew_,x,lam);
      r    = residual_1d(*xnew_);
      while ( r > zero && cnt < maxit ) {
        lam2  = lam;
        r2    = r;
        s     = std::max(r2/r-one,c1);
        dlam += dlam/s;
        lam  -= dlam;
        update_primal_1d(*xnew_,x,lam);
        r     = residual_1d(*xnew_);
      }
      lam1 = lam;
      r1   = r;
    }

    // Secant phase
    s     = one - r1/r2;
    dlam /= s;
    lam   = lam2-dlam;
    update_primal_1d(*xnew_,x,lam);
    r     = residual_1d(*xnew_);
    cnt   = 0;
    while ( std::abs(r) > ctol && std::abs(dlam) > ctol && cnt < maxit ) {
      if ( r > zero ) {
        if ( s <= two ) {
          lam2 = lam;
          r2   = r;
          s    = one-r1/r2;
          dlam = (lam2-lam1)/s;
          lam  = lam2-dlam;
        }
        else {
          s    = std::max(r2/r-one,c1);
          dlam = (lam2-lam)/s;
          lam3 = std::max(lam-dlam,c2*lam1+c3*lam);
          lam2 = lam;
          r2   = r;
          lam  = lam3;
          s    = (lam2-lam1)/(lam2-lam);
        }
      }
      else {
        if ( s >= two ) {
          lam1 = lam;
          r1   = r;
          s    = one-r1/r2;
          dlam = (lam2-lam1)/s;
          lam  = lam2-dlam;
        }
        else {
          s    = std::max(r1/r-one,c1);
          dlam = (lam-lam1)/s;
          lam3 = std::min(lam+dlam,c2*lam2+c3*lam);
          lam1 = lam;
          r1   = r;
          lam  = lam3;
          s    = (lam2-lam1)/(lam2-lam);
        }
	cnt++;
      }
      update_primal_1d(*xnew_,x,lam);
      r = residual_1d(*xnew_);
    }
    if (std::abs(r) > ctol) {
      //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
      std::cout << ">>> ROL::PolyhedralProjection::project : Projection failed!  rnorm = " << std::abs(r) << "  rtol = " << ctol << std::endl;
    }
    // Return projection
    x.set(*xnew_);
  }

  // r = A y + b
  // y = P(x + A^T lam)
  Real residual_nd(Vector<Real> &r, const Vector<Real> &y) const {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    con_->value(r,y,tol);
    return r.norm();
  }

  // Jv = inv(A DP(y) A^T) v
  // y  = P(x + A^T lam)
  class Jacobian_nd : public LinearOperator<Real> {
    private:
      const Ptr<Constraint<Real>>      con_;
      const Ptr<BoundConstraint<Real>> bnd_;
      const Ptr<const Vector<Real>>    y_;
      const Ptr<Vector<Real>> xdual_, xprim_;
      const Real alpha_;
    public:
      Jacobian_nd(const Ptr<Constraint<Real>>      &con,
                  const Ptr<BoundConstraint<Real>> &bnd,
                  const Ptr<const Vector<Real>>    &y,
                  const Ptr<Vector<Real>>          &xdual,
                  const Ptr<Vector<Real>>          &xprim,
                  const Real                        alpha = 1e-4)
        : con_(con), bnd_(bnd), y_(y), xdual_(xdual), xprim_(xprim), alpha_(alpha) {}
      void apply(Vector<Real> &Jx, const Vector<Real> &x, Real &tol) const {
        con_->applyAdjointJacobian(*xdual_,x,*y_,tol);
        xprim_->set(xdual_->dual());
        bnd_->pruneActive(*xprim_,*y_);
        con_->applyJacobian(Jx,*xprim_,*y_,tol);
        // This is a hack to make the Jacobian invertible
        Jx.axpy(alpha_,x);
      }
  };
  class Precond_nd : public LinearOperator<Real> {
  public:
    void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      Hv.set(v.dual());
    }
    void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      Hv.set(v.dual());
    }
  };
  void solve_newton_system_nd(Vector<Real> &s, const Vector<Real> &r, const Vector<Real> &y,
                              const Real mu, const Real rho) const {
    int iter(0), flag(0);
    Ptr<Precond_nd>  M = makePtr<Precond_nd>();
    Ptr<Jacobian_nd> J = makePtr<Jacobian_nd>(con_,bnd_,makePtrFromRef(y),xdual_,xprim_, mu);
    krylov_->run(s,*J,r,*M,iter,flag);
  }

  void update_primal_nd(Vector<Real> &y, const Vector<Real> &x, const Vector<Real> &lam) const {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    y.set(x);
    con_->applyAdjointJacobian(*xdual_,lam,x,tol);
    y.plus(xdual_->dual());
    bnd_->project(y);
  }

  // Inexact Newton method for monotone dual optimality system
  // motivated by Solodov and Svaiter 1999
  void project_nd(Vector<Real> &x, Vector<Real> &lam, Vector<Real> &dlam) const {
    const Real atol(1e-4*std::sqrt(ROL_EPSILON<Real>())), rtol(1e-2);
    const Real one(1), decr(1e-4), stol(std::sqrt(ROL_EPSILON<Real>()));
    const Real factor(0.5), c1(1e-4), c2(1e-2), half(0.5);
    update_primal_nd(*xnew_,x,lam);
    Real rnorm = residual_nd(*res_,*xnew_);
    const Real ctol = std::min(atol,rtol*rnorm);
    Real alpha(1), tmp(0), mu(0), rho(1), dd(0);
    int cnt = 0, maxit = 1000;
    for (cnt = 0; cnt < maxit; ++cnt) {
      // Compute Newton step
      mu  = c1*std::max(rnorm,std::sqrt(rnorm));
      rho = std::min(half,c2*std::min(std::sqrt(rnorm),rnorm)); // Unused
      solve_newton_system_nd(dlam,*res_,*xnew_,mu,rho);
      lnew_->set(lam); lnew_->axpy(-alpha, dlam);
      update_primal_nd(*xnew_,x,*lnew_);
      //tmp = residual_nd(*res_,*xnew_);
      /* Begin Solodov and Svaiter */
      rnorm = residual_nd(*res_,*xnew_);
      tmp   = dlam.dot(res_->dual());
      dd    = dlam.dot(dlam);
      /* End Solodov and Svaiter */
      // Perform backtracking line search
      //while ( tmp > (one-decr*alpha)*rnorm && alpha > stol ) {
      /* Begin Solodov and Svaiter */
      while ( tmp < decr*(one-rho)*mu*dd && alpha > stol ) {
      /* End Solodov and Svaiter */
        alpha *= factor;
        lnew_->set(lam); lnew_->axpy(-alpha, dlam);
        update_primal_nd(*xnew_,x,*lnew_);
        //tmp = residual_nd(*res_,*xnew_);
        /* Begin Solodov and Svaiter */
        rnorm = residual_nd(*res_,*xnew_);
        tmp   = dlam.dot(res_->dual());
        /* End Solodov and Svaiter */
      }
      // Update iterate
      lam.set(*lnew_);
      //rnorm = tmp;
      /* Begin Solodov and Svaiter */
      //lam.axpy(-alpha*tmp/(rnorm*rnorm),res_->dual());
      //update_primal_nd(*xnew_,x,lam);
      //rnorm = residual_nd(*res_,*xnew_);
      /* End Solodov and Svaiter */
      //std::cout << "  cnt = " << cnt << "  rnorm = " << rnorm << "  alpha = " << alpha << "  mu = " << mu << "  rho = " << rho << std::endl;
      if (rnorm <= ctol) { // = covers the case of identically zero residual
        break;
      }
      alpha = one;
    }
    if (rnorm > ctol) {
      //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
      std::cout << ">>> ROL::PolyhedralProjection::project : Projection failed!  rnorm = " << rnorm << "  rtol = " << ctol << std::endl;
    }
    x.set(*xnew_);
  }

}; // class PolyhedralProjection

} // namespace ROL

#endif
