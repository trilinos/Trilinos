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


#ifndef ROL_SEMISMOOTHNEWTONPROJECTION_DEF_H
#define ROL_SEMISMOOTHNEWTONPROJECTION_DEF_H

namespace ROL {

template<typename Real>
SemismoothNewtonProjection<Real>::SemismoothNewtonProjection(const Vector<Real>               &xprim,
                                                             const Vector<Real>               &xdual,
                                                             const Ptr<BoundConstraint<Real>> &bnd,
                                                             const Ptr<Constraint<Real>>      &con,
                                                             const Vector<Real>               &mul,
                                                             const Vector<Real>               &res)
  : PolyhedralProjection<Real>(xprim,xdual,bnd,con,mul,res),
    DEFAULT_atol_      (std::sqrt(ROL_EPSILON<Real>()*std::sqrt(ROL_EPSILON<Real>()))),
    DEFAULT_rtol_      (std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_stol_      (std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_decr_      (1e-4),
    DEFAULT_factor_    (0.5),
    DEFAULT_regscale_  (1e-4),
    DEFAULT_errscale_  (1e-2),
    DEFAULT_maxit_     (5000),
    DEFAULT_lstype_    (0),
    DEFAULT_verbosity_ (0),
    DEFAULT_useproj_   (false),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    stol_      (DEFAULT_stol_),
    decr_      (DEFAULT_decr_),
    factor_    (DEFAULT_factor_),
    regscale_  (DEFAULT_regscale_),
    errscale_  (DEFAULT_errscale_),
    maxit_     (DEFAULT_maxit_),
    lstype_    (DEFAULT_lstype_),
    verbosity_ (DEFAULT_verbosity_),
    useproj_   (DEFAULT_useproj_) {
  dim_   = mul.dimension();
  xnew_  = xprim.clone();
  lnew_  = mul.clone();
  dlam_  = mul.clone();
  
  ParameterList list;
  list.sublist("General").sublist("Krylov").set("Type",               "CG");
  list.sublist("General").sublist("Krylov").set("Absolute Tolerance", 1e-6);
  list.sublist("General").sublist("Krylov").set("Relative Tolerance", 1e-4);
  list.sublist("General").sublist("Krylov").set("Iteration Limit",    dim_);
  list.sublist("General").set("Inexact Hessian-Times-A-Vector",      false);
  krylov_ = KrylovFactory<Real>(list);

  // Set tolerance
  Real resl = residual(*res_,*bnd_->getLowerBound());
  Real resu = residual(*res_,*bnd_->getUpperBound());
  Real res0 = std::max(resl,resu);
  if (res0 < atol_) {
    res0 = static_cast<Real>(1);
  }
  ctol_ = std::min(atol_,rtol_*res0);
}

template<typename Real>
SemismoothNewtonProjection<Real>::SemismoothNewtonProjection(const Vector<Real>               &xprim,
                                                             const Vector<Real>               &xdual,
                                                             const Ptr<BoundConstraint<Real>> &bnd,
                                                             const Ptr<Constraint<Real>>      &con,
                                                             const Vector<Real>               &mul,
                                                             const Vector<Real>               &res,
                                                             ParameterList                    &list)
  : SemismoothNewtonProjection<Real>(xprim,xdual,bnd,con,mul,res) {
  ParameterList &ppl = list.sublist("General").sublist("Polyhedral Projection");
  atol_      = ppl.get("Absolute Tolerance",                                              DEFAULT_atol_);
  rtol_      = ppl.get("Relative Tolerance",                                              DEFAULT_rtol_);
  stol_      = ppl.sublist("Semismooth Newton").get("Step Tolerance",                     DEFAULT_stol_);
  decr_      = ppl.sublist("Semismooth Newton").get("Sufficient Decrease Tolerance",      DEFAULT_decr_);
  factor_    = ppl.sublist("Semismooth Newton").get("Backtracking Rate",                  DEFAULT_factor_);
  regscale_  = ppl.sublist("Semismooth Newton").get("Regularization Scale",               DEFAULT_regscale_);
  errscale_  = ppl.sublist("Semismooth Newton").get("Relative Error Scale",               DEFAULT_errscale_);
  maxit_     = ppl.get("Iteration Limit",                                                 DEFAULT_maxit_);
  lstype_    = ppl.sublist("Semismooth Newton").get("Line Search Type",                   DEFAULT_lstype_);
  verbosity_ = list.sublist("General").get("Output Level",                                DEFAULT_verbosity_);
  useproj_   = ppl.sublist("Semismooth Newton").get("Project onto Separating Hyperplane", DEFAULT_useproj_);
  
  ParameterList klist;
  klist.sublist("General").sublist("Krylov") = ppl.sublist("Semismooth Newton").sublist("Krylov");
  klist.sublist("General").set("Inexact Hessian-Times-A-Vector", false);
  krylov_ = KrylovFactory<Real>(klist);
}

template<typename Real>
void SemismoothNewtonProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    project_ssn(x, *mul_, *dlam_, stream);
  }
}

template<typename Real>
Real SemismoothNewtonProjection<Real>::residual(Vector<Real> &r, const Vector<Real> &y) const {
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  con_->update(y,UpdateType::Temp);
  con_->value(r,y,tol);
  return r.norm();
}

template<typename Real>
void SemismoothNewtonProjection<Real>::solve_newton_system(Vector<Real>       &s,
                                                           const Vector<Real> &r,
                                                           const Vector<Real> &y,
                                                           const Real          mu,
                                                           const Real          rho,
                                                           int                &iter,
                                                           int                &flag) const {
  Ptr<Precond>  M = makePtr<Precond>();
  Ptr<Jacobian> J = makePtr<Jacobian>(con_,bnd_,makePtrFromRef(y),xdual_,xprim_,mu);
  krylov_->run(s,*J,r,*M,iter,flag);
}

template<typename Real>
void SemismoothNewtonProjection<Real>::update_primal(Vector<Real>       &y,
                                                     const Vector<Real> &x,
                                                     const Vector<Real> &lam) const {
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  y.set(x);
  con_->update(x,UpdateType::Temp);
  con_->applyAdjointJacobian(*xdual_,lam,x,tol);
  y.plus(xdual_->dual());
  bnd_->project(y);
}

template<typename Real>
void SemismoothNewtonProjection<Real>::project_ssn(Vector<Real> &x,
                                                   Vector<Real> &lam,
                                                   Vector<Real> &dlam,
                                                   std::ostream &stream) const {
  const Real zero(0), half(0.5), one(1);
  // Compute initial residual
  update_primal(*xnew_,x,lam);
  Real rnorm = residual(*res_,*xnew_);
  if (rnorm == zero) {
    x.set(*xnew_);
    return;
  }
  Real alpha(1), tmp(0), mu(0), rho(1), dd(0);
  int iter(0), flag(0);
  std::ios_base::fmtflags streamFlags(stream.flags());
  if (verbosity_ > 2) {
    stream << std::endl;
    stream << std::scientific << std::setprecision(6);
    stream << " Polyhedral Projection using Dual Semismooth Newton" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "rnorm";
    stream << std::setw(15) << std::left << "alpha";
    stream << std::setw(15) << std::left << "mu";
    stream << std::setw(15) << std::left << "rho";
    stream << std::setw(15) << std::left << "rtol";
    stream << std::setw(8)  << std::left << "kiter";
    stream << std::setw(8)  << std::left << "kflag";
    stream << std::endl;
  }
  for (int cnt = 0; cnt < maxit_; ++cnt) {
    // Compute Newton step
    mu  = regscale_*std::max(rnorm,std::sqrt(rnorm));
    rho = std::min(half,errscale_*std::min(std::sqrt(rnorm),rnorm));
    solve_newton_system(dlam,*res_,*xnew_,mu,rho,iter,flag);
    lnew_->set(lam); lnew_->axpy(-alpha, dlam);
    update_primal(*xnew_,x,*lnew_);
    // Perform line search
    if (lstype_ == 1) { // Usual line search for nonlinear equations
      tmp = residual(*res_,*xnew_);
      while ( tmp > (one-decr_*alpha)*rnorm && alpha > stol_ ) {
        alpha *= factor_;
        lnew_->set(lam); lnew_->axpy(-alpha, dlam);
        update_primal(*xnew_,x,*lnew_);
        tmp = residual(*res_,*xnew_);
      }
      rnorm = tmp;
    }
    else { // Default Solodov and Svaiter line search
      rnorm = residual(*res_,*xnew_);
      //tmp   = dlam.dot(res_->dual());
      tmp   = dlam.apply(*res_);
      dd    = dlam.dot(dlam);
      while ( tmp < decr_*(one-rho)*mu*dd && alpha > stol_ ) {
        alpha *= factor_;
        lnew_->set(lam); lnew_->axpy(-alpha, dlam);
        update_primal(*xnew_,x,*lnew_);
        rnorm = residual(*res_,*xnew_);
        //tmp   = dlam.dot(res_->dual());
        tmp   = dlam.apply(*res_);
      }
    }
    // Update iterate
    lam.set(*lnew_);
    // Project onto separating hyperplane
    if (useproj_) {
      lam.axpy(-alpha*tmp/(rnorm*rnorm),res_->dual());
      update_primal(*xnew_,x,lam);
      rnorm = residual(*res_,*xnew_);
    }
    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << rnorm;
      stream << std::setw(15) << std::left << alpha;
      stream << std::setw(15) << std::left << mu;
      stream << std::setw(15) << std::left << rho;
      stream << std::setw(15) << std::left << ctol_;
      stream << std::setw(8)  << std::left << iter;
      stream << std::setw(8)  << std::left << flag;
      stream << std::endl;
    }
    if (rnorm <= ctol_) break;
    alpha = one;
  }
  if (verbosity_ > 2) {
    stream << std::endl;
  }
  if (rnorm > ctol_) {
    //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
    stream << ">>> ROL::PolyhedralProjection::project : Projection may be inaccurate!  rnorm = ";
    stream << rnorm << "  rtol = " << ctol_ << std::endl;
  }
  x.set(*xnew_);
  stream.flags(streamFlags);
}

} // namespace ROL

#endif
