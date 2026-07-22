// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BRENTSPROJECTION_DEF_H
#define ROL_BRENTSPROJECTION_DEF_H

namespace ROL {

template<typename Real>
BrentsProjection<Real>::BrentsProjection(const Vector<Real>               &xprim,
                                           const Vector<Real>               &xdual,
                                           const Ptr<BoundConstraint<Real>> &bnd,
                                           const Ptr<Constraint<Real>>      &con,
                                           const Vector<Real>               &mul,
                                           const Vector<Real>               &res)
  : PolyhedralProjection<Real>(xprim,xdual,bnd,con,mul,res),
    DEFAULT_atol_      (std::sqrt(ROL_EPSILON<Real>()*std::sqrt(ROL_EPSILON<Real>()))),
    DEFAULT_rtol_      (std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_ltol_      (ROL_EPSILON<Real>()),
    DEFAULT_maxit_     (5000),
    DEFAULT_verbosity_ (0),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    ltol_      (DEFAULT_ltol_),
    maxit_     (DEFAULT_maxit_),
    verbosity_ (DEFAULT_verbosity_) {
  initialize(xprim,xdual,bnd,con,mul,res);
}

template<typename Real>
BrentsProjection<Real>::BrentsProjection(const Vector<Real>               &xprim,
                                           const Vector<Real>               &xdual,
                                           const Ptr<BoundConstraint<Real>> &bnd,
                                           const Ptr<Constraint<Real>>      &con,
                                           const Vector<Real>               &mul,
                                           const Vector<Real>               &res,
                                           ParameterList                    &list)
  : PolyhedralProjection<Real>(xprim,xdual,bnd,con,mul,res),
    DEFAULT_atol_      (std::sqrt(ROL_EPSILON<Real>()*std::sqrt(ROL_EPSILON<Real>()))),
    DEFAULT_rtol_      (std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_ltol_      (ROL_EPSILON<Real>()),
    DEFAULT_maxit_     (5000),
    DEFAULT_verbosity_ (0),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    ltol_      (DEFAULT_ltol_),
    maxit_     (DEFAULT_maxit_),
    verbosity_ (DEFAULT_verbosity_) {
  atol_      = list.sublist("General").sublist("Polyhedral Projection").get("Absolute Tolerance",   DEFAULT_atol_);
  rtol_      = list.sublist("General").sublist("Polyhedral Projection").get("Relative Tolerance",   DEFAULT_rtol_);
  ltol_      = list.sublist("General").sublist("Polyhedral Projection").get("Multiplier Tolerance", DEFAULT_ltol_);
  maxit_     = list.sublist("General").sublist("Polyhedral Projection").get("Iteration Limit",      DEFAULT_maxit_);
  verbosity_ = list.sublist("General").get("Output Level", DEFAULT_verbosity_);
  initialize(xprim,xdual,bnd,con,mul,res);
}

template<typename Real>
void BrentsProjection<Real>::initialize(const Vector<Real>               &xprim,
                                         const Vector<Real>               &xdual,
                                         const Ptr<BoundConstraint<Real>> &bnd,
                                         const Ptr<Constraint<Real>>      &con,
                                         const Vector<Real>               &mul,
                                         const Vector<Real>               &res) {
  dim_ = mul.dimension();
  ROL_TEST_FOR_EXCEPTION(dim_!=1,std::logic_error,
    ">>> ROL::BrentsProjection : The range of the linear constraint must be one dimensional!");
  xnew_  = xprim.clone();
  Px_    = xprim.clone();
  mul1_  = static_cast<Real>(0);
  dlam1_ = static_cast<Real>(2);
  // con.value(x) = xprim_->dot(x) + b_
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  xprim_->zero();
  con_->update(*xprim_,UpdateType::Temp);
  con_->value(*res_,*xprim_,tol);
  b_ = res_->dot(*res_->basis(0));
  mul_->setScalar(static_cast<Real>(1));
  con_->applyAdjointJacobian(*xdual_,*mul_,xprim,tol);
  xprim_->set(xdual_->dual());
  cdot_ = xprim_->dot(*xprim_);
  // Set tolerance
  //xnew_->zero();
  //bnd_->project(*xnew_);
  //Real res0 = std::abs(residual(*xnew_));
  Real resl = ROL_INF<Real>(), resu = ROL_INF<Real>();
  if (bnd_->isLowerActivated()) resl = residual(*bnd_->getLowerBound());
  if (bnd_->isUpperActivated()) resu = residual(*bnd_->getUpperBound());
  Real res0 = std::max(resl,resu);
  if (res0 < atol_) res0 = static_cast<Real>(1);
  ctol_ = std::min(atol_,rtol_*res0);
}

template<typename Real>
void BrentsProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    mul1_  = -residual(x)/cdot_;
    //mul1_  = static_cast<Real>(0);
    dlam1_ = static_cast<Real>(2);
    //dlam1_ = static_cast<Real>(1)+std::abs(mul1_);
    project_df(x, mul1_, dlam1_, stream);
    mul_->setScalar(mul1_);
  }
}

template<typename Real>
Real BrentsProjection<Real>::residual(const Vector<Real> &x) const {
  return xprim_->dot(x) + b_;
}

template<typename Real>
void BrentsProjection<Real>::update_primal(Vector<Real> &y, const Vector<Real> &x, const Real lam) const {
  y.set(x);
  y.axpy(lam,*xprim_);
  bnd_->project(y);
}

template<typename Real>
void BrentsProjection<Real>::project_df(Vector<Real> &x, Real &lam, Real &dlam, std::ostream &stream) const {
  const Real zero(0), one(1), c1(0.1);
  Real lamLower(0), lamUpper(0), res(0), resLower(0), resUpper(0), s(0);
  Real rtol = ctol_;
  int cnt(0);
  // Compute initial residual
  update_primal(*xnew_,x,lam);
  res = residual(*xnew_);
  if (res == zero) {
    x.set(*xnew_);
    return;
  }
  std::ios_base::fmtflags streamFlags(stream.flags());
  if (verbosity_ > 2) {
    stream << std::scientific << std::setprecision(6);
    stream << std::endl;
    stream << " Polyhedral Projection using Brents' Algorithm" << std::endl;
    stream << "  Bracketing Phase" << std::endl;
  }
  // Bracketing phase
  if ( res < zero ) {
    lamLower = lam;
    resLower = res;
    lam     += dlam;
    update_primal(*xnew_,x,lam);
    res      = residual(*xnew_);
    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << "iter";
      stream << std::setw(15) << std::left << "lam";
      stream << std::setw(15) << std::left << "res";
      stream << std::setw(15) << std::left << "lower lam";
      stream << std::setw(15) << std::left << "lower res";
      stream << std::endl;
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << lam;
      stream << std::setw(15) << std::left << res;
      stream << std::setw(15) << std::left << lamLower;
      stream << std::setw(15) << std::left << resLower;
      stream << std::endl;
    }
    while ( res < zero && std::abs(res) > rtol && cnt < maxit_ ) {
      s         = std::max(resLower/res-one,c1);
      dlam     += dlam/s;
      lamLower  = lam;
      resLower  = res;
      lam      += dlam;
      update_primal(*xnew_,x,lam);
      res       = residual(*xnew_);
      cnt++;
      if (verbosity_ > 2) {
        stream << "  ";
        stream << std::setw(6)  << std::left << cnt;
        stream << std::setw(15) << std::left << lam;
        stream << std::setw(15) << std::left << res;
        stream << std::setw(15) << std::left << lamLower;
        stream << std::setw(15) << std::left << resLower;
        stream << std::endl;
      }
    }
    lamUpper = lam;
    resUpper = res;
  }
  else {
    lamUpper = lam;
    resUpper = res;
    lam     -= dlam;
    update_primal(*xnew_,x,lam);
    res      = residual(*xnew_);
    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << "iter";
      stream << std::setw(15) << std::left << "lam";
      stream << std::setw(15) << std::left << "res";
      stream << std::setw(15) << std::left << "upper lam";
      stream << std::setw(15) << std::left << "upper res";
      stream << std::endl;
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << lam;
      stream << std::setw(15) << std::left << res;
      stream << std::setw(15) << std::left << lamUpper;
      stream << std::setw(15) << std::left << resUpper;
      stream << std::endl;
    }
    while ( res > zero && std::abs(res) > rtol && cnt < maxit_ ) {
      s         = std::max(resUpper/res-one,c1);
      dlam     += dlam/s;
      lamUpper  = lam;
      resUpper  = res;
      lam      -= dlam;
      update_primal(*xnew_,x,lam);
      res       = residual(*xnew_);
      cnt++;
      if (verbosity_ > 2) {
        stream << "  ";
        stream << std::setw(6)  << std::left << cnt;
        stream << std::setw(15) << std::left << lam;
        stream << std::setw(15) << std::left << res;
        stream << std::setw(15) << std::left << lamUpper;
        stream << std::setw(15) << std::left << resUpper;
        stream << std::endl;
      }
    }
    lamLower = lam;
    resLower = res;
  }
  if (verbosity_ > 2) {
    stream << "  Bracket: ";
    stream << std::setw(15) << std::left << lamLower;
    stream << std::setw(15) << std::left << lamUpper;
    stream << std::endl;
  }

  // Secant phase
  //rtol = ctol_*std::max(one,std::min(std::abs(resLower),std::abs(resUpper)));
  cnt  = 0;
  if (verbosity_ > 2) {
    stream << std::endl;
    stream << "  Brents' Phase" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "rtol";
    stream << std::setw(15) << std::left << "lam";
    stream << std::setw(15) << std::left << "res";
    stream << std::setw(15) << std::left << "lam low";
    stream << std::setw(15) << std::left << "res low";
    stream << std::setw(15) << std::left << "lam up";
    stream << std::setw(15) << std::left << "res up";
    stream << std::endl;
  }
  const Real half(0.5), two(2), three(3);
  const Real eps(ROL_EPSILON<Real>()), tol0(rtol); // tol0(1e1*eps);
  Real d1(1), d2(1), tol(1);
  Real p(0), q(0), r(0), m(0);
  lam = lamUpper; res = resUpper;
  update_primal(*xnew_,x,lamUpper);
  for (cnt = 0; cnt < maxit_; cnt++) {
    if ((resUpper > zero && res > zero) || (resUpper <= zero && res <= zero)) {
      lam = lamLower; res = resLower;
      d1 = lamUpper-lamLower; d2 = d1;
    }
    if (std::abs(res) < std::abs(resUpper)) {
      lamLower = lamUpper; lamUpper = lam; lam = lamLower;
      resLower = resUpper; resUpper = res; res = resLower;
    }
    tol = two*eps*std::abs(lamUpper) + half*tol0;
    m = half*(lam - lamUpper);
    if (std::abs(m) <= tol || std::abs(resUpper) <= rtol) break;
    if (std::abs(d2) < tol || std::abs(resLower) <= std::abs(resUpper)) {
      d1 = m; d2 = d1;
    }
    else {
      s = resUpper/resLower;
      if (lamLower == lam) {
        p = two*m*s;
        q = one-s;
      }
      else {
        q = resLower/res;
        r = resUpper/res;
        p = s*(two*m*q*(q-r)-(lamUpper-lamLower)*(r-one));
        q = (q-one)*(r-one)*(s-one);
      }
      if (p > zero) q = -q;
      else          p = -p;
      if (two*p < three*m*q-std::abs(tol*q) && p < std::abs(half*d2*q)) {
        d2 = d1; d1 = p/q;
      }
      else {
        d1 = m; d2 = d1;
      }
    }
    lamLower = lamUpper; resLower = resUpper;
    if (std::abs(d1) > tol) lamUpper += d1;
    else if (m > zero)      lamUpper += tol;
    else                    lamUpper -= tol;
    update_primal(*xnew_,x,lamUpper);
    resUpper = residual(*xnew_);

    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << rtol;
      stream << std::setw(15) << std::left << lam;
      stream << std::setw(15) << std::left << res;
      stream << std::setw(15) << std::left << lamLower;
      stream << std::setw(15) << std::left << resLower;
      stream << std::setw(15) << std::left << lamUpper;
      stream << std::setw(15) << std::left << resUpper;
      stream << std::endl;
    }
  }
  if (verbosity_ > 2) {
    if (cnt < maxit_) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << rtol;
      stream << std::setw(15) << std::left << lam;
      stream << std::setw(15) << std::left << res;
      stream << std::setw(15) << std::left << lamLower;
      stream << std::setw(15) << std::left << resLower;
      stream << std::setw(15) << std::left << lamUpper;
      stream << std::setw(15) << std::left << resUpper;
      stream << std::endl;
    }
    stream << std::endl;
  }
  // Return projection
  res = resUpper;
  x.set(*xnew_);
  if (std::abs(res) > rtol ) {
    //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
    stream << ">>> ROL::PolyhedralProjection::project : Projection may be inaccurate!  rnorm = ";
    stream << std::abs(res) << "  rtol = " << rtol << std::endl;
  }
  stream.flags(streamFlags);
}

} // namespace ROL

#endif
