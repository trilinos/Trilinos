// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DAIFLETCHERPROJECTION_DEF_H
#define ROL_DAIFLETCHERPROJECTION_DEF_H

namespace ROL {

template<typename Real>
DaiFletcherProjection<Real>::DaiFletcherProjection(const Vector<Real>               &xprim,
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
DaiFletcherProjection<Real>::DaiFletcherProjection(const Vector<Real>               &xprim,
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
void DaiFletcherProjection<Real>::initialize(const Vector<Real>               &xprim,
                                             const Vector<Real>               &xdual,
                                             const Ptr<BoundConstraint<Real>> &bnd,
                                             const Ptr<Constraint<Real>>      &con,
                                             const Vector<Real>               &mul,
                                             const Vector<Real>               &res) {
  dim_ = mul.dimension();
  ROL_TEST_FOR_EXCEPTION(dim_!=1,std::logic_error,
    ">>> ROL::DaiFletcherProjection : The range of the linear constraint must be one dimensional!");
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
void DaiFletcherProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    Px_->set(x); bnd_->project(*Px_);
    mul1_  = -residual(*Px_)/cdot_;
    //mul1_  = -residual(x)/cdot_;
    //mul1_  = static_cast<Real>(0);
    dlam1_ = static_cast<Real>(2);
    //dlam1_ = static_cast<Real>(1)+std::abs(mul1_);
    project_df(x, mul1_, dlam1_, stream);
    mul_->setScalar(mul1_);
  }
}

template<typename Real>
Real DaiFletcherProjection<Real>::residual(const Vector<Real> &x) const {
  return xprim_->dot(x) + b_;
}

template<typename Real>
void DaiFletcherProjection<Real>::update_primal(Vector<Real> &y, const Vector<Real> &x, const Real lam) const {
  y.set(x);
  y.axpy(lam,*xprim_);
  bnd_->project(y);
}

template<typename Real>
void DaiFletcherProjection<Real>::project_df(Vector<Real> &x, Real &lam, Real &dlam, std::ostream &stream) const {
  const Real zero(0), one(1), two(2), c1(0.1), c2(0.75), c3(0.25);
  Real lamLower(0), lamUpper(0), lamNew(0), res(0), resLower(0), resUpper(0), s(0);
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
    stream << " Polyhedral Projection using the Dai-Fletcher Algorithm" << std::endl;
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
  rtol = ctol_*std::max(one,std::min(std::abs(resLower),std::abs(resUpper)));
  //s    = one - resLower / resUpper;
  //dlam = (lamUpper - lamLower) / s;
  //lam  = lamUpper - dlam;
  s    = (resUpper - resLower) / resUpper;
  lam  = (resUpper * lamLower - resLower * lamUpper) / (resUpper - resLower);
  dlam = lamUpper - lam;
  update_primal(*xnew_,x,lam);
  res  = residual(*xnew_);
  cnt  = 0;
  if (verbosity_ > 2) {
    stream << std::endl;
    stream << "  Secant Phase" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "lam";
    stream << std::setw(15) << std::left << "res";
    stream << std::setw(15) << std::left << "stepsize";
    stream << std::setw(15) << std::left << "rtol";
    stream << std::setw(15) << std::left << "lbnd";
    stream << std::setw(15) << std::left << "lres";
    stream << std::setw(15) << std::left << "ubnd";
    stream << std::setw(15) << std::left << "ures";
    stream << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << cnt;
    stream << std::setw(15) << std::left << lam;
    stream << std::setw(15) << std::left << res;
    stream << std::setw(15) << std::left << dlam;
    stream << std::setw(15) << std::left << rtol;
    stream << std::setw(15) << std::left << lamLower;
    stream << std::setw(15) << std::left << resLower;
    stream << std::setw(15) << std::left << lamUpper;
    stream << std::setw(15) << std::left << resUpper;
    stream << std::endl;
  }
  for (cnt = 1; cnt < maxit_; cnt++) {
    // Exit if residual or bracket length are sufficiently small
    if ( std::abs(res) <= rtol ||
         std::abs(lamUpper-lamLower) < ltol_*std::max(std::abs(lamUpper),std::abs(lamLower)) ) {
      break;
    }

    if ( res > zero ) {
      if ( s <= two ) {
        lamUpper = lam;
        resUpper = res;
        //s        = one - resLower / resUpper;
        //dlam     = (lamUpper - lamLower) / s;
        //lam      = lamUpper - dlam;
        s        = (resUpper - resLower) / resUpper;
        lam      = (lamLower * resUpper - lamUpper * resLower) / (resUpper - resLower);
        dlam     = lamUpper - lam;
      }
      else {
        //s        = std::max(resUpper / res - one, c1);
        //dlam     = (lamUpper - lam) / s;
        //lamNew   = std::max(lam - dlam, c2*lamLower + c3*lam);
        if (resUpper <= (c1+one)*res) {
          dlam   = (lamUpper - lam) / c1;
          lamNew = std::max(lam - dlam, c2*lamLower + c3*lam);
        }
        else {
          lamNew = std::max((lam * resUpper - lamUpper * res) / (resUpper - res),
                            c2*lamLower + c3*lam);
          dlam   = lam - lamNew;
        }
        lamUpper = lam;
        resUpper = res;
        lam      = lamNew;
        s        = (lamUpper - lamLower) / (lamUpper - lam);
      }
    }
    else {
      if ( s >= two ) {
        lamLower = lam;
        resLower = res;
        //s        = one - resLower / resUpper;
        //dlam     = (lamUpper - lamLower) / s;
        //lam      = lamUpper - dlam;
        s        = (resUpper - resLower) / resUpper;
        lam      = (lamLower * resUpper - lamUpper * resLower) / (resUpper - resLower);
        dlam     = lamUpper - lam;
      }
      else {
        //s        = std::max(resLower / res - one, c1);
        //dlam     = (lam + lamLower) / s;
        //lamNew   = std::min(lam + dlam, c2*lamUpper + c3*lam);
        if (resLower >= (c1+one)*res) {
          dlam   = (lam - lamLower) / c1;
          lamNew = std::max(lam + dlam, c2*lamUpper + c3*lam);
        }
        else {
          lamNew = std::max((lamLower * res - lam * resLower) / (res - resLower),
                            c2*lamUpper + c3*lam);
          dlam   = lamNew - lamLower;
        }
        lamLower = lam;
        resLower = res;
        lam      = lamNew;
        s        = (lamUpper - lamLower) / (lamUpper - lam);
      }
    }
    update_primal(*xnew_,x,lam);
    res = residual(*xnew_);

    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << lam;
      stream << std::setw(15) << std::left << res;
      stream << std::setw(15) << std::left << dlam;
      stream << std::setw(15) << std::left << rtol;
      stream << std::setw(15) << std::left << lamLower;
      stream << std::setw(15) << std::left << resLower;
      stream << std::setw(15) << std::left << lamUpper;
      stream << std::setw(15) << std::left << resUpper;
      stream << std::endl;
    }
  }
  if (verbosity_ > 2) {
    stream << std::endl;
  }
  // Return projection
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
