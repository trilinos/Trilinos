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
    DEFAULT_atol_      (1e-4*std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_rtol_      (1e-2),
    DEFAULT_ltol_      (ROL_EPSILON<Real>()),
    DEFAULT_maxit_     (5000),
    DEFAULT_verbosity_ (0),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    ltol_      (DEFAULT_ltol_),
    maxit_     (DEFAULT_maxit_),
    verbosity_ (DEFAULT_verbosity_) {
  dim_ = mul.dimension();
  if (dim_ != 1) {
    throw Exception::NotImplemented(">>> ROL::DaiFletcherProjection : The range of the linear constraint must be one dimensional!");
  }
  xnew_  = xprim.clone();
  mul1_  = static_cast<Real>(0);
  dlam1_ = static_cast<Real>(2);
  // con.value(x) = xprim_->dot(x) + b_
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  xprim_->zero();
  con_->value(*res_,*xprim_,tol);
  b_ = res_->dot(*res_->basis(0));
  mul_->setScalar(static_cast<Real>(1));
  con_->applyAdjointJacobian(*xdual_,*mul_,xprim,tol);
  xprim_->set(xdual_->dual());
  cdot_ = xprim_->dot(*xprim_);
}

template<typename Real>
DaiFletcherProjection<Real>::DaiFletcherProjection(const Vector<Real>               &xprim,
                                                   const Vector<Real>               &xdual,
                                                   const Ptr<BoundConstraint<Real>> &bnd,
                                                   const Ptr<Constraint<Real>>      &con,
                                                   const Vector<Real>               &mul,
                                                   const Vector<Real>               &res,
                                                   ParameterList                    &list)
  : DaiFletcherProjection<Real>(xprim,xdual,bnd,con,mul,res) {
  atol_      = list.sublist("General").sublist("Polyhedral Projection").get("Absolute Tolerance",   DEFAULT_atol_);
  rtol_      = list.sublist("General").sublist("Polyhedral Projection").get("Relative Tolerance",   DEFAULT_rtol_);
  ltol_      = list.sublist("General").sublist("Polyhedral Projection").get("Multiplier Tolerance", DEFAULT_ltol_);
  maxit_     = list.sublist("General").sublist("Polyhedral Projection").get("Iteration Limit",      DEFAULT_maxit_);
  verbosity_ = list.sublist("General").get("Output Level", DEFAULT_verbosity_);
}

template<typename Real>
void DaiFletcherProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    mul1_  = -residual(x)/cdot_;
    //mul1_  = static_cast<Real>(0);
    dlam1_ = static_cast<Real>(2);
    //dlam1_ = static_cast<Real>(1)+std::abs(mul1_);
    project_df(x, mul1_, dlam1_, stream);
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
  Real lam1(0), lam2(0), lam3(0), r(0), r0(0), r1(0), r2(0), s(0);
  int cnt(0);
  // Set residual tolerance to min(atol,rtol*abs(<c,P(y)>+b))
  update_primal(*xnew_,x,lam1);
  r0 = residual(*xnew_);
  if (r0 == zero) {
    lam = lam1;
    x.set(*xnew_);
    return;
  }
  // Compute initial residual
  update_primal(*xnew_,x,lam);
  r = residual(*xnew_);
  if (r == zero) {
    x.set(*xnew_);
    return;
  }
  const Real ctol = std::min(atol_,rtol_*std::max(std::abs(r),std::abs(r0)));
  std::ios_base::fmtflags streamFlags(stream.flags());
  if (verbosity_ > 1) {
    stream << std::scientific << std::setprecision(6);
    stream << std::endl;
    stream << " Polyhedral Projection using the Dai-Fletcher Algorithm" << std::endl;
    stream << "  Bracketing Phase" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "Lower lam";
    stream << std::setw(15) << std::left << "Lower res";
    stream << std::setw(15) << std::left << "Upper lam";
    stream << std::setw(15) << std::left << "Upper res";
    stream << std::endl;
  }
  // Bracketing phase
  if ( r < zero ) {
    lam1 = lam;
    r1   = r;
    lam += dlam;
    update_primal(*xnew_,x,lam);
    r    = residual(*xnew_);
    if (verbosity_ > 1) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << lam1;
      stream << std::setw(15) << std::left << r1;
      stream << std::setw(15) << std::left << lam2;
      stream << std::setw(15) << std::left << r2;
      stream << std::endl;
    }
    while ( r < zero && std::abs(r) > ctol && cnt < maxit_ ) {
      lam1  = lam;
      r1    = r;
      s     = std::max(r1/r-one,c1);
      dlam += dlam/s;
      lam  += dlam;
      update_primal(*xnew_,x,lam);
      r     = residual(*xnew_);
      cnt++;
      if (verbosity_ > 1) {
        stream << "  ";
        stream << std::setw(6)  << std::left << cnt;
        stream << std::setw(15) << std::left << lam1;
        stream << std::setw(15) << std::left << r1;
        stream << std::setw(15) << std::left << lam2;
        stream << std::setw(15) << std::left << r2;
        stream << std::endl;
      }
    }
    lam2 = lam;
    r2   = r;
  }
  else {
    lam2 = lam;
    r2   = r;
    lam -= dlam;
    update_primal(*xnew_,x,lam);
    r    = residual(*xnew_);
    if (verbosity_ > 1) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << lam1;
      stream << std::setw(15) << std::left << r1;
      stream << std::setw(15) << std::left << lam2;
      stream << std::setw(15) << std::left << r2;
      stream << std::endl;
    }
    while ( r > zero && std::abs(r) > ctol && cnt < maxit_ ) {
      lam2  = lam;
      r2    = r;
      s     = std::max(r2/r-one,c1);
      dlam += dlam/s;
      lam  -= dlam;
      update_primal(*xnew_,x,lam);
      r     = residual(*xnew_);
      cnt++;
      if (verbosity_ > 1) {
        stream << "  ";
        stream << std::setw(6)  << std::left << cnt;
        stream << std::setw(15) << std::left << lam1;
        stream << std::setw(15) << std::left << r1;
        stream << std::setw(15) << std::left << lam2;
        stream << std::setw(15) << std::left << r2;
        stream << std::endl;
      }
    }
    lam1 = lam;
    r1   = r;
  }

  // Secant phase
  s     = one - r1/r2;
  dlam /= s;
  lam   = lam2-dlam;
  update_primal(*xnew_,x,lam);
  r     = residual(*xnew_);
  cnt   = 0;
  if (verbosity_ > 1) {
    stream << std::endl;
    stream << "  Secant Phase" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "res";
    stream << std::setw(15) << std::left << "stepsize";
    stream << std::setw(15) << std::left << "residual tol";
    stream << std::endl;
  }
  while ( std::abs(r) > ctol && std::abs(dlam) > ltol_ && cnt < maxit_ ) {
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
    }
    update_primal(*xnew_,x,lam);
    r = residual(*xnew_);
    if (verbosity_ > 1) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << r;
      stream << std::setw(15) << std::left << dlam;
      stream << std::setw(15) << std::left << ctol;
      stream << std::endl;
    }
    cnt++;
  }
  if (verbosity_ > 1) {
    stream << std::endl;
  }
  // Return projection
  x.set(*xnew_);
  if (std::abs(r) > ctol) {
    //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
    stream << ">>> ROL::PolyhedralProjection::project : Projection may be inaccurate!  rnorm = ";
    stream << std::abs(r) << "  rtol = " << ctol << std::endl;
  }
  stream.flags(streamFlags);
}

} // namespace ROL

#endif
