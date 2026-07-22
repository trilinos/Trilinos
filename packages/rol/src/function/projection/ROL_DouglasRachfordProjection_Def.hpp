// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DOUGLASRACHFORDPROJECTION_DEF_H
#define ROL_DOUGLASRACHFORDPROJECTION_DEF_H

namespace ROL {

template<typename Real>
DouglasRachfordProjection<Real>::DouglasRachfordProjection(const Vector<Real>               &xprim,
                                                           const Vector<Real>               &xdual,
                                                           const Ptr<BoundConstraint<Real>> &bnd,
                                                           const Ptr<Constraint<Real>>      &con,
                                                           const Vector<Real>               &mul,
                                                           const Vector<Real>               &res)
  : PolyhedralProjection<Real>(xprim,xdual,bnd,con,mul,res),
    DEFAULT_atol_      (1e-2*std::sqrt(ROL_EPSILON<Real>()*std::sqrt(ROL_EPSILON<Real>()))),
    DEFAULT_rtol_      (1e-2*std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_maxit_     (10000),
    DEFAULT_verbosity_ (0),
    DEFAULT_alpha1_    (0.5),
    DEFAULT_gamma_     (10.0),
    DEFAULT_t0_        (1.9),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    maxit_     (DEFAULT_maxit_),
    verbosity_ (DEFAULT_verbosity_),
    alpha1_    (DEFAULT_alpha1_),
    alpha2_    (1.0-alpha1_),
    gamma_     (DEFAULT_gamma_),
    t0_        (DEFAULT_t0_) {
  dim_ = mul.dimension();
  tmp_ = xprim.clone();
  y_   = xprim.clone();
  q_   = xprim.clone();
  p_   = xprim.clone();
  z_   = xdual.clone();
  if (dim_ == 1) {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    xprim_->zero();
    con_->update(*xprim_,UpdateType::Temp);
    con_->value(*res_,*xprim_,tol);
    b_ = res_->dot(*res_->basis(0));
    mul_->setScalar(static_cast<Real>(1));
    con_->applyAdjointJacobian(*z_,*mul_,xprim,tol);
    xprim_->set(z_->dual());
    cdot_ = xprim_->dot(*xprim_);
  }
  z_->zero();
}

template<typename Real>
DouglasRachfordProjection<Real>::DouglasRachfordProjection(const Vector<Real>               &xprim,
                                                           const Vector<Real>               &xdual,
                                                           const Ptr<BoundConstraint<Real>> &bnd,
                                                           const Ptr<Constraint<Real>>      &con,
                                                           const Vector<Real>               &mul,
                                                           const Vector<Real>               &res,
                                                           ParameterList                    &list)
  : DouglasRachfordProjection<Real>(xprim,xdual,bnd,con,mul,res) {
  atol_      = list.sublist("General").sublist("Polyhedral Projection").get("Absolute Tolerance", DEFAULT_atol_);
  rtol_      = list.sublist("General").sublist("Polyhedral Projection").get("Relative Tolerance", DEFAULT_rtol_);
  maxit_     = list.sublist("General").sublist("Polyhedral Projection").get("Iteration Limit",    DEFAULT_maxit_);
  verbosity_ = list.sublist("General").get("Output Level", DEFAULT_verbosity_);
  alpha1_    = list.sublist("General").sublist("Polyhedral Projection").sublist("Douglas-Rachford").get("Constraint Weight", DEFAULT_alpha1_);
  alpha2_    = static_cast<Real>(1)-alpha1_;
  gamma_     = list.sublist("General").sublist("Polyhedral Projection").sublist("Douglas-Rachford").get("Penalty Parameter", DEFAULT_gamma_);
  t0_        = list.sublist("General").sublist("Polyhedral Projection").sublist("Douglas-Rachford").get("Relaxation Parameter", DEFAULT_t0_);
}

template<typename Real>
void DouglasRachfordProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    project_DouglasRachford(x, stream);
  }
}

template<typename Real>
Real DouglasRachfordProjection<Real>::residual_1d(const Vector<Real> &x) const {
  return xprim_->dot(x) + b_;
}

template<typename Real>
void DouglasRachfordProjection<Real>::residual_nd(Vector<Real> &r, const Vector<Real> &y) const {
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  con_->update(y,UpdateType::Temp);
  con_->value(r,y,tol);
}

template<typename Real>
void DouglasRachfordProjection<Real>::project_bnd(Vector<Real> &x, const Vector<Real> &y) const {
  x.set(y);
  bnd_->project(x);
}

template<typename Real>
void DouglasRachfordProjection<Real>::project_con(Vector<Real> &x, const Vector<Real> &y) const {
  if (dim_ == 1) {
    Real rhs = residual_1d(y);
    Real lam = -rhs/cdot_;
    x.set(y);
    x.axpy(lam,*xprim_);
  }
  else {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    residual_nd(*res_,y);
    con_->solveAugmentedSystem(x,*mul_,*z_,*res_,y,tol);
    x.scale(static_cast<Real>(-1));
    x.plus(y);
  }
}

template<typename Real>
void DouglasRachfordProjection<Real>::project_DouglasRachford(Vector<Real> &x, std::ostream &stream) const {
  const Real one(1), two(2), xnorm(x.norm()), ctol(std::min(atol_,rtol_*xnorm));
  Real rnorm(0);
  p_->zero(); q_->zero(); y_->set(x);
  std::ios_base::fmtflags streamFlags(stream.flags());
  if (verbosity_ > 2) {
    stream << std::scientific << std::setprecision(6);
    stream << std::endl;
    stream << " Polyhedral Projection using Douglas Rachford Splitting" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "error";
    stream << std::setw(15) << std::left << "tol";
    stream << std::endl;
  }
  for (int cnt=0; cnt < maxit_; ++cnt) {
    // Constraint projection
    tmp_->set(*y_);
    tmp_->axpy(alpha1_*gamma_,x);
    tmp_->scale(one/(alpha1_*gamma_+one));
    project_con(*p_,*tmp_);
    // Bounds projection
    tmp_->zero();
    tmp_->axpy(two,*p_);
    tmp_->axpy(-one,*y_);
    tmp_->axpy(alpha2_*gamma_,x);
    tmp_->scale(one/(alpha2_*gamma_+one));
    project_bnd(*q_,*tmp_);
    // Dual update
    tmp_->set(*q_);
    tmp_->axpy(-one,*p_);
    rnorm = tmp_->norm();
    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << rnorm;
      stream << std::setw(15) << std::left << ctol;
      stream << std::endl;
    }
    if (rnorm <= ctol) break;
    y_->axpy(t0_,*tmp_);
  }
  if (verbosity_ > 2) stream << std::endl;
  x.set(*q_);
  if (rnorm > ctol) {
    //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
    stream << ">>> ROL::PolyhedralProjection::project : Projection may be inaccurate!  rnorm = ";
    stream << rnorm << "  rtol = " << ctol << std::endl;
  }
  stream.flags(streamFlags);
}

} // namespace ROL

#endif
