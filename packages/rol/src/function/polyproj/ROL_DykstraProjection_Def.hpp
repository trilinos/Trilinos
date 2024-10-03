// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DYKSTRAPROJECTION_DEF_H
#define ROL_DYKSTRAPROJECTION_DEF_H

namespace ROL {

template<typename Real>
DykstraProjection<Real>::DykstraProjection(const Vector<Real>               &xprim,
                                           const Vector<Real>               &xdual,
                                           const Ptr<BoundConstraint<Real>> &bnd,
                                           const Ptr<Constraint<Real>>      &con,
                                           const Vector<Real>               &mul,
                                           const Vector<Real>               &res)
  : PolyhedralProjection<Real>(xprim,xdual,bnd,con,mul,res),
    DEFAULT_atol_      (std::sqrt(ROL_EPSILON<Real>()*std::sqrt(ROL_EPSILON<Real>()))),
    DEFAULT_rtol_      (std::sqrt(ROL_EPSILON<Real>())),
    DEFAULT_maxit_     (10000),
    DEFAULT_verbosity_ (0),
    atol_      (DEFAULT_atol_),
    rtol_      (DEFAULT_rtol_),
    maxit_     (DEFAULT_maxit_),
    verbosity_ (DEFAULT_verbosity_) {
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
DykstraProjection<Real>::DykstraProjection(const Vector<Real>               &xprim,
                                           const Vector<Real>               &xdual,
                                           const Ptr<BoundConstraint<Real>> &bnd,
                                           const Ptr<Constraint<Real>>      &con,
                                           const Vector<Real>               &mul,
                                           const Vector<Real>               &res,
                                           ParameterList                    &list)
  : DykstraProjection<Real>(xprim,xdual,bnd,con,mul,res) {
  atol_      = list.sublist("General").sublist("Polyhedral Projection").get("Absolute Tolerance", DEFAULT_atol_);
  rtol_      = list.sublist("General").sublist("Polyhedral Projection").get("Relative Tolerance", DEFAULT_rtol_);
  maxit_     = list.sublist("General").sublist("Polyhedral Projection").get("Iteration Limit",    DEFAULT_maxit_);
  verbosity_ = list.sublist("General").get("Output Level", DEFAULT_verbosity_);
}

template<typename Real>
void DykstraProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    project_Dykstra(x, stream);
  }
}

template<typename Real>
Real DykstraProjection<Real>::residual_1d(const Vector<Real> &x) const {
  return xprim_->dot(x) + b_;
}

template<typename Real>
void DykstraProjection<Real>::residual_nd(Vector<Real> &r, const Vector<Real> &y) const {
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  con_->update(y,UpdateType::Temp);
  con_->value(r,y,tol);
}

template<typename Real>
void DykstraProjection<Real>::project_bnd(Vector<Real> &x, const Vector<Real> &y) const {
  x.set(y);
  bnd_->project(x);
}

template<typename Real>
void DykstraProjection<Real>::project_con(Vector<Real> &x, const Vector<Real> &y) const {
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
void DykstraProjection<Real>::project_Dykstra(Vector<Real> &x, std::ostream &stream) const {
  const Real one(1), xnorm(x.norm()), ctol(std::min(atol_,rtol_*xnorm));
  Real norm1(0), norm2(0), rnorm(0);
  p_->zero(); q_->zero();
  std::ios_base::fmtflags streamFlags(stream.flags());
  if (verbosity_ > 2) {
    stream << std::scientific << std::setprecision(6);
    stream << std::endl;
    stream << " Polyhedral Projection using Dykstra's Algorithm" << std::endl;
    stream << "  ";
    stream << std::setw(6)  << std::left << "iter";
    stream << std::setw(15) << std::left << "con norm";
    stream << std::setw(15) << std::left << "bnd norm";
    stream << std::setw(15) << std::left << "error";
    stream << std::setw(15) << std::left << "tol";
    stream << std::endl;
  }
  for (int cnt=0; cnt < maxit_; ++cnt) {
    // Constraint projection
    tmp_->set(x);   tmp_->plus(*p_);
    project_con(*y_,*tmp_);
    p_->set(*tmp_); p_->axpy(-one,*y_);
    // compute error between pnew and pold
    tmp_->set(x);   tmp_->axpy(-one,*y_);
    norm1 = tmp_->norm();
    // Bounds projection
    tmp_->set(*y_); tmp_->plus(*q_);
    project_bnd(x,*tmp_);
    q_->set(*tmp_); q_->axpy(-one,x);
    // compute error between qnew and qold
    tmp_->set(x);   tmp_->axpy(-one,*y_);
    norm2 = tmp_->norm();
    // stopping condition based on Birgin/Raydan paper
    // Robust Stopping Criteria for Dykstra's Algorithm
    // SISC Vol. 26, No. 4, 2005
    rnorm = std::sqrt(norm1*norm1 + norm2*norm2);
    if (verbosity_ > 2) {
      stream << "  ";
      stream << std::setw(6)  << std::left << cnt;
      stream << std::setw(15) << std::left << norm1;
      stream << std::setw(15) << std::left << norm2;
      stream << std::setw(15) << std::left << rnorm;
      stream << std::setw(15) << std::left << ctol;
      stream << std::endl;
    }
    if (rnorm <= ctol) break;
  }
  if (verbosity_ > 2) {
    stream << std::endl;
  }
  if (rnorm > ctol) {
    //throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : Projection failed!");
    stream << ">>> ROL::PolyhedralProjection::project : Projection may be inaccurate!  rnorm = ";
    stream << rnorm << "  rtol = " << ctol << std::endl;
  }
  stream.flags(streamFlags);
}

} // namespace ROL

#endif
