// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_05.cpp
    \brief Shows how to compute the objective function,
           \f[
              J(u,z) = \frac{1}{2} \int_0^T \int_0^1 (u(t,x)-\bar{u}(x))^2\,\mathrm{d}x\,\mathrm{d}t
                       \frac{\alpha}{2} \int_0^T z(t)^2\,\mathrm{d}t.
           \f]
*/

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_DynamicObjective.hpp"

template<class Real>
class Objective_ParabolicControl : public ROL::DynamicObjective<Real> {
  
  typedef typename std::vector<Real>::size_type uint;

private:
  Real alpha_;
  Real theta_;
  uint nx_;
  Real dx_;
  Real dt_;
  int  type_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) const {
    const Real two(2), four(4), six(6);
    Real ip(0);
    Real c = ((x.size()==nx_) ? four : two);
    for (uint i=0; i<x.size(); i++) {
      if ( i == 0 ) {
        ip += dx_/six*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == x.size()-1 ) {
        ip += dx_/six*(x[i-1] + c*x[i])*y[i];
      }
      else {
        ip += dx_/six*(x[i-1] + four*x[i] + x[i+1])*y[i];
      }
    }
    return ip;
  }

  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) const {
    const Real two(2), four(4), six(6);
    Mu.resize(u.size());
    Real c = ((u.size()==nx_) ? four : two);
    for (uint i=0; i<u.size(); i++) {
      if ( i == 0 ) {
        Mu[i] = dx_/six*(c*u[i] + u[i+1]);
      }
      else if ( i == u.size()-1 ) {
        Mu[i] = dx_/six*(u[i-1] + c*u[i]);
      }
      else {
        Mu[i] = dx_/six*(u[i-1] + four*u[i] + u[i+1]);
      }
    }
  }

  Real evaluate_target(Real x) const {
    const Real zero(0), half(0.5), eight(8), pi(M_PI);
    Real val(0);
    switch (type_) {
      case 1:  val = (x<half ? half : zero);                 break;
      case 2:  val = half;                                   break;
      case 3:  val = half*std::abs(std::sin(eight*pi*x));    break;
      case 4:  val = half*std::exp(-half*(x-half)*(x-half)); break;
      default: val = half;                                   break;
    }
    return val;
  }

  ROL::Ptr<const std::vector<Real>> getVector( const ROL::Vector<Real>& x ) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real>> getVector( ROL::Vector<Real>& x ) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();  
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_ParabolicControl(ROL::ParameterList &pl) {
    const Real one(1);
    uint nt   = pl.get("Temporal Discretization",  100);
    Real  T   = pl.get("End Time",                 1.0);
    type_     = pl.get("Target Type",                1);
    alpha_    = pl.get("Control Penalty",         1e-4);
    nx_       = pl.get("Spatial Discretization",   128); 
    theta_    = pl.get("Integration Factor",       0.5);
    dx_       = one/(static_cast<Real>(nx_)-one);
    dt_       = T/(static_cast<Real>(nt)-one);
  }

  Real value( const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real half(0.5), one(1);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    // Compute Norm Squared of State
    Real target(0);
    std::vector<Real> uOld(nx_), uNew(nx_);
    for (uint n = 0; n < nx_; n++) {
      target = evaluate_target(static_cast<Real>(n)*dx_);
      uOld[n] = (*uop)[n] - target;
      uNew[n] = (*unp)[n] - target;
    }
    Real uoval = dot(uOld,uOld);
    Real unval = dot(uNew,uNew);
    // Add Norm Squared of Control
    Real zval = dot(*zp, *zp);
    return half * dt_ * (theta_*uoval + (one-theta_)*unval + alpha_ * zval);
  }

  void gradient_uo( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<std::vector<Real>>        gp = getVector(g);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    std::vector<Real> uDiff(nx_);
    Real target(0);
    for (uint n = 0; n < nx_; n++) {
      target = evaluate_target(static_cast<Real>(n)*dx_);
      uDiff[n] = (*uop)[n] - target;
    }
    apply_mass(*gp, uDiff);
    g.scale(theta_*dt_);
  }

  void gradient_un( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    ROL::Ptr<std::vector<Real>>        gp = getVector(g);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    std::vector<Real> uDiff(nx_);
    Real target(0);
    for (uint n = 0; n < nx_; n++) {
      target = evaluate_target(static_cast<Real>(n)*dx_);
      uDiff[n] = (*unp)[n] - target;
    }
    apply_mass(*gp, uDiff);
    g.scale((one-theta_)*dt_);
  }

  void gradient_z( ROL::Vector<Real> &g,
             const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
             const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<std::vector<Real>>       gp = getVector(g);
    ROL::Ptr<const std::vector<Real>> zp = getVector(z);
    apply_mass(*gp, *zp);
    g.scale(dt_*alpha_);
  }

  void hessVec_uo_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv); 
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    apply_mass(*hvp, *vp);
    hv.scale(theta_*dt_);
  }

  void hessVec_uo_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_uo_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_un_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_un_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv); 
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    apply_mass(*hvp, *vp);
    hv.scale((one-theta_)*dt_);
  }

  void hessVec_un_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<std::vector<Real>>      hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>> vp = getVector(v);
    apply_mass(*hvp, *vp);
    hv.scale(dt_*alpha_);
  }
};
