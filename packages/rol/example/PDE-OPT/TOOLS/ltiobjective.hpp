// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDE_LTIOBJECTIVE_HPP
#define PDE_LTIOBJECTIVE_HPP

#include "ROL_DynamicObjective.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_LinearCombinationObjective_SimOpt.hpp"
#include "integralobjective.hpp"
#include "qoi.hpp"
#include "assembler.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template <class Real>
class LTI_Objective : public ROL::DynamicObjective<Real> {
private:
  const ROL::Ptr<ROL::Objective_SimOpt<Real>> obj_, objT_;
  Real theta_, T_;
  mutable ROL::Ptr<ROL::Vector<Real>> zdual_, udual_;
  mutable bool isInit_;

  void initialize(const ROL::Vector<Real> &z, const ROL::Vector<Real> &u) const {
    if (!isInit_) {
      zdual_ = z.dual().clone();
      udual_ = u.dual().clone();
      isInit_ = true;
    }
  }

  bool integratedObjective(void) const {
    return obj_ != ROL::nullPtr;
  }

  bool finalTimeObjective(const Real T) const {
    return objT_ != ROL::nullPtr && T_ == T;
  }

public:
  LTI_Objective(ROL::ParameterList                    &parlist,
          const ROL::Ptr<ROL::Objective_SimOpt<Real>> &obj,
          const bool finalTime = false)
    : obj_ (!finalTime ? obj : ROL::nullPtr),
      objT_( finalTime ? obj : ROL::nullPtr),
      isInit_(false) {
    theta_ = parlist.sublist("Time Discretization").get("Theta",    1.0);
    T_     = parlist.sublist("Time Discretization").get("End Time", 1.0);
  }

  LTI_Objective(ROL::ParameterList                    &parlist,
          const ROL::Ptr<ROL::Objective_SimOpt<Real>> &obj  = ROL::nullPtr, // Integrated objective
          const ROL::Ptr<ROL::Objective_SimOpt<Real>> &objT = ROL::nullPtr) // Final time objective
    : obj_(obj), objT_(objT), isInit_(false) {
    theta_ = parlist.sublist("Time Discretization").get("Theta",    1.0);
    T_     = parlist.sublist("Time Discretization").get("End Time", 1.0);
  }

  Real value( const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    // Integrated Objective
    Real valo(0), valn(0);
    if (integratedObjective()) {
      obj_->update(uo,z);
      valo = obj_->value(uo,z,tol);
      obj_->update(un,z);
      valn = obj_->value(un,z,tol);
    }
    // Final Time Objective
    Real valT(0);
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      valT = objT_->value(un,z,tol);
    }
    return dt*((one-theta_)*valo + theta_*valn) + valT;
  }

  void gradient_uo( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    if (integratedObjective()) {
      const Real one(1);
      Real timeOld = ts.t[0], timeNew = ts.t[1];
      Real dt = timeNew - timeOld;
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      obj_->update(uo,z);
      obj_->gradient_1(g,uo,z,tol);
      g.scale(dt*(one-theta_));
    }
    else {
      g.zero();
    }
  }

  void gradient_un( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    g.zero();
    if (integratedObjective()) {
      obj_->update(un,z);
      obj_->gradient_1(*udual_,un,z,tol);
      g.axpy(dt*theta_, *udual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->gradient_1(*udual_,un,z,tol);
      g.plus(*udual_);
    }
  }

  void gradient_z( ROL::Vector<Real> &g,
             const ROL::Vector<Real> &uo,
             const ROL::Vector<Real> &un,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    g.zero();
    if (integratedObjective()) {
      obj_->update(uo,z);
      obj_->gradient_2(*zdual_,uo,z,tol);
      g.axpy(dt*(one-theta_), *zdual_);
      obj_->update(un,z);
      obj_->gradient_2(*zdual_,un,z,tol);
      g.axpy(dt*theta_,*zdual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->gradient_2(*zdual_,un,z,tol);
      g.plus(*zdual_);
    }
  }

  void hessVec_uo_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    if (integratedObjective()) {
      const Real one(1);
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      Real timeOld = ts.t[0], timeNew = ts.t[1];
      Real dt = timeNew - timeOld;
      obj_->update(uo,z);
      obj_->hessVec_11(*udual_,v,uo,z,tol);
      hv.axpy(dt*(one-theta_), *udual_);
    }
    else {
      hv.zero();
    }
  }

  void hessVec_uo_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    hv.zero();
  }

  void hessVec_uo_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    if (integratedObjective()) {
      const Real one(1);
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      Real timeOld = ts.t[0], timeNew = ts.t[1];
      Real dt = timeNew - timeOld;
      obj_->update(uo,z);
      obj_->hessVec_12(hv,v,uo,z,tol);
      hv.scale(dt*(one-theta_));
    }
    else {
      hv.zero();
    }
  }

  void hessVec_un_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    hv.zero();
  }

  void hessVec_un_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    if (integratedObjective()) {
      obj_->update(un,z);
      obj_->hessVec_11(*udual_,v,un,z,tol);
      hv.axpy(dt*theta_, *udual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->hessVec_11(*udual_,v,un,z,tol);
      hv.plus(*udual_);
    }
  }

  void hessVec_un_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    if (integratedObjective()) {
      obj_->update(un,z);
      obj_->hessVec_12(*udual_,v,un,z,tol);
      hv.axpy(dt*theta_, *udual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->hessVec_12(*udual_,v,un,z,tol);
      hv.plus(*udual_);
    }
  }

  void hessVec_z_uo( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    if (integratedObjective()) {
      const Real one(1);
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      Real timeOld = ts.t[0], timeNew = ts.t[1];
      Real dt = timeNew - timeOld;
      obj_->update(uo,z);
      obj_->hessVec_21(hv,v,uo,z,tol);
      hv.scale(dt*(one-theta_));
    }
    else {
      hv.zero();
    }
  }

  void hessVec_z_un( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    if (integratedObjective()) {
      obj_->update(un,z);
      obj_->hessVec_21(*zdual_,v,un,z,tol);
      hv.axpy(dt*theta_,*zdual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->hessVec_21(*zdual_,v,un,z,tol);
      hv.plus(*zdual_);
    }
  }

  void hessVec_z_z( ROL::Vector<Real> &hv,
              const ROL::Vector<Real> &v,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    if (integratedObjective()) {
      obj_->update(uo,z);
      obj_->hessVec_22(*zdual_,v,uo,z,tol);
      hv.axpy(dt*(one-theta_),*zdual_);
      obj_->update(un,z);
      obj_->hessVec_22(*zdual_,v,un,z,tol);
      hv.axpy(dt*theta_,*zdual_);
    }
    if (finalTimeObjective(timeNew)) {
      objT_->update(un,z);
      objT_->hessVec_22(*zdual_,v,un,z,tol);
      hv.plus(*zdual_);
    }
  }
}; // class LTI_Objective

#endif
