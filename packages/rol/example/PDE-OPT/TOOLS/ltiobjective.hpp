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
  const ROL::Ptr<ROL::Objective_SimOpt<Real>> obj_;
  Real theta_;
  mutable ROL::Ptr<ROL::Vector<Real>> zdual_;

public:
  LTI_Objective(const ROL::Ptr<ROL::Objective_SimOpt<Real>> &obj,
                const ROL::Vector<Real> &z,
                ROL::ParameterList &parlist) :obj_(obj) {
    // Get time discretization parameters
    theta_  = parlist.sublist("Time Discretization").get("Theta",                1.0);
    zdual_ = z.dual().clone();
  }

  Real value( const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    obj_->update(uo,z);
    Real valo = obj_->value(uo,z,tol);
    obj_->update(un,z);
    Real valn = obj_->value(un,z,tol);
    return dt*((one-theta_)*valo + theta_*valn);
  }

  void gradient_uo( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    obj_->update(uo,z);
    obj_->gradient_1(g,uo,z,tol);
    g.scale(dt*(one-theta_));
  }

  void gradient_un( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(un,z);
    obj_->gradient_1(g,un,z,tol);
    g.scale(dt*theta_);
  }

  void gradient_z( ROL::Vector<Real> &g,
             const ROL::Vector<Real> &uo,
             const ROL::Vector<Real> &un,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z);
    obj_->gradient_2(g,uo,z,tol);
    g.scale(dt*(one-theta_));
    obj_->update(un,z);
    obj_->gradient_2(*zdual_,un,z,tol);
    g.axpy(dt*theta_,*zdual_);
  }

  void hessVec_uo_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z);
    obj_->hessVec_11(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
  }

  void hessVec_uo_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_uo_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z);
    obj_->hessVec_12(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
  }

  void hessVec_un_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_un_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(un,z);
    obj_->hessVec_11(hv,v,un,z,tol);
    hv.scale(dt*theta_);
  }

  void hessVec_un_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(un,z);
    obj_->hessVec_12(hv,v,un,z,tol);
    hv.scale(dt*theta_);
  }

  void hessVec_z_uo( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z);
    obj_->hessVec_21(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
  }

  void hessVec_z_un( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(un,z);
    obj_->hessVec_21(hv,v,un,z,tol);
    hv.scale(dt*theta_);
  }

  void hessVec_z_z( ROL::Vector<Real> &hv,
              const ROL::Vector<Real> &v,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z);
    obj_->hessVec_22(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
    obj_->update(un,z);
    obj_->hessVec_22(*zdual_,v,un,z,tol);
    hv.axpy(dt*theta_,*zdual_);
  }
}; // class LTI_Objective

#endif
