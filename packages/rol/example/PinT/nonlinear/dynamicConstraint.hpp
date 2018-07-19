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

/*! \file  dynamicConstraint.cpp
    \brief Shows how to solve arbitrary nonlinear constraint
           \f[
              e^{u_{n} - u_{n-1} - \delta_t z_n} = 1.
           \f]
*/

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_DynamicConstraint.hpp"

template<class Real>
class Constraint_Nonlinear : public ROL::DynamicConstraint<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Real dt_;

  ROL::Ptr<const std::vector<Real>> getVector( const ROL::Vector<Real>& x ) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real>> getVector( ROL::Vector<Real>& x ) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

public:

  Constraint_Nonlinear(ROL::ParameterList &pl) {
    const Real one(1);
    uint nt   = pl.get("Temporal Discretization", 100);
    Real  T   = pl.get("End Time",                1.0);
    dt_       = T/(static_cast<Real>(nt)-one);
  }

  void value(ROL::Vector<Real>    &c,    const ROL::Vector<Real> &uold,
       const ROL::Vector<Real>    &unew, const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) const {
    const Real one(1);
    ROL::Ptr<std::vector<Real>>        cp = getVector(c);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*cp)[0] = std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) - one;
  }

  void solve(ROL::Vector<Real>    &c,    const ROL::Vector<Real> &uold,
             ROL::Vector<Real>    &unew, const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<std::vector<Real>>       unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*unp)[0] = (*uop)[0] + dt_*(*zp)[0];
    value(c, uold, unew, z, ts);
  }

  void applyJacobian_uo(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                  const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*jvp)[0] = -std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*vp)[0];
  }

  void applyJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                  const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*jvp)[0] = std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*vp)[0];
  }

  void applyJacobian_z(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                 const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                 const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*jvp)[0] = -dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*vp)[0];
  }

  void applyAdjointJacobian_uo(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_uo(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_z(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                        const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                        const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_z(jv, v, uold, unew, z, ts);
  }

  void applyInverseJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*jvp)[0] = (*vp)[0] / std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]);
  }

  void applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyInverseJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointHessian_un_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_un_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = -std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_un_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = -dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_uo_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = -std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_uo_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_uo_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_z_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = -dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_z_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_z_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    (*ahwvp)[0] = dt_ * dt_ * std::exp((*unp)[0] - (*uop)[0] - dt_*(*zp)[0]) * (*wp)[0] * (*vp)[0];
  }
};
