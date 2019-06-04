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

/*! \file  example_05.cpp
    \brief Shows how to compute the objective function,
           \f[
              \frac{\delta_t}{2} (u_{n-1}-1)^2 (u_n-1)^2 (z_n-1)^2.
           \f]
*/

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_DynamicObjective.hpp"

template<class Real>
class Objective_Nonlinear : public ROL::DynamicObjective<Real> {
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

  Objective_Nonlinear(ROL::ParameterList &pl) {
    const Real one(1);
    uint nt   = pl.get("Temporal Discretization",  100);
    Real  T   = pl.get("End Time",                 1.0);
    dt_       = T/(static_cast<Real>(nt)-one);
  }

  Real value( const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real half(0.5), one(1), two(2);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    return dt_ * half * std::pow(((*zp)[0]-one) * ((*uop)[0]-one) * ((*unp)[0]-one), two);
  }

  void gradient_uo( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>        gp = getVector(g);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*gp)[0] = dt_ * std::pow(((*zp)[0]-one) * ((*unp)[0]-one), two) * ((*uop)[0]-one);
  }

  void gradient_un( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
              const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>        gp = getVector(g);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*gp)[0] = dt_ * std::pow(((*zp)[0]-one) * ((*uop)[0]-one), two) * ((*unp)[0]-one);
  }

  void gradient_z( ROL::Vector<Real> &g,
             const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
             const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>        gp = getVector(g);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*gp)[0] = dt_ * std::pow(((*uop)[0]-one) * ((*unp)[0]-one), two) * ((*zp)[0]-one);
  }

  void hessVec_uo_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = dt_ * std::pow(((*zp)[0]-one) * ((*unp)[0]-one), two) * (*vp)[0];
  }

  void hessVec_uo_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*zp)[0]-one), two) * ((*uop)[0]-one) * ((*unp)[0]-one) * (*vp)[0];
  }

  void hessVec_uo_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*unp)[0]-one), two) * ((*uop)[0]-one) * ((*zp)[0]-one) * (*vp)[0];
  }

  void hessVec_un_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*zp)[0]-one), two) * ((*uop)[0]-one) * ((*unp)[0]-one) * (*vp)[0];
  }

  void hessVec_un_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = dt_ * std::pow(((*zp)[0]-one) * ((*uop)[0]-one), two) * (*vp)[0];
  }

  void hessVec_un_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*uop)[0]-one), two) * ((*unp)[0]-one) * ((*zp)[0]-one) * (*vp)[0];
  }

  void hessVec_z_uo( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*unp)[0]-one), two) * ((*uop)[0]-one) * ((*zp)[0]-one) * (*vp)[0];
  }

  void hessVec_z_un( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = two * dt_ * std::pow(((*uop)[0]-one), two) * ((*unp)[0]-one) * ((*zp)[0]-one) * (*vp)[0];
  }

  void hessVec_z_z( ROL::Vector<Real> &hv,   const ROL::Vector<Real>    &v,
               const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
               const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1), two(2);
    ROL::Ptr<std::vector<Real>>       hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    (*hvp)[0] = dt_ * std::pow(((*uop)[0]-one) * ((*unp)[0]-one), two) * (*vp)[0];
  }
};
