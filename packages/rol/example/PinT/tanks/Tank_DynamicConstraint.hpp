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


#pragma once
#ifndef TANK_DYNAMICCONSTRAINT_HPP
#define TANK_DYNAMICCONSTRAINT_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "TankState.hpp"
#include <utility>


/** \class Tank_DynamicConstraint based on the new DynamicConstraint interface
    \brief Compute time-step for the coupled tank network
*/

namespace details {

using namespace ROL;

template<typename Real>
class Tank_DynamicConstraint : public DynamicConstraint<Real> {

  using State   = TankStateVector<Real>;
  using Control = TankControlVector<Real>;
  using size_type = typename std::vector<Real>::size_type;

private:

  ROL::Ptr<TankState<Real>> tankState_;

  State&   to_state  ( Vector<Real>& x ) { return static_cast<State&>(x); }
  Control& to_control( Vector<Real>& x ) { return static_cast<Control&>(x); }
  const State&   to_state  ( const Vector<Real>& x ) const { return static_cast<const State&>(x);   }
  const Control& to_control( const Vector<Real>& x ) const { return static_cast<const Control&>(x); }

  State   zero_state_;
  Control zero_ctrl_;

public: 

  Tank_DynamicConstraint( const ROL::Ptr<TankState<Real>>& tankState, 
                  Teuchos::ParameterList& pl ) : tankState_(tankState), 
    rows_(pl.get("Number of Rows", 3)), cols_(pl.get("Number of Columns", 3)),
    zero_state_(rows_,cols_,"Zero State"), zero_ctrl_(rows_,cols_,"Zero Control") {}

  void value( Vector& c, const Vector& u_old, const Vector& u_new, 
              const Vector& z, TimeStamp<Real>& ts ) override {

    auto& c_state  = to_state(c); 
    auto& uo_state = to_state(u_old);
    auto& un_state = to_state(u_new);
    auto& z_ctrl   = to_control(z);
    c.zero();
    tankState_->value( c_state, uo_state, un_state, z_ctrl ) ;
  }

  void solve( Vector& c, const Vector& u_old, Vector& u_new, 
              const Vector& z, TimeStamp<Real>& ts ) override {
  
    u_new.zero();  
    auto& c_state  = to_state(c);      
    auto& un_state = to_state(u_new);
    auto& uo_state = to_state(u_old);
    auto& z_ctrl   = to_control(z);
    tankState_->solve( c_state, un_state, uo_state, z_ctrl );
  }

  void applyJacobian_uo( Vector& jv, const Vector& v_old,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z,  TimeStamp<Real>& ts) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vo_state = to_state(v_old);
    tankState_->applyJacobian_1_old(jv_state,vo_state);
  }

  void applyAdjointJacobian_uo( Vector &ajv_old, const Vector &dualv,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, TimeStamp<Real>& ts) override { 
    ajv_old.zero();
    auto& ajv_state = to_state(ajv_old);
    auto& dv_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_1_old( ajv_state, dv_state );
  }

  //----------------------------------------------------------------------------

  void applyJacobian_un( Vector& jv, const Vector& v_new,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, TimeStamp<Real>& ts ) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vn_state = to_state(v_new);
    tankState_->applyJacobian_1_new(jv_state,vn_state);
  }
 
  void applyAdjointJacobian_un( Vector& ajv_new, const Vector &dualv,
                                   const Vector &u_old, const Vector& u_new,
                                   const Vector &z, TimeStamp<Real>& ts) override {
    ajv_new.zero();
    auto& ajv_state = to_state(ajv_new);
    auto& dv_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_1_new( ajv_state, dv_state );
  }

  void applyInverseJacobian_un( Vector &ijv, const Vector &v_new,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, TimeStamp<Real>& ts ) override {
    ijv.zero();
    auto& ijv_state = to_state(ijv);      
    auto& vn_state  = to_state(v_new);
    tankState_->applyInverseJacobian_1_new(ijv_state,vn_state);
 }

  void applyInverseAdjointJacobian_un( Vector& iajv, const Vector& v_new,
                                          const Vector& u_old, const Vector& u_new,
                                          const Vector& z, TimeStamp<Real>& ts) override {
    iajv.zero();
    auto& iajv_state = to_state(iajv);      
    auto& vn_state  = to_state(v_new);
    tankState_->applyInverseAdjointJacobian_1_new(iajv_state,vn_state);
  }


  //----------------------------------------------------------------------------

  void applyJacobian_z( Vector &jv, const Vector &v,
                        const Vector &u_old, const Vector &u_new,
                        const Vector &z, Real &tol TimeStamp<Real>& ts) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& v_ctrl   = to_control(v);
    tankState_->applyJacobian_2( jv_state, v_ctrl );
}

  void applyAdjointJacobian_z( Vector& ajv, const Vector &dualv,
                                    const Vector &u_old, const Vector& u_new,
                                    const Vector &z, TimeStamp<Real>& ts ) override {
    ajv.zero();
    auto& ajv_ctrl = to_control(ajv);
    auto& v_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_2( ajv_ctrl, v_state );
  }

  //-----------------------------------------------------------------------------

  void applyAdjointHessian_un_un( Vector& ahwv_old, const Vector& w,
                                   const Vector& v_new, const Vector& u_old, 
                                   const Vector& u_new, const Vector& z,
                                   TimeStamp<Real>& ts ) override { ahwv_old.zero(); }


};

using details::Tank_DynamicConstraint;

} // namespace details


#endif // TANK_DYNAMICCONSTRAINT_HPP

