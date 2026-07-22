// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKCONSTRAINT_HPP
#define TANKCONSTRAINT_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ROL_Constraint_TimeSimOpt.hpp"
#include "ROL_StdVector.hpp"
#include "TankState.hpp"
#include <utility>

/** \class TankConstraint
    \brief Compute time-step for the coupled tank network
*/
namespace details {

using std::vector;

template<typename Real>
class TankConstraint : public ROL::Constraint_TimeSimOpt<Real> {

  using Vector    = ROL::Vector<Real>;
  using StdVector = ROL::StdVector<Real>;
 
  using StateVector = TankStateVector<Real>;
  using ControlVector = TankControlVector<Real>;
 
  using size_type = typename std::vector<Real>::size_type;

private:
  ROL::Ptr<TankState<Real>> tankState_;

  StateVector&   to_state  ( Vector& x ) { return static_cast<StateVector&>(x); }
  ControlVector& to_control( Vector& x ) { return static_cast<ControlVector&>(x); }
  const StateVector&   to_state  ( const Vector& x ) const { return static_cast<const StateVector&>(x);   }
  const ControlVector& to_control( const Vector& x ) const { return static_cast<const ControlVector&>(x); }

  size_type rows_, cols_;

  StateVector   zero_state_;
  ControlVector zero_ctrl_;

  bool has_applyJacobian_1_old_;
  bool has_applyJacobian_1_new_;
  bool has_applyJacobian_2_;

public:
  TankConstraint( const ROL::Ptr<TankState<Real>>& tankState, 
                  ROL::ParameterList& pl ) : tankState_(tankState), 
    rows_(pl.get("Number of Rows", 3)), cols_(pl.get("Number of Columns", 3)),
    zero_state_(rows_,cols_,"Zero State"), zero_ctrl_(rows_,cols_,"Zero Control"),
    has_applyJacobian_1_old_( pl.get("Has applyJacobian_1_old", false ) ),
    has_applyJacobian_1_new_( pl.get("Has applyJacobian_1_new", false ) ),
    has_applyJacobian_2_( pl.get("Has applyJacobian_2", false ) ) {
  }

  void print_tankstate_parameters( std::ostream& os ) { tankState_->print_members(os); }

  void value( Vector& c, const Vector& u_old, const Vector& u_new, 
              const Vector& z, Real& tol ) override {

    auto& c_state  = to_state(c); 
    auto& uo_state = to_state(u_old);
    auto& un_state = to_state(u_new);
    auto& z_ctrl   = to_control(z);
    tankState_->value( c_state, uo_state, un_state, z_ctrl ) ;
  }


  void solve( Vector& c, const Vector& u_old, Vector& u_new, 
              const Vector& z, Real& tol ) override {
  
    u_new.zero();  
    auto& c_state  = to_state(c);      
    auto& un_state = to_state(u_new);
    auto& uo_state = to_state(u_old);
    auto& z_ctrl   = to_control(z);
    tankState_->solve( c_state, un_state, uo_state, z_ctrl );
  }

  //----------------------------------------------------------------------------

  void applyJacobian_1_old( Vector& jv, const Vector& v_old,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, Real& tol ) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vo_state = to_state(v_old);
    tankState_->applyJacobian_1_old(jv_state,vo_state);
  }

  void applyAdjointJacobian_1_old( Vector &ajv_old, const Vector &dualv,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, Real& tol) override { 
    ajv_old.zero();
    auto& ajv_state = to_state(ajv_old);
    auto& dv_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_1_old( ajv_state, dv_state );
  }

  //----------------------------------------------------------------------------

  void applyJacobian_1_new( Vector& jv, const Vector& v_new,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, Real& tol ) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vn_state = to_state(v_new);
    tankState_->applyJacobian_1_new(jv_state,vn_state);
  }
 
  void applyAdjointJacobian_1_new( Vector& ajv_new, const Vector &dualv,
                                   const Vector &u_old, const Vector& u_new,
                                   const Vector &z, Real& tol) override {
    ajv_new.zero();
    auto& ajv_state = to_state(ajv_new);
    auto& dv_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_1_new( ajv_state, dv_state );
  }

  void applyInverseJacobian_1_new( Vector &ijv, const Vector &v_new,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, Real& tol ) override {
    ijv.zero();
    auto& ijv_state = to_state(ijv);      
    auto& vn_state  = to_state(v_new);
    tankState_->applyInverseJacobian_1_new(ijv_state,vn_state);
 }

  void applyInverseAdjointJacobian_1_new( Vector& iajv, const Vector& v_new,
                                          const Vector& u_old, const Vector& u_new,
                                          const Vector& z, Real& tol) override {
    iajv.zero();
    auto& iajv_state = to_state(iajv);      
    auto& vn_state  = to_state(v_new);
    tankState_->applyInverseAdjointJacobian_1_new(iajv_state,vn_state);
  }

  //----------------------------------------------------------------------------

  void applyJacobian_2( Vector &jv, const Vector &v,
                        const Vector &u_old, const Vector &u_new,
                        const Vector &z, Real &tol ) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& v_ctrl   = to_control(v);
    tankState_->applyJacobian_2( jv_state, v_ctrl );
}

  void applyAdjointJacobian_2_time( Vector& ajv, const Vector &dualv,
                                    const Vector &u_old, const Vector& u_new,
                                    const Vector &z, Real &tol) override {
    ajv.zero();
    auto& ajv_ctrl = to_control(ajv);
    auto& v_state  = to_state(dualv);
    tankState_->applyAdjointJacobian_2( ajv_ctrl, v_state );
  }


  //-----------------------------------------------------------------------------

  void applyAdjointHessian_11_old( Vector& ahwv_old, const Vector& w,
                                   const Vector& v_new, const Vector& u_old, 
                                   const Vector& u_new, const Vector& z,
                                   Real &tol) override { ahwv_old.zero(); }

  void applyAdjointHessian_11_new( Vector& ahwv_new, const Vector &w,
                                   const Vector &v_new, const Vector &u_old, 
                                   const Vector &u_new, const Vector &z,
                                   Real &tol) override { ahwv_new.zero(); }


}; // class TankConstraint

} // namespace details

using details::TankConstraint;

#endif // TANKCONSTRAINT_HPP

