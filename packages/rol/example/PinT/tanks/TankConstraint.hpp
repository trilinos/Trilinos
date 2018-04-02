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

private:
  TankState<Real> tankState_;

  vector<Real>& getVector( Vector& x ) const {
    return *( dynamic_cast<StdVector&>(x).getVector() );
  }

  const vector<Real>& getVector( const Vector& x ) const {
    return *( dynamic_cast<const StdVector&>(x).getVector() );
  }

public:
  TankConstraint( Teuchos::ParameterList& pl ) :
    tankState_(pl) {}

  void value( Vector& c, const Vector& u_old, const Vector& u_new, 
              const Vector& z, Real& tol ) override {

    c.zero();

    StateVector& c_state  = dynamic_cast<StateVector&>(c);      

    const StateVector& uo_state = dynamic_cast<const StateVector&>(u_old);
    const StateVector& un_state = dynamic_cast<const StateVector&>(u_new);
    const ControlVector& z_ctrl = dynamic_cast<const ControlVector&>(z);

    tankState_.value( c_state, uo_state, un_state, z_ctrl ) ;
  }

  void solve( Vector& c, const Vector& u_old, Vector& u_new, 
              const Vector& z, Real& tol ) override {

    u_new.zero();

    StateVector& c_state  = dynamic_cast<StateVector&>(c);      
    StateVector& un_state = dynamic_cast<StateVector&>(u_new);

    const StateVector& uo_state = dynamic_cast<const StateVector&>(u_old);
    const ControlVector& z_ctrl = dynamic_cast<const ControlVector&>(z);
 
    tankState_.solve( c_state, un_state, uo_state, z_ctrl );
  }

  //----------------------------------------------------------------------------

  void applyJacobian_1_old( Vector& jv, const Vector& v_old,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, Real& tol ) override {
//     tankState_.applyJacobian_1_old( getVector(jv), getVector(v_old) );
  }

  void applyAdjointJacobian_1_old( Vector &ajv_old, const Vector &dualv,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, Real& tol) override { 
  }

  //----------------------------------------------------------------------------

  void applyJacobian_1_new( Vector& jv, const Vector& v_new,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, Real& tol ) override {
 //    tankState_.applyJacobian_1_new( getVector(jv), getVector(v_new) );
  }
 
  void applyAdjointJacobian_1_new( Vector& ajv_new, const Vector &dualv,
                                   const Vector &u_old, const Vector& u_new,
                                   const Vector &z, Real& tol) override {

  }

  void applyInverseJacobian_1_new( Vector &ijv, const Vector &v_new,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, Real& tol ) override {
  }

  void applyInverseAdjointJacobian_1_new( Vector& iajv, const Vector& v_new,
                                          const Vector& u_old, const Vector& u_new,
                                          const Vector& z, Real& tol) override {

  }

  //----------------------------------------------------------------------------

  void applyJacobian_2( Vector &jv, const Vector &v_new,
                        const Vector &u_old, const Vector &u_new,
                        const Vector &z, Real &tol ) override {
//    tankState_.applyJacobian_2( getVector(jv), getVector(v_new) );
  }

  void applyAdjointJacobian_2_time( Vector& ajv, const Vector &dualv,
                                    const Vector &u_old, const Vector& u_new,
                                    const Vector &z, Real &tol) override {
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

