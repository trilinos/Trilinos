// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKSTATE_HPP
#define TANKSTATE_HPP

#include <memory>
#include "Teuchos_ParameterList.hpp"
#include "LowerBandedMatrix.hpp"
#include "TankVector.hpp"

namespace details {

using namespace std;

template<typename Real> 
class TankState {
public:
  using Vector        = vector<Real>;
  using Matrix        = LowerBandedMatrix<Real>;
  using StateVector   = TankStateVector<Real>;
  using ControlVector = TankControlVector<Real>;

  using size_type  = typename vector<Real>::size_type;
   
  TankState( ROL::ParameterList& pl );

  void solve( StateVector& c, StateVector& u_new, 
              const StateVector& u_old, const ControlVector& z ) const;

  void value( StateVector& c, const StateVector& u_old, 
              const StateVector& u_new, const ControlVector& z ) const;
 
  void applyJacobian_1_old( StateVector& jv, const StateVector& v_old ) const;
  void applyJacobian_1_new( StateVector& jv, const StateVector& v_old ) const;
  void applyInverseJacobian_1_new( StateVector& ijv, const StateVector& v_new ) const;
  void applyAdjointJacobian_1_old( StateVector& ajv, const StateVector& v_old ) const;
  void applyAdjointJacobian_1_new( StateVector& ajv, const StateVector& v_new ) const;
  void applyInverseAdjointJacobian_1_new( StateVector& iajv, const StateVector& v_new ) const;
  void applyJacobian_2( StateVector& jv, const ControlVector& v ) const;
  void applyAdjointJacobian_2( ControlVector& ajv, const StateVector& v ) const;

  void print_members( ostream& os ) const;

private:
  size_type rows_;             // Number of tank rows
  size_type cols_;             // Number of tank columns

  //---------- Physical Parameters ---------------------------------------------
  Real Cv_;                    // Valve constant 
  Real rho_;                   // Density of fluid
  Real h0_;                    // Initial fluid level
  Real H_;                     // Height of tanks
  Real A_;                     // Cross-sectional area of tanks
  Real g_;                     // Gravity constant
  
  //---------- Time Discretization ---------------------------------------------
  Real T_;                     // Total time
  Real theta_;                 // Implicit/Explicit splitting factor
  int  Nt_;                    // Number of Time steps
  Real dt_;                    // Time step size

  size_type    Ntanks_;        // Total number of tanks 

  ControlVector p_;            // Passthrough coefficients
  mutable StateVector scratch_;

  Real         kappa_;         // Cv*rho*g
  Real         beta_;          // dt/A
  Real         alphaL_;        // kappa*(theta-1)
  Real         alphaR_;        // kappa*theta

  shared_ptr<Matrix> L_, R_, S_;

  //--------- Subvector addressing ---------------------------------------------
  size_type  h_, Qout_, Qin_,  z_;

}; // class TankState

} // namespace details

using details::TankState;

#include "TankState_Impl.hpp"

#endif // TANKSTATE_HPP

