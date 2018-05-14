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
#ifndef TANKS_DYNAMICCONSTRAINT_HPP
#define TANKS_DYNAMICCONSTRAINT_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ROL_DynamicConstraint.hpp"

#include "Tanks_StateVector.hpp"
#include "Tanks_ControlVector.hpp"
#include "LowerBandedMatrix.hpp"

#include <utility>


/** \class Tanks_DynamicConstraint based on the new DynamicConstraint interface
    \brief Compute time-step for the coupled tank network
*/

namespace Tanks {

template<typename Real>
class DynamicConstraint : public ROL::DynamicConstraint<Real> {

  using State   = StateVector<Real>;
  using Control = ControlVector<Real>;

  using size_type = typename State::size_type;

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

  State   zero_state_;
  Control zero_ctrl_;

public: 

  DynamicConstraint( ROL::ParameterList& pl ) :
  // ----------- Begin Initializer List ----------------//
  rows_   ( pl.get( "Number of Rows",        3      ) ),
  cols_   ( pl.get( "Number of Columns",     3      ) ),
  Cv_     ( pl.get( "Valve Constant",        1.0e-2 ) ),
  rho_    ( pl.get( "Density of Fluid",      1.0e3  ) ),
  h0_     ( pl.get( "Initial Fluid Level",   2.0    ) ),
  H_      ( pl.get( "Height of Tank",        10.0   ) ),
  A_      ( pl.get( "Cross-sectional Area",  10.0   ) ),
  g_      ( pl.get( "Gravity Constant",      9.8    ) ),
  T_      ( pl.get( "Total Time",            20.0   ) ),
  theta_  ( pl.get( "Theta",                 0.5    ) ),
  Nt_     ( pl.get( "Number of Time Steps",  100    ) ),
  //----------------------------------------------------//
  dt_( T_/Nt_ ),
  Ntanks_(rows_*cols_),
  p_(rows_,cols_), scratch_(rows_,cols_),
  kappa_( Cv_*rho_*g_ ), 
  beta_( dt_/A_ ),
  alphaL_( 0.5*kappa_*(1-theta_)*beta_ ),
  alphaR_( 0.5*kappa_*theta_*beta_ ),
  h_(0), Qout_(2*Ntanks_), Qin_(Ntanks_), z_(0),
  zero_state_(rows_,cols_,"Zero State"), 
  zero_ctrl_(rows_,cols_,"Zero Control") {
  // ------------- End Initializer List ----------------//
  
    auto ptrows = ROL::getArrayFromStringParameter<int>( pl, "Pass-Through Rows"    );
    auto ptcols = ROL::getArrayFromStringParameter<int>( pl, "Pass-Through Columns" );
  
    vector<size_type> band_index{0, 1, cols_};
    vector<size_type> shift_index{1, cols_};
  
    p_.setScalar(1.0);
  
    for( size_type i=0; i<static_cast<size_type>(ptrows.size()); ++i ) {
      p_( ptrows.at(i), ptcols.at(i) ) = 0.0;
    }
  
    L_ = make_shared<TankLevelMatrix<Real>>( rows_, cols_,  alphaL_, *(p_.getVector()) );
    R_ = make_shared<TankLevelMatrix<Real>>( rows_, cols_, -alphaR_, *(p_.getVector()) );
    S_ = make_shared<SplitterMatrix<Real>>( rows_, cols_ );  
  }

  static ROL::Ptr<DynamicConstraint> create( ROL::ParameterList& pl ) {
    return ROL::makePtr<DynamicConstraint>( pl );
  }

  void value( Vector& c, const Vector& u_old, const Vector& u_new, 
              const Vector& z, TimeStamp<Real>& ts ) override {

    auto& c_state  = to_state(c); 
    auto& uo_state = to_state(u_old);
    auto& un_state = to_state(u_new);
    auto& z_ctrl   = to_control(z);
    c.zero();

    c_state.axpy( -beta_, z_ctrl, h_ );
    c_state.hadamard( p_, h_ );
    c_state.set( un_state, Qin_, Qin_ );
    c_state.axpy( -1.0, z_ctrl, Qin_ );
    c_state.set( un_state, Qout_, Qout_ );
    c_state.axpy( -kappa_, un_state, Qout_, h_ );
  
    L_->apply( c_state, un_state,  1.0, h_, h_ ); 
    R_->apply( c_state, uo_state, -1.0, h_, h_ );
    S_->apply( c_state, un_state, -1.0, Qin_, Qout_ );  
  }


  void solve( Vector& c, const Vector& u_old, Vector& u_new, 
              const Vector& z, TimeStamp<Real>& ts ) override {
  
    u_new.zero();  
    auto& c_state  = to_state(c);      
    auto& un_state = to_state(u_new);
    auto& uo_state = to_state(u_old);
    auto& z_ctrl   = to_control(z);

    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        u_new.h(i,j) =  u_old.h(i,j) + p_(i,j)*(beta_*z(i,j)-2.0*alphaR_*u_old.h(i,j));
  
        if( i>0 ) u_new.h(i,j) += p_(i,j)*(alphaL_*u_new.h(i-1,j) + alphaR_*u_old.h(i-1,j));
        if( j>0 ) u_new.h(i,j) += p_(i,j)*(alphaL_*u_new.h(i,j-1) + alphaR_*u_old.h(i,j-1));
  
        u_new.h(i,j) /= (1.0+2.0*alphaL_*p_(i,j));
        u_new.Qout(i,j) = kappa_*u_new.h(i,j);   
        u_new.Qin(i,j)  = 1.0*z(i,j);
       }
    }
    S_->apply( u_new, u_new, 1.0, Ntanks_, 2*Ntanks_ );

  }

  void applyJacobian_uo( Vector& jv, const Vector& v_old,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z,  TimeStamp<Real>& ts) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vo_state = to_state(v_old);
    R_->apply( jv_state, vo_state, -1.0, h_, h_);      
  }

  void applyAdjointJacobian_uo( Vector &ajv_old, const Vector &dualv,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, TimeStamp<Real>& ts) override { 
    ajv_old.zero();
    auto& ajv_state = to_state(ajv_old);
    auto& dv_state  = to_state(dualv);
    R_->applyTranspose( ajv_state, dv_state, -1.0, h_, h_);  
  }

  //----------------------------------------------------------------------------

  void applyJacobian_un( Vector& jv, const Vector& v_new,
                            const Vector& u_old, const Vector& u_new,
                            const Vector& z, TimeStamp<Real>& ts ) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vn_state = to_state(v_new);
    L_->apply( jv_state, vn_state, 1.0, h_, h_ );  
    jv.set( vn_state, Qout_, Qout_ );
    jv.axpy( -kappa_, vn_state, Qout_, h_ );
    jv.set( vn_state, Qin_, Qin_ );
    S_->apply( jv, vn_state, -1.0, Qin_, Qout_ );   
  }
 
  void applyAdjointJacobian_un( Vector& ajv_new, const Vector &dualv,
                                   const Vector &u_old, const Vector& u_new,
                                   const Vector &z, TimeStamp<Real>& ts) override {
    ajv_new.zero();
    auto& ajv_state = to_state(ajv_new);
    auto& dv_state  = to_state(dualv);
    L_->applyTranspose( jv_state, vn_state 1.0, h_, h_ );  
    jv.axpy(-kappa_, vn_state, h_, Qout_ );
    jv.set( vn_state, Qout_, Qout_ );
    jv.set( vn_state, Qin_, Qin_ );
    S_->applyTranspose( jv_state, vn_state, -1.0, Qout_, Qin_ );
  }

  void applyInverseJacobian_un( Vector &ijv, const Vector &v_new,
                                   const Vector &u_old, const Vector &u_new,
                                   const Vector &z, TimeStamp<Real>& ts ) override {
    ijv.zero();
    auto& ijv_state = to_state(ijv);      
    auto& vn_state  = to_state(v_new);
    L_->solve( ijv_state, vn_state, 1.0, h_, h_ );  
    ijv.set( v_new, Qout_, Qout_ );
    ijv.axpy( kappa_, ijv_state, Qout_, h_ );
    ijv.set( vn_state, Qin_, Qin_ );
    S_->apply( ijv_state, ijv_state, 1.0, Qin_, Qout_ );
 }

  void applyInverseAdjointJacobian_un( Vector& iajv, const Vector& v_new,
                                          const Vector& u_old, const Vector& u_new,
                                          const Vector& z, TimeStamp<Real>& ts) override {
    iajv.zero();
    auto& iajv_state = to_state(iajv);      
    auto& vn_state  = to_state(v_new);
    iajv_state.set( vn_state, Qin_, Qin_ );
    iajv_state.set( vn_state, Qout_, Qout_ );
    S_->applyTranspose( iajv_state, vn_state, 1.0, Qout_, Qin_ );
    scratch_.zero();
    scratch_.set( vn_state, h_, h_ );
    scratch_.axpy( kappa_, iajv_state, h_, Qout_ );
    L_->solveTranspose( iajv_state, scratch_, 1.0, h_, h_ );  

  }


  //----------------------------------------------------------------------------

  void applyJacobian_z( Vector &jv, const Vector &v,
                        const Vector &u_old, const Vector &u_new,
                        const Vector &z, Real &tol TimeStamp<Real>& ts) override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& v_ctrl   = to_control(v);
    jv_state.axpy(-beta_,v_ctrl, h_);
    jv_state.hadamard( p_, h_ );
    jv_state.axpy(-1.0, v_ctrl, Qin_ );    
}

  void applyAdjointJacobian_z( Vector& ajv, const Vector &dualv,
                                    const Vector &u_old, const Vector& u_new,
                                    const Vector &z, TimeStamp<Real>& ts ) override {
    ajv.zero();
    auto& ajv_ctrl = to_control(ajv);
    auto& v_state  = to_state(dualv);
    for( size_type i=0; i<rows_; ++i ) 
      for( size_type j=0; j<cols_; ++j ) 
        ajv_ctrl(i,j) = -beta_*p_(i,j)*v_state.h(i,j)-v_state.Qin(i,j);
    }

};

} // namespace Tanks


#endif // TANKS_DYNAMICCONSTRAINT_HPP

