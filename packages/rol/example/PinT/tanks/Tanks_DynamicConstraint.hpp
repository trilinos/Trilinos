// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKS_DYNAMICCONSTRAINT_HPP
#define TANKS_DYNAMICCONSTRAINT_HPP

#include "ROL_ParameterList.hpp"
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

  using State     = StateVector<Real>;
  using Control   = ControlVector<Real>;
  using Matrix    = LowerBandedMatrix<Real>;
  using V         = ROL::Vector<Real>;
  using TS        = ROL::TimeStamp<Real>;
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

  Control p_;                  // Passthrough coefficients
  mutable State scratch_;

  Real         kappa_;         // Cv*rho*g
  Real         beta_;          // dt/A
  Real         alphaL_;        // kappa*(theta-1)
  Real         alphaR_;        // kappa*theta

  //--------- Subvector addressing ---------------------------------------------
  size_type  h_, Qout_, Qin_,  z_;

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

  void value( V& c, const V& uo, const V& un, 
              const V& z, const TS& ts ) const override {

    auto& c_state  = to_state(c); 
    auto& uo_state = to_state(uo);
    auto& un_state = to_state(un);
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

  void solve( V& c, const V& uo, V& un, 
              const V& z, const TS& ts ) override {
  
    un.zero();  
    auto& un_state = to_state(un);
    auto& uo_state = to_state(uo);
    auto& z_ctrl   = to_control(z);

    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        un_state.h(i,j) =  uo_state.h(i,j) + p_(i,j)*(beta_*z_ctrl(i,j)-2.0*alphaR_*uo_state.h(i,j));
  
        if( i>0 ) un_state.h(i,j) += p_(i,j)*(alphaL_*un_state.h(i-1,j) + alphaR_*uo_state.h(i-1,j));
        if( j>0 ) un_state.h(i,j) += p_(i,j)*(alphaL_*un_state.h(i,j-1) + alphaR_*uo_state.h(i,j-1));
  
        un_state.h(i,j) /= (1.0+2.0*alphaL_*p_(i,j));
        un_state.Qout(i,j) = kappa_*un_state.h(i,j);   
        un_state.Qin(i,j)  = 1.0*z_ctrl(i,j);
       }
    }
    S_->apply( un_state, un_state, 1.0, Ntanks_, 2*Ntanks_ );
    value(c,uo,un,z,ts);
  }

  void applyJacobian_uo( V& jv, const V& v,
                         const V& uo, const V& un,
                         const V& z,  const TS& ts) const override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vo_state = to_state(v);
    R_->apply( jv_state, vo_state, -1.0, h_, h_);      
  }

  void applyAdjointJacobian_uo( V &ajv, const V& v,
                                const V &uo, const V &un,
                                const V &z, const TS& ts) const override { 
    ajv.zero();
    auto& ajv_state = to_state(ajv);
    auto& v_state   = to_state(v);
    R_->applyTranspose( ajv_state, v_state, -1.0, h_, h_);  
  }

  //----------------------------------------------------------------------------

  void applyJacobian_un( V& jv, const V& vn,
                         const V& uo, const V& un,
                         const V& z, const TS& ts ) const override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& vn_state = to_state(vn);
    L_->apply( jv_state, vn_state, 1.0, h_, h_ );  
    jv_state.set( vn_state, Qout_, Qout_ );
    jv_state.axpy( -kappa_, vn_state, Qout_, h_ );
    jv_state.set( vn_state, Qin_, Qin_ );
    S_->apply( jv_state, vn_state, -1.0, Qin_, Qout_ );   
  }
 
  void applyAdjointJacobian_un( V& ajv, const V &v,
                                const V &uo, const V& un,
                                const V &z, const TS& ts ) const override {
    ajv.zero();
    auto& ajv_state = to_state(ajv);
    auto& vn_state  = to_state(v);
    L_->applyTranspose( ajv_state, vn_state, 1.0, h_, h_ );  
    ajv_state.axpy(-kappa_, vn_state, h_, Qout_ );
    ajv_state.set( vn_state, Qout_, Qout_ );
    ajv_state.set( vn_state, Qin_, Qin_ );   
    S_->applyTranspose( ajv_state, vn_state, -1.0, Qout_, Qin_ );
  }

  void applyInverseJacobian_un( V &ijv, const V &v,
                                const V &uo, const V &un,
                                const V &z, const TS& ts ) const override {
    ijv.zero();
    auto& ijv_state = to_state(ijv);      
    auto& vn_state  = to_state(v);
    L_->solve( ijv_state, vn_state, 1.0, h_, h_ );  
    ijv_state.set( vn_state, Qout_, Qout_ );
    ijv_state.axpy( kappa_, ijv_state, Qout_, h_ );
    ijv_state.set( vn_state, Qin_, Qin_ );
    S_->apply( ijv_state, ijv_state, 1.0, Qin_, Qout_ );
 }

  void applyInverseAdjointJacobian_un( V& iajv, const V& v,
                                       const V& uo, const V& un,
                                       const V& z, const TS& ts) const override {
    iajv.zero();
    auto& iajv_state = to_state(iajv);      
    auto& vn_state  = to_state(v);
    iajv_state.set( vn_state, Qin_, Qin_ );
    iajv_state.set( vn_state, Qout_, Qout_ );
    S_->applyTranspose( iajv_state, vn_state, 1.0, Qout_, Qin_ );
    scratch_.zero();
    scratch_.set( vn_state, h_, h_ );
    scratch_.axpy( kappa_, iajv_state, h_, Qout_ );
    L_->solveTranspose( iajv_state, scratch_, 1.0, h_, h_ );  

  }


  //----------------------------------------------------------------------------

  void applyJacobian_z( V &jv, const V &v,
                        const V &uo, const V &un,
                        const V &z, const TS& ts) const override {
    jv.zero();
    auto& jv_state = to_state(jv);
    auto& v_ctrl   = to_control(v);
    jv_state.axpy(-beta_,v_ctrl, h_);
    jv_state.hadamard( p_, h_ );
    jv_state.axpy(-1.0, v_ctrl, Qin_ );    
  }

  void applyAdjointJacobian_z( V& ajv, const V& v,
                               const V &uo, const V& un,
                               const V &z, const TS& ts ) const override {
    ajv.zero();
    auto& ajv_ctrl = to_control(ajv);
    auto& v_state  = to_state(v);
    for( size_type i=0; i<rows_; ++i ) 
      for( size_type j=0; j<cols_; ++j ) 
        ajv_ctrl(i,j) = -beta_*p_(i,j)*v_state.h(i,j)-v_state.Qin(i,j);
    }
 
  

}; // Tanks::DynamicConstraint

} // namespace Tanks


#endif // TANKS_DYNAMICCONSTRAINT_HPP

