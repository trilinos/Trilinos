// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKSTATE_IMPL_HPP
#define TANKSTATE_IMPL_HPP

namespace details {

using namespace std;

template<typename Real>
TankState<Real>::TankState( ROL::ParameterList& pl ) :
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
  h_(0), Qout_(2*Ntanks_), Qin_(Ntanks_), z_(0) {
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

} // end Constructor

template<typename Real>
void TankState<Real>::solve( StateVector& c, StateVector& u_new, 
                             const StateVector& u_old, const ControlVector& z ) const {
  
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

template<typename Real>
void TankState<Real>::value( StateVector& c, const StateVector& u_old, 
                             const StateVector& u_new, const ControlVector& z ) const {

  c.axpy( -beta_, z, h_ );
  c.hadamard( p_, h_ );
  c.set( u_new, Qin_, Qin_ );
  c.axpy( -1.0, z, Qin_ );
  c.set( u_new, Qout_, Qout_ );
  c.axpy( -kappa_, u_new, Qout_, h_ );

  L_->apply( c, u_new,  1.0, h_, h_ ); 
  R_->apply( c, u_old, -1.0, h_, h_ );
  S_->apply( c, u_new, -1.0, Qin_, Qout_ );

}



template<typename Real>
void TankState<Real>::applyJacobian_1_old( StateVector& jv, const StateVector& v_old ) const {
  R_->apply( jv, v_old, -1.0, h_, h_);  
}

template<typename Real>
void TankState<Real>::applyAdjointJacobian_1_old( StateVector& jv, const StateVector& v_old ) const {
  R_->applyTranspose( jv, v_old, -1.0, h_, h_);  
}

template<typename Real>
void TankState<Real>::applyJacobian_1_new( StateVector& jv, const StateVector& v_new ) const {

  L_->apply( jv, v_new, 1.0, h_, h_ );  
  jv.set( v_new, Qout_, Qout_ );
  jv.axpy( -kappa_, v_new, Qout_, h_ );
  jv.set( v_new, Qin_, Qin_ );
  S_->apply( jv, v_new, -1.0, Qin_, Qout_ );
}

template<typename Real>
void TankState<Real>::applyInverseJacobian_1_new( StateVector& ijv, const StateVector& v_new ) const {

  L_->solve( ijv, v_new, 1.0, h_, h_ );  
  ijv.set( v_new, Qout_, Qout_ );
  ijv.axpy( kappa_, ijv, Qout_, h_ );
  ijv.set( v_new, Qin_, Qin_ );
  S_->apply( ijv, ijv, 1.0, Qin_, Qout_ );
}

template<typename Real>
void TankState<Real>::applyAdjointJacobian_1_new( StateVector& jv, const StateVector& v_new ) const {

  L_->applyTranspose( jv, v_new, 1.0, h_, h_ );  
  jv.axpy(-kappa_, v_new, h_, Qout_ );
  jv.set( v_new, Qout_, Qout_ );
  jv.set( v_new, Qin_, Qin_ );
  S_->applyTranspose( jv, v_new, -1.0, Qout_, Qin_ );
}

template<typename Real>
void TankState<Real>::applyInverseAdjointJacobian_1_new( StateVector& iajv, const StateVector& v_new ) const {

  iajv.set( v_new, Qin_, Qin_ );
  iajv.set( v_new, Qout_, Qout_ );
  S_->applyTranspose( iajv, v_new, 1.0, Qout_, Qin_ );
  scratch_.zero();
  scratch_.set( v_new, h_, h_ );
  scratch_.axpy( kappa_, iajv, h_, Qout_ );
  L_->solveTranspose( iajv, scratch_, 1.0, h_, h_ );  
}

template<typename Real>
void TankState<Real>::applyJacobian_2( StateVector& jv, const ControlVector &v ) const {
  jv.axpy(-beta_,v, h_);
  jv.hadamard( p_, h_ );
  jv.axpy(-1.0, v, Qin_ );
}

template<typename Real>
void TankState<Real>::applyAdjointJacobian_2( ControlVector& ajv, const StateVector &v ) const {
  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      ajv(i,j) = -beta_*p_(i,j)*v.h(i,j)-v.Qin(i,j);
    }
  }
}



template<typename Real>
void TankState<Real>::print_members( ostream& os ) const {
 
  os << endl;
  os << "+---------------------------------+" << endl;
  os << "|     Tank System Parameters      |" << endl;
  os << "+---------------------------------+" << endl;
  os << "  Number of rows       = " << rows_    << endl;          
  os << "  Number of columns    = " << cols_    << endl;          
  os << "  Valve Constant       = " << Cv_      << endl; 
  os << "  Density of Fluid     = " << rho_     << endl; 
  os << "  Initial Fluid Level  = " << h0_      << endl; 
  os << "  Height of Tank       = " << H_       << endl; 
  os << "  Cross-sectional Area = " << A_       << endl; 
  os << "  Gravity Constant     = " << g_       << endl; 
  os << "  Total Time           = " << T_       << endl; 
  os << "  Theta                = " << theta_   << endl; 
  os << "  Number of Time Steps = " << Nt_      << endl; 
  os << "  Ntanks_              = " << Ntanks_  << endl;  
  os << "  kappa_               = " << kappa_   << endl;
  os << "  alphaL_              = " << alphaL_  << endl;
  os << "  alphaR_              = " << alphaR_  << endl;

  os << "\nPass-through tanks (indicated by 0)" << endl;
  for( size_type i=0; i<rows_; ++i ) {
    os << "  ";
    for( size_type j=0; j<cols_; ++j ) {
      size_type k = cols_*i+j;
      os << p_[k] << " ";
    }
    os << endl;
  }

//  os << "\nLHS Matrix" << endl;
//  L_->print(os);
//  os << "\nRHS Matrix" << endl;
//  R_->print(os);

}



} // namespace details


#endif // TANKSTATE_IMPL_HPP

