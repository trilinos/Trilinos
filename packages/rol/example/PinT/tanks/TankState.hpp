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
#ifndef TANKSTATE_HPP
#define TANKSTATE_HPP

#include <memory>

#include "Teuchos_ParameterList.hpp"
#include "LowerBandedMatrix.hpp"


namespace details {

using namespace std;

template<typename Real> 
class TankState {
private:

  using Matrix    = LowerBandedMatrix<Real>;
  using size_type = typename vector<Real>::size_type;

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
  Real dt_;                    // Time step size
  Real theta_;                 // Implicit/Explicit splitting factor
  int  Nt_;                    // Number of Time steps

  size_type    Ntanks_;        // Total number of tanks 
  vector<Real> p_;             // Passthrough coefficients
 
  Real         coeff1_;        // Cv*rho*g
  Real         kappa_;         // Cv*rho*g*dt/2A

  shared_ptr<Matrix> L_, R_;

public:
   
  TankState( Teuchos::ParameterList& pl ) :
    // ----------- Begin Initializer List ----------------//
    rows_   ( pl.get( "Number of Rows",          3      ) ),
    cols_   ( pl.get( "Number of Columns",       3      ) ),
    Cv_     ( pl.get( "Valve Constant",          1.0e-2 ) ),
    rho_    ( pl.get( "Density of Fluid",        1.0e3  ) ),
    h0_     ( pl.get( "Initial Fluid Level",     2.0    ) ),
    H_      ( pl.get( "Height of Tank",          10.0   ) ),
    A_      ( pl.get( "Cross-sectional Area",    10.0   ) ),
    g_      ( pl.get( "Gravity Constant",        9.8    ) ),
    T_      ( pl.get( "Total Time",              20.0   ) ),
    theta_  ( pl.get( "Theta",                   0.5    ) ),
    Nt_     ( pl.get( "Number of Time Steps",    100    ) ),
    //----------------------------------------------------//
    Ntanks_(rows_*cols_),  p_(Ntanks_,1.0), 
    coeff1_( Cv_*rho_*g_ ), 
    kappa_( 0.5*coeff1_*dt_/A_ ) {
    // ------------- End Initializer List ----------------//

    auto ptrows = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Rows"    );
    auto ptcols = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Columns" );

    vector<Real> w(Ntanks_,1.0);

    vector<size_type> band_index{0, 1, cols_};
    
    for( size_type j=0; j<cols_; ++j ) w[j]       = 0.0;
    for( size_type i=0; i<rows_; ++i ) w[i*cols_] = 0.0;

    for( size_type i=0; i<ptrows.size(); ++i ) {
      size_type k = cols_*ptrows[i]+ptcols[i];
      p_[k] = 0.0;
    }

    vector<Real> band_L0(Ntanks_);       vector<Real> band_R0(Ntanks_);
    vector<Real> band_L1(Ntanks_-1);     vector<Real> band_R1(Ntanks_-1);
    vector<Real> band_Lc(Ntanks_-cols_); vector<Real> band_Rc(Ntanks_-cols_);

    Real alpha_L = (1-theta_)*kappa_;    Real alpha_R = theta_*kappa_;

    band_L0[0] = 1.0-2.0*alpha_L*p_[0];
    band_R0[0] = 1.0-2.0*alpha_R*p_[0];

    for( size_type l=1; l<cols_; ++l ) {
      band_L0[l]   = 1.0-2.0*alpha_L*p_[l];
      band_R0[l]   = 1.0-2.0*alpha_R*p_[l];

      if( l>=1 ) {
        
        band_L1[l-1] = alpha_L*w[l]*p_[l] * (l%cols_!=0);
        band_R1[l-1] = alpha_R*w[l]*p_[l] * (l%cols_!=0);

        if( l>=cols_ ) {
          band_Lc[l-cols_] = alpha_L*w[l]*p_[l];
          band_Rc[l-cols_] = alpha_R*w[l]*p_[l];
        }
      } 
    } // end for
  } // end Constructor

  // Subvector Accessor Methods
  Real& h( vector<Real>& x, size_type r, size_type c ) const {
    return x[cols_*r+c];
  }

  const Real& h( const vector<Real>& x, size_type r, size_type c ) const {
    return x[cols_*r+c];
  }

  Real& Qout( vector<Real>& x, size_type r, size_type c ) const {
    return x[Ntanks_+cols_*r+c];
  }

  const Real& Qout( const vector<Real>& x, size_type r, size_type c ) const {
    return x[Ntanks_+cols_*r+c];
  }

  Real& Qin( vector<Real>& x, size_type r, size_type c ) const {
    return x[2*Ntanks_+cols_*r+c];
  }

  const Real& Qin( const vector<Real>& x, size_type r, size_type c ) const {
    return x[2*Ntanks_+cols_*r+c];
  }

  // Given a vector u=(h,Qin,Qout) and f, update Qin and Qout from h
  void compute_flow( vector<Real>& u, const vector<Real>& f ) const {
    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        size_type l = cols_*i+j;
        Qout(u,i,j) = coeff1_*h(u,i,j);   
        Qin(u,i,j)  = f(i,j);
        if( i>0 ) Qin(u,i,j) += 0.5*Qout(u,i-1,j);
        if( j>0 ) Qin(u,i,j) += 0.5*Qout(u,i,j-1);
      }
    }
  }

  // This requires the input f to already be theta-weighted
  void residual_level( vector<Real>& c, 
                       const vector<Real>& u_new, 
                       const vector<Real>& u_old, 
                       const vector<Real>& f ) const {
    // c += L*u_new
    L_->apply(c,u_new,1.0,0,Ntanks_);

    // c -= R*u_new
    R_->apply(c,u_old,-1.0,0,Ntanks_);

    for( size_type l=0; l<Ntanks_; ++l )  c[l] += dt_*p_[l]*f[l]/A_;
   
  }
  

  void solve_level( vector<Real>& c, 
                    vector<Real>& u_new, 
                    const vector<Real>& u_old, 
                    const vector<Real>& f ) const {
    // c += R*u_new
    R_->apply(c,u_old,1.0,0,Ntanks_);

    for( size_type l=0; l<Ntanks_; ++l )  c[l] += dt_*p_[l]*f[l]/A_;
 
    L_->solve(u_new,c,1.0,0,Ntanks_);
       
  }

  
 

}; // class TankState


} // namespace details

using details::TankState;

#endif // TANKSTATE_HPP

