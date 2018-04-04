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
#ifndef TANKSTATE_IMPL_HPP
#define TANKSTATE_IMPL_HPP

namespace details {

using namespace std;

template<typename Real>
TankState<Real>::TankState( Teuchos::ParameterList& pl ) :
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
  dt_( T_/Nt_ ),
  Ntanks_(rows_*cols_),
  p_(Ntanks_,1.0), 
  kappa_( Cv_*rho_*g_ ), 
  beta_( dt_/A_ ),
//  betaL_( (1-theta_)*dt_/A_ ),
//  betaR_( theta_*dt_/A_ ),
  alphaL_( 0.5*kappa_*(1-theta_)*beta_ ),
  alphaR_( 0.5*kappa_*theta_*beta_ ) {
  // ------------- End Initializer List ----------------//

  auto ptrows = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Rows"    );
  auto ptcols = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Columns" );

  vector<size_type> band_index{0, 1, cols_};
  vector<size_type> shift_index{1, cols_};
  
  for( size_type i=0; i<static_cast<size_type>(ptrows.size()); ++i ) {
    p_.at( cols_*ptrows.at(i)+ptcols.at(i) ) = 0.0;
  }
 
//  vector<Real> band_D1(Ntanks_-1); 
//  vector<Real> band_Dc(Ntanks_-cols_);  
//
//  vector<Real> band_L0( Ntanks_ );       vector<Real> band_R0( Ntanks_ );
//  vector<Real> band_L1( Ntanks_-1 );     vector<Real> band_R1( Ntanks_-1 );
//  vector<Real> band_Lc( Ntanks_-cols_ ); vector<Real> band_Rc( Ntanks_-cols_ );
//
//  for( size_type l=0; l<Ntanks_; ++l ) {
//    band_L0.at(l) = 1.0 + 2.0*alphaL_*p_.at(l);
//    band_R0.at(l) = 1.0 - 2.0*alphaR_*p_.at(l);
//
//    if( l>=1 ) {
//      band_D1.at(l-1) = 0.5*(l%cols_!=0);
//      band_L1.at(l-1) = -alphaL_*p_.at(l)*(l%cols_!=0);
//      band_R1.at(l-1) =  alphaR_*p_.at(l)*(l%cols_!=0);
//
//      if( l>=cols_ ) {
//        band_Dc.at(l-cols_) = 0.5*(l%cols_!=0);
//        band_Lc.at(l-cols_) = -alphaL_*p_.at(l);
//        band_Rc.at(l-cols_) =  alphaR_*p_.at(l);
//      }
//    } 
//  } // end for
// 
//  vector<vector<Real>> lbands{ band_L0, band_L1, band_Lc };
//  vector<vector<Real>> rbands{ band_R0, band_R1, band_Rc };
//  vector<vector<Real>> dbands{ band_D1, band_Dc };

  L_ = make_shared<TankLevelMatrix<Real>>( rows_, cols_, alphaL_, p_ );
//  R_ = make_shared<Matrix>( band_index, rbands );
  S_ = make_shared<SplitterMatrix<Real>>( rows_, cols_ );
//  print_members(cout);

} // end Constructor

template<typename Real>
void TankState<Real>::solve( StateVector& c, StateVector& u_new, 
                             const StateVector& u_old, const ControlVector& z ) const {

  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      u_new.h(i,j) =  u_old.h(i,j) + p(i,j)*(beta_*z(i,j)-2.0*alphaR_*u_old.h(i,j));

      if( i>0 ) u_new.h(i,j) += p(i,j)*(alphaL_*u_new.h(i-1,j) + alphaR_*u_old.h(i-1,j));
      if( j>0 ) u_new.h(i,j) += p(i,j)*(alphaL_*u_new.h(i,j-1) + alphaR_*u_old.h(i,j-1));
      u_new.h(i,j) /= (1.0+2.0*alphaL_*p(i,j));

      u_new.Qout(i,j) = kappa_*u_new.h(i,j);   
      u_new.Qin(i,j)  = 1.0*z(i,j);
      if( i>0 ) u_new.Qin(i,j) += 0.5*u_new.Qout(i-1,j);
      if( j>0 ) u_new.Qin(i,j) += 0.5*u_new.Qout(i,j-1);
     }
  }
} 

template<typename Real>
void TankState<Real>::value( StateVector& c, const StateVector& u_old, 
                             const StateVector& u_new, const ControlVector& z ) const {

  L_->apply( *(c.getVector()), *(u_new.getVector()), 1.0, 0, 0 ); 

  for(size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      c.h(i,j) += // (1.0+2.0*alphaL_*p(i,j))*u_new.h(i,j) 
                   - u_old.h(i,j) 
                   - p(i,j)*( beta_*z(i,j)-2.0*alphaR_*u_old.h(i,j) );
      c.Qout(i,j) = u_new.Qout(i,j) - kappa_*u_new.h(i,j);
      c.Qin(i,j)  = u_new.Qin(i,j)  - z(i,j);
 
      if( i>0 ) {
        c.h(i,j)   -=// p(i,j)*(alphaL_*u_new.h(i-1,j) + 
                        p(i,j)*alphaR_*u_old.h(i-1,j);
//        c.Qin(i,j) -= 0.5*u_new.Qout(i-1,j);
      }
      if( j>0 ) {
        c.h(i,j)   -= //p(i,j)*(alphaL_*u_new.h(i,j-1) + 
                      p(i,j)*alphaR_*u_old.h(i,j-1);
//        c.Qin(i,j) -= 0.5*u_new.Qout(i,j-1);
      }
    }
  }
  S_->apply( *(c.getVector()), *(u_new.getVector()), -1.0, Ntanks_, 2*Ntanks_ );
}

template<typename Real>
void TankState<Real>::applyJacobian_1_old( StateVector& jv, const StateVector& v_old ) const {
  R_->apply( *(jv.getVector()), *(v_old.getVector()), -1.0, 0, Ntanks_ );  
}

template<typename Real>
void TankState<Real>::applyJacobian_1_new( StateVector& jv, const StateVector& v_new ) const {

  L_->apply( *(jv.getVector()), *(v_new.getVector()), 1.0, 0, Ntanks_ );  
 
  for(size_type i=0; i<rows_; ++i ) {                                  
    for( size_type j=0; j<cols_; ++j ) {        
                       
      jv.Qout(i,j) = v_new.Qout(i,j) - kappa_*v_new.h(i,j);                                  
      jv.Qin(i,j)  = v_new.Qin(i,j);                                   
                                                                     
      if( i>0 ) jv.Qin(i,j) -= 0.5*v_new.Qout(i-1,j);                           
      if( j>0 ) jv.Qin(i,j) -= 0.5*v_new.Qout(i,j-1);                           
    }
  }
}

template<typename Real>
void TankState<Real>::applyInverseJacobian_1_new( StateVector& ijv, const StateVector& v_new ) const {
  L_->solve( *(ijv.getVector()), *(v_new.getVector()), 1.0, 0, Ntanks_ );  
  for(size_type i=0; i<rows_; ++i ) {                                  
    for( size_type j=0; j<cols_; ++j ) {                               
      ijv.Qout(i,j) = v_new.Qout(i,j) + kappa_*ijv.h(i,j);
      ijv.Qin(i,j)  = v_new.Qin(i,j);
      if( i>0 ) {                                                       
        ijv.Qin(i,j) += 0.5*ijv.Qout(i-1,j);                           
      }                                                                 
      if( j>0 ) {                                                       
        ijv.Qin(i,j) += 0.5*ijv.Qout(i,j-1);                           
      }
    }
  }
}

template<typename Real>
void TankState<Real>::applyJacobian_2( StateVector& jv, const ControlVector &v_new ) const {
  for(size_type i=0; i<rows_; ++i ) {                                  
    for( size_type j=0; j<cols_; ++j ) {                               
      jv.h(i,j) =  (1.0+2.0*alphaL_*p(i,j))*v_new.h(i,j);              
      jv.Qout(i,j) = v_new.Qout(i,j);                                  
      jv.Qin(i,j)  = v_new.Qin(i,j);                                   
                                                                       
     if( i>0 ) {                                                       
       jv.h(i,j) -= p(i,j)*alphaL_*v_new.h(i-1,j);                     
       jv.Qin(i,j) -= 0.5*v_new.Qout(i-1,j);                           
     }                                                                 
     if( j>0 ) {                                                       
       jv.h(i,j) -= p(i,j)*alphaL_*v_new.h(i,j-1);                     
       jv.Qin(i,j) -= 0.5*v_new.Qout(i,j-1);                           
     }
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

