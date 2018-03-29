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
  coeff1_( Cv_*rho_*g_ ), 
  kappa_( 0.5*coeff1_*dt_/A_ ),
  betaL_( (theta_-1)*dt_/A_ ),
  betaR_( theta_*dt_/A_ ) {
  // ------------- End Initializer List ----------------//

  auto ptrows = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Rows"    );
  auto ptcols = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Columns" );

  vector<Real> w(Ntanks_,1.0);

  vector<size_type> band_index{0, 1, cols_};
  
  for( size_type j=0; j<cols_; ++j ) w.at(j)       = 0.0;
  for( size_type i=0; i<rows_; ++i ) w.at(i*cols_) = 0.0;

  for( size_type i=0; i<static_cast<size_type>(ptrows.size()); ++i ) {
    size_type k = cols_*ptrows.at(i)+ptcols.at(i);
    p_[k] = 0.0;
  }

  vector<Real> band_L0(Ntanks_);       vector<Real> band_R0(Ntanks_);
  vector<Real> band_L1(Ntanks_-1);     vector<Real> band_R1(Ntanks_-1);
  vector<Real> band_Lc(Ntanks_-cols_); vector<Real> band_Rc(Ntanks_-cols_);

  Real alpha_L = (1-theta_)*kappa_;    Real alpha_R = theta_*kappa_;

  band_L0[0] = 1.0-2.0*alpha_L*p_[0];
  band_R0[0] = 1.0-2.0*alpha_R*p_[0];

  for( size_type l=1; l<Ntanks_; ++l ) {
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
 
  vector<vector<Real>> lbands{band_L0,band_L1,band_Lc};
  vector<vector<Real>> rbands{band_R0,band_R1,band_Rc};

  L_ = make_shared<Matrix>( band_index, lbands );
  R_ = make_shared<Matrix>( band_index, rbands );

  print_members(cout);

} // end Constructor

template<typename Real>
void TankState<Real>::print_members( ostream& os ) const {

  os << "Number of rows       = " << rows_   << endl;          
  os << "Number of columns    = " << cols_   << endl;          
  os << "Valve Constant       = " << Cv_     << endl; 
  os << "Density of Fluid     = " << rho_    << endl; 
  os << "Initial Fluid Level  = " << h0_     << endl; 
  os << "Height of Tank       = " << H_      << endl; 
  os << "Cross-sectional Area = " << A_      << endl; 
  os << "Gravity Constant     = " << g_      << endl; 
  os << "Total Time           = " << T_      << endl; 
  os << "Theta                = " << theta_  << endl; 
  os << "Number of Time Steps = " << Nt_     << endl; 
  os << "Ntanks_              = " << Ntanks_ << endl;  
  os << "coeff1_              = " << coeff1_ << endl;
  os << "kappa_               = " << kappa_  << endl;
  os << "betaL_               = " << betaL_  << endl;
  os << "betaR_               = " << betaR_  << endl;

  os << "\nPass Throughs" << endl;
  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      size_type k = cols_*i+j;
      os << p_[k] << " ";
    }
    os << endl;
  }

  os << "\nLHS Matrix" << endl;
  L_->print(os);
  os << "\nRHS Matrix" << endl;
  R_->print(os);

}

template<typename Real>
void TankState<Real>::value( vector<Real>& c, const vector<Real>& u_old, 
                             const vector<Real>& u_new, const vector<Real>& z ) const {

  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {

      size_type l = cols_*i+j;

      auto& h_val    = h(c,i,j);    auto& h_new    = h(u_new,i,j);    auto& h_old    = h(u_old,i,j);
      auto& Qout_val = Qout(c,i,j); auto& Qout_new = Qout(u_new,i,j); auto& Qout_old = Qout(u_old,i,j);
      auto& Qin_val  = Qin(c,i,j);  auto& Qin_new  = Qin(u_new,i,j);  auto& Qin_old  = Qin(u_old,i,j);

      h_val    = h_new - h_old - p_.at(l)*( betaL_*(Qin_new-Qout_new) + betaR_*(Qin_old-Qout_old) );
      Qout_val = Qout_new - coeff1_*h_new;
      Qin_val  = Qin_new - dt_*z.at(l);

      if( i>0 ) Qin_val -= 0.5*Qout(u_new,i-1,j);
      if( j>0 ) Qin_val -= 0.5*Qout(u_new,i,j-1);

    }
      cout << endl;
  }  

}

template<typename Real>
void TankState<Real>::applyJacobian_1_old( vector<Real>& jv, const vector<Real>& v_old ) const {

  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      size_type l = cols_*i+j;

      auto& h_jv    = h(jv,i,j);    auto& h_vo    = h(v_old,i,j);    
      auto& Qout_jv = Qout(jv,i,j); auto& Qout_vo = Qout(v_old,i,j); 
      auto& Qin_jv  = Qin(jv,i,j);  auto& Qin_vo  = Qin(v_old,i,j);  

      h_jv    = - h_vo - p_.at(l)*( betaR_*(Qin_vo-Qout_vo) );
      Qout_jv = 0; 
      Qin_jv  = 0;
    }
  }  
}

template<typename Real>
void TankState<Real>::applyJacobian_1_new( vector<Real>& jv, const vector<Real>& v_new ) const {

  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      size_type l = cols_*i+j;

      auto& h_jv    = h(jv,i,j);    auto& h_vn    = h(v_new,i,j);    
      auto& Qout_jv = Qout(jv,i,j); auto& Qout_vn = Qout(v_new,i,j); 
      auto& Qin_jv  = Qin(jv,i,j);  auto& Qin_vn  = Qin(v_new,i,j);  

      h_jv    = h_vn - p_.at(l)*( betaL_*(Qin_vn-Qout_vn) );
      Qout_jv = Qout_vn - coeff1_*h_vn;
      Qin_jv  = Qin_vn;

      if( i>0 ) Qin_jv -= 0.5*Qout(v_new,i-1,j);
      if( j>0 ) Qin_jv -= 0.5*Qout(v_new,i,j-1);
    }
  }  
}

template<typename Real>
void TankState<Real>::applyJacobian_2( vector<Real> &jv, const vector<Real> &v_new ) const {

  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      size_type l = cols_*i+j;

      auto& h_jv    = h(jv,i,j);    
      auto& Qout_jv = Qout(jv,i,j); 
      auto& Qin_jv  = Qin(jv,i,j);  

      h_jv    = 0;
      Qout_jv = 0;
      Qin_jv  = -v_new.at(l);

      if( i>0 ) Qin_jv -= 0.5*Qout(v_new,i-1,j);
      if( j>0 ) Qin_jv -= 0.5*Qout(v_new,i,j-1);
    }
  }  
}



template<typename Real>
void TankState<Real>::compute_flow(       vector<Real>& u, 
                                    const vector<Real>& z ) const {
  for( size_type i=0; i<rows_; ++i ) {
    for( size_type j=0; j<cols_; ++j ) {
      size_type l = cols_*i+j;
      Qout(u,i,j) = coeff1_*h(u,i,j);   
      Qin(u,i,j)  = dt_*z.at(l);
      if( i>0 ) Qin(u,i,j) += 0.5*Qout(u,i-1,j);
      if( j>0 ) Qin(u,i,j) += 0.5*Qout(u,i,j-1);
    }
  }
  cout << "u_new = ";
  for( auto e: u ) cout << e << " ";
  cout << endl;
} // compute_flow



template<typename Real>
void TankState<Real>::solve_level(       vector<Real>& c, 
                                         vector<Real>& u_new, 
                                   const vector<Real>& u_old, 
                                   const vector<Real>& z     ) const {
  // c += R*u_new
  R_->apply(c,u_old,1.0,0,Ntanks_);

  for( size_type l=0; l<Ntanks_; ++l )  c.at(l) += dt_*p_.at(l)*z.at(l)/A_;

  cout << "c = ";
  for( auto e: c ) cout << e << " ";
  L_->solve(u_new,c,1.0,0,Ntanks_);
  cout << "\nu_new = ";
  for( auto e: u_new ) cout << e << " ";
  cout << endl;
} // solve_level

} // namespace details


#endif // TANKSTATE_IMPL_HPP

