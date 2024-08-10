// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_MINRES_HPP
#define ROL_MINRES_HPP

#include <array>
#include "ROL_Krylov.hpp"
#include "ROL_VectorClone.hpp"

namespace ROL {

  /** \class ROL::MINRES
    \brief Implements the MINRES algorithm for solving symmetric indefinite
    systems 
   */


namespace details {

using namespace std;

template<typename Real>
class MINRES : public Krylov<Real> {

  using V  = Vector<Real>;
  using OP = LinearOperator<Real>;

private:

  // Givens rotation matrix elements
  Real resnorm_;
  int maxiter_;  
  bool useInexact_;
  array<Real,4> H_;
  array<Real,2> rhs_;

  VectorCloneMap<Real> clones_;      

  void givens( Real& c, Real& s, Real& r, Real a, Real b ) const {

    Real zero(0), one(1);

    if( b == zero ) {
      c = ( a >= zero ? one : -one );
      s = zero;
      r = abs(a);
    }
    else if( a == zero ) {
      c = zero;
      s = ( b >= zero ? one : -one );
      r = abs(b);
    }
    else if( abs(a) > abs(b) ) {
      auto t = b/a;
      auto u = copysign(sqrt(one+t*t),a);
      c = one/u;
      s = c*t;
      r = a*u;
    }
    else {
      auto t = a/b;
      auto u = copysign(sqrt(one+t*t),b);
      s = 1/u;
      c = s*t;
      r = b*u;
    }
  } // givens()

public:

  MINRES( Real absTol = 1.e-4, Real relTol = 1.e-2, unsigned maxit = 100, bool useInexact = false) :
    Krylov<Real>(absTol,relTol,maxit), useInexact_(useInexact),
    clones_("v_prev","v_curr","v_next","w_prev","w_curr","w_next") { }

  // Note: Preconditioner is not implemented
  virtual Real run( V &x, OP &A, const V &b, OP &M, int &iter, int &flag ) override {

    auto v_prev = clones_( x, "v_prev" );  v_prev->zero(); 
    auto v_curr = clones_( x, "v_curr" );  v_curr->set(b);
    auto v_next = clones_( x, "v_next" );  v_next->zero();
    auto w_prev = clones_( x, "w_prev" );  w_prev->zero();
    auto w_curr = clones_( x, "w_curr" );  w_curr->zero();
    auto w_next = clones_( x, "w_next" );  w_next->zero();

    Real c_prev{0}, s_prev{0}, c_curr{0}, s_curr{0}, c_next{0}, s_next{0};

    resnorm_ = v_curr->norm();    
    Real rtol = min(Krylov<Real>::getAbsoluteTolerance(),Krylov<Real>::getRelativeTolerance()*resnorm_);
    Real itol = sqrt(ROL_EPSILON<Real>());

    for( auto &e: H_ ) e = 0;

    rhs_[0] = resnorm_; rhs_[1] = 0;

    v_curr->scale(1.0/resnorm_);

    for( iter=0;  iter < (int)Krylov<Real>::getMaximumIteration(); iter++) {
      if ( useInexact_ ) {
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * resnorm_);
      }

      if( resnorm_ < rtol ) break;

      A.apply( *v_next, *v_curr, itol );

      if( iter>0 ) v_next->axpy(-H_[1],*v_prev);    

      H_[2] = v_next->dot(*v_curr);

      v_next->axpy(-H_[2],*v_curr);

      H_[3] = v_next->norm();

      v_next->scale(1.0/H_[3]);

      // Rotation on H_[0] and H_[1]. 
      if( iter>1 ) {
        H_[0]  = s_prev*H_[1];
        H_[1] *= c_prev;
      }     

      // Rotation on H_[1] and H_[2]
      if( iter>0 ) {
        auto tmp = c_curr*H_[2]-s_curr*H_[1];
        H_[1] = c_curr*H_[1] + s_curr*H_[2];
        H_[2] = tmp;
      }

      // Form rotation coefficients
      givens( c_next, s_next, H_[2], H_[2], H_[3] );
      rhs_[1]  = -s_next*rhs_[0];
      rhs_[0] *=  c_next;

      w_next->set( *v_curr );

      if( iter>0 )  w_next->axpy( -H_[1], *w_curr );
      if( iter>1 )  w_next->axpy( -H_[0], *w_prev );

      w_next->scale(1.0/H_[2]);

      x.axpy( rhs_[0], *w_next );

      v_prev->set( *v_curr ); 
      v_curr->set( *v_next ); 
      w_prev->set( *w_curr ); 
      w_curr->set( *w_next ); 

      c_prev = c_curr;
      c_curr = c_next;
      s_prev = s_curr;
      s_curr = s_next;

      rhs_[0] = rhs_[1];

      H_[1] = H_[3];

      resnorm_ = abs( rhs_[1] );

    } // for (iter)

    if ( iter == (int)Krylov<Real>::getMaximumIteration() ) flag = 1;
    else iter++;

    return resnorm_;
  } // run()

}; // class MINRES

} // namespace details


using details::MINRES;


} // namespace ROL


#endif // ROL_MINRES_HPP

