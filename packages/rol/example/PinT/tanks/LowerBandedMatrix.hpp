// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef LOWER_BANDED_MATRIX_HPP
#define LOWER_BANDED_MATRIX_HPP

#include<iostream>
#include<iomanip>

#include "ROL_StdVector.hpp"

namespace details {
using namespace std;

/* \class LowerBandedMatrix
   \brief Implements linear solve and multiplication by
          a three-banded matrix and its transpose when the matrix
          has the specific form \f]

         \f[ a_{ii} = 1 + \alpha*p_i,\; 
             a_{i,i-1}=\beta (1-\delta_{i,c}) p_i,\; 
             a_{i,i-c}=\beta p_i \f]

 */

template<typename Real>
class LowerBandedMatrix {
  
  using size_type = typename vector<Real>::size_type;

private:

  size_type c_;  
  size_type n_;
  Real alpha_;
  Real beta_;


  Real coeff( size_type i ) const { return (1.0-static_cast<Real>(i%c_==0)); }

public:   
  LowerBandedMatrix(){}
  LowerBandedMatrix( size_type r, size_type c, Real alpha, Real beta ) :
    c_{c}, n_{r*c}, alpha_{alpha}, beta_{beta} {}
  virtual ~LowerBandedMatrix() {}

  // Compute b += \kappa*A*x 
  // with the elements of b and x being offset by the arguments Ax_begin and
  // x_begin, if the vectors are larger than A
  virtual void apply( vector<Real>& b, const vector<Real> &x, Real kappa=1.0, 
                      size_type bi0=0, size_type xi0=0 ) const {

    b.at(bi0)   += kappa*x.at(xi0)*alpha_;

    for( size_type i=1; i<c_; ++i ) 
      b.at(bi0+i) += kappa*(x.at(xi0+i)*alpha_ + beta_*x.at(xi0+i-1));

    for( size_type i=c_; i<n_; ++i ) 
      b.at(bi0+i) += kappa*(x.at(xi0+i)*alpha_ + 
                             beta_*(coeff(i)*x.at(xi0+i-1) + x.at(xi0+i-c_)));
  } 

  void apply( ROL::StdVector<Real>& b, const ROL::StdVector<Real> &x, Real kappa=1.0, 
                      size_type bi0=0, size_type xi0=0 ) const {
    apply( *(b.getVector()), *(x.getVector()), kappa, bi0, xi0 ); 
  }


  virtual void solve( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                      size_type xi0=0, size_type bi0=0 ) const {
    x.at(xi0) += b.at(bi0)*kappa/alpha_;
    
    for( size_type i=1; i<c_; ++i ) 
      x.at(xi0+i) += ( b.at(bi0+i) - beta_*x.at(xi0+i-1) )*kappa/alpha_;

    for( size_type i=c_; i<n_; ++i ) 
      x.at(xi0+i) += ( b.at(bi0+i) - beta_*(coeff(i)*x.at(xi0+i-1) 
                                            +x.at(xi0+i-c_)))*kappa/alpha_;
  }

  void solve( ROL::StdVector<Real>& x, const ROL::StdVector<Real>& b, Real kappa=1.0,
                      size_type xi0=0, size_type bi0=0 ) const { 
    solve( *(x.getVector()), *(b.getVector()), kappa, xi0, bi0 );
  }


  virtual void applyTranspose( vector<Real>&b, const vector<Real>& x, Real kappa=1.0,
                               size_type bi0=0, size_type xi0=0 ) const {

    b.at(bi0+n_-1) += kappa*x.at(xi0+n_-1)*alpha_;
    for( size_type j=1; j<c_; ++j ) { 
      size_type i=n_-j-1;
      b.at(bi0+i) += kappa*(x.at(xi0+i)*alpha_ + beta_*x.at(xi0+i+1)); 
    }

    for( size_type j=c_; j<n_; ++j ) {
      size_type i=n_-j-1;
      b.at(bi0+i) += kappa*(x.at(xi0+i)*alpha_ + 
                            beta_*(coeff(i+1)*x.at(xi0+i+1) + x.at(xi0+i+c_))); 
    }
  }

  void applyTranspose( ROL::StdVector<Real>&b, const ROL::StdVector<Real>& x, Real kappa=1.0,
                               size_type bi0=0, size_type xi0=0 ) const {
    applyTranspose( *(b.getVector()), *(x.getVector()), kappa, bi0, xi0 );
  }


  virtual void solveTranspose( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                               size_type xi0=0, size_type bi0=0 ) const {

    x.at(xi0+n_-1) += b.at(bi0+n_-1)*kappa/alpha_;
    for( size_type j=1; j<c_; ++j ) { 
      size_type i=n_-j-1;
      x.at(xi0+i) += ( b.at(bi0+i) - beta_*x.at(xi0+i+1) )*kappa/alpha_;
    }
    for( size_type j=c_; j<n_; ++j ) {
      size_type i=n_-j-1;
      x.at(xi0+i) += ( b.at(bi0+i) - beta_*( coeff(i+1)*x.at(xi0+i+1) +
                       x.at(xi0+i+c_) ) )*kappa/alpha_;
    }
  }

  void solveTranspose( ROL::StdVector<Real>& x, const ROL::StdVector<Real>& b, Real kappa=1.0,
                               size_type xi0=0, size_type bi0=0 ) const {
    solveTranspose( *(x.getVector()), *(b.getVector()), kappa, xi0, bi0 );
  }

}; // class LowerBandedMatrix


/** Matrix multiplication for a matrix with three lower bands
    that have the special structure:

    \f[   a_{i,i-1}  = 0.5 (1-\delta_{i,c}) \;
          a_{i,i-c}  = 0.5 
    \f]
*/
template<typename Real> 
class SplitterMatrix : public LowerBandedMatrix<Real> {

  using size_type = typename vector<Real>::size_type;

private:
  size_type c_;  
  size_type n_;

  Real coeff( size_type i ) const { return (1.0-static_cast<Real>(i%c_==0)); }

public:
  SplitterMatrix( size_type r, size_type c ) : c_{c}, n_{r*c} {}

  void apply( vector<Real>& b, const vector<Real> &x, Real kappa=1.0, 
                      size_type bi0=0, size_type xi0=0 ) const override {

    for( size_type i=1; i<c_; ++i ) 
      b.at(bi0+i) += 0.5*kappa*x.at(xi0+i-1)*coeff(i);

    for( size_type i=c_; i<n_; ++i ) 
      b.at(bi0+i) += 0.5*kappa*( coeff(i)*x.at(xi0+i-1) + x.at(xi0+i-c_) );
  } 

  void solve( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                      size_type xi0=0, size_type bi0=0 ) const override {
    throw ROL::Exception::NotImplemented("SplitterMatrix is singular.");
  }

  void applyTranspose( vector<Real>&b, const vector<Real>& x, Real kappa=1.0,
                               size_type bi0=0, size_type xi0=0 ) const override {

    for( size_type j=1; j<c_; ++j ) { 
      size_type i=n_-j-1;
      b.at(bi0+i) += 0.5*kappa*x.at(xi0+i+1); 
    }
    for( size_type j=c_; j<n_; ++j ) {
      size_type i=n_-j-1;
      b.at(bi0+i) += 0.5*kappa*(coeff(i+1)*x.at(xi0+i+1) + x.at(xi0+i+c_)); 
    }
  }

  void solveTranspose( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                               size_type xi0=0, size_type bi0=0 ) const override {
    throw ROL::Exception::NotImplemented("SplitterMatrix is singular.");
  }
};


/** Linear solve and matrix multiplication for a matrix with three lower bands
    that have the special structure:

    \f[   a_{i,i} = 1 + 2\alpha p_{i}, \;
          a_{i,i-1} = -\alpha p_i(1-\delta_{c\%i,0}),\;
          a_{i,i-c}  = -alpha p_i    \f]
*/

template<typename Real>
class TankLevelMatrix : public LowerBandedMatrix<Real> {
  using size_type = typename vector<Real>::size_type;

private:

  size_type c_;  
  size_type n_;
  vector<Real> alpha_;

  Real coeff( size_type i ) const { return (1.0-static_cast<Real>(i%c_==0)); }

  Real diag( size_type i ) const { return 1.0+2.0*alpha_.at(i); }

public:   

  TankLevelMatrix( size_type r, size_type c, Real alpha, const vector<Real>& p ) :
    c_{c}, n_{r*c}, alpha_(n_) {
    for( size_type i=0; i<n_; ++i ) alpha_.at(i) = alpha*p.at(i);
  }

  // Compute b += \kappa*A*x 
  // with the elements of b and x being offset by the arguments Ax_begin and
  // x_begin, if the vectors are larger than A
  virtual void apply( vector<Real>& b, const vector<Real> &x, Real kappa=1.0, 
                      size_type bi0=0, size_type xi0=0 ) const {
    for( size_type i=0; i<n_; ++i ) {
      auto value = diag(i)*x.at(xi0+i);
      if( i>=1 ) value -= alpha_.at(i)*coeff(i)*x.at(xi0+i-1);
      if( i>=c_) value -= alpha_.at(i)*x.at(xi0+i-c_);
      b.at(bi0+i) += kappa*value;
    }
  } 

  virtual void solve( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                      size_type xi0=0, size_type bi0=0 ) const {
    for( size_type i=0; i<n_; ++i ) {
      auto value = b.at(bi0+i);
      if( i>=1 ) value  += alpha_.at(i)*coeff(i)*x.at(xi0+i-1);
      if( i>=c_ ) value += alpha_.at(i)*x.at(xi0+i-c_);
      x.at(xi0+i) += kappa*value/diag(i);
    }
  }

  virtual void applyTranspose( vector<Real>&b, const vector<Real>& x, Real kappa=1.0,
                               size_type bi0=0, size_type xi0=0 ) const {

    for( size_type j=0; j<n_; ++j ) {
      size_type i = n_-j-1;
      auto value = diag(i)*x.at(xi0+i);
      if( i<n_-1  ) value -= alpha_.at(i+1)*x.at(xi0+i+1)*coeff(i+1);
      if( i<n_-c_ ) value -= alpha_.at(i+c_)*x.at(xi0+i+c_);
      b.at(bi0+i) += kappa*value;
    }
  }

  virtual void solveTranspose( vector<Real>& x, const vector<Real>& b, Real kappa=1.0,
                               size_type xi0=0, size_type bi0=0 ) const {
    for( size_type j=0; j<n_; ++j ) {
      size_type i = n_-j-1;
       auto value = b.at(bi0+i);
       if( i<n_-1 )  value += alpha_.at(i+1)*coeff(i+1)*x.at(xi0+i+1);
       if( i<n_-c_ ) value += alpha_.at(i+c_)*x.at(xi0+i+c_);
       x.at(xi0+i) += kappa*value/diag(i);
      }
    }

};


} // namespace details

using details::LowerBandedMatrix;
using details::SplitterMatrix; 
using details::TankLevelMatrix;


#endif // LOWER_BANDED_MATRIX_HPP

