// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDTRIDIAGONALOPERATOR_H
#define ROL_STDTRIDIAGONALOPERATOR_H

#include "ROL_StdLinearOperator.hpp"

/** @ingroup func_group
    \class ROL::StdTridiagonalOperator
    \brief Provides the std::vector implementation to apply a linear operator,
      which encapsulates a tridiagonal matrix

    \f[ T = \begin{matrix} 
    a_0 & b_0 \\
    c_0 & a_1 & b_1 \\
        & c_1 & \ddots & \ddots \\
        &     & \ddots & a_{n-2} & b_{n-2} \\
        &     &        & c_{n-2} & a_{n-1}
    \end{pmatrix} \f]
    ---
*/


namespace ROL {

template <class Real>
class StdTridiagonalOperator : public StdLinearOperator<Real> {
 
  template <typename T> using vector = std::vector<T>;
  
private:

  const ROL::Ptr<const vector<Real> > a_; // Diagonal
  const ROL::Ptr<const vector<Real> > b_; // Superdiagonal
  const ROL::Ptr<const vector<Real> > c_; // Subdiagonal

  mutable vector<Real> dl_;
  mutable vector<Real> d_;
  mutable vector<Real> du_;
  mutable vector<Real> du2_;

  mutable vector<int>  ipiv_;

  int N_;
  
  ROL::LAPACK<int,Real>  lapack_;

  void copy(void) const {
    for(int i=0;i<N_-1;++i) {
      dl_[i]  = (*c_)[i];
      d_[i]   = (*a_)[i];
      du_[i]  = (*b_)[i]; 
    }  
    d_[N_-1] = (*a_)[N_-1];
    du2_.assign(N_-2,0.0);
  }

public:

  StdTridiagonalOperator( const ROL::Ptr<const vector<Real> > &a,
                          const ROL::Ptr<const vector<Real> > &b,
                          const ROL::Ptr<const vector<Real> > &c ) : 
    StdLinearOperator<Real>(), a_(a), b_(b), c_(c) { 

    N_ = a_->size();
    
    dl_.resize(N_-1);
    d_.resize(N_);
    du_.resize(N_-1);
    du2_.resize(N_-2);
    ipiv_.resize(N_);   
  }

  StdTridiagonalOperator( const ROL::Ptr<const vector<Real> > &a,
                          const ROL::Ptr<const vector<Real> > &b ) {
    StrdTridiagonalOperator(a,b,b);
  }
   

  virtual ~StdTridiagonalOperator() {}
 
  using StdLinearOperator<Real>::apply;
  using StdLinearOperator<Real>::applyInverse;
  using StdLinearOperator<Real>::applyAdjoint;
  using StdLinearOperator<Real>::applyAdjointInverse;

 
  virtual void apply( vector<Real> &Hv, const vector<Real> &v, Real &tol ) const {
    Hv[0] = (*a_)[0]*v[0] + (*b_)[0]*v[1];

    for( int i=1; i<N_-1; ++i ) {
      Hv[i] = (*c_)[i-1]*v[i-1] + (*a_)[i]*v[i] + (*b_)[i]*v[i+1];
    }
          
    Hv[N_-1] = (*a_)[N_-1]*v[N_-1] + (*c_)[N_-2]*v[N_-2];

  }

  virtual void applyAdjoint( vector<Real> &Hv, const vector<Real> &v, Real &tol ) const {
    Hv[0] = (*a_)[0]*v[0] + (*c_)[0]*v[1];

    for( int i=1; i<N_-1; ++i ) {
      Hv[i] = (*b_)[i-1]*v[i-1] + (*a_)[i]*v[i] + (*c_)[i]*v[i+1];
    }
          
    Hv[N_-1] = (*a_)[N_-1]*v[N_-1] + (*b_)[N_-2]*v[N_-2];
  }

  virtual void applyInverse( vector<Real> &Hv, const vector<Real> &v, Real &tol ) const {
    Hv = v;
    copy();
    int INFO;
    lapack_.GTTRF(N_,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&INFO); 
    lapack_.GTTRS('N',N_,1,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&Hv[0],N_,&INFO);
}

  virtual void applyAdjointInverse( vector<Real> &Hv, const vector<Real> &v, Real &tol ) const {
    Hv = v;
    copy();
    int INFO;
    lapack_.GTTRF(N_,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&INFO); 
    lapack_.GTTRS('T',N_,1,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&Hv[0],N_,&INFO);
  }

}; // class StdTridiagonalOperator

} // namespace ROL

#endif
