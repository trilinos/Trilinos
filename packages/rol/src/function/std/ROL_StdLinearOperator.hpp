// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDLINEAROPERATOR_H
#define ROL_STDLINEAROPERATOR_H

#include "ROL_LinearOperator.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_LAPACK.hpp"


/** @ingroup func_group
    \class ROL::StdLinearOperator
    \brief Provides the std::vector implementation to apply a linear operator,
      which is a std::vector representation of column-stacked matrix

     Currently, this interface requires that the underlying matrix be square
    ---
*/


namespace ROL {

template <class Real>
class StdLinearOperator : public LinearOperator<Real> {
 
  typedef StdVector<Real> SV;

  typedef std::vector<Real> vector;

private:

  ROL::Ptr<std::vector<Real> > A_;
  int N_;
  int INFO_;
  
  mutable vector           PLU_;
  mutable std::vector<int> ipiv_; 

  ROL::LAPACK<int,Real>  lapack_;

public:

  StdLinearOperator() {}

  StdLinearOperator( ROL::Ptr<std::vector<Real> > &A ) : A_(A) { 
    int N2 = A_->size();
    N_ = (std::round(std::sqrt(N2)));
    bool isSquare = N_*N_ == N2;
    ROL_TEST_FOR_EXCEPTION( !isSquare, std::invalid_argument,
      "Error: vector representation of matrix must have a square "
      "number of elements.");
    ipiv_.resize(N_);   
  }

  virtual ~StdLinearOperator() {}
  
  using LinearOperator<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();
    update(*xp,flag,iter);   
  }

  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}

  // Matrix multiplication
  using LinearOperator<Real>::apply;
  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
        
    ROL::Ptr<vector> Hvp = dynamic_cast<SV&>(Hv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    apply(*Hvp,*vp,tol);
  }

  virtual void apply( std::vector<Real> &Hv, const std::vector<Real> &v, Real &tol ) const {
    for( int i=0; i<N_; ++i ) {
      Hv[i] = Real(0);
      for( int j=0; j<N_; ++j ) {
        Hv.at(i) += A_->at(N_*j+i)*v.at(j);
      }
    }
  }

  // Matrix multiplication with transpose
  using LinearOperator<Real>::applyAdjoint;
  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
        
    ROL::Ptr<vector> Hvp = dynamic_cast<SV&>(Hv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    applyAdjoint(*Hvp,*vp,tol);
  }

  virtual void applyAdjoint( std::vector<Real> &Hv, const std::vector<Real> &v, Real &tol ) const {
    for( int i=0; i<N_; ++i ) {
      Hv[i] = Real(0);
      for( int j=0; j<N_; ++j ) {
        Hv.at(i) += A_->at(N_*i+j)*v.at(j);
      }
    }
  }
  

  // Solve the system

  using LinearOperator<Real>::applyInverse;
  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    
    ROL::Ptr<vector> Hvp = dynamic_cast<SV&>(Hv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    applyInverse(*Hvp,*vp,tol);
  }

  virtual void applyInverse( std::vector<Real> &Hv, const std::vector<Real> &v, Real &tol ) const {

    const int LDA = N_;
    const int LDB = N_;
    int INFO;
    int NRHS = 1;
 
    Hv = v;
    PLU_  = *A_;

    // Do LU factorization
    lapack_.GETRF(N_,N_,&PLU_[0],LDA,&ipiv_[0],&INFO);

    ROL_TEST_FOR_EXCEPTION(INFO>0,std::logic_error,"Error in StdLinearOperator::applyInverse(): "
      "Zero diagonal element encountered in matrix factor U(" << INFO << "," << INFO << ").");

    ROL_TEST_FOR_EXCEPTION(INFO<0,std::logic_error,"Error in StdLinearOperator::applyInverse(): "
      "Illegal value encountered in element " << -INFO << " when performing LU factorization.");    

    // Solve factored system
    lapack_.GETRS('N',N_,NRHS,&PLU_[0],LDA,&ipiv_[0],&Hv[0],LDB,&INFO);

    ROL_TEST_FOR_EXCEPTION(INFO<0,std::logic_error,"Error in StdLinearOperator::applyInverse(): "
      "Illegal value encountered in element " << -INFO << " when solving the factorized system. "); 
  
  }

  // Solve the system with transposed matrix

  using LinearOperator<Real>::applyAdjointInverse;
  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    
    ROL::Ptr<vector> Hvp = dynamic_cast<SV&>(Hv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    applyAdjointInverse(*Hvp,*vp,tol);
  }

  virtual void applyAdjointInverse( std::vector<Real> &Hv, const std::vector<Real> &v, Real &tol ) const {

    const int LDA = N_;
    const int LDB = N_;
    int INFO;
    int NRHS = 1;
 
    Hv = v;
    PLU_  = *A_;

    // Do LU factorization
    lapack_.GETRF(N_,N_,&PLU_[0],LDA,&ipiv_[0],&INFO);

    ROL_TEST_FOR_EXCEPTION(INFO>0,std::logic_error,"Error in StdLinearOperator::applyAdjointInverse(): "
      "Zero diagonal element encountered in matrix factor U(" << INFO << "," << INFO << ").");

    ROL_TEST_FOR_EXCEPTION(INFO<0,std::logic_error,"Error in StdLinearOperator::applyAdjointInverse(): "
      "Illegal value encountered in element " << -INFO << " when performing LU factorization.");    

    // Solve factored system
    lapack_.GETRS('T',N_,NRHS,&PLU_[0],LDA,&ipiv_[0],&Hv[0],LDB,&INFO);

    ROL_TEST_FOR_EXCEPTION(INFO<0,std::logic_error,"Error in StdLinearOperator::applyAdjointInverse(): "
      "Illegal value encountered in element " << -INFO << " when solving the factorized system. "); 
  
  }

}; // class LinearOperator

} // namespace ROL

#endif
