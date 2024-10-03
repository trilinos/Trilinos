// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LANCZOS_H
#define ROL_LANCZOS_H

#include "ROL_Krylov.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_LAPACK.hpp"

namespace ROL {

/** \class ROL::Lanczos 
    \brief Interface for computing the Lanczos vectors
           and approximate solutions to symmetric indefinite
           linear systems

*/


template<class Real>
class Lanczos { 

  template <typename T> using ROL::Ptr = ROL::Ptr<T>;
  template <typename T> using vector = std::vector<T>;

  template typename vector<Real> size_type uint;

  typedef Vector<Real>         V;
  typedef LinearOperator<Real> OP;

  typedef ROL::ParameterList PL;

private:

  ROL::LAPACK<int,Real> lapack_;

  vector<ROL::Ptr<V> > Q_;     // Orthogonal basis
  vector<Real>    alpha_; // Diagonal recursion coefficients
  vector<Real>    beta_;  // Sub/super-diagonal recursion coefficients
  
  // Temporary vectors for factorizations, linear solves, and eigenvalue calculations
  vector<Real> dl_;    
  vector<Real> d_;
  vector<Real> du_;
  vector<Real> du2_;      
  vector<Real> y_;        // Arnoldi expansion coefficients

  vector<Real> work_;     // Scratch space for eigenvalue decomposition
  vector<int>  ipiv_;     // Pivots for LU

  ROL::Ptr<V> u_;              // An auxilliary vector

  Real max_beta_;          // maximum beta encountered
  Real tol_beta_;          // relative smallest beta allowed

  Real tol_ortho_;        // Maximum orthogonality loss tolerance

  int maxit_;             // Maximum number of vectors to store

  int k_;                 // current iterate number


  // Allocte memory for Arnoldi vectors and recurions coefficients
  void allocate( void ) {

    u_ = b.clone();

    alpha_.reserve(maxit_);
    beta_.reserve(maxit_);

    dl_.reserve(maxit_);
    d_.reserve(maxit_);
    du_.reserve(maxit_);
    du2_.reserve(maxit_);

    work_.reserve(4*maxit_);

    ipiv_.reserve(maxit_);
 
    y_.reserve(maxit_);
 
    alpha_.reserve(maxit_);
    beta_.reserve(maxit_);

    for( uint j=0; j<maxit_; ++j ) {
        Q_.push_back(b.clone());
    }  
  }

public:

  enum class FLAG_ITERATE : unsigned {
    ITERATE_SUCCESS = 0,
    ITERATE_SMALL_BETA,   // Beta too small to continue
    ITERATE_MAX_REACHED,  // Reached maximum number of iterations
    ITERATE_ORTHO_TOL,    // Reached maximim orthogonality loss
    ITERATE_LAST
  };

  enum class FLAG_SOLVE : unsigned {
    SOLVE_SUCCESS = 0,
    SOLVE_ILLEGAL_VALUE,
    SOLVE_SINGULAR_U,
    SOLVE_LAST
  };
 

  Lanczos( ROL::ParameterList &PL ) {
    PL &krylovList = parlist.sublist("General").sublist("Krylov");
    PL &lanczosList = krylovList.sublist("Lanczos");

    Real tol_default = std::sqrt(ROL_EPSILON<Real>());

    maxit_     = krylovList_.get("Iteration Limit",10);
    tol_beta_  = lanczosList.get("Beta Relative Tolerance", tol_default);
    tol_ortho_ = lanczosList.get("Orthogonality Tolerance", tol_default); 

  }

  void initialize( const V& b ) {
    allocate();
    reset(b);
     
  }

  void initialize( const V &x0, const V &b, const LO &A, Real &tol ) {
    allocate();
    reset(x0,b,A,tol);

  }


  void reset( const V &b ) {
    k_ = 0;
    max_beta_ = 0;
    Q_[0]->set(b);
    beta_[0] = Q_[0]->norm();
    max_beta_ = std::max(max_beta_,beta_[0]);
    Q_[0]->scale(1.0/beta_[0]);
  }

  void reset( const V &x0, const V &b, const LO &A, Real &tol ) {
    k_ = 0;
    max_beta_ = 0;
    Q_[0]->set(b);
    A.apply(*u_,x0,tol);
    Q_[0]->axpy(-1.0,*u_);
    beta_[0] = Q_[0]->norm();
    max_beta_ = std::max(max_beta_,beta_[0]);
    Q_[0]->scale(1.0/beta_[0]);     
  }

  FLAG_ITERATE iterate( const OP &A, Real &tol ) {

    if( k_ == maxit_ ) {
      return ITERATE_MAX_REACHED;
    }

    A.apply(*u_,*(Q_[k]),tol);
    Real delta;

    if( k_>0 ) { 
      u_->axpy(-beta_[k],V_[k_-1]);
    }
    alpha_[k] = u_->dot(*(V_[k]));
    u_->axpy(alpha_[k],V_[k_]);

    if( k_>0 ) {
      delta = u_->dot(*(V_[k-1]));
      u_->axpy(-delta,*(V_[k-1]));  
    }
    delta = u_->dot(*(V_[k]));
    alpha_[k] += delta;

    if( k_ < maxit_-1 ) {
      u_->axpy(-delta,*(V_[k_]));
      beta_[k+1] = u_->norm();
      max_beta_ = std::max(max_beta_,beta_[k+1]);
      if(beta_[k+1] < tol_beta_*max_beta_) {
        return ITERATE_SMALL_BETA;
      }

      V_[k+1]->set(*u_);
      V_[k+1]->scale(1.0/beta_[k+1]); 

      // Check orthogonality 
      Real dotprod = V_[k+1]->dot(*(V_[0]));

      if( std::sqrt(dotprod) > tol_ortho_ ) {
        return ITERATE_ORTHO_TOL;
      }  
    }

    ++k_; 
    return ITERATE_SUCCESS;
  } 
  
  // Compute the eigenvalues of the tridiagonal matrix T
  void eigenvalues( std::vector<Real> &E ) {

    std::vector<Real> Z(1,0); // Don't compute eigenvectors

    int INFO;
    int LDZ = 0;
    const char COMPZ = 'N':
    d_  = alpha_;
    du_ = beta_;

    lapack_->STEQR(COMPZ,k_,&d_[0],&du_[0],&Z[0],LDZ,&work_[0],&INFO);
 
    if( INFO < 0 ) {
      return SOLVE_ILLEGAL_VALUE;
    } 
    else if( INFO > 0 ) {
      return SOLVE_SINGULAR_U;
    }
    
  } 

  FLAG_SOLVE solve( V &x, Real tau=0 ) {

    const char TRANS = 'N';
    const int NRHS = 1;
    int INFO;
    
    // Fill arrays
    for(uint j=0;j<k_;++j) { 
      d_[j]  = alpha_[j]+tau;
    }

    dl_ = beta_;
    du_ = beta_;
    du2_.assign(maxit_,0);

    // Do Tridiagonal LU factorization 
    lapack_->GTTRF(k_,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&INFO);

    // Solve the factorized system for Arnoldi expansion coefficients
    lapack_->GTTRS(TRANS,k_,1,&dl[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&y_[0],k_,&INFO);       
    
  }



}; // class LanczosFactorization
}  // namespace ROL 

#endif // ROL_LANCZOS_H
