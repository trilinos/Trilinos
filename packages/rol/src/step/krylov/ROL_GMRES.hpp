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

#ifndef ROL_GMRES_H
#define ROL_GMRES_H

/** \class ROL::GMRES
    \brief Preconditioned GMRES solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_LAPACK.hpp"
#include "ROL_LinearAlgebra.hpp"


namespace ROL {

template<class Real>
class GMRES : public Krylov<Real> {

  typedef LA::Matrix<Real> SDMatrix;
  typedef LA::Vector<Real> SDVector;

private:
 
  ROL::Ptr<Vector<Real> > r_;
  ROL::Ptr<Vector<Real> > z_;
  ROL::Ptr<Vector<Real> > w_;
   
  ROL::Ptr<SDMatrix> H_;      // quasi-Hessenberg matrix
  ROL::Ptr<SDVector> cs_;     // Givens Rotations cosine components
  ROL::Ptr<SDVector> sn_;     // Givens Rotations sine components
  ROL::Ptr<SDVector> s_;      
  ROL::Ptr<SDVector> y_;      
  ROL::Ptr<SDVector> cnorm_;   

  ROL::Ptr<std::vector<Real> > res_;
  
  bool isInitialized_;
  bool useInexact_;
  bool useInitialGuess_;    // If false, inital x will be ignored and zero vec used
 
  ROL::LAPACK<int,Real> lapack_;

public:
  
  GMRES( ROL::ParameterList &parlist ) : Krylov<Real>(parlist), isInitialized_(false) {

    
    
    using std::vector; 

    Real zero(0);

    ROL::ParameterList &gList = parlist.sublist("General");
    ROL::ParameterList &kList = gList.sublist("Krylov");
    
    useInexact_      = gList.get("Inexact Hessian-Times-A-Vector",false);
    useInitialGuess_ = kList.get("Use Initial Guess",false);
    int maxit = Krylov<Real>::getMaximumIteration();

    H_     = ROL::makePtr<SDMatrix>( maxit+1, maxit );
    cs_    = ROL::makePtr<SDVector>( maxit );
    sn_    = ROL::makePtr<SDVector>( maxit );
    s_     = ROL::makePtr<SDVector>( maxit+1 ); 
    y_     = ROL::makePtr<SDVector>( maxit+1 );
    cnorm_ = ROL::makePtr<SDVector>( maxit );   
    res_   = ROL::makePtr<std::vector<Real>>(maxit+1,zero);
       
  }
 
  Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b,
            LinearOperator<Real> &M, int &iter, int &flag ) {

    Real absTol = Krylov<Real>::getAbsoluteTolerance();
    Real relTol = Krylov<Real>::getRelativeTolerance();
    int maxit = Krylov<Real>::getMaximumIteration();
 
    flag = 0;

    Real zero(0), one(1);

    if ( !isInitialized_ ) {
      r_  = b.clone();
      w_  = b.clone();
      z_  = x.clone();

      isInitialized_ = true;
    }

    Real itol  = std::sqrt(ROL_EPSILON<Real>());

    // Compute initial residual
    if(useInitialGuess_) {

      A.apply(*r_,x,itol);
      r_->scale(-one);
      r_->plus(b);       // r = b-Ax
 
    }
    else {
      x.zero();
      r_->set(b);
    }

    Real temp  = 0;

    std::vector<ROL::Ptr<Vector<Real > > > V;
    std::vector<ROL::Ptr<Vector<Real > > > Z;

    (*res_)[0] = r_->norm();

    // This should be a tolerance check
    Real rtol = std::min(absTol,relTol*(*res_)[0]);
    if ((*res_)[0] <= rtol) {
      iter = 0;
      flag = 0;
      return (*res_)[0];
    }

    V.push_back(b.clone());
    (V[0])->set(*r_);
    (V[0])->scale(one/(*res_)[0]);

    (*s_)(0) = (*res_)[0];

    for( iter=0; iter<maxit; ++iter ) {

      // std::cout << "iter = " << iter << "  rnorm = " << (*res_)[iter] << "  xnorm = " << x.norm() << "  bnorm = " << b.norm() << std::endl;

      if( useInexact_ ) {
        itol = rtol/(maxit*(*res_)[iter]);
      }

      Z.push_back(x.clone());

      // Apply right preconditioner
      M.applyInverse(*(Z[iter]),*(V[iter]),itol);

      // Apply operator
      A.apply(*w_,*(Z[iter]),itol);

      // Evaluate coefficients and orthogonalize using Gram-Schmidt
      for( int k=0; k<=iter; ++k ) {
        (*H_)(k,iter) = w_->dot(*(V[k]));
        w_->axpy( -(*H_)(k,iter), *(V[k]) );
      }

      (*H_)(iter+1,iter) = w_->norm();

      V.push_back( b.clone() );
      (V[iter+1])->set(*w_);
      (V[iter+1])->scale(one/((*H_)(iter+1,iter)));

      // Apply Givens rotations
      for( int k=0; k<=iter-1; ++k ) {
        temp            =  (*cs_)(k)*(*H_)(k,iter) + (*sn_)(k)*(*H_)(k+1,iter);
        (*H_)(k+1,iter) = -(*sn_)(k)*(*H_)(k,iter) + (*cs_)(k)*(*H_)(k+1,iter); 
        (*H_)(k,iter)   = temp;
      }

      // Form i-th rotation matrix
      if( (*H_)(iter+1,iter) == zero ) {
        (*cs_)(iter) = one;
        (*sn_)(iter) = zero;
      }
      else if ( std::abs((*H_)(iter+1,iter)) > std::abs((*H_)(iter,iter)) ) {
        temp = (*H_)(iter,iter) / (*H_)(iter+1,iter);
        (*sn_)(iter) = one / std::sqrt( one + temp*temp );
        (*cs_)(iter) = temp*(*sn_)(iter);
      }
      else {
        temp = (*H_)(iter+1,iter) / (*H_)(iter,iter);
        (*cs_)(iter) = one / std::sqrt( one + temp*temp );
        (*sn_)(iter) = temp*(*cs_)(iter);
      }

      // Approximate residual norm
      temp               = (*cs_)(iter)*(*s_)(iter);
      (*s_)(iter+1)      = -(*sn_)(iter)*(*s_)(iter);
      (*s_)(iter)        = temp;
      (*H_)(iter,iter)   = (*cs_)(iter)*(*H_)(iter,iter) + (*sn_)(iter)*(*H_)(iter+1,iter);
      (*H_)(iter+1,iter) = zero;
      (*res_)[iter+1]    = std::abs((*s_)(iter+1));
 
      // Update solution approximation.
      const char uplo = 'U';
      const char trans = 'N';
      const char diag = 'N';
      const char normin = 'N';
      Real scaling = zero;
      int info = 0;
      *y_ = *s_;
      lapack_.LATRS(uplo, trans, diag, normin, iter+1, H_->values(), maxit+1, y_->values(), &scaling, cnorm_->values(), &info);

      z_->zero();

      for( int k=0; k<=iter;++k ) {
        z_->axpy((*y_)(k),*(Z[k]));
      }

      if( (*res_)[iter+1] <= rtol ) {
        // Update solution vector
        x.plus(*z_);
        break;
      }
    } // loop over iter

    if(iter == maxit) {
      flag = 1;
      x.plus(*z_);
      return (*res_)[iter];
    }

    return (*res_)[iter+1];
  }

}; // class GMRES

} // namespace ROL

#endif // ROL_GMRES_H

