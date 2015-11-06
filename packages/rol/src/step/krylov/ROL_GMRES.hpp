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

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

namespace ROL {

template<class Real>
class GMRES : public Krylov<Real> {

  typedef Teuchos::SerialDenseMatrix<int, Real> SDMatrix;
  typedef Teuchos::SerialDenseVector<int, Real> SDVector;

  typedef typename std::vector<Real>::size_type uint;  

private:

 
  Teuchos::RCP<Vector<Real> > r_;
  Teuchos::RCP<Vector<Real> > z_;
  Teuchos::RCP<Vector<Real> > w_;
   
  Teuchos::RCP<SDMatrix> H_;      // quasi-Hessenberg matrix
  Teuchos::RCP<SDVector> cs_;     // Givens Rotations cosine components
  Teuchos::RCP<SDVector> sn_;     // Givens Rotations sine components
  Teuchos::RCP<SDVector> s_;      
  Teuchos::RCP<SDVector> y_;      
  Teuchos::RCP<SDVector> cnorm_;   

  Teuchos::RCP<std::vector<Real> > res_;
  
  bool isInitialized_;
  bool useInexact_;
  bool useInitialGuess_;
  int maxit_; 
  Real absTol_;
  Real relTol_;
 
  Teuchos::LAPACK<int,Real> lapack_;

public:
  
  GMRES( Teuchos::ParameterList &parlist ) : isInitialized_(false) {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::vector; 

    Teuchos::ParameterList &klist = parlist.sublist("General").sublist("Krylov");
    
    maxit_           = klist.get("Iteration Limit",50);
    absTol_          = klist.get("Absolute Tolerance", 1.e-4);
    relTol_          = klist.get("Relative Tolerance", 1.e-2);
    useInitialGuess_ = klist.get("Use Initial Guess",false);


    H_     = rcp( new SDMatrix( maxit_+1, maxit_ ) );
    cs_    = rcp( new SDVector( maxit_ ) );
    sn_    = rcp( new SDVector( maxit_ ) );
    s_     = rcp( new SDVector( maxit_+1 ) ); 
    y_     = rcp( new SDVector( maxit_+1 ) );
    cnorm_ = rcp( new SDVector( maxit_ ) );   
    res_   = rcp( new std::vector<Real>(maxit_+1,0.0) );
       
  }
 
  void run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b,
            LinearOperator<Real> &M, int &iter, int &flag ) {

    using Teuchos::RCP;
 
    flag = 0; 

    Real zero = 0.0;
    Real one =  1.0;

    if ( !isInitialized_ ) {
      r_  = b.clone();
      w_  = b.clone();
      z_  = x.clone();

      isInitialized_ = true;
    }

    Real itol  = std::sqrt(ROL_EPSILON); 

    // Compute initial residual
    if(useInitialGuess_) {
    
      A.apply(*r_,x,itol);
      r_->scale(-1.0);
      r_->plus(b);       // r = b-Ax
 
    }
    else {
      x.zero();
      r_->set(b);
    }

    Real temp  = 0;

    std::vector<RCP<Vector<Real > > > V;
    std::vector<RCP<Vector<Real > > > Z;

    (*res_)[0] = r_->norm();
     
    Real rtol  = std::min(absTol_,relTol_,(*res_)[0]);

    V.push_back(b.clone());
    (V[0])->set(*r_);
    (V[0])->scale(one/(*res_)[0]);    

    (*s_)(0) = (*res_)[0];

    uint m = static_cast<uint>(maxit_);

    for( uint i=0; i<m; ++i ) {

      if( useInexact_ ) {
        itol = rtol/(maxit_*(*res_)(i)); 
      }

      Z.push_back(x.clone());

      // Apply right preconditioner
      M.apply(*(Z[i]),*(V[i]),itol);

      // Apply operator
      A.apply(*w_,*(Z[i]),itol);

      // Evaluate coefficients and orthogonalize using Gram-Schmidt
      for( uint k=0; k<=i; ++k ) {
        (*H_)(k,i) = w_->dot(*(V[k]));
        w_->axpy( -(*H_)(k,i), *(V[k]) );
      } 
     
      (*H_)(i+1,i) = w_->norm();

      V.push_back( b.clone() );
      (V[i+1])->set(*w_);
      (V[i+1])->scale(one/((*H_)(i+1,i)));

      // Apply Givens rotations
      for( uint k=0; k<=i-1; ++k ) {
        temp         =  (*cs_)(k)*(*H_)(k,i) + (*sn_)(k)*(*H_)(k+1,i);
        (*H_)(k+1,i) = -(*sn_)(k)*(*H_)(k,i) + (*cs_)(k)*(*H_)(k+1,i); 
        (*H_)(k,i)   = temp;
      } 

      // Form i-th rotation matrix
      if( (*H_)(i+1,i) == zero ) {
        (*cs_)(i) = one;
        (*sn_)(i) = zero;
      }
      else if ( std::abs((*H_)(i+1,i)) > std::abs((*H_)(i,i)) ) { 
        temp = (*H_)(i,i) / (*H_)(i+1,i);
        (*sn_)(i) = one / std::sqrt( one + temp*temp );
        (*cs_)(i) = temp*(*sn_)(i); 
      }
      else {
        temp = (*H_)(i+1,i) / (*H_)(i,i);
        (*cs_)(i) = one / std::sqrt( one + temp*temp );
        (*sn_)(i) = temp*(cs_)(i);  
      }
     
      // Approximate residual norm
      temp         = (*cs_)(i)*(*sn_)(i);
      (*s_)(i+1)   = -(*sn_)(i)*(*s_)(i);
      (*s_)(i)     = temp;
      (*H_)(i,i)   = (*cs_)(i)*(*H_)(i,i) + (*sn_)*(*H_)(i+1,i);
      (*H_)(i+1,i) = zero;
      (*res_)[i+1] = std::abs((*s_)(i+1));
  
      if( (*res_)[i+1] < rtol ) {
        break;
      }

      // Update solution approximation.
      const char uplo = 'U';
      const char trans = 'N';
      const char diag = 'N';
      const char normin = 'N';
      Real scaling = zero;
      int info = 0;
      *y_ = *s_;
      lapack_.LATRS(uplo, trans, diag, normin, i+1, H_->values(), m+1, y_->values(), &scaling, cnorm_->values(), &info);

      z_->zero();

      for( uint k=0; k<=i;++k ) {
        z_->axpy((*y_)(k),*(Z[k]));
      }

      if( (*res_)[i+1] <= rtol ) {
        // Update solution vector
        x.plus(*z_);  
        break;
      }

      if(i == maxit_) {
        flag = 1;
      }
      iter = static_cast<int>(i);

    } // loop over i

  }  


}; // class GMRES

} // namespace ROL

#endif // ROL_GMRES_H

