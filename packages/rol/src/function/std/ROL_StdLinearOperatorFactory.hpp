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


#ifndef ROL_STDLINEAROPERATORFACTORY_H
#define ROL_STDLINEAROPERATORFACTORY_H

#include "ROL_StdLinearOperator.hpp"

#include <ctime>
#include <cstdlib>
#include <string>

/** @ingroup func_group
    \class ROL::StdLinearOperatorFactory
    \brief Creates StdLinearOperator objects which wrap random
           matrices of the desired size and property
   ---
*/

namespace ROL {

template<class Real> 
class StdLinearOperatorFactory {

  template <typename T> using RCP = Teuchos::RCP<T>;

  typedef LinearOperator<Real>    OP;
  typedef StdLinearOperator<Real> StdOP;

  typedef std::vector<Real>       vector;

private:


  Teuchos::BLAS<int,Real>     blas_;
  Teuchos::LAPACK<int,Real>   lapack_;

  // Fill x with uniformly-distributed random values from [lower,upper]
  void randomize( vector &x, Real lower=0.0, Real upper=1.0 ) {
    int N = x.size(); 
    for( int i=0; i<N; ++i ) {
      x[i] = lower+(upper-lower)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    }
  }

  void diagonal( vector &D, const vector& d ) {
    int N = d.size(); 
    int N2 = N*N;
    D.reserve(N2);
    D.assign(N2,0.0);
    for( int i=0; i<N; ++i ) {
      D[(N+1)*i] = d[i];
    }
  }

  // C = A*B with optional transposes
  void multiply( vector &C, const vector &A, const vector &B, bool transA=false, bool transB=false ) {
    int N2 = A.size();
    int N = (std::round(std::sqrt(N2)));
    bool isSquare = N*N == N2;
    TEUCHOS_TEST_FOR_EXCEPTION( !isSquare, std::invalid_argument,
      "Error: vector representation A of matrix must have a square "
      "number of elements."); 
    TEUCHOS_TEST_FOR_EXCEPTION( B.size() != N2, std::invalid_argument,
      "Error: vectors A and B must have the same length.");
    TEUCHOS_TEST_FOR_EXCEPTION( C.size() != N2, std::invalid_argument,
      "Error: vectors A and C must have the same length.");

    char tra = transA : 'T' : 'N';
    char trb = transB : 'T' : 'N';

    blas_.GEMM(tra,trb,N,N,N,1.0,&A[0],&B[0],N,0.0,N);

  }

  // Orthogonalize the columns of A
  void orthogonalize( vector &A ) {
    int N2 = A.size();
    int N = (std::round(std::sqrt(N2)));
    bool isSquare = N*N == N2;
    TEUCHOS_TEST_FOR_EXCEPTION( !isSquare, std::invalid_argument,
      "Error: vector representation of matrix must have a square "
      "number of elements."); 

    vector TAU(N,0.0);
    
    int LDA = N;
    Real LWORK1, LWORK2; 
    int INFO = 0;

    // Query workspace 
    lapack_.GEQRF(N,N,&A[0],LDA,&TAU[0],&LWORK1,-1,&INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(INFO,std::logic_error,"LAPACK GEQRF LWORK query failed.");
    
    lapack_.ORGQR(N,N,N,&A[0],LDA,&TAU[0],&LWORK2,-1,&INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(INFO,std::logic_error,"LAPACK ORGQR LWORK query failed.");

    const int LWORK = std::max(std::abs(LWORK1),std::abs(LWORK2));

    vector WORK(LWORK);
 
    // Factor the input matrix
    lapack_.GEQRF(N,N,&A[0],LDA,&TAU[0],LWORK,&INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(INFO,std::logic_error,"LAPACK GEQRF failed with INFO = " << INFO );
    
    // Overwrite the input matrix with the orthogonal matrix Q
    lapack_.ORGQR(N,N,N,&A[0],LDA,&TAU[0],&WORK[0],LWORK,&INFO);    
    TEUCHOS_TEST_FOR_EXCEPTION(INFO,std::logic_error,"LAPACK ORGQR failed with INFO = " << INFO );

  }  


public:

  enum class EMatrixType unsigned {
    MATRIX_SPD = 0,    // Symmetric Positive Definite matrix
    MATRIX_SYMMETRIC,  // Symmetric Indefinite matrix
    MATRIX_UNITARY,    // Unitary matrix
    MATRIX_SINGULAR_S, // Singular symmetric matrix
    MATRIX_SINGULAR_N, // Singular nonsymmetric matrix
    MATRIX_DEFAULT,    // Random nonsymmetric matrix
    MATRIX_LAST 
  }; 

  ~StdLinearOperatorFactor(void) {}

  EMatrixType StringToEMatrixType( satd::string s ) {
    s = removeStringFormat(s);
    for ( EMatrixType mt = MATRIX_SPD; mt < MATRIX_LAST; ++ ) {
      if ( !s.compare(removeStringFormat(EStepToString(mt))) ) {
        return mt;
      }
  }

  std::string EMatrixTypeToString( EMatrixType mt ) {
    std::string retString;
    switch( mt ) {
      case MATRIX_SPD:        retString = "Symmetric Positive Definite";  break;
      case MATRIX_SYMMETRIC:  retString = "Symmetric Indefinite";         break;
      case MATRIX_UNITARY:    retString = "Unitary";                      break;
      case MATRIX_SINGULAR_S: retString = "Singular Symmetric";           break;
      case MATRIX_SINGILAR_N: retString = "Singular Nonsymmetric";        break;
      case MATRIX_DEFAULT:    retString = "Default";                      break;
      case MATRIX_LAST:       retString = "Last (dummy)";                 break; 
    }
    return retString;    
  }

  Teuchos::RCP<LinearOperator<Real> > getOperator( int size, const std::string &type="" ) const {
    EMatrixType emt = StringToEMatrixType(type);

    using Teuchos::rcp;

    int n2 = size*size;

    RCP<vector> Ap = rcp( new vector(n2) );

    switch( emt ) {
      case MATRIX_SPD: {
        
        vector d(size);
        randomize(d,1.0,2.0);

        // A = D
        diagonal(*Ap,d);
        
        vector Q(n2);
        randomize(Q,-1.0,1.0);        
        orthogonalize(Q);
 
        // A = D*Q
        multiply(*Ap,*Ap,Q);
        
        // A = Q'*D*Q
        multiply(*Ap,Q,*Ap,true);

      }
      break;

      case MATRIX_SYMMETRIC: {
        vector d(size);
        randomize(d);

        // A = D
        diagonal(*Ap,d);
        
        vector Q(n2);
        randomize(Q,-1.0,1.0);        
        orthogonalize(Q);
 
        // A = D*Q
        multiply(*Ap,*Ap,Q);
        
        // A = Q'*D*Q
        multiply(*Ap,Q,*Ap,true);
      }
      break; 

      case MATRIX_UNITARY: {
        randomize(*Ap);
        orthogonalize(*Ap);  
      } 
  
      case MATRIX_SINGULAR_S: {
        vector d(size);
        randomize(d);
  
        d[0] = 0;

        // A = D
        diagonal(*Ap,d);
        
        vector Q(n2);
        randomize(Q,-1.0,1.0);        
        orthogonalize(Q);
 
        // A = D*Q
        multiply(*Ap,*Ap,Q);
        
        // A = Q'*D*Q
        multiply(*Ap,Q,*Ap,true);
      

      case MATRIX_SINGULAR_N: {

        vector d(size);
        randomize(d,0.0,1.0);
  
        d[0] = 0;

        // A = D
        diagonal(*Ap,d);

        vector V(n2);
        randomize(V,-1.0,1.0);        
        orthogonalize(V);

        // A = D*V'
        multiply(*Ap,*Ap,Q,false,true);
       
        vector U(n2);
        randomize(U,-1.0,1.0);        
        orthogonalize(U);

        // A = U*D*V'
        multiply(*Ap,U,*Ap);
 
      }

      case MATRIX_DEFAULT: 
      default: {
        randomize(*Ap);
      }   

    }
    return rcp( new StdOP(Ap) );
  }

  

}; // class StdLinearOperatorFactory

} // namespace ROL


#endif // ROL_STDLINEAROPERATORFACTORY_H
