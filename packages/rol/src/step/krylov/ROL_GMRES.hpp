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

namespace ROL {

template<class Real>
class GMRES : public Krylov<Real> {

  typedef Teuchos::SerialDenseMatrix<int, Real> SDMatrix;
  typedef Teuchos::SerialDenseVector<int, Real> SDVector;

  typedef typename std::vector<Real>::size_type size_type;  

private:

 
  Teuchos::RCP<Vector<Real> > r_;
  Teuchos::RCP<Vector<Real> > z_;
  Teuchos::RCP<Vector<Real> > w_;

  Teuchos::RCP<std::vector<Teuchos::RCP<Vector<Real> > > > V_;
  Teuchos::RCP<std::vector<Teuchos::RCP<Vector<Real> > > > Z_;

  Teuchos::RCP<SDMatrix> H_;      // quasi-Hessenberg matrix
  Teuchos::RCP<SDVector> cs_;     // Givens Rotations cosine components
  Teuchos::RCP<SDVector> sn_;     // Givens Rotations sine components
  Teuchos::RCP<SDVector> s_;      
  Teuchos::RCP<SDVector> y_;      
  Teuchos::RCP<SDVector> cnorm_;   

  
  bool isInitialized_;
  int maxit_; 
  Real absTol_;
  Real relTol_;
 
  Teuchos::LAPACK<int,Real> lapack_;

public:
  
  GMRES( Teuchos::ParameterList &parlist ) : isInitialized_(false) {

    using Teuchos::rcp; 

    Teuchos::ParameterList &klist = parlist.sublist("General").sublist("Krylov");
    
    maxit_  = klist.get("Iteration Limit",50);
    absTol_ = klist.get("Absolute Tolerance", 1.e-4);
    relTol_ = klist.get("Relative Tolerance", 1.e-2);

    H_     = rcp( new SDMatrix( maxit_+1, maxit_ ) );
    cs_    = rcp( new SDVector( maxit_ ) );
    sn_    = rcp( new SDVector( maxit_ ) );
    s_     = rcp( new SDVector( maxit_+1 ) ); 
    y_     = rcp( new SDVector( maxit_+1 ) );
    cnorm_ = rcp( new SDVector( maxit_ ) );   
           
  }
 
  void run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
            int &iter, int &flag ) {

    if ( !isInitialized_ ) {
      r_  = b.clone();
      w_  = b.clone();
      z_  = x.clone();
      isInitialized_ = true;
    }

    
    
  }  


}; // class GMRES

} // namespace ROL

#endif // ROL_GMRES_H

