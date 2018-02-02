
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

#pragma once

#include "XROL_Vector.hpp"

namespace XROL {

namespace Teuchos { class ParameterList; }
template<class> LinearOperator;

namespace Krylov {

template<class X, class Y=X>
class ConjugateGradients : public Solver<ConjugateGradients<X,Y>,X,Y> {
public: 

  ConjugateGradients( Teuchos::ParameterList& parlist ) : iter_(0), isInitialized_(false) {

    auto &cglist = parlist.sublist("Krylov").sublist("Conjugate Gradients");
    zeroInitialGuess_ = cglist.get("Zero Initial Guess", true);
    maxit_ = cglist.get("Maximum Iterations",10);

  }

  std::unique_ptr<Output> run( X& x, LinearOperator<X,Y> &A, const Y& b ) {

    static_assert(is_same<X,Y>::value,"Preconditioner is required when x and b "
                  "are elements of different vector spaces.");

    initialize(x,b);

    auto rho0 = dot(*r_,*r_);
   

    set(*p_,*r_);

    // Loop until maxit, negative dot product, or convergence
    for( iter_ = 0; iter_ < maxit_; ++iter ) {
      
      A.apply(*q_,*p_);
      auto alpha = rho0/dot(*p_,*q_);
      axpy(*x_, alpha,*p_);
      axpy(*x_,-alpha,*q_);
   

    }



    }

  }

  std::unique_ptr<Output> run( X& x, LinearOperator<X,Y> &A, 
                               const Y& b, LinearOperator<Y,X> &M ) {
    initialize(x,b);
    
    // Loop until maxit, negative dot product, or convergence
    for( iter_ = 0; iter_ < maxit_; ++iter ) {


    }


  }
  


private:

  void initialize(const X& x, const Y& b) {
    iter_ = 0;
   
    if( !isInitialized_ ) {
      v_ = clone(x);
      p_ = clone(x);   
      r_ = clone(b);   
      q_ = clone(b);
      isInitialized_ = true;
    }
 
    if( zeroInitialGuess_ ) {
      fill(x,0.0);
      set(*r_,b);
      A.apply(*v_,x);
      axpy(*r_,-1.0,*v);  // r = b-Ax
    }
    else {
      set(*r_,b);  // r = b
    }
  }

  std::unique_ptr<X> x_;
  std::unique_ptr<X> p_;
  std::unique_ptr<Y> q_; 
  std::unique_ptr<Y> r_;

  std::size_t iter_;
  std::size_t maxit_;
 
  bool isInitialized_;
  bool zeroInitialGuess_;
  magnitude_t<X> tol_;


};

} // namespace Krylov

} // namespace XROL

