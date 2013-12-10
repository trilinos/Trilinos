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

#ifndef ROL_SECANT_H
#define ROL_SECANT_H

/** \class ROL::Secant
    \brief Provides interface for and implements limited-memory secant operators.
*/

#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
struct SecantState {
  std::vector<Teuchos::RCP<Vector<Real> > > iterDiff; // Step Storage
  std::vector<Teuchos::RCP<Vector<Real> > > gradDiff; // Gradient Storage
  std::vector<Real>                         product;  // Step-Gradient Inner Product Storage
  std::vector<Real>                         product2; // Step-Gradient Inner Product Storage
  int storage;                                        // Storage Size
  int current;                                        // Current Storage Size
  int iter;                                           // Current Optimization Iteration
};

template<class Real>
class Secant {
private:

  Teuchos::RCP<SecantState<Real> > state_; // Secant State

public:

  virtual ~Secant() {}

  // Constructor
  Secant( int M = 10 ) {
    state_ = Teuchos::rcp( new SecantState<Real> ); 
    state_->storage = M;
    state_->current = -1;
    state_->iter    = 0;
  }

  Teuchos::RCP<SecantState<Real> >& get_state() { return this->state_; }

  // Update Secant Approximation
  virtual void update( const Vector<Real> &grad, const Vector<Real> &gp, const Vector<Real> &s, 
                       const Real snorm, const int iter ) {
    this->state_->iter = iter;
    Teuchos::RCP<Vector<Real> > gradDiff = grad.clone();
    gradDiff->set(grad);
    gradDiff->axpy(-1.0,gp);

    Real sy = s.dot(*gradDiff);
    if (sy > ROL_EPSILON*snorm*snorm) {
      if (this->state_->current < this->state_->storage-1) {
        this->state_->current++;                                      // Increment Storage
      }
      else {
        this->state_->iterDiff.erase(this->state_->iterDiff.begin()); // Remove first element of s list 
        this->state_->gradDiff.erase(this->state_->gradDiff.begin()); // Remove first element of y list
        this->state_->product.erase(this->state_->product.begin());   // Remove first element of rho list
      }
      this->state_->iterDiff.push_back(s.clone()); 
      this->state_->iterDiff[this->state_->current]->set(s);          // s=x_{k+1}-x_k
      this->state_->gradDiff.push_back(s.clone()); 
      this->state_->gradDiff[this->state_->current]->set(*gradDiff);  // y=g_{k+1}-g_k
      this->state_->product.push_back(sy);                            // ys=1/rho  
    }
  }

  // Apply Secant Approximate Inverse Hessian
  virtual void applyH( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x ) = 0;

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x ) {
    Hv.set(v);
    if (this->state_->iter != 0 && this->state_->current != -1) {
      Real yy = this->state_->gradDiff[this->state_->current]->dot(*(this->state_->gradDiff[this->state_->current]));
      Hv.scale(this->state_->product[this->state_->current]/yy);
    }
  }

  // Apply Secant Approximate Hessian
  virtual void applyB( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x ) = 0;

  // Apply Initial Secant Approximate Hessian 
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x ) {
    Bv.set(v);
    if (this->state_->iter != 0 && this->state_->current != -1) {
      Real yy = this->state_->gradDiff[this->state_->current]->dot(*(this->state_->gradDiff[this->state_->current]));
      Bv.scale(yy/this->state_->product[this->state_->current]);
    }
  }

  // Test Secant Approximations 
  void test( const Vector<Real> &x, const Vector<Real> &s ) {
    Teuchos::RCP<Vector<Real> > vec  = x.clone();
    Teuchos::RCP<Vector<Real> > Hvec = x.clone();
    Teuchos::RCP<Vector<Real> > Bvec = x.clone();
  
    // Print BHv -> Should be v
    vec->set(s);
    this->applyH(*Hvec,*vec,x);
    this->applyB(*Bvec,*Hvec,x);
    vec->axpy(-1.0,*Bvec);
    std::cout << " ||BHv-v|| = " << vec->norm() << "\n";
  
    // Print HBv -> Should be v
    vec->set(s);
    this->applyB(*Bvec,*vec,x);
    this->applyH(*Hvec,*Bvec,x);
    vec->axpy(-1.0,*Hvec);
    std::cout << " ||HBv-v|| = " << vec->norm() << "\n";
  }

};

}

#include "ROL_lBFGS.hpp"
#include "ROL_lDFP.hpp"
#include "ROL_lSR1.hpp"
#include "ROL_BarzilaiBorwein.hpp"

namespace ROL {
  template<class Real>
  inline Teuchos::RCP<Secant<Real> > getSecant( ESecant esec = SECANT_LBFGS, int L = 10, int BBtype = 1 ) {
    switch (esec) {
      case SECANT_LBFGS:           return Teuchos::rcp( new lBFGS<Real>(L) );
      case SECANT_LDFP:            return Teuchos::rcp( new lDFP<Real>(L) );
      case SECANT_LSR1:            return Teuchos::rcp( new lSR1<Real>(L) );
      case SECANT_BARZILAIBORWEIN: return Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
      default:                     return Teuchos::null; 
    }
  }
}

#endif
