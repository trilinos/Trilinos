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

#include "ROL_ParameterList.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
struct SecantState {
  ROL::Ptr<Vector<Real> >               iterate;
  std::vector<ROL::Ptr<Vector<Real> > > iterDiff; // Step Storage
  std::vector<ROL::Ptr<Vector<Real> > > gradDiff; // Gradient Storage
  std::vector<Real>                         product;  // Step-Gradient Inner Product Storage
  std::vector<Real>                         product2; // Step-Gradient Inner Product Storage
  int storage;                                        // Storage Size
  int current;                                        // Current Storage Size
  int iter;                                           // Current Optimization Iteration
};

template<class Real>
class Secant : public LinearOperator<Real> {
private:

  ROL::Ptr<SecantState<Real> > state_; // Secant State
  bool isInitialized_;

public:

  virtual ~Secant() {}

  // Constructor
  Secant( int M = 10 ) : isInitialized_(false) {
    state_ = ROL::makePtr<SecantState<Real>>(); 
    state_->storage = M;
    state_->current = -1;
    state_->iter    = 0;
  }

  ROL::Ptr<SecantState<Real> >& get_state() { return state_; }
  const ROL::Ptr<SecantState<Real> >& get_state() const { return state_; }

  // Update Secant Approximation
  virtual void updateStorage( const Vector<Real> &x,  const Vector<Real> &grad,
                              const Vector<Real> &gp, const Vector<Real> &s,
                              const Real snorm,       const int iter ) {
    Real one(1);
    if ( !isInitialized_ ) {
      state_->iterate = x.clone();
      isInitialized_ = true;
    }
    state_->iterate->set(x);
    state_->iter = iter;
    ROL::Ptr<Vector<Real> > gradDiff = grad.clone();
    gradDiff->set(grad);
    gradDiff->axpy(-one,gp);

    Real sy = s.dot(gradDiff->dual());
    if (sy > ROL_EPSILON<Real>()*snorm*snorm) {
      if (state_->current < state_->storage-1) {
        state_->current++;                                // Increment Storage
      }
      else {
        state_->iterDiff.erase(state_->iterDiff.begin()); // Remove first element of s list 
        state_->gradDiff.erase(state_->gradDiff.begin()); // Remove first element of y list
        state_->product.erase(state_->product.begin());   // Remove first element of rho list
      }
      state_->iterDiff.push_back(s.clone()); 
      state_->iterDiff[state_->current]->set(s);          // s=x_{k+1}-x_k
      state_->gradDiff.push_back(grad.clone()); 
      state_->gradDiff[state_->current]->set(*gradDiff);  // y=g_{k+1}-g_k
      state_->product.push_back(sy);                      // ys=1/rho  
    }
  }

  // Apply Secant Approximate Inverse Hessian
  virtual void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const = 0;

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v ) const {
    Hv.set(v.dual());
    if (state_->iter != 0 && state_->current != -1) {
      Real yy = state_->gradDiff[state_->current]->dot(*(state_->gradDiff[state_->current]));
      Hv.scale(state_->product[state_->current]/yy);
    }
  }

  // Apply Secant Approximate Hessian
  virtual void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const = 0;

  // Apply Initial Secant Approximate Hessian 
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v ) const {
    Bv.set(v.dual());
    if (state_->iter != 0 && state_->current != -1) {
      Real yy = state_->gradDiff[state_->current]->dot(*(state_->gradDiff[state_->current]));
      Bv.scale(yy/state_->product[state_->current]);
    }
  }

  // Test Secant Approximations 
  void test( const Vector<Real> &x, const Vector<Real> &s ) const {
    ROL::Ptr<Vector<Real> > vec  = x.clone();
    ROL::Ptr<Vector<Real> > Hvec = x.clone();
    ROL::Ptr<Vector<Real> > Bvec = x.clone();
    Real one(1);
  
    // Print BHv -> Should be v
    vec->set(s);
    applyH(*Hvec,*vec);
    applyB(*Bvec,*Hvec);
    vec->axpy(-one,*Bvec);
    std::cout << " ||BHv-v|| = " << vec->norm() << "\n";
  
    // Print HBv -> Should be v
    vec->set(s);
    applyB(*Bvec,*vec);
    applyH(*Hvec,*Bvec);
    vec->axpy(-one,*Hvec);
    std::cout << " ||HBv-v|| = " << vec->norm() << "\n";
  }

  void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
    applyB(Hv,v);
  }

  void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
    applyH(Hv,v);
  }

};

}

#include "ROL_SecantFactory.hpp"

#endif
