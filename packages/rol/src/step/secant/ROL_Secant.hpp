// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

enum ESecantMode {
  SECANTMODE_FORWARD = 0,
  SECANTMODE_INVERSE,
  SECANTMODE_BOTH
};

template<class Real>
struct SecantState {
  Ptr<Vector<Real>>              iterate;
  std::vector<Ptr<Vector<Real>>> iterDiff; // Step Storage
  std::vector<Ptr<Vector<Real>>> gradDiff; // Gradient Storage
  std::vector<Real>              product;  // Step-Gradient Inner Product Storage
  std::vector<Real>              product2; // Step-Gradient Inner Product Storage
  int storage;                             // Storage Size
  int current;                             // Current Storage Size
  int iter;                                // Current Optimization Iteration
  ESecantMode mode;                        // Intended application mode

  SecantState(int M, ESecantMode sm) : storage(M), current(-1), iter(0), mode(sm) {}
};

template<class Real>
class Secant : public LinearOperator<Real> {
protected:

  const Ptr<SecantState<Real>> state_; // Secant State
  Ptr<Vector<Real>> y_;
  bool useDefaultScaling_;
  Real Bscaling_;

private:

  bool isInitialized_;

public:

  virtual ~Secant() {}

  // Constructor
  Secant( int M = 10, bool useDefaultScaling = true, Real Bscaling = Real(1), ESecantMode mode = SECANTMODE_BOTH )
    : state_(makePtr<SecantState<Real>>(M,mode)),
      useDefaultScaling_(useDefaultScaling), Bscaling_(Bscaling),
      isInitialized_(false) {}

  Ptr<SecantState<Real>>& get_state() { return state_; }
  const Ptr<SecantState<Real>>& get_state() const { return state_; }

  // Update Secant Approximation
  virtual void updateStorage( const Vector<Real> &x,  const Vector<Real> &grad,
                              const Vector<Real> &gp, const Vector<Real> &s,
                              const Real snorm,       const int iter ) {
    const Real one(1);
    if ( !isInitialized_ ) {
      state_->iterate = x.clone();
      y_              = grad.clone();
      isInitialized_  = true;
    }
    state_->iterate->set(x);
    state_->iter = iter;
    y_->set(grad);
    y_->axpy(-one,gp);

    //Real sy = s.dot(y_->dual());
    Real sy = s.apply(*y_);
    if (sy > ROL_EPSILON<Real>()*snorm*snorm) {
      if (state_->current < state_->storage-1) {
        state_->current++;                                // Increment Storage
        state_->iterDiff.push_back(s.clone());            // Create new memory
        state_->gradDiff.push_back(grad.clone());         // Create new memory
      }
      else {
        state_->iterDiff.push_back(state_->iterDiff[0]);  // Move first element to the last
        state_->gradDiff.push_back(state_->gradDiff[0]);  // Move first element to the last
        state_->iterDiff.erase(state_->iterDiff.begin()); // Remove first element of s list 
        state_->gradDiff.erase(state_->gradDiff.begin()); // Remove first element of y list
        state_->product.erase(state_->product.begin());   // Remove first element of rho list
      }
      state_->iterDiff[state_->current]->set(s);          // s=x_{k+1}-x_k
      state_->gradDiff[state_->current]->set(*y_);        // y=g_{k+1}-g_k
      state_->product.push_back(sy);                      // ys=1/rho  
    }
  }

  // Apply Secant Approximate Inverse Hessian
  virtual void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const = 0;

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v ) const {
    Hv.set(v.dual());
    if (useDefaultScaling_) {
      if (state_->iter != 0 && state_->current != -1) {
        Real yy = state_->gradDiff[state_->current]->dot(*(state_->gradDiff[state_->current]));
        Hv.scale(state_->product[state_->current]/yy);
      }
    }
    else {
      Hv.scale(static_cast<Real>(1)/Bscaling_);
    }
  }

  // Apply Secant Approximate Hessian
  virtual void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const = 0;

  // Apply Initial Secant Approximate Hessian 
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v ) const {
    Bv.set(v.dual());
    if (useDefaultScaling_) {
      if (state_->iter != 0 && state_->current != -1) {
        Real yy = state_->gradDiff[state_->current]->dot(*(state_->gradDiff[state_->current]));
        Bv.scale(yy/state_->product[state_->current]);
      }
    }
    else {
      Bv.scale(Bscaling_);
    }
  }

  // Test Secant Approximations 
  void test(std::ostream &stream = std::cout ) const {
    if (isInitialized_) {
      Ptr<Vector<Real>> v  = state_->iterate->clone();
      Ptr<Vector<Real>> Hv = state_->iterate->clone();
      Ptr<Vector<Real>> Bv = state_->iterate->dual().clone();
      const Real one(1);
  
      // Print BHv -> Should be v
      v->randomize(-one,one);
      applyH(*Hv,*v);
      applyB(*Bv,*Hv);
      v->axpy(-one,*Bv);
      stream << " ||BHv-v|| = " << v->norm() << std::endl;
  
      // Print HBv -> Should be v
      v->randomize(-one,one);
      applyB(*Bv,*v);
      applyH(*Hv,*Bv);
      v->axpy(-one,*Hv);
      stream << " ||HBv-v|| = " << v->norm() << std::endl;
    }
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
