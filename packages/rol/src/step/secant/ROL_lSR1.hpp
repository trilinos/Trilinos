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

#ifndef ROL_LSR1_H
#define ROL_LSR1_H

/** \class ROL::lSR1
    \brief Provides definitions for limited-memory SR1 operators.
*/

#include "ROL_Secant.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class lSR1 : public Secant<Real> {
private:

  mutable bool updateIterate_;
  bool isInitialized_;

public:
  lSR1(int M) : Secant<Real>(M), isInitialized_(false) {
    updateIterate_ = true;
  }

  // Update Secant Approximation
  void updateStorage( const Vector<Real> &x,  const Vector<Real> &grad,
                      const Vector<Real> &gp, const Vector<Real> &s,
                      const Real snorm,       const int iter ) {
    Real one(1);
    // Get Generic Secant State
    ROL::Ptr<SecantState<Real> >& state = Secant<Real>::get_state();
    if ( !isInitialized_ ) {
      state->iterate = x.clone();
      isInitialized_ = true;
    }

    state->iterate->set(x);
    state->iter = iter;
    ROL::Ptr<Vector<Real> > gradDiff = grad.clone();
    gradDiff->set(grad);
    gradDiff->axpy(-one,gp);

    Real sy = s.dot(gradDiff->dual());
    if (updateIterate_ || state->current == -1) {
      if (state->current < state->storage-1) {
        state->current++;                               // Increment Storage
      }
      else {
        state->iterDiff.erase(state->iterDiff.begin()); // Remove first element of s list 
        state->gradDiff.erase(state->gradDiff.begin()); // Remove first element of y list
        state->product.erase(state->product.begin());   // Remove first element of rho list
      }
      state->iterDiff.push_back(s.clone());
      state->iterDiff[state->current]->set(s);          // s=x_{k+1}-x_k
      state->gradDiff.push_back(grad.clone());
      state->gradDiff[state->current]->set(*gradDiff);  // y=g_{k+1}-g_k
      state->product.push_back(sy);                     // ys=1/rho  
    }
    updateIterate_ = true;
  }

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v ) const {
    Hv.set(v.dual());
  }

  // Apply lSR1 Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const {
    // Get Generic Secant State
    const ROL::Ptr<SecantState<Real> >& state = Secant<Real>::get_state();

    // Apply initial Hessian approximation to v
    applyH0(Hv,v);

    std::vector<ROL::Ptr<Vector<Real> > > a(state->current+1);
    std::vector<ROL::Ptr<Vector<Real> > > b(state->current+1);
    Real byi(0), byj(0), bv(0), normbi(0), normyi(0), one(1);
    for (int i = 0; i <= state->current; i++) {
      // Compute Hy
      a[i] = Hv.clone();
      applyH0(*(a[i]),*(state->gradDiff[i]));
      for (int j = 0; j < i; j++) {
        byj = b[j]->dot((state->gradDiff[j])->dual());
        byi = b[j]->dot((state->gradDiff[i])->dual());
        a[i]->axpy(byi/byj,*(b[j]));
      }
      // Compute s - Hy
      b[i] = Hv.clone();
      b[i]->set(*(state->iterDiff[i]));
      b[i]->axpy(-one,*(a[i]));

      // Compute Hv
      byi    = b[i]->dot((state->gradDiff[i])->dual());
      normbi = b[i]->norm();
      normyi = (state->gradDiff[i])->norm();
      if ( i == state->current && std::abs(byi) < sqrt(ROL_EPSILON<Real>())*normbi*normyi ) {
        updateIterate_ = false;
      }
      else {
        updateIterate_ = true;
        bv  = b[i]->dot(v.dual());
        Hv.axpy(bv/byi,*(b[i]));
      }
    }
  }

  // Apply Initial Secant Approximate Hessian
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v ) const {
    Bv.set(v.dual());
  }

  // Apply lSR1 Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const {
    // Get Generic Secant State
    const ROL::Ptr<SecantState<Real> >& state = Secant<Real>::get_state();

    // Apply initial Hessian approximation to v
    applyB0(Bv,v);

    std::vector<ROL::Ptr<Vector<Real> > > a(state->current+1);
    std::vector<ROL::Ptr<Vector<Real> > > b(state->current+1);
    Real bsi(0), bsj(0), bv(0), normbi(0), normsi(0), one(1);
    for (int i = 0; i <= state->current; i++) {
      // Compute Hy
      a[i] = Bv.clone();
      applyB0(*(a[i]),*(state->iterDiff[i]));
      for (int j = 0; j < i; j++) {
        bsj = (state->iterDiff[j])->dot(b[j]->dual());
        bsi = (state->iterDiff[i])->dot(b[j]->dual());
        a[i]->axpy(bsi/bsj,*(b[j]));
      }
      // Compute s - Hy
      b[i] = Bv.clone();
      b[i]->set(*(state->gradDiff[i]));
      b[i]->axpy(-one,*(a[i]));

      // Compute Hv
      bsi    = (state->iterDiff[i])->dot(b[i]->dual());
      normbi = b[i]->norm();
      normsi = (state->iterDiff[i])->norm();
      if ( i == state->current && std::abs(bsi) < sqrt(ROL_EPSILON<Real>())*normbi*normsi ) {
        updateIterate_ = false;
      }
      else {
        updateIterate_ = true;
        bv  = b[i]->dot(v.dual());
        Bv.axpy(bv/bsi,*(b[i]));
      }
    }
  }
};

}

#endif
