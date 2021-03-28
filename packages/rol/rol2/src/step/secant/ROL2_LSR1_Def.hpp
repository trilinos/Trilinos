// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL2) Package
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
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
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
#ifndef ROL2_LSR1_DEF_HPP
#define ROL2_LSR1_DEF_HPP

namespace ROL2 {

// Update lSR1 Approximation
template<class Real>
void lSR1<Real>::updateStorage( const Vector<Real>& x,  
                                const Vector<Real>& grad,
                                const Vector<Real>& gp, 
                                const Vector<Real>& s,
                                      Real          snorm,      
                                      int           iter )  {

  auto& state = Secant<Real>::getState();
  const Real one(1);
  Real tol = default_tolerance<Real>();

  if ( !isInitialized_ ) {
    state.iterate_ = x.clone();
    y_              = grad.clone();
    if (state.mode_ == Secant<Real>::Mode::Forward) {
      Bs_ = grad.clone(); dual_ = grad.clone();
    }
    else if (state.mode_ == Secant<Real>::Mode::Inverse) {
      Hy_ = x.clone();    prim_ = x.clone();
    }
    else {
      Bs_ = grad.clone(); dual_ = grad.clone();
      Hy_ = x.clone();    prim_ = x.clone();
    }
    isInitialized_ = true;
  } // endif( !isInitialized_ )

  // Update iterate
  state.iter_ = iter;
  state.iterate_->set(x);

  // Compute gradient difference
  y_->set(grad);
  y_->axpy(-one,gp);

  Real dotF(ROL_INF<Real>), tolF(0), dotI(ROL_INF<Real>), tolI(0);

  if( non_inverse_ ) {
    // Compute y - Bs and <s, y - Bs>
    applyB(*Bs_,s);
    Bs_->scale(-one);
    Bs_->plus(*y_);
    dotF = s.apply(*Bs_);
    tolF = tol*snorm*Bs_->norm();
  }

  if( non_forward_ ) {
    // Compute s - Hy and <y, s - Hy>
    applyH(*Hy_,*y_);
    Hy_->scale(-one);
    Hy_->plus(s);
    dotI = y_->apply(*Hy_);
    tolI = tol*y_->norm()*Hy_->norm();
  }

  if( std::abs(dotF) > tolF && std::abs(dotI) > tolI ) {

    if (state.current_ < state.storage_-1) {
      state.current_++;
      if( non_forward_ ) state.iterDiff_.push_back(x.clone());     // Create new memory
      if( non_inverse_ ) state.gradDiff_.push_back(grad.clone());  // Create new memory
    }
    else {
      if( non_forward_ ) {
        state.iterDiff_.push_back(state.iterDiff_[0]);  // Move first element to the last
        state.iterDiff_.erase(state.iterDiff_.begin()); // Remove first element of s list
        state.product2_.erase(state.product2_.begin()); // Remove first element of rho list
      }
      if( non_inverse_ ) { 
        state.gradDiff_.push_back(state.gradDiff_[0]);  // Move first element to the last
        state.gradDiff_.erase(state.gradDiff_.begin()); // Remove first element of y list
        state.product_.erase(state.product_.begin());   // Remove first element of rho list
      }
    }

    if( non_forward_ ) {
      state.iterDiff_[state.current_]->set(*Hy_);       // s_k - H_k y_k
      state.product2_.push_back(dotI);                   // (s_k - H_k y_k)' y_k
    } 

    if( non_inverse_ ) {
      state.gradDiff_[state.current_]->set(*Bs_);       // y_k - B_k s_k
      state.product_.push_back(dotF);                    // (y_k - B_k s_k)' s_k
    }

    if (useDefaultScaling_) Bscaling_ = s.apply(*y_)/(snorm*snorm);
  }
} // lSR1<Real>::updateStorage


// Apply lSR1 Approximate Inverse Hessian
template<class Real>
void lSR1<Real>::applyH(       Vector<Real>& Hv, 
                         const Vector<Real>& v ) const {

  auto& state = Secant<Real>::getState();

  if (state.mode_ == Secant<Real>::Mode::Inverse || 
      state.mode_ == Secant<Real>::Mode::Both) {
    // Apply initial Hessian approximation to v
    H0called_ = false;
    applyH0(Hv,v);
    // Apply rank one updates
    if (state.current_ > -1) {
      Real prod(0);
      if (!H0called_) prim_->set(v.dual());
      for (int i = 0; i <= state.current_; ++i) {
        prod = state.iterDiff_[i]->dot(*prim_);
        Hv.axpy(prod/state.product2_[i],*state.iterDiff_[i]);
      }
    }
  }
  else throw std::logic_error(">>> ROL::lSR1::applyH : Not supported in forward mode!");
} // lSR1<Real>::applyH



// Apply Initial lSR1 Approximate Inverse Hessian
template<class Real>
void lSR1<Real>::applyH0(       Vector<Real>& Hv, 
                          const Vector<Real>& v ) const  { 

  auto& state = Secant<Real>::getState();

  if (state.current_ > -1) {
    prim_->set(v.dual());
    Hv.set(*prim_);
    H0called_ = true;
  }
  else Hv.set(v.dual());
  Hv.scale(static_cast<Real>(1)/Bscaling_);
} // lSR1<Real>::applyH0




// Apply lSR1 Approximate Hessian
template<class Real>
void lSR1<Real>::applyB(       Vector<Real>& Bv, 
                         const Vector<Real>& v ) const {

  auto& state = Secant<Real>::getState();

  if (state.current_ > -1) {
    dual_->set(v.dual());
    Bv.set(*dual_);
    B0called_ = true;
  }
  else Bv.set(v.dual());

  Bv.scale(Bscaling_);
}

// Apply Initial lSR1 Approximate Hessian 
template<class Real>
void lSR1<Real>::applyB0(       Vector<Real>& Bv, 
                          const Vector<Real>& v ) const {
  auto& state = Secant<Real>::getState();

  if( state.mode_ == Secant<Real>::Mode::Forward || 
      state.mode_ == Secant<Real>::Mode::Both) {
    // Apply initial Hessian approximation to v
    B0called_ = false;
    applyB0(Bv,v);
    // Apply rank one updates
    if (state.current_ > -1) {
      Real prod(0);
      if (!B0called_) dual_->set(v.dual());
      for (int i = 0; i <= state.current_; ++i) {
        prod = state.gradDiff_[i]->dot(*dual_);
        Bv.axpy(prod/state.product_[i],*state.gradDiff_[i]);
      }
    }
  }
  else throw std::logic_error(">>> ROL::lSR1::applyB : Not supported in inverse mode!");
} // lSR1<Real>::applyB0

} // namespace ROL2

#endif // ROL2_LSR1_DEF_HPP
