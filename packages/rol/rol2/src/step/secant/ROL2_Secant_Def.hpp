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
#ifndef ROL2_SECANT_DEF_H
#define ROL2_SECANT_DEF_H

/** \class ROL2::Secant
    \brief Provides interface for and implements limited-memory secant operators.
*/

namespace ROL2 {

template<class Real>
Secant<Real>::Secant( int                M, 
                      bool               useDefaultScaling, 
                      Real               Bscaling, 
                      Secant<Real>::Mode mode )
 : state_( makePtr<State>(M,mode) ),
   useDefaultScaling_(useDefaultScaling), 
   Bscaling_(Bscaling),
   mode_(mode) {}

template<class Real>
void Secant<Real>::updateStorage( const Vector<Real>& x,  
                                  const Vector<Real>& grad,
                                  const Vector<Real>& gp, 
                                  const Vector<Real>& s,
                                        Real          snorm,       
                                        int           iter ) {
  const Real one(1);
  if ( !isInitialized_ ) {
    state_->iterate_ = x.clone();
    y_              = grad.clone();
    isInitialized_  = true;
  }
  state_->iterate_->set(x);
  state_->iter_ = iter;
  y_->set(grad);
  y_->axpy(-one,gp);

  //Real sy = s.dot(y_->dual());
  Real sy = s.apply(*y_);
  if (sy > ROL_EPSILON<Real>*snorm*snorm) {
    if (state_->current_ < state_->storage_-1) {
      state_->current_++;                                // Increment Storage
      state_->iterDiff_.push_back(s.clone());            // Create new memory
      state_->gradDiff_.push_back(grad.clone());         // Create new memory
    }
    else {
      state_->iterDiff_.push_back(state_->iterDiff_[0]);  // Move first element to the last
      state_->gradDiff_.push_back(state_->gradDiff_[0]);  // Move first element to the last
      state_->iterDiff_.erase(state_->iterDiff_.begin()); // Remove first element of s list 
      state_->gradDiff_.erase(state_->gradDiff_.begin()); // Remove first element of y list
      state_->product_.erase(state_->product_.begin());   // Remove first element of rho list
    }
    state_->iterDiff_[state_->current_]->set(s);          // s=x_{k+1}-x_k
    state_->gradDiff_[state_->current_]->set(*y_);        // y=g_{k+1}-g_k
    state_->product_.push_back(sy);                      // ys=1/rho  
  }
}

template<class Real>
void Secant<Real>::applyH0(       Vector<Real>& Hv,
                            const Vector<Real>& v ) const {
  Hv.set(v.dual());
  if (useDefaultScaling_) {
    if (state_->iter_ != 0 && state_->current_ != -1) {
      Real yy = state_->gradDiff_[state_->current_]->dot(*(state_->gradDiff_[state_->current_]));
      Hv.scale(state_->product_[state_->current_]/yy);
    }
  }
  else Hv.scale(static_cast<Real>(1)/Bscaling_);
}

template<class Real>
void Secant<Real>::applyB0(       Vector<Real>& Bv, 
                            const Vector<Real>& v ) const {
  Bv.set(v.dual());
  if (useDefaultScaling_) {
    if (state_->iter_ != 0 && state_->current_ != -1) {
      Real yy = state_->gradDiff_[state_->current_]->dot(*(state_->gradDiff_[state_->current_]));
      Bv.scale(yy/state_->product_[state_->current_]);
    }
  }
  else Bv.scale(Bscaling_);
}

template<class Real>
void Secant<Real>::test( std::ostream &os ) const {
  if (isInitialized_) {
    auto v  = state_->iterate_->clone();
    auto Hv = state_->iterate_->clone();
    auto Bv = state_->iterate_->dual().clone();
    const Real one(1);

    // Print BHv -> Should be v
    v->randomize(-one,one);
    applyH(*Hv,*v);
    applyB(*Bv,*Hv);
    v->axpy(-one,*Bv);
    os << " ||BHv-v|| = " << v->norm() << std::endl;

    // Print HBv -> Should be v
    v->randomize(-one,one);
    applyB(*Bv,*v);
    applyH(*Hv,*Bv);
    v->axpy(-one,*Hv);
    os << " ||HBv-v|| = " << v->norm() << std::endl;
  }
}

template<class Real>
Ptr<Secant<Real>> 
Secant<Real>::create( Secant<Real>::Type type,
                      int                L,
                      int                BBtype ) {
  switch( type ) {
    case Type::LBFGS:           return makePtr<lBFGS<Real>>(L);
    case Type::LDFP:            return makePtr<lDFP<Real>>(L);
    case Type::LSR1:            return makePtr<lSR1<Real>>(L);
    case Type::BarzilaiBorwein: return makePtr<BarzilaiBorwein<Real>>(BBtype);
    default:                    return nullPtr; // Should we throw an exception here?
  }
}


template<class Real>
Ptr<Secant<Real>> 
Secant<Real>::create( ParameterList&     parlist,
                      Secant<Real>::Mode mode ) {

  auto& slist = parlist.sublist("General").sublist("Secant");
  auto sstr = slist.get("Type","Limited-Memory BFGS");
  int L    = slist.get("Maximum Storage",10);
  int BB   = slist.get("Barzilai-Borwein",1);
  bool uds = slist.get("Use Default Scaling",true);
  Real s   = slist.get("Initial Hessian Scale",1.0);

  switch( type_dict[sstr] ) {
    case Type::LBFGS:           return makePtr<lBFGS<Real>>(L,uds,s);
    case Type::LDFP:            return makePtr<lDFP<Real>>(L,uds,s);
    case Type::LSR1:            return makePtr<lSR1<Real>>(L,uds,s,mode);
    case Type::BarzilaiBorwein: return makePtr<BarzilaiBorwein<Real>>(BB);
    default:                    return nullPtr; // Should we throw an exception here?
  }
}



} // namespace ROL2

#endif // ROL2_SECANT_DEF_H
