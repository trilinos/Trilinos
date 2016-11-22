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

#ifndef ROL_LINEARCOMBINATIONOBJECTIVE_H
#define ROL_LINEARCOMBINATIONOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "Teuchos_RCP.hpp"

namespace ROL {

template <class Real>
class LinearCombinationObjective : public Objective<Real> {
private:
  const std::vector<Real> weights_;
  std::vector<Teuchos::RCP<Objective<Real> > > obj_;
  size_t size_;

  Teuchos::RCP<Vector<Real> > xdual_;
  bool initialized_;

public:
  LinearCombinationObjective(const std::vector<Real> &weights,
                             const std::vector<Teuchos::RCP<Objective<Real> > > &obj)
    : Objective<Real>(), weights_(weights), obj_(obj), size_(weights.size()),
      xdual_(Teuchos::null), initialized_(false) {}

  Real value( const Vector<Real> &x, Real &tol ) {
    Real val = 0.;
    for (size_t i = 0; i < size_; i++) {
      val += weights_[i]*obj_[i]->value(x,tol);
    } 
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    if (!initialized_) {
      xdual_ = g.clone();
      initialized_ = true;
    }
    g.zero();
    for (size_t i = 0; i < size_; i++) {
      obj_[i]->gradient(*xdual_,x,tol);
      g.axpy(weights_[i],*xdual_);
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    if (!initialized_) {
      xdual_ = hv.clone();
      initialized_ = true;
    }
    hv.zero();
    for (size_t i = 0; i < size_; i++) {
      obj_[i]->hessVec(*xdual_,v,x,tol);
      hv.axpy(weights_[i],*xdual_);
    }
  }

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) {
    Objective<Real>::setParameter(param);
    for (size_t i = 0; i < size_; ++i) {
      obj_[i]->setParameter(param);
    }
  }
}; // class LinearCombinationObjective

} // namespace ROL

#endif
