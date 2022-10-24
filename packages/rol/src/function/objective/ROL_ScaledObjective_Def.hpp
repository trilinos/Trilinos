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

#ifndef ROL_SCALED_OBJECTIVE_DEF_HPP
#define ROL_SCALED_OBJECTIVE_DEF_HPP

namespace ROL {

template<typename Real>
void ScaledObjective<Real>::update(const Vector<Real> &x, UpdateType type, int iter) {
  obj_->update(x,type,iter);
}

template<typename Real>
void ScaledObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  obj_->setParameter(param);
}

template<typename Real>
Real ScaledObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return scale_ * obj_->value(x,tol);
}

template<typename Real>
void ScaledObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  obj_->gradient(g,x,tol);
  g.scale(scale_);
}

template<typename Real>
void ScaledObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->hessVec(hv,v,x,tol);
  hv.scale(scale_);
}

template<typename Real>
void ScaledObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->invHessVec(hv,v,x,tol);
  hv.scale(static_cast<Real>(1)/scale_);
}

template<typename Real>
void ScaledObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->precond(Pv,v,x,tol);
  Pv.scale(static_cast<Real>(1)/scale_);
}

} // End ROL Namespace

#endif
