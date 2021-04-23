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

#ifndef ROL_QUADRATIC_OBJECTIVE_DEF_H
#define ROL_QUADRATIC_OBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
QuadraticObjective<Real>::QuadraticObjective(const Ptr<const LinearOperator<Real>> &H,
                                             const Ptr<const Vector<Real>>         &g,
                                             Real                                   c)
  : H_(H), g_(g), c_(c) {
  tmp_ = g_->clone();
}

template<typename Real>
Real QuadraticObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  H_->apply(*tmp_,x,tol);
  tmp_->scale(static_cast<Real>(0.5));
  tmp_->plus(*g_);
  //return x.dot(tmp_->dual()) + c_;
  return x.apply(*tmp_) + c_;
}

template<typename Real>
void QuadraticObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  H_->apply(g,x,tol);
  g.plus(*g_);
}

template<typename Real>
void QuadraticObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  H_->apply(hv,v,tol);
}

template<typename Real>
void QuadraticObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  H_->applyInverse(hv,v,tol);
}

} // namespace ROL

#endif
