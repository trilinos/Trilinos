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


#ifndef ROL_POLYHEDRALPROJECTION_DEF_H
#define ROL_POLYHEDRALPROJECTION_DEF_H

namespace ROL {

template<typename Real>
PolyhedralProjection<Real>::PolyhedralProjection(const Ptr<BoundConstraint<Real>> &bnd)
  : bnd_(bnd), con_(nullPtr) {}

template<typename Real>
PolyhedralProjection<Real>::PolyhedralProjection(const Vector<Real>               &xprim,
                                                 const Vector<Real>               &xdual,
                                                 const Ptr<BoundConstraint<Real>> &bnd,
                                                 const Ptr<Constraint<Real>>      &con,
                                                 const Vector<Real>               &mul,
                                                 const Vector<Real>               &res)
  : bnd_(bnd), con_(con) {
  xprim_ = xprim.clone();
  xdual_ = xdual.clone();
  mul_   = mul.clone();
  res_   = res.clone();
}

template<typename Real>
void PolyhedralProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : No projection implemented!");
  }
}

template<typename Real>
const Ptr<Constraint<Real>> PolyhedralProjection<Real>::getLinearConstraint(void) const {
  return con_;
}

template<typename Real>
const Ptr<Vector<Real>> PolyhedralProjection<Real>::getMultiplier(void) const {
  return mul_;
}

template<typename Real>
const Ptr<Vector<Real>> PolyhedralProjection<Real>::getResidual(void) const {
  return res_;
}

} // namespace ROL

#endif
