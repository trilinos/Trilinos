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

#ifndef ROL_STDOBJECTIVE_DEF_H
#define ROL_STDOBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
void StdObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  update(*(xs.getVector()),flag,iter);
}

template<typename Real>
void StdObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  update(*(xs.getVector()),type,iter);
}

template<typename Real>
Real StdObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  return value(*(xs.getVector()),tol);
}

template<typename Real>
void StdObjective<Real>::gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
  const unsigned size = x.size();
  std::vector<Real> y; y.assign(x.begin(),x.end());
  const Real cbrteps = std::cbrt(ROL::ROL_EPSILON<Real>()), one(1);
  Real h(1), xi(0);
  const Real val = value(x,tol);
  for (unsigned i = 0; i < size; ++i) {
    xi   = x[i];
    h    = cbrteps * std::max(std::abs(xi),one) * sgn(xi);
    y[i] = xi + h;
    h    = y[i] - xi;
    update(y);
    g[i] = (value(y,tol) - val)/h;
    y[i] = xi;
  }
  update(x);
}

template<typename Real>
void StdObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  StdVector<Real> gs = dynamic_cast<StdVector<Real>&>(g);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  gradient(*(gs.getVector()),*(xs.getVector()),tol);
}

template<typename Real>
Real StdObjective<Real>::dirDeriv( const std::vector<Real> &x, const std::vector<Real> &d, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
    ">>> ERROR (ROL::StdObjective): dirDeriv not implemented!");
}

template<typename Real>
Real StdObjective<Real>::dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  const StdVector<Real> ds = dynamic_cast<const StdVector<Real>&>(d);
  try {
    return dirDeriv(*(xs.getVector()),*(ds.getVector()),tol);
  }
  catch (std::exception &e) {
    return Objective<Real>::dirDeriv(x,d,tol);
  }
}

template<typename Real>
void StdObjective<Real>::hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
    ">>> ERROR (ROL::StdObjective): hessVec not implemented!");
}

template<typename Real>
void StdObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  try {
    StdVector<Real> hvs = dynamic_cast<StdVector<Real>&>(hv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    hessVec(*(hvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
  }
  catch (std::exception &e) {
    Objective<Real>::hessVec(hv,v,x,tol);
  }
}

template<typename Real>
void StdObjective<Real>::invHessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
    ">>> ERROR (ROL::StdObjective): invHessVec not implemented!");
}

template<typename Real>
void StdObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  StdVector<Real> hvs = dynamic_cast<StdVector<Real>&>(hv);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  invHessVec(*(hvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
}

template<typename Real>
void StdObjective<Real>::precond( std::vector<Real> &Pv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
  Pv.assign(v.begin(),v.end());
}

template<typename Real>
void StdObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  StdVector<Real> Pvs = dynamic_cast<StdVector<Real>&>(Pv);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  precond(*(Pvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
}

template<typename Real>
Real StdObjective<Real>::sgn(Real x) const {
  const Real zero(0), one(1);
  return (x < zero ? -one : one);
}

} // namespace ROL

#endif
