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

#ifndef ROL_OED_PROBABILITY_CONSTRAINT_DEF_HPP
#define ROL_OED_PROBABILITY_CONSTRAINT_DEF_HPP

namespace ROL {
namespace OED {

/***************************************************************************/
/* Begin Accessor Functions                                                */
/***************************************************************************/
template<typename Real>
std::vector<Real>& ProbabilityConstraint<Real>::getData(Vector<Real> &x) const {
  return *dynamic_cast<StdVector<Real>&>(x).getVector();
}

template<typename Real>
const std::vector<Real>& ProbabilityConstraint<Real>::getConstData(const Vector<Real> &x) const {
  return *dynamic_cast<const StdVector<Real>&>(x).getVector();
}

template<typename Real>
void ProbabilityConstraint<Real>::sumAll(Real *input, Real *output, int size, const Vector<Real> &x) const {
  dynamic_cast<const ProbabilityVector<Real>&>(x).getBatchManager()->sumAll(input,output,size);
}
/***************************************************************************/
/* End Accessor Functions                                                  */
/***************************************************************************/

template<typename Real>
ProbabilityConstraint<Real>::ProbabilityConstraint(const Vector<Real> &p,
                      bool useScale,
                      Real scale)
  : useScale_(useScale), scale_(scale) {
  if (useScale_ && scale_ < static_cast<Real>(0)) {
    Real N(p.dimension());
    scale_ = static_cast<Real>(1)/std::sqrt(N);
    //scale_ = static_cast<Real>(1)/std::pow(N,2.0/3.0);
  }
}

template<typename Real>
void ProbabilityConstraint<Real>::value(Vector<Real> &c,
           const Vector<Real> &x,
           Real &tol) {
  c.zero();
  std::vector<Real>       &cdata = getData(c);
  const std::vector<Real> &xdata = getConstData(x);
  Real mval(0), gval(0); 
  for (const auto &xi : xdata) mval += xi;
  sumAll(&mval,&gval,1,x);
  cdata[0] = gval-static_cast<Real>(1);
  if (useScale_) c.scale(scale_);
}

template<typename Real>
void ProbabilityConstraint<Real>::applyJacobian(Vector<Real> &jv,
                   const Vector<Real> &v,
                   const Vector<Real> &x,
                   Real &tol) {
  jv.zero();
  std::vector<Real>       &jdata = getData(jv);
  const std::vector<Real> &vdata = getConstData(v);
  Real mval(0), gval(0); 
  for (const auto &vi : vdata) mval += vi;
  sumAll(&mval,&gval,1,x);
  jdata[0] = gval;
  if (useScale_) jv.scale(scale_);
}

template<typename Real>
void ProbabilityConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                          const Vector<Real> &v,
                          const Vector<Real> &x,
                          Real &tol) {
  const std::vector<Real> &vdata = getConstData(v);
  ajv.setScalar(vdata[0]);
  if (useScale_) ajv.scale(scale_);
}

template<typename Real>
void ProbabilityConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                         const Vector<Real> &u,
                         const Vector<Real> &v,
                         const Vector<Real> &x,
                         Real &tol) {
  ahuv.zero();
}

} // End OED Namespace
} // End ROLNamespace

#endif
