// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
