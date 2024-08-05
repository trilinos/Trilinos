// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_L1_PENALTY_DEF_HPP
#define ROL_OED_L1_PENALTY_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
L1Penalty<Real>::L1Penalty()
  : ProfiledClass<Real,std::string>("OED::L1Penalty") {}

template<typename Real>
std::vector<Real>& L1Penalty<Real>::getData(Vector<Real> &x) const {
  return *dynamic_cast<StdVector<Real>&>(x).getVector();
}

template<typename Real>
const std::vector<Real>& L1Penalty<Real>::getConstData(const Vector<Real> &x) const {
  return *dynamic_cast<const StdVector<Real>&>(x).getVector();
}

template<typename Real>
void L1Penalty<Real>::sumAll(Real *input, Real *output, int size, const Vector<Real> &x) const {
  dynamic_cast<const ProbabilityVector<Real>&>(x).getBatchManager()->sumAll(input,output,size);
}

template<typename Real>
Real L1Penalty<Real>::value( const Vector<Real> &x, Real &tol ) {
  startTimer("value");
  const std::vector<Real> &xdata = getConstData(x);
  Real mval(0), gval(0); 
  for (const auto &xi : xdata) mval += xi;
  sumAll(&mval,&gval,1,x);
  stopTimer("value");
  return gval;
}

template<typename Real>
void L1Penalty<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  startTimer("gradient");
  g.setScalar(static_cast<Real>(1));
  stopTimer("gradient");
}

template<typename Real>
void L1Penalty<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  startTimer("hessVec");
  hv.zero();
  stopTimer("hessVec");
}

} // End OED Namespace
} // End ROL Namespace

#endif
