// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_DOUBLE_WELL_PENALTY_DEF_HPP
#define ROL_OED_DOUBLE_WELL_PENALTY_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
DoubleWellPenalty<Real>::DoubleWellPenalty()
  : ProfiledClass<Real,std::string>("OED::DoubleWellPenalty") {}

template<typename Real>
std::vector<Real>& DoubleWellPenalty<Real>::getData(Vector<Real> &x) const {
  return *dynamic_cast<StdVector<Real>&>(x).getVector();
}

template<typename Real>
const std::vector<Real>& DoubleWellPenalty<Real>::getConstData(const Vector<Real> &x) const {
  return *dynamic_cast<const StdVector<Real>&>(x).getVector();
}

template<typename Real>
void DoubleWellPenalty<Real>::sumAll(Real *input, Real *output, int size, const Vector<Real> &x) const {
  dynamic_cast<const ProbabilityVector<Real>&>(x).getBatchManager()->sumAll(input,output,size);
}

template<typename Real>
Real DoubleWellPenalty<Real>::value( const Vector<Real> &x, Real &tol ) {
  startTimer("value");
  const Real one(1);
  const std::vector<Real> &xdata = getConstData(x);
  Real mval(0), gval(0); 
  for (const auto &xi : xdata) mval += xi*xi*(one-xi)*(one-xi);
  sumAll(&mval,&gval,1,x);
  stopTimer("value");
  return gval;
}

template<typename Real>
void DoubleWellPenalty<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  startTimer("gradient");
  const Real one(1), two(2), three(3);
  Real xi(0);
  std::vector<Real> &gdata = getData(g);
  const std::vector<Real> &xdata = getConstData(x);
  for (int i = 0; i < static_cast<int>(gdata.size()); ++i) {
    xi       = xdata[i];
    gdata[i] = two*xi*(two*xi*xi-three*xi+one);
  }
  stopTimer("gradient");
}

template<typename Real>
void DoubleWellPenalty<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  startTimer("hessVec");
  const Real one(1), two(2), three(3);
  Real xi(0), vi(0);
  std::vector<Real> &hvdata = getData(hv);
  const std::vector<Real> &xdata = getConstData(x);
  const std::vector<Real> &vdata = getConstData(v);
  for (int i = 0; i < static_cast<int>(hvdata.size()); ++i) {
    xi        = xdata[i];
    vi        = vdata[i];
    hvdata[i] = (three*(two*xi-one)*(two*xi-one)-one)*vi;
  }
  stopTimer("hessVec");
}

} // End OED Namespace
} // End ROL Namespace

#endif
