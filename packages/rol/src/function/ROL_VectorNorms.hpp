// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTORNORMS_H
#define ROL_VECTORNORMS_H

#include "ROL_Vector.hpp"

namespace ROL {

template<class Real>
Real normL1( const Vector<Real> &x ) {
 
  ROL::Ptr<Vector<Real> > xabs = x.clone();
  xabs->set(x);

  xabs->applyUnary(Elementwise::AbsoluteValue<Real>());
  return xabs->reduce(Elementwise::ReductionSum<Real>());
}

template<class Real, class Exponent>
Real normLp( const Vector<Real> &x, Exponent p ) {

  ROL::Ptr<Vector<Real> > xabsp = x.clone();
  xabsp->set(x);
  xabsp->applyUnary(Elementwise::AbsoluteValue<Real>());
  xabsp->applyUnary(Elementwise::Power<Real>(p)); 
  Real sum = xabsp->reduce(Elementwise::ReductionSum<Real>());
  return std::pow(sum,1.0/p);
}

template<class Real>
Real normLinf( const Vector<Real> &x ) {
 
  ROL::Ptr<Vector<Real> > xabs = x.clone();
  xabs->set(x);

  xabs->applyUnary(Elementwise::AbsoluteValue<Real>());
  return xabs->reduce(Elementwise::ReductionMax<Real>());
}



} // namespace ROL

#endif // ROL_VECTORNORMS_H
