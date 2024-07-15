// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_FINITEDIFFERENCEDEF_HPP
#define ROL_FINITEDIFFERENCEDEF_HPP

#include <iomanip>

namespace ROL {

namespace details {

using namespace std;
using ::ROL::Finite_Difference_Arrays::shifts;
using ::ROL::Finite_Difference_Arrays::weights;

template<typename Real>
FiniteDifference<Real>::FiniteDifference( const int order ) : order_(order),
  workspace_(makePtr<VectorWorkspace<Real>>()) {}

template<typename Real>
FiniteDifference<Real>::FiniteDifference( const int order,
                                          const Ptr<VectorWorkspace<Real>>& workspace ) : 
  order_(order), workspace_(workspace) {}

template<typename Real>
Real FiniteDifference<Real>::operator()( f_scalar_t<Real>& f_value,
                                         f_update_t<Real>& f_update,
                                         const Vector<Real>& v, 
                                         const Vector<Real>& x,
                                         const Real h ) const {
  f_update(x);

  Real f  = f_value(x); 
  Real fx = weights[order_-1][0] * f;
  auto xc = workspace_->copy(x);

  for(int j=0; j<order_; ++j) {
    xc->axpy( h*shifts[order_-1][j], v );

    // Only evaluate at shifts where the weight is nonzero  
    if( weights[order_-1][j+1] != 0 ) {
      f_update(*xc);
      fx += weights[order_-1][j+1] * f_value(*xc);
    }
  }
  return fx/h;
}

template<typename Real>
void FiniteDifference<Real>::operator()( f_vector_t<Real>& f_value, 
                                         f_update_t<Real>& f_update,
                                         V& Jv, 
                                         const V& v, 
                                         const V& x,
                                         const Real h ) const {
  auto xc  = workspace_->copy(x);
  auto Jvc = workspace_->clone(Jv);

  f_update(x);
  f_value(Jv,x);
  Jv.scale(weights[order_-1][0]);

  for(int j=0; j<order_; ++j) {

    xc->axpy( h*shifts[order_-1][j], v );

    if( weights[order_-1][j+1] != 0 ) {
      f_update(*xc);
      f_value(*Jvc,*xc);  
      Jv.axpy(weights[order_-1][j+1],*Jvc);
    }
  }
  Jv.scale(1.0/h);
}

} // namespace details

} // namespace ROL

#endif // ROL_FINITEDIFFERENCEDEF_HPP

