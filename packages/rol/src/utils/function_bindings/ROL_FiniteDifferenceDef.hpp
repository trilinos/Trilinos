// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (width_14) Sandia Corporation
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

