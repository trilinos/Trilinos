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


#pragma once
#ifndef ROL_CONSTRAINT_CHECKINTERFACE_HPP
#define ROL_CONSTRAINT_CHECKINTERFACE_HPP

#include "ROL_FunctionBindings.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {
namespace details {


using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class Constraint_CheckInterface {
private:
  Constraint<Real>& con_;
  Real tol_;

public:
  using V = Vector<Real>;

  Constraint_CheckInterface( Constraint<Real>& con ) : 
    con_(con), tol_(sqrt(ROL_EPSILON<Real>())) {}
   
  f_update_t<Real> update() {
    return bind( &Constraint<Real>::update, &con_, ph::_1, true, 0 );
  }

  f_vector_t<Real> value() {
    return bind( &Constraint<Real>::value, &con_, ph::_1, ph::_2, tol_);
  }


  f_dderiv_t<Real> jacobian() {
    return bind( &Constraint<Real>::applyJacobian, &con_, ph::_1, ph::_2, ph::_3, tol_);
  }

  // Provide a vector in the dual constraint space
  f_dderiv_t<Real> adjointJacobian( ) {
    return bind( static_cast<void (Constraint<Real>::*)
                              ( V&, const V&, const V&, Real& )>
               (&Constraint<Real>::applyAdjointJacobian), 
                &con_, ph::_1, ph::_2, ph::_3, tol_);
  }

  f_dderiv_t<Real> adjointHessian( const V& l ) {
    return bind( &Constraint<Real>::applyAdjointHessian, &con_, ph::_1, cref(l), ph::_2, ph::_3, tol_);
  }


}; // Constraint_CheckInterface

} // namespace details 

using details::Constraint_CheckInterface;

template<typename Real>
Constraint_CheckInterface<Real> make_check( Constraint<Real>& con ) {
  return Constraint_CheckInterface<Real>(con);
}



} // namespace ROL


#endif // ROL_CONSTRAINT_CHECKINTERFACE_HPP

