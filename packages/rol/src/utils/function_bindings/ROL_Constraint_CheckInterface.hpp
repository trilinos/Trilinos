// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    return bind( (void(Constraint<Real>::*)(const Vector<Real>&,bool,int))&Constraint<Real>::update, &con_, ph::_1, true, 0 );
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

