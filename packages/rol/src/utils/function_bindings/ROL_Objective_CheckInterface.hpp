// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_OBJECTIVE_CHECKINTERFACE_HPP
#define ROL_OBJECTIVE_CHECKINTERFACE_HPP

#include "ROL_Objective.hpp"
#include "ROL_FunctionBindings.hpp"

namespace ROL {
namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class Objective_CheckInterface {
private:
  using V = Vector<Real>;
  Objective<Real>& obj_;
  Real tol_;

public:

  Objective_CheckInterface( Objective<Real>& obj ) : 
    obj_(obj), tol_(sqrt(ROL_EPSILON<Real>())) {}

  f_update_t<Real> update() {
    return bind( (void(Objective<Real>::*)(const Vector<Real>&,bool,int))&Objective<Real>::update, &obj_, ph::_1, true, 0 );
  }

  f_scalar_t<Real> value() {
    return bind( &Objective<Real>::value, &obj_, ph::_1, tol_);
  }

  f_vector_t<Real> gradient() {
    return bind( &Objective<Real>::gradient, &obj_, ph::_1, ph::_2, tol_);
  }

  f_dderiv_t<Real> hessVec() {
    return bind( &Objective<Real>::hessVec, &obj_, ph::_1, ph::_2, ph::_3, tol_);
  }

}; // Objective_CheckInterface

} // namespace details

using details::Objective_CheckInterface;
template<typename Real>
Objective_CheckInterface<Real> make_check( Objective<Real>& obj ) {
  return Objective_CheckInterface<Real>(obj);
}

} // namespace ROL


#endif // ROL_OBJECTIVE_CHECKINTERFACE_HPP

