// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_OBJECTIVE_SIMOPT_CHECKINTERFACEDEF_HPP
#define ROL_OBJECTIVE_SIMOPT_CHECKINTERFACEDEF_HPP


#include "ROL_Objective.hpp"
#include <functional>

namespace ROL {
namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class Objective_SimOpt_CheckInterface {
private:
  using V = Vector<Real>;
  Objective_SimOpt<Real>& obj_;
  Real tol_;

public:

  Objective_SimOpt_CheckInterface( Objective<Real>& obj ) : 
    obj_(obj), tol_(sqrt(ROL_EPSILON<Real>())) {}
   
  // Takes a Vector_SimOpt
  f_update_t<Real> update() {
    return bind( &Objective_SimOpt<Real>::update, &obj_, ph::_1, true, 0 );
  }

  // Takes a Vector_SimOpt
  f_scalar_t<Real> value() {
    return bind( &Objective_SimOpt<Real>::value, &obj_, ph::_1, tol_);
  }

  f_vector_t<Real> gradient_1( const V& z ) {
    return bind( &Objective_SimOpt<Real>::gradient_1, &obj_, ph::_1, ph::_2, cref(z), tol_);
  }

  f_vector_t<Real> gradient_2( const V& u ) {
    return bind( &Objective_SimOpt<Real>::gradient_2, &obj_, ph::_1, cref(u), ph::_2, tol_);
  }

  f_dderiv_t<Real> hessVec_11( const  ) {
    return bind( &Objective_SimOpt<Real>::hessVec, &obj_, ph::_1, ph::_2, ph::_3, tol_);
  }

}; // Objective_CheckInterface

} // namespace details

using details::Objective_SimOpt_CheckInterface;
template<typename Real>
Objective_SimOpt_CheckInterface<Real> make_check( Objective_SimOpt<Real>& obj ) {
  return Objective_SimOpt_CheckInterface<Real>(obj);
}

} // namespace ROL



#endif // ROL_OBJECTIVE_SIMOPT_CHECKINTERFACEDEF_HPP

