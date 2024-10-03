// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP
#define ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP

#include <functional>

#include "ROL_DynamicObjective.hpp"
#include "ROL_FunctionBindings.hpp"


namespace ROL {
namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class DynamicObjective_CheckInterface {
private:

  using V   = Vector<Real>;
  using Obj = DynamicObjective<Real>;

  Obj& obj_;
  Real tol_;
  TimeStamp<Real> ts_;

public:

  DynamicObjective_CheckInterface( Obj& obj ) : obj_(obj) {
    ts_.t.resize(2);
    ts_.t.at(0) = 0.0;
    ts_.t.at(1) = 1.0;
    ts_.k = 0;    
  }

  DynamicObjective_CheckInterface( Obj& obj, TimeStamp<Real>& ts ) : 
    obj_(obj), ts_(ts) { }


  f_update_t<Real> update_uo( const V& un, const V& z ) {
    return bind( &Obj::update, &obj_, ph::_1, cref(un), cref(z), ts_ );
  }

  f_update_t<Real> update_un( const V& uo, const V& z ) {
    return bind( &Obj::update, &obj_, cref(uo), ph::_1, cref(z), ts_ );
  }

  f_update_t<Real> update_z( const V& uo, const V& un ) {
    return bind( &Obj::update, &obj_, cref(uo), cref(un), ph::_1, ts_ );
  }

  //----------------------------------------------------------------------------

  f_scalar_t<Real> value_uo( const V& un, const V& z ) {
    return bind( &Obj::value, &obj_, ph::_1, cref(un), cref(z), ts_ );
  }

  f_scalar_t<Real> value_un( const V& uo, const V& z ) {
    return bind( &Obj::value, &obj_, cref(uo), ph::_1, cref(z), ts_ );
  }

  f_scalar_t<Real> value_z( const V& uo, const V& un ) {
    return bind( &Obj::value, &obj_, cref(uo), cref(un), ph::_1, ts_ );
  }

  //----------------------------------------------------------------------------

  f_vector_t<Real> gradient_uo( const V& un, const V& z ) {
    return bind( &Obj::gradient_uo, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> gradient_un( const V& uo, const V& z ) {
    return bind( &Obj::gradient_un, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }
  
  f_vector_t<Real> gradient_z( const V& uo, const V& un ) {
    return bind( &Obj::gradient_z, &obj_, ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  // For hessian checks
  f_vector_t<Real> gradient_uo_uo( const V& un, const V& z ) {
    return bind( &Obj::gradient_uo, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> gradient_uo_un( const V& uo, const V& z ) {
    return bind( &Obj::gradient_uo, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  f_vector_t<Real> gradient_uo_z( const V& uo, const V& un ) {
    return bind( &Obj::gradient_uo, &obj_, ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  f_vector_t<Real> gradient_un_uo( const V& un, const V& z ) {
    return bind( &Obj::gradient_un, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> gradient_un_un( const V& uo, const V& z ) {
    return bind( &Obj::gradient_un, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  f_vector_t<Real> gradient_un_z( const V& uo, const V& un ) {
    return bind( &Obj::gradient_un, &obj_, ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  f_vector_t<Real> gradient_z_uo( const V& un, const V& z ) {
    return bind( &Obj::gradient_z, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> gradient_z_un( const V& uo, const V& z ) {
    return bind( &Obj::gradient_z, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  f_vector_t<Real> gradient_z_z( const V& uo, const V& un ) {
    return bind( &Obj::gradient_z, &obj_, ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> hessVec_uo_uo( const V& un, const V& z ) {
    return bind( &Obj::hessVec_uo_uo, &obj_, ph::_1, ph::_2, ph::_3, cref(un), cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_uo_un( const V& uo, const V& z ) {
    return bind( &Obj::hessVec_uo_un, &obj_, ph::_1, ph::_2, cref(uo), ph::_3, cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_uo_z( const V& uo, const V& un ) {
    return bind( &Obj::hessVec_uo_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3, ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> hessVec_un_uo( const V& un, const V& z ) {
    return bind( &Obj::hessVec_un_uo, &obj_, ph::_1, ph::_2, ph::_3, cref(un), cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_un_un( const V& uo, const V& z ) {
    return bind( &Obj::hessVec_un_un, &obj_, ph::_1, ph::_2, cref(uo), ph::_3, cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_un_z( const V& uo, const V& un ) {
    return bind( &Obj::hessVec_un_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3, ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> hessVec_z_uo( const V& un, const V& z  ) {
    return bind( &Obj::hessVec_z_uo, &obj_, ph::_1, ph::_2, ph::_3, cref(un), cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_z_un( const V& uo, const V& z ) {
    return bind( &Obj::hessVec_z_un, &obj_, ph::_1, ph::_2, cref(uo), ph::_3, cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_z_z( const V& uo, const V& un ) {
    return bind( &Obj::hessVec_z_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3, ts_ );
  }


}; // DynamicObjective_CheckInterface


} // namespace details 

using details::DynamicObjective_CheckInterface;

template<typename Real>
DynamicObjective_CheckInterface<Real> make_check( DynamicObjective<Real>& obj ) {
  return DynamicObjective_CheckInterface<Real>( obj );
}

template<typename Real>
DynamicObjective_CheckInterface<Real> make_check( DynamicObjective<Real>& obj, 
                                                  TimeStamp<Real>& ts ) {
  return DynamicObjective_CheckInterface<Real>( obj, ts );
}


} // namespace ROL


#endif // ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP

