// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICOBJECTIVECHECK_HPP
#define ROL_DYNAMICOBJECTIVECHECK_HPP

#include "ROL_DynamicObjective_CheckInterface.hpp"
#include "ROL_ValidateFunction.hpp"
#include <string>

// TODO: Add symmetry check for diagonal Hessian blocks and adjoint consistency
//       for off-diagonal block pairs

namespace ROL {

template<typename Real>
struct DynamicObjectiveCheck {

  static void check( DynamicObjective<Real>& obj,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     const std::vector<std::string>& methods ) {

    auto obj_check = make_check( obj );
    check( obj_check, validator, uo, un, z, methods );
  }

  static void check( DynamicObjective<Real>& obj,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     TimeStamp<Real>& timeStamp,
                     const std::vector<std::string>& methods ) {

    auto obj_check = make_check( obj, timeStamp );
    check( obj_check, validator, uo, un, z, methods );
  }

  static void check( DynamicObjective_CheckInterface<Real>& obj_check, 
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     const std::vector<std::string>& methods ) {

    auto gu = uo.dual().clone();
    auto gz = z.dual().clone();
    auto vu = uo.clone();
    auto vz = z.clone();

    vu->randomize();
    vz->randomize();


    //-------------------------------------------------------------------------
    // Check gradient components
    if( std::find(methods.begin(),methods.end(),"gradient_uo") != methods.end() ) {
      auto value  = obj_check.value_uo( un, z );
      auto grad   = obj_check.gradient_uo( un, z );
      auto update = obj_check.update_uo( un, z );
      validator.derivative_check( value, grad, update, *gu, *vu, uo, "grad_uo'*dir" );
    }
    if( std::find(methods.begin(),methods.end(),"gradient_un") != methods.end() ) {
      auto value  = obj_check.value_un( uo, z );
      auto grad   = obj_check.gradient_un( uo, z );
      auto update = obj_check.update_un( uo, z );
      validator.derivative_check( value, grad, update, *gu, *vu, un, "grad_un'*dir" );
    }
    if( std::find(methods.begin(),methods.end(),"gradient_z") != methods.end() ) {
      auto value  = obj_check.value_z( uo, un );
      auto grad   = obj_check.gradient_z( uo, un );
      auto update = obj_check.update_z( uo, un );
      validator.derivative_check( value, grad, update, *gz, *vz,  z, "grad_z'*dir" );
    }

    //-------------------------------------------------------------------------
    // Check Hessian components
    if( std::find(methods.begin(),methods.end(),"hessVec_uo_uo") != methods.end() ) {
      auto grad    = obj_check.gradient_uo_uo( un, z ); 
      auto hessVec = obj_check.hessVec_uo_uo( un, z );
      auto update  = obj_check.update_uo( un, z );
      validator.derivative_check( grad, hessVec, update, *gu, *vu, uo, "norm(H_uo_uo*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_un") != methods.end() ) {
      auto grad    = obj_check.gradient_uo_un( uo, z ); 
      auto hessVec = obj_check.hessVec_uo_un( uo, z );
      auto update  = obj_check.update_un( uo, z );
      validator.derivative_check( grad, hessVec, update, *gu, *vu, un, "norm(H_uo_un*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_z") != methods.end() ) {
      auto grad    = obj_check.gradient_uo_z( uo, un ); 
      auto hessVec = obj_check.hessVec_uo_z( uo, un );
      auto update  = obj_check.update_z( uo, un );
      validator.derivative_check( grad, hessVec, update, *gu, *vz,  z, "norm(H_uo_z*vec)" );
    }



    if( std::find(methods.begin(),methods.end(),"hessVec_un_uo") != methods.end() ) {
      auto grad    = obj_check.gradient_un_uo( un, z ); 
      auto hessVec = obj_check.hessVec_un_uo( un, z );
      auto update  = obj_check.update_uo( un, z );
      validator.derivative_check( grad, hessVec, update, *gu, *vu, uo, "norm(H_un_uo*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_un_un") != methods.end() ) {
      auto grad    = obj_check.gradient_un_un( uo, z );
      auto hessVec = obj_check.hessVec_un_un( uo, z );
      auto update  = obj_check.update_un( uo, z );
      validator.derivative_check( grad, hessVec, update, *gu, *vu, un, "norm(H_un_un*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_un_z") != methods.end() ) {
      auto grad    = obj_check.gradient_un_z( uo, un );
      auto hessVec = obj_check.hessVec_un_z( uo, un );
      auto update  = obj_check.update_z( uo, un );
      validator.derivative_check( grad, hessVec, update, *gu, *vz,  z, "norm(H_un_z*vec)" );
    }



    if( std::find(methods.begin(),methods.end(),"hessVec_z_uo") != methods.end() ) {
      auto grad    = obj_check.gradient_z_uo( un, z );
      auto hessVec = obj_check.hessVec_z_uo( un, z );
      auto update  = obj_check.update_uo( un, z );
      validator.derivative_check( grad, hessVec, update, *gz, *vu, uo, "norm(H_z_uo*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_z_un") != methods.end() ) {
      auto grad    = obj_check.gradient_z_un( uo, z );
      auto hessVec = obj_check.hessVec_z_un( uo, z );
      auto update  = obj_check.update_un( uo, z );
      validator.derivative_check( grad, hessVec, update, *gz, *vu, un, "norm(H_z_un*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_z_z") != methods.end() ) {
      auto grad    = obj_check.gradient_z_z( uo, un );
      auto hessVec = obj_check.hessVec_z_z( uo, un );
      auto update  = obj_check.update_z( uo, un );
      auto H = obj_check.hessVec_z_z(uo,un);
      validator.derivative_check( grad, hessVec, update, *gz, *vz,  z, "norm(H_z_z*vec)" );
    }
  }

  static void check( DynamicObjective<Real>& obj,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z ) {
    std::vector<std::string> methods = {"gradient_uo",
                                        "gradient_un",
                                        "gradient_z",
                                        "hessVec_uo_uo",
                                        "hessVec_uo_un",
                                        "hessVec_uo_z",
                                        "hessVec_un_uo",
                                        "hessVec_un_un",
                                        "hessVec_un_z",
                                        "hessVec_z_uo",
                                        "hessVec_z_un",
                                        "hessVec_z_z"};
    check(obj, validator, uo, un, z, methods);
  }

  static void check( DynamicObjective<Real>& obj,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     TimeStamp<Real> &ts ) {
    std::vector<std::string> methods = {"gradient_uo",
                                        "gradient_un",
                                        "gradient_z",
                                        "hessVec_uo_uo",
                                        "hessVec_uo_un",
                                        "hessVec_uo_z",
                                        "hessVec_un_uo",
                                        "hessVec_un_un",
                                        "hessVec_un_z",
                                        "hessVec_z_uo",
                                        "hessVec_z_un",
                                        "hessVec_z_z"};
    check(obj, validator, uo, un, z, ts, methods);
  }
}; // DynamicObjectiveCheck

} // namespace ROL

#endif  // ROL_DYNAMICOBJECTIVECHECK_HPP
