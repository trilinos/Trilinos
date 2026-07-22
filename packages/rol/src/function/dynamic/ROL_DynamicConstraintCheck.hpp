// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICCONSTRAINTCHECK_HPP
#define ROL_DYNAMICCONSTRAINTCHECK_HPP

#include "ROL_DynamicConstraint_CheckInterface.hpp"
#include "ROL_ValidateFunction.hpp"
#include <string>


namespace ROL {


template<typename Real>
struct DynamicConstraintCheck {

   static void check( DynamicConstraint<Real>& con,
                      ValidateFunction<Real>& validator,
                      const Vector<Real>& uo,
                      const Vector<Real>& un,
                      const Vector<Real>& z,
                      const std::vector<std::string>& methods ) {

    auto con_check = make_check( con );
    check( con_check, validator, uo, un, z, methods );
  }

  static void check( DynamicConstraint<Real>& con,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     TimeStamp<Real>& timeStamp,
                     const std::vector<std::string>& methods ) {
 
    auto con_check = make_check( con, timeStamp );
    check( con_check, validator, uo, un, z, methods );
  }

  static void check( DynamicConstraint_CheckInterface<Real>& con_check, 
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     const std::vector<std::string>& methods ) {

    auto c  = uo.clone();
    auto vu = uo.clone();
    auto vz = z.clone();
    auto l  = uo.dual().clone();
     
    c->randomize();
    vu->randomize();
    vz->randomize();
    l->randomize();

    auto update_uo = con_check.update_uo( un, z );
    auto update_un = con_check.update_un( uo, z );
    auto update_z  = con_check.update_z( un, uo );
  
    auto value_uo = con_check.value_uo( un, z );
    auto value_un = con_check.value_un( uo, z );
    auto value_z  = con_check.value_z( uo, un );
    
    //-------------------------------------------------------------------------
    // Check Jacobian components
    if( std::find(methods.begin(),methods.end(),"applyJacobian_uo") != methods.end() ) {
      auto J = con_check.jacobian_uo( un, z );
      validator.derivative_check( value_uo, J, update_uo, *c, *vu, uo, "norm(J_uo*vec)" );
    } //else os << "\napplyJacobian_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyJacobian_un") != methods.end() ) {
      auto J = con_check.jacobian_un( uo, z );
      validator.derivative_check( value_un, J, update_un, *c, *vu, un, "norm(J_un*vec)" );
    } //else os << "\napplyJacobian_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyJacobian_z") != methods.end() ) {
      auto J = con_check.jacobian_z( uo, un );
      validator.derivative_check( value_z,  J, update_z,  *c, *vz, z,  "norm(J_z*vec)"  );
    } //else os << "\napplyJacobian_z not implemented.\n";


    //-------------------------------------------------------------------------
    // Check Adjoint Jacobian component consistencies
    if( std::find(methods.begin(),methods.end(),"applyAdjointJacobian_uo") != methods.end() ) {
      auto J = con_check.jacobian_uo( un, z );
      auto aJ = con_check.adjointJacobian_uo( un, z );
      validator.adjoint_consistency_check( J, aJ, update_uo, *c, *vu, uo, 
                                           "Jacobian with respect to uo", "J_uo");
    } //else os << "\napplyAdjointJacobian_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointJacobian_un") != methods.end() ) {
      auto J  = con_check.jacobian_un( uo, z );
      auto aJ = con_check.adjointJacobian_un( uo, z );
      validator.adjoint_consistency_check( J, aJ, update_un, *c, *vu, un, 
                                           "Jacobian with respect to un", "J_un");
    } //else os << "\napplyAdjointJacobian_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointJacobian_z") != methods.end() ) {
      auto J  = con_check.jacobian_z(  uo, un );
      auto aJ = con_check.adjointJacobian_z(  uo, un );
      validator.adjoint_consistency_check( J,  aJ, update_z, *vz, *c, z,  
                                           "Jacobian with respect to z", "J_z");
    } //else os << "\napplyAdjointJacobian_z not implemented.\n";

 
    //-------------------------------------------------------------------------
    // Check inverses
    if( std::find(methods.begin(),methods.end(),"solve") != methods.end() ) {
      auto S = con_check.solve_un( uo, z );
      validator.solve_check( S, value_un, update_un, *c, un, "Dynamic Constraint");
    } //else os << "\nsolve not implemented.\n";


    if( std::find(methods.begin(),methods.end(),"applyInverseJacobian_un") != methods.end() ) {
      auto J  = con_check.jacobian_un( uo, z );
      auto iJ = con_check.inverseJacobian_un( uo, z );
      validator.inverse_check( J, iJ, update_un, *vu, un, 
                               "Jacobian with respect to un", "J_un");
    } //else os << "\napplyInverseJacobian_un not implemented.\n";


    if( std::find(methods.begin(),methods.end(),"applyInverseAdjointJacobian_un") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_un( uo, z );
      auto iaJ = con_check.inverseAdjointJacobian_un( uo, z );
      validator.inverse_check( aJ, iaJ, update_un, *vu, un, 
                               "adjoint Jacobian with respect to un", "aJ_un");
    } //else os << "\napplyInverseAdjointJacobian_un not implemented.\n";


    //-------------------------------------------------------------------------
    // Check Adjoint Hessian components
    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_uo_uo") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_uo_uo( un, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_uo_uo( un, z, *l );
      validator.derivative_check( aJl, aH, update_uo, *c, *vu, uo, "H_uo_uo");
    } //else os << "\napplyAdjointHessian_uo_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_uo_un") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_un_uo( un, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_uo_un( un, z, *l );
      validator.derivative_check( aJl, aH, update_uo, *c, *vu, uo, "H_uo_un");
    } //else os << "\napplyAdjointHessian_uo_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_uo_z") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_z_uo( un, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_uo_z( un, z, *l );
      validator.derivative_check( aJl, aH, update_uo, *vz, *vu, uo, "H_uo_z");
 } //else os << "\napplyAdjointHessian_uo_z not implemented.\n";


    
    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_un_uo") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_uo_un( uo, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_un_uo( uo, z, *l );
      validator.derivative_check( aJl, aH, update_un, *c, *vu, un, "H_un_uo");
    } //else os << "\napplyAdjointHessian_un_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_un_un") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_un_un( uo, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_un_un( uo, z, *l );
      validator.derivative_check( aJl, aH, update_un, *c, *vu, un, "H_un_un");
    } //else os << "\napplyAdjointHessian_un_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_un_z") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_z_un( uo, z );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_un_z( un, z, *l );
      validator.derivative_check( aJl, aH, update_un, *vz, *vu, un, "H_un_z");
    } //else os << "\napplyAdjointHessian_uo_uo not implemented.\n";



    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_z_uo") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_uo_z( uo, un );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_z_uo( uo, un, *l );
      validator.derivative_check( aJl, aH, update_z, *c, *vz, z, "H_z_uo");
    } //else os << "\napplyAdjointHessian_z_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_z_un") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_un_z( uo, un );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_z_un( uo, un, *l );
      validator.derivative_check( aJl, aH, update_z, *c, *vz, z, "H_z_un");
    } //else os << "\napplyAdjointHessian_z_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyAdjointHessian_z_z") != methods.end() ) {
      auto aJ  = con_check.adjointJacobian_z_z( uo, un );
      auto aJl = fix_direction( aJ, *l );
      auto aH  = con_check.adjointHessian_z_z( uo, un, *l );
      validator.derivative_check( aJl, aH, update_z, *vz, *vz, z, "H_z_z");
    } //else os << "\napplyAdjointHessian_z_z not implemented.\n";

  } // check()

  static void check( DynamicConstraint<Real>& con,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z) {
    std::vector<std::string> methods = {"applyJacobian_uo",
                                        "applyJacobian_un",
                                        "applyJacobian_z",
                                        "applyAdjointJacobian_uo",
                                        "applyAdjointJacobian_un",
                                        "applyAdjointJacobian_z",
                                        "solve",
                                        "applyInverseJacobian_un",
                                        "applyInverseAdjointJacobian_un",
                                        "applyAdjointHessian_uo_uo",
                                        "applyAdjointHessian_uo_un",
                                        "applyAdjointHessian_uo_z",
                                        "applyAdjointHessian_un_uo",
                                        "applyAdjointHessian_un_un",
                                        "applyAdjointHessian_un_z",
                                        "applyAdjointHessian_z_uo",
                                        "applyAdjointHessian_z_un",
                                        "applyAdjointHessian_z_z"};
    check(con, validator, uo, un, z, methods);
  }

  static void check( DynamicConstraint<Real>& con,
                     ValidateFunction<Real>& validator,
                     const Vector<Real>& uo,
                     const Vector<Real>& un,
                     const Vector<Real>& z,
                     TimeStamp<Real> &ts) {
    std::vector<std::string> methods = {"applyJacobian_uo",
                                        "applyJacobian_un",
                                        "applyJacobian_z",
                                        "applyAdjointJacobian_uo",
                                        "applyAdjointJacobian_un",
                                        "applyAdjointJacobian_z",
                                        "solve",
                                        "applyInverseJacobian_un",
                                        "applyInverseAdjointJacobian_un",
                                        "applyAdjointHessian_uo_uo",
                                        "applyAdjointHessian_uo_un",
                                        "applyAdjointHessian_uo_z",
                                        "applyAdjointHessian_un_uo",
                                        "applyAdjointHessian_un_un",
                                        "applyAdjointHessian_un_z",
                                        "applyAdjointHessian_z_uo",
                                        "applyAdjointHessian_z_un",
                                        "applyAdjointHessian_z_z"};
    check(con, validator, uo, un, z, ts, methods);
  }

  
};

} // namespace ROL

#endif // ROL_DYNAMICCONSTRAINTCHECK_HPP


