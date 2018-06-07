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
#ifndef ROL_DYNAMICCONSTRAINTCHECK_HPP
#define ROL_DYNAMICCONSTRAINTCHECK_HPP

#include "ROL_DynamicConstraint_CheckInterface.hpp"
#include "ROL_ValidateFunction.hpp"
#include "ROL_RandomVector.hpp"
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
 
    auto c  = uo.clone();
    auto vu = uo.clone();
    auto vz = z.clone();
    auto l  = uo.dual().clone();
     
    RandomizeVector( *c  );
    RandomizeVector( *vu );
    RandomizeVector( *vz );
    RandomizeVector( *l  );

    std::ostream& os = validator.getStream();

    auto con_check = make_check( con );

    auto update_uo = con_check.update_uo();
    auto update_un = con_check.update_un();
    auto update_z  = con_check.update_z();
  
    auto value_uo = con_check.value_uo( un, z  );
    auto value_un = con_check.value_un( uo, z  );
    auto value_z  = con_check.value_z( uo, un );
    
    //-------------------------------------------------------------------------
    // Check Jacobian components
    if( std::find(methods.begin(),methods.end(),"applyJacobian_uo") != methods.end() ) {
      auto J = con_check.jacobian_uo( un, z  );
      validator.derivative_check( value_uo, J, update_uo, *c, *vu, uo, "norm(J_uo*vec)" );
    } else os << "\napplyJacobian_uo not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyJacobian_un") != methods.end() ) {
      auto J = con_check.jacobian_un( uo, z  );
      validator.derivative_check( value_un, J, update_un, *c, *vu, un, "norm(J_un*vec)" );
    } else os << "\napplyJacobian_un not implemented.\n";

    if( std::find(methods.begin(),methods.end(),"applyJacobian_z") != methods.end() ) {
      auto J = con_check.jacobian_z( uo, un );
      validator.derivative_check( value_z,  J,  update_z,  *c, *vz, z,  "norm(J_z*vec)" );
    } else os << "\napplyJacobian_z not implemented.\n";

//
//    //-------------------------------------------------------------------------
//    // Check Adjoint Jacobian component consistencies
//    if( &Derived::applyAdjointJacobian_uo != &Base::applyAdjointJacobian_uo ) {
//      auto J = con_check.jacobian_uo( un, z );
//      auto aJ = con_check.adjointJacobian_uo( un, z  );
//      validator.adjoint_consistency_check( J, aJ, update_uo, *c, *vu, uo, 
//                                           "Jacobian with respect to uo", "J_uo");
//    } else os << "\napplyAdjointJacobian_uo not implemented.\n";
//
//    if( &Derived::applyAdjointJacobian_un != &Base::applyAdjointJacobian_un ) {
//      auto J  = con_check.jacobian_un( uo, z );
//      auto aJ = con_check.adjointJacobian_un( uo, z );
//      validator.adjoint_consistency_check( J, aJ, update_un, *c, *vu, un, 
//                                           "Jacobian with respect to un", "J_un");
//    } else os << "\napplyAdjointJacobian_un not implemented.\n";
//
//    if( &Derived::applyAdjointJacobian_z != &Base::applyAdjointJacobian_z ) {
//      auto J  = con_check.jacobian_z(  uo, un );
//      auto aJ = con_check.adjointJacobian_z(  uo, un );
//      validator.adjoint_consistency_check( J,  aJ, update_z, *vz, *c, z,  
//                                           "Jacobian with respect to z", "J_z");
//    } else os << "\napplyAdjointJacobian_z not implemented.\n";
//
// 
//    //-------------------------------------------------------------------------
//    // Check inverses
//    if( &Derived::applyInverseJacobian_un != &Base::applyInverseJacobian_un ) {
//      auto J  = con_check.jacobian_uo( un, z );
//      auto iJ = con_check.inverseJacobian_un( uo, z );
//      validator.inverse_check( J, iJ, update_un, *vu, un, 
//                               "Jacobian with respect to un", "J_un");
//    } else os << "\napplyInverseJacobian_un not implemented.\n";
//
//    if( &Derived::applyInverseAdjointJacobian_un != &Base::applyInverseAdjointJacobian_un ) {
//      auto aJ  = con_check.adjointJacobian_un( uo, z );
//      auto iaJ = con_check.inverseAdjointJacobian_un( uo, z );
//      validator.inverse_check( aJ, iaJ, update_un, *vu, un, 
//                               "adjoint Jacobian with respect to un", "J_un");
//    } else os << "\napplyInverseAdjointJacobian_un not implemented.\n";
//
//
//    //-------------------------------------------------------------------------
//    // Check Adjoint Hessian components
//    if( &Derived::applyAdjointHessian_uo_uo != &Base::applyAdjointHessian_uo_uo ) {
//      auto aJ  = con_check.adjointJacobian_uo( un, z );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_uo_uo( un, z, *l );
//      validator.derivative_check( aJl, aJ, update_uo, *c, *vu, uo, "H_uo_uo");
//    } else os << "\napplyAdjointHessian_uo_uo not implemented.\n";
//
//    if( &Derived::applyAdjointHessian_uo_un != &Base::applyAdjointHessian_uo_un ) {
//      auto aJ  = con_check.adjointJacobian_uo( un, z );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_uo_un( uo, z, *l );
//      validator.derivative_check( aJl, aJ, update_un, *c, *vu, un, "H_uo_un");
//    } else os << "\napplyAdjointHessian_uo_un not implemented.\n";
//
//    if( &Derived::applyAdjointHessian_uo_z != &Base::applyAdjointHessian_uo_z ) {
//      auto aJ  = con_check.adjointJacobian_uo( un, z );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_uo_uo( uo, un, *l );
//      validator.derivative_check( aJl, aJ, update_uo, *c, *vz, z, "H_uo_z");
//    } else os << "\napplyAdjointHessian_uo_z not implemented.\n";
//    
//
//    if( &Derived::applyAdjointHessian_un_uo != &Base::applyAdjointHessian_un_uo ) {
//      auto aJ  = con_check.adjointJacobian_un( un, z );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_un_uo( un, z, *l );
//      validator.derivative_check( aJl, aJ, update_uo, *c, *vu, uo, "H_un_uo");
//    } else os << "\napplyAdjointHessian_un_uo not implemented.\n";
//
//    if( &Derived::applyAdjointHessian_un_un != &Base::applyAdjointHessian_un_un ) {
//      auto aJ  = con_check.adjointJacobian_un( uo, z );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_un_un( uo, z, *l );
//      validator.derivative_check( aJl, aJ, update_un, *c, *vu, un, "H_un_un");
//    } else os << "\napplyAdjointHessian_un_un not implemented.\n";
//
//
//    if( &Derived::applyAdjointHessian_un_z != &Base::applyAdjointHessian_un_z ) {
//      auto aJ  = con_check.adjointJacobian_un( uo, un );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_un_z( uo, un, *l );
//      validator.derivative_check( aJl, aJ, update_z, *c, *vz, z, "H_un_z");
//    } else os << "\napplyAdjointHessian_uo_uo not implemented.\n";
//
//
//    if( &Derived::applyAdjointHessian_z_uo != &Base::applyAdjointHessian_z_uo ) {
//      auto aJ  = con_check.adjointJacobian_z( uo, un );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_z_un( uo, z, *l );
//      validator.derivative_check( aJl, aJ, update_un, *c, *vu, un, "H_z_un");
//    } else os << "\napplyAdjointHessian_z_uo not implemented.\n";
//
//    if( &Derived::applyAdjointHessian_z_un != &Base::applyAdjointHessian_z_uo ) {
//      auto aJ  = con_check.adjointJacobian_z( uo, un );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_z_un( uo, z, *l );
//      validator.derivative_check( aJl, aJ, update_un, *c, *vu, un, "H_z_un");
//    } else os << "\napplyAdjointHessian_z_un not implemented.\n";
//
//    if( &Derived::applyAdjointHessian_z_z != &Base::applyAdjointHessian_z_z ) {
//      auto aJ  = con_check.adjointJacobian_z( uo, un );
//      auto aJl = fix_direction( aJ, *l );
//      auto aH  = con_check.adjointHessian_z_z( uo, un, *l );
//      validator.derivative_check( aJl, aJ, update_z, *c, *vz, z, "H_z_un");
//    } else os << "\napplyAdjointHessian_z_z not implemented.\n";

  } // check()

  
};

} // namespace ROL

#endif // ROL_DYNAMICCONSTRAINTCHECK_HPP


