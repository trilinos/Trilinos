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
#ifndef ROL_DYNAMICOBJECTIVECHECK_HPP
#define ROL_DYNAMICOBJECTIVECHECK_HPP

#include "ROL_DynamicObjective_CheckInterface.hpp"
#include "ROL_ValidateFunction.hpp"
#include "ROL_RandomVector.hpp"
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

    auto gu = uo.dual().clone();
    auto gz = z.dual().clone();
    auto vu = uo.clone();
    auto vz = z.clone();
    
    RandomizeVector(*vu);
    RandomizeVector(*vz);

    auto obj_check = make_check( obj );

    auto update_uo = obj_check.update_uo( un, z );
    auto update_un = obj_check.update_un( uo, z );
    auto update_z  = obj_check.update_z( un, uo );
  
    auto value_uo = obj_check.value_uo( un, z );
    auto value_un = obj_check.value_un( uo, z );
    auto value_z  = obj_check.value_z( uo, un );

    auto grad_uo = obj_check.gradient_uo( uo, z  );
    auto grad_un = obj_check.gradient_un( un, z  );
    auto grad_z  = obj_check.gradient_z(  uo, un );

    //-------------------------------------------------------------------------
    // Check gradient components
    if( std::find(methods.begin(),methods.end(),"gradient_uo") != methods.end() ) {
      validator.derivative_check( value_uo, grad_uo, update_uo, *gu, *vu, uo, "grad_uo'*dir" );
    }             
    if( std::find(methods.begin(),methods.end(),"gradient_un") != methods.end() ) {
      validator.derivative_check( value_un, grad_un, update_un, *gu, *vu, un, "grad_un'*dir" );
    }             
    if( std::find(methods.begin(),methods.end(),"gradient_z") != methods.end() ) {
      validator.derivative_check( value_z, grad_z, update_z, *gz, *vz, z, "grad_z'*dir" );
    }             

    //-------------------------------------------------------------------------
    // Check Hessian components
    if( std::find(methods.begin(),methods.end(),"hessVec_uo_uo") != methods.end() ) {
      auto H = obj_check.hessVec_uo_uo(un,z);
      validator.derivative_check( grad_uo, H, update_uo, *gu, *vu, uo, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_un") != methods.end() ) {
      auto H = obj_check.hessVec_uo_un(uo,z);
      validator.derivative_check( grad_uo, H, update_un, *gu, *vu, un, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_z") != methods.end() ) {
      auto H = obj_check.hessVec_uo_z(uo,un);
      validator.derivative_check( grad_uo, H, update_z, *gz, *vz, z, "norm(Hess*vec)" );
    }



    if( std::find(methods.begin(),methods.end(),"hessVec_uo_uo") != methods.end() ) {
      auto H = obj_check.hessVec_un_uo(un,z);
      validator.derivative_check( grad_un, H, update_uo, *gu, *vu, uo, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_un") != methods.end() ) {
      auto H = obj_check.hessVec_un_un(uo,z);
      validator.derivative_check( grad_un, H, update_un, *gu, *vu, un, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_z") != methods.end() ) {
      auto H = obj_check.hessVec_un_z(uo,un);
      validator.derivative_check( grad_un, H, update_z, *gz, *vz, z, "norm(Hess*vec)" );
    }



    if( std::find(methods.begin(),methods.end(),"hessVec_uo_uo") != methods.end() ) {
      auto H = obj_check.hessVec_z_uo(un,z);
      validator.derivative_check( grad_z, H, update_uo, *gu, *vu, uo, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_un") != methods.end() ) {
      auto H = obj_check.hessVec_z_un(uo,z);
      validator.derivative_check( grad_z, H, update_un, *gu, *vu, un, "norm(Hess*vec)" );
    }

    if( std::find(methods.begin(),methods.end(),"hessVec_uo_z") != methods.end() ) {
      auto H = obj_check.hessVec_z_z(uo,un);
      validator.derivative_check( grad_z, H, update_z, *gz, *vz, z, "norm(Hess*vec)" );
    }
  }
}; // DynamicObjectiveCheck


} // namespace ROL

#endif  // ROL_DYNAMICOBJECTIVECHECK_HPP
