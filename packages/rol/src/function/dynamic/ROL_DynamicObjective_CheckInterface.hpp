// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Objtract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following Objditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of Objditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of Objditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// Objtributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// ObjTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR ObjSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN ObjTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Objtact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#pragma once
#ifndef ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP
#define ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP

#include <functional>

#include "ROL_DynamicObjective.hpp"


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

  DynamicObjective_CheckInterface( Obj& obj ) : obj_(obj) {}

  f_update_t<Real> update_uo() {
    return bind( &Obj::update_uo, &obj_, ph::_1 );
  }

  f_update_t<Real> update_un() {
    return bind( &Obj::update_un, &obj_, ph::_1 );
  }

  f_update_t<Real> update_z() {
    return bind( &Obj::update_z, &obj_, ph::_1 );
  }

  //----------------------------------------------------------------------------

  f_vector_t<Real> value_uo( const V& un, const V& z ) {
    return bind( &Obj::value, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> value_un( const V& uo, const V& z ) {
    return bind( &Obj::value, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  f_vector_t<Real> value_z( const V& uo, const V& un ) {
    return bind( &Obj::value, &obj_, ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  //----------------------------------------------------------------------------

  f_vector_t<Real> gradient_uo( const V& un, const V& z ) {
    return bind( &Obj::gradient_uo, &obj_, ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> gradient_un( const V& uo, const V& z ) {
    return bind( &Obj::gradient_un, &obj_, ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }
  
  f_vector_t<Real> gradient_uo( const V& un, const V& z ) {
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
    return bind( &Obj::hessVec_uo_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3 ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> hessVec_un_uo( const V& un, const V& z ) {
    return bind( &Obj::hessVec_un_uo, &obj_, ph::_1, ph::_2, ph::_3, cref(un), cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_un_un( const V& uo, const V& z ) {
    return bind( &Obj::hessVec_un_un, &obj_, ph::_1, ph::_2, cref(uo), ph::_3, cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_un_z( const V& uo, const V& un ) {
    return bind( &Obj::hessVec_un_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3 ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> hessVec_z_uo( const V& un, const V&  ) {
    return bind( &Obj::hessVec_z_uo, &obj_, ph::_1, ph::_2, ph::_3, cref(un), cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_z_un( const V& uo, const V& z ) {
    return bind( &Obj::hessVec_z_un, &obj_, ph::_1, ph::_2, cref(uo), ph::_3, cref(z), ts_ );
  }

  f_dderiv_t<Real> hessVec_z_z( const V& uo, const V& un ) {
    return bind( &Obj::hessVec_z_z, &obj_, ph::_1, ph::_2, cref(uo), cref(un), ph::_3 ts_ );
  }


}; // DynamicObjective_CheckInterface


} // namespace details 

using details::DynamicObjective_CheckInterface;

template<typename Real>
DynamicConstraint_CheckInterface<Real> make_check( DynamicObjective<Real>& obj ) {
  return DynamicConstraint_CheckInterface<Real>( obj );
}


} // namespace ROL


#endif // ROL_DYNAMICOBJECTIVE_CHECKINTERFACE_HPP

