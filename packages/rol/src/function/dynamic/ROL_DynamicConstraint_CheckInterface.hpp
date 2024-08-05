// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICCONSTRAINT_CHECKINTERFACE_HPP
#define ROL_DYNAMICCONSTRAINT_CHECKINTERFACE_HPP

#include "ROL_FunctionBindings.hpp"
#include "ROL_DynamicConstraint.hpp"

namespace ROL {
namespace details {

using namespace std;
namespace ph = std::placeholders;

template<typename Real>
class DynamicConstraint_CheckInterface {
private:

  using V  = Vector<Real>;
  using Con = DynamicConstraint<Real>;

  Con& con_;
  Real tol_;
  TimeStamp<Real> ts_;

public:

  DynamicConstraint_CheckInterface( Con& con ) : con_(con) { 
    ts_.t.resize(2);
    ts_.t.at(0) = 0.01;
    ts_.t.at(1) = 0.02345;
    ts_.k = 0;
  }

  DynamicConstraint_CheckInterface( Con& con, TimeStamp<Real>& ts ) : con_(con), ts_(ts) { 
  }
 
  f_update_t<Real> update_uo( const V& un, const V& z ) {
    return bind( &Con::update, &con_, ph::_1, cref(un), cref(z), ts_ );
  }

  f_update_t<Real> update_un( const V& uo, const V& z ) {
    return bind( &Con::update, &con_, cref(uo), ph::_1, cref(z), ts_ );
  }

  f_update_t<Real> update_z( const V& uo, const V& un ) {
    return bind( &Con::update, &con_, cref(uo), cref(un), ph::_1, ts_ );
  }

  //----------------------------------------------------------------------------

  f_vector_t<Real> value_uo( const V& un, const V& z ) {
    return bind( &Con::value, &con_, 
                 ph::_1, ph::_2, cref(un), cref(z), ts_ );
  }

  f_vector_t<Real> value_un( const V& uo, const V& z ) {
    return bind( &Con::value, &con_, 
                 ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  f_vector_t<Real> value_z( const V& uo, const V& un ) {
    return bind( &Con::value, &con_, 
                 ph::_1, cref(uo), cref(un), ph::_2, ts_ );
  }

  f_solve_t<Real> solve_un( const V& uo, const V& z ) {
    return bind( &Con::solve, &con_,
                 ph::_1, cref(uo), ph::_2, cref(z), ts_ );
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> jacobian_uo( const V& un, const V& z ) {
    return bind( &Con::applyJacobian_uo, &con_, ph::_1, ph::_2, ph::_3, 
                 cref(un), cref(z), ts_ ); 
  }

  f_dderiv_t<Real> jacobian_un( const V& uo, const V& z ) {
    return bind( &Con::applyJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }

  f_dderiv_t<Real> inverseJacobian_un( const V& uo, const V& z  ) {
    return bind( &Con::applyInverseJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }

  f_dderiv_t<Real> jacobian_z( const V& uo, const V& un ) {
    return bind( &Con::applyJacobian_z, &con_, ph::_1, ph::_2, cref(uo), 
                 cref(un), ph::_3, ts_ ); 
  }

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> adjointJacobian_uo( const V& un, const V& z ) {
    return bind( &Con::applyAdjointJacobian_uo, &con_, ph::_1, ph::_2, ph::_3,
                 cref(un), cref(z), ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_un( const V& uo, const V& z ) {
    return bind( &Con::applyAdjointJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }

  f_dderiv_t<Real> inverseAdjointJacobian_un( const V& uo, const V& z  ) {
    return bind( &Con::applyInverseAdjointJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }
 
  f_dderiv_t<Real> adjointJacobian_z( const V& uo, const V& un ) {
    return bind( &Con::applyAdjointJacobian_z, &con_, ph::_1, ph::_2, cref(uo), 
                 cref(un), ph::_3, ts_ ); 
  }
 
  //----------------------------------------------------------------------------

  f_dderiv_t<Real> adjointJacobian_uo_uo( const V& un, const V& z ) {
    return bind( &Con::applyAdjointJacobian_uo, &con_, ph::_1, ph::_2, ph::_3,
                 cref(un), cref(z), ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_uo_un( const V& uo, const V& z ) {
    return bind( &Con::applyAdjointJacobian_uo, &con_, ph::_1, ph::_2, cref(uo),
                 ph::_3, cref(z), ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_uo_z( const V& uo, const V& un ) {
    return bind( &Con::applyAdjointJacobian_uo, &con_, ph::_1, ph::_2, cref(uo), 
                 cref(un), ph::_3, ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_un_uo( const V& un, const V& z ) {
    return bind( &Con::applyAdjointJacobian_un, &con_, ph::_1, ph::_2, ph::_3,
                 cref(un), cref(z), ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_un_un( const V& uo, const V& z ) {
    return bind( &Con::applyAdjointJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }

  f_dderiv_t<Real> adjointJacobian_un_z( const V& uo, const V& un ) {
    return bind( &Con::applyAdjointJacobian_un, &con_, ph::_1, ph::_2, cref(uo), 
                 cref(un), ph::_3, ts_ ); 
  }
 
  f_dderiv_t<Real> adjointJacobian_z_uo( const V& un, const V& z ) {
    return bind( &Con::applyAdjointJacobian_z, &con_, ph::_1, ph::_2, ph::_3,
                 cref(un), cref(z), ts_ ); 
  }
 
  f_dderiv_t<Real> adjointJacobian_z_un( const V& uo, const V& z ) {
    return bind( &Con::applyAdjointJacobian_z, &con_, ph::_1, ph::_2, cref(uo), 
                 ph::_3, cref(z), ts_ ); 
  }
 
  f_dderiv_t<Real> adjointJacobian_z_z( const V& uo, const V& un ) {
    return bind( &Con::applyAdjointJacobian_z, &con_, ph::_1, ph::_2, cref(uo), 
                 cref(un), ph::_3, ts_ ); 
  }
 
  //----------------------------------------------------------------------------

  f_dderiv_t<Real> adjointHessian_un_un( const V& uo, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_un_un, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), ph::_3, cref(z), ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_un_uo( const V& uo, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_un_uo, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), ph::_3, cref(z), ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_un_z( const V& uo, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_un_z, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), ph::_3, cref(z), ts_ );
  } 

  //----------------------------------------------------------------------------
  
  f_dderiv_t<Real> adjointHessian_uo_un( const V& un, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_uo_un, &con_, ph::_1, cref(l), ph::_2, 
                 ph::_3, cref(un), cref(z), ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_uo_uo( const V& un, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_uo_uo, &con_, ph::_1, cref(l), ph::_2, 
                 ph::_3, cref(un), cref(z), ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_uo_z( const V& un, const V& z, const V& l ) {
    return bind( &Con::applyAdjointHessian_uo_z, &con_, ph::_1, cref(l), ph::_2, 
                 ph::_3, cref(un), cref(z), ts_ );
  } 

  //----------------------------------------------------------------------------

  f_dderiv_t<Real> adjointHessian_z_un( const V& uo, const V& un, const V& l ) {
    return bind( &Con::applyAdjointHessian_z_un, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), cref(un), ph::_3, ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_z_uo( const V& uo, const V& un, const V& l ) {
    return bind( &Con::applyAdjointHessian_z_uo, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), cref(un), ph::_3, ts_ );
  } 

  f_dderiv_t<Real> adjointHessian_z_z( const V& uo, const V& un, const V& l ) {
    return bind( &Con::applyAdjointHessian_z_z, &con_, ph::_1, cref(l), ph::_2, 
                 cref(uo), cref(un), ph::_3, ts_ );
  } 

}; // class DynamicConstraint_CheckInterface

} // namespace details 

using details::DynamicConstraint_CheckInterface;

template<typename Real>
DynamicConstraint_CheckInterface<Real> make_check( DynamicConstraint<Real>& con ) {
  return DynamicConstraint_CheckInterface<Real>(con);
}

template<typename Real>
DynamicConstraint_CheckInterface<Real> make_check( DynamicConstraint<Real>& con, 
                                                   TimeStamp<Real>& timeStamp ) {
  return DynamicConstraint_CheckInterface<Real>(con,timeStamp);
}

} // namespace ROL


#endif // ROL_DYNAMICCONSTRAINT_CHECKINTERFACE_HPP


