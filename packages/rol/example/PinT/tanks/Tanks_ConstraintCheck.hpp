// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKS_CONSTRAINTCHECK_HPP
#define TANKS_CONSTRAINTCHECK_HPP

#include <ostream>

#include "ROL_ValidateFunction.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_DynamicConstraint_CheckInterface.hpp"
#include "ROL_SerialConstraint.hpp"

namespace Tanks {

template<typename Real>
void check( DynamicConstraint<Real>& con, 
            ROL::BoundConstraint<Real>& ubnd, 
            ROL::BoundConstraint<Real>& zbnd,
            ROL::ValidateFunction<Real>& validator ) { 
          
  auto uo = ubnd.getLowerBound()->clone();
  auto un = ubnd.getLowerBound()->clone();
  auto vu = uo->clone();
  auto z  = zbnd.getLowerBound()->clone();
  auto vz = z->clone();
  auto c  = uo->clone();
  auto l  = c->dual().clone();

  ROL::RandomizeFeasibleVector(*un, ubnd);
  ROL::RandomizeFeasibleVector(*uo, ubnd);
  ROL::RandomizeVector( *vu );
  ROL::RandomizeVector( *z  );
  ROL::RandomizeVector( *vz );
  ROL::RandomizeVector( *c  );
  ROL::RandomizeVector( *l  );

  auto con_check = ROL::make_check( con );
  auto up_uo = con_check.update_uo(*un,*z);
  auto up_un = con_check.update_un(*uo,*z);
  auto up_z  = con_check.update_z(*uo,*un);

  auto val_uo = con_check.value_uo( *un, *z  );
  auto val_un = con_check.value_un( *uo, *z  );
  auto val_z  = con_check.value_z(  *uo, *un );

  auto J_uo = con_check.jacobian_uo( *un, *z  );
  auto J_un = con_check.jacobian_un( *uo, *z  );
  auto J_z  = con_check.jacobian_z(  *uo, *un );

  auto iJ_un = con_check.inverseJacobian_un( *uo, *z );

  auto aJ_uo = con_check.adjointJacobian_uo( *un, *z  );
  auto aJ_un = con_check.adjointJacobian_un( *uo, *z  );
  auto aJ_z  = con_check.adjointJacobian_z(  *uo, *un );
  
  auto aJ_uol = fix_direction( aJ_uo, *l );
  auto aJ_unl = fix_direction( aJ_un, *l );
  auto aJ_zl  = fix_direction( aJ_z,  *l );
  
  auto iaJ_un = con_check.inverseAdjointJacobian_un( *uo, *z );

  // Check uo,un,z Jacobians
  validator.derivative_check( val_uo, J_uo, up_uo, *c, *vu, *uo, "norm(J_uo*vec)" );
  validator.derivative_check( val_un, J_un, up_un, *c, *vu, *un, "norm(J_un*vec)" );
  validator.derivative_check( val_z,  J_z,  up_z,  *c, *vz, *z,  "norm(J_z*vec)" );

  // Check un inverse Jacobian
  validator.inverse_check( J_un, iJ_un, up_un, *vu, *un, "Jacobian with respect to un", "J_un");
 
  // Check Jacobian adjoint consistencies
  validator.adjoint_consistency_check( J_uo, aJ_uo, up_uo, *c, *vu, *uo, "Jacobian with respect to uo", "J_uo");
  validator.adjoint_consistency_check( J_un, aJ_un, up_un, *c, *vu, *un, "Jacobian with respect to un", "J_un");
  validator.adjoint_consistency_check( J_z,  aJ_z,  up_z,  *vz, *c, *z,  "Jacobian with respect to z", "J_z");

} // check(ROL::DynamicConstraint)



template<typename Real>
void check( ROL::SerialConstraint<Real>& con, 
            ROL::BoundConstraint<Real>& ubnd, // For a single time step
            ROL::BoundConstraint<Real>& zbnd,
            std::ostream& outStream ) { 

    auto Nt = con.numTimeSteps();
    auto up = ubnd.getLowerBound();
    auto zp = zbnd.getLowerBound();

    auto u_serial  = ROL::PartitionedVector<Real>::create( *up, Nt );
    auto vu_serial = ROL::PartitionedVector<Real>::create( *up, Nt );
    auto c_serial  = ROL::PartitionedVector<Real>::create( *up, Nt );
    auto z_serial  = ROL::PartitionedVector<Real>::create( *zp, Nt );
    auto vz_serial = ROL::PartitionedVector<Real>::create( *zp, Nt );

    for( size_type i=0; i<Nt; ++i ) {
      ROL::RandomizeFeasibleVector( *(u_serial->get(i)), ubnd );
      ROL::RandomizeVector( *( c_serial->get(i)) );
      ROL::RandomizeVector( *(vu_serial->get(i)) );
      ROL::RandomizeVector( *( z_serial->get(i)) );
      ROL::RandomizeVector( *(vz_serial->get(i)) );
    }

    con.checkApplyJacobian_1(               *u_serial,  *z_serial,  *vu_serial, *c_serial, true, outStream );
    con.checkApplyJacobian_2(               *u_serial,  *z_serial,  *vz_serial, *c_serial, true, outStream );
    con.checkInverseJacobian_1( *c_serial, *vu_serial,  *u_serial,   *z_serial,            true, outStream );
    con.checkAdjointConsistencyJacobian_1(  *c_serial,  *vu_serial,  *u_serial, *z_serial, true, outStream );
    con.checkAdjointConsistencyJacobian_2(  *c_serial,  *vz_serial,  *u_serial, *z_serial, true, outStream );
    con.checkInverseAdjointJacobian_1(     * c_serial,  *vu_serial,  *u_serial, *z_serial, true, outStream );
} // check(Tanks::SerialConstraint)

} // namespace Tanks


#endif // TANKS_CONSTRAINTCHECK_HPP
