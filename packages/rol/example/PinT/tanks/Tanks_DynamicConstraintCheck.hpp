#pragma once
#ifndef TANKS_DYNAMICCONSTRAINTCHECK_HPP
#define TANKS_DYNAMICCONSTRAINTCHECK_HPP

#include <ostream>

#include "ROL_ValidateFunction.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_DynamicConstraint_CheckInterface.hpp"

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
  auto up_uo = con_check.update_uo();
  auto up_un = con_check.update_un();
  auto up_z  = con_check.update_z();

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

} // check

} // namespace Tanks


#endif // TANKS_DYNAMICCONSTRAINTCHECK_HPP
