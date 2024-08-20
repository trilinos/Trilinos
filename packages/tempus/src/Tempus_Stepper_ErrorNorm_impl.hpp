//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Stepper_ErrorNorm_impl_hpp
#define Tempus_Stepper_ErrorNorm_impl_hpp

#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolve_decl.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

#include "Tempus_NumericalUtils.hpp"

namespace Tempus {

template <class Scalar>
Stepper_ErrorNorm<Scalar>::Stepper_ErrorNorm()
  : relTol_(1.0e-4), abssTol_(1.0e-4)
{
}

template <class Scalar>
Stepper_ErrorNorm<Scalar>::Stepper_ErrorNorm(const Scalar relTol,
                                             const Scalar absTol)
  : relTol_(relTol), abssTol_(absTol)
{
}

template <class Scalar>
Scalar Stepper_ErrorNorm<Scalar>::computeWRMSNorm(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &xNext,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &err)
{
  if (errorWeightVector_ == Teuchos::null)
    errorWeightVector_ = Thyra::createMember(x->space());

  if (u_ == Teuchos::null) u_ = Thyra::createMember(x->space());

  if (uNext_ == Teuchos::null) uNext_ = Thyra::createMember(x->space());

  // Compute: Atol + max(|u^n|, |u^{n+1}| ) * Rtol
  Thyra::abs(*x, u_.ptr());
  Thyra::abs(*xNext, uNext_.ptr());
  Thyra::pair_wise_max_update(relTol_, *u_, uNext_.ptr());
  if (!approxZero(abssTol_)) {
    Thyra::add_scalar(abssTol_, uNext_.ptr());
  }
  else {
    Scalar absTol = Thyra::norm_2(*uNext_) * numericalTol<Scalar>();
    if (approxZero(absTol)) absTol = numericalTol<Scalar>();
    Thyra::add_scalar(absTol, uNext_.ptr());
  }

  Thyra::assign(errorWeightVector_.ptr(),
                Teuchos::ScalarTraits<Scalar>::zero());
  Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *err, *uNext_,
                         errorWeightVector_.ptr());

  const auto space_dim = err->space()->dim();
  Scalar err_norm      = std::abs(Thyra::norm(*errorWeightVector_) / space_dim);
  return err_norm;
}

template <class Scalar>
Scalar Stepper_ErrorNorm<Scalar>::errorNorm(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x)
{
  if (scratchVector_ == Teuchos::null)
    scratchVector_ = Thyra::createMember(x->space());

  Thyra::assign(scratchVector_.ptr(), *x);  // | U |
  Thyra::abs(*x, scratchVector_.ptr());
  Thyra::Vt_S(scratchVector_.ptr(), relTol_);
  if (!approxZero(abssTol_)) {
    Thyra::Vp_S(scratchVector_.ptr(), abssTol_);
  }
  else {
    Scalar absTol = Thyra::norm_2(*scratchVector_) * numericalTol<Scalar>();
    if (approxZero(absTol)) absTol = numericalTol<Scalar>();
    Thyra::add_scalar(absTol, scratchVector_.ptr());
  }
  Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *x, *scratchVector_,
                         scratchVector_.ptr());
  Scalar err = Thyra::norm_inf(*scratchVector_);
  return err;
}

}  // namespace Tempus

#endif
