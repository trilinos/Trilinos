// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_Stepper_ErrorNorm_impl_hpp
#define Tempus_Stepper_ErrorNorm_impl_hpp

#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolve_decl.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorStdOps_decl.hpp"
namespace Tempus {

template<class Scalar>
Stepper_ErrorNorm<Scalar>::Stepper_ErrorNorm() : relTol_(1.0e-4), abssTol_(1.0e-4)
{}

template<class Scalar>
Stepper_ErrorNorm<Scalar>::Stepper_ErrorNorm(const Scalar relTol, const Scalar absTol) :
  relTol_(relTol), abssTol_(absTol)
{}

template<class Scalar>
Scalar Stepper_ErrorNorm<Scalar>::
computeErrorNorm(const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x, 
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &err)
{
  if (errorWeightVector_ == Teuchos::null)
    errorWeightVector_ = Thyra::createMember(x->space());

  Thyra::assign(errorWeightVector_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *err, *x, errorWeightVector_.ptr());

  const auto space_dim = err->space()->dim();
  Scalar err_norm = std::abs( Thyra::norm(*errorWeightVector_) / space_dim);
  return err_norm;
}


template<class Scalar>
Scalar Stepper_ErrorNorm<Scalar>::
errorNorm(const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x)
{
  if (scratchVector_ == Teuchos::null) 
    scratchVector_  = Thyra::createMember(x->space());

  Thyra::assign(scratchVector_.ptr(), *x); // | U |
  Thyra::abs(*x, scratchVector_.ptr());
  Thyra::Vt_S(scratchVector_.ptr(), relTol_);
  Thyra::Vp_S(scratchVector_.ptr(), abssTol_);
  Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *x, *scratchVector_, scratchVector_.ptr());
  Scalar err = Thyra::norm_inf(*scratchVector_);
  return err;
  
}


}

#endif
