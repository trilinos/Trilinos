//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "Tempus_NumericalUtils.hpp"
#include "Tempus_Stepper_ErrorNorm.hpp"

#include "../TestModels/DahlquistTestModel.hpp"

namespace Tempus_Unit_Test {

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Stepper_ErrorNorm, computeWRMSNorm)
{
  int N = 3;
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > xSpace =
      Thyra::defaultSpmdVectorSpace<double>(N);
  auto x     = Thyra::createMember(xSpace);
  auto xNext = Thyra::createMember(xSpace);
  auto eVec  = Thyra::createMember(xSpace);

  auto tol   = Tempus::numericalTol<double>();
  auto eNorm = Teuchos::rcp(new Tempus::Stepper_ErrorNorm<double>(tol, 0.0));

  Thyra::assign(x.ptr(), 0.0);
  Thyra::assign(xNext.ptr(), 0.0);
  Thyra::assign(eVec.ptr(), 0.0);
  double rmsNorm = eNorm->computeWRMSNorm(x, xNext, eVec);
  TEST_FLOATING_EQUALITY(rmsNorm, 0.0, 1.0e-12);

  Thyra::assign(x.ptr(), 1.0);
  Thyra::assign(xNext.ptr(), 1.0 + tol);
  Thyra::assign(eVec.ptr(), tol);
  rmsNorm = eNorm->computeWRMSNorm(x, xNext, eVec);
  TEST_FLOATING_EQUALITY(rmsNorm, 1.0 / std::sqrt(N), 1.0e-12);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Stepper_ErrorNorm, errorNorm)
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > xSpace =
      Thyra::defaultSpmdVectorSpace<double>(3);
  auto x = Thyra::createMember(xSpace);

  auto tol   = Tempus::numericalTol<double>();
  auto eNorm = Teuchos::rcp(new Tempus::Stepper_ErrorNorm<double>(tol, tol));

  Thyra::assign(x.ptr(), 0.0);
  double norm = eNorm->errorNorm(x);
  TEST_FLOATING_EQUALITY(norm, tol, 1.0e-12);

  Thyra::assign(x.ptr(), 1.0);
  norm = eNorm->errorNorm(x);
  TEST_FLOATING_EQUALITY(norm, 1.0 / (2.0 * tol), 1.0e-12);
}

}  // namespace Tempus_Unit_Test
