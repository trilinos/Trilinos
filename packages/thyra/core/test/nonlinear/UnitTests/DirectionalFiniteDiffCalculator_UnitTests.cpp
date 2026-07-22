// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"


namespace Thyra {


using Teuchos::null;
using Teuchos::as;
using Teuchos::getParameter;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DirectionalFiniteDiffCalculator, defaultConstruct, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  const RCP<DirectionalFiniteDiffCalculator<Scalar> > fdCalc =
    directionalFiniteDiffCalculator<Scalar>();

  TEST_EQUALITY_CONST( fdCalc->fd_method_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO );
  TEST_EQUALITY_CONST( fdCalc->fd_step_select_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size(), as<ScalarMag>(-1.0) );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size_min(), as<ScalarMag>(-1.0) );

  const RCP<const ParameterList> validPL = fdCalc->getValidParameters();
  TEST_ASSERT(nonnull(validPL));
  TEST_EQUALITY_CONST(getParameter<std::string>(*validPL, "FD Method"),
    "order-one");
  TEST_EQUALITY_CONST(getParameter<std::string>(*validPL, "FD Step Select Type"),
    "Absolute");
  TEST_EQUALITY_CONST(getParameter<double>(*validPL, "FD Step Length"),
    as<double>(-1.0));

  // ToDo: Add more tests!

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( DirectionalFiniteDiffCalculator, defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DirectionalFiniteDiffCalculator, plConstruct, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  const RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("FD Method", "order-two-central");
  pl->set("FD Step Select Type", "Relative");
  pl->set("FD Step Length", as<double>(1e-5));

  const RCP<DirectionalFiniteDiffCalculator<Scalar> > fdCalc =
    directionalFiniteDiffCalculator<Scalar>(pl);

  TEST_EQUALITY_CONST( fdCalc->fd_method_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL );
  TEST_EQUALITY_CONST( fdCalc->fd_step_select_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_STEP_RELATIVE );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size(), as<ScalarMag>(1e-5) );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size_min(), as<ScalarMag>(-1.0) );

  // ToDo: Add more tests!

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( DirectionalFiniteDiffCalculator, plConstruct )


} // namespace Thyra
