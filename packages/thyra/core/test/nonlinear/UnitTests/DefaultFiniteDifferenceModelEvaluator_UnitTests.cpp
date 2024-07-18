// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Thyra::DefaultFiniteDifferenceModelEvaluator;
using Thyra::defaultFiniteDifferenceModelEvaluator;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultFiniteDifferenceModelEvaluator, defaultConstruct, Scalar )
{

  RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> > fdModel =
    defaultFiniteDifferenceModelEvaluator<Scalar>();

  TEST_EQUALITY( fdModel->get_direcFiniteDiffCalculator(), null );
  // ToDo: Add more tests!

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( DefaultFiniteDifferenceModelEvaluator, defaultConstruct )


} // namespace
