// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_ModelEvaluatorBase.hpp"
#include "Thyra_DummyTestModelEvaluator.hpp"

namespace Thyra {

  // In TEUCHOS_DEBUG builds, an exception is thrown in
  // ModelEvaluatorDefaultBase runtime checking if Np() on the inArgs
  // and outArgs is not the same. Lazy instantiation normally sets
  // Np() on the model evaluator and underlying inArgs and outArgs. If
  // this parameter changes in the model evaluator, the code is
  // expected to manually call a new initialization. In some cases we
  // want the derived ME to allow for lazy instatiaion since
  // initialization could be costly. So the function resetDefaultBase
  // resets the lazyInitializaion flag such that DefaultBase will
  // automatically re-initialize at the right time and only once.
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(ModelEvaluatorDefaultBase,resetDefaultBase,Scalar)
  {
    const Ordinal x_size = 4;
    const Ordinal p_size = 2;
    const RCP<Thyra::DummyTestModelEvaluator<Scalar> > model = 
      Thyra::dummyTestModelEvaluator<Scalar>(x_size,Teuchos::tuple<Ordinal>(p_size,p_size));

    TEST_EQUALITY(model->createInArgs().Np(),2);
    TEST_EQUALITY(model->createOutArgs().Np(),2);
    TEST_EQUALITY(model->Np(),2);

    // Lazy instantiation is already called so would normally not
    // reset the ModelEvaluatorDefaultBase to Np() == 3 and doesn't
    // seem to resize the outArgs structures correctly.  This tests
    // that by adding resetDefaultBase() in set_p() function of the
    // model evaluator results in a new lazyInstantion that makes Np
    // consistent across all objects.
    model->change_p_size_incorrectly(3);
    TEST_EQUALITY(model->createInArgs().Np(),3);
    TEST_EQUALITY(model->createOutArgs().Np(),2);
    TEST_EQUALITY(model->Np(),2);

    model->change_p_size_correctly(3);
    TEST_EQUALITY(model->createInArgs().Np(),3);
    TEST_EQUALITY(model->createOutArgs().Np(),3);
    TEST_EQUALITY(model->Np(),3);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(ModelEvaluatorDefaultBase,resetDefaultBase)

}
