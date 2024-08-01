// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_Simple2DModelEvaluator.hpp"
#include "Thyra_DummyTestModelEvaluator.hpp"
#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


using Teuchos::null;
typedef ModelEvaluatorBase MEB;


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( StepSize_checker, step_size , Scalar )
{
  RCP<DummyTestModelEvaluator<Scalar> > model = dummyTestModelEvaluator<Scalar>();
  TEST_ASSERT(model != null);
  TEST_EQUALITY(model->Np(), 0);
  TEST_EQUALITY(model->Ng(), 0);
  TEST_ASSERT(model->get_x_space() != null);
  TEST_EQUALITY(model->get_x_space()->dim(), 2);
  TEST_ASSERT(model->get_f_space() != null);
  TEST_EQUALITY(model->get_f_space()->dim(), 2);
  // ToDo: Test getNominalValues()
  TEST_ASSERT(model->create_W_op() != null);
  TEST_ASSERT(model->get_W_factory() != null);
  TEST_ASSERT(model->create_W_prec() != null);
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_step_size));
  TEST_EQUALITY(inArgs.Np(), 0);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  StepSize_checker, step_size)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( StepSize_checker, stageNumber , Scalar )
{
  RCP<DummyTestModelEvaluator<Scalar> > model = dummyTestModelEvaluator<Scalar>();
  TEST_ASSERT(model != null);
  TEST_EQUALITY(model->Np(), 0);
  TEST_EQUALITY(model->Ng(), 0);
  TEST_ASSERT(model->get_x_space() != null);
  TEST_EQUALITY(model->get_x_space()->dim(), 2);
  TEST_ASSERT(model->get_f_space() != null);
  TEST_EQUALITY(model->get_f_space()->dim(), 2);
  // ToDo: Test getNominalValues()
  TEST_ASSERT(model->create_W_op() != null);
  TEST_ASSERT(model->get_W_factory() != null);
  TEST_ASSERT(model->create_W_prec() != null);
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_stage_number));
  TEST_EQUALITY(inArgs.Np(), 0);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  StepSize_checker, stageNumber)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( StepSize_checker, stepSizeAndStageNumber , Scalar )
{
  RCP<DummyTestModelEvaluator<Scalar> > model = dummyTestModelEvaluator<Scalar>();
  TEST_ASSERT(model != null);
  TEST_EQUALITY(model->Np(), 0);
  TEST_EQUALITY(model->Ng(), 0);
  TEST_ASSERT(model->get_x_space() != null);
  TEST_EQUALITY(model->get_x_space()->dim(), 2);
  TEST_ASSERT(model->get_f_space() != null);
  TEST_EQUALITY(model->get_f_space()->dim(), 2);
  // ToDo: Test getNominalValues()
  TEST_ASSERT(model->create_W_op() != null);
  TEST_ASSERT(model->get_W_factory() != null);
  TEST_ASSERT(model->create_W_prec() != null);
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_step_size));
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_stage_number));
  TEST_EQUALITY(inArgs.Np(), 0);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  StepSize_checker, stepSizeAndStageNumber)

} // namespace Thyra
