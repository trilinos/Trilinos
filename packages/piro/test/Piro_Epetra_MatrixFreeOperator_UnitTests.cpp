// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_MatrixFreeDecorator.hpp"

#include "Piro_Test_EpetraSupport.hpp"
#include "MockModelEval_A.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

using namespace Teuchos;
using namespace Piro;
using namespace Piro::Test;

RCP<EpetraExt::ModelEvaluator> epetraModelNew()
{
#ifdef HAVE_MPI
  const MPI_Comm comm = MPI_COMM_WORLD;
#else /*HAVE_MPI*/
  const int comm = 0;
#endif /*HAVE_MPI*/
  return rcp(new MockModelEval_A(comm));
}

EpetraExt::ModelEvaluator::InArgs createStaticNominalInArgs(const EpetraExt::ModelEvaluator &model)
{
  EpetraExt::ModelEvaluator::InArgs result = model.createInArgs();

  if (nonnull(model.get_x_init())) {
    result.set_x(model.get_x_init());
  }

  const int parameterCount = result.Np();
  for (int l = 0; l < parameterCount; ++l) {
    if (nonnull(model.get_p_init(l))) {
      result.set_p(l, model.get_p_init(l));
    }
  }

  return result;
}

EpetraExt::ModelEvaluator::InArgs createDynamicNominalInArgs(const EpetraExt::ModelEvaluator &model)
{
  EpetraExt::ModelEvaluator::InArgs result = createStaticNominalInArgs(model);

  if (nonnull(model.get_x_dot_init())) {
    result.set_x_dot(model.get_x_init());
  }

  return result;
}

// Floating point tolerance
const double tol = 2.0e-6;

TEUCHOS_UNIT_TEST(Epetra_MatrixFreeOperator, Spaces)
{
  const RCP<EpetraExt::ModelEvaluator> model = epetraModelNew();
  const RCP<Epetra::MatrixFreeOperator> op(new Epetra::MatrixFreeOperator(model));

  TEST_ASSERT(op->OperatorDomainMap().SameAs(*model->get_x_map()));
  TEST_ASSERT(op->OperatorRangeMap().SameAs(*model->get_f_map()));
}

TEUCHOS_UNIT_TEST(Epetra_MatrixFreeOperator, Static)
{
  const RCP<EpetraExt::ModelEvaluator> model = epetraModelNew();

  const EpetraExt::ModelEvaluator::InArgs modelInArgs = createStaticNominalInArgs(*model);

  const RCP<Epetra_Operator> expectedJacobian = model->create_W();
  TEST_ASSERT(nonnull(expectedJacobian));
  const RCP<Epetra_Vector> f_base = vectorNew(*model->get_f_map());
  {
    EpetraExt::ModelEvaluator::OutArgs modelOutArgs = model->createOutArgs();
    modelOutArgs.set_W(expectedJacobian);
    modelOutArgs.set_f(f_base);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const RCP<Epetra::MatrixFreeOperator> jacobian(new Epetra::MatrixFreeOperator(model));
  jacobian->setBase(modelInArgs, f_base, /*haveXdot =*/ false, false);

  TEST_EQUALITY(
      jacobian->OperatorDomainMap().NumGlobalElements(),
      expectedJacobian->OperatorDomainMap().NumGlobalElements());

  const int colCount = expectedJacobian->OperatorDomainMap().NumGlobalElements();
  for (int i = 0; i < colCount; ++i) {
    const Array<double> expected = arrayFromLinOp(*expectedJacobian, i);
    const Array<double> actual = arrayFromLinOp(*jacobian, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Epetra_MatrixFreeOperator, Dynamic)
{
  const RCP<EpetraExt::ModelEvaluator> model = epetraModelNew();

  EpetraExt::ModelEvaluator::InArgs modelInArgs = createDynamicNominalInArgs(*model);
  modelInArgs.set_alpha(2.0);
  modelInArgs.set_beta(0.5);

  const RCP<Epetra_Operator> expectedJacobian = model->create_W();
  const RCP<Epetra_Vector> f_base = vectorNew(*model->get_f_map());
  {
    EpetraExt::ModelEvaluator::OutArgs modelOutArgs = model->createOutArgs();
    modelOutArgs.set_W(expectedJacobian);
    modelOutArgs.set_f(f_base);
    model->evalModel(modelInArgs, modelOutArgs);
  }

  const RCP<Epetra::MatrixFreeOperator> jacobian(new Epetra::MatrixFreeOperator(model));
  jacobian->setBase(modelInArgs, f_base, /*haveXdot =*/ true, false);

  TEST_EQUALITY(
      jacobian->OperatorDomainMap().NumGlobalElements(),
      expectedJacobian->OperatorDomainMap().NumGlobalElements());

  const int colCount = expectedJacobian->OperatorDomainMap().NumGlobalElements();
  for (int i = 0; i < colCount; ++i) {
    const Array<double> expected = arrayFromLinOp(*expectedJacobian, i);
    const Array<double> actual = arrayFromLinOp(*jacobian, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}
