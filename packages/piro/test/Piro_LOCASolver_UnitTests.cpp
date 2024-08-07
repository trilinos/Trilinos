// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_LOCASolver.hpp"

#include "MockModelEval_A.hpp"

#include "Piro_Test_WeakenedModelEvaluator.hpp"
#include "Piro_Test_ThyraSupport.hpp"
#include "Piro_Test_MockObserver.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include <stdexcept>

using namespace Teuchos;
using namespace Piro;
using namespace Piro::Test;

namespace Thyra {
  typedef ModelEvaluatorBase MEB;
} // namespace Thyra

// Setup support

const RCP<EpetraExt::ModelEvaluator> epetraModelNew()
{
#ifdef HAVE_MPI
  const MPI_Comm comm = MPI_COMM_WORLD;
#else /*HAVE_MPI*/
  const int comm = 0;
#endif /*HAVE_MPI*/
  return rcp(new MockModelEval_A(comm));
}

const RCP<Thyra::ModelEvaluator<double> > thyraModelNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory(new Thyra::AmesosLinearOpWithSolveFactory);
  return epetraModelEvaluator(epetraModel, lowsFactory);
}

const RCP<Thyra::ModelEvaluator<double> > solverNew(
    const RCP<Thyra::ModelEvaluator<double> > &thyraModel,
    const RCP<Thyra::ModelEvaluator<double> > &thyraAdjointModel = Teuchos::null,
    const RCP<Piro::ObserverBase<double> > &observer = Teuchos::null)
{
  const RCP<ParameterList> piroParams(new ParameterList("Piro Parameters"));
  updateParametersFromXmlFile("input_Solve_LOCA_1.xml", piroParams.ptr());
  return observedLocaSolver<double>(piroParams, thyraModel, thyraAdjointModel, observer);
}

const RCP<Thyra::ModelEvaluator<double> > solverNew(
    const RCP<EpetraExt::ModelEvaluator> &epetraModel,
    const RCP<Piro::ObserverBase<double> > &observer = Teuchos::null)
{
  return solverNew(thyraModelNew(epetraModel), Teuchos::null, observer);
}

// Floating point tolerance
const double tol = 1.0e-8;

// Tests

TEUCHOS_UNIT_TEST(Piro_LOCASolver, Spaces)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  TEST_ASSERT(solver->Np() == 1);
  TEST_ASSERT(solver->Ng() == 2);

  const int parameterIndex = 0;
  const int responseIndex = 0;
  const int solutionResponseIndex = solver->Ng() - 1;
  TEST_ASSERT(nonnull(solver->get_p_space(parameterIndex)));
  TEST_ASSERT(nonnull(solver->get_g_space(responseIndex)));
  TEST_ASSERT(nonnull(solver->get_g_space(solutionResponseIndex)));

  // TODO
  //TEST_THROW(solver->get_x_space(), std::exception);
  //TEST_THROW(solver->get_f_space(), std::exception);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, Solution)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionObserver)
{
  const RCP<MockObserver<double> > observer(new MockObserver<double>);
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew(), observer);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*observer->lastSolution());
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionForMissingParameterValues)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->createInArgs();
  const int parameterIndex = 0;
  inArgs.set_p(parameterIndex, null);

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionForAlternateParameterValues)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  {
    const int parameterIndex = 0;
    const RCP<Thyra::VectorBase<double> > p_in = Thyra::createMember(*solver->get_p_space(0));
    TEST_EQUALITY(p_in->space()->dim(), 2);
    Thyra::set_ele(0, 1.0, p_in.ptr());
    Thyra::set_ele(1, 0.0, p_in.ptr());
    inArgs.set_p(parameterIndex, p_in);
  }

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 1.0, 2.0, 3.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, Response)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(8.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, ResponseForMissingParameterValues)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->createInArgs();
  const int parameterIndex = 0;
  inArgs.set_p(parameterIndex, null);

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(8.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, ResponseForAlternateParameterValues)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  {
    const int parameterIndex = 0;
    const RCP<Thyra::VectorBase<double> > p_in = Thyra::createMember(*solver->get_p_space(0));
    TEST_EQUALITY(p_in->space()->dim(), 2);
    Thyra::set_ele(0, 1.0, p_in.ptr());
    Thyra::set_ele(1, 0.0, p_in.ptr());
    inArgs.set_p(parameterIndex, p_in);
  }

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(18.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionSensitivityMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionSensitivityOp)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));

  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromLinOp(*dxdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionSensitivityOp_NoDfDpMv)
{
  // Disable support for MultiVector-based DfDp derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDfDpMv(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));

  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  solver->evalModel(inArgs, outArgs);

  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromLinOp(*dxdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvGrad)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp->col(parameterIndex));
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityOp)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromLinOp(*dgdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvJac_NoDgDxMv)
{
  // Disable support for MultiVector-based DgDx derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDxMv(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvGrad_NoDgDpMvJac)
{
  // Disable support for Jacobian-oriented MultiVector DgDx derivative
  // (Only gradient layout is available)
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMvJac(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityOp_NoDgDpMv)
{
  // Disable support for MultiVector-based DgDp derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMv(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromLinOp(*dgdp, i);
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvJacWithResponseSensitivityMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int parameterIndex = 0;

  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  const int responseIndex = 0;
  const Thyra::MEB::Derivative<double> dgdp_deriv =
    Thyra::create_DgDp_mv(*solver, responseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dgdp = dgdp_deriv.getMultiVector();
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromVector(*dxdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      const Array<double> actual = arrayFromVector(*dgdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvGradWithSolutionSensitivityOp_NoDgDpMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMvJac(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (LINOP layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const RCP<Thyra::LinearOpBase<double> > dxdp =
    solver->create_DgDp_op(solutionResponseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dxdp));
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp);

  // Request response sensitivity (MV_GRAD layout)
  const int responseIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromLinOp(*dxdp, i);
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Verifiy response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    const Array<double> actual = arrayFromVector(*dgdp);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvGradWithSolutionSensitivityMvJac_NoDgDpMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
    rcp(new WeakenedModelEvaluator_NoDgDpMvJac(thyraModelNew(epetraModelNew())));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (MV_JAC layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  // Request response sensitivity (MV_GRAD layout)
  const int responseIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromVector(*dxdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Verifiy response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    const Array<double> actual = arrayFromVector(*dgdp);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityOpWithSolutionSensitivityMvJac)
{
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (MV_JAC layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  // Request response sensitivity (LINOP layout)
  const int responseIndex = 0;
  const RCP<Thyra::LinearOpBase<double> > dgdp =
    solver->create_DgDp_op(responseIndex, parameterIndex);
  TEST_ASSERT(nonnull(dgdp));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  const Array<Array<double> > expected = tuple(
      Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
      Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
  TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
  TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
    const Array<double> actual = arrayFromVector(*dxdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
  }

  // Verify response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    for (int i = 0; i < expected.size(); ++i) {
      const Array<double> actual = arrayFromLinOp(*dgdp, i);
      TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
    }
  }
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SensitivityMvGradWithSolutionSensitivityMvJac_NoAdjointW_NoDgDpMvJac)
{
  // Disable support for Jacobian adjoint solve
  // (Only forward solve is available)
  const RCP<Thyra::ModelEvaluator<double> > weakenedModel =
      rcp(new WeakenedModelEvaluator_NoDgDpMvJac(rcp(new WeakenedModelEvaluator_NoAdjointW(thyraModelNew(epetraModelNew())))));
  const RCP<Thyra::ModelEvaluator<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();

  const int parameterIndex = 0;

  // Request solution sensitivity (MV_JAC layout)
  const int solutionResponseIndex = solver->Ng() - 1;
  const Thyra::MEB::Derivative<double> dxdp_deriv =
    Thyra::create_DgDp_mv(*solver, solutionResponseIndex, parameterIndex, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  const RCP<Thyra::MultiVectorBase<double> > dxdp = dxdp_deriv.getMultiVector();
  outArgs.set_DgDp(solutionResponseIndex, parameterIndex, dxdp_deriv);

  // Request response sensitivity (MV_GRAD layout)
  const int responseIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  // Verify solution sensitivity
  {
    const Array<Array<double> > expected = tuple(
        Array<double>(tuple(0.5, 0.0, 0.0, 0.0)),
        Array<double>(tuple(0.0, 1.0, 1.0, 1.0)));
    TEST_EQUALITY(dxdp->domain()->dim(), expected.size());
    for (int i = 0; i < expected.size(); ++i) {
      TEST_EQUALITY(dxdp->range()->dim(), expected[i].size());
      const Array<double> actual = arrayFromVector(*dxdp->col(i));
      TEST_COMPARE_FLOATING_ARRAYS(actual, expected[i], tol);
    }
  }

  // Verifiy response sensitivity
  {
    const Array<double> expected = tuple(2.0, -8.0);
    const Array<double> actual = arrayFromVector(*dgdp);
    TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
  }
}

#endif /*HAVE_PIRO_NOX*/
