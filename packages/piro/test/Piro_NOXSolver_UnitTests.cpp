/*
// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER
*/

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_NOXSolver.hpp"

#include "MockModelEval_A.hpp"

#include "Piro_Test_WeakenedModelEvaluator.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

using namespace Teuchos;
using namespace Piro;

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

const RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyraModelNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory(new Thyra::AmesosLinearOpWithSolveFactory);
  return epetraModelEvaluator(epetraModel, lowsFactory);
}

const RCP<NOXSolver<double> > solverNew(const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel)
{
  const RCP<Teuchos::ParameterList> piroParams(new Teuchos::ParameterList("Piro Parameters"));
  return rcp(new NOXSolver<double>(piroParams, thyraModel));
}

const RCP<NOXSolver<double> > solverNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  return solverNew(thyraModelNew(epetraModel));
}

// Testing support

Array<double> arrayFromVector(const Thyra::VectorBase<double> &v)
{
  const Thyra::ConstDetachedVectorView<double> view(v, Thyra::Range1D(), true);
  return Array<double>(view.values(), view.values() + view.subDim());
}

const double tol = 1.0e-8;

TEUCHOS_UNIT_TEST(Piro_NOXSolver, Solution)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int solutionResponseIndex = solver->Ng() - 1;
  outArgs.set_g(solutionResponseIndex, Thyra::createMember(*solver->get_g_space(solutionResponseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(solutionResponseIndex));
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, Response)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  outArgs.set_g(responseIndex, Thyra::createMember(*solver->get_g_space(responseIndex)));

  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*outArgs.get_g(responseIndex));
  const Array<double> expected = tuple(8.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvJac)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::MultiVectorBase<double> > dgdp =
    Thyra::createMembers(solver->get_g_space(responseIndex), solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvGrad)
{
  const RCP<NOXSolver<double> > solver = solverNew(epetraModelNew());

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(responseIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvJac_NoDgDxMv)
{
  // Disable support for MultiVector-based DgDx derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new Test::WeakenedModelEvaluator_NoDgDxMv(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::MultiVectorBase<double> > dgdp =
    Thyra::createMembers(solver->get_g_space(responseIndex), solver->get_p_space(parameterIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_JACOBIAN_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  TEST_EQUALITY(dgdp->domain()->dim(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const Array<double> actual = arrayFromVector(*dgdp->col(i));
    TEST_COMPARE_FLOATING_ARRAYS(actual, arrayView(&expected[i], 1), tol);
  }
}

TEUCHOS_UNIT_TEST(Piro_NOXSolver, SensitivityMvGrad_NoDgDpMvJac)
{
  // Disable support for MultiVector-based DgDx derivative
  // (Only LinOp form is available)
  const RCP<Thyra::ModelEvaluatorDefaultBase<double> > weakenedModel =
      rcp(new Test::WeakenedModelEvaluator_NoDgDpMvJac(thyraModelNew(epetraModelNew())));
  const RCP<NOXSolver<double> > solver = solverNew(weakenedModel);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();

  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  const int responseIndex = 0;
  const int parameterIndex = 0;
  const RCP<Thyra::VectorBase<double> > dgdp =
    Thyra::createMember(solver->get_p_space(responseIndex));
  const Thyra::MEB::Derivative<double> dgdp_deriv(dgdp, Thyra::MEB::DERIV_MV_GRADIENT_FORM);
  outArgs.set_DgDp(responseIndex, parameterIndex, dgdp_deriv);

  solver->evalModel(inArgs, outArgs);

  const Array<double> expected = tuple(2.0, -8.0);
  const Array<double> actual = arrayFromVector(*dgdp);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

#endif /*Piro_ENABLE_NOX*/
