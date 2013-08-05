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
#include "Piro_LOCASolver.hpp"

#include "MockModelEval_A.hpp"

#include "Piro_Test_ThyraSupport.hpp"
#include "Piro_Test_MockObserver.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

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

const RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyraModelNew(const RCP<EpetraExt::ModelEvaluator> &epetraModel)
{
  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory(new Thyra::AmesosLinearOpWithSolveFactory);
  return epetraModelEvaluator(epetraModel, lowsFactory);
}

const RCP<LOCASolver<double> > solverNew(
    const RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    const RCP<Piro::ObserverBase<double> > &observer = Teuchos::null)
{
  const RCP<ParameterList> piroParams(new ParameterList("Piro Parameters"));
  updateParametersFromXmlFile("input_Solve_LOCA_1.xml", piroParams.ptr());
  return observedLocaSolver<double>(piroParams, thyraModel, observer);
}

const RCP<LOCASolver<double> > solverNew(
    const RCP<EpetraExt::ModelEvaluator> &epetraModel,
    const RCP<Piro::ObserverBase<double> > &observer = Teuchos::null)
{
  return solverNew(thyraModelNew(epetraModel), observer);
}

// Floating point tolerance
const double tol = 1.0e-8;

// Tests

TEUCHOS_UNIT_TEST(Piro_LOCASolver, Spaces)
{
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew(), observer);

  const Thyra::MEB::InArgs<double> inArgs = solver->getNominalValues();
  Thyra::MEB::OutArgs<double> outArgs = solver->createOutArgs();
  solver->evalModel(inArgs, outArgs);

  const Array<double> actual = arrayFromVector(*observer->lastSolution());
  const Array<double> expected = tuple(1.0, 2.0, 3.0, 4.0);
  TEST_COMPARE_FLOATING_ARRAYS(actual, expected, tol);
}

TEUCHOS_UNIT_TEST(Piro_LOCASolver, SolutionForAlternateParameterValues)
{
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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
  const RCP<LOCASolver<double> > solver = solverNew(epetraModelNew());

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

#endif /*Piro_ENABLE_NOX*/
