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

RCP<Epetra_Vector> vectorFromLinOp(const Epetra_Operator &op, int col)
{
  const RCP<Epetra_Vector> result = vectorNew(op.OperatorRangeMap());

  const RCP<Epetra_Vector> rhs = vectorNew(op.OperatorDomainMap());
  const double value = 1.0;
  const int ierr = rhs->ReplaceGlobalValues(1, &value, &col);
  // Some processes might not hold the entry
  TEUCHOS_ASSERT(ierr == 0 || ierr == 1);
  {
    const Epetra_Comm &comm = rhs->Comm();
    int ierrLoc = ierr, ierrSum;
    comm.SumAll(&ierrLoc, &ierrSum, 1);
    // At least one process holds the entry
    TEUCHOS_ASSERT(ierrSum < comm.NumProc());
  }

  op.Apply(*rhs, *result);
  return result;
}

Array<double> arrayFromLinOp(const Epetra_Operator &op, int col)
{
  const RCP<const Epetra_Vector> v = vectorFromLinOp(op, col);
  return arrayFromVector(*v);
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
  jacobian->setBase(modelInArgs, f_base, /*haveXdot =*/ false);

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
  jacobian->setBase(modelInArgs, f_base, /*haveXdot =*/ true);

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
