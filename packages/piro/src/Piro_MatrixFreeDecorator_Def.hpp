// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_MatrixFreeDecorator.hpp"

#include "Piro_MatrixFreeLinearOp.hpp"

#include "Thyra_VectorStdOps.hpp"

template <typename Scalar>
Piro::MatrixFreeDecorator<Scalar>::MatrixFreeDecorator(
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model, 
  double lambda_) :
  Thyra::ModelEvaluatorDelegatorBase<Scalar>(model), 
  lambda(lambda_) 
{
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::MatrixFreeDecorator<Scalar>::create_W_op() const
{
  return Teuchos::rcp(new MatrixFreeLinearOp<Scalar>(this->getUnderlyingModel(), lambda));
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::MatrixFreeDecorator<Scalar>::createOutArgsImpl() const
{
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs =
    this->getUnderlyingModel()->createOutArgs();

  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> result = modelOutArgs;
  result.setModelEvalDescription(this->description());

  const bool modelSupportsResidual = modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, modelSupportsResidual);

  return result;
}

template <typename Scalar>
void
Piro::MatrixFreeDecorator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;

  // Forward all input arguments
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs =
    this->getUnderlyingModel()->createInArgs();
  modelInArgs.setArgs(inArgs, /*ignoreUnsupported =*/ false);

  // Forward all supported output arguments, except Jacobian operator
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs =
    this->getUnderlyingModel()->createOutArgs();
  modelOutArgs.setArgs(outArgs, /*ignoreUnsupported =*/ true);
  if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op)) {
    modelOutArgs.set_W_op(Teuchos::null);
  }

  const bool computeResidual = Teuchos::nonnull(outArgs.get_f());
  const bool computeJacobian = Teuchos::nonnull(outArgs.get_W_op());

  if (computeJacobian && !computeResidual) {
    // Residual required to compute Jacobian
    // Must be computed even if not requested by user
    modelOutArgs.set_f(Thyra::createMember(this->getUnderlyingModel()->get_f_space()));
  }

  this->getUnderlyingModel()->evalModel(inArgs, modelOutArgs);

  if (computeJacobian) {
    const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
    const RCP<MatrixFreeLinearOp<Scalar> > W_actual =
      Teuchos::rcp_dynamic_cast<MatrixFreeLinearOp<Scalar> >(W_out);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(W_actual));

    RCP<const Thyra::VectorBase<Scalar> > f_base;
    {
      const RCP<const Thyra::VectorBase<Scalar> > f_out = modelOutArgs.get_f();
      if (computeResidual) {
        // Make a deep copy to assert exclusive ownership of *f_out
        // to avoid aliasing or dangling pointer
        f_base = f_out->clone_v();
      } else {
        // A shallow copy is safe since we created and therefore are the exclusive owner of *f_out
        f_base = f_out;
      }
    }

    // Deep copy of the input parameters since we are not guaranteed to own the referenced objects
    Thyra::ModelEvaluatorBase::InArgs<Scalar> clonedModelInArgs =
      this->getUnderlyingModel()->createInArgs();
    clonedModelInArgs.setArgs(modelInArgs, /*ignoreUnsupported =*/ false, /*cloneObjects =*/ true);

    // Pass references to operator for shallow copies
    W_actual->setBase(clonedModelInArgs, f_base);
  }
}
