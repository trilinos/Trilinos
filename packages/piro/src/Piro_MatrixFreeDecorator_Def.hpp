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

#include "Piro_MatrixFreeDecorator.hpp"

#include "Piro_MatrixFreeLinearOp.hpp"

#include "Thyra_VectorStdOps.hpp"

template <typename Scalar>
Piro::MatrixFreeDecorator<Scalar>::MatrixFreeDecorator(
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model) :
  Thyra::ModelEvaluatorDelegatorBase<Scalar>(model)
{
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::MatrixFreeDecorator<Scalar>::create_W_op() const
{
  return Teuchos::rcp(new MatrixFreeLinearOp<Scalar>(this->getUnderlyingModel()));
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
