// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP

#include "Thyra_DefaultFiniteDifferenceModelEvaluator_decl.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


// Constructors/initializers/accessors/utilities


template<class Scalar>
ScaledModelEvaluator<Scalar>::ScaledModelEvaluator()
{}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string ScaledModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::ScaledModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ScaledModelEvaluator<Scalar>::
set_f_scaling(const RCP<const Thyra::VectorBase<Scalar> >& f_scaling)
{
  f_scaling_ = f_scaling;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
void ScaledModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::ScaledModelEvaluator",inArgs,outArgs
    );

  thyraModel->evalModel(inArgs, outArgs);

  if (nonnull(f_scaling_)) {

    const RCP<VectorBase<Scalar> > f = outArgs.get_f();
    if (nonnull(f)) {
      ele_wise_scale(*f_scaling_, f.ptr());
    }
    
    const RCP<LinearOpBase<Scalar> > W_op = outArgs.get_W_op();
    if (nonnull(W_op)) {
      const RCP<ScaledLinearOpBase<Scalar> > W_scaled =
        rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(W_op, true);
      W_scaled->scaleLeft(*f_scaling_);
    }

  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
 
}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
