//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_ImplicitAdjointModelEvaluator_hpp
#define Thyra_ImplicitAdjointModelEvaluator_hpp

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAdjointLinearOpWithSolve.hpp"
#include "Thyra_AdjointLinearOpWithSolveFactory.hpp"
#include "Thyra_AdjointPreconditionerFactory.hpp"

namespace Thyra {

/** \brief An implementation of AdjointModelEvaluatorBase that creates an
 implicit adjoint from the supplied model evaluator. */
/**
 */
template <typename Scalar>
class ImplicitAdjointModelEvaluator
  : public ModelEvaluatorDelegatorBase<Scalar> {
 public:
  //! Constructor
  ImplicitAdjointModelEvaluator(const RCP<const ModelEvaluator<Scalar> >& model)
    : ModelEvaluatorDelegatorBase<Scalar>(model)
  {
  }

  //! Constructor
  ImplicitAdjointModelEvaluator(const RCP<ModelEvaluator<Scalar> >& model)
    : ModelEvaluatorDelegatorBase<Scalar>(model)
  {
  }

  //! Destructor
  virtual ~ImplicitAdjointModelEvaluator() = default;

  //! Create adjoint solver
  RCP<LinearOpWithSolveBase<Scalar> > create_W() const
  {
    return nonconstAdjointLows(this->getUnderlyingModel()->create_W());
  }

  //! Create adjoint op
  RCP<LinearOpBase<Scalar> > create_W_op() const
  {
    return nonconstAdjoint(this->getUnderlyingModel()->create_W_op());
  }

  //! Create adjoint preconditioner
  RCP<PreconditionerBase<Scalar> > create_W_prec() const
  {
    return nonconstAdjointPreconditioner(
        this->getUnderlyingModel()->create_W_prec());
  }

  //! Get adjoint solver factory
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const
  {
    return adjointLinearOpWithSolveFactory(
        this->getUnderlyingModel()->get_W_factory());
  }

 private:
  void evalModelImpl(const ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                     const ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::OutArgs<Scalar> model_outArgs =
        this->getUnderlyingModel()->createOutArgs();

    if (model_outArgs.supports(MEB::OUT_ARG_W) &&
        outArgs.get_W() != Teuchos::null) {
      RCP<DefaultAdjointLinearOpWithSolve<Scalar> > adjoint_op =
          Teuchos::rcp_dynamic_cast<DefaultAdjointLinearOpWithSolve<Scalar> >(
              outArgs.get_W(), true);
      model_outArgs.set_W(adjoint_op->getNonconstOp());
    }

    if (model_outArgs.supports(MEB::OUT_ARG_W_op) &&
        outArgs.get_W_op() != Teuchos::null) {
      RCP<DefaultScaledAdjointLinearOp<Scalar> > adjoint_op =
          Teuchos::rcp_dynamic_cast<DefaultScaledAdjointLinearOp<Scalar> >(
              outArgs.get_W_op(), true);
      model_outArgs.set_W_op(adjoint_op->getNonconstOp());
    }

    if (model_outArgs.supports(MEB::OUT_ARG_W_prec) &&
        outArgs.get_W_prec() != Teuchos::null) {
      RCP<AdjointPreconditioner<Scalar> > adjoint_op =
          Teuchos::rcp_dynamic_cast<AdjointPreconditioner<Scalar> >(
              outArgs.get_W_prec(), true);
      model_outArgs.set_W_prec(adjoint_op->getNonconstPreconditioner());
    }

    this->getUnderlyingModel()->evalModel(inArgs, model_outArgs);
  }
};

template <typename Scalar>
RCP<ImplicitAdjointModelEvaluator<Scalar> > implicitAdjointModelEvaluator(
    const RCP<const ModelEvaluator<Scalar> >& model)
{
  return Teuchos::rcp(new ImplicitAdjointModelEvaluator<Scalar>(model));
}

template <typename Scalar>
RCP<ImplicitAdjointModelEvaluator<Scalar> > implicitAdjointModelEvaluator(
    const RCP<ModelEvaluator<Scalar> >& model)
{
  return Teuchos::rcp(new ImplicitAdjointModelEvaluator<Scalar>(model));
}

}  // namespace Thyra

#endif
