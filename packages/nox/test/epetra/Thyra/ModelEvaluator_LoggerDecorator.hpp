// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_MODEL_EVALUATOR_LOGGER_DECORATOR_HPP
#define NOX_MODEL_EVALUATOR_LOGGER_DECORATOR_HPP

namespace NOX_TEST {


  /** \brief Model Evaluator Decorator unit testing class that logs
      calls to evalModel() to count the number of calls made to
      evaluate the residual and if each call is exact or approximate
      (used in a Jacobian approximation in JFNK perturbation).  This
      is used to unit test the JFNK support in NOX.
  */
  template<typename ScalarT>
  class ModelEvaluatorLoggerDecorator : public ::Thyra::ModelEvaluatorDelegatorBase<ScalarT> {

  public:

    ModelEvaluatorLoggerDecorator(const Teuchos::RCP< ::Thyra::ModelEvaluator<ScalarT> > &model) :
      ::Thyra::ModelEvaluatorDelegatorBase<ScalarT>(model),
      num_calls_exact_(0),
      num_calls_approx_deriv_(0),
      num_calls_very_approx_deriv_(0)
    { }

    void evalModelImpl(const ::Thyra::ModelEvaluatorBase::InArgs<ScalarT> &inArgs,
               const ::Thyra::ModelEvaluatorBase::OutArgs<ScalarT> &outArgs) const
    {
      ::Thyra::ModelEvaluatorBase::InArgs<ScalarT> modelInArgs = this->getUnderlyingModel()->createInArgs();
      modelInArgs.setArgs(inArgs,false);

      ::Thyra::ModelEvaluatorBase::OutArgs<ScalarT> modelOutArgs = this->getUnderlyingModel()->createOutArgs();
      modelOutArgs.setArgs(outArgs,false);

      ::Thyra::ModelEvaluatorBase::Evaluation< ::Thyra::VectorBase<ScalarT> > eval = outArgs.get_f();
      if (nonnull(eval)) {
    if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_EXACT) {
      ++num_calls_exact_;
      //std::cout << "EXACT!" << std::endl;
    }
    else if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_APPROX_DERIV) {
      ++num_calls_approx_deriv_;
      //std::cout << "APPROX DERIV!" << std::endl;
    }
    else if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_VERY_APPROX_DERIV) {
      ++num_calls_very_approx_deriv_;
      //std::cout << "VERY APPROX DERIV!" << std::endl;
    }
      }

      this->getUnderlyingModel()->evalModel(modelInArgs,modelOutArgs);
    }

    int getNumExactCalls() const
    { return num_calls_exact_; }

    int getNumApproxDerivCalls() const
    { return num_calls_approx_deriv_; }

    int getNumVeryApproxDerivCalls() const
    { return num_calls_very_approx_deriv_; }

  private:

    void reset()
    {
      num_calls_exact_ = 0;
      num_calls_approx_deriv_ = 0;
      num_calls_very_approx_deriv_ = 0;
    }

    mutable int num_calls_exact_;
    mutable int num_calls_approx_deriv_;
    mutable int num_calls_very_approx_deriv_;

  };

}

#endif
