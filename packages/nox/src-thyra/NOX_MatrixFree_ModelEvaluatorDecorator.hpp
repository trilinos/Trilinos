#ifndef NOX_MATRIX_FREE_MODEL_EVALUATOR_DECORATOR_HPP
#define NOX_MATRIX_FREE_MODEL_EVALUATOR_DECORATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

namespace NOX {
  
  /** \brief Model Evaluator Decorator class that adds support for the evaluation of a matrix-free W_op.

     The NOX solvers evaluate W_op as part of a Newton-based solve.
     In terms of a Jacobian-free Newton-Krylov method, we need to wrap
     the incoming model evaluator to say that it supports a ficticous
     W_op (otherwise all solvers would require built-in knowledge of
     the Matrix-Free algorithm).  This wrapper adds support for that.
  */
  template<typename ScalarT>
  class MatrixFreeModelEvaluatorDecorator : public ::Thyra::ModelEvaluatorDelegatorBase<ScalarT> {
    
  public:
    
    MatrixFreeModelEvaluatorDecorator(const Teuchos::RCP< ::Thyra::ModelEvaluator<ScalarT> > &model) :
      ::Thyra::ModelEvaluatorDelegatorBase<ScalarT>(model)
    { }

    ::Thyra::ModelEvaluatorBase::OutArgs<ScalarT> createOutArgsImpl() const
    {
      ::Thyra::ModelEvaluatorBase::OutArgsSetup<ScalarT> outArgs(this->getUnderlyingModel()->createOutArgs());
      if (outArgs.supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_op))
	  return outArgs;

      outArgs.setSupports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_op,true);
      return outArgs;
    }

    void evalModelImpl(const ::Thyra::ModelEvaluatorBase::InArgs<ScalarT> &inArgs,
		       const ::Thyra::ModelEvaluatorBase::OutArgs<ScalarT> &outArgs) const
    {
      ::Thyra::ModelEvaluatorBase::OutArgs<ScalarT> outArgsWithoutW_op = this->getUnderlyingModel()->createOutArgs();
      
      outArgsWithoutW_op.setArgs(outArgs, true);
      
      if (outArgsWithoutW_op.supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_op))
	outArgsWithoutW_op.set_W_op(Teuchos::null);

      this->getUnderlyingModel()->evalModel(inArgs,outArgsWithoutW_op);
    }

  };
}

#endif
