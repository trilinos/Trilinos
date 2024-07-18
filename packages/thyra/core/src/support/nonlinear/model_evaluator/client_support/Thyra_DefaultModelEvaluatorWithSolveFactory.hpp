// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MODEL_EVALUATOR_WITH_SOLVE_FACTORY_HPP
#define THYRA_DEFAULT_MODEL_EVALUATOR_WITH_SOLVE_FACTORY_HPP


#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Teuchos_Time.hpp"


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and uses a compatible
 * LinearOpWithSolveFactory object to create a LinearOpWithSolveBase version
 * of W.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultModelEvaluatorWithSolveFactory
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultModelEvaluatorWithSolveFactory();

  /** \brief . */
  DefaultModelEvaluatorWithSolveFactory(
    const RCP<ModelEvaluator<Scalar> > &thyraModel,
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &W_factory
    );

  /** \brief . */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel,
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &W_factory
    );

  /** \brief . */
  void uninitialize(
    RCP<ModelEvaluator<Scalar> > *thyraModel = NULL,
    RCP<LinearOpWithSolveFactoryBase<Scalar> > *W_factory = NULL
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> > create_W() const;

  //@}

  /** \name Public functions overridden from ModelEvaluatorDelegatorBase. */
  //@{

  /** \brief . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  RCP<LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

};


// /////////////////////////////////
// Implementations


// Constructors/initializers/accessors/utilities


template<class Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::DefaultModelEvaluatorWithSolveFactory()
{}


template<class Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::DefaultModelEvaluatorWithSolveFactory(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &W_factory
  )
{
  initialize(thyraModel,W_factory);
}


template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &W_factory
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  W_factory_ = W_factory;
}


template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::uninitialize(
  RCP<ModelEvaluator<Scalar> > *thyraModel,
  RCP<LinearOpWithSolveFactoryBase<Scalar> > *W_factory
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  if(W_factory) *W_factory = W_factory_;
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
  W_factory_ = Teuchos::null;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultModelEvaluatorWithSolveFactory<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultModelEvaluatorWithSolveFactory{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << ",W_factory=";
  if(W_factory_.get())
    oss << "\'"<<W_factory_->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from ModelEvaluator.


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_W() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    W_factory_.get()==NULL, std::logic_error
    ,"Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>::create_W(): "
    "Error, the client did not set a LinearOpWithSolveFactoryBase object for W!"
    );
  W_factory_->setOStream(this->getOStream());
  W_factory_->setVerbLevel(this->getVerbLevel());
  return W_factory_->createOp();
}


// Overridden from ModelEvaluatorDelegatorBase.


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_W_factory() const
{
  return W_factory_;
}


// Private functions overridden from ModelEvaluatorDefaultBase.


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::createOutArgsImpl() const
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(wrappedOutArgs.Np(),wrappedOutArgs.Ng());
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setSupports(MEB::OUT_ARG_W,
    wrappedOutArgs.supports(MEB::OUT_ARG_W_op)&&W_factory_.get()!=NULL);
  return outArgs;
}


template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultModelEvaluatorWithSolveFactory",inArgs,outArgs
    );

  Teuchos::Time timer("");

  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveFactoryBase<Scalar> >
    VOTSLOWSF;
  VOTSLOWSF W_factory_outputTempState(W_factory_,out,verbLevel);

  // InArgs

  MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();

  wrappedInArgs.setArgs(inArgs,true);

  // OutArgs

  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();

  wrappedOutArgs.setArgs(outArgs,true);

  RCP<LinearOpWithSolveBase<Scalar> > W;
  RCP<const LinearOpBase<Scalar> > fwdW;
  if( outArgs.supports(MEB::OUT_ARG_W) && (W = outArgs.get_W()).get() ) {
    Thyra::uninitializeOp<Scalar>(*W_factory_, W.ptr(), outArg(fwdW));

    {
      // Handle this case later if we need to!
      const bool both_W_and_W_op_requested = nonnull(outArgs.get_W_op());
      TEUCHOS_TEST_FOR_EXCEPT(both_W_and_W_op_requested);
    }

    RCP<LinearOpBase<Scalar> > nonconst_fwdW;
    if(fwdW.get()) {
      nonconst_fwdW = rcp_const_cast<LinearOpBase<Scalar> >(fwdW);
    }
    else {
      nonconst_fwdW = thyraModel->create_W_op();
      fwdW = nonconst_fwdW;
    }

    wrappedOutArgs.set_W_op(nonconst_fwdW);
  }

  // Do the evaluation

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nEvaluating the output functions on model \'"
         << thyraModel->description() << "\' ...\n";
  timer.start(true);

  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);

  timer.stop();
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to evaluate underlying model = "
                   << timer.totalElapsedTime()<<" sec\n";

  // Postprocess arguments

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nPost processing the output objects ...\n";
  timer.start(true);

  if( W.get() ) {
    Thyra::initializeOp<Scalar>(*W_factory_, fwdW, W.ptr());
    W->setVerbLevel(this->getVerbLevel());
    W->setOStream(this->getOStream());
  }

  timer.stop();
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to process output objects = "
                   << timer.totalElapsedTime()<<" sec\n";

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Thyra


#endif // THYRA_DEFAULT_MODEL_EVALUATOR_WITH_SOLVE_FACTORY_HPP
