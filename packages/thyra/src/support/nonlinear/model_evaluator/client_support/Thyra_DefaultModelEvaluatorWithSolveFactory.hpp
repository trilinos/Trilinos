// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
 */
template<class Scalar>
class DefaultModelEvaluatorWithSolveFactory : virtual public ModelEvaluatorDelegatorBase<Scalar> {
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultModelEvaluatorWithSolveFactory();

  /** \brief . */
  DefaultModelEvaluatorWithSolveFactory(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  &W_factory
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  &W_factory
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel = NULL
    ,Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  *W_factory   = NULL
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >   W_factory_;
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::DefaultModelEvaluatorWithSolveFactory()
{}

template<class Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::DefaultModelEvaluatorWithSolveFactory(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  &W_factory
  )
{
  initialize(thyraModel,W_factory);
}

template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  &W_factory
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  W_factory_ = W_factory;
}

template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel
  ,Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  *W_factory
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  if(W_factory) *W_factory = W_factory_;
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
  W_factory_ = Teuchos::null;
}

// Overridden from ModelEvaulator.

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_W() const
{
  TEST_FOR_EXCEPTION(
    W_factory_.get()==NULL, std::logic_error
    ,"Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>::create_W(): "
    "Error, the client did not set a LinearOpWithSolveFactoryBase object for W!"
    );
  W_factory_->setOStream(this->getOStream());
  W_factory_->setVerbLevel(this->getVerbLevel());
  return W_factory_->createOp();
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(wrappedOutArgs.Np(),wrappedOutArgs.Ng());
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setSupports(MEB::OUT_ARG_W,wrappedOutArgs.supports(MEB::OUT_ARG_W_op)&&W_factory_.get()!=NULL);
  return outArgs;
}

template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultModelEvaluatorWithSolveFactory",inArgs,outArgs
    );

  Teuchos::Time timer("");

  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveFactoryBase<Scalar> > VOTSLOWSF;
  VOTSLOWSF W_factory_outputTempState(W_factory_,out,verbLevel);
  
  // InArgs

  MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();

  wrappedInArgs.setArgs(inArgs,true);

  // OutArgs

  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();

  wrappedOutArgs.setArgs(outArgs,true);
  
  RefCountPtr<LinearOpWithSolveBase<Scalar> > W;
  RefCountPtr<LinearOpBase<Scalar> >          W_op;
  RefCountPtr<const LinearOpBase<Scalar> >    fwdW;
  RefCountPtr<LinearOpBase<Scalar> >          nonconst_fwdW;
  if( outArgs.supports(MEB::OUT_ARG_W) && (W = outArgs.get_W()).get() ) {
    Thyra::uninitializeOp<Scalar>(*W_factory_,&*W,&fwdW);
    if(fwdW.get()) {
      nonconst_fwdW = rcp_const_cast<LinearOpBase<Scalar> >(fwdW);
    }
    else {
      nonconst_fwdW = thyraModel->create_W_op();
      fwdW = nonconst_fwdW;
    }
  }
  if( outArgs.supports(MEB::OUT_ARG_W_op) && (W_op = outArgs.get_W_op()).get() ) {
    if( W_op.get() && !nonconst_fwdW.get() )
      nonconst_fwdW = rcp_const_cast<LinearOpBase<Scalar> >(fwdW);
  }
  if(nonconst_fwdW.get()) {
    wrappedOutArgs.set_W_op(nonconst_fwdW);
  }

  // Do the evaluation
  
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEvaluating the output functions on model \'" << thyraModel->description() << "\'  ...\n";
  timer.start(true);
  
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  
  timer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to evaluate underlying model = "<<timer.totalElapsedTime()<<" sec\n";

  // Postprocess arguments
  
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nPost processing the output objects ...\n";
  timer.start(true);
  
  if( W.get() ) {
    Thyra::initializeOp<Scalar>(*W_factory_,fwdW,&*W);
    W->setVerbLevel(this->getVerbLevel());
    W->setOStream(this->getOStream());
  }
  
  if( W_op.get() ) {
    TEST_FOR_EXCEPT(true); // Handle this case later if we need to!
  }

  timer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to process output objects = "<<timer.totalElapsedTime()<<" sec\n";

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultModelEvaluatorWithSolveFactory<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
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

} // namespace Thyra

#endif // THYRA_DEFAULT_MODEL_EVALUATOR_WITH_SOLVE_FACTORY_HPP
