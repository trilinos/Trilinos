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

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and uses a compatible
 * LinearOpWithSolveFactory object to create a LinearOpWithSolveBase version
 * of W.
 *
 * ToDo: Finish documentation!
 *
 */
template<class Scalar>
class DefaultModelEvaluatorWithSolveFactory : virtual public ModelEvaluator<Scalar> {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

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
  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > getUnderlyingModel() const;

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel = NULL
    ,Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  *W_factory   = NULL
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
	int Np() const;
  /** \brief . */
	int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
	Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
	Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_init() const;
  /** \brief . */
	Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_init(int l) const;
  /** \brief . */
  ScalarMag get_t_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_lower_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_upper_bounds() const;
  /** \brief . */
	Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_lower_bounds(int l) const;
  /** \brief . */
	Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_upper_bounds(int l) const;
  /** \brief . */
  ScalarMag get_t_lower_bound() const;
  /** \brief . */
  ScalarMag get_t_upper_bound() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DfDp_op(int l) const;
  /** \brief . */
  ModelEvaluatorBase::DerivativeMultiVector<Scalar> create_DfDp_mv(int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  /** \brief . */
  ModelEvaluatorBase::DerivativeMultiVector<Scalar> create_DgDx_mv(int j, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const;
  /** \brief . */
  ModelEvaluatorBase::DerivativeMultiVector<Scalar> create_DgDp_mv( int j, int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation ) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                      thyraModel_;
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >        W_factory_;
  
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
  thyraModel_ = thyraModel;
  W_factory_ = W_factory;
}

template<class Scalar>
Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::getUnderlyingModel() const
{
  return thyraModel_;
}

template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel
  ,Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<Scalar> >  *W_factory
  )
{
  if(thyraModel) *thyraModel = thyraModel_;
  if(W_factory) *W_factory = W_factory_;
  thyraModel_ = Teuchos::null;
  W_factory_ = Teuchos::null;
}

// Overridden from ModelEvaulator.

template<class Scalar>
int DefaultModelEvaluatorWithSolveFactory<Scalar>::Np() const
{
  return thyraModel_->Np();
}

template<class Scalar>
int DefaultModelEvaluatorWithSolveFactory<Scalar>::Ng() const
{
  return thyraModel_->Ng();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_x_space() const
{
  return thyraModel_->get_x_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_f_space() const
{
  return thyraModel_->get_f_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_p_space(int l) const
{
  return thyraModel_->get_p_space(l);
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_g_space(int j) const
{
  return thyraModel_->get_g_space(j);
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_x_init() const
{
  return thyraModel_->get_x_init();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_p_init(int l) const
{
  return thyraModel_->get_p_init(l);
}

template<class Scalar>
typename DefaultModelEvaluatorWithSolveFactory<Scalar>::ScalarMag
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_t_init() const
{
  return thyraModel_->get_t_init();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_x_lower_bounds() const
{
  return thyraModel_->get_x_lower_bounds();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_x_upper_bounds() const
{
  return thyraModel_->get_x_upper_bounds();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_p_lower_bounds(int l) const
{
  return thyraModel_->get_p_lower_bounds(l);
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_p_upper_bounds(int l) const
{
  return thyraModel_->get_p_upper_bounds(l);
}

template<class Scalar>
typename DefaultModelEvaluatorWithSolveFactory<Scalar>::ScalarMag
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_t_lower_bound() const
{
  return thyraModel_->get_t_lower_bound();
}

template<class Scalar>
typename DefaultModelEvaluatorWithSolveFactory<Scalar>::ScalarMag
DefaultModelEvaluatorWithSolveFactory<Scalar>::get_t_upper_bound() const
{
  return thyraModel_->get_t_upper_bound();
}

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
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_W_op() const
{
  return thyraModel_->create_W_op();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DfDp_op(int l) const
{
  return thyraModel_->create_DfDp_op(l);
}

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DfDp_mv(int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation) const
{
  return thyraModel_->create_DfDp_mv(l,orientation);
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DgDx_op(int j) const
{
  return thyraModel_->create_DgDx_op(j);
}

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DgDx_mv(int j, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation) const
{
  return thyraModel_->create_DgDx_mv(j,orientation);
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DgDp_op( int j, int l ) const
{
  return thyraModel_->create_DgDp_op(j,l);
}

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::create_DgDp_mv( int j, int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation ) const
{
  return thyraModel_->create_DgDp_mv(j,l,orientation);
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::createInArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const MEB::InArgs<Scalar> wrappedInArgs = thyraModel_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(wrappedInArgs.Np());
  inArgs.setSupports(wrappedInArgs);
  return inArgs;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultModelEvaluatorWithSolveFactory<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel_->createOutArgs();
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

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>::evalModel(...) ...\n";

  typedef Teuchos::VerboseObjectTempState<ModelEvaluatorBase> VOTSME;
  VOTSME thyraModel_outputTempState(thyraModel_,out,verbLevel);

  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveFactoryBase<Scalar> > VOTSLOWSF;
  VOTSLOWSF W_factory_outputTempState(W_factory_,out,verbLevel);
  
  // InArgs

  MEB::InArgs<Scalar> wrappedInArgs = thyraModel_->createInArgs();

  wrappedInArgs.setArgs(inArgs,true);

  // OutArgs

  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel_->createOutArgs();

  wrappedOutArgs.setArgs(outArgs,true);
  
  RefCountPtr<LinearOpWithSolveBase<Scalar> > W;
  RefCountPtr<LinearOpBase<Scalar> >          W_op;
  RefCountPtr<const LinearOpBase<Scalar> >    fwdW;
  RefCountPtr<LinearOpBase<Scalar> >          nonconst_fwdW;
  if( outArgs.supports(MEB::OUT_ARG_W) && (W = outArgs.get_W()).get() ) {
    W_factory_->uninitializeOp(&*W,&fwdW);
    if(fwdW.get()) {
      nonconst_fwdW = rcp_const_cast<LinearOpBase<Scalar> >(fwdW);
    }
    else {
      nonconst_fwdW = thyraModel_->create_W_op();
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
    *out << "\nEvaluating the output functions on model \'" << thyraModel_->description() << "\'  ...\n";
  timer.start(true);
  
  thyraModel_->evalModel(wrappedInArgs,wrappedOutArgs);
  
  timer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *OSTab(out).getOStream() << "\nTime to evaluate underlying model = "<<timer.totalElapsedTime()<<" sec\n";

  // Postprocess arguments
  
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nPost processing the output objects ...\n";
  timer.start(true);
  
  if( W.get() ) {
    W_factory_->initializeOp(fwdW,&*W);
    W->setVerbLevel(this->getVerbLevel());
    W->setOStream(this->getOStream());
  }
  
  if( W_op.get() ) {
    TEST_FOR_EXCEPT(true); // Handle this case later if we need to!
  }

  timer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *OSTab(out).getOStream() << "\nTime to process output objects = "<<timer.totalElapsedTime()<<" sec\n";


  totalTimer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal evaluation time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>::evalModel(...) ...\n";
  
}

template<class Scalar>
void DefaultModelEvaluatorWithSolveFactory<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  thyraModel_->reportFinalPoint(finalPoint,wasSolved);
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultModelEvaluatorWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "Thyra::DefaultModelEvaluatorWithSolveFactory{";
  oss << "thyraModel=";
  if(thyraModel_.get())
    oss << "\'"<<thyraModel_->description()<<"\'";
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
