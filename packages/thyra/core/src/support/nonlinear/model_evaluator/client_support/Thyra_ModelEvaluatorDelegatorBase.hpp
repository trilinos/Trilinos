// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP
#define THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP


#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


namespace Thyra {


/** \brief This is a base class that delegetes almost all function to a
 * wrapped model evaluator object.
 *
 * This class makes it easy to write many different types of decorator
 * subclasses by only requiring those subclasses to override just the behavior
 * that they need to overide and leave the rest alone.  Note that whenever the
 * definition of <tt>ModelEvaluator</tt> changes, this class will have to be
 * updated.  However, the decorator subclasses that derive from this base
 * class might be able to remain unchanged.  Note that changes in only the
 * <tt>ModelEvaluatorBase::InArgs</tt> and
 * <tt>ModelEvaluatorBase::OutArgs</tt> classes should not require any changes
 * here.
 * 
 * The only functions that a client must override in order to create a
 * concrete subcalss is the <tt>evalModel()</tt> function.  All other
 * functions have implementations here that simply delegate to the model
 * evaluator object returned from <tt>getUnderlyingModel()</tt>.  However,
 * most decorator classes will need to override at least one other function.
 *
 * This class provides a default implemntation of the machinary to store and
 * access the wrapped model evaluator object.  A subclass can choose to ignore
 * this and override the functions <tt>isUnderlyingModelConst()<tt>,
 * <tt>getConstUnderlyingModel()</tt>, and <tt>getUnderlyingModel()</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class ModelEvaluatorDelegatorBase
  : virtual public ModelEvaluatorDefaultBase<Scalar>
{
public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Constructs to uninitialized. */
  ModelEvaluatorDelegatorBase();

  /** \brief Calls <tt>initialize()</tt>. */
  ModelEvaluatorDelegatorBase(
    const RCP<ModelEvaluator<Scalar> > &model
    );

  /** \brief Calls <tt>initialize()</tt>. */
  ModelEvaluatorDelegatorBase(
    const RCP<const ModelEvaluator<Scalar> > &model
    );

  /** \brief Initialize given a non-const model evaluator. */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &model
    );

  /** \brief Initialize given a const model evaluator. */
  void initialize(
    const RCP<const ModelEvaluator<Scalar> > &model
    );

  /** \brief Uninitialize. */
  void uninitialize();

  //@}

  /** \name Virtual functions that can overriden */
  //@{

  /** \brief . */
  virtual bool isUnderlyingModelConst() const;

  /** \brief . */
  virtual RCP<ModelEvaluator<Scalar> > getNonconstUnderlyingModel();

  /** \brief . */
  virtual RCP<const ModelEvaluator<Scalar> > getUnderlyingModel() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

  //@}

protected:

  /** \name Producted utility functions to be called by subclasses */
  //@{

  /** \brief Set a valid parameter for reading the local verbosity level. */
  void setLocalVerbosityLevelValidatedParameter(
    ParameterList *paramList
    ) const;

  /** \brief Read the local verbosity level parameter. */
  Teuchos::EVerbosityLevel readLocalVerbosityLevelValidatedParameter(
    ParameterList &paramList 
    ) const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase */
  //@{

  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DfDp_op_impl(int l) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDx_dot_op_impl(int j) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDx_op_impl(int j) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDp_op_impl( int j, int l ) const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  //@}

private: // Data members

  Teuchos::ConstNonconstObjectContainer<ModelEvaluator<Scalar> > model_;

  static
  RCP<
    Teuchos::StringToIntegralParameterEntryValidator<
      Teuchos::EVerbosityLevel
      >
    > LocalVerbosityLevel_validator_;
  static const std::string LocalVerbosityLevel_name_;
  static const Teuchos::EVerbosityLevel LocalVerbosityLevel_enum_default_;
  static const std::string LocalVerbosityLevel_default_;
  
};


#define THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_LOCALVERBLEVEL_BEGIN(CLASS_NAME,INARGS,OUTARGS,UNDERLYINGMODEL,LOCALVERBLEVEL) \
  \
  using Teuchos::includesVerbLevel; \
  using Teuchos::RCP; \
  using Teuchos::EVerbosityLevel; \
  const std::string blahblah_classNameStr \
    = std::string(CLASS_NAME)+"<"+Teuchos::ScalarTraits<Scalar>::name()+">"; \
  const std::string blahblah_classFuncNameStr \
    = blahblah_classNameStr+"::evalModel(...)"; \
  TEUCHOS_FUNC_TIME_MONITOR(blahblah_classFuncNameStr); \
  \
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &blahblah_outArgs = (OUTARGS); \
  \
  Teuchos::Time totalTimer(""); \
  totalTimer.start(true); \
  \
  const RCP<Teuchos::FancyOStream> out = this->getOStream(); \
  const EVerbosityLevel verbLevel = this->getVerbLevel(); \
  const EVerbosityLevel localVerbLevelInput = (LOCALVERBLEVEL); \
  const EVerbosityLevel localVerbLevel = \
    ( localVerbLevelInput==Teuchos::VERB_DEFAULT ? verbLevel : localVerbLevelInput ); \
  Teuchos::OSTab tab(out); \
  if(out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_LOW)) \
    *out << "\nEntering " << blahblah_classFuncNameStr << " ...\n"; \
  \
  if(out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_MEDIUM)) \
    *out \
      << "\ninArgs =\n" << Teuchos::describe((INARGS),localVerbLevel) \
      << "\noutArgs on input =\n" << Teuchos::describe((OUTARGS),Teuchos::VERB_LOW); \
  \
  const RCP<const Thyra::ModelEvaluator<Scalar> > \
    thyraModel = (UNDERLYINGMODEL); \
  \
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME; \
  VOTSME thyraModel_outputTempState(thyraModel,out,verbLevel)


#define THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_LOCALVERBLEVEL_BEGIN(CLASS_NAME,INARGS,OUTARGS,LOCALVERBLEVEL) \
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_LOCALVERBLEVEL_BEGIN(CLASS_NAME,INARGS,OUTARGS,this->getUnderlyingModel(),LOCALVERBLEVEL)


#define THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(CLASS_NAME,INARGS,OUTARGS,UNDERLYINGMODEL) \
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_LOCALVERBLEVEL_BEGIN(CLASS_NAME,INARGS,OUTARGS,UNDERLYINGMODEL,Teuchos::VERB_DEFAULT)


#define THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(CLASS_NAME,INARGS,OUTARGS) \
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(CLASS_NAME,INARGS,OUTARGS,this->getUnderlyingModel())


#define THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END() \
  if(out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_MEDIUM)) \
    *out \
      << "\noutArgs on output =\n" << Teuchos::describe(blahblah_outArgs,localVerbLevel); \
  \
  totalTimer.stop(); \
  if(out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_LOW)) \
    *out \
      << "\nTotal evaluation time = "<<totalTimer.totalElapsedTime()<<" sec\n" \
      << "\nLeaving " << blahblah_classFuncNameStr << " ...\n"


// /////////////////////////////////
// Implementations


// Static class data members


template<class Scalar>
RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
  Teuchos::EVerbosityLevel
  >
>
ModelEvaluatorDelegatorBase<Scalar>::LocalVerbosityLevel_validator_;

template<class Scalar>
const std::string
ModelEvaluatorDelegatorBase<Scalar>::LocalVerbosityLevel_name_
= "Local Verbosity Level";

template<class Scalar>
const Teuchos::EVerbosityLevel
ModelEvaluatorDelegatorBase<Scalar>::LocalVerbosityLevel_enum_default_
= Teuchos::VERB_DEFAULT;

template<class Scalar>
const std::string
ModelEvaluatorDelegatorBase<Scalar>::LocalVerbosityLevel_default_
= getVerbosityLevelParameterValueName(
  ModelEvaluatorDelegatorBase<Scalar>::LocalVerbosityLevel_enum_default_
  );


// Constructors/initializers


template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase()
{}


template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase(
  const RCP<ModelEvaluator<Scalar> >   &model
  )
{
  this->initialize(model);
}


template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase(
  const RCP<const ModelEvaluator<Scalar> >   &model
  )
{
  this->initialize(model);
}


template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> >   &model
  )
{
  model_.initialize(model);
}


template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::initialize(
  const RCP<const ModelEvaluator<Scalar> >   &model
  )
{
  model_.initialize(model);
}


template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::uninitialize()
{
  model_.uninitialize();
}


// Virtual functions that can overriden


template<class Scalar>
bool ModelEvaluatorDelegatorBase<Scalar>::isUnderlyingModelConst() const
{
  return model_.isConst();
}


template<class Scalar>
RCP<ModelEvaluator<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::getNonconstUnderlyingModel()
{
  return model_.getNonconstObj();
}


template<class Scalar>
RCP<const ModelEvaluator<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::getUnderlyingModel() const
{
  return model_.getConstObj();
}


// Overridden from ModelEvaulator.


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_x_space() const
{
  return getUnderlyingModel()->get_x_space();
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_f_space() const
{
  return getUnderlyingModel()->get_f_space();
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_p_space(int l) const
{
  return getUnderlyingModel()->get_p_space(l);
}


template<class Scalar>
RCP<const Teuchos::Array<std::string> >
ModelEvaluatorDelegatorBase<Scalar>::get_p_names(int l) const
{
  return getUnderlyingModel()->get_p_names(l);
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_g_space(int j) const
{
  return getUnderlyingModel()->get_g_space(j);
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getNominalValues() const
{
  return getUnderlyingModel()->getNominalValues();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getLowerBounds() const
{
  return getUnderlyingModel()->getLowerBounds();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getUpperBounds() const
{
  return getUnderlyingModel()->getUpperBounds();
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_W() const
{
  return getUnderlyingModel()->create_W();
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_W_op() const
{
  return getUnderlyingModel()->create_W_op();
}


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_W_factory() const
{
  return getUnderlyingModel()->get_W_factory();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::createInArgs() const
{
  ModelEvaluatorBase::InArgsSetup<Scalar> inArgs = getUnderlyingModel()->createInArgs();
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}


template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  getNonconstUnderlyingModel()->reportFinalPoint(finalPoint,wasSolved);
}


// protected


// Producted utility functions to be called by subclasses


template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::setLocalVerbosityLevelValidatedParameter(
    ParameterList *paramList
    ) const
{
  TEST_FOR_EXCEPT(0==paramList);
  if (is_null(LocalVerbosityLevel_validator_))
    LocalVerbosityLevel_validator_ =
      Teuchos::verbosityLevelParameterEntryValidator(
        LocalVerbosityLevel_name_
        );
  paramList->set(
    LocalVerbosityLevel_name_, LocalVerbosityLevel_default_,
    "Overriding verbosity level for this model evaluator object.\n"
    "This level will not propagate to nested model evaluator objects\n"
    "The value of \"default\" result in the object verbosity level being\n"
    "used instead.",
    LocalVerbosityLevel_validator_
    );
}


template<class Scalar>
Teuchos::EVerbosityLevel
ModelEvaluatorDelegatorBase<Scalar>::readLocalVerbosityLevelValidatedParameter(
  ParameterList &paramList 
  ) const
{
  return LocalVerbosityLevel_validator_->getIntegralValue(
    paramList, LocalVerbosityLevel_name_, LocalVerbosityLevel_default_ );
}


// private


// Producted functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DfDp_op_impl(int l) const
{
  return getUnderlyingModel()->create_DfDp_op(l);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DgDx_dot_op_impl(
  int j
  ) const
{
  return getUnderlyingModel()->create_DgDx_dot_op(j);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DgDx_op_impl(
  int j
  ) const
{
  return getUnderlyingModel()->create_DgDx_op(j);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DgDp_op_impl(
  int j, int l
  ) const
{
  return getUnderlyingModel()->create_DgDp_op(j,l);
}


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::createOutArgsImpl() const
{
  ModelEvaluatorBase::OutArgsSetup<Scalar>
    outArgs = getUnderlyingModel()->createOutArgs();
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}


} // namespace Thyra

#endif // THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP
