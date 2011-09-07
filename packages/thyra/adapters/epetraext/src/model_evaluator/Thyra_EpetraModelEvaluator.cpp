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

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp" // Gives verbose macros!
#include "EpetraExt_ModelEvaluatorScalingTools.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_Time.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace {


const std::string StateFunctionScaling_name = "State Function Scaling";
Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::EpetraModelEvaluator::EStateFunctionScaling
    >
  >
stateFunctionScalingValidator;
const std::string StateFunctionScaling_default = "None";


// Extract out the Epetra_RowMatrix from the set W in an Epetra OutArgs object
Teuchos::RCP<Epetra_RowMatrix>
get_Epetra_RowMatrix(
  const EpetraExt::ModelEvaluator::OutArgs &epetraOutArgs
  )
{
  using Teuchos::RCP;
  const RCP<Epetra_Operator>
    eW = epetraOutArgs.get_W();
  const RCP<Epetra_RowMatrix>
    ermW = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(eW,false);
  TEST_FOR_EXCEPTION(
    is_null(ermW), std::logic_error,
    "Thyra::EpetraModelEvaluator::evalModel(...): Error, if\n"
    "scaling is turned on, the underlying Epetra_Operator created\n"
    "an initialized by the underlying epetra model evaluator\n"
    "\"" << epetraOutArgs.modelEvalDescription() << "\"\n"
    "must support the Epetra_RowMatrix interface through a dynamic cast.\n"
    "The concrete type " << Teuchos::typeName(*eW) << " does not support\n"
    "Epetra_RowMatrix!"
    );
  return ermW;
}


Teuchos::RCP<Epetra_Operator>
create_and_assert_W( 
  const EpetraExt::ModelEvaluator &epetraModel
  )
{
  using Teuchos::RCP;
  RCP<Epetra_Operator>
    eW = epetraModel.create_W();
  TEST_FOR_EXCEPTION(
    is_null(eW), std::logic_error,
    "Error, the call to create_W() returned null on the "
    "EpetraExt::ModelEvaluator object "
    "\"" << epetraModel.description() << "\".  This may mean that "
    "the underlying model does not support more than one copy of "
    "W at one time!" );
  return eW;
}


} // namespace


namespace Thyra {


// Constructors/initializers/accessors.


EpetraModelEvaluator::EpetraModelEvaluator()
  :nominalValuesAndBoundsAreUpdated_(false), stateFunctionScaling_(STATE_FUNC_SCALING_NONE),
   currentInArgsOutArgs_(false), finalPointWasSolved_(false)
{}


EpetraModelEvaluator::EpetraModelEvaluator(
  const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
  const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
  )
  :nominalValuesAndBoundsAreUpdated_(false), stateFunctionScaling_(STATE_FUNC_SCALING_NONE),
   currentInArgsOutArgs_(false), finalPointWasSolved_(false)
{
  initialize(epetraModel,W_factory);
}


void EpetraModelEvaluator::initialize(
  const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
  const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
  )
{
  using Teuchos::implicit_cast;
  typedef ModelEvaluatorBase MEB;
  //
  epetraModel_ = epetraModel;
  //
  W_factory_ = W_factory;
  //
  x_map_ = epetraModel_->get_x_map();
  f_map_ = epetraModel_->get_f_map();
  if (!is_null(x_map_)) {
    x_space_ = create_VectorSpace(x_map_);
    f_space_ = create_VectorSpace(f_map_);
  }
  //
  EpetraExt::ModelEvaluator::InArgs inArgs = epetraModel_->createInArgs();
  p_map_.resize(inArgs.Np()); p_space_.resize(inArgs.Np());
  p_map_is_local_.resize(inArgs.Np(),false);
  for( int l = 0; l < implicit_cast<int>(p_space_.size()); ++l ) {
    RCP<const Epetra_Map>
      p_map_l = ( p_map_[l] = epetraModel_->get_p_map(l) );
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      is_null(p_map_l), std::logic_error,
      "Error, the the map p["<<l<<"] for the model \""
      <<epetraModel->description()<<"\" can not be null!");
#endif

    p_map_is_local_[l] = !p_map_l->DistributedGlobal();
    p_space_[l] = create_VectorSpace(p_map_l);
  }
  //
  EpetraExt::ModelEvaluator::OutArgs outArgs = epetraModel_->createOutArgs();
  g_map_.resize(outArgs.Ng()); g_space_.resize(outArgs.Ng());
  g_map_is_local_.resize(outArgs.Ng(),false);
  for( int j = 0; j < implicit_cast<int>(g_space_.size()); ++j ) {
    RCP<const Epetra_Map>
      g_map_j = ( g_map_[j] = epetraModel_->get_g_map(j) );
    g_map_is_local_[j] = !g_map_j->DistributedGlobal();
    g_space_[j] = create_VectorSpace( g_map_j );
  }
  //
  epetraInArgsScaling_ = epetraModel_->createInArgs();
  epetraOutArgsScaling_ = epetraModel_->createOutArgs();
  nominalValuesAndBoundsAreUpdated_ = false;
  finalPointWasSolved_ = false;
  stateFunctionScalingVec_ = Teuchos::null; // Must set new scaling!
  //
  currentInArgsOutArgs_ = false;
}


RCP<const EpetraExt::ModelEvaluator>
EpetraModelEvaluator::getEpetraModel() const
{
  return epetraModel_;
}


void EpetraModelEvaluator::setNominalValues(
  const ModelEvaluatorBase::InArgs<double>& nominalValues
 )
{
  nominalValues_.setArgs(nominalValues);
  // Note: These must be the scaled values so we don't need to scale!
}


void EpetraModelEvaluator::setStateVariableScalingVec(
  const RCP<const Epetra_Vector> &stateVariableScalingVec
  )
{
  typedef ModelEvaluatorBase MEB;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !this->createInArgs().supports(MEB::IN_ARG_x) );
#endif  
  stateVariableScalingVec_ = stateVariableScalingVec.assert_not_null();
  invStateVariableScalingVec_ = Teuchos::null;
  nominalValuesAndBoundsAreUpdated_ = false;
}


RCP<const Epetra_Vector>
EpetraModelEvaluator::getStateVariableScalingVec() const
{
  return stateVariableScalingVec_;
}


RCP<const Epetra_Vector>
EpetraModelEvaluator::getStateVariableInvScalingVec() const
{
  updateNominalValuesAndBounds();
  return invStateVariableScalingVec_;
}


void EpetraModelEvaluator::setStateFunctionScalingVec(
  const RCP<const Epetra_Vector> &stateFunctionScalingVec
  )
{
  stateFunctionScalingVec_ = stateFunctionScalingVec;
}


RCP<const Epetra_Vector>
EpetraModelEvaluator::getStateFunctionScalingVec() const
{
  return stateFunctionScalingVec_;
}


void EpetraModelEvaluator::uninitialize(
  RCP<const EpetraExt::ModelEvaluator> *epetraModel,
  RCP<LinearOpWithSolveFactoryBase<double> > *W_factory
  )
{
  if(epetraModel) *epetraModel = epetraModel_;
  if(W_factory) *W_factory = W_factory_;
  epetraModel_ = Teuchos::null;
  W_factory_ = Teuchos::null;
  stateFunctionScalingVec_ = Teuchos::null;
  stateVariableScalingVec_ = Teuchos::null;
  invStateVariableScalingVec_ = Teuchos::null;
  currentInArgsOutArgs_ = false;
}


const ModelEvaluatorBase::InArgs<double>&
EpetraModelEvaluator::getFinalPoint() const
{
  return finalPoint_;
}


bool EpetraModelEvaluator::finalPointWasSolved() const
{
  return finalPointWasSolved_;
}


// Public functions overridden from Teuchos::Describable


std::string EpetraModelEvaluator::description() const
{
  std::ostringstream oss;
  oss << "Thyra::EpetraModelEvaluator{";
  oss << "epetraModel=";
  if(epetraModel_.get())
    oss << "\'"<<epetraModel_->description()<<"\'";
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


// Overridden from Teuchos::ParameterListAcceptor


void EpetraModelEvaluator::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters(),0); // Just validate my params
  paramList_ = paramList;
  const EStateFunctionScaling stateFunctionScaling_old = stateFunctionScaling_; 
  stateFunctionScaling_ = stateFunctionScalingValidator->getIntegralValue(
    *paramList_, StateFunctionScaling_name, StateFunctionScaling_default
    );
  if( stateFunctionScaling_ != stateFunctionScaling_old )
    stateFunctionScalingVec_ = Teuchos::null;
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}


RCP<Teuchos::ParameterList>
EpetraModelEvaluator::getNonconstParameterList()
{
  return paramList_;
}


RCP<Teuchos::ParameterList>
EpetraModelEvaluator::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


RCP<const Teuchos::ParameterList>
EpetraModelEvaluator::getParameterList() const
{
  return paramList_;
}


RCP<const Teuchos::ParameterList>
EpetraModelEvaluator::getValidParameters() const
{
  using Teuchos::rcp;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::tuple;
  using Teuchos::rcp_implicit_cast;
  typedef Teuchos::ParameterEntryValidator PEV;
  static RCP<const Teuchos::ParameterList> validPL;
  if(is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    stateFunctionScalingValidator = rcp(
      new StringToIntegralParameterEntryValidator<EStateFunctionScaling>(
        tuple<std::string>(
          "None",
          "Row Sum"
          ),
        tuple<std::string>(
          "Do not scale the state function f(...) in this class.",

          "Scale the state function f(...) and all its derivatives\n"
          "using the row sum scaling from the initial Jacobian\n"
          "W=d(f)/d(x).  Note, this only works with Epetra_CrsMatrix\n"
          "currently."
          ),
        tuple<EStateFunctionScaling>(
          STATE_FUNC_SCALING_NONE,
          STATE_FUNC_SCALING_ROW_SUM
          ),
        StateFunctionScaling_name
        )
      );
    pl->set(StateFunctionScaling_name,StateFunctionScaling_default,
      "Determines if and how the state function f(...) and all of its\n"
      "derivatives are scaled.  The scaling is done explicitly so there should\n"
      "be no impact on the meaning of inner products or tolerances for\n"
      "linear solves.",
      rcp_implicit_cast<const PEV>(stateFunctionScalingValidator)
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from ModelEvaulator.


int EpetraModelEvaluator::Np() const
{
  return p_space_.size();
}


int EpetraModelEvaluator::Ng() const
{
  return g_space_.size();
}


RCP<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_x_space() const
{
  return x_space_;
}


RCP<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_f_space() const
{
  return f_space_;
}


RCP<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_p_space(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, this->Np() );
#endif
  return p_space_[l];
}


RCP<const Teuchos::Array<std::string> >
EpetraModelEvaluator::get_p_names(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, this->Np() );
#endif
  return epetraModel_->get_p_names(l);
}


RCP<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_g_space(int j) const
{
  TEST_FOR_EXCEPT( ! ( 0 <= j && j < this->Ng() ) );
  return g_space_[j];
}


ModelEvaluatorBase::InArgs<double>
EpetraModelEvaluator::getNominalValues() const
{
  updateNominalValuesAndBounds();
  return nominalValues_;
}


ModelEvaluatorBase::InArgs<double>
EpetraModelEvaluator::getLowerBounds() const
{
  updateNominalValuesAndBounds();
  return lowerBounds_;
}


ModelEvaluatorBase::InArgs<double>
EpetraModelEvaluator::getUpperBounds() const
{
  updateNominalValuesAndBounds();
  return upperBounds_;
}


RCP<LinearOpWithSolveBase<double> >
EpetraModelEvaluator::create_W() const
{
  TEST_FOR_EXCEPTION(
    W_factory_.get()==NULL, std::logic_error
    ,"Thyra::EpetraModelEvaluator::create_W(): "
    "Error, the client did not set a LinearOpWithSolveFactoryBase"
    " object for W!"
    );
  W_factory_->setOStream(this->getOStream());
  W_factory_->setVerbLevel(this->getVerbLevel());
  return W_factory_->createOp();
}


RCP<LinearOpBase<double> >
EpetraModelEvaluator::create_W_op() const
{
  return this->create_epetra_W_op();
}


RCP<const LinearOpWithSolveFactoryBase<double> >
EpetraModelEvaluator::get_W_factory() const
{
  return W_factory_;
}


ModelEvaluatorBase::InArgs<double> EpetraModelEvaluator::createInArgs() const
{
  if (!currentInArgsOutArgs_)
    updateInArgsOutArgs();
  return prototypeInArgs_;
}


void EpetraModelEvaluator::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<double> &finalPoint,
  const bool wasSolved
  )
{
  finalPoint_ = this->createInArgs();
  finalPoint_.setArgs(finalPoint);
  finalPointWasSolved_ = wasSolved;
}


// Private functions overridden from ModelEvaulatorDefaultBase


RCP<LinearOpBase<double> >
EpetraModelEvaluator::create_DfDp_op_impl(int l) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


RCP<LinearOpBase<double> >
EpetraModelEvaluator::create_DgDx_dot_op_impl(int j) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


RCP<LinearOpBase<double> >
EpetraModelEvaluator::create_DgDx_op_impl(int j) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


RCP<LinearOpBase<double> >
EpetraModelEvaluator::create_DgDp_op_impl( int j, int l ) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


ModelEvaluatorBase::OutArgs<double>
EpetraModelEvaluator::createOutArgsImpl() const
{
  if (!currentInArgsOutArgs_)
    updateInArgsOutArgs();
  return prototypeOutArgs_;
}


void EpetraModelEvaluator::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double>& inArgs_in,
  const ModelEvaluatorBase::OutArgs<double>& outArgs
  ) const
{

  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::includesVerbLevel;
  typedef EpetraExt::ModelEvaluator EME;

  //
  // A) Initial setup
  //

  // Make sure that we are fully initialized!
  this->updateNominalValuesAndBounds();

  // Make sure we grab the initial guess first!
  InArgs<double> inArgs = this->getNominalValues();
  // Now, copy the parameters from the input inArgs_in object to the inArgs
  // object.  Any input objects that are set in inArgs_in will overwrite those
  // in inArgs that will already contain the nominal values.  This will insure
  // that all input parameters are set and those that are not set by the
  // client will be at their nominal values (as defined by the underlying
  // EpetraExt::ModelEvaluator object).  The full set of Thyra input arguments
  // must be set before these can be translated into Epetra input arguments.
  inArgs.setArgs(inArgs_in);

  // Print the header and the values of the inArgs and outArgs objects!
  typedef double Scalar; // Needed for below macro!
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Thyra::EpetraModelEvaluator",inArgs,outArgs,Teuchos::null
    );

  // State function Scaling
  const bool firstTimeStateFuncScaling
    = (
      stateFunctionScaling_ != STATE_FUNC_SCALING_NONE
      && is_null(stateFunctionScalingVec_)
      );
  
  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveFactoryBase<double> > VOTSLOWSF;
  VOTSLOWSF W_factory_outputTempState(W_factory_,out,verbLevel);

  Teuchos::Time timer("");

  //
  // B) Prepressess the InArgs and OutArgs in preparation to call
  // the underlying EpetraExt::ModelEvaluator
  //

  //
  // B.1) Translate InArgs from Thyra to Epetra objects
  //
  
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nSetting-up/creating input arguments ...\n";
  timer.start(true);

  // Unwrap input Thyra objects to get Epetra objects
  EME::InArgs epetraScaledInArgs = epetraModel_->createInArgs();
  convertInArgsFromThyraToEpetra( inArgs, &epetraScaledInArgs );

  // Unscale the input Epetra objects which will be passed to the underlying
  // EpetraExt::ModelEvaluator object.
  EME::InArgs epetraInArgs = epetraModel_->createInArgs();
  EpetraExt::unscaleModelVars(
    epetraScaledInArgs, epetraInArgsScaling_, &epetraInArgs,
    out.get(), verbLevel
    );

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o() << "\nTime to setup InArgs = "<<timer.totalElapsedTime()<<" sec\n";

  //
  // B.2) Convert from Thyra to Epetra OutArgs
  //
  
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nSetting-up/creating output arguments ...\n";
  timer.start(true);
  
  // The unscaled Epetra OutArgs that will be passed to the
  // underlying EpetraExt::ModelEvaluator object
  EME::OutArgs epetraUnscaledOutArgs = epetraModel_->createOutArgs();

  // Various objects that are needed later (see documentation in
  // the function convertOutArgsFromThyraToEpetra(...)
  RCP<LinearOpWithSolveBase<double> > W;
  RCP<LinearOpBase<double> > W_op;
  RCP<const LinearOpBase<double> > fwdW;
  RCP<EpetraLinearOp> efwdW;
  RCP<Epetra_Operator> eW;
  
  // Convert from Thyra to Epetra OutArgs and grap some of the intermediate
  // objects accessed along the way that are needed later.
  convertOutArgsFromThyraToEpetra(
    outArgs,
    &epetraUnscaledOutArgs,
    &W, &W_op, &fwdW, &efwdW, &eW
    );

  //
  // B.3) Setup OutArgs to computing scaling if needed
  //

  if (firstTimeStateFuncScaling) {
    preEvalScalingSetup(&epetraInArgs,&epetraUnscaledOutArgs,out,verbLevel);
  }

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o()
      << "\nTime to setup OutArgs = "
      << timer.totalElapsedTime() <<" sec\n";

  //
  // C) Evaluate the underlying EpetraExt model to compute the Epetra outputs
  //
  
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nEvaluating the Epetra output functions ...\n";
  timer.start(true);

  epetraModel_->evalModel(epetraInArgs,epetraUnscaledOutArgs);

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o()
      << "\nTime to evaluate Epetra output functions = "
      << timer.totalElapsedTime() <<" sec\n";

  //
  // D) Postprocess the output objects
  //

  //
  // D.1) Compute the scaling factors if needed
  //
  
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nCompute scale factors if needed ...\n";
  timer.start(true);

  if (firstTimeStateFuncScaling) {
    postEvalScalingSetup(epetraUnscaledOutArgs,out,verbLevel);
  }

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o()
      << "\nTime to compute scale factors = "
      << timer.totalElapsedTime() <<" sec\n";

  //
  // D.2) Scale the output Epetra objects
  //

  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nScale the output objects ...\n";
  timer.start(true);

  EME::OutArgs epetraOutArgs = epetraModel_->createOutArgs();
  bool allFuncsWhereScaled = false;
  EpetraExt::scaleModelFuncs(
    epetraUnscaledOutArgs, epetraInArgsScaling_, epetraOutArgsScaling_,
    &epetraOutArgs, &allFuncsWhereScaled,
    out.get(), verbLevel
    );
  TEST_FOR_EXCEPTION(
    !allFuncsWhereScaled, std::logic_error,
    "Error, we can not currently handle epetra output objects that could not be"
    " scaled.  Special code will have to be added to handle this (i.e. using"
    " implicit diagonal and multiplied linear operators to implicitly do"
    " the scaling."
    );

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o()
      << "\nTime to scale the output objects = "
      << timer.totalElapsedTime() << " sec\n";

  //
  // D.3) Convert any Epetra objects to Thyra OutArgs objects that still need to
  // be converted
  //

  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nFinish processing and wrapping the output objects ...\n";
  timer.start(true);

  finishConvertingOutArgsFromEpetraToThyra(
    epetraOutArgs, W, W_op, fwdW, efwdW, eW,
    outArgs
    );

  timer.stop();
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    OSTab(out).o()
      << "\nTime to finish processing and wrapping the output objects = "
      << timer.totalElapsedTime() <<" sec\n";

  //
  // E) Print footer to end the function
  //

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


// private


void EpetraModelEvaluator::convertInArgsFromEpetraToThyra(
  const EpetraExt::ModelEvaluator::InArgs &epetraInArgs,
  ModelEvaluatorBase::InArgs<double> *inArgs
  ) const
{
  
  using Teuchos::implicit_cast;
  typedef ModelEvaluatorBase MEB;

  TEST_FOR_EXCEPT(!inArgs);

  if(inArgs->supports(MEB::IN_ARG_x)) {
    inArgs->set_x( create_Vector( epetraInArgs.get_x(), x_space_ ) );
  }
  
  if(inArgs->supports(MEB::IN_ARG_x_dot)) {
    inArgs->set_x_dot( create_Vector( epetraInArgs.get_x_dot(), x_space_ ) );
  }

  const int l_Np = inArgs->Np();
  for( int l = 0; l < l_Np; ++l ) {
    inArgs->set_p( l, create_Vector( epetraInArgs.get_p(l), p_space_[l] ) );
  }
  
  if(inArgs->supports(MEB::IN_ARG_t)) {
    inArgs->set_t(epetraInArgs.get_t());
  }
  
}


void EpetraModelEvaluator::convertInArgsFromThyraToEpetra(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  EpetraExt::ModelEvaluator::InArgs *epetraInArgs
  ) const
{

  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  using Teuchos::Polynomial;
#endif // HAVE_THYRA_ME_POLYNOMIAL


  TEST_FOR_EXCEPT(0==epetraInArgs);

  RCP<const VectorBase<double> > x_dot;
  if( inArgs.supports(IN_ARG_x_dot) && (x_dot = inArgs.get_x_dot()).get() ) {
    RCP<const Epetra_Vector> e_x_dot = get_Epetra_Vector(*x_map_,x_dot);
    epetraInArgs->set_x_dot(e_x_dot);
  }

  RCP<const VectorBase<double> > x;
  if( inArgs.supports(IN_ARG_x) && (x = inArgs.get_x()).get() ) {
    RCP<const Epetra_Vector> e_x = get_Epetra_Vector(*x_map_,x);
    epetraInArgs->set_x(e_x);
  }

  RCP<const VectorBase<double> > p_l;
  for(int l = 0;  l < inArgs.Np(); ++l ) {
    p_l = inArgs.get_p(l);
    if(p_l.get()) epetraInArgs->set_p(l,get_Epetra_Vector(*p_map_[l],p_l));
  }

#ifdef HAVE_THYRA_ME_POLYNOMIAL

  RCP<const Polynomial< VectorBase<double> > > x_dot_poly;
  RCP<Epetra_Vector> epetra_ptr;
  if(
    inArgs.supports(IN_ARG_x_dot_poly)
    && (x_dot_poly = inArgs.get_x_dot_poly()).get()
    )
  {
    RCP<Polynomial<Epetra_Vector> > epetra_x_dot_poly = 
      rcp(new Polynomial<Epetra_Vector>(x_dot_poly->degree()));
    for (unsigned int i=0; i<=x_dot_poly->degree(); i++) {
      epetra_ptr = rcp_const_cast<Epetra_Vector>(
        get_Epetra_Vector(*x_map_, x_dot_poly->getCoefficient(i)) );
      epetra_x_dot_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraInArgs->set_x_dot_poly(epetra_x_dot_poly);
  }
  
  RCP<const Polynomial< VectorBase<double> > > x_poly;
  if(
    inArgs.supports(IN_ARG_x_poly)
    && (x_poly = inArgs.get_x_poly()).get()
    )
  {
    RCP<Polynomial<Epetra_Vector> > epetra_x_poly = 
      rcp(new Polynomial<Epetra_Vector>(x_poly->degree()));
    for (unsigned int i=0; i<=x_poly->degree(); i++) {
      epetra_ptr = rcp_const_cast<Epetra_Vector>(
        get_Epetra_Vector(*x_map_, x_poly->getCoefficient(i)) );
      epetra_x_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraInArgs->set_x_poly(epetra_x_poly);
  }

#endif // HAVE_THYRA_ME_POLYNOMIAL

  if( inArgs.supports(IN_ARG_t) )
    epetraInArgs->set_t(inArgs.get_t());
  
  if( inArgs.supports(IN_ARG_alpha) )
    epetraInArgs->set_alpha(inArgs.get_alpha());
  
  if( inArgs.supports(IN_ARG_beta) )
    epetraInArgs->set_beta(inArgs.get_beta());

}


void EpetraModelEvaluator::convertOutArgsFromThyraToEpetra(
  const ModelEvaluatorBase::OutArgs<double> &outArgs,
  EpetraExt::ModelEvaluator::OutArgs *epetraUnscaledOutArgs_inout,
  RCP<LinearOpWithSolveBase<double> > *W_out,
  RCP<LinearOpBase<double> > *W_op_out,
  RCP<const LinearOpBase<double> > *fwdW_out,
  RCP<EpetraLinearOp> *efwdW_out,
  RCP<Epetra_Operator> *eW_out
  ) const
{

  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::implicit_cast;
  using Thyra::get_Epetra_Vector;
  typedef EpetraExt::ModelEvaluator EME;

  // Assert input
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(epetraUnscaledOutArgs_inout);
  TEUCHOS_ASSERT(W_out);
  TEUCHOS_ASSERT(W_op_out);
  TEUCHOS_ASSERT(fwdW_out);
  TEUCHOS_ASSERT(efwdW_out);
  TEUCHOS_ASSERT(eW_out);
#endif

  // Create easy to use references
  EpetraExt::ModelEvaluator::OutArgs &epetraUnscaledOutArgs = *epetraUnscaledOutArgs_inout;
  RCP<LinearOpWithSolveBase<double> > &W = *W_out;
  RCP<LinearOpBase<double> > &W_op = *W_op_out;
  RCP<const LinearOpBase<double> > &fwdW = *fwdW_out;
  RCP<EpetraLinearOp> &efwdW = *efwdW_out;
  RCP<Epetra_Operator> &eW = *eW_out;

  // f
  { 
    RCP<VectorBase<double> > f;
    if( outArgs.supports(OUT_ARG_f) && (f = outArgs.get_f()).get() )
      epetraUnscaledOutArgs.set_f(get_Epetra_Vector(*f_map_,f));
  }
    
  // g
  {
    RCP<VectorBase<double> > g_j;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      g_j = outArgs.get_g(j);
      if(g_j.get()) epetraUnscaledOutArgs.set_g(j,get_Epetra_Vector(*g_map_[j],g_j));
    }
  }
  
  // W and W_op
  {

    if( outArgs.supports(OUT_ARG_W) && (W = outArgs.get_W()).get() ) {
      Thyra::uninitializeOp<double>(*W_factory_, W.ptr(), Teuchos::outArg(fwdW));
      if(fwdW.get()) {
        efwdW = rcp_const_cast<EpetraLinearOp>(
          rcp_dynamic_cast<const EpetraLinearOp>(fwdW,true));
      }
      else {
        efwdW = this->create_epetra_W_op();
        fwdW = efwdW;
      }
    }
    if( outArgs.supports(OUT_ARG_W_op) && (W_op = outArgs.get_W_op()).get() ) {
      if( W_op.get() && !efwdW.get() )
        efwdW = rcp_const_cast<EpetraLinearOp>(
          rcp_dynamic_cast<const EpetraLinearOp>(W_op,true));
    }
    if(efwdW.get()) {
      // By the time we get here, if we have an object in efwdW, then it
      // should already be embeadded with an underlying Epetra_Operator object
      // that was allocated by the EpetraExt::ModelEvaluator object.
      // Therefore, we should just have to grab this object and be on our way.
      eW = efwdW->epetra_op();
      epetraUnscaledOutArgs.set_W(eW);
    }
    // NOTE: Above, if both W and W_op are set and have been through at least
    // one prior evaluation (and therefore have Epetra_Operator objects embedded
    // in them), then we will use the Epetra_Operator embedded in W to pass to
    // the EpetraExt::ModelEvaluator object and ignore the Epetra_Operator
    // object in W_op.  In the standard use case, these will be the same
    // Epetra_Operator objects.  However, it is possible that the client could
    // use this interface in such a way that these would have different
    // Epetra_Operator objects embedded in them.  In this (very unlikely) case,
    // the Epetra_Operator embedded in W_op will be discarded!  This might be
    // surprising to a client but it is very unlikely that this will ever be a
    // problem, but the issue is duly noted here!  Only dangerous programming
    // use of this interface would cause any problem.
    
    // Note: The following derivative objects update in place!

  }

  // DfDp
  {
    Derivative<double> DfDp_l;
    for(int l = 0;  l < outArgs.Np(); ++l ) {
      if( !outArgs.supports(OUT_ARG_DfDp,l).none()
        && !(DfDp_l = outArgs.get_DfDp(l)).isEmpty() )
      {
        epetraUnscaledOutArgs.set_DfDp(l,convert(DfDp_l,f_map_,p_map_[l]));
      }
    }
  }

  // DgDx_dot
  {
    Derivative<double> DgDx_dot_j;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      if( !outArgs.supports(OUT_ARG_DgDx_dot,j).none()
        && !(DgDx_dot_j = outArgs.get_DgDx_dot(j)).isEmpty() )
      {
        epetraUnscaledOutArgs.set_DgDx_dot(j,convert(DgDx_dot_j,g_map_[j],x_map_));
      }
    }
  }

  // DgDx
  {
    Derivative<double> DgDx_j;
    for(int j = 0;  j < outArgs.Ng(); ++j ) {
      if( !outArgs.supports(OUT_ARG_DgDx,j).none()
        && !(DgDx_j = outArgs.get_DgDx(j)).isEmpty() )
      {
        epetraUnscaledOutArgs.set_DgDx(j,convert(DgDx_j,g_map_[j],x_map_));
      }
    }
  }

  // DgDp
  {
    DerivativeSupport DgDp_j_l_support;
    Derivative<double> DgDp_j_l;
    for (int j = 0;  j < outArgs.Ng(); ++j ) {
      for (int l = 0;  l < outArgs.Np(); ++l ) {
        if (!(DgDp_j_l_support = outArgs.supports(OUT_ARG_DgDp,j,l)).none()
          && !(DgDp_j_l = outArgs.get_DgDp(j,l)).isEmpty() )
        {
          epetraUnscaledOutArgs.set_DgDp(j,l,convert(DgDp_j_l,g_map_[j],p_map_[l]));
        }
      }
    }
  }

#ifdef HAVE_THYRA_ME_POLYNOMIAL

  // f_poly
  RCP<const Teuchos::Polynomial< VectorBase<double> > > f_poly;
  if( outArgs.supports(OUT_ARG_f_poly) && (f_poly = outArgs.get_f_poly()).get() )
  {
    RCP<Teuchos::Polynomial<Epetra_Vector> > epetra_f_poly = 
      Teuchos::rcp(new Teuchos::Polynomial<Epetra_Vector>(f_poly->degree()));
    for (unsigned int i=0; i<=f_poly->degree(); i++) {
      RCP<Epetra_Vector> epetra_ptr
        = Teuchos::rcp_const_cast<Epetra_Vector>(get_Epetra_Vector(*f_map_,
            f_poly->getCoefficient(i)));
      epetra_f_poly->setCoefficientPtr(i,epetra_ptr);
    }
    epetraUnscaledOutArgs.set_f_poly(epetra_f_poly);
  }

#endif // HAVE_THYRA_ME_POLYNOMIAL

}


void EpetraModelEvaluator::preEvalScalingSetup(
  EpetraExt::ModelEvaluator::InArgs *epetraInArgs_inout,
  EpetraExt::ModelEvaluator::OutArgs *epetraUnscaledOutArgs_inout,
  const RCP<Teuchos::FancyOStream> &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  
  typedef EpetraExt::ModelEvaluator EME;
  
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(epetraInArgs_inout);
  TEUCHOS_ASSERT(epetraUnscaledOutArgs_inout);
#endif

  EpetraExt::ModelEvaluator::InArgs
    &epetraInArgs = *epetraInArgs_inout;
  EpetraExt::ModelEvaluator::OutArgs
    &epetraUnscaledOutArgs = *epetraUnscaledOutArgs_inout;

  if (
    ( stateFunctionScaling_ == STATE_FUNC_SCALING_ROW_SUM )
    &&
    (
      epetraUnscaledOutArgs.supports(EME::OUT_ARG_f) 
      &&
      epetraUnscaledOutArgs.funcOrDerivesAreSet(EME::OUT_ARG_f)
      )
    &&
    (
      epetraUnscaledOutArgs.supports(EME::OUT_ARG_W)
      &&
      is_null(epetraUnscaledOutArgs.get_W())
      )
    )
  {
    // This is the first pass through with scaling turned on and the client
    // turned on automatic scaling but did not ask for W.  We must compute W
    // in order to compute the scale factors so we must allocate a temporary W
    // just to compute the scale factors and then throw it away.  If the
    // client wants to evaluate W at the same point, then it should have
    // passed W in but that is not our problem here.  The ModelEvaluator
    // relies on the client to set up the calls to allow for efficient
    // evaluation.

    if(out.get() && verbLevel >= Teuchos::VERB_LOW)
      *out
        << "\nCreating a temporary Epetra W to compute scale factors"
        << " for f(...) ...\n";
    epetraUnscaledOutArgs.set_W(create_and_assert_W(*epetraModel_));
    if( epetraInArgs.supports(EME::IN_ARG_beta) )
      epetraInArgs.set_beta(1.0);
    if( epetraInArgs.supports(EME::IN_ARG_alpha) )
      epetraInArgs.set_alpha(0.0);
  }
  
}


void EpetraModelEvaluator::postEvalScalingSetup(
  const EpetraExt::ModelEvaluator::OutArgs &epetraUnscaledOutArgs,
  const RCP<Teuchos::FancyOStream> &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::OSTab;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::includesVerbLevel;

  // Compute the scale factors for the state function f(...)
  switch(stateFunctionScaling_) {

    case STATE_FUNC_SCALING_ROW_SUM: {

      // Compute the inverse row-sum scaling from W

      const RCP<Epetra_RowMatrix>
        ermW = get_Epetra_RowMatrix(epetraUnscaledOutArgs);
      // Note: Above, we get the Epetra W object directly from the Epetra
      // OutArgs object since this might be a temporary matrix just to
      // compute scaling factors.  In this case, the stack funtion variable
      // eW might be empty!

      RCP<Epetra_Vector>
        invRowSums = rcp(new Epetra_Vector(ermW->OperatorRangeMap()));
      // Above: From the documentation is seems that the RangeMap should be
      // okay but who knows for sure!

      ermW->InvRowSums(*invRowSums);

      if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW)) {
        *out
          << "\nComputed inverse row sum scaling from W that"
          " will be used to scale f(...) and its derivatives:\n";
        double minVal = 0, maxVal = 0, avgVal = 0;
        invRowSums->MinValue(&minVal);
        invRowSums->MaxValue(&maxVal);
        invRowSums->MeanValue(&avgVal);
        OSTab tab(out);
        *out
          << "min(invRowSums) = " << minVal << "\n"
          << "max(invRowSums) = " << maxVal << "\n"
          << "avg(invRowSums) = " << avgVal << "\n";
      }

      stateFunctionScalingVec_ = invRowSums;

      break;

    }

    default:
      TEST_FOR_EXCEPT("Should never get here!");

  }

  epetraOutArgsScaling_ = epetraModel_->createOutArgs();

  epetraOutArgsScaling_.set_f(
    rcp_const_cast<Epetra_Vector>(stateFunctionScalingVec_) );

}


void EpetraModelEvaluator::finishConvertingOutArgsFromEpetraToThyra(
  const EpetraExt::ModelEvaluator::OutArgs &epetraOutArgs,
  RCP<LinearOpWithSolveBase<double> > &W,
  RCP<LinearOpBase<double> > &W_op,
  RCP<const LinearOpBase<double> > &fwdW,
  RCP<EpetraLinearOp> &efwdW,
  RCP<Epetra_Operator> &eW,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  typedef EpetraExt::ModelEvaluator EME;

  if(efwdW.get()) {
    efwdW->setFullyInitialized(true); 
    // NOTE: Above will directly update W_op also if W.get()==NULL!
  }
  
  if( W.get() ) {
    Thyra::initializeOp<double>(*W_factory_, fwdW, W.ptr());
    W->setOStream(this->getOStream());
  }

  if( W_op.get() ) {
    if( W_op.shares_resource(efwdW) ) {
      // W_op was already updated above since *efwdW is the same object as *W_op
    }
    else {
      rcp_dynamic_cast<EpetraLinearOp>(W_op,true)->setFullyInitialized(true);
    }
  }

}


void EpetraModelEvaluator::updateNominalValuesAndBounds() const
{

  using Teuchos::rcp;
  using Teuchos::implicit_cast;
  typedef ModelEvaluatorBase MEB;
  typedef EpetraExt::ModelEvaluator EME;

  if( !nominalValuesAndBoundsAreUpdated_ ) {

    // Gather the nominal values and bounds into Epetra InArgs objects

    EME::InArgs epetraOrigNominalValues;
    EpetraExt::gatherModelNominalValues(
      *epetraModel_, &epetraOrigNominalValues );

    EME::InArgs epetraOrigLowerBounds;
    EME::InArgs epetraOrigUpperBounds;
    EpetraExt::gatherModelBounds(
      *epetraModel_, &epetraOrigLowerBounds, &epetraOrigUpperBounds );

    // Set up Epetra InArgs scaling object

    epetraInArgsScaling_ = epetraModel_->createInArgs();

    if( !is_null(stateVariableScalingVec_) ) {
      invStateVariableScalingVec_
        = EpetraExt::createInverseModelScalingVector(stateVariableScalingVec_);
      if( epetraOrigNominalValues.supports(EME::IN_ARG_x_dot) ) {
        epetraInArgsScaling_.set_x_dot(invStateVariableScalingVec_);
      }
      if( epetraOrigNominalValues.supports(EME::IN_ARG_x) ) {
        epetraInArgsScaling_.set_x(invStateVariableScalingVec_);
      }
    }
    
    // Scale the original variables and bounds

    EME::InArgs epetraScaledNominalValues = epetraModel_->createInArgs();
    EpetraExt::scaleModelVars(
      epetraOrigNominalValues, epetraInArgsScaling_, &epetraScaledNominalValues
      );

    EME::InArgs epetraScaledLowerBounds = epetraModel_->createInArgs();
    EME::InArgs epetraScaledUpperBounds = epetraModel_->createInArgs();
    EpetraExt::scaleModelBounds(
      epetraOrigLowerBounds, epetraOrigUpperBounds, epetraModel_->getInfBound(),
      epetraInArgsScaling_,
      &epetraScaledLowerBounds, &epetraScaledUpperBounds
      );

    // Wrap the scaled epetra InArgs objects as Thyra InArgs objects!

    nominalValues_ = this->createInArgs();
    lowerBounds_ = this->createInArgs();
    upperBounds_ = this->createInArgs();
    convertInArgsFromEpetraToThyra(epetraScaledNominalValues, &nominalValues_);
    convertInArgsFromEpetraToThyra(epetraScaledLowerBounds, &lowerBounds_);
    convertInArgsFromEpetraToThyra(epetraScaledUpperBounds, &upperBounds_);

    nominalValuesAndBoundsAreUpdated_ = true;

  }
  else {

    // The nominal values and bounds should already be updated an should have
    // the currect scaling!

  }

}


void EpetraModelEvaluator::updateInArgsOutArgs() const
{

  typedef EpetraExt::ModelEvaluator EME;

  const EpetraExt::ModelEvaluator &epetraModel = *epetraModel_;
  EME::InArgs  epetraInArgs  = epetraModel.createInArgs();
  EME::OutArgs epetraOutArgs = epetraModel.createOutArgs();
  const int l_Np = epetraOutArgs.Np();
  const int l_Ng = epetraOutArgs.Ng();

  //
  // InArgs
  //

  InArgsSetup<double> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(epetraInArgs.Np());
  inArgs.setSupports(IN_ARG_x_dot, epetraInArgs.supports(EME::IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, epetraInArgs.supports(EME::IN_ARG_x));
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  inArgs.setSupports(IN_ARG_x_dot_poly,
    epetraInArgs.supports(EME::IN_ARG_x_dot_poly));
  inArgs.setSupports(IN_ARG_x_poly, epetraInArgs.supports(EME::IN_ARG_x_poly));
#endif // HAVE_THYRA_ME_POLYNOMIAL
  inArgs.setSupports(IN_ARG_t, epetraInArgs.supports(EME::IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, epetraInArgs.supports(EME::IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, epetraInArgs.supports(EME::IN_ARG_beta));
  prototypeInArgs_ = inArgs;

  //
  // OutArgs
  //

  OutArgsSetup<double> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(l_Np, l_Ng);
  // f
  outArgs.setSupports(OUT_ARG_f, epetraOutArgs.supports(EME::OUT_ARG_f));
  if (outArgs.supports(OUT_ARG_f)) {
    // W
    outArgs.setSupports(
      OUT_ARG_W, epetraOutArgs.supports(EME::OUT_ARG_W)&&!is_null(W_factory_));
    outArgs.setSupports(OUT_ARG_W_op,  epetraOutArgs.supports(EME::OUT_ARG_W));
    outArgs.set_W_properties(convert(epetraOutArgs.get_W_properties()));
    // DfDp
    for(int l=0; l<l_Np; ++l) {
      outArgs.setSupports(OUT_ARG_DfDp, l,
        convert(epetraOutArgs.supports(EME::OUT_ARG_DfDp, l)));
      if(!outArgs.supports(OUT_ARG_DfDp, l).none())
        outArgs.set_DfDp_properties(l,
          convert(epetraOutArgs.get_DfDp_properties(l)));
    }
  }
  // DgDx_dot and DgDx
  for(int j=0; j<l_Ng; ++j) {
    if (inArgs.supports(IN_ARG_x_dot))
      outArgs.setSupports(OUT_ARG_DgDx_dot, j,
        convert(epetraOutArgs.supports(EME::OUT_ARG_DgDx_dot, j)));
    if(!outArgs.supports(OUT_ARG_DgDx_dot, j).none())
      outArgs.set_DgDx_dot_properties(j,
        convert(epetraOutArgs.get_DgDx_dot_properties(j)));
    if (inArgs.supports(IN_ARG_x))
      outArgs.setSupports(OUT_ARG_DgDx, j,
        convert(epetraOutArgs.supports(EME::OUT_ARG_DgDx, j)));
    if(!outArgs.supports(OUT_ARG_DgDx, j).none())
      outArgs.set_DgDx_properties(j,
        convert(epetraOutArgs.get_DgDx_properties(j)));
  }
  // DgDp
  for(int j=0; j < l_Ng; ++j) for(int l=0; l < l_Np; ++l) {
    const EME::DerivativeSupport epetra_DgDp_j_l_support =
      epetraOutArgs.supports(EME::OUT_ARG_DgDp, j, l);
    outArgs.setSupports(OUT_ARG_DgDp, j, l,
      convert(epetra_DgDp_j_l_support));
    if(!outArgs.supports(OUT_ARG_DgDp, j, l).none())
      outArgs.set_DgDp_properties(j, l,
        convert(epetraOutArgs.get_DgDp_properties(j, l)));
  }
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  outArgs.setSupports(OUT_ARG_f_poly,
    epetraOutArgs.supports(EME::OUT_ARG_f_poly));
#endif // HAVE_THYRA_ME_POLYNOMIAL
  prototypeOutArgs_ = outArgs;

  // We are current!
  currentInArgsOutArgs_ = true;

}


RCP<EpetraLinearOp>
EpetraModelEvaluator::create_epetra_W_op() const
{
  return Thyra::partialNonconstEpetraLinearOp(
    this->get_f_space(), this->get_x_space(),
    create_and_assert_W(*epetraModel_)
    );
}


} // namespace Thyra


//
// Non-member utility functions
//


Teuchos::RCP<Thyra::EpetraModelEvaluator>
Thyra::epetraModelEvaluator(
  const RCP<const EpetraExt::ModelEvaluator> &epetraModel,
  const RCP<LinearOpWithSolveFactoryBase<double> > &W_factory
  )
{
  return Teuchos::rcp(new EpetraModelEvaluator(epetraModel,W_factory));
}


Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation
Thyra::convert(
  const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation &mvOrientation
  )
{
  switch(mvOrientation) {
    case EpetraExt::ModelEvaluator::DERIV_MV_BY_COL :
      return ModelEvaluatorBase::DERIV_MV_BY_COL;
    case EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW :
      return ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DERIV_MV_BY_COL; // Should never be called!
}


EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation
Thyra::convert(
  const ModelEvaluatorBase::EDerivativeMultiVectorOrientation &mvOrientation
  )
{
  switch(mvOrientation) {
    case ModelEvaluatorBase::DERIV_MV_BY_COL :
      return EpetraExt::ModelEvaluator::DERIV_MV_BY_COL;
    case ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW :
      return EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return EpetraExt::ModelEvaluator::DERIV_MV_BY_COL; // Should never be called!
}


Thyra::ModelEvaluatorBase::DerivativeProperties
Thyra::convert(
  const EpetraExt::ModelEvaluator::DerivativeProperties &derivativeProperties
  )
{
  ModelEvaluatorBase::EDerivativeLinearity linearity;
  switch(derivativeProperties.linearity) {
    case EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN;
      break;
    case EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_CONST;
      break;
    case  EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST:
      linearity = ModelEvaluatorBase::DERIV_LINEARITY_NONCONST;
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
  ModelEvaluatorBase::ERankStatus rank;
  switch(derivativeProperties.rank) {
    case EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN:
      rank = ModelEvaluatorBase::DERIV_RANK_UNKNOWN;
      break;
    case EpetraExt::ModelEvaluator::DERIV_RANK_FULL:
      rank = ModelEvaluatorBase::DERIV_RANK_FULL;
      break;
    case EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT:
      rank = ModelEvaluatorBase::DERIV_RANK_DEFICIENT;
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DerivativeProperties(
    linearity,rank,derivativeProperties.supportsAdjoint);
}


Thyra::ModelEvaluatorBase::DerivativeSupport
Thyra::convert(
  const EpetraExt::ModelEvaluator::DerivativeSupport &derivativeSupport
  )
{
  ModelEvaluatorBase::DerivativeSupport ds;
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_LINEAR_OP))
    ds.plus(ModelEvaluatorBase::DERIV_LINEAR_OP);
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_MV_BY_COL))
    ds.plus(ModelEvaluatorBase::DERIV_MV_BY_COL);
  if(derivativeSupport.supports(EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW))
    ds.plus(ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW);
  return ds;
}


EpetraExt::ModelEvaluator::Derivative
Thyra::convert(
  const ModelEvaluatorBase::Derivative<double> &derivative,
  const RCP<const Epetra_Map> &fnc_map,
  const RCP<const Epetra_Map> &var_map
  )
{
  typedef ModelEvaluatorBase MEB;
  if(derivative.getLinearOp().get()) {
    return EpetraExt::ModelEvaluator::Derivative(
      Teuchos::rcp_const_cast<Epetra_Operator>(
        Teuchos::dyn_cast<const EpetraLinearOp>(*derivative.getLinearOp()).epetra_op()
        )
      );
  }
  else if(derivative.getDerivativeMultiVector().getMultiVector().get()) {
    return EpetraExt::ModelEvaluator::Derivative(
      EpetraExt::ModelEvaluator::DerivativeMultiVector(
        get_Epetra_MultiVector(
          ( derivative.getDerivativeMultiVector().getOrientation() == MEB::DERIV_MV_BY_COL
            ? *fnc_map
            : *var_map
            )
          ,derivative.getDerivativeMultiVector().getMultiVector()
          )
        ,convert(derivative.getDerivativeMultiVector().getOrientation())
        )
      );
  }
  return EpetraExt::ModelEvaluator::Derivative();
}
