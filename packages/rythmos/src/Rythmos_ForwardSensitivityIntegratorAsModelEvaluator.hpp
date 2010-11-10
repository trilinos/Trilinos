//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_FORWARD_SENSITIVITY_INTEGRATOR_AS_MODEL_EVALUATOR_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_INTEGRATOR_AS_MODEL_EVALUATOR_HPP


#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_InterpolationBufferHelpers.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_ForwardResponseSensitivityComputer.hpp"
#include "Rythmos_extractStateAndSens.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


namespace ForwardSensitivityIntegratorAsModelEvaluatorTypes {
/** \brief . */
enum EResponseType {
  /** \brief The response functions should be summed and returned. */
  RESPONSE_TYPE_SUM,
  /** \brief The response functions should be blocked up and returned. */
  RESPONSE_TYPE_BLOCK
};
} // namespace ForwardSensitivityIntegratorAsModelEvaluatorTypes


/** \brief Concrete <tt>Thyra::ModelEvaluator</tt> subclass that turns a
 * forward ODE/DAE with an observation into a parameterized evaluation of
 * <tt>p -> g</tt> with forward sensitivities <tt>DgDp</tt>.
 *
 * The form of the parameterized state equation is:

 \verbatim

   f(x_dot(t),x(t),{p_l},t) = 0, over t = [t0,tf]

   x(t0) = x_init(p)

 \endverbatim

 * The full response function takes one of two possible forms:
 * 
 * <ul>
 *
 * <li>Summation form (<tt>responseType==RESPONSE_TYPE_SUM</tt>):

 \verbatim

   d(x_dot,x,p) = sum( beta_k * g_k(x_dot(t_k),x(t_k),p,t_k), k=0...NO-1 )

 \endverbatim

 * <li>Blocked form (<tt>responseType==RESPONSE_TYPE_BLOCKED</tt>)

 \verbatim

   d(x_dot,x,p)[k] = g_k(x_dot(t_k),x(t_k),p,t_k), k=0...NO-1

 \endverbatim
 
 * </ul>
 *
 * The individual reduced response functions are:

 \verbatim

   g_hat_k(p) =  g_k(x_dot(t_k),x(t_k),p,t_k), k=0...NO-1
   
 \endverbatim

 * Above, the reduced response function <tt>g_hat_k(p)</tt> is the evaluation
 * of the full response function at the state solution.
 *
 * The first derivative of the individual reduced response functions is:

 \verbatim

   d(g_hat_k)/d(p) =
      beta_k
      *
      (
        d(g_k)/d(x_dot) * d(x_dot)/d(p)
        + d(g_k)/d(x) * d(x)/d(p)
        + d(g_k)/d(p)
        ),

    for k=0...NO-1

 \endverbatim

 * Above, <tt>S_dot = d(x_dot)/d(p)</tt> and <tt>S = d(x_dot)/d(p)</tt> are
 * computed by solving the forward sensitivity equations (see
 * <tt>ForwardSensitivityStepper</tt>).
 *
 * This class implements the reduced response function <tt>p -> d_hat(p)</tt>
 * by integrating the state+sens equations and assembling the reduced response
 * functions.
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class ForwardSensitivityIntegratorAsModelEvaluator
  : virtual public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
  , virtual public Teuchos::ParameterListAcceptor
{
public:

  //@}

  /** \name Constructors, Initialization, Misc. */
  //@{

  /** \brief . */
  ForwardSensitivityIntegratorAsModelEvaluator();

  /** \brief Initalized with all objects needed to define evaluation 
   *
   * \param stateStepper [inout,persisting] The state stepper.  On input, this
   * object must be fully initialized with everything it needs except for the
   * initial condition.
   *
   * \param stateIntegrator [inout,persisting] The integrator that will be used
   * to integrate the state equations.  This integrator need not be set up
   * with a stepper as that will be done internally.  It just needs to have
   * its parameters set to know which algorithm to run.  This integrator will
   * be cloned to create the integrator for the state + sens equations if
   * <tt>stateAndSensIntegrator==null</tt>.
   *
   * \param stateAndSensStepper [inout,persisting] The state + sens stepper.
   * On input, this object must be fully initialized with everything it needs
   * except for the initial condition.
   *
   * \param stateAndSensIntegrator [inout,persisting] The integrator that will
   * be used to integrate the state + sens equations.  This integrator need
   * not be set up with a stepper as that will be done internally.  It just
   * needs to have its parameters set to know which algorithm to run.  The
   * <tt>stateIntegrator</tt> will be cloned for this purpose if
   * <tt>stateAndSensIntegrator==null</tt>.
   *
   * \param stateAndSensInitCond [in,persisting] The initial condition for the
   * state + sens system.  This object should also include any parameters or
   * other data that is needed to define the initial condition.  The
   * underlying parameter subvector will be set and reset again and again in
   * order to satisfy the evaluation <tt>p -> g(p)</tt> defined by this
   * interface.
   *
   * \param responseTimes [in] Array giving the sorted time points for each
   * response function <tt>g_k(...)</tt>.  Note that the final time for the
   * integrator will be taken to be <tt>finalTime = responseTimes.back()</tt>.
   *
   * \param responseFuncs [in,persisting] Array giving the response functions
   * for each time point <tt>responseTimes[k]</tt>.  It is assumed that the
   * response function will be <tt>g(0)</tt> in this model and that the same
   * parameter subvector used in the state function is used here also.
   *
   * \param responseFuncBasePoints [in,persisting] Array giving the base point
   * of the response functions.
   *
   * \param responseType [in] Determines how the response functions are
   * treated and how they are returned.  If
   * <tt>responseType==RESPONSE_TYPE_SUM</tt>, then the response functions are
   * summed and the single summed response function is returned from
   * <tt>evalModel()</tt>.  If <tt>responseType==RESPONSE_TYPE_BLOCK</tt>, then
   * the response functions are blocked up and returned as a product vector
   * from <tt>evalModel()</tt>.
   *
   * ToDo: Finish Documentation!
   */
  void initialize(
    const RCP<StepperBase<Scalar> > &stateStepper,
    const RCP<IntegratorBase<Scalar> > &stateIntegrator,
    const RCP<ForwardSensitivityStepper<Scalar> > &stateAndSensStepper,
    const RCP<IntegratorBase<Scalar> > &stateAndSensIntegrator,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateAndSensInitCond,
    const Array<Scalar> &responseTimes,
    const Array<RCP<const Thyra::ModelEvaluator<Scalar> > > &responseFuncs,
    const Array<Thyra::ModelEvaluatorBase::InArgs<Scalar> > &responseFuncBasePoints,
    const ForwardSensitivityIntegratorAsModelEvaluatorTypes::EResponseType responseType
    );

  /** \brief . */
  const Array<Scalar>& getResponseTimes() const;

  //@}


  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief .
   *
   * Note that <tt>observationTargetIO()</tt> and <tt>parameterBaseIO()</tt>
   * must be set before calling this function in order to use the parameter
   * sublist to read in the vectors <tt>observationTarget()</tt> and
   * <tt>parameterBase()</tt>.
   */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief .
   *
   * Note that <tt>observationTargetIO()</tt> and <tt>parameterBaseIO()</tt>
   * must be set before calling this function in order to have the sublists
   * added that will allow the vectors <tt>observationTarget()</tt> and
   * <tt>parameterBase()</tt> to be read in latter when the parameter list is
   * set..
   */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  // //////////////////////
  // Private types

  typedef Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > SpaceArray_t;
  

  // //////////////////////
  // Private data members

  RCP<Teuchos::ParameterList>  paramList_;
  bool dumpSensitivities_;

  RCP<StepperBase<Scalar> > stateStepper_;
  RCP<IntegratorBase<Scalar> > stateIntegrator_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> stateInitCond_;

  RCP<ForwardSensitivityStepper<Scalar> > stateAndSensStepper_;
  RCP<IntegratorBase<Scalar> > stateAndSensIntegrator_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> stateAndSensInitCond_;

  Array<Scalar> responseTimes_;
  Array<RCP<const Thyra::ModelEvaluator<Scalar> > > responseFuncs_;
  mutable Array<Thyra::ModelEvaluatorBase::InArgs<Scalar> > responseFuncInArgs_;
  ForwardSensitivityIntegratorAsModelEvaluatorTypes::EResponseType responseType_;
  bool response_func_supports_x_dot_;
  bool response_func_supports_D_x_dot_;
  bool response_func_supports_D_p_;

  int p_index_;
  int g_index_;
  Scalar finalTime_;

  int Np_;
  int Ng_;

  SpaceArray_t g_space_;
  SpaceArray_t p_space_;

  static const std::string dumpSensitivities_name_;
  static const bool dumpSensitivities_default_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
RCP<ForwardSensitivityIntegratorAsModelEvaluator<Scalar> >
forwardSensitivityIntegratorAsModelEvaluator(
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<IntegratorBase<Scalar> > &stateIntegrator,
  const RCP<ForwardSensitivityStepper<Scalar> > &stateAndSensStepper,
  const RCP<IntegratorBase<Scalar> > &stateAndSensIntegrator,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateAndSensInitCond,
  const Array<Scalar> &responseTimes,
  const Array<RCP<const Thyra::ModelEvaluator<Scalar> > > &responseFuncs,
  const Array<Thyra::ModelEvaluatorBase::InArgs<Scalar> > &responseFuncBasePoints,
  const ForwardSensitivityIntegratorAsModelEvaluatorTypes::EResponseType responseType
  )
{
  using Teuchos::rcp;
  RCP<ForwardSensitivityIntegratorAsModelEvaluator<Scalar> >
    sensIntegratorAsModelEvaluator = rcp(new ForwardSensitivityIntegratorAsModelEvaluator<Scalar>());
  sensIntegratorAsModelEvaluator->initialize(
    stateStepper, stateIntegrator,
    stateAndSensStepper, stateAndSensIntegrator,
    stateAndSensInitCond,
    responseTimes, responseFuncs,
    responseFuncBasePoints, responseType
    );
  return sensIntegratorAsModelEvaluator;
}


// ////////////////////////
// Definitions


// Static members


template<class Scalar>
const std::string
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::dumpSensitivities_name_
= "Dump Sensitivities";

template<class Scalar>
const bool
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::dumpSensitivities_default_
= false;


// Constructors, Initialization, Misc.


template<class Scalar>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::ForwardSensitivityIntegratorAsModelEvaluator()
  :dumpSensitivities_(dumpSensitivities_default_),
   responseType_(ForwardSensitivityIntegratorAsModelEvaluatorTypes::RESPONSE_TYPE_SUM),
   response_func_supports_x_dot_(false),
   response_func_supports_D_x_dot_(false),
   response_func_supports_D_p_(false),
   p_index_(-1),
   g_index_(-1),
   finalTime_(-std::numeric_limits<Scalar>::max()),
   Np_(0),
   Ng_(0)
{}


template<class Scalar>
void ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::initialize(
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<IntegratorBase<Scalar> > &stateIntegrator,
  const RCP<ForwardSensitivityStepper<Scalar> > &stateAndSensStepper,
  const RCP<IntegratorBase<Scalar> > &stateAndSensIntegrator,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateAndSensInitCond,
  const Array<Scalar> &responseTimes,
  const Array<RCP<const Thyra::ModelEvaluator<Scalar> > > &responseFuncs,
  const Array<Thyra::ModelEvaluatorBase::InArgs<Scalar> > &responseFuncBasePoints,
  const ForwardSensitivityIntegratorAsModelEvaluatorTypes::EResponseType responseType
  )
{

  using Teuchos::as;
  typedef Thyra::ModelEvaluatorBase MEB;
  namespace FSIAMET = ForwardSensitivityIntegratorAsModelEvaluatorTypes;

  //
  // A) Validate and set input
  //

#ifdef RYTHMOS_DEBUG
  const int numResponseTimes = responseTimes.size();

  TEST_FOR_EXCEPT(is_null(stateStepper));
  TEST_FOR_EXCEPT(is_null(stateIntegrator));
  TEST_FOR_EXCEPT(is_null(stateAndSensStepper));
  TEST_FOR_EXCEPT(is_null(stateAndSensInitCond.get_x()));
  TEST_FOR_EXCEPT(is_null(stateAndSensInitCond.get_x_dot()));
  TEST_FOR_EXCEPT( !( numResponseTimes > 0 ) );
  assertTimePointsAreSorted(responseTimes);
  TEST_FOR_EXCEPT( as<int>(responseFuncs.size()) != numResponseTimes );
  TEST_FOR_EXCEPT( as<int>(responseFuncBasePoints.size()) != numResponseTimes );
  // ToDo: Assert that all of the observation models have the same response
  // function spaces so that they can be added together!
#endif // RYTHMOS_DEBUG

  stateStepper_ = stateStepper;
  stateAndSensStepper_ = stateAndSensStepper;
  stateAndSensInitCond_ = stateAndSensInitCond;
  stateIntegrator_ = stateIntegrator;
  if (!is_null(stateAndSensIntegrator))
    stateAndSensIntegrator_ = stateAndSensIntegrator;
  else
    stateAndSensIntegrator_ = stateIntegrator_->cloneIntegrator().assert_not_null();
  responseTimes_ = responseTimes;
  responseFuncs_ = responseFuncs;
  responseFuncInArgs_ = responseFuncBasePoints;
  responseType_ = responseType;
  
  finalTime_ = responseTimes_.back();
  
  //
  // B) Set the initial condition for the state-only problem
  //

  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    stateModel = stateStepper_->getModel();

  // Get base point pased in with no x or x_dot
  MEB::InArgs<Scalar>
    basePoint_no_x = stateAndSensInitCond_;
  basePoint_no_x.set_x(Teuchos::null);
  basePoint_no_x.set_x_dot(Teuchos::null);

  // Create an empty InArgs for the state model/stepper
  stateInitCond_ = stateModel->createInArgs();

  // Set the base point (except x, x_dot).
  stateInitCond_.setArgs(basePoint_no_x);

  // Set x and x_dot
  const RCP<const Thyra::ProductVectorBase<Scalar> >
    x_bar_init = Thyra::productVectorBase<Scalar>(
      stateAndSensInitCond_.get_x()
      ),
    x_bar_dot_init = Thyra::productVectorBase<Scalar>(
      stateAndSensInitCond_.get_x_dot()
      );
  stateInitCond_.set_x(x_bar_init->getVectorBlock(0));
  stateInitCond_.set_x_dot(x_bar_dot_init->getVectorBlock(0));

  //
  // C) Set up the info for this model evaluators interface
  //

  Np_ = 1;
  p_index_ = getParameterIndex(*stateAndSensStepper_);
  p_space_.clear();
  p_space_.push_back(stateModel->get_p_space(p_index_));
  
  Ng_ = 1;
  g_index_ = 0; // ToDo: Accept this from input!
  g_space_.clear();

  if (responseType_ == FSIAMET::RESPONSE_TYPE_SUM) {
    g_space_.push_back(responseFuncs[0]->get_g_space(0));
  }
  else if (responseType_ == FSIAMET::RESPONSE_TYPE_BLOCK) {
    g_space_.push_back(
      Thyra::productVectorSpace(
        responseFuncs[0]->get_g_space(g_index_), responseFuncs.size()
        )
      );
  }

  MEB::InArgs<Scalar>
    responseInArgs = responseFuncs[0]->createInArgs();
  response_func_supports_x_dot_ =
    responseInArgs.supports(MEB::IN_ARG_x_dot);
  MEB::OutArgs<Scalar>
    responseOutArgs = responseFuncs[0]->createOutArgs();
  response_func_supports_D_x_dot_ =
    !responseOutArgs.supports(MEB::OUT_ARG_DgDx_dot,g_index_).none();
  response_func_supports_D_p_ =
    !responseOutArgs.supports(MEB::OUT_ARG_DgDp,g_index_,p_index_).none();
  
}


template<class Scalar>
const Array<Scalar>&
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::getResponseTimes() const
{
  return responseTimes_;
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(0==paramList.get());
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  dumpSensitivities_ = paramList_->get(
    dumpSensitivities_name_, dumpSensitivities_default_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef RYTHMOS_DEBUG
  paramList_->validateParameters(*getValidParameters());
#endif // RYTHMOS_DEBUG

}


template<class Scalar>
RCP<Teuchos::ParameterList>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(dumpSensitivities_name_, dumpSensitivities_default_,
      "Set to true to have the individual sensitivities printed to\n"
      "the output stream."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
int ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::Np() const
{
  return Np_;
}


template<class Scalar>
int ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::Ng() const
{
  return Ng_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::get_p_space(int l) const
{
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_[l];
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::get_g_space(int j) const
{
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_[j];
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(Np_);
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeSupport DS;
  namespace FSIAMET = ForwardSensitivityIntegratorAsModelEvaluatorTypes;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np_,Ng_);
  outArgs.setSupports(MEB::OUT_ARG_DgDp, 0 ,0, MEB::DERIV_MV_BY_COL);
  return outArgs;
}


template<class Scalar>
void ForwardSensitivityIntegratorAsModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::describe;
  using Teuchos::OSTab;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::VerboseObjectTempState<InterpolationBufferBase<Scalar> > VOTSSB;
  typedef Thyra::ModelEvaluatorBase MEB;
  namespace FSIAMET = ForwardSensitivityIntegratorAsModelEvaluatorTypes;

  //
  // Stream output and other basic setup
  //
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::ForwardSensitivityIntegratorAsModelEvaluator", inArgs, outArgs, null
    );
  VOTSSB stateIntegrator_outputTempState(
    stateIntegrator_,out,incrVerbLevel(verbLevel,-1));
  VOTSSB stateAndSensIntegrator_outputTempState(
    stateAndSensIntegrator_,out,incrVerbLevel(verbLevel,-1));

  const bool trace =
    out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_LOW);
  const bool print_norms =
    out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_MEDIUM);
  const bool print_x =
    out.get() && includesVerbLevel(localVerbLevel,Teuchos::VERB_EXTREME);

  const RCP<const Thyra::ModelEvaluator<Scalar> >
    stateModel = stateStepper_->getModel();
  
  //
  // A) Process OutArgs first to see what functions we will be computing
  //
  
  RCP<Thyra::VectorBase<Scalar> >
    d_hat = outArgs.get_g(0);

  MEB::Derivative<Scalar>
    D_d_hat_D_p = outArgs.get_DgDp(0,p_index_);
    
  // Integrate state+sens or just state?
  const bool integrateStateAndSens = !D_d_hat_D_p.isEmpty();

  // Shortcut exit if no output was requested
  if ( is_null(d_hat) && D_d_hat_D_p.isEmpty() ) {
    if (trace)
      *out << "\nSkipping evaluation since no outputs were requested!\n";
    return;
  }
    
  //
  // B) Process the InArgs knowing if we will integrate just the state or the
  // state+sens.
  //
    
  const RCP<const Thyra::VectorBase<Scalar> >
    p = inArgs.get_p(0).assert_not_null();
    
  if (integrateStateAndSens) {
    if (trace)
      *out << "\nComputing D_d_hat_d_p by integrating state+sens ...\n";
    stateAndSensInitCond_.set_p(p_index_,p);
  }
  else {
    if (trace)
      *out << "\nComputing just d_hat by integrating the state ...\n";
    stateInitCond_.set_p(p_index_,p);
  }
    
  //
  // C) Setup the stepper and the integrator for state+sens or only state
  //

  RCP<IntegratorBase<Scalar> > integrator;
  if (integrateStateAndSens) {
    stateAndSensStepper_->setInitialCondition(stateAndSensInitCond_);
    stateAndSensIntegrator_->setStepper(stateAndSensStepper_,finalTime_);
    integrator = stateAndSensIntegrator_;
  }
  else {
    stateStepper_->setInitialCondition(stateInitCond_);
    stateIntegrator_->setStepper(stateStepper_,finalTime_);
    integrator = stateIntegrator_;
  }
    
  //
  // D) Setup for computing and processing the individual response functions
  //

  ForwardResponseSensitivityComputer<Scalar>
    forwardResponseSensitivityComputer;
  forwardResponseSensitivityComputer.setOStream(out);
  forwardResponseSensitivityComputer.setVerbLevel(localVerbLevel);
  forwardResponseSensitivityComputer.dumpSensitivities(dumpSensitivities_);
  forwardResponseSensitivityComputer.setResponseFunction(
    responseFuncs_[0],
    responseFuncs_[0]->createInArgs(), // Will replace this for each point!
    p_index_, g_index_
    );

  // D.a) Create storage for the individual response function ouptuts g_k
  // and its derivaitves

  RCP<Thyra::VectorBase<Scalar> > g_k;
  RCP<Thyra::MultiVectorBase<Scalar> > D_g_hat_D_p_k;

  if (!is_null(d_hat)) {
    g_k = forwardResponseSensitivityComputer.create_g_hat();
  }
  if (!D_d_hat_D_p.isEmpty()) {
    D_g_hat_D_p_k = forwardResponseSensitivityComputer.create_D_g_hat_D_p();
  }

  // D.b) Zero out d_hat and D_d_hat_D_p if we are doing a summation type of
  // evaluation
  if (responseType_ == FSIAMET::RESPONSE_TYPE_SUM) {
    if (!is_null(d_hat)) {
      assign( d_hat.ptr(), ST::zero() );
    }
    if (!D_d_hat_D_p.isEmpty()) {
      assign( D_d_hat_D_p.getMultiVector().ptr(), ST::zero() );
    }
  }

  // D.c) Get product vector and multi-vector interfaces if
  // we are just blocking up response functions
  RCP<Thyra::ProductVectorBase<Scalar> > prod_d_hat;
  RCP<Thyra::ProductMultiVectorBase<Scalar> > prod_D_d_hat_D_p_trans;
  if (responseType_ == FSIAMET::RESPONSE_TYPE_BLOCK) {
    prod_d_hat = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
      d_hat, true);
    prod_D_d_hat_D_p_trans = rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(
      D_d_hat_D_p.getMultiVector(), true);
  }
    
  //
  // E) Run the integrator at the response time points and assemble
  //    the response function and/or the response function
  //    derivative
  //

  if (trace) *out << "\nStarting integration and assembly loop ...\n";
    
  const int numResponseFunctions = responseTimes_.size();

  for (int k = 0; k < numResponseFunctions; ++k ) {
      
    OSTab tab(out);
      
    const Scalar t = responseTimes_[k];
      
    //
    // E.1) Integrate up to t and get x_bar and x_bar_dot
    //
    // Note, x_bar and x_bar_dot may be the state or the state+sens!
      
    if (trace)
      *out << "\nIntegrating state (or state+sens) for response k = "
           << k << " for t = " << t << " ...\n";
      
    RCP<const Thyra::VectorBase<Scalar> > x_bar, x_bar_dot;

    {
#ifdef ENABLE_RYTHMOS_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR(
        "Rythmos:ForwardSensitivityIntegratorAsModelEvaluator::evalModel: integrate"
        );
#endif
      get_fwd_x_and_x_dot( *integrator, t, &x_bar, &x_bar_dot );
    }
      
    if (print_norms) {
      *out << "\n||x_bar||inf = " << norms_inf(*x_bar) << endl;
      *out << "\n||x_bar_dot||inf = " << norms_inf(*x_bar_dot) << endl;
    }
    if (print_x) {
      *out << "\nx_bar = " << *x_bar << endl;
      *out << "\nx_bar_dot = " << *x_bar_dot << endl;
    }
      
    // Extra underlying state and sensitivities
    RCP<const Thyra::VectorBase<Scalar> > x, x_dot;
    RCP<const Thyra::MultiVectorBase<Scalar> > S, S_dot;
    if (integrateStateAndSens) {
      extractStateAndSens( x_bar, x_bar_dot, &x, &S, &x_dot, &S_dot );
    }
    else {
      x = x_bar;
      x_dot = x_bar_dot;
    }
      
    //
    // E.2) Evaluate the response function g_k and its derivatives at this
    // time point
    //
      
    if (trace)
      *out << "\nEvaluating response function k = " << k << " at time t = " << t << " ...\n";
      
    RCP<const Thyra::ModelEvaluator<Scalar> > responseFunc = responseFuncs_[k];
      
    MEB::InArgs<Scalar> responseInArgs = responseFunc->createInArgs();
    responseInArgs.setArgs(responseFuncInArgs_[k]);
    responseInArgs.set_p(p_index_,p);

    forwardResponseSensitivityComputer.resetResponseFunction(
      responseFunc, responseInArgs
      );

    if (!is_null(D_g_hat_D_p_k)) {
      forwardResponseSensitivityComputer.computeResponseAndSensitivity(
        x_dot.get(), S_dot.get(), *x, *S, t, g_k.get(), &*D_g_hat_D_p_k
        );
    }
    else {
      forwardResponseSensitivityComputer.computeResponse(
        x_dot.get(), *x, t, g_k.get()
        );
    }

    //
    // E.3) Assemble the final response functions and derivatives
    //

    // E.3.a) Assemble d_hat from g_k
    if (!is_null(d_hat)) {
      if (responseType_ == FSIAMET::RESPONSE_TYPE_SUM) {
        if (trace) *out << "\nd_hat += g_k ...\n";
        Vp_V( d_hat.ptr(), *g_k );
      }
      else if (responseType_ == FSIAMET::RESPONSE_TYPE_BLOCK) {
        if (trace) *out << "\nd_hat["<<k<<"] = g_k ...\n";
        assign( prod_d_hat->getNonconstVectorBlock(k).ptr(), *g_k );
      }
    }

    // E.3.b) Assemble D_d_hat_Dp from D_g_hat_D_p_k
    if (!D_d_hat_D_p.isEmpty()) {
      if (responseType_ == FSIAMET::RESPONSE_TYPE_SUM) {
        if (trace) *out << "\nD_d_hat_D_p += D_g_hat_D_p_k ...\n";
        Vp_V( D_d_hat_D_p.getMultiVector().ptr(), *D_g_hat_D_p_k );
        if (dumpSensitivities_ || print_x)
          *out << "\nD_d_hat_D_p = "
               << Teuchos::describe(*D_d_hat_D_p.getMultiVector(),Teuchos::VERB_EXTREME);
      }
      else if (responseType_ == FSIAMET::RESPONSE_TYPE_BLOCK) {
        if (trace) *out << "\nD_d_hat_D_p["<<k<<"] = D_g_hat_D_p_k ...\n";
        assign(
          prod_D_d_hat_D_p_trans->getNonconstMultiVectorBlock(k).ptr(),
          *D_g_hat_D_p_k
          );
      }
    }

  }
    
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_SENSITIVITY_INTEGRATOR_AS_MODEL_EVALUATOR_HPP
