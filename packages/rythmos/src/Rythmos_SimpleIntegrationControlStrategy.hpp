
#ifndef RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_HPP
#define RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_HPP


#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


/** \brief Base class for strategy objects that control integration by
 * selecting step sizes for a stepper.
 *
 * ToDo: Finish Implementation!
 */
template<class Scalar>
class SimpleIntegrationControlStrategy
  : virtual public IntegrationControlStrategyBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief Constructors/Initializers. */
  //@{

  /** \brief . */
  SimpleIntegrationControlStrategy();

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \brief Overridden from IntegrationControlStrategyBase */
  //@{

  /** \brief . */
  RCP<IntegrationControlStrategyBase<Scalar> >
  cloneIntegrationControlStrategy() const;

  /** \brief . */
  void resetIntegrationControlStrategy(
    const TimeRange<Scalar> &integrationTimeDomain
    );

  /** \brief . */
  StepControlInfo<Scalar>
  getNextStepControlInfo(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfoLast,
    const int timeStepIter
    );

  //@}

private:

  bool takeVariableSteps_;
  Scalar max_dt_;
  int numTimeSteps_;
  Scalar fixed_dt_;

  TimeRange<Scalar> integrationTimeDomain_;

  static const std::string takeVariableSteps_name_;
  static const bool takeVariableSteps_default_;

  static const std::string max_dt_name_;
  static const double max_dt_default_;

  static const std::string numTimeSteps_name_;
  static const int numTimeSteps_default_;

  static const std::string fixed_dt_name_;
  static const double fixed_dt_default_;

};


/** \brief .
 *
 * \relates SimpleIntegrationControlStrategy
 */
template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy()
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrationControl = Teuchos::rcp(new SimpleIntegrationControlStrategy<Scalar>());
  return integrationControl;
}


/** \brief .
 *
 * \relates SimpleIntegrationControlStrategy
 */
template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy( const RCP<ParameterList> &paramList )
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrationControl = Teuchos::rcp(new SimpleIntegrationControlStrategy<Scalar>());
  integrationControl->setParameterList(paramList);
  return integrationControl;
}


//
// Implementation
//


// Static members


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::takeVariableSteps_name_
= "Take Variable Steps";

template<class Scalar> 
const bool
SimpleIntegrationControlStrategy<Scalar>::takeVariableSteps_default_
= true;


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::max_dt_name_
= "Max dt";

template<class Scalar> 
const double
SimpleIntegrationControlStrategy<Scalar>::max_dt_default_
= std::numeric_limits<Scalar>::max();


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::numTimeSteps_name_
= "Number of Time Steps";

template<class Scalar> 
const int
SimpleIntegrationControlStrategy<Scalar>::numTimeSteps_default_
= -1;


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::fixed_dt_name_
= "Fixed dt";

template<class Scalar> 
const double
SimpleIntegrationControlStrategy<Scalar>::fixed_dt_default_
= -1.0;


// Constructors/Initializers


template<class Scalar>
SimpleIntegrationControlStrategy<Scalar>::SimpleIntegrationControlStrategy()
  :takeVariableSteps_(takeVariableSteps_default_),
   max_dt_(max_dt_default_),
   numTimeSteps_(numTimeSteps_default_),
   fixed_dt_(fixed_dt_default_)
{}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void SimpleIntegrationControlStrategy<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  using Teuchos::as;
  using Teuchos::get;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);
  takeVariableSteps_ = paramList->get(
    takeVariableSteps_name_, takeVariableSteps_ );
  if (!takeVariableSteps_) {
    numTimeSteps_ = paramList->get(numTimeSteps_name_,numTimeSteps_);
    fixed_dt_ = paramList->get(fixed_dt_name_,fixed_dt_);
    TEST_FOR_EXCEPTION(
      numTimeSteps_ < 0 && fixed_dt_ <= ST::zero(), std::logic_error,
      "Error, when taking fixed steps, the user must set the parameters "
      "\""<<numTimeSteps_name_<<"\" > 0 or \""<<fixed_dt_name_<<"\" > 0.0!" );
  }
  else {
    max_dt_ = paramList->get(max_dt_name_,max_dt_);
  }
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const ParameterList>
SimpleIntegrationControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(
      takeVariableSteps_name_, takeVariableSteps_default_,
      "Take variable time steps or fixed time steps.\n"
      "If set to false, then the parameter \"" + fixed_dt_name_ + "\"\n"
      "or \"" + numTimeSteps_name_ + "\" must be set!"
      );
    pl->set(
      max_dt_name_, max_dt_default_,
      "Gives the max size of the variable time steps.  This is only read and used if\n"
      "\"" + takeVariableSteps_name_ + "\" is set to true."
      );
    pl->set(
      numTimeSteps_name_, numTimeSteps_default_,
      "Gives the number of fixed time steps.  The actual step size gets computed\n"
      "on the fly given the size of the time domain.\n"
      "This is only read and used if \"" + takeVariableSteps_name_ + "\" is set to false\n"
      "and \"" + fixed_dt_name_ + "\" is set to < 0.0."
      );
    pl->set(
      fixed_dt_name_, fixed_dt_default_,
      "Gives the size of the fixed time steps.  This is only read and used if\n"
      "\"" + takeVariableSteps_name_ + "\" is set to false."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from IntegrationControlStrategyBase


template<class Scalar>
RCP<IntegrationControlStrategyBase<Scalar> >
SimpleIntegrationControlStrategy<Scalar>::cloneIntegrationControlStrategy() const
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrCtrlStry = simpleIntegrationControlStrategy<Scalar>();
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList))
    integrCtrlStry->setParameterList(Teuchos::parameterList(*paramList));
  integrCtrlStry->takeVariableSteps_ = takeVariableSteps_;
  integrCtrlStry->max_dt_ = max_dt_;
  integrCtrlStry->numTimeSteps_ = numTimeSteps_;
  integrCtrlStry->fixed_dt_ = fixed_dt_;
  return integrCtrlStry;
}


template<class Scalar>
void
SimpleIntegrationControlStrategy<Scalar>::resetIntegrationControlStrategy(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain.length() > ST::zero());
#endif
  integrationTimeDomain_ = integrationTimeDomain;
  if (takeVariableSteps_) {
    if (max_dt_ < ST::zero()) {
      max_dt_ = integrationTimeDomain_.length();
    }
  }
  else {
    if (fixed_dt_ < ST::zero()) {
#ifdef HAVE_RYTHMOS_DEBUG
      TEUCHOS_ASSERT(numTimeSteps_ > 0);
#endif
      fixed_dt_ = integrationTimeDomain_.length()/numTimeSteps_;
    }
  }
}


template<class Scalar>
StepControlInfo<Scalar>
SimpleIntegrationControlStrategy<Scalar>::getNextStepControlInfo(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfoLast,
  const int timeStepIter
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain_.length() > ST::zero());
#endif
  
  StepControlInfo<Scalar> trialStepCtrlInfo;
  
  if (takeVariableSteps_) {
    trialStepCtrlInfo.stepType = STEP_TYPE_VARIABLE;
    trialStepCtrlInfo.stepSize = max_dt_;
  }
  else {
    trialStepCtrlInfo.stepType = STEP_TYPE_FIXED;
    trialStepCtrlInfo.stepSize = fixed_dt_;
  }
  
  return trialStepCtrlInfo;
  
}


} // namespace Rythmos


#endif // RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_HPP
