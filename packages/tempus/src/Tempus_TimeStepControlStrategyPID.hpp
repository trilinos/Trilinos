// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategy_PID_hpp
#define Tempus_TimeStepControlStrategy_PID_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"


namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 */
template<class Scalar>
class TimeStepControlStrategyPID
  : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Constructor
  TimeStepControlStrategyPID(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null){
     this->setParameterList(pList);
  }

  /// Destructor
  virtual ~TimeStepControlStrategyPID(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(const TimeStepControl<Scalar> tsc,
    Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
    Status & integratorStatus) override
  {

     // Take first step with initial time step provided
     if (!firstSuccessfulStep_){
        firstSuccessfulStep_ = true;
        return;
     }

     Teuchos::RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
     Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
     const Scalar errorAbs = metaData->getErrorAbs();
     volatile const int iStep = metaData->getIStep();
     volatile const Scalar errorRel = metaData->getErrorRel();
     int order = metaData->getOrder();
     Scalar dt = metaData->getDt();
     bool printChanges = solutionHistory->getVerbLevel() !=
        Teuchos::as<int>(Teuchos::VERB_NONE);

     Teuchos::RCP<Teuchos::FancyOStream> out = tsc.getOStream();
     Teuchos::OSTab ostab(out,1,"getNextTimeStep");

     // For now, only time step is changed (not order)
     /*
     auto changeOrder = [] (int order_old, int order_new, std::string reason) {
        std::stringstream message;
        message << "     (order = " << std::setw(2) << order_old
           <<       ", new = " << std::setw(2) << order_new
           << ")  " << reason << std::endl;
        return message.str();
     };
     */

     // update errors
     errNm2_ = errNm1_;
     errNm1_ = errN_;
     errN_ = errorRel;

     Scalar k1 = Teuchos::as<Scalar>(-k1_ / order);
     Scalar k2 = Teuchos::as<Scalar>(k2_ / order);
     Scalar k3 = Teuchos::as<Scalar>(-k3_ / order);

     k1 = std::pow(errN_, k1);
     k2 = std::pow(errNm1_, k2);
     k3 = std::pow(errNm2_, k3);
     Scalar beta = safetyFactor_*k1*k2*k3;
     beta = std::max(facMin_, beta);
     beta = std::min(facMax_, beta);

     // new (optimal) suggested time step
     dt = beta * dt;

     if (workingState->getSolutionStatus() == Status::PASSED or
         workingState->getSolutionStatus() == Status::WORKING) {
        if(lastStepRejected_){
           dt = std::min(dt, metaData->getDt());
        } else {
           facMax_ = tscsPL_->get<Scalar>("Maximum Safety Factor");
        }
        lastStepRejected_ = false;
     } else {
        facMax_ = 0.90;
        lastStepRejected_ = true;
     }

     // update dt
     metaData->setDt(dt);
  }

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pList){

     if (pList == Teuchos::null) {
        // Create default parameters if null, otherwise keep current parameters.
        if (tscsPL_ == Teuchos::null) {
           tscsPL_ = Teuchos::parameterList("TimeStepControlStrategyPID");
           *tscsPL_= *(this->getValidParameters());
        }
     } else {
        tscsPL_ = pList;
     }
     tscsPL_->validateParametersAndSetDefaults(*this->getValidParameters());

     errN_ = Scalar(1.0);
     errNm1_ = Scalar(1.0);
     errNm2_ = Scalar(1.0);
     k1_ = tscsPL_->get<Scalar>("K1");
     k2_ = tscsPL_->get<Scalar>("K2");
     k3_ = tscsPL_->get<Scalar>("K3");
     safetyFactor_ = tscsPL_->get<Scalar>("Safety Factor");
     facMax_ = tscsPL_->get<Scalar>("Maximum Safety Factor");
     facMin_ = tscsPL_->get<Scalar>("Minimum Safety Factor");

     TEUCHOS_TEST_FOR_EXCEPTION(safetyFactor_ <= 0.0, std::out_of_range,
     "Error - Invalid value of Safety Factory= " << safetyFactor_ << "!  \n"
     << "Safety Factor must be > 0.0.\n");

     TEUCHOS_TEST_FOR_EXCEPTION(facMax_ <= 0.0, std::out_of_range,
     "Error - Invalid value of Maximum Safety Factory= " << facMax_ << "!  \n"
     << "Maximum Safety Factor must be > 0.0.\n");

     TEUCHOS_TEST_FOR_EXCEPTION(facMax_<= 0.0, std::out_of_range,
     "Error - Invalid value of Minimum Safety Factory= " << facMin_ << "!  \n"
     << "Minimum Safety Factor must be > 0.0.\n");

  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
     Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

     pl->set<std::string>("Name","PID");
     pl->set<Scalar>("K1" , 0.58, "");
     pl->set<Scalar>("K2" , 0.21, "");
     pl->set<Scalar>("K3" , 0.10, "");
     pl->set<Scalar>("Safety Factor" , 0.90, "Safety Factor");
     pl->set<Scalar>("Maximum Safety Factor" , 5.0, "Maximum Safety Factor");
     pl->set<Scalar>("Minimum Safety Factor" , 0.5, "Minimum Safety Factor");
     return pl;
  }

  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() {
     return tscsPL_;
  }

  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() {
     Teuchos::RCP<Teuchos::ParameterList> temp_plist = tscsPL_;
     tscsPL_ = Teuchos::null;
     return(temp_plist);
  }
  //@}

private:
    Teuchos::RCP<Teuchos::ParameterList> tscsPL_;
    Scalar k1_;
    Scalar k2_;
    Scalar k3_;
    Scalar errN_;
    Scalar errNm1_;
    Scalar errNm2_;
    Scalar safetyFactor_;
    Scalar facMax_;
    Scalar facMin_;
    bool firstSuccessfulStep_ = false;
    bool lastStepRejected_ = false;

};
} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_PID_hpp
