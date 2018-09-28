// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategy_IntegralController_hpp
#define Tempus_TimeStepControlStrategy_IntegralController_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"


namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 *
 * Gustaf Soderlind.
 * Automatic control and adaptive time-stepping.
 * Numerical Algorithms, 31(1):281–310, Dec 2002.
 *
 * The step size is chosen based on "Controller Type":
 *
 * PID = Proportional-Integral-Derivative Controller
 * \f[
 *      (\Delta t)_{n+1} =
 *      (\Delta t)_n \left( \epsilon_n ^{-k_1 / p} \epsilon_{n-1}^{k_2 / p} \epsilon_{n-2}^{-k_3 / p} \right)
 * \f]
 *
 * PI = Proportional-Integral Controller
 * \f[
 *      (\Delta t)_{n+1} =
 *      (\Delta t)_n \left( \epsilon_n ^{-k_1 / p} \epsilon_{n-1}^{k_2 / p} \right)
 * \f]
 *
 * I = Integral Controller
 * \f[
 *      (\Delta t)_{n+1} =
 *      (\Delta t)_n \left( \epsilon_n ^{-k_1 / p} \right)
 * \f]
 *
 * where \f$\epsilon_n \f$ is the error at time step \f$n\f$.
 *
 * Appropriate for Explicit Methods
 */
template<class Scalar>
class TimeStepControlStrategyIntegralController
  : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Constructor
  TimeStepControlStrategyIntegralController(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null){
     this->setParameterList(pList);
  }

  /// Destructor
  virtual ~TimeStepControlStrategyIntegralController(){}

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
     const Scalar errorRel = metaData->getErrorRel();
     Scalar beta = 1.0;

     // assumes the embedded solution is the low order solution
     int order = metaData->getOrder() - 1;
     Scalar dt = metaData->getDt();
     //bool printChanges = solutionHistory->getVerbLevel() !=
        //Teuchos::as<int>(Teuchos::VERB_NONE);

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

     if (controller_ == "I")
        beta = safetyFactor_*k1;
     else if (controller_ == "PI")
        beta = safetyFactor_ *k1*k2;
     else // (controller_ == "PID")
        beta = safetyFactor_*k1*k2*k3;

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
        facMax_ = tscsPL_->get<Scalar>("Safety Factor After Step Rejection");
        lastStepRejected_ = true;
     }

     // update dt
     metaData->setDt(dt);
  }

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
  void setParameterList(
    const Teuchos::RCP<Teuchos::ParameterList> & pList) override
  {

     if (pList == Teuchos::null) {
        // Create default parameters if null, otherwise keep current parameters.
        if (tscsPL_ == Teuchos::null) {
           tscsPL_ = Teuchos::parameterList("TimeStepControlStrategyIntegralController");
           *tscsPL_= *(this->getValidParameters());
        }
     } else {
        tscsPL_ = pList;
     }
     tscsPL_->validateParametersAndSetDefaults(*this->getValidParameters());

     errN_ = Scalar(1.0);
     errNm1_ = Scalar(1.0);
     errNm2_ = Scalar(1.0);
     k1_ = tscsPL_->get<Scalar>("KI");
     k2_ = tscsPL_->get<Scalar>("KP");
     k3_ = tscsPL_->get<Scalar>("KD");
     safetyFactor_ = tscsPL_->get<Scalar>("Safety Factor");
     facMax_ = tscsPL_->get<Scalar>("Maximum Safety Factor");
     facMin_ = tscsPL_->get<Scalar>("Minimum Safety Factor");
     controller_ = tscsPL_->get<std::string>("Controller Type");

     TEUCHOS_TEST_FOR_EXCEPTION(safetyFactor_ <= 0.0, std::out_of_range,
     "Error - Invalid value of Safety Factory= " << safetyFactor_ << "!  \n"
     << "Safety Factor must be > 0.0.\n");

     TEUCHOS_TEST_FOR_EXCEPTION(facMax_ <= 0.0, std::out_of_range,
     "Error - Invalid value of Maximum Safety Factory= " << facMax_ << "!  \n"
     << "Maximum Safety Factor must be > 0.0.\n");

     TEUCHOS_TEST_FOR_EXCEPTION(facMax_<= 0.0, std::out_of_range,
     "Error - Invalid value of Minimum Safety Factory= " << facMin_ << "!  \n"
     << "Minimum Safety Factor must be > 0.0.\n");

     TEUCHOS_TEST_FOR_EXCEPTION(((controller_ != "I") and
                                 (controller_ != "PI") and
                                 (controller_ != "PID")), std::invalid_argument,
     "Error - Invalid choice of Controller Type = " << controller_ << "!  \n"
     << "Valid Choice are ['I', 'PI', 'PID'].\n");

  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override
  {
     Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

     pl->set<std::string>("Name","Integral Controller","Integral Controller");
     pl->set<std::string>("Controller Type","PID",
                          "Proportional-Integral-Derivative");
     pl->set<Scalar>("KI" , 0.58, "Integral gain");
     pl->set<Scalar>("KP" , 0.21, "Proportional gain");
     pl->set<Scalar>("KD" , 0.10, "Derivative gain");
     pl->set<Scalar>("Safety Factor" , 0.90, "Safety Factor");
     pl->set<Scalar>("Maximum Safety Factor" , 5.0, "Maximum Safety Factor");
     pl->set<Scalar>("Minimum Safety Factor" , 0.5, "Minimum Safety Factor");
     pl->set<Scalar>("Safety Factor After Step Rejection" , 0.9,
                     "Safety Factor Following Step Rejection");
     return pl;
  }

  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override {
     return tscsPL_;
  }

  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() override {
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
    std::string controller_;

};
} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_IntegralController_hpp
