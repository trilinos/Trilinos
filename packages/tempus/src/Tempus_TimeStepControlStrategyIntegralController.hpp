//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeStepControlStrategy_IntegralController_hpp
#define Tempus_TimeStepControlStrategy_IntegralController_hpp

#include "Tempus_config.hpp"
#include "Tempus_NumericalUtils.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
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
 *      (\Delta t)_n \left( \epsilon_n ^{-k_I / p}
 *           \epsilon_{n-1}^{k_P / p} \epsilon_{n-2}^{-k_D / p} \right)
 * \f]
 *
 * PI = Proportional-Integral Controller
 * \f[
 *      (\Delta t)_{n+1} =
 *      (\Delta t)_n \left( \epsilon_n ^{-k_I / p}
 *           \epsilon_{n-1}^{k_P / p} \right)
 * \f]
 *
 * I = Integral Controller
 * \f[
 *      (\Delta t)_{n+1} =
 *      (\Delta t)_n \left( \epsilon_n ^{-k_I / p} \right)
 * \f]
 *
 * where \f$\epsilon_n \f$ is the error at time step \f$n\f$
 * and \f$p\f$ is the order of the embedded solution, which
 * is assumed to be the low order solution (i.e., the time
 * step order minus one).
 *
 * Appropriate for Explicit Methods
 */
template <class Scalar>
class TimeStepControlStrategyIntegralController
  : virtual public TimeStepControlStrategy<Scalar> {
 public:
  /// Default Constructor
  TimeStepControlStrategyIntegralController()
    : controller_("PID"),
      KI_(0.58),
      KP_(0.21),
      KD_(0.10),
      safetyFactor_(0.90),
      safetyFactorAfterReject_(0.9),
      facMax_(5.0),
      facMin_(0.5)
  {
    facMaxINPUT_ = facMax_;

    this->setStrategyType("Integral Controller");
    this->setStepType("Variable");
    this->setName("Integral Controller");
    this->initialize();
  }

  /// Full Constructor
  TimeStepControlStrategyIntegralController(
      std::string controller, Scalar KI, Scalar KP, Scalar KD,
      Scalar safetyFactor, Scalar safetyFactorAfterReject, Scalar facMax,
      Scalar facMin, std::string name = "Integral Controller")
    : controller_(controller),
      KI_(KI),
      KP_(KP),
      KD_(KD),
      safetyFactor_(safetyFactor),
      safetyFactorAfterReject_(safetyFactorAfterReject),
      facMax_(facMax),
      facMin_(facMin)
  {
    facMaxINPUT_ = facMax_;

    this->setStrategyType("Integral Controller");
    this->setStepType("Variable");
    this->setName(name);
    this->initialize();
  }

  /// Destructor
  virtual ~TimeStepControlStrategyIntegralController() {}

  /** \brief Set the time step size.*/
  virtual void setNextTimeStep(
      const TimeStepControl<Scalar> &tsc,
      Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
      Status & /* integratorStatus */) override
  {
    using Teuchos::RCP;

    this->checkInitialized();

    // Take first step with initial time step provided
    if (!firstSuccessfulStep_) {
      firstSuccessfulStep_ = true;
      return;
    }

    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    Scalar beta = 1.0;

    // assumes the embedded solution is the low order solution
    int order = workingState->getOrder() - 1;
    Scalar dt = workingState->getTimeStep();

    // Get the relative errors.
    Scalar errN   = workingState->getErrorRel();
    Scalar errNm1 = workingState->getErrorRelNm1();
    Scalar errNm2 = workingState->getErrorRelNm2();

    if (errN < numericalTol<Scalar>()) errN = 1.0;
    if (errNm1 < numericalTol<Scalar>()) errNm1 = 1.0;
    if (errNm2 < numericalTol<Scalar>()) errNm2 = 1.0;

    Scalar k1 = Teuchos::as<Scalar>(-KI_ / order);
    Scalar k2 = Teuchos::as<Scalar>(KP_ / order);
    Scalar k3 = Teuchos::as<Scalar>(-KD_ / order);

    k1 = std::pow(errN, k1);
    k2 = std::pow(errNm1, k2);
    k3 = std::pow(errNm2, k3);

    if (controller_ == "I")
      beta = safetyFactor_ * k1;
    else if (controller_ == "PI")
      beta = safetyFactor_ * k1 * k2;
    else  // (controller_ == "PID")
      beta = safetyFactor_ * k1 * k2 * k3;

    beta = std::max(facMin_, beta);
    beta = std::min(facMax_, beta);

    // new (optimal) suggested time step
    dt = beta * dt;

    if (workingState->getSolutionStatus() == Status::PASSED ||
        workingState->getSolutionStatus() == Status::WORKING) {
      if (lastStepRejected_) {
        dt = std::min(dt, workingState->getTimeStep());
      }
      else {
        facMax_ = facMaxINPUT_;
      }
      lastStepRejected_ = false;
    }
    else {
      facMax_           = safetyFactorAfterReject_;
      lastStepRejected_ = true;
    }

    // update dt
    workingState->setTimeStep(dt);
    workingState->setTime(solutionHistory->getCurrentState()->getTime() + dt);
  }

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override
  {
    return "Tempus::TimeStepControlStrategyIntegralController";
  }

  void describe(Teuchos::FancyOStream &out,
                const Teuchos::EVerbosityLevel verbLevel) const override
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, this->description());
    l_out->setOutputToRootOnly(0);

    *l_out << "\n--- " << this->description() << " ---" << std::endl;

    if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
      *l_out << "  Strategy Type                      = "
             << this->getStrategyType() << std::endl
             << "  Step Type                          = " << this->getStepType()
             << std::endl
             << "  Controller Type                    = " << getController()
             << std::endl
             << "  KI                                 = " << getKI()
             << std::endl
             << "  KP                                 = " << getKP()
             << std::endl
             << "  KD                                 = " << getKD()
             << std::endl
             << "  Safety Factor                      = " << getSafetyFactor()
             << std::endl
             << "  Safety Factor After Step Rejection = "
             << getSafetyFactorAfterReject() << std::endl
             << "  Maximum Safety Factor (INPUT)      = " << facMaxINPUT_
             << std::endl
             << "  Maximum Safety Factor              = " << getFacMax()
             << std::endl
             << "  Minimum Safety Factor              = " << getFacMin()
             << std::endl;
      *l_out << std::string(this->description().length() + 8, '-') << std::endl;
    }
  }
  //@}

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::parameterList("Time Step Control Strategy");

    pl->set<std::string>("Strategy Type", this->getStrategyType(),
                         "Integral Controller");
    pl->set<std::string>("Controller Type", getController(),
                         "Proportional-Integral-Derivative");
    pl->set<Scalar>("KI", getKI(), "Integral gain");
    pl->set<Scalar>("KP", getKP(), "Proportional gain");
    pl->set<Scalar>("KD", getKD(), "Derivative gain");
    pl->set<Scalar>("Safety Factor", getSafetyFactor(), "Safety Factor");
    pl->set<Scalar>("Safety Factor After Step Rejection",
                    getSafetyFactorAfterReject(),
                    "Safety Factor Following Step Rejection");
    pl->set<Scalar>("Maximum Safety Factor", getFacMax(),
                    "Maximum Safety Factor");
    pl->set<Scalar>("Minimum Safety Factor", getFacMin(),
                    "Minimum Safety Factor");
    return pl;
  }

  virtual void initialize() const override
  {
    TEUCHOS_TEST_FOR_EXCEPTION(safetyFactor_ <= 0.0, std::out_of_range,
                               "Error - Invalid value of Safety Factory= "
                                   << safetyFactor_ << "!  \n"
                                   << "Safety Factor must be > 0.0.\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        facMax_ <= 0.0, std::out_of_range,
        "Error - Invalid value of Maximum Safety Factory= "
            << facMax_ << "!  \n"
            << "Maximum Safety Factor must be > 0.0.\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        facMax_ <= 0.0, std::out_of_range,
        "Error - Invalid value of Minimum Safety Factory= "
            << facMin_ << "!  \n"
            << "Minimum Safety Factor must be > 0.0.\n");

    TEUCHOS_TEST_FOR_EXCEPTION(((controller_ != "I") && (controller_ != "PI") &&
                                (controller_ != "PID")),
                               std::invalid_argument,
                               "Error - Invalid choice of Controller Type = "
                                   << controller_ << "!  \n"
                                   << "Valid Choice are ['I', 'PI', 'PID'].\n");

    this->isInitialized_ = true;  // Only place where this is set to true!
  }

  virtual std::string getController() const { return controller_; }
  virtual Scalar getKI() const { return KI_; }
  virtual Scalar getKP() const { return KP_; }
  virtual Scalar getKD() const { return KD_; }
  virtual Scalar getSafetyFactor() const { return safetyFactor_; }
  virtual Scalar getSafetyFactorAfterReject() const
  {
    return safetyFactorAfterReject_;
  }
  virtual Scalar getFacMax() const { return facMax_; }
  virtual Scalar getFacMin() const { return facMin_; }

  virtual void setController(std::string c)
  {
    controller_          = c;
    this->isInitialized_ = false;
  }
  virtual void setKI(Scalar k)
  {
    KI_                  = k;
    this->isInitialized_ = false;
  }
  virtual void setKP(Scalar k)
  {
    KP_                  = k;
    this->isInitialized_ = false;
  }
  virtual void setKD(Scalar k)
  {
    KD_                  = k;
    this->isInitialized_ = false;
  }
  virtual void setSafetyFactor(Scalar f)
  {
    safetyFactor_        = f;
    this->isInitialized_ = false;
  }
  virtual void setSafetyFactorAfterReject(Scalar f)
  {
    safetyFactorAfterReject_ = f;
    this->isInitialized_     = false;
  }
  virtual void setFacMax(Scalar f)
  {
    facMax_              = f;
    facMaxINPUT_         = f;
    this->isInitialized_ = false;
  }
  virtual void setFacMin(Scalar f)
  {
    facMin_              = f;
    this->isInitialized_ = false;
  }

 private:
  std::string controller_;          ///< Control type ['I', 'PI', 'PID']
  Scalar KI_;                       ///< Integral gain
  Scalar KP_;                       ///< Proportional gain
  Scalar KD_;                       ///< Derivative gain
  Scalar safetyFactor_;             ///< Safety Factor
  Scalar safetyFactorAfterReject_;  ///< Safety Factor Following Step Rejection
  Scalar facMaxINPUT_;              ///< Maximum Safety Factor from input
  Scalar facMax_;                   ///< Maximum Safety Factor
  Scalar facMin_;                   ///< Minimum Safety Factor
  bool firstSuccessfulStep_ = false;
  bool lastStepRejected_    = false;
};

// Nonmember constructor.
template <class Scalar>
Teuchos::RCP<TimeStepControlStrategyIntegralController<Scalar> >
createTimeStepControlStrategyIntegralController(
    const Teuchos::RCP<Teuchos::ParameterList> pList,
    std::string name = "Integral Controller")
{
  using Teuchos::rcp;
  auto tscs = rcp(new TimeStepControlStrategyIntegralController<Scalar>());
  if (pList == Teuchos::null || pList->numParams() == 0) return tscs;

  TEUCHOS_TEST_FOR_EXCEPTION(
      pList->get<std::string>("Strategy Type") != "Integral Controller",
      std::logic_error,
      "Error - Strategy Type != 'Integral Controller'.  (='" +
          pList->get<std::string>("Strategy Type") + "')\n");

  pList->validateParametersAndSetDefaults(*tscs->getValidParameters());

  tscs->setController(pList->get<std::string>("Controller Type"));
  tscs->setKI(pList->get<Scalar>("KI"));
  tscs->setKP(pList->get<Scalar>("KP"));
  tscs->setKD(pList->get<Scalar>("KD"));
  tscs->setSafetyFactor(pList->get<Scalar>("Safety Factor"));
  tscs->setSafetyFactorAfterReject(
      pList->get<Scalar>("Safety Factor After Step Rejection"));
  tscs->setFacMax(pList->get<Scalar>("Maximum Safety Factor"));
  tscs->setFacMin(pList->get<Scalar>("Minimum Safety Factor"));

  tscs->setName(name);
  tscs->initialize();

  return tscs;
}

/// Nonmember function to return ParameterList with default values.
template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
getTimeStepControlStrategyIntegralControllerPL()
{
  auto t = rcp(new Tempus::TimeStepControlStrategyIntegralController<Scalar>());
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      t->getValidParameters());
}

}  // namespace Tempus
#endif  // Tempus_TimeStepControlStrategy_IntegralController_hpp
