//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeStepControlStrategy_BasicVS_hpp
#define Tempus_TimeStepControlStrategy_BasicVS_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"

namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 *  This TimeStepControlStrategy primarily tries to maintain a
 *  certain level of change in the solution ill-respective of the
 *  error involved, e.g., the solution should change between 1% and
 *  3% (\f$\eta_{min}=0.01\f$ and \f$\eta_{max}=0.03\f$) every
 *  time step.  The relative solution change is measured by
 *  \f[
 *    \eta_{n-1} = \frac{|| x_{n-1} - x_{n-2} ||}{ || x_{n-2} || + \epsilon }
 *  \f]
 *  where \f$\epsilon\f$ is a small constant to ensure that \f$\eta_{n-1}\f$
 *  remains finite.  The user can select the desired relative
 *  change in the solution by choosing a range for \f$\eta_{n-1}\f$
 *  \f[
 *    \eta_{min} < \eta_{n-1} < \eta_{max}
 *  \f]
 *  If the solution change is outside this range, an amplification
 *  (\f$\rho\f$) or reduction factor (\f$\sigma\f$) is applied to
 *  the timestep to bring the solution change back into the desired
 *  range.  This can be written as
 *  \f[
 *    \Delta t_n = \left\{
 *      \begin{array}{rll}
 *        \sigma \Delta t_{n-1} & \mbox{if $\eta_{n-1} > \eta_{max}$}
 *                              & \mbox{where $0 < \sigma < 1$}             \\
 *        \rho   \Delta t_{n-1} & \mbox{else if $\eta_{n-1} < \eta_{min}$}
 *                              & \mbox{where $\rho > 1$}                   \\
 *               \Delta t_{n-1} &
 *                        \mbox{else if $\eta_{min}<\eta_{n-1}<\eta_{max}$} \\
 *      \end{array}
 *    \right.
 *  \f]
 *  In the full implementation, several other mechanisms can amplify
 *  or reduce the timestep.
 *  - Stepper fails
 *  - Maximum Absolute error, \f$e_{abs}^{max}\f$
 *  - Maximum Relative error, \f$e_{rel}^{max}\f$
 *  \f[
 *    \Delta t_n = \left\{
 *      \begin{array}{rll}
 *        \sigma \Delta t_{n-1} & \mbox{if Stepper fails}
 *                              & \mbox{where $0 < \sigma < 1$}            \\
 *        \rho   \Delta t_{n-1} & \mbox{else if $\eta_{n-1} < \eta_{min}$}
 *                              & \mbox{where $\rho > 1$}                  \\
 *        \sigma \Delta t_{n-1} & \mbox{else if $\eta_{n-1} > \eta_{max}$}
 *                              & \mbox{where $0 < \sigma < 1$}            \\
 *        \sigma \Delta t_{n-1} & \mbox{else if $e_{abs} > e_{abs}^{max}$}
 *                              & \mbox{where $0 < \sigma < 1$}            \\
 *        \sigma \Delta t_{n-1} & \mbox{else if $e_{rel} > e_{rel}^{max}$}
 *                              & \mbox{where $0 < \sigma < 1$}            \\
 *        \rho   \Delta t_{n-1} & \mbox{else if $p < p_{min}$}
 *                              & \mbox{where $\rho > 1$}                  \\
 *        \sigma \Delta t_{n-1} & \mbox{else if $p > p_{max}$}
 *                              & \mbox{where $0 < \sigma < 1$}            \\
 *               \Delta t_{n-1} & \mbox{else} &                            \\
 *      \end{array}
 *    \right.
 *  \f]
 *
 *  Note
 *  - Only one amplification or reduction is applied each timestep.
 *  - The priority is specified by the order of list.
 *  - The timestep, \f$\Delta t_n\f$, is still constrained to the
 *    maximum and minimum timestep size.
 *    \f$\Delta t_{min} < \Delta t_n < \Delta t_{max}\f$
 *  - If \f$ \eta_{min} < \eta_n < \eta_{max}\f$, the timestep
 *    is unchanged, i.e., constant timestep size.
 *  - To have constant timesteps, set \f$\eta_{min}=0\f$ and
 *    \f$\eta_{max}=10^{16}\f$.  These are the defaults.
 *  - From (Denner, 2014), amplification factor, \f$\rho\f$, is
 *    required to be less than 1.91 for stability (\f$\rho < 1.91\f$).
 *  - Denner (2014) suggests that \f$\eta_{min} = 0.1*\eta_{max}\f$
 *    and the numerical tests confirm this for their problems.
 *
 *  #### References
 *  Section 2.2.1 / Algorithm 2.4 of A. Denner, "Experiments on
 *  Temporal Variable Step BDF2 Algorithms", Masters Thesis,
 *  U Wisconsin-Madison, 2014.
 *
 */
template <class Scalar>
class TimeStepControlStrategyBasicVS
  : virtual public TimeStepControlStrategy<Scalar> {
 public:
  /// Default Constructor
  TimeStepControlStrategyBasicVS()
    : rho_(1.75), sigma_(0.5), minEta_(0.0), maxEta_(1.0e+16)
  {
    this->setStrategyType("Basic VS");
    this->setStepType("Variable");
    this->setName("Basic VS");
    this->initialize();
  }

  /// Full Constructor
  TimeStepControlStrategyBasicVS(Scalar rho, Scalar sigma, Scalar minEta,
                                 Scalar maxEta, std::string name = "Basic VS")
    : rho_(rho), sigma_(sigma), minEta_(minEta), maxEta_(maxEta)
  {
    this->setStrategyType("Basic VS");
    this->setStepType("Variable");
    this->setName(name);
    this->initialize();
  }

  /// Destructor
  virtual ~TimeStepControlStrategyBasicVS() {}

  /** \brief Set the time step size.*/
  virtual void setNextTimeStep(
      const TimeStepControl<Scalar> &tsc,
      Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
      Status & /* integratorStatus */) override
  {
    using Teuchos::RCP;

    this->checkInitialized();

    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    const Scalar errorAbs = workingState->getErrorAbs();
    const Scalar errorRel = workingState->getErrorRel();
    const int iStep       = workingState->getIndex();
    Scalar dt             = workingState->getTimeStep();

    Scalar rho   = getAmplFactor();
    Scalar sigma = getReductFactor();
    Scalar eta   = solutionHistory->getCurrentState()->getDxNormL2Rel();
    if (iStep == 1) eta = getMinEta();  // For first step use initial dt.

    // General rule: only increase/decrease dt once for any given reason.
    if (workingState->getSolutionStatus() == Status::FAILED) {
      tsc.printDtChanges(iStep, dt, dt * sigma,
                         "Stepper failure - Decreasing dt.");
      dt *= sigma;
    }
    else {                      // Stepper passed
      if (eta < getMinEta()) {  // increase dt
        tsc.printDtChanges(iStep, dt, dt * rho,
                           "Change too small (" + std::to_string(eta) + " < " +
                               std::to_string(getMinEta()) +
                               ").  Increasing dt.");
        dt *= rho;
      }
      else if (eta > getMaxEta()) {  // reduce dt
        tsc.printDtChanges(iStep, dt, dt * sigma,
                           "Change too large (" + std::to_string(eta) + " > " +
                               std::to_string(getMaxEta()) +
                               ").  Decreasing dt.");
        dt *= sigma;
      }
      else if (errorAbs > tsc.getMaxAbsError()) {  // reduce dt
        tsc.printDtChanges(
            iStep, dt, dt * sigma,
            "Absolute error is too large (" + std::to_string(errorAbs) + " > " +
                std::to_string(tsc.getMaxAbsError()) + ").  Decreasing dt.");
        dt *= sigma;
      }
      else if (errorRel > tsc.getMaxRelError()) {  // reduce dt
        tsc.printDtChanges(
            iStep, dt, dt * sigma,
            "Relative error is too large (" + std::to_string(errorRel) + " > " +
                std::to_string(tsc.getMaxRelError()) + ").  Decreasing dt.");
        dt *= sigma;
      }
    }

    if (dt < tsc.getMinTimeStep()) {  // decreased below minimum dt
      tsc.printDtChanges(iStep, dt, tsc.getMinTimeStep(),
                         "dt is too small.  Resetting to minimum dt.");
      dt = tsc.getMinTimeStep();
    }
    if (dt > tsc.getMaxTimeStep()) {  // increased above maximum dt
      tsc.printDtChanges(iStep, dt, tsc.getMaxTimeStep(),
                         "dt is too large.  Resetting to maximum dt.");
      dt = tsc.getMaxTimeStep();
    }

    workingState->setTimeStep(dt);
    workingState->setTime(solutionHistory->getCurrentState()->getTime() + dt);
    workingState->setComputeNorms(true);
  }

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override
  {
    return "Tempus::TimeStepControlStrategyBasicVS";
  }

  void describe(Teuchos::FancyOStream &out,
                const Teuchos::EVerbosityLevel verbLevel) const override
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, this->description());
    l_out->setOutputToRootOnly(0);

    *l_out << "\n--- " << this->description() << " ---" << std::endl;

    if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
      *l_out << "  StrategyType                      = "
             << this->getStrategyType() << std::endl
             << "  Step Type                         = " << this->getStepType()
             << std::endl
             << "  Amplification Factor              = " << getAmplFactor()
             << std::endl
             << "  Reduction Factor                  = " << getReductFactor()
             << std::endl
             << "  Minimum Value Monitoring Function = " << getMinEta()
             << std::endl
             << "  Maximum Value Monitoring Function = " << getMaxEta()
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

    pl->set<std::string>("Strategy Type", this->getStrategyType(), "Basic VS");
    pl->set<double>("Amplification Factor", getAmplFactor(),
                    "Amplification factor");
    pl->set<double>("Reduction Factor", getReductFactor(), "Reduction factor");
    pl->set<double>("Minimum Value Monitoring Function", getMinEta(),
                    "Min value eta");
    pl->set<double>("Maximum Value Monitoring Function", getMaxEta(),
                    "Max value eta");
    return pl;
  }

  virtual void initialize() const override
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        getAmplFactor() <= 1.0, std::out_of_range,
        "Error - Invalid value of Amplification Factor = "
            << getAmplFactor() << "!  \n"
            << "Amplification Factor must be > 1.0.\n");

    TEUCHOS_TEST_FOR_EXCEPTION(getReductFactor() >= 1.0, std::out_of_range,
                               "Error - Invalid value of Reduction Factor = "
                                   << getReductFactor() << "!  \n"
                                   << "Reduction Factor must be < 1.0.\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        getMinEta() > getMaxEta(), std::out_of_range,
        "Error - Invalid values of 'Minimum Value Monitoring Function' = "
            << getMinEta()
            << "\n and 'Maximum Value Monitoring Function' = " << getMaxEta()
            << "! \n Mininum Value cannot be > Maximum Value! \n");

    this->isInitialized_ = true;  // Only place where this is set to true!
  }

  virtual Scalar getAmplFactor() const { return rho_; }
  virtual Scalar getReductFactor() const { return sigma_; }
  virtual Scalar getMinEta() const { return minEta_; }
  virtual Scalar getMaxEta() const { return maxEta_; }

  virtual void setAmplFactor(Scalar rho)
  {
    rho_                 = rho;
    this->isInitialized_ = false;
  }
  virtual void setReductFactor(Scalar sigma)
  {
    sigma_               = sigma;
    this->isInitialized_ = false;
  }
  virtual void setMinEta(Scalar minEta)
  {
    minEta_              = minEta;
    this->isInitialized_ = false;
  }
  virtual void setMaxEta(Scalar maxEta)
  {
    maxEta_              = maxEta;
    this->isInitialized_ = false;
  }

 private:
  Scalar rho_;     ///< Amplification Factor
  Scalar sigma_;   ///< Reduction Factor
  Scalar minEta_;  ///< Minimum Value Monitoring Function
  Scalar maxEta_;  ///< Maximum Value Monitoring Function
};

/// Nonmember constructor.
template <class Scalar>
Teuchos::RCP<TimeStepControlStrategyBasicVS<Scalar> >
createTimeStepControlStrategyBasicVS(
    const Teuchos::RCP<Teuchos::ParameterList> &pList,
    std::string name = "Basic VS")
{
  auto tscs = Teuchos::rcp(new TimeStepControlStrategyBasicVS<Scalar>());
  if (pList == Teuchos::null || pList->numParams() == 0) return tscs;

  TEUCHOS_TEST_FOR_EXCEPTION(
      pList->get<std::string>("Strategy Type") != "Basic VS", std::logic_error,
      "Error - Strategy Type != 'Basic VS'.  (='" +
          pList->get<std::string>("Strategy Type") + "')\n");

  pList->validateParametersAndSetDefaults(*tscs->getValidParameters());

  tscs->setAmplFactor(pList->get<double>("Amplification Factor"));
  tscs->setReductFactor(pList->get<double>("Reduction Factor"));
  tscs->setMinEta(pList->get<double>("Minimum Value Monitoring Function"));
  tscs->setMaxEta(pList->get<double>("Maximum Value Monitoring Function"));

  tscs->setName(name);
  tscs->initialize();

  return tscs;
}

/// Nonmember function to return ParameterList with default values.
template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> getTimeStepControlStrategyBasicVS_PL()
{
  auto t = rcp(new Tempus::TimeStepControlStrategyBasicVS<Scalar>());
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      t->getValidParameters());
}

}  // namespace Tempus
#endif  // Tempus_TimeStepControlStrategy_BasicVS_hpp
